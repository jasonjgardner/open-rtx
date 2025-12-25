/* MIT License
 *
 * Copyright (c) 2025 OpenRTX Contributors
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// =============================================================================
// Diffuse Ray Generation for Raytraced Global Illumination
// Uses hardware-accelerated raytracing with deterministic directions
// =============================================================================

#include "Include/Renderer.hlsl"
#include "Include/Util.hlsl"
#include "Include/OpenRTX.hlsl"
#include "Include/GI.hlsl"

// =============================================================================
// RAYTRACED DIFFUSE INDIRECT
// =============================================================================

// Trace single-bounce indirect ray and compute lighting at hit
float3 traceIndirectBounce(
    float3 origin,
    float3 direction,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float3 ambientColor)
{
    float3 radiance = 0.0;

    // Trace indirect ray using hardware RT
    RayQuery<RAY_FLAG_NONE> indirectQuery;
    RayDesc indirectRay;
    indirectRay.Origin = origin;
    indirectRay.Direction = direction;
    indirectRay.TMin = 0.001;
    indirectRay.TMax = GI_MAX_DISTANCE;

    indirectQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, indirectRay);

    while (indirectQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(indirectQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            indirectQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (indirectQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        HitInfo hitInfo = GetCommittedHitInfo(indirectQuery);
        ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
        GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
        SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

        float hitDistance = hitInfo.rayT;

        // Calculate direct lighting at hit point
        float NdotL = saturate(dot(surfaceInfo.normal, sunDir));

        // Trace shadow ray at hit point
        float sunVisibility = 1.0;
        if (NdotL > 0.0)
        {
            sunVisibility = traceSunShadow(surfaceInfo.position + surfaceInfo.normal * 0.02, sunDir);
        }

        // Direct sun at hit
        radiance += surfaceInfo.color * sunColor * sunIntensity * NdotL * sunVisibility;

        // Ambient at hit
        radiance += surfaceInfo.color * ambientColor;

        // Emissive at hit (boosted for GI)
        radiance += surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY * INDIRECT_EMISSIVE_BOOST;

        // Sky approximation at hit
        float skyVis = saturate(surfaceInfo.normal.y * 0.5 + 0.5);
        radiance += surfaceInfo.color * skyColor * skyVis * 0.2;

        // Distance falloff
        float falloff = 1.0 / (1.0 + hitDistance * hitDistance * GI_DISTANCE_FALLOFF);
        radiance *= falloff;
    }
    else
    {
        // Hit sky
        float skyBrightness = saturate(direction.y * 0.5 + 0.5);
        radiance = skyColor * skyBrightness * 0.5;
    }

    return radiance;
}

// Calculate raytraced diffuse GI using deterministic sample directions
float3 calculateRaytracedDiffuseGI(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float3 ambientColor,
    int numSamples,
    float time)
{
    float3 indirectLight = 0.0;
    float3 biasedPos = position + normal * 0.02;

    // Temporal offset for sample rotation
    int frameOffset = int(time * 60.0) % 8;

    [unroll]
    for (int i = 0; i < numSamples && i < 8; i++)
    {
        // Use deterministic hemisphere directions with temporal rotation
        int sampleIdx = (i * 2 + frameOffset) % 16;
        float3 rayDir = getGIRayDirection(normal, sampleIdx, numSamples);

        // Trace indirect bounce
        float3 bounceLight = traceIndirectBounce(
            biasedPos, rayDir,
            sunDir, sunColor, sunIntensity,
            skyColor, ambientColor);

        indirectLight += bounceLight;
    }

    indirectLight /= float(min(numSamples, 8));

    // Apply surface albedo (diffuse transfer)
    return albedo * indirectLight;
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.DiffuseRayGenCombinedInline.group_size}})]
void DiffuseRayGenCombinedInline(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID)
{
#if !ENABLE_ADVANCED_GI || !GI_DIFFUSE_BOUNCES
    return;
#endif

    uint2 pixelCoord = dispatchThreadID.xy;

    // Bounds check
    if (any(pixelCoord >= g_view.renderResolution))
        return;

    // Trace primary ray to get surface data
    RayDesc primaryRay;
    primaryRay.Direction = rayDirFromNDC(getNDCjittered(pixelCoord));
    primaryRay.Origin = g_view.viewOriginSteveSpace;
    primaryRay.TMin = 0;
    primaryRay.TMax = 10000;

    RayQuery<RAY_FLAG_NONE> primaryQuery;
    primaryQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, primaryRay);

    while (primaryQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(primaryQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            primaryQuery.CommitNonOpaqueTriangleHit();
        }
    }

    // Only process valid surface hits
    if (primaryQuery.CommittedStatus() != COMMITTED_TRIANGLE_HIT)
        return;

    HitInfo hitInfo = GetCommittedHitInfo(primaryQuery);
    ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
    GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
    SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

    // Skip emissive surfaces (they provide GI, don't receive it)
    if (surfaceInfo.emissive > 0.5)
        return;

    // Skip highly metallic surfaces (specular only)
    if (surfaceInfo.metalness > 0.9)
        return;

    // Get lighting parameters
    float3 sunDir = getTrueDirectionToSun();
    float3 sunColor = g_view.sunColour;
    float sunIntensity = g_view.sunMeshIntensity;
    float3 skyColor = g_view.skyColor;
    float3 ambientColor = g_view.constantAmbient;

    // Determine sample count based on distance (LOD)
    float distanceFromCamera = hitInfo.rayT;
    float distanceLOD = saturate(distanceFromCamera / 128.0);

    // Reduce samples at distance for performance
    int numSamples = int(lerp(float(GI_DIFFUSE_SAMPLES), 1.0, distanceLOD * distanceLOD));
    numSamples = max(1, numSamples);

    // Calculate raytraced diffuse GI
    float3 diffuseGI = calculateRaytracedDiffuseGI(
        surfaceInfo.position,
        surfaceInfo.normal,
        surfaceInfo.color,
        sunDir, sunColor, sunIntensity,
        skyColor, ambientColor,
        numSamples,
        g_view.time);

    // Apply strength
    diffuseGI *= GI_DIFFUSE_STRENGTH;

    // Clamp to prevent fireflies
    diffuseGI = min(diffuseGI, 10.0);

    // Output would go to diffuse GI buffer for denoising
    // denoisingInputs[bufferIndex][pixelCoord] = float4(diffuseGI, 1.0);
}
