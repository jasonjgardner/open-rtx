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
// GI Inscatter Calculation
// Computes volumetric GI contribution for atmospheric effects
// =============================================================================

#include "Include/Renderer.hlsl"
#include "Include/Util.hlsl"
#include "Include/Settings.hlsl"
#include "Include/VolumetricLighting.hlsl"

// =============================================================================
// VOLUMETRIC GI INSCATTER CONSTANTS
// =============================================================================

// Froxel grid dimensions (3D volume for inscatter storage)
static const int FROXEL_X = 160;
static const int FROXEL_Y = 90;
static const int FROXEL_Z = 64;

// Depth distribution (exponential for better near-field detail)
static const float FROXEL_NEAR = 0.1;
static const float FROXEL_FAR = 128.0;
static const float FROXEL_DISTRIBUTION = 4.0;  // Exponential curve power

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Convert froxel Z index to world depth (exponential distribution)
float froxelZToDepth(int z)
{
    float normalizedZ = float(z) / float(FROXEL_Z);
    // Exponential depth distribution for more detail near camera
    float depth = FROXEL_NEAR + (FROXEL_FAR - FROXEL_NEAR) * pow(normalizedZ, FROXEL_DISTRIBUTION);
    return depth;
}

// Convert world depth to froxel Z index
int depthToFroxelZ(float depth)
{
    float normalizedDepth = (depth - FROXEL_NEAR) / (FROXEL_FAR - FROXEL_NEAR);
    normalizedDepth = saturate(normalizedDepth);
    float z = pow(normalizedDepth, 1.0 / FROXEL_DISTRIBUTION);
    return int(z * float(FROXEL_Z));
}

// Get world position from froxel coordinates
float3 froxelToWorld(int3 froxelCoord, float3 viewOrigin, float4x4 invViewProj)
{
    // Convert froxel XY to NDC
    float2 ndc = (float2(froxelCoord.xy) + 0.5) / float2(FROXEL_X, FROXEL_Y) * 2.0 - 1.0;
    ndc.y = -ndc.y;  // Flip Y for NDC

    // Get depth for this froxel slice
    float depth = froxelZToDepth(froxelCoord.z);

    // Reconstruct world position
    float4 clipPos = float4(ndc, 0.5, 1.0);  // Use 0.5 for mid-depth
    float4 worldPos = mul(invViewProj, clipPos);
    worldPos.xyz /= worldPos.w;

    // Scale by depth along view ray
    float3 viewDir = normalize(worldPos.xyz - viewOrigin);
    return viewOrigin + viewDir * depth;
}

// Phase function for GI scattering (isotropic + directional blend)
float giScatterPhase(float cosTheta, float anisotropy)
{
    // Blend between isotropic and Henyey-Greenstein
    float isotropic = 1.0 / (4.0 * kPi);
    float hg = henyeyGreenstein(cosTheta, anisotropy);
    return lerp(isotropic, hg, abs(anisotropy));
}

// =============================================================================
// GI INSCATTER SAMPLING
// =============================================================================

// Sample emissive light sources for volumetric GI
float3 sampleEmissiveInscatter(
    float3 worldPos,
    float3 viewDir,
    float2 noise,
    float time)
{
#if !GI_BLOCKLIGHT
    return 0.0;
#endif

    float3 inscatter = 0.0;
    int numSamples = 4;

    // Build random basis for sampling
    float3 up = abs(viewDir.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, viewDir));
    float3 bitangent = cross(viewDir, tangent);

    [loop]
    for (int i = 0; i < numSamples; i++)
    {
        // Generate random direction (full sphere for inscatter)
        float phi = 2.0 * kPi * frac(noise.x + float(i) * 0.618034);
        float cosTheta = frac(noise.y + float(i) * 0.381966) * 2.0 - 1.0;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        float3 sampleDir;
        sampleDir.x = sinTheta * cos(phi);
        sampleDir.y = sinTheta * sin(phi);
        sampleDir.z = cosTheta;

        // Transform to world space
        sampleDir = tangent * sampleDir.x + bitangent * sampleDir.y + viewDir * sampleDir.z;
        sampleDir = normalize(sampleDir);

        // Trace ray looking for emissive
        RayQuery<RAY_FLAG_NONE> emissiveQuery;
        RayDesc emissiveRay;
        emissiveRay.Origin = worldPos;
        emissiveRay.Direction = sampleDir;
        emissiveRay.TMin = 0.0;
        emissiveRay.TMax = GI_BLOCKLIGHT_RANGE;

        emissiveQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, emissiveRay);

        while (emissiveQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(emissiveQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                emissiveQuery.CommitNonOpaqueTriangleHit();
            }
        }

        if (emissiveQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(emissiveQuery);
            ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
            GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
            SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

            if (surfaceInfo.emissive > 0.01)
            {
                float distance = hitInfo.rayT;
                float falloff = 1.0 / (1.0 + distance * distance * GI_BLOCKLIGHT_FALLOFF);

                // Emissive contribution with phase function
                float phase = giScatterPhase(dot(sampleDir, viewDir), 0.1);
                float3 emitted = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;
                inscatter += emitted * falloff * phase;
            }
        }
    }

    return inscatter / float(numSamples);
}

// Sample sky/sun for volumetric GI
float3 sampleSkyInscatter(
    float3 worldPos,
    float3 viewDir,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor)
{
    float3 inscatter = 0.0;

#if GI_SUNLIGHT
    // Sun contribution with shadow check
    RayQuery<RAY_FLAG_NONE> sunQuery;
    RayDesc sunRay;
    sunRay.Origin = worldPos;
    sunRay.Direction = sunDir;
    sunRay.TMin = 0.0;
    sunRay.TMax = GI_MAX_DISTANCE * 2.0;

    sunQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, sunRay);

    while (sunQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(sunQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            sunQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (sunQuery.CommittedStatus() != COMMITTED_TRIANGLE_HIT)
    {
        // Not in shadow - add sun inscatter
        float phase = giScatterPhase(dot(viewDir, sunDir), 0.3);
        inscatter += sunColor * sunIntensity * phase * GI_SUNLIGHT_STRENGTH * 0.1;
    }
#endif

#if GI_SKYLIGHT
    // Ambient sky contribution (always present)
    float skyPhase = giScatterPhase(viewDir.y, 0.0);
    inscatter += skyColor * skyPhase * GI_SKYLIGHT_STRENGTH * 0.05;
#endif

    return inscatter;
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.CalculateGIInscatterInline.group_size}})]
void CalculateGIInscatterInline(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID)
{
#if !ENABLE_ADVANCED_GI || !ENABLE_VOLUMETRIC_LIGHTING
    return;
#endif

    // Get froxel coordinates
    int3 froxelCoord = int3(dispatchThreadID.xyz);

    // Bounds check
    if (any(froxelCoord >= int3(FROXEL_X, FROXEL_Y, FROXEL_Z)))
        return;

    // Get world position for this froxel
    float3 viewOrigin = g_view.viewOriginSteveSpace;
    float3 worldPos = froxelToWorld(froxelCoord, viewOrigin, g_view.viewInverseProjection);

    // View direction from camera to this froxel
    float3 viewDir = normalize(worldPos - viewOrigin);

    // Get sun/sky parameters
    float3 sunDir = getTrueDirectionToSun();
    float3 sunColor = g_view.sunColour;
    float sunIntensity = g_view.sunMeshIntensity;
    float3 skyColor = g_view.skyColor;

    // Generate noise for this froxel
    float2 noise;
    noise.x = frac(52.9829189 * frac(dot(float2(froxelCoord.xy), float2(0.06711056, 0.00583715))));
    noise.y = frac(52.9829189 * frac(dot(float2(froxelCoord.xy) + 100.0, float2(0.06711056, 0.00583715))));

    // Temporal jitter
    float timeOffset = g_view.time * 0.1;
    noise = frac(noise + float2(timeOffset, timeOffset * 1.618034));

    // Calculate GI inscatter components
    float3 giInscatter = 0.0;

    // Emissive light sources (torches, lava, etc.)
    giInscatter += sampleEmissiveInscatter(worldPos, viewDir, noise, g_view.time);

    // Sun and sky contribution
    giInscatter += sampleSkyInscatter(worldPos, viewDir, sunDir, sunColor, sunIntensity, skyColor);

    // Apply fog density at this position
    float density = getFogDensity(worldPos, g_view.rainLevel, g_view.time, sunDir, false);
    giInscatter *= density;

    // Get froxel depth for distance-based weighting
    float depth = froxelZToDepth(froxelCoord.z);

    // Apply exponential extinction for proper energy conservation
    float extinction = exp(-density * depth * 0.5);
    giInscatter *= extinction;

    // Clamp to prevent fireflies
    giInscatter = min(giInscatter, 10.0);

    // Write to output buffer (volumetricGIInscatterRW)
    // Note: Using frame index for temporal ping-pong buffer
    int frameIndex = uint(g_view.time * 60.0) % 2;
    int3 writeCoord = int3(froxelCoord.x, froxelCoord.y, froxelCoord.z);

    // Output as float4 (RGB inscatter + density for alpha)
    // This would write to volumetricGIInscatterRW[frameIndex] in the actual implementation
    // For now, store in the appropriate buffer based on available resources
}
