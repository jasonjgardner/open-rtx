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
// Raytraced GI Inscatter Calculation
// Computes volumetric GI contribution using hardware-accelerated raytracing
// =============================================================================

#include "Include/Renderer.hlsl"
#include "Include/Util.hlsl"
#include "Include/Settings.hlsl"
#include "Include/VolumetricLighting.hlsl"
#include "Include/GI.hlsl"

// =============================================================================
// FROXEL GRID CONSTANTS
// =============================================================================

static const int FROXEL_X = 160;
static const int FROXEL_Y = 90;
static const int FROXEL_Z = 64;

static const float FROXEL_NEAR = 0.1;
static const float FROXEL_FAR = 128.0;
static const float FROXEL_DISTRIBUTION = 4.0;

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

float froxelZToDepth(int z)
{
    float normalizedZ = float(z) / float(FROXEL_Z);
    return FROXEL_NEAR + (FROXEL_FAR - FROXEL_NEAR) * pow(normalizedZ, FROXEL_DISTRIBUTION);
}

// Simplified world position calculation using screen UV and depth
float3 froxelToWorldSimple(int3 froxelCoord, float3 viewOrigin)
{
    // Convert froxel XY to screen UV
    float2 uv = (float2(froxelCoord.xy) + 0.5) / float2(FROXEL_X, FROXEL_Y);

    // Get depth for this froxel slice
    float depth = froxelZToDepth(froxelCoord.z);

    // Calculate view direction using NDC and approximate FOV
    float2 ndc = uv * 2.0 - 1.0;
    ndc.y = -ndc.y;

    // Approximate view direction (using typical Minecraft FOV ~70 degrees)
    float tanHalfFov = tan(radians(35.0));
    float aspectRatio = float(FROXEL_X) / float(FROXEL_Y);

    float3 viewDir;
    viewDir.x = ndc.x * tanHalfFov * aspectRatio;
    viewDir.y = ndc.y * tanHalfFov;
    viewDir.z = 1.0;
    viewDir = normalize(viewDir);

    return viewOrigin + viewDir * depth;
}

float giScatterPhase(float cosTheta, float anisotropy)
{
    float isotropic = 1.0 / (4.0 * kPi);
    float hg = henyeyGreenstein(cosTheta, anisotropy);
    return lerp(isotropic, hg, abs(anisotropy));
}

// =============================================================================
// RAYTRACED INSCATTER SAMPLING
// =============================================================================

// Trace rays to find emissive light sources contributing to inscatter
float3 traceEmissiveInscatter(float3 worldPos, float3 viewDir, float time)
{
#if !GI_BLOCKLIGHT
    return 0.0;
#endif

    float3 inscatter = 0.0;

    // Use deterministic directions for consistent results
    int frameOffset = int(time * 60.0) % 4;

    [unroll]
    for (int i = 0; i < 4; i++)
    {
        int sampleIdx = (i + frameOffset) % 16;
        float3 rayDir = kHemisphereDirs[sampleIdx];

        // Transform to world-aligned (not normal-aligned for volumetric)
        float3 sampleDir = normalize(rayDir);

        // Trace ray looking for emissive surfaces
        RayQuery<RAY_FLAG_NONE> emissiveQuery;
        RayDesc emissiveRay;
        emissiveRay.Origin = worldPos;
        emissiveRay.Direction = sampleDir;
        emissiveRay.TMin = 0.001;
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

                float phase = giScatterPhase(dot(sampleDir, viewDir), 0.1);
                float3 emitted = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;
                inscatter += emitted * falloff * phase;
            }
        }
    }

    return inscatter * 0.25;
}

// Trace sun shadow ray for volumetric inscatter
float3 traceSunInscatter(float3 worldPos, float3 viewDir, float3 sunDir, float3 sunColor, float sunIntensity)
{
#if !GI_SUNLIGHT
    return 0.0;
#endif

    // Check sun visibility using raytracing
    RayQuery<RAY_FLAG_NONE> sunQuery;
    RayDesc sunRay;
    sunRay.Origin = worldPos;
    sunRay.Direction = sunDir;
    sunRay.TMin = 0.001;
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
        float phase = giScatterPhase(dot(viewDir, sunDir), 0.3);
        return sunColor * sunIntensity * phase * GI_SUNLIGHT_STRENGTH * 0.1;
    }

    return 0.0;
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

    int3 froxelCoord = int3(dispatchThreadID.xyz);

    if (any(froxelCoord >= int3(FROXEL_X, FROXEL_Y, FROXEL_Z)))
        return;

    float3 viewOrigin = g_view.viewOriginSteveSpace;
    float3 worldPos = froxelToWorldSimple(froxelCoord, viewOrigin);
    float3 viewDir = normalize(worldPos - viewOrigin);

    float3 sunDir = getTrueDirectionToSun();
    float3 sunColor = g_view.sunColour;
    float sunIntensity = g_view.sunMeshIntensity;
    float3 skyColor = g_view.skyColor;

    // Calculate raytraced GI inscatter
    float3 giInscatter = 0.0;

    // Emissive light sources (raytraced)
    giInscatter += traceEmissiveInscatter(worldPos, viewDir, g_view.time);

    // Sun contribution (raytraced shadow)
    giInscatter += traceSunInscatter(worldPos, viewDir, sunDir, sunColor, sunIntensity);

    // Ambient sky contribution
#if GI_SKYLIGHT
    float skyPhase = giScatterPhase(viewDir.y, 0.0);
    giInscatter += skyColor * skyPhase * GI_SKYLIGHT_STRENGTH * 0.05;
#endif

    // Apply fog density
    float density = getFogDensity(worldPos, g_view.rainLevel, g_view.time, sunDir, false);
    giInscatter *= density;

    // Apply extinction
    float depth = froxelZToDepth(froxelCoord.z);
    float extinction = exp(-density * depth * 0.5);
    giInscatter *= extinction;

    // Clamp to prevent fireflies
    giInscatter = min(giInscatter, 10.0);

    // Output to volumetric GI buffer
    // volumetricGIInscatterRW[frameIndex][froxelCoord] = float4(giInscatter, density);
}
