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
// OpenRTX Global Illumination System
// Hardware-accelerated raytraced global illumination using DXR
// Uses deterministic ray directions and single-bounce indirect lighting
// =============================================================================

#ifndef __OPENRTX_GI_HLSL__
#define __OPENRTX_GI_HLSL__

#include "Settings.hlsl"

// =============================================================================
// GI STRUCTURES
// =============================================================================

// Material properties for GI calculations
struct GIMaterial
{
    float3 albedo;          // Base color
    float roughness;        // Surface roughness (0 = mirror, 1 = diffuse)
    float metalness;        // Metalness (0 = dielectric, 1 = metal)
    float3 emission;        // Emissive light
    float subsurface;       // Subsurface scattering amount
};

// Accumulated GI result
struct GIResult
{
    float3 diffuse;         // Diffuse indirect lighting
    float3 specular;        // Specular indirect lighting (glossy reflections)
    float3 emission;        // Accumulated emissive light from nearby sources
    float ao;               // Raytraced ambient occlusion factor
    float visibility;       // Overall visibility (for shadows)
};

// Raytraced hit information
struct RTHitInfo
{
    bool hit;               // Did the ray hit something?
    float distance;         // Distance to hit
    float3 position;        // World position of hit
    float3 normal;          // Surface normal at hit
    float3 albedo;          // Surface albedo color
    float3 emission;        // Emissive light from surface
};

// =============================================================================
// DETERMINISTIC RAY DIRECTIONS
// =============================================================================

// Pre-computed hemisphere sample directions for consistent raytracing
// These are distributed using Fibonacci spiral for uniform coverage
static const float3 kHemisphereDirs[16] = {
    float3(0.0000, 1.0000, 0.0000),   // Up
    float3(0.5257, 0.8507, 0.0000),   // Upper ring
    float3(0.1625, 0.8507, 0.5000),
    float3(-0.4253, 0.8507, 0.3090),
    float3(-0.4253, 0.8507, -0.3090),
    float3(0.1625, 0.8507, -0.5000),
    float3(0.6882, 0.5257, 0.5000),   // Middle ring
    float3(-0.2629, 0.5257, 0.8090),
    float3(-0.8507, 0.5257, 0.0000),
    float3(-0.2629, 0.5257, -0.8090),
    float3(0.6882, 0.5257, -0.5000),
    float3(0.5878, 0.3090, 0.7500),   // Lower ring
    float3(-0.4755, 0.3090, 0.8236),
    float3(-0.9511, 0.3090, 0.0000),
    float3(-0.4755, 0.3090, -0.8236),
    float3(0.5878, 0.3090, -0.7500)
};

// Transform hemisphere direction to align with surface normal
float3 alignToNormal(float3 dir, float3 normal)
{
    // Build tangent space from normal
    float3 up = abs(normal.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, normal));
    float3 bitangent = cross(normal, tangent);

    // Transform direction from Y-up to normal-aligned
    return normalize(tangent * dir.x + normal * dir.y + bitangent * dir.z);
}

// Get deterministic ray direction for sample index
float3 getGIRayDirection(float3 normal, int sampleIndex, int totalSamples)
{
    int idx = sampleIndex % 16;
    float3 localDir = kHemisphereDirs[idx];
    return alignToNormal(localDir, normal);
}

// =============================================================================
// TEMPORAL JITTER (for noise reduction across frames)
// =============================================================================

// Get frame-based jitter offset for temporal accumulation
float2 getTemporalJitter(float time)
{
    // 8-frame Halton sequence rotation
    int frame = int(time * 60.0) % 8;
    static const float2 kHaltonSequence[8] = {
        float2(0.5, 0.333),
        float2(0.25, 0.666),
        float2(0.75, 0.111),
        float2(0.125, 0.444),
        float2(0.625, 0.777),
        float2(0.375, 0.222),
        float2(0.875, 0.555),
        float2(0.0625, 0.888)
    };
    return kHaltonSequence[frame];
}

// R2 quasi-random sequence for sample distribution
float2 getR2Sequence(int n)
{
    const float g = 1.32471795724474602596;
    const float a1 = 1.0 / g;
    const float a2 = 1.0 / (g * g);
    return frac(float2(a1, a2) * float(n) + 0.5);
}

// Initialize empty GI result
GIResult initGIResult()
{
    GIResult result;
    result.diffuse = 0.0;
    result.specular = 0.0;
    result.emission = 0.0;
    result.ao = 1.0;
    result.visibility = 1.0;
    return result;
}

// =============================================================================
// RAYTRACED AMBIENT OCCLUSION (RTAO)
// =============================================================================

// Trace a single occlusion ray using hardware RT
float traceOcclusionRay(float3 origin, float3 direction, float maxDistance)
{
    RayQuery<RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH> aoQuery;
    RayDesc aoRay;
    aoRay.Origin = origin;
    aoRay.Direction = direction;
    aoRay.TMin = 0.001;
    aoRay.TMax = maxDistance;

    aoQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, aoRay);

    while (aoQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(aoQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            aoQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (aoQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        float hitDist = aoQuery.CommittedRayT();
        // Smooth distance-weighted occlusion
        return saturate(hitDist / maxDistance);
    }

    return 1.0;  // No occlusion
}

// Calculate raytraced ambient occlusion using deterministic ray directions
float calculateRTAO(float3 position, float3 normal, int samples, float time)
{
#if !GI_AMBIENT_OCCLUSION
    return 1.0;
#endif

    float visibility = 0.0;
    float3 biasedPos = position + normal * 0.02;

    // Use temporal rotation for sample distribution
    int frameOffset = int(time * 60.0) % 4;

    [unroll]
    for (int i = 0; i < samples && i < 8; i++)
    {
        int sampleIdx = (i + frameOffset) % 16;
        float3 rayDir = getGIRayDirection(normal, sampleIdx, samples);

        // Trace occlusion ray
        float sampleVisibility = traceOcclusionRay(biasedPos, rayDir, GI_AO_RADIUS);
        visibility += sampleVisibility;
    }

    visibility /= float(min(samples, 8));
    return lerp(1.0, visibility, GI_AO_STRENGTH);
}

// =============================================================================
// RAYTRACED INDIRECT LIGHTING
// =============================================================================

// Trace a single indirect ray and return hit information
RTHitInfo traceIndirectRay(float3 origin, float3 direction, float maxDistance)
{
    RTHitInfo result;
    result.hit = false;
    result.distance = maxDistance;
    result.position = origin + direction * maxDistance;
    result.normal = -direction;
    result.albedo = 0.0;
    result.emission = 0.0;

    RayQuery<RAY_FLAG_NONE> indirectQuery;
    RayDesc indirectRay;
    indirectRay.Origin = origin;
    indirectRay.Direction = direction;
    indirectRay.TMin = 0.001;
    indirectRay.TMax = maxDistance;

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

        result.hit = true;
        result.distance = hitInfo.rayT;
        result.position = surfaceInfo.position;
        result.normal = surfaceInfo.normal;
        result.albedo = surfaceInfo.color;
        result.emission = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;
    }

    return result;
}

// Check sun visibility from a point (shadow ray)
float traceSunShadow(float3 position, float3 sunDir)
{
    RayQuery<RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH> shadowQuery;
    RayDesc shadowRay;
    shadowRay.Origin = position;
    shadowRay.Direction = sunDir;
    shadowRay.TMin = 0.001;
    shadowRay.TMax = GI_MAX_DISTANCE * 2.0;

    shadowQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, shadowRay);

    while (shadowQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(shadowQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            shadowQuery.CommitNonOpaqueTriangleHit();
        }
    }

    return (shadowQuery.CommittedStatus() != COMMITTED_TRIANGLE_HIT) ? 1.0 : 0.0;
}

// =============================================================================
// SINGLE-BOUNCE INDIRECT LIGHTING
// =============================================================================

// Calculate indirect diffuse lighting using raytracing
// Uses single bounce for performance, with direct lighting at hit points
float3 calculateIndirectDiffuse(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    int samples,
    float time)
{
#if !GI_DIFFUSE_BOUNCES
    return 0.0;
#endif

    float3 indirectLight = 0.0;
    float3 biasedPos = position + normal * 0.02;

    // Temporal rotation for better coverage
    int frameOffset = int(time * 60.0) % 8;

    [unroll]
    for (int i = 0; i < samples && i < 8; i++)
    {
        int sampleIdx = (i * 2 + frameOffset) % 16;
        float3 rayDir = getGIRayDirection(normal, sampleIdx, samples);

        // Trace indirect ray
        RTHitInfo hit = traceIndirectRay(biasedPos, rayDir, GI_MAX_DISTANCE);

        float3 sampleLight = 0.0;

        if (hit.hit)
        {
            // Calculate lighting at hit point
            float hitNdotL = saturate(dot(hit.normal, sunDir));

            // Check sun visibility at hit point
            float sunShadow = traceSunShadow(hit.position + hit.normal * 0.02, sunDir);

            // Direct sun contribution at hit surface
            sampleLight += hit.albedo * sunColor * sunIntensity * hitNdotL * sunShadow;

            // Emissive contribution (light sources)
            sampleLight += hit.emission * INDIRECT_EMISSIVE_BOOST;

            // Ambient sky approximation at hit
            float skyVisibility = saturate(hit.normal.y * 0.5 + 0.5);
            sampleLight += hit.albedo * skyColor * skyVisibility * 0.2;

            // Distance falloff
            float falloff = 1.0 / (1.0 + hit.distance * hit.distance * GI_DISTANCE_FALLOFF);
            sampleLight *= falloff;
        }
        else
        {
            // Hit sky - sample sky color
            float skyBrightness = saturate(rayDir.y * 0.5 + 0.5);
            sampleLight = skyColor * skyBrightness * 0.3;
        }

        indirectLight += sampleLight;
    }

    indirectLight /= float(min(samples, 8));

    // Apply surface albedo (diffuse transfer)
    return albedo * indirectLight * GI_DIFFUSE_STRENGTH;
}

// =============================================================================
// EMISSIVE LIGHT GATHERING
// =============================================================================

// Gather light from nearby emissive surfaces using raytracing
float3 gatherEmissiveLight(
    float3 position,
    float3 normal,
    float3 albedo,
    int samples,
    float time)
{
#if !GI_BLOCKLIGHT
    return 0.0;
#endif

    float3 emissiveLight = 0.0;
    float3 biasedPos = position + normal * 0.02;

    // Temporal rotation
    int frameOffset = int(time * 60.0) % 4;

    [unroll]
    for (int i = 0; i < samples && i < 8; i++)
    {
        int sampleIdx = (i + frameOffset * 2) % 16;
        float3 rayDir = getGIRayDirection(normal, sampleIdx, samples);

        // Trace ray looking for emissive surfaces
        RTHitInfo hit = traceIndirectRay(biasedPos, rayDir, GI_BLOCKLIGHT_RANGE);

        if (hit.hit)
        {
            // Check for emissive
            float emissiveStrength = length(hit.emission);
            if (emissiveStrength > 0.01)
            {
                // Distance-based falloff
                float falloff = 1.0 / (1.0 + hit.distance * hit.distance * GI_BLOCKLIGHT_FALLOFF);
                emissiveLight += hit.emission * falloff;
            }
        }
    }

    emissiveLight /= float(min(samples, 8));

    return albedo * emissiveLight * GI_BLOCKLIGHT_STRENGTH;
}

// =============================================================================
// SKY VISIBILITY
// =============================================================================

// Calculate sky visibility factor using raytracing
float3 calculateSkyVisibility(
    float3 position,
    float3 normal,
    float3 skyColor,
    int samples,
    float time)
{
#if !GI_SKYLIGHT
    return skyColor * 0.5;  // Fallback ambient
#endif

    float visibleSamples = 0.0;
    float3 biasedPos = position + normal * 0.02;

    // Check upward hemisphere for sky visibility
    int frameOffset = int(time * 60.0) % 4;

    [unroll]
    for (int i = 0; i < samples && i < 4; i++)
    {
        int sampleIdx = (i + frameOffset) % 6;  // Use upper hemisphere samples
        float3 rayDir = getGIRayDirection(normal, sampleIdx, samples);

        // Only check if direction has positive sky component
        if (rayDir.y > 0.1)
        {
            float vis = traceOcclusionRay(biasedPos, rayDir, GI_MAX_DISTANCE);
            visibleSamples += vis;
        }
        else
        {
            visibleSamples += 0.5;  // Partial visibility for side directions
        }
    }

    float skyFactor = visibleSamples / float(min(samples, 4));
    return skyColor * skyFactor * GI_SKYLIGHT_STRENGTH;
}

// =============================================================================
// RAYTRACED SPECULAR (Single-bounce glossy reflections)
// =============================================================================

// Calculate specular indirect using raytracing with reflection direction
float3 calculateSpecularIndirect(
    float3 position,
    float3 normal,
    float3 viewDir,
    float3 albedo,
    float roughness,
    float metalness,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float time)
{
#if !GI_SPECULAR
    return 0.0;
#endif

    // Skip very rough surfaces
    if (roughness > GI_SPECULAR_MAX_ROUGHNESS)
        return 0.0;

    // Calculate reflection direction
    float3 reflectDir = reflect(-viewDir, normal);

    // Add roughness-based cone spread using temporal jitter
    if (roughness > 0.05)
    {
        float2 jitter = getTemporalJitter(time) * 2.0 - 1.0;
        float3 tangent = normalize(cross(normal, abs(normal.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0)));
        float3 bitangent = cross(normal, tangent);

        float coneAngle = roughness * roughness * 0.5;
        reflectDir = normalize(reflectDir + (tangent * jitter.x + bitangent * jitter.y) * coneAngle);

        // Ensure reflection stays above surface
        if (dot(reflectDir, normal) < 0.0)
            reflectDir = reflect(-viewDir, normal);
    }

    float3 biasedPos = position + normal * 0.02;

    // Trace reflection ray
    RTHitInfo hit = traceIndirectRay(biasedPos, reflectDir, GI_MAX_DISTANCE);

    float3 specularLight = 0.0;

    if (hit.hit)
    {
        // Calculate lighting at hit point
        float hitNdotL = saturate(dot(hit.normal, sunDir));
        float sunShadow = traceSunShadow(hit.position + hit.normal * 0.02, sunDir);

        specularLight += hit.albedo * sunColor * sunIntensity * hitNdotL * sunShadow;
        specularLight += hit.emission;
        specularLight += hit.albedo * skyColor * 0.15;
    }
    else
    {
        // Hit sky
        float skyBrightness = saturate(reflectDir.y * 0.5 + 0.5);
        specularLight = skyColor * skyBrightness;
    }

    // Apply Fresnel
    float3 f0 = lerp(float3(0.04, 0.04, 0.04), albedo, metalness);
    float NdotV = saturate(dot(normal, viewDir));
    float3 fresnel = f0 + (1.0 - f0) * pow(1.0 - NdotV, 5.0);

    // Roughness fade
    float roughnessFade = 1.0 - saturate(roughness / GI_SPECULAR_MAX_ROUGHNESS);

    return specularLight * fresnel * GI_SPECULAR_STRENGTH * roughnessFade;
}

// =============================================================================
// MAIN RAYTRACED GI CALCULATION
// =============================================================================

// Calculate complete raytraced global illumination for a surface
// Uses hardware-accelerated raytracing with deterministic directions
GIResult calculateGI(
    float3 position,
    float3 normal,
    float3 viewDir,
    GIMaterial material,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float2 pixelCoord,
    float time,
    int quality)  // 0 = low, 1 = medium, 2 = high
{
    GIResult result = initGIResult();

#if !ENABLE_ADVANCED_GI
    return result;
#endif

    // Adjust sample counts based on quality
    int aoSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_AO_SAMPLES);
    int diffuseSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_DIFFUSE_SAMPLES);
    int skylightSamples = quality == 0 ? 1 : (quality == 1 ? 2 : GI_SKYLIGHT_SAMPLES);
    int blocklightSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_BLOCKLIGHT_SAMPLES);

    // Raytraced Ambient Occlusion
#if GI_AMBIENT_OCCLUSION
    result.ao = calculateRTAO(position, normal, aoSamples, time);
#endif

    // Sky visibility (raytraced)
#if GI_SKYLIGHT
    float3 skyGI = calculateSkyVisibility(position, normal, skyColor, skylightSamples, time);
    result.diffuse += material.albedo * skyGI;
#endif

    // Emissive light gathering (raytraced)
#if GI_BLOCKLIGHT
    float3 emissiveGI = gatherEmissiveLight(position, normal, material.albedo, blocklightSamples, time);
    result.emission += emissiveGI;
#endif

    // Single-bounce indirect diffuse (raytraced)
#if GI_DIFFUSE_BOUNCES
    float3 indirectDiffuse = calculateIndirectDiffuse(
        position, normal, material.albedo,
        sunDir, sunColor, sunIntensity, skyColor,
        diffuseSamples, time);
    result.diffuse += indirectDiffuse;
#endif

    // Specular indirect (raytraced reflection)
#if GI_SPECULAR
    if (material.roughness < GI_SPECULAR_MAX_ROUGHNESS)
    {
        result.specular = calculateSpecularIndirect(
            position, normal, viewDir,
            material.albedo, material.roughness, material.metalness,
            sunDir, sunColor, sunIntensity, skyColor,
            time);
    }
#endif

    // Apply AO to indirect lighting
    result.diffuse *= result.ao;
    result.emission *= lerp(1.0, result.ao, 0.5);

    return result;
}

// Simplified GI for shadow rays (just visibility)
float calculateGIShadow(float3 position, float3 direction, float maxDistance)
{
    RayQuery<RAY_FLAG_NONE> shadowQuery;
    RayDesc shadowRay;
    shadowRay.Origin = position;
    shadowRay.Direction = direction;
    shadowRay.TMin = 0.001;
    shadowRay.TMax = maxDistance;

    shadowQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, shadowRay);

    while (shadowQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(shadowQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            shadowQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (shadowQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        return 0.0;  // Occluded
    }

    return 1.0;  // Visible
}

// =============================================================================
// GI INTEGRATION HELPERS
// =============================================================================

// Apply GI result to surface lighting
float3 applyGI(float3 directLighting, GIResult gi, float giIntensity)
{
    float3 indirect = gi.diffuse + gi.specular + gi.emission;
    return directLighting + indirect * giIntensity;
}

// Combine GI with existing emissive sampling (for backwards compatibility)
float3 combineWithEmissiveGI(float3 existingEmissiveGI, GIResult gi)
{
    // Blend new GI system with existing emissive GI
    // Weight towards new system when advanced GI is enabled
#if ENABLE_ADVANCED_GI
    return gi.emission + gi.diffuse + existingEmissiveGI * 0.3;
#else
    return existingEmissiveGI;
#endif
}

#endif // __OPENRTX_GI_HLSL__
