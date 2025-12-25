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
// Path-traced global illumination with multi-bounce support
// Inspired by SEUS PTGI techniques
// =============================================================================

#ifndef __OPENRTX_GI_HLSL__
#define __OPENRTX_GI_HLSL__

#include "Settings.hlsl"

// =============================================================================
// GI STRUCTURES
// =============================================================================

// Ray traversal data for efficient voxel-based ray marching
struct GIRayData
{
    float3 origin;          // Ray origin in world space
    float3 direction;       // Normalized ray direction
    float3 invDirection;    // 1.0 / direction for efficient intersection
    float maxDistance;      // Maximum trace distance
};

// GI hit result containing illumination data
struct GIHitResult
{
    bool hit;               // Did the ray hit something?
    float distance;         // Distance to hit
    float3 position;        // World position of hit
    float3 normal;          // Surface normal at hit
    float3 albedo;          // Surface albedo color
    float3 emission;        // Emissive light from surface
    float roughness;        // Surface roughness
    float metalness;        // Surface metalness
    float ao;               // Ambient occlusion at hit
};

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
    float3 emission;        // Accumulated emissive light
    float ao;               // Screen-space AO factor
    float visibility;       // Overall visibility (for shadows)
};

// =============================================================================
// GI UTILITY FUNCTIONS
// =============================================================================

// Initialize GI ray from origin and direction
GIRayData initGIRay(float3 origin, float3 direction, float maxDist)
{
    GIRayData ray;
    ray.origin = origin;
    ray.direction = normalize(direction);
    ray.invDirection = 1.0 / (ray.direction + sign(ray.direction) * 1e-7);
    ray.maxDistance = maxDist;
    return ray;
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

// Blue noise sampling for temporal stability
// Based on interleaved gradient noise with temporal offset
float3 getGIBlueNoise(float2 pixelCoord, float time, int sampleIndex)
{
    // Interleaved gradient noise (fast, good distribution)
    float2 offset = float2(sampleIndex * 7.0 + time * 50.0, sampleIndex * 11.0 + time * 30.0);
    float2 coord = pixelCoord + offset;

    float3 noise;
    noise.x = frac(52.9829189 * frac(dot(coord, float2(0.06711056, 0.00583715))));
    noise.y = frac(52.9829189 * frac(dot(coord + 100.0, float2(0.06711056, 0.00583715))));
    noise.z = frac(52.9829189 * frac(dot(coord + 200.0, float2(0.06711056, 0.00583715))));

    return noise;
}

// R2 quasi-random sequence for better sample distribution
float2 getR2Sequence(int n)
{
    // Plastic constant derived from tribonacci
    const float g = 1.32471795724474602596;
    const float a1 = 1.0 / g;
    const float a2 = 1.0 / (g * g);

    return frac(float2(a1, a2) * float(n) + 0.5);
}

// Cosine-weighted hemisphere sampling
float3 sampleCosineHemisphere(float2 xi, float3 normal)
{
    // Build tangent space
    float3 tangent = normalize(cross(normal, abs(normal.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0)));
    float3 bitangent = cross(normal, tangent);

    // Cosine-weighted sampling
    float r = sqrt(xi.x);
    float theta = 2.0 * kPi * xi.y;

    float3 localDir;
    localDir.x = r * cos(theta);
    localDir.y = r * sin(theta);
    localDir.z = sqrt(max(0.0, 1.0 - xi.x));

    // Transform to world space
    return normalize(tangent * localDir.x + bitangent * localDir.y + normal * localDir.z);
}

// GGX importance sampling for specular GI
float3 sampleGGXHemisphere(float2 xi, float3 normal, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;

    float phi = 2.0 * kPi * xi.x;
    float cosTheta = sqrt((1.0 - xi.y) / (1.0 + (a2 - 1.0) * xi.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    // Microfacet normal in tangent space
    float3 H;
    H.x = sinTheta * cos(phi);
    H.y = sinTheta * sin(phi);
    H.z = cosTheta;

    // Build tangent space
    float3 up = abs(normal.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, normal));
    float3 bitangent = cross(normal, tangent);

    // Transform to world space
    return normalize(tangent * H.x + bitangent * H.y + normal * H.z);
}

// =============================================================================
// VISIBILITY AND OCCLUSION
// =============================================================================

// Calculate ambient occlusion factor using cone tracing approximation
float calculateGIAO(float3 position, float3 normal, float2 pixelCoord, float time, int samples)
{
#if !GI_AMBIENT_OCCLUSION
    return 1.0;
#endif

    float occlusion = 0.0;
    float maxOcclusionDist = GI_AO_RADIUS;

    [loop]
    for (int i = 0; i < samples; i++)
    {
        // Get sample direction in hemisphere
        float2 xi = getGIBlueNoise(pixelCoord, time, i).xy;
        xi = frac(xi + getR2Sequence(i));

        float3 sampleDir = sampleCosineHemisphere(xi, normal);

        // Trace short ray for AO
        RayQuery<RAY_FLAG_NONE> aoQuery;
        RayDesc aoRay;
        aoRay.Origin = position + normal * 0.01;
        aoRay.Direction = sampleDir;
        aoRay.TMin = 0.0;
        aoRay.TMax = maxOcclusionDist;

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
            // Distance-weighted occlusion
            float occlusionFactor = 1.0 - saturate(hitDist / maxOcclusionDist);
            occlusion += occlusionFactor * occlusionFactor;
        }
    }

    occlusion /= float(samples);
    return saturate(1.0 - occlusion * GI_AO_STRENGTH);
}

// =============================================================================
// LIGHT SAMPLING
// =============================================================================

// Sample sunlight GI contribution
float3 calculateSunlightGI(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity)
{
#if !GI_SUNLIGHT
    return 0.0;
#endif

    // Direct sun contribution with hemisphere visibility
    float NdotL = saturate(dot(normal, sunDir));

    if (NdotL <= 0.0)
        return 0.0;

    // Simple shadow check
    RayQuery<RAY_FLAG_NONE> shadowQuery;
    RayDesc shadowRay;
    shadowRay.Origin = position + normal * 0.01;
    shadowRay.Direction = sunDir;
    shadowRay.TMin = 0.0;
    shadowRay.TMax = GI_MAX_DISTANCE;

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
        // In shadow - return minimal ambient
        return albedo * sunColor * 0.05 * GI_SUNLIGHT_STRENGTH;
    }

    // Lit by sun
    return albedo * sunColor * sunIntensity * NdotL * GI_SUNLIGHT_STRENGTH;
}

// Sample skylight GI contribution (ambient sky hemisphere)
float3 calculateSkylightGI(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 skyColor,
    float2 pixelCoord,
    float time)
{
#if !GI_SKYLIGHT
    return 0.0;
#endif

    float3 skylightAccum = 0.0;
    int samples = GI_SKYLIGHT_SAMPLES;

    [loop]
    for (int i = 0; i < samples; i++)
    {
        // Sample sky hemisphere
        float2 xi = getGIBlueNoise(pixelCoord, time, i + 100).xy;
        xi = frac(xi + getR2Sequence(i));

        float3 sampleDir = sampleCosineHemisphere(xi, normal);

        // Check if this direction can see the sky
        RayQuery<RAY_FLAG_NONE> skyQuery;
        RayDesc skyRay;
        skyRay.Origin = position + normal * 0.01;
        skyRay.Direction = sampleDir;
        skyRay.TMin = 0.0;
        skyRay.TMax = GI_MAX_DISTANCE;

        skyQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, skyRay);

        while (skyQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(skyQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                skyQuery.CommitNonOpaqueTriangleHit();
            }
        }

        if (skyQuery.CommittedStatus() != COMMITTED_TRIANGLE_HIT)
        {
            // Sees sky - add sky contribution
            // Weight by hemisphere direction (more from above)
            float skyWeight = saturate(sampleDir.y * 0.5 + 0.5);
            skylightAccum += skyColor * skyWeight;
        }
    }

    skylightAccum /= float(samples);
    return albedo * skylightAccum * GI_SKYLIGHT_STRENGTH;
}

// Sample blocklight (emissive) GI contribution
float3 calculateBlocklightGI(
    float3 position,
    float3 normal,
    float3 albedo,
    float2 pixelCoord,
    float time,
    int samples)
{
#if !GI_BLOCKLIGHT
    return 0.0;
#endif

    float3 blocklightAccum = 0.0;
    float totalWeight = 0.0;

    [loop]
    for (int i = 0; i < samples; i++)
    {
        // Generate random direction in hemisphere
        float3 noise = getGIBlueNoise(pixelCoord, time, i);
        float2 xi = frac(noise.xy + getR2Sequence(i));

        float3 sampleDir = sampleCosineHemisphere(xi, normal);

        // Trace ray looking for emissive surfaces
        RayQuery<RAY_FLAG_NONE> emissiveQuery;
        RayDesc emissiveRay;
        emissiveRay.Origin = position + normal * 0.01;
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

            // Check if hit surface is emissive
            if (surfaceInfo.emissive > 0.01)
            {
                float distance = hitInfo.rayT;

                // Physically-based inverse square falloff with range limit
                float falloff = 1.0 / (1.0 + distance * distance * GI_BLOCKLIGHT_FALLOFF);

                // Emissive surface contributes light
                float3 emittedLight = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;

                // Weight by angle (cosine already from hemisphere sampling)
                float weight = falloff;
                blocklightAccum += emittedLight * weight;
                totalWeight += weight;
            }
            else
            {
                // Non-emissive hit - could do indirect bounce here
                // For now, add tiny ambient contribution
                float distance = hitInfo.rayT;
                float falloff = 1.0 / (1.0 + distance * distance * 0.1);
                blocklightAccum += surfaceInfo.color * 0.01 * falloff;
                totalWeight += 0.1;
            }
        }
    }

    if (totalWeight > 0.0)
    {
        blocklightAccum /= totalWeight;
    }

    return albedo * blocklightAccum * GI_BLOCKLIGHT_STRENGTH;
}

// =============================================================================
// DIFFUSE GI (Multi-bounce)
// =============================================================================

// Single bounce diffuse GI
float3 traceDiffuseBounce(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float2 pixelCoord,
    float time,
    int bounceIndex)
{
    // Generate sample direction
    float3 noise = getGIBlueNoise(pixelCoord, time, bounceIndex);
    float2 xi = frac(noise.xy + getR2Sequence(bounceIndex));

    float3 sampleDir = sampleCosineHemisphere(xi, normal);

    // Trace ray
    RayQuery<RAY_FLAG_NONE> bounceQuery;
    RayDesc bounceRay;
    bounceRay.Origin = position + normal * 0.01;
    bounceRay.Direction = sampleDir;
    bounceRay.TMin = 0.0;
    bounceRay.TMax = GI_MAX_DISTANCE;

    bounceQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, bounceRay);

    while (bounceQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(bounceQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            bounceQuery.CommitNonOpaqueTriangleHit();
        }
    }

    float3 bounceLight = 0.0;

    if (bounceQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        HitInfo hitInfo = GetCommittedHitInfo(bounceQuery);
        ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
        GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
        SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

        float distance = hitInfo.rayT;

        // Get hit surface lighting
        float3 hitNormal = surfaceInfo.normal;
        float NdotL = saturate(dot(hitNormal, sunDir));

        // Direct light at hit point
        float3 hitDirect = surfaceInfo.color * sunColor * sunIntensity * NdotL;

        // Emissive contribution (with indirect boost)
        float3 hitEmissive = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY * INDIRECT_EMISSIVE_BOOST;

        // Ambient at hit point
        float3 hitAmbient = surfaceInfo.color * skyColor * 0.2;

        // Combine hit surface lighting
        bounceLight = hitDirect + hitEmissive + hitAmbient;

        // Distance falloff
        float falloff = 1.0 / (1.0 + distance * distance * GI_DISTANCE_FALLOFF);
        bounceLight *= falloff;
    }
    else
    {
        // Hit sky - sample sky color in that direction
        float skyBrightness = saturate(sampleDir.y * 0.5 + 0.5);
        bounceLight = skyColor * skyBrightness * 0.5;
    }

    // Apply surface albedo (diffuse transfer)
    return albedo * bounceLight;
}

// Full diffuse GI with multiple bounces
float3 calculateDiffuseGI(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float2 pixelCoord,
    float time)
{
#if !GI_DIFFUSE_BOUNCES
    return 0.0;
#endif

    float3 diffuseGI = 0.0;
    int samples = GI_DIFFUSE_SAMPLES;

    // First bounce
    [loop]
    for (int i = 0; i < samples; i++)
    {
        diffuseGI += traceDiffuseBounce(
            position, normal, albedo,
            sunDir, sunColor, sunIntensity, skyColor,
            pixelCoord, time, i);
    }

    diffuseGI /= float(samples);
    return diffuseGI * GI_DIFFUSE_STRENGTH;
}

// =============================================================================
// SPECULAR GI (Glossy Reflections)
// =============================================================================

float3 calculateSpecularGI(
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
    float2 pixelCoord,
    float time)
{
#if !GI_SPECULAR
    return 0.0;
#endif

    // Skip very rough surfaces
    if (roughness > GI_SPECULAR_MAX_ROUGHNESS)
        return 0.0;

    float3 specularGI = 0.0;
    int samples = GI_SPECULAR_SAMPLES;

    // Calculate F0 for Fresnel
    float3 f0 = lerp(float3(0.04, 0.04, 0.04), albedo, metalness);

    [loop]
    for (int i = 0; i < samples; i++)
    {
        // Generate sample direction using GGX importance sampling
        float3 noise = getGIBlueNoise(pixelCoord, time, i + 200);
        float2 xi = frac(noise.xy + getR2Sequence(i));

        float3 H = sampleGGXHemisphere(xi, normal, max(roughness, 0.04));
        float3 reflectDir = reflect(-viewDir, H);

        // Skip if reflection direction goes into surface
        if (dot(reflectDir, normal) <= 0.0)
            continue;

        // Trace reflection ray
        RayQuery<RAY_FLAG_NONE> specQuery;
        RayDesc specRay;
        specRay.Origin = position + normal * 0.01;
        specRay.Direction = reflectDir;
        specRay.TMin = 0.0;
        specRay.TMax = GI_MAX_DISTANCE;

        specQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, specRay);

        while (specQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(specQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                specQuery.CommitNonOpaqueTriangleHit();
            }
        }

        float3 sampleLight = 0.0;

        if (specQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(specQuery);
            ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
            GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
            SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

            // Shade reflected surface
            float3 hitNormal = surfaceInfo.normal;
            float NdotL = saturate(dot(hitNormal, sunDir));

            sampleLight = surfaceInfo.color * sunColor * sunIntensity * NdotL;
            sampleLight += surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;
            sampleLight += surfaceInfo.color * skyColor * 0.1;
        }
        else
        {
            // Hit sky
            float skyBrightness = saturate(reflectDir.y * 0.5 + 0.5);
            sampleLight = skyColor * skyBrightness;
        }

        // Apply Fresnel
        float NdotV = saturate(dot(normal, viewDir));
        float VdotH = saturate(dot(viewDir, H));
        float3 fresnel = f0 + (1.0 - f0) * pow(1.0 - VdotH, 5.0);

        specularGI += sampleLight * fresnel;
    }

    specularGI /= float(samples);

    // Roughness-based fade
    float roughnessFade = 1.0 - saturate(roughness / GI_SPECULAR_MAX_ROUGHNESS);

    return specularGI * GI_SPECULAR_STRENGTH * roughnessFade;
}

// =============================================================================
// MAIN GI CALCULATION
// =============================================================================

// Calculate complete global illumination for a surface
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
    // Fallback to simple emissive GI only
    return result;
#endif

    // Adjust sample counts based on quality
    int aoSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_AO_SAMPLES);
    int diffuseSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_DIFFUSE_SAMPLES);
    int skylightSamples = quality == 0 ? 1 : (quality == 1 ? 2 : GI_SKYLIGHT_SAMPLES);
    int blocklightSamples = quality == 0 ? 2 : (quality == 1 ? 4 : GI_BLOCKLIGHT_SAMPLES);

    // Ambient Occlusion
#if GI_AMBIENT_OCCLUSION
    result.ao = calculateGIAO(position, normal, pixelCoord, time, aoSamples);
#endif

    // Sunlight GI (direct + shadowed)
#if GI_SUNLIGHT
    float3 sunGI = calculateSunlightGI(position, normal, material.albedo, sunDir, sunColor, sunIntensity);
    result.diffuse += sunGI;
#endif

    // Skylight GI (ambient hemisphere)
#if GI_SKYLIGHT
    float3 skyGI = calculateSkylightGI(position, normal, material.albedo, skyColor, pixelCoord, time);
    result.diffuse += skyGI;
#endif

    // Blocklight GI (emissive sources)
#if GI_BLOCKLIGHT
    float3 blockGI = calculateBlocklightGI(position, normal, material.albedo, pixelCoord, time, blocklightSamples);
    result.emission += blockGI;
#endif

    // Diffuse bounce GI
#if GI_DIFFUSE_BOUNCES
    float3 diffuseBounce = calculateDiffuseGI(
        position, normal, material.albedo,
        sunDir, sunColor, sunIntensity, skyColor,
        pixelCoord, time);
    result.diffuse += diffuseBounce;
#endif

    // Specular GI (glossy reflections)
#if GI_SPECULAR
    if (material.roughness < GI_SPECULAR_MAX_ROUGHNESS)
    {
        result.specular = calculateSpecularGI(
            position, normal, viewDir,
            material.albedo, material.roughness, material.metalness,
            sunDir, sunColor, sunIntensity, skyColor,
            pixelCoord, time);
    }
#endif

    // Apply AO to indirect lighting
    result.diffuse *= result.ao;
    result.emission *= lerp(1.0, result.ao, 0.5);  // Partial AO on emission

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
