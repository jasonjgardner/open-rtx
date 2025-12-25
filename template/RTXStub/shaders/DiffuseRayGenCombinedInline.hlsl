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
// Diffuse Ray Generation for Global Illumination
// Generates and traces diffuse bounce rays for indirect lighting
// =============================================================================

#include "Include/Renderer.hlsl"
#include "Include/Util.hlsl"
#include "Include/OpenRTX.hlsl"
#include "Include/GI.hlsl"

// =============================================================================
// DIFFUSE RAY GENERATION STRUCTURES
// =============================================================================

// G-Buffer data for current pixel
struct GBufferData
{
    float3 position;        // World position
    float3 normal;          // Surface normal
    float3 albedo;          // Surface albedo
    float roughness;        // Surface roughness
    float metalness;        // Surface metalness
    float depth;            // Linear depth
    bool valid;             // Is this a valid surface?
};

// Diffuse ray hit result
struct DiffuseHitResult
{
    float3 radiance;        // Incoming radiance from hit
    float distance;         // Hit distance
    float3 hitNormal;       // Normal at hit point
    float3 hitAlbedo;       // Albedo at hit point
    bool hit;               // Did we hit anything?
};

// =============================================================================
// NOISE AND SAMPLING UTILITIES
// =============================================================================

// High-quality temporal blue noise
float2 getTemporalBlueNoise(uint2 pixelCoord, float time, uint sampleIndex)
{
    // Base interleaved gradient noise
    float2 baseNoise;
    float2 coord = float2(pixelCoord) + float2(sampleIndex * 17, sampleIndex * 31);
    baseNoise.x = frac(52.9829189 * frac(dot(coord, float2(0.06711056, 0.00583715))));
    baseNoise.y = frac(52.9829189 * frac(dot(coord + 100.0, float2(0.06711056, 0.00583715))));

    // Temporal offset using golden ratio
    float temporalOffset = frac(time * 7.0 + float(sampleIndex) * 0.618034);
    baseNoise = frac(baseNoise + temporalOffset);

    // R2 sequence for better sample distribution
    float2 r2 = getR2Sequence(sampleIndex);
    return frac(baseNoise + r2);
}

// Stratified sampling with jitter
float2 getStratifiedSample(uint sampleIndex, uint totalSamples, float2 jitter)
{
    uint gridSize = uint(sqrt(float(totalSamples)));
    uint x = sampleIndex % gridSize;
    uint y = sampleIndex / gridSize;

    float2 cell = float2(x, y) / float(gridSize);
    float2 cellJitter = jitter / float(gridSize);

    return cell + cellJitter;
}

// =============================================================================
// G-BUFFER READING (would read from primary pass output)
// =============================================================================

// Reconstruct world position from depth
float3 reconstructWorldPosition(uint2 pixelCoord, float depth, float4x4 invViewProj)
{
    float2 ndc = (float2(pixelCoord) + 0.5) * g_view.recipRenderResolution * 2.0 - 1.0;
    ndc.y = -ndc.y;

    float4 clipPos = float4(ndc, depth, 1.0);
    float4 worldPos = mul(invViewProj, clipPos);
    return worldPos.xyz / worldPos.w;
}

// =============================================================================
// DIFFUSE RAY TRACING
// =============================================================================

// Trace a single diffuse bounce ray
DiffuseHitResult traceDiffuseRay(
    float3 origin,
    float3 direction,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float3 ambientColor)
{
    DiffuseHitResult result;
    result.radiance = 0.0;
    result.distance = 0.0;
    result.hitNormal = float3(0, 1, 0);
    result.hitAlbedo = 0.0;
    result.hit = false;

    // Trace ray
    RayQuery<RAY_FLAG_NONE> diffuseQuery;
    RayDesc diffuseRay;
    diffuseRay.Origin = origin;
    diffuseRay.Direction = direction;
    diffuseRay.TMin = 0.001;
    diffuseRay.TMax = GI_MAX_DISTANCE;

    diffuseQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, diffuseRay);

    while (diffuseQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(diffuseQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            diffuseQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (diffuseQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        HitInfo hitInfo = GetCommittedHitInfo(diffuseQuery);
        ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
        GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
        SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

        result.hit = true;
        result.distance = hitInfo.rayT;
        result.hitNormal = surfaceInfo.normal;
        result.hitAlbedo = surfaceInfo.color;

        // Calculate lighting at hit point
        float NdotL = saturate(dot(surfaceInfo.normal, sunDir));

        // Simple shadow check for hit point
        float shadow = 1.0;
        if (NdotL > 0.0)
        {
            RayQuery<RAY_FLAG_NONE> shadowQuery;
            RayDesc shadowRay;
            shadowRay.Origin = surfaceInfo.position + surfaceInfo.normal * 0.01;
            shadowRay.Direction = sunDir;
            shadowRay.TMin = 0.0;
            shadowRay.TMax = GI_MAX_DISTANCE;

            shadowQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, shadowRay);

            while (shadowQuery.Proceed())
            {
                HitInfo shadowHit = GetCandidateHitInfo(shadowQuery);
                if (AlphaTestHitLogic(shadowHit))
                {
                    shadowQuery.CommitNonOpaqueTriangleHit();
                }
            }

            if (shadowQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
            {
                shadow = 0.0;
            }
        }

        // Direct lighting at hit point
        float3 directLight = surfaceInfo.color * sunColor * sunIntensity * NdotL * shadow;

        // Ambient lighting
        float3 ambient = surfaceInfo.color * ambientColor;

        // Emissive contribution (with indirect boost for GI)
        float3 emissive = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY * INDIRECT_EMISSIVE_BOOST;

        // Sky visibility approximation (hemisphere above normal)
        float skyVisibility = saturate(surfaceInfo.normal.y * 0.5 + 0.5);
        float3 skylight = surfaceInfo.color * skyColor * skyVisibility * 0.3;

        // Combine all lighting
        result.radiance = directLight + ambient + emissive + skylight;

        // Distance-based falloff to prevent distant surfaces from contributing too much
        float distanceFalloff = 1.0 / (1.0 + result.distance * result.distance * GI_DISTANCE_FALLOFF);
        result.radiance *= distanceFalloff;
    }
    else
    {
        // Hit sky - sample sky color in ray direction
        float skyBrightness = saturate(direction.y * 0.5 + 0.5);
        result.radiance = skyColor * skyBrightness;

        // Add sun contribution if looking toward sun
        float sunDot = saturate(dot(direction, sunDir));
        if (sunDot > 0.99)
        {
            result.radiance += sunColor * sunIntensity * pow(sunDot, 64.0);
        }
    }

    return result;
}

// =============================================================================
// MULTI-SAMPLE DIFFUSE GI
// =============================================================================

// Calculate diffuse GI with multiple samples and temporal accumulation
float3 calculateDiffuseGIMultiSample(
    float3 position,
    float3 normal,
    float3 albedo,
    float3 sunDir,
    float3 sunColor,
    float sunIntensity,
    float3 skyColor,
    float3 ambientColor,
    uint2 pixelCoord,
    float time,
    uint numSamples)
{
    float3 totalRadiance = 0.0;
    float totalWeight = 0.0;

    // Build tangent space
    float3 tangent = normalize(cross(normal, abs(normal.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0)));
    float3 bitangent = cross(normal, tangent);

    [loop]
    for (uint i = 0; i < numSamples; i++)
    {
        // Get sample direction using cosine-weighted hemisphere
        float2 xi = getTemporalBlueNoise(pixelCoord, time, i);

        // Cosine-weighted sampling
        float r = sqrt(xi.x);
        float theta = 2.0 * kPi * xi.y;

        float3 localDir;
        localDir.x = r * cos(theta);
        localDir.y = r * sin(theta);
        localDir.z = sqrt(max(0.0, 1.0 - xi.x));

        // Transform to world space
        float3 sampleDir = normalize(tangent * localDir.x + bitangent * localDir.y + normal * localDir.z);

        // Trace ray
        DiffuseHitResult hitResult = traceDiffuseRay(
            position + normal * 0.01,
            sampleDir,
            sunDir, sunColor, sunIntensity,
            skyColor, ambientColor);

        // PDF for cosine-weighted sampling is cos(theta)/pi
        // But since we're importance sampling, the weight is 1
        float weight = 1.0;

        totalRadiance += hitResult.radiance * weight;
        totalWeight += weight;
    }

    if (totalWeight > 0.0)
    {
        totalRadiance /= totalWeight;
    }

    // Apply surface albedo (diffuse reflection)
    return albedo * totalRadiance;
}

// =============================================================================
// SPHERICAL HARMONICS ENCODING (for efficient storage and filtering)
// =============================================================================

// L1 spherical harmonics bands (4 coefficients)
struct SH_L1
{
    float4 coeffs[3];  // RGB for each SH band
};

// Encode radiance into L1 spherical harmonics
SH_L1 encodeSH_L1(float3 direction, float3 radiance)
{
    SH_L1 sh;

    // L0 band (constant)
    float Y00 = 0.282095;  // 1/(2*sqrt(pi))

    // L1 bands (linear)
    float Y1n1 = 0.488603 * direction.y;  // sqrt(3/(4*pi))
    float Y10 = 0.488603 * direction.z;
    float Y11 = 0.488603 * direction.x;

    // Encode each color channel
    sh.coeffs[0] = float4(radiance.r * Y00, radiance.r * Y1n1, radiance.r * Y10, radiance.r * Y11);
    sh.coeffs[1] = float4(radiance.g * Y00, radiance.g * Y1n1, radiance.g * Y10, radiance.g * Y11);
    sh.coeffs[2] = float4(radiance.b * Y00, radiance.b * Y1n1, radiance.b * Y10, radiance.b * Y11);

    return sh;
}

// Decode radiance from L1 spherical harmonics for given direction
float3 decodeSH_L1(SH_L1 sh, float3 direction)
{
    float Y00 = 0.282095;
    float Y1n1 = 0.488603 * direction.y;
    float Y10 = 0.488603 * direction.z;
    float Y11 = 0.488603 * direction.x;

    float3 radiance;
    radiance.r = sh.coeffs[0].x * Y00 + sh.coeffs[0].y * Y1n1 + sh.coeffs[0].z * Y10 + sh.coeffs[0].w * Y11;
    radiance.g = sh.coeffs[1].x * Y00 + sh.coeffs[1].y * Y1n1 + sh.coeffs[1].z * Y10 + sh.coeffs[1].w * Y11;
    radiance.b = sh.coeffs[2].x * Y00 + sh.coeffs[2].y * Y1n1 + sh.coeffs[2].z * Y10 + sh.coeffs[2].w * Y11;

    return max(radiance, 0.0);
}

// Accumulate SH sample
void accumulateSH_L1(inout SH_L1 sh, float3 direction, float3 radiance)
{
    SH_L1 sample = encodeSH_L1(direction, radiance);
    sh.coeffs[0] += sample.coeffs[0];
    sh.coeffs[1] += sample.coeffs[1];
    sh.coeffs[2] += sample.coeffs[2];
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

    // Read G-Buffer data (would be from primary pass outputs)
    // For now, we reconstruct from available data
    float2 uv = (float2(pixelCoord) + 0.5) * g_view.recipRenderResolution;

    // Trace a primary ray to get surface data
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

    // Skip emissive surfaces (they don't need GI, they provide it)
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

    // Calculate number of samples based on quality and distance
    float distanceFromCamera = hitInfo.rayT;
    float distanceLOD = saturate(distanceFromCamera / 128.0);

    // Reduce samples at distance
    uint numSamples = uint(lerp(float(GI_DIFFUSE_SAMPLES), 1.0, distanceLOD * distanceLOD));
    numSamples = max(1, numSamples);

    // Calculate diffuse GI
    float3 diffuseGI = calculateDiffuseGIMultiSample(
        surfaceInfo.position,
        surfaceInfo.normal,
        surfaceInfo.color,
        sunDir, sunColor, sunIntensity,
        skyColor, ambientColor,
        pixelCoord,
        g_view.time,
        numSamples);

    // Apply GI strength
    diffuseGI *= GI_DIFFUSE_STRENGTH;

    // Clamp to prevent fireflies
    diffuseGI = min(diffuseGI, 10.0);

    // Output would go to diffuse GI buffer for denoising
    // denoisingInputs[bufferIndex][pixelCoord] = float4(diffuseGI, 1.0);
}
