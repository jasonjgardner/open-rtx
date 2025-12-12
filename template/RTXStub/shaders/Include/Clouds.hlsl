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
// OpenRTX Volumetric Clouds
// Advanced volumetric cloud rendering with multiple noise types and lighting
// =============================================================================

#ifndef __OPENRTX_CLOUDS_HLSL__
#define __OPENRTX_CLOUDS_HLSL__

#include "Settings.hlsl"

// =============================================================================
// CLOUD CONTEXT STRUCTURE
// =============================================================================

struct CloudContext
{
    float3 cameraPosition;
    float3 sunDirection;
    float3 lightColor;
    float3 ambientColor;
    float time;
    float rainStrength;
    float shadowFade;      // Day/night transition (1 = day, 0 = night)
    float sunVisibility;   // Sun visibility (1 = visible, 0 = below horizon)
    float moonVisibility;  // Moon visibility (1 = visible, 0 = below horizon)
    float fogDensity;
};

// =============================================================================
// NOISE FUNCTIONS
// =============================================================================

// Simple 2D hash function
float cloudHash2D(float2 p)
{
    return frac(sin(dot(p, float2(12.9898, 4.1414))) * 43758.5453);
}

// Simple 3D hash function
float cloudHash3D(float3 p)
{
    p = frac(p * 0.3183099 + 0.1);
    p *= 17.0;
    return frac(p.x * p.y * p.z * (p.x + p.y + p.z));
}

// Smooth noise 2D (replaces texture lookup)
float cloudNoise2D(float2 p)
{
    float2 i = floor(p);
    float2 f = frac(p);
    f = f * f * (3.0 - 2.0 * f);

    float a = cloudHash2D(i);
    float b = cloudHash2D(i + float2(1.0, 0.0));
    float c = cloudHash2D(i + float2(0.0, 1.0));
    float d = cloudHash2D(i + float2(1.0, 1.0));

    return lerp(lerp(a, b, f.x), lerp(c, d, f.x), f.y);
}

// Smooth noise 3D
float cloudNoise3D(float3 p)
{
    float3 i = floor(p);
    float3 f = frac(p);
    f = f * f * (3.0 - 2.0 * f);

    float n000 = cloudHash3D(i);
    float n100 = cloudHash3D(i + float3(1, 0, 0));
    float n010 = cloudHash3D(i + float3(0, 1, 0));
    float n110 = cloudHash3D(i + float3(1, 1, 0));
    float n001 = cloudHash3D(i + float3(0, 0, 1));
    float n101 = cloudHash3D(i + float3(1, 0, 1));
    float n011 = cloudHash3D(i + float3(0, 1, 1));
    float n111 = cloudHash3D(i + float3(1, 1, 1));

    float4 low = float4(n000, n100, n010, n110);
    float4 high = float4(n001, n101, n011, n111);
    float4 mix1 = lerp(low, high, f.z);
    float2 mix2 = lerp(mix1.xy, mix1.zw, f.y);
    return lerp(mix2.x, mix2.y, f.x);
}

// FBM noise for clouds
float cloudFBM(float2 p, int octaves)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    [loop]
    for (int i = 0; i < octaves; i++)
    {
        value += amplitude * cloudNoise2D(p * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value;
}

// Worley noise (cellular)
float cloudWorley(float2 p)
{
    float2 i = floor(p);
    float2 f = frac(p);

    float minDist = 1.0;

    [unroll]
    for (int y = -1; y <= 1; y++)
    {
        [unroll]
        for (int x = -1; x <= 1; x++)
        {
            float2 neighbor = float2(x, y);
            float2 cellPoint = cloudHash2D(i + neighbor).xx * 0.5 + 0.5;
            cellPoint = 0.5 + 0.5 * sin(6.2831 * cellPoint);
            float2 diff = neighbor + cellPoint - f;
            float dist = length(diff);
            minDist = min(minDist, dist);
        }
    }

    return minDist;
}

// =============================================================================
// CLOUD BASE SAMPLING FUNCTIONS
// =============================================================================

// Perlin-based cloud sampling
float cloudSampleBasePerlin(float2 coord)
{
    float noise = cloudFBM(coord * 4.0, 4);
    return noise;
}

// Worley-based cloud sampling
float cloudSampleBaseWorley(float2 coord)
{
    float noise = cloudWorley(coord * 8.0);
    noise = pow(1.0 - noise, 2.0) * 0.5 + 0.25;
    return noise;
}

// Blocky (Minecraft-style) cloud sampling
float cloudSampleBaseBlocky(float2 coord, float rainStrength)
{
    float noiseRes = 512.0;

    coord = coord * noiseRes - 0.5;

    float2 flr = floor(coord);
    float2 frc = coord - flr;

    // Blocky interpolation
    frc = saturate(frc * 5.0 - 2.0);
    frc = frc * frc * (3.0 - 2.0 * frc);

    coord = (flr + frc + 0.5) / noiseRes;

    float noiseBase = cloudNoise2D(coord * 32.0);
    noiseBase = (1.0 - noiseBase) * 4.0;

    // Rain variation
    float noiseRain = cloudNoise2D((coord + float2(0.5, 0.0)) * 32.0);
    noiseRain = (1.0 - noiseRain) * 4.0 * smoothstep(0.0, 0.5, rainStrength);

    noiseBase = min(noiseBase + noiseRain, 1.0);

    return noiseBase;
}

// Detail noise sampling
float cloudSampleDetail(float2 coord, float cloudGradient)
{
    float detailZ = floor(cloudGradient * float(CLOUD_THICKNESS)) * 0.04;
    float detailFrac = frac(cloudGradient * float(CLOUD_THICKNESS));

    float noiseDetailLow = cloudNoise2D(coord + detailZ);
    float noiseDetailHigh = cloudNoise2D(coord + detailZ + 0.04);

    return lerp(noiseDetailLow, noiseDetailHigh, detailFrac);
}

// =============================================================================
// CLOUD COVERAGE FUNCTIONS
// =============================================================================

// Default coverage (for Perlin and Worley)
float cloudCoverageDefault(float cloudGradient)
{
    float coverage = abs(cloudGradient - 0.125);

    coverage *= cloudGradient > 0.125 ? (2.14 - CLOUD_AMOUNT * 0.1) : 8.0;
    coverage = coverage * coverage * 4.0;

    return coverage;
}

// Blocky coverage
float cloudCoverageBlocky(float cloudGradient)
{
    float coverage = abs(cloudGradient - 0.5) * 2.0;

    coverage *= coverage;
    coverage *= coverage;

    return coverage;
}

// Apply density function
float cloudApplyDensity(float noise, float rainStrength)
{
    noise *= CLOUD_DENSITY * 0.125;
    noise *= (1.0 - 0.75 * rainStrength);
    noise = noise / sqrt(noise * noise + 0.5);

    return noise;
}

// =============================================================================
// CLOUD COMBINATION FUNCTIONS
// =============================================================================

// Combine for default (Perlin/Worley) clouds
float cloudCombineDefault(float noiseBase, float noiseDetail, float noiseCoverage,
                          float sunCoverage, float rainStrength)
{
    float noise = lerp(noiseBase, noiseDetail, 0.0476 * CLOUD_DETAIL) * 21.0;

    noise = lerp(noise - noiseCoverage, 21.0 - noiseCoverage * 2.5, 0.33 * rainStrength);
    noise = max(noise - (sunCoverage * 3.0 + CLOUD_AMOUNT), 0.0);

    noise = cloudApplyDensity(noise, rainStrength);

    return noise;
}

// Combine for blocky clouds
float cloudCombineBlocky(float noiseBase, float noiseCoverage, float rainStrength)
{
    float noise = (noiseBase - noiseCoverage) * 2.0;
    noise = max(noise, 0.0);
    noise = cloudApplyDensity(noise, rainStrength);

    return noise;
}

// =============================================================================
// MAIN CLOUD SAMPLING FUNCTION
// =============================================================================

float cloudSample(float2 coord, float2 wind, float cloudGradient, float sunCoverage,
                  float rainStrength, float dither)
{
    coord *= 0.004 * CLOUD_STRETCH;

#if CLOUD_BASE_TYPE == 0
    // Perlin-based clouds
    float2 baseCoord = coord * 0.25 + wind;
    float2 detailCoord = coord * 0.5 - wind * 2.0;

    float noiseBase = cloudSampleBasePerlin(baseCoord);
    float noiseDetail = cloudSampleDetail(detailCoord, cloudGradient);
    float noiseCoverage = cloudCoverageDefault(cloudGradient);

    return cloudCombineDefault(noiseBase, noiseDetail, noiseCoverage, sunCoverage, rainStrength);

#elif CLOUD_BASE_TYPE == 1
    // Worley-based clouds
    float2 baseCoord = coord * 0.5 + wind * 2.0;
    float2 detailCoord = coord * 0.5 - wind * 2.0;

    float noiseBase = cloudSampleBaseWorley(baseCoord);
    float noiseDetail = cloudSampleDetail(detailCoord, cloudGradient);
    float noiseCoverage = cloudCoverageDefault(cloudGradient);

    return cloudCombineDefault(noiseBase, noiseDetail, noiseCoverage, sunCoverage, rainStrength);

#else
    // Blocky (Minecraft-style) clouds
    float2 baseCoord = coord * 0.125 + wind * 0.5;

    float noiseBase = cloudSampleBaseBlocky(baseCoord, rainStrength);
    float noiseCoverage = cloudCoverageBlocky(cloudGradient);

    return cloudCombineBlocky(noiseBase, noiseCoverage, rainStrength);
#endif
}

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

float invLerp(float v, float l, float h)
{
    // Guard against division by zero
    float range = h - l;
    return saturate((v - l) / max(range, 0.0001));
}

// Get cloud mask based on fog
float getCloudMask(float3 viewPos, float z, float fogDensity)
{
    if (z == 1.0) return 1.0;

    float fogFactor = length(viewPos);
    float fogFar = 256.0;  // Render distance approximation
    float vanillaDensity = 0.2 * max(fogDensity, 0.5);

    float vanillaFog = 1.0 - (fogFar - fogFactor) / (vanillaDensity * fogFar);
    vanillaFog = saturate(vanillaFog);

    return vanillaFog;
}

// =============================================================================
// CLOUD OUTPUT STRUCTURE
// =============================================================================

struct CloudOutput
{
    float3 color;
    float transmittance;
    float depth;
};

// =============================================================================
// SKYBOX CLOUD RENDERING
// =============================================================================

CloudOutput drawCloudSkybox(float3 viewPos, float z, float dither, CloudContext ctx, bool fadeFaster)
{
    CloudOutput output;
    output.color = float3(0.0, 0.0, 0.0);
    output.transmittance = 1.0;
    output.depth = 10000.0;

#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    float cloudMask = getCloudMask(viewPos, z, ctx.fogDensity);
    if (cloudMask == 0.0) return output;

    // Temporal dithering
    dither = frac(dither + ctx.time * 0.618);

    int samples = CLOUD_THICKNESS * 2;

    float cloud = 0.0;
    float cloudLighting = 0.0;

    float sampleStep = 1.0 / float(samples);
    float currentStep = dither * sampleStep;

    float3 nViewPos = normalize(viewPos);
    float VoU = dot(nViewPos, float3(0, 1, 0));  // Up vector
    float VoL = dot(nViewPos, ctx.sunDirection);

    // Sun coverage for cloud reveal effect
#if ENABLE_CLOUD_REVEAL
    float sunCoverage = lerp(abs(VoL), max(VoL, 0.0), ctx.shadowFade);
    sunCoverage = pow(saturate(sunCoverage * 2.0 - 1.0), 12.0) * (1.0 - ctx.rainStrength);
#else
    float sunCoverage = 0.0;
#endif

    float convertedHeight = (max(CLOUD_HEIGHT, 128.0) - 72.0) / CLOUD_SCALE;

    float2 wind = float2(
        ctx.time * CLOUD_SPEED * 0.0005,
        sin(ctx.time * CLOUD_SPEED * 0.001) * 0.005
    ) * 0.667;

    float3 cloudColor = float3(0.0, 0.0, 0.0);

    if (VoU > 0.025)
    {
        float3 wpos = normalize(viewPos);
        // Guard against division by near-zero y component
        float safeWposY = max(abs(wpos.y), 0.001) * sign(wpos.y + 0.0001);

        float halfVoL = lerp(abs(VoL) * 0.8, VoL, ctx.shadowFade) * 0.5 + 0.5;
        float halfVoLSqr = halfVoL * halfVoL;
        float scattering = pow(halfVoL, 6.0);
        float noiseLightFactor = (2.0 - 1.5 * VoL * ctx.shadowFade) * CLOUD_DENSITY * 0.5;

        [loop]
        for (int i = 0; i < samples; i++)
        {
            if (cloud > 0.99) break;

#if CLOUD_BASE_TYPE == 2
            float planeY = convertedHeight + currentStep * float(CLOUD_THICKNESS) * 0.5;
#else
            float planeY = convertedHeight + currentStep * float(CLOUD_THICKNESS);
#endif

            float3 planeCoord = wpos * (planeY / safeWposY);
            float2 cloudCoord = ctx.cameraPosition.xz * 0.0625 * 16.0 / CLOUD_SCALE + planeCoord.xz;

            float noise = cloudSample(cloudCoord, wind, currentStep, sunCoverage, ctx.rainStrength, dither);

            // Use max() to avoid undefined behavior when base is 0 with non-integer exponent
            float sampleLighting = pow(max(currentStep, 0.0001), 1.125 * halfVoLSqr + 0.875) * 0.8 + 0.2;
            sampleLighting *= 1.0 - pow(max(noise, 0.0001), max(noiseLightFactor, 0.0001));

            cloudLighting = lerp(cloudLighting, sampleLighting, noise * (1.0 - cloud * cloud));
            cloud = lerp(cloud, 1.0, noise);

            currentStep += sampleStep;
        }

        cloudLighting = lerp(cloudLighting, 1.0, (1.0 - cloud * cloud) * scattering * 0.5);
        cloudLighting *= (1.0 - 0.9 * ctx.rainStrength);

        cloudColor = lerp(
            ctx.ambientColor * (0.3 * ctx.sunVisibility + 0.5),
            ctx.lightColor * (0.85 + 1.15 * scattering),
            cloudLighting
        );
        cloudColor *= 1.0 - 0.4 * ctx.rainStrength;

        cloud *= saturate(1.0 - exp(-16.0 / max(ctx.fogDensity, 0.5) * VoU + 0.5));

        if (fadeFaster)
        {
            cloud *= 1.0 - pow(1.0 - VoU, 4.0);
        }
    }

    cloudColor *= CLOUD_BRIGHTNESS * (0.5 - 0.25 * (1.0 - ctx.sunVisibility) * (1.0 - ctx.rainStrength));

    // Height-based visibility
    cloudColor *= saturate((ctx.cameraPosition.y + 70.0) / 8.0);

    cloud *= cloud;
    cloud *= cloudMask;
    cloud *= CLOUD_OPACITY;

    output.color = cloudColor;
    output.transmittance = 1.0 - cloud;

    return output;
}

// =============================================================================
// VOLUMETRIC CLOUD RENDERING
// =============================================================================

CloudOutput drawCloudVolumetric(float3 viewPos, float3 cameraPos, float z, float dither,
                                CloudContext ctx, inout float cloudViewLength, bool fadeFaster)
{
    CloudOutput output;
    output.color = float3(0.0, 0.0, 0.0);
    output.transmittance = 1.0;
    output.depth = 10000.0;

#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    float cloudMask = getCloudMask(viewPos, z, ctx.fogDensity);

    // Temporal dithering
    dither = frac(dither + ctx.time * 0.618);

    float3 nViewPos = normalize(viewPos);
    float3 worldPos = viewPos;  // Already in world space
    float3 nWorldPos = normalize(worldPos);
    float viewLength = length(viewPos);

    float cloudThickness = float(CLOUD_THICKNESS);

    int maxSamples = CLOUD_MARCH_STEPS;

#if CLOUD_BASE_TYPE == 2
    cloudThickness *= 0.5;
    maxSamples = min(maxSamples * 2, 64);
#endif

    float lowerY = CLOUD_HEIGHT;
    float upperY = lowerY + cloudThickness * CLOUD_SCALE;

    // Guard against division by near-zero y component
    float safeNWorldPosY = (abs(nWorldPos.y) < 0.001) ? 0.001 * sign(nWorldPos.y + 0.0001) : nWorldPos.y;
    float lowerPlane = (lowerY - cameraPos.y) / safeNWorldPosY;
    float upperPlane = (upperY - cameraPos.y) / safeNWorldPosY;

    float nearestPlane = max(min(lowerPlane, upperPlane), 0.0);
    float furthestPlane = max(lowerPlane, upperPlane);

    float maxCloudViewLength = cloudViewLength;

    if (furthestPlane < 0.0) return output;

    float planeDifference = furthestPlane - nearestPlane;

    float3 startPos = cameraPos + nearestPlane * nWorldPos;

    float heightRange = max((upperY - lowerY) * 0.5, 0.001);
    float lengthScaling = abs(cameraPos.y - (upperY + lowerY) * 0.5) / heightRange;
    lengthScaling = saturate((lengthScaling - 1.0) * cloudThickness * 0.125);

    float sampleLength = cloudThickness * CLOUD_SCALE / 2.0;
    float lengthDivisor = (4.0 * safeNWorldPosY * safeNWorldPosY) * lengthScaling + 1.0;
    sampleLength /= max(lengthDivisor, 0.001);

    float3 sampleStep = nWorldPos * sampleLength;
    float safeSampleLength = max(sampleLength, 0.001);
    int samples = int(min(planeDifference / safeSampleLength, float(maxSamples)) + 1.0);

    float3 samplePos = startPos + sampleStep * dither;
    float sampleTotalLength = nearestPlane + sampleLength * dither;

    float2 wind = float2(
        ctx.time * CLOUD_SPEED * 0.0005,
        sin(ctx.time * CLOUD_SPEED * 0.001) * 0.005
    ) * 0.667;

    float cloud = 0.0;
    float cloudFaded = 0.0;
    float cloudLighting = 0.0;

    float VoU = dot(nViewPos, float3(0, 1, 0));
    float VoL = dot(nViewPos, ctx.sunDirection);
    float VoS = dot(nViewPos, ctx.sunDirection);  // Same as VoL for sun

    // Sun coverage for cloud reveal
#if ENABLE_CLOUD_REVEAL
    float sunCoverage = lerp(abs(VoL), max(VoL, 0.0), ctx.shadowFade);
    sunCoverage = pow(saturate(sunCoverage * 2.0 - 1.0), 12.0) * (1.0 - ctx.rainStrength);
#else
    float sunCoverage = 0.0;
#endif

    float halfVoL = lerp(abs(VoL) * 0.8, VoL, ctx.shadowFade) * 0.5 + 0.5;
    float halfVoLSqr = halfVoL * halfVoL;
    float halfVoS = VoS * 0.5 + 0.5;
    float halfVoSSqr = halfVoS * halfVoS;

    float scattering = pow(halfVoL, 6.0);
    float noiseLightFactor = (2.0 - 1.5 * VoL * ctx.shadowFade) * CLOUD_DENSITY * 0.5;

    float viewLengthSoftMin = viewLength - sampleLength * 0.5;
    float viewLengthSoftMax = viewLength + sampleLength * 0.5;

    float distanceFade = 1.0;
    float fadeStart = 32.0 / max(ctx.fogDensity, 0.5);
    float fadeEnd = (fadeFaster ? 80.0 : 240.0) / max(ctx.fogDensity, 0.5);

    float xzNormalizeFactor = 10.0 / max(abs(CLOUD_HEIGHT - 72.0), 56.0);

    [loop]
    for (int i = 0; i < samples; i++)
    {
        if (cloud > 0.99) break;
        if (sampleTotalLength > viewLengthSoftMax && cloudMask == 0.0) break;

        float cloudGradient = invLerp(samplePos.y, lowerY, upperY);
        float xzNormalizedDistance = length(samplePos.xz - cameraPos.xz) * xzNormalizeFactor;
        float2 cloudCoord = samplePos.xz / CLOUD_SCALE;

        float noise = cloudSample(cloudCoord, wind, cloudGradient, sunCoverage, ctx.rainStrength, dither);
        noise *= step(lowerY, samplePos.y) * step(samplePos.y, upperY);

        // Use max() to avoid undefined behavior when base is 0 with non-integer exponent
        float sampleLighting = pow(max(cloudGradient, 0.0001), 1.125 * halfVoLSqr + 0.875) * 0.8 + 0.2;
        sampleLighting *= 1.0 - pow(max(noise, 0.0001), max(noiseLightFactor, 0.0001));

        float sampleFade = invLerp(xzNormalizedDistance, fadeEnd, fadeStart);
        distanceFade *= lerp(1.0, sampleFade, noise * (1.0 - cloud));
        noise *= step(xzNormalizedDistance, fadeEnd);

        cloudLighting = lerp(cloudLighting, sampleLighting, noise * (1.0 - cloud * cloud));

        cloud = lerp(cloud, 1.0, noise);

        if (z < 1.0)
        {
            float softFade = invLerp(sampleTotalLength, viewLengthSoftMax, viewLengthSoftMin);
            noise *= lerp(softFade, 1.0, cloudMask * cloudMask);
        }

        cloudFaded = lerp(cloudFaded, 1.0, noise);

        if (cloudViewLength == maxCloudViewLength && cloud > 0.5)
        {
            cloudViewLength = sampleTotalLength;
        }

        samplePos += sampleStep;
        sampleTotalLength += sampleLength;
    }

    cloudFaded *= distanceFade;

    cloudLighting = lerp(cloudLighting, 1.0, (1.0 - cloud * cloud) * scattering * 0.5);
    cloudLighting *= (1.0 - 0.9 * ctx.rainStrength);

    float horizonSunFactor = ctx.sunVisibility * (1.0 - ctx.sunVisibility) * halfVoSSqr * 0.5;
    // Use component-wise pow with max() to avoid undefined behavior with zero/negative values
    float horizonExp = 4.0 - 3.0 * ctx.sunVisibility;
    float3 horizonSunCol = float3(
        pow(max(ctx.lightColor.r, 0.0001), horizonExp),
        pow(max(ctx.lightColor.g, 0.0001), horizonExp),
        pow(max(ctx.lightColor.b, 0.0001), horizonExp)
    );

    float3 cloudAmbientCol = ctx.ambientColor * (0.3 * ctx.sunVisibility + 0.5);
    float3 cloudLightCol = lerp(ctx.lightColor, horizonSunCol, horizonSunFactor);
    cloudLightCol *= 0.85 + 1.15 * scattering;

    float3 cloudColor = lerp(cloudAmbientCol, cloudLightCol, cloudLighting);

    cloudColor *= 1.0 - 0.4 * ctx.rainStrength;
    cloudColor *= CLOUD_BRIGHTNESS * (0.5 - 0.25 * (1.0 - ctx.sunVisibility) * (1.0 - ctx.rainStrength));

    cloudFaded *= cloudFaded * CLOUD_OPACITY;

    if (cloudFaded < dither)
    {
        cloudViewLength = maxCloudViewLength;
    }

    // Height-based visibility
    cloudColor *= saturate((cameraPos.y + 70.0) / 8.0);

    output.color = cloudColor;
    output.transmittance = 1.0 - cloudFaded;
    output.depth = cloudViewLength;

    return output;
}

// =============================================================================
// STAR RENDERING
// =============================================================================

void drawStars(inout float3 color, float3 viewPos, CloudContext ctx)
{
#if !ENABLE_STARS
    return;
#endif

    float3 wpos = viewPos * 100.0;
    // Guard against division by zero
    float denom = wpos.y + length(wpos.xz);
    float3 planeCoord = wpos / max(denom, 0.001);
    float2 wind = float2(ctx.time, 0.0);
    float2 coord = planeCoord.xz * 0.4 + ctx.cameraPosition.xz * 0.0001 + wind * 0.00125;
    coord = floor(coord * 1024.0) / 1024.0;

    float3 nViewPos = normalize(viewPos);
    float VoU = saturate(dot(nViewPos, float3(0, 1, 0)));
    float VoL = dot(nViewPos, ctx.sunDirection);
    float multiplier = sqrt(sqrt(VoU)) * STAR_BRIGHTNESS * (1.0 - ctx.rainStrength) * ctx.moonVisibility;

    float star = 1.0;
    if (VoU > 0.0)
    {
        star *= cloudHash2D(coord);
        star *= cloudHash2D(coord + 0.10);
        star *= cloudHash2D(coord + 0.23);
    }

    // Adjust threshold based on density setting
    float threshold = 0.8125 + (1.0 - STAR_DENSITY) * 0.15;
    star = saturate(star - threshold) * multiplier;

    // Height-based visibility
    star *= saturate((ctx.cameraPosition.y + 70.0) / 8.0);

    // Fade near sun/moon
    float moonFade = smoothstep(-0.997, -0.992, VoL);
    star *= moonFade;

    // Night sky tint - use component-wise pow to avoid HLSL issues
    float3 nightColor = float3(
        pow(0.6, 0.8),
        pow(0.7, 0.8),
        pow(1.0, 0.8)
    );
    color += star * nightColor;
}

// =============================================================================
// AURORA RENDERING
// =============================================================================

float auroraSample(float2 coord, float2 wind, float VoU)
{
    float noise = cloudNoise2D(coord * 0.0625 + wind * 0.25) * 3.0;
    noise += cloudNoise2D(coord * 0.03125 + wind * 0.15) * 3.0;

    noise = max(1.0 - 4.0 * (0.5 * VoU + 0.5) * abs(noise - 3.0), 0.0);

    return noise;
}

float3 drawAurora(float3 viewPos, float dither, CloudContext ctx)
{
#if AURORA_MODE == 0
    return float3(0.0, 0.0, 0.0);
#endif

    // Temporal dithering
    dither = frac(dither + ctx.time * 0.618);

    float sampleStep = 1.0 / float(AURORA_SAMPLES);
    float currentStep = dither * sampleStep;

    float VoU = dot(normalize(viewPos), float3(0, 1, 0));

    float visibility = ctx.moonVisibility * (1.0 - ctx.rainStrength) * (1.0 - ctx.rainStrength);

    float2 wind = float2(
        ctx.time * CLOUD_SPEED * 0.000125,
        sin(ctx.time * CLOUD_SPEED * 0.05) * 0.00025
    );

    float3 aurora = float3(0.0, 0.0, 0.0);

    if (VoU > 0.0 && visibility > 0.0)
    {
        float3 wpos = normalize(viewPos);
        // Guard against division by near-zero y component
        float safeWposY = max(wpos.y, 0.001);

        [loop]
        for (int i = 0; i < AURORA_SAMPLES; i++)
        {
            float3 planeCoord = wpos * ((8.0 + currentStep * 7.0) / safeWposY) * 0.004;

            float2 coord = ctx.cameraPosition.xz * 0.00004 + planeCoord.xz;
            coord += float2(coord.y, -coord.x) * 0.3;

            float noise = auroraSample(coord, wind, VoU);

            if (noise > 0.0)
            {
                noise *= cloudNoise2D(coord * 0.125 + wind * 0.25);
                noise *= 0.5 * cloudNoise2D(coord + wind * 16.0) + 0.75;
                noise = noise * noise * 3.0 * sampleStep;
                noise *= max(sqrt(1.0 - length(planeCoord.xz) * 3.75), 0.0);

                float3 auroraColor = lerp(AURORA_LOW_COLOR, AURORA_HIGH_COLOR, pow(max(currentStep, 0.0001), 0.4));
                aurora += noise * auroraColor * exp2(-6.0 * float(i) * sampleStep);
            }
            currentStep += sampleStep;
        }
    }

    // Height-based visibility
    visibility *= saturate((ctx.cameraPosition.y + 70.0) / 8.0);

    return aurora * visibility * AURORA_INTENSITY;
}

// =============================================================================
// CIRRUS CLOUDS
// =============================================================================

float3 renderCirrusClouds(float3 rayDir, float3 sunDir, float3 sunColor, float time)
{
    if (rayDir.y <= 0.01)
        return float3(0.0, 0.0, 0.0);

    // Project onto high-altitude cirrus plane
    float t = CIRRUS_HEIGHT / rayDir.y;
    float2 hitPos = rayDir.xz * t;

    // Animate with wind
    float2 animUV = hitPos * 0.0002 + CLOUD_WIND_DIRECTION * time * 0.001;

    // Multi-octave noise for wispy pattern
    float density = 0.0;
    density += cloudNoise2D(animUV * 1.0) * 0.5;
    density += cloudNoise2D(animUV * 2.0 + 100.0) * 0.25;
    density += cloudNoise2D(animUV * 4.0 + 200.0) * 0.125;

    // Threshold for wispy appearance
    density = smoothstep(0.35, 0.65, density) * 0.3;

    if (density < 0.001)
        return float3(0.0, 0.0, 0.0);

    // Simple Henyey-Greenstein phase function
    float cosTheta = dot(rayDir, sunDir);
    float g = 0.5;
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    float phase = (1.0 - g2) / (4.0 * 3.14159265 * pow(max(denom, 0.0001), 1.5));

    float3 cirrusColor = lerp(float3(1.0, 1.0, 1.0), sunColor, 0.3);

    return cirrusColor * phase * density * CIRRUS_OPACITY;
}

// =============================================================================
// MAIN CLOUD RENDERING FUNCTION
// =============================================================================

CloudOutput renderVolumetricClouds(float3 rayOrigin, float3 rayDir, float3 sunDir,
                                   float3 sunColor, float time, float maxDistance)
{
    CloudContext ctx;
    ctx.cameraPosition = rayOrigin;
    ctx.sunDirection = sunDir;
    ctx.lightColor = sunColor;
    ctx.ambientColor = sunColor * 0.3;
    ctx.time = time;
    ctx.rainStrength = 0.0;
    ctx.shadowFade = saturate(sunDir.y * 2.0 + 0.5);
    ctx.sunVisibility = saturate(sunDir.y + 0.1);
    ctx.moonVisibility = saturate(-sunDir.y + 0.1);
    ctx.fogDensity = 1.0;

    float dither = frac(52.9829189 * frac(dot(rayDir.xz, float2(0.06711056, 0.00583715))));

#if CLOUD_MODE == 0
    // Skybox mode
    return drawCloudSkybox(rayDir * 1000.0, 1.0, dither, ctx, false);
#else
    // Volumetric mode
    float cloudViewLength = maxDistance;
    return drawCloudVolumetric(rayDir * maxDistance, rayOrigin, 1.0, dither, ctx, cloudViewLength, false);
#endif
}

// Full context version
CloudOutput renderVolumetricCloudsWithContext(float3 rayOrigin, float3 rayDir,
                                               CloudContext ctx, float maxDistance)
{
    float dither = frac(52.9829189 * frac(dot(rayDir.xz, float2(0.06711056, 0.00583715))));

#if CLOUD_MODE == 0
    // Skybox mode
    return drawCloudSkybox(rayDir * 1000.0, 1.0, dither, ctx, false);
#else
    // Volumetric mode
    float cloudViewLength = maxDistance;
    return drawCloudVolumetric(rayDir * maxDistance, rayOrigin, 1.0, dither, ctx, cloudViewLength, false);
#endif
}

// =============================================================================
// CLOUD SHADOWS
// =============================================================================

float sampleCloudShadow(float3 worldPos, float3 sunDir, float time)
{
    float tToCloud = (CLOUD_HEIGHT - worldPos.y) / max(0.001, sunDir.y);

    if (tToCloud < 0.0)
        return 1.0;

    float3 cloudPos = worldPos + sunDir * tToCloud;

    // Sample cloud density at the shadow position
    float2 wind = float2(
        time * CLOUD_SPEED * 0.0005,
        sin(time * CLOUD_SPEED * 0.001) * 0.005
    ) * 0.667;

    float2 cloudCoord = cloudPos.xz / CLOUD_SCALE;
    float cloudGradient = 0.5;  // Middle of cloud layer
    float density = cloudSample(cloudCoord * 0.004 * CLOUD_STRETCH, wind, cloudGradient, 0.0, 0.0, 0.0);

    return 1.0 - saturate(density * CLOUD_SHADOW_STRENGTH * 5.0);
}

#endif // __OPENRTX_CLOUDS_HLSL__
