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
// Implements real-time ray-marched volumetric clouds with multiple layers
// Based on: "The Real-time Volumetric Cloudscapes of Horizon Zero Dawn"
// =============================================================================

#ifndef __OPENRTX_CLOUDS_HLSL__
#define __OPENRTX_CLOUDS_HLSL__

#include "Settings.hlsl"

// =============================================================================
// PHASE FUNCTIONS (use Sky.hlsl's phaseMieHG - already included via OpenRTX.hlsl)
// =============================================================================

// phaseMieHG is defined in Sky.hlsl which is included before Clouds.hlsl

// =============================================================================
// NOISE FUNCTIONS
// =============================================================================

// Hash functions for procedural noise
float hashCloud(float3 p)
{
    p = frac(p * 0.3183099 + 0.1);
    p *= 17.0;
    return frac(p.x * p.y * p.z * (p.x + p.y + p.z));
}

float hashCloud2(float2 p)
{
    return frac(sin(dot(p, float2(127.1, 311.7))) * 43758.5453);
}

// 3D Value noise
float valueNoise3D(float3 p)
{
    float3 pi = floor(p);
    float3 pf = frac(p);

    // Smoothstep interpolation
    float3 w = pf * pf * (3.0 - 2.0 * pf);

    float n000 = hashCloud(pi + float3(0, 0, 0));
    float n001 = hashCloud(pi + float3(0, 0, 1));
    float n010 = hashCloud(pi + float3(0, 1, 0));
    float n011 = hashCloud(pi + float3(0, 1, 1));
    float n100 = hashCloud(pi + float3(1, 0, 0));
    float n101 = hashCloud(pi + float3(1, 0, 1));
    float n110 = hashCloud(pi + float3(1, 1, 0));
    float n111 = hashCloud(pi + float3(1, 1, 1));

    float n00 = lerp(n000, n100, w.x);
    float n01 = lerp(n001, n101, w.x);
    float n10 = lerp(n010, n110, w.x);
    float n11 = lerp(n011, n111, w.x);

    float n0 = lerp(n00, n10, w.y);
    float n1 = lerp(n01, n11, w.y);

    return lerp(n0, n1, w.z);
}

// Worley/Cellular noise for cloud detail
float worleyNoise3D(float3 p)
{
    float3 pi = floor(p);
    float3 pf = frac(p);

    float minDist = 1.0;

    for (int x = -1; x <= 1; x++)
    {
        for (int y = -1; y <= 1; y++)
        {
            for (int z = -1; z <= 1; z++)
            {
                float3 offset = float3(x, y, z);
                float3 cellId = pi + offset;

                float3 cellCenter = offset + float3(
                    hashCloud(cellId),
                    hashCloud(cellId + 100.0),
                    hashCloud(cellId + 200.0));

                float dist = length(cellCenter - pf);
                minDist = min(minDist, dist);
            }
        }
    }

    return minDist;
}

// Fractal Brownian Motion
float fbm(float3 p, int octaves)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    [unroll]
    for (int i = 0; i < octaves; i++)
    {
        value += amplitude * valueNoise3D(p * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value;
}

// Worley FBM for erosion
float worleyFBM(float3 p, int octaves)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    [unroll]
    for (int i = 0; i < octaves; i++)
    {
        value += amplitude * (1.0 - worleyNoise3D(p * frequency));
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value;
}

// =============================================================================
// CLOUD DENSITY FUNCTIONS
// =============================================================================

// Height-based density gradient for cumulus clouds
float cloudHeightGradient(float height, float cloudType)
{
    // cloudType: 0 = stratus, 0.5 = stratocumulus, 1 = cumulus
    float4 gradientParams;

    // Stratus: flat, low
    float4 stratusParams = float4(0.02, 0.05, 0.1, 0.3);
    // Stratocumulus: puffy but wide
    float4 stratocumulusParams = float4(0.0, 0.1, 0.3, 0.6);
    // Cumulus: tall, puffy
    float4 cumulusParams = float4(0.0, 0.1, 0.7, 1.0);

    gradientParams = lerp(lerp(stratusParams, stratocumulusParams, cloudType * 2.0),
                          cumulusParams, saturate(cloudType * 2.0 - 1.0));

    float gradient = smoothstep(gradientParams.x, gradientParams.y, height) *
                     smoothstep(gradientParams.w, gradientParams.z, height);

    return gradient;
}

// Remap function
float remap(float value, float low1, float high1, float low2, float high2)
{
    return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}

// Sample cloud density at a point
float sampleCloudDensity(float3 worldPos, float time, bool detailed)
{
    // Convert world position to cloud space
    float cloudThickness = CLOUD_MAX_HEIGHT - CLOUD_MIN_HEIGHT;
    float heightFraction = (worldPos.y - CLOUD_MIN_HEIGHT) / cloudThickness;

    // Check bounds
    if (heightFraction < 0.0 || heightFraction > 1.0)
        return 0.0;

    // Wind animation
    float2 windOffset = CLOUD_WIND_DIRECTION * CLOUD_WIND_SPEED * time;
    float3 animatedPos = worldPos + float3(windOffset.x, 0.0, windOffset.y);

    // Base cloud shape (large scale)
    float3 baseCoords = animatedPos * 0.0003;
    float baseNoise = fbm(baseCoords, 4);

    // Coverage map (simplified - in production would use weather texture)
    float2 coverageCoords = animatedPos.xz * 0.00005 + time * 0.001;
    float coverage = valueNoise3D(float3(coverageCoords, time * 0.1));
    coverage = smoothstep(0.3, 0.7, coverage) * CLOUD_COVERAGE + (1.0 - CLOUD_COVERAGE) * 0.5;

    // Cloud type from coverage (simplified)
    float cloudType = saturate(coverage * 2.0);

    // Height gradient
    float heightGradient = cloudHeightGradient(heightFraction, cloudType);

    // Combine base shape
    float baseDensity = remap(baseNoise, 1.0 - coverage, 1.0, 0.0, 1.0);
    baseDensity *= heightGradient;
    baseDensity = saturate(baseDensity);

    // Early exit for low density
    if (baseDensity < 0.01)
        return 0.0;

    // Detail erosion (only when needed)
    if (detailed && baseDensity > 0.01)
    {
        float3 detailCoords = animatedPos * 0.003;

        // High frequency erosion
        float detailNoise = worleyFBM(detailCoords, 3);
        float detailErosion = detailNoise * CLOUD_DETAIL_STRENGTH;

        // Height-based detail
        float detailModifier = lerp(detailErosion, 1.0 - detailErosion, saturate(heightFraction * 5.0));

        baseDensity = remap(baseDensity, detailModifier, 1.0, 0.0, 1.0);
        baseDensity = saturate(baseDensity);
    }

    return baseDensity * CLOUD_DENSITY;
}

// =============================================================================
// PHASE FUNCTIONS FOR CLOUDS
// =============================================================================

// Dual-lobe Henyey-Greenstein phase function
float cloudPhaseFunction(float cosTheta)
{
    float g1 = CLOUD_PHASE_G1;
    float g2 = CLOUD_PHASE_G2;
    float blend = CLOUD_PHASE_BLEND;

    float phase1 = phaseMieHG(cosTheta, g1);
    float phase2 = phaseMieHG(cosTheta, g2);

    return lerp(phase2, phase1, blend);
}

// =============================================================================
// LIGHTING FUNCTIONS
// =============================================================================

// Beer-Lambert absorption
float beerLambert(float density)
{
    return exp(-density);
}

// Beer-Powder function for dark edges
float beerPowder(float density, float cosTheta)
{
    float beer = beerLambert(density);
    float powder = 1.0 - exp(-density * 2.0);

    // Mix based on angle to sun
    float depthProbability = lerp(powder, 1.0, saturate(cosTheta * 0.5 + 0.5));

    return beer * depthProbability * 2.0;
}

// Sample light transmittance through cloud (cone sample)
float sampleLightTransmittance(float3 pos, float3 lightDir, float time, int numSamples)
{
    float stepSize = (CLOUD_MAX_HEIGHT - pos.y) / float(numSamples);
    float density = 0.0;

    float3 currentPos = pos;

    [unroll]
    for (int i = 0; i < numSamples; i++)
    {
        currentPos += lightDir * stepSize;

        // Add some cone spread for softer shadows
        float3 offset = float3(
            hashCloud(currentPos + float(i)) - 0.5,
            0.0,
            hashCloud(currentPos + float(i) + 100.0) - 0.5) * stepSize * 0.5;

        density += sampleCloudDensity(currentPos + offset, time, false) * stepSize;
    }

    return beerLambert(density * CLOUD_SCATTERING_COEFFICIENT);
}

// =============================================================================
// RAY-CLOUD INTERSECTION
// =============================================================================

// Get intersection with cloud layer
float2 rayCloudLayerIntersect(float3 rayOrigin, float3 rayDir)
{
    // Intersect with two horizontal planes
    float tMin = (CLOUD_MIN_HEIGHT - rayOrigin.y) / rayDir.y;
    float tMax = (CLOUD_MAX_HEIGHT - rayOrigin.y) / rayDir.y;

    // Sort
    if (tMin > tMax)
    {
        float temp = tMin;
        tMin = tMax;
        tMax = temp;
    }

    // Clamp to positive
    tMin = max(0.0, tMin);
    tMax = max(0.0, tMax);

    return float2(tMin, tMax);
}

// =============================================================================
// MAIN CLOUD RENDERING
// =============================================================================

struct CloudOutput
{
    float3 color;
    float transmittance;
    float depth;
};

CloudOutput renderVolumetricClouds(float3 rayOrigin, float3 rayDir, float3 sunDir, float3 sunColor, float time, float maxDistance)
{
    CloudOutput output;
    output.color = 0.0;
    output.transmittance = 1.0;
    output.depth = maxDistance;

    // Check if clouds are enabled
#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    // Get intersection with cloud layer
    float2 intersection = rayCloudLayerIntersect(rayOrigin, rayDir);
    float tMin = intersection.x;
    float tMax = min(intersection.y, maxDistance);

    // No intersection
    if (tMax <= tMin)
        return output;

    // Ray marching parameters
    float rayLength = tMax - tMin;
    int numSteps = CLOUD_MARCH_STEPS;
    float baseStepSize = rayLength / float(numSteps);

    // Phase function
    float cosTheta = dot(rayDir, sunDir);
    float phase = cloudPhaseFunction(cosTheta);

    // Ambient light (from sky)
    float3 ambientColor = sunColor * 0.15;

    // March through cloud layer
    float3 currentPos = rayOrigin + rayDir * tMin;
    float currentT = tMin;
    float accumDensity = 0.0;

    for (int i = 0; i < numSteps && output.transmittance > 0.01; i++)
    {
        // Adaptive step size based on density
        float stepSize = baseStepSize;

        // Sample cloud density
        float density = sampleCloudDensity(currentPos, time, true);

        if (density > 0.001)
        {
            // Record first hit depth
            if (output.depth >= maxDistance)
                output.depth = currentT;

            // Light sampling
            float lightTransmittance = sampleLightTransmittance(currentPos, sunDir, time, CLOUD_LIGHT_MARCH_STEPS);

            // Beer-Powder for scattered light
            float scatterAmount = beerPowder(density * stepSize, cosTheta);

            // Combine lighting
            float3 lightColor = sunColor * lightTransmittance * phase;
            float3 scatteredLight = (lightColor + ambientColor) * scatterAmount;

            // Powder effect for dark edges
            float powder = 1.0 - exp(-density * stepSize * 2.0);
            scatteredLight *= lerp(1.0, powder, CLOUD_POWDER_STRENGTH);

            // Accumulate
            float transmittance = beerLambert(density * stepSize * CLOUD_SCATTERING_COEFFICIENT);
            output.color += scatteredLight * output.transmittance;
            output.transmittance *= transmittance;

            // Smaller steps in dense regions
            stepSize *= lerp(1.0, 0.5, density);
        }
        else
        {
            // Larger steps in empty regions
            stepSize *= 2.0;
        }

        // Advance
        stepSize = clamp(stepSize, baseStepSize * 0.25, baseStepSize * 4.0);
        currentT += stepSize;
        currentPos += rayDir * stepSize;

        // Check bounds
        if (currentT > tMax)
            break;
    }

    return output;
}

// =============================================================================
// CIRRUS CLOUD LAYER
// =============================================================================

float sampleCirrusCloud(float2 uv, float time)
{
    // Animated coordinates
    float2 animUV = uv + CLOUD_WIND_DIRECTION * time * 0.0001;

    // Multi-layer noise
    float cirrus = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    for (int i = 0; i < 3; i++)
    {
        float2 sampleUV = animUV * frequency * 0.0001;
        cirrus += amplitude * valueNoise3D(float3(sampleUV, float(i) * 0.5));
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    // Sharp threshold
    cirrus = smoothstep(0.4, 0.6, cirrus);

    return cirrus * 0.3; // Thin cirrus
}

float3 renderCirrusClouds(float3 rayDir, float3 sunDir, float3 sunColor, float time)
{
    if (rayDir.y <= 0.01)
        return 0.0;

    // Project onto cirrus plane
    float t = (CIRRUS_HEIGHT - 0.0) / rayDir.y;
    float2 hitPos = rayDir.xz * t;

    // Sample cirrus
    float density = sampleCirrusCloud(hitPos, time);

    if (density < 0.001)
        return 0.0;

    // Simple lighting
    float cosTheta = dot(rayDir, sunDir);
    float phase = phaseMieHG(cosTheta, 0.5);

    float3 color = sunColor * phase * density;

    // Forward scattering
    float forwardScatter = pow(saturate(cosTheta), 8.0);
    color += sunColor * forwardScatter * density * 0.5;

    return color;
}

// =============================================================================
// CLOUD SHADOWS
// =============================================================================

// Sample cloud shadow at a world position
float sampleCloudShadow(float3 worldPos, float3 sunDir, float time)
{
    // Project position up to cloud layer along sun direction
    float tToCloud = (CLOUD_MIN_HEIGHT - worldPos.y) / max(0.001, sunDir.y);

    if (tToCloud < 0.0)
        return 1.0; // Above clouds

    float3 cloudSamplePos = worldPos + sunDir * tToCloud;

    // Sample cloud density
    float density = sampleCloudDensity(cloudSamplePos, time, false);

    // Soft shadow
    return 1.0 - saturate(density * 10.0);
}

#endif // __OPENRTX_CLOUDS_HLSL__
