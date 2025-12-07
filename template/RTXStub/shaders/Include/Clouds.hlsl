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
// OpenRTX Minecraft-Style Volumetric Clouds
// Grid-based clouds matching Minecraft's blocky aesthetic with volumetric softening
// =============================================================================

#ifndef __OPENRTX_CLOUDS_HLSL__
#define __OPENRTX_CLOUDS_HLSL__

#include "Settings.hlsl"

// =============================================================================
// MINECRAFT CLOUD CONSTANTS
// =============================================================================

// Minecraft cloud cell size (12x12 blocks in XZ, 4 blocks thick)
static const float CLOUD_CELL_SIZE = 12.0;
static const float CLOUD_THICKNESS = 4.0;

// Cloud layer Y position (Minecraft uses Y=192 in modern versions)
static const float MC_CLOUD_HEIGHT = 192.0;

// Edge softening radius (in blocks) - keeps clouds mostly blocky
static const float CLOUD_SOFT_EDGE = 1.5;

// =============================================================================
// HASH FUNCTIONS
// =============================================================================

// Integer hash for deterministic cloud pattern
uint hashInt(uint x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

// 2D cell hash - determines if a cloud cell exists
float cellHash(int2 cell)
{
    uint h = hashInt(uint(cell.x) + hashInt(uint(cell.y) * 15731U));
    return float(h) / float(0xFFFFFFFFU);
}

// Simple 3D hash for internal variation
float hash3D(float3 p)
{
    p = frac(p * 0.3183099 + 0.1);
    p *= 17.0;
    return frac(p.x * p.y * p.z * (p.x + p.y + p.z));
}

// =============================================================================
// MINECRAFT-STYLE CLOUD PATTERN
// =============================================================================

// Check if a cell contains a cloud (mimics Minecraft's cloud pattern)
bool hasCloud(int2 cell, float time)
{
    // Wind animation - shift cells over time
    float2 windOffset = CLOUD_WIND_DIRECTION * CLOUD_WIND_SPEED * time * 0.01;
    int2 shiftedCell = cell - int2(floor(windOffset));

    // Deterministic pattern based on cell position
    float h = cellHash(shiftedCell);

    // Cloud coverage threshold
    return h < CLOUD_COVERAGE;
}

// Get soft density within a cloud cell (1.0 at center, 0.0 at edges)
float getCloudCellDensity(float3 localPos, float cellSize, float thickness)
{
    // localPos is position within the cell (0 to cellSize in XZ, 0 to thickness in Y)

    // Calculate distance from edges
    float3 halfSize = float3(cellSize * 0.5, thickness * 0.5, cellSize * 0.5);
    float3 center = halfSize;
    float3 dist = abs(localPos - center);

    // Soft edges only - keep interior solid
    float3 edgeDist = halfSize - dist;

    // Apply softening only near edges
    float softX = saturate(edgeDist.x / CLOUD_SOFT_EDGE);
    float softY = saturate(edgeDist.y / CLOUD_SOFT_EDGE);
    float softZ = saturate(edgeDist.z / CLOUD_SOFT_EDGE);

    // Smooth the edge falloff
    softX = smoothstep(0.0, 1.0, softX);
    softY = smoothstep(0.0, 1.0, softY);
    softZ = smoothstep(0.0, 1.0, softZ);

    return softX * softY * softZ;
}

// Sample cloud density at world position
float sampleMinecraftCloudDensity(float3 worldPos, float time)
{
    // Check if within cloud layer height
    float cloudBottom = MC_CLOUD_HEIGHT;
    float cloudTop = MC_CLOUD_HEIGHT + CLOUD_THICKNESS;

    if (worldPos.y < cloudBottom || worldPos.y > cloudTop)
        return 0.0;

    // Apply wind offset to world position
    float2 windOffset = CLOUD_WIND_DIRECTION * CLOUD_WIND_SPEED * time * 0.01;
    float3 animatedPos = worldPos;
    animatedPos.xz += windOffset;

    // Get cell coordinates
    int2 cell = int2(floor(animatedPos.xz / CLOUD_CELL_SIZE));

    // Check if this cell has a cloud
    if (!hasCloud(cell, 0.0))  // time=0 since we already applied wind to position
        return 0.0;

    // Get local position within cell
    float3 localPos;
    localPos.xz = fmod(animatedPos.xz, CLOUD_CELL_SIZE);
    if (localPos.x < 0.0) localPos.x += CLOUD_CELL_SIZE;
    if (localPos.z < 0.0) localPos.z += CLOUD_CELL_SIZE;
    localPos.y = worldPos.y - cloudBottom;

    // Get density with soft edges
    float density = getCloudCellDensity(localPos, CLOUD_CELL_SIZE, CLOUD_THICKNESS);

    return density * CLOUD_DENSITY;
}

// =============================================================================
// PHASE FUNCTIONS
// =============================================================================

// Henyey-Greenstein phase function
float cloudPhaseHG(float cosTheta, float g)
{
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return (1.0 - g2) / (4.0 * 3.14159265 * pow(max(denom, 0.0001), 1.5));
}

// Dual-lobe phase function for clouds
float cloudPhaseFunction(float cosTheta)
{
    float forward = cloudPhaseHG(cosTheta, 0.8);   // Strong forward scattering
    float back = cloudPhaseHG(cosTheta, -0.3);     // Weak back scattering
    return lerp(back, forward, 0.7);
}

// =============================================================================
// LIGHTING
// =============================================================================

// Beer-Lambert absorption
float beerLambert(float opticalDepth)
{
    return exp(-opticalDepth);
}

// Sample light transmittance toward sun
float sampleLightTransmittance(float3 pos, float3 lightDir, float time, int numSamples)
{
    float cloudTop = MC_CLOUD_HEIGHT + CLOUD_THICKNESS;
    float maxDist = (cloudTop - pos.y) / max(lightDir.y, 0.001);
    maxDist = min(maxDist, CLOUD_CELL_SIZE * 3.0);  // Limit sample distance

    float stepSize = maxDist / float(numSamples);
    float density = 0.0;
    float3 samplePos = pos;

    [loop]
    for (int i = 0; i < numSamples; i++)
    {
        samplePos += lightDir * stepSize;
        density += sampleMinecraftCloudDensity(samplePos, time) * stepSize;
    }

    return beerLambert(density * CLOUD_SCATTERING_COEFFICIENT);
}

// =============================================================================
// RAY-BOX INTERSECTION
// =============================================================================

// Intersect ray with cloud layer bounds
float2 rayCloudLayerIntersect(float3 rayOrigin, float3 rayDir)
{
    float cloudBottom = MC_CLOUD_HEIGHT;
    float cloudTop = MC_CLOUD_HEIGHT + CLOUD_THICKNESS;

    float tMin = (cloudBottom - rayOrigin.y) / rayDir.y;
    float tMax = (cloudTop - rayOrigin.y) / rayDir.y;

    if (tMin > tMax)
    {
        float temp = tMin;
        tMin = tMax;
        tMax = temp;
    }

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

#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    // Intersect with cloud layer
    float2 intersection = rayCloudLayerIntersect(rayOrigin, rayDir);
    float tMin = intersection.x;
    float tMax = min(intersection.y, maxDistance);

    if (tMax <= tMin || tMin > 20000.0)
        return output;

    // Ray march parameters - fewer steps since clouds are blocky
    float rayLength = tMax - tMin;
    int numSteps = min(CLOUD_MARCH_STEPS, int(rayLength / 2.0) + 8);
    float stepSize = rayLength / float(numSteps);

    // Jitter to reduce banding
    float2 jitterSeed = rayDir.xz * 1000.0 + time * 10.0;
    float jitter = frac(52.9829189 * frac(dot(jitterSeed, float2(0.06711056, 0.00583715))));

    // Phase function
    float cosTheta = dot(rayDir, sunDir);
    float phase = cloudPhaseFunction(cosTheta);

    // Time of day scaling
    float dayFactor = saturate(sunDir.y * 2.0 + 0.5);
    float timeScale = max(0.05, dayFactor);

    // Ambient from sky
    float3 ambientColor = sunColor * 0.2 * timeScale;

    // March through cloud layer
    float3 currentPos = rayOrigin + rayDir * (tMin + jitter * stepSize);
    float currentT = tMin;

    [loop]
    for (int i = 0; i < numSteps; i++)
    {
        if (output.transmittance < 0.01)
            break;

        float density = sampleMinecraftCloudDensity(currentPos, time);

        if (density > 0.001)
        {
            // Record first hit
            if (output.depth >= maxDistance)
                output.depth = currentT;

            // Light sampling
            float lightTrans = sampleLightTransmittance(currentPos, sunDir, time, CLOUD_LIGHT_MARCH_STEPS);

            // Direct + scattered light
            float3 directLight = sunColor * phase * lightTrans;

            // Multi-scatter approximation
            float3 scatterLight = sunColor * 0.25 * pow(lightTrans, 0.3) * 0.4;

            // Silver lining when backlit
            float silverLining = pow(saturate(-cosTheta), 4.0) * lightTrans * 0.3;

            float3 totalLight = directLight + scatterLight + ambientColor + sunColor * silverLining;

            // Accumulate
            float absorption = density * stepSize * CLOUD_SCATTERING_COEFFICIENT;
            float transmittance = beerLambert(absorption);

            output.color += totalLight * (1.0 - transmittance) * output.transmittance;
            output.transmittance *= transmittance;
        }

        currentT += stepSize;
        currentPos += rayDir * stepSize;

        if (currentT > tMax)
            break;
    }

    return output;
}

// =============================================================================
// CLOUD SHADOWS
// =============================================================================

float sampleCloudShadow(float3 worldPos, float3 sunDir, float time)
{
    // Project to cloud layer
    float tToCloud = (MC_CLOUD_HEIGHT - worldPos.y) / max(0.001, sunDir.y);

    if (tToCloud < 0.0)
        return 1.0;

    float3 cloudPos = worldPos + sunDir * tToCloud;
    float density = sampleMinecraftCloudDensity(cloudPos, time);

    return 1.0 - saturate(density * 5.0);
}

// =============================================================================
// CIRRUS CLOUDS (High-altitude wispy layer - 2D, not volumetric)
// =============================================================================

// Simple 2D noise for cirrus
float cirrusNoise2D(float2 p)
{
    float2 pi = floor(p);
    float2 pf = frac(p);
    float2 w = pf * pf * (3.0 - 2.0 * pf);

    float n00 = hash3D(float3(pi, 0.0));
    float n01 = hash3D(float3(pi + float2(0, 1), 0.0));
    float n10 = hash3D(float3(pi + float2(1, 0), 0.0));
    float n11 = hash3D(float3(pi + float2(1, 1), 0.0));

    return lerp(lerp(n00, n10, w.x), lerp(n01, n11, w.x), w.y);
}

float3 renderCirrusClouds(float3 rayDir, float3 sunDir, float3 sunColor, float time)
{
    if (rayDir.y <= 0.01)
        return 0.0;

    // Project onto high-altitude cirrus plane
    float cirrusHeight = 400.0;  // High above blocky clouds
    float t = cirrusHeight / rayDir.y;
    float2 hitPos = rayDir.xz * t;

    // Animate with wind
    float2 animUV = hitPos * 0.0002 + CLOUD_WIND_DIRECTION * time * 0.001;

    // Multi-octave noise for wispy pattern
    float density = 0.0;
    density += cirrusNoise2D(animUV * 1.0) * 0.5;
    density += cirrusNoise2D(animUV * 2.0 + 100.0) * 0.25;
    density += cirrusNoise2D(animUV * 4.0 + 200.0) * 0.125;

    // Threshold for wispy appearance
    density = smoothstep(0.35, 0.65, density) * 0.3;

    if (density < 0.001)
        return 0.0;

    // Simple lighting
    float cosTheta = dot(rayDir, sunDir);
    float phase = cloudPhaseHG(cosTheta, 0.5);

    float3 cirrusColor = lerp(float3(1.0, 1.0, 1.0), sunColor, 0.3);

    return cirrusColor * phase * density * CIRRUS_OPACITY;
}

#endif // __OPENRTX_CLOUDS_HLSL__
