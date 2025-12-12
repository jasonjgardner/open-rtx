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
// Physically-inspired raymarched volumetric clouds with tileable 3D noise
// =============================================================================

#ifndef __OPENRTX_CLOUDS_HLSL__
#define __OPENRTX_CLOUDS_HLSL__

#include "Settings.hlsl"

// =============================================================================
// CLOUD VOLUME CONSTANTS
// =============================================================================

// Cloud layer bounds (world-space Y coordinates)
static const float kCloudLayerBottom = CLOUD_MIN_HEIGHT;
static const float kCloudLayerTop = CLOUD_MAX_HEIGHT;
static const float kCloudLayerThickness = kCloudLayerTop - kCloudLayerBottom;

// =============================================================================
// HASH FUNCTIONS FOR PROCEDURAL NOISE
// =============================================================================

// High-quality 3D hash function (avoids sin() precision issues)
float hash31(float3 p)
{
    p = frac(p * float3(0.1031, 0.1030, 0.0973));
    p += dot(p, p.yxz + 33.33);
    return frac((p.x + p.y) * p.z);
}

// 3D hash returning float3 - uses integer-like operations for better quality
float3 hash33(float3 p)
{
    p = frac(p * float3(0.1031, 0.1030, 0.0973));
    p += dot(p, p.yxz + 33.33);
    return frac((p.xxy + p.yxx) * p.zyx);
}

// Integer hash for consistent random
uint hashInt(uint x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

// =============================================================================
// TILEABLE 3D NOISE FUNCTIONS
// =============================================================================

// Smooth interpolation for noise (quintic for C2 continuity)
float3 quinticInterp(float3 x)
{
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

// Helper: wrap coordinates for tiling with proper negative handling
float3 wrapCoords(float3 p, float period)
{
    return p - floor(p / period) * period;
}

// Gradient noise (Perlin-like) - base shape noise
// Tileable with period for seamless texturing
float gradientNoise3D(float3 p, float period)
{
    // Proper tiling that handles negatives correctly
    float3 pTiled = wrapCoords(p, period);

    float3 i = floor(pTiled);
    float3 f = frac(pTiled);

    // Smooth interpolation
    float3 u = quinticInterp(f);

    // Gradient vectors at corners (using improved hash)
    float3 ga = hash33(wrapCoords(i + float3(0, 0, 0), period)) * 2.0 - 1.0;
    float3 gb = hash33(wrapCoords(i + float3(1, 0, 0), period)) * 2.0 - 1.0;
    float3 gc = hash33(wrapCoords(i + float3(0, 1, 0), period)) * 2.0 - 1.0;
    float3 gd = hash33(wrapCoords(i + float3(1, 1, 0), period)) * 2.0 - 1.0;
    float3 ge = hash33(wrapCoords(i + float3(0, 0, 1), period)) * 2.0 - 1.0;
    float3 gf = hash33(wrapCoords(i + float3(1, 0, 1), period)) * 2.0 - 1.0;
    float3 gg = hash33(wrapCoords(i + float3(0, 1, 1), period)) * 2.0 - 1.0;
    float3 gh = hash33(wrapCoords(i + float3(1, 1, 1), period)) * 2.0 - 1.0;

    // Normalize gradients for consistency
    ga = normalize(ga + 0.0001);
    gb = normalize(gb + 0.0001);
    gc = normalize(gc + 0.0001);
    gd = normalize(gd + 0.0001);
    ge = normalize(ge + 0.0001);
    gf = normalize(gf + 0.0001);
    gg = normalize(gg + 0.0001);
    gh = normalize(gh + 0.0001);

    // Distance vectors
    float3 pa = f - float3(0, 0, 0);
    float3 pb = f - float3(1, 0, 0);
    float3 pc = f - float3(0, 1, 0);
    float3 pd = f - float3(1, 1, 0);
    float3 pe = f - float3(0, 0, 1);
    float3 pf = f - float3(1, 0, 1);
    float3 pg = f - float3(0, 1, 1);
    float3 ph = f - float3(1, 1, 1);

    // Dot products
    float va = dot(ga, pa);
    float vb = dot(gb, pb);
    float vc = dot(gc, pc);
    float vd = dot(gd, pd);
    float ve = dot(ge, pe);
    float vf = dot(gf, pf);
    float vg = dot(gg, pg);
    float vh = dot(gh, ph);

    // Trilinear interpolation
    return lerp(lerp(lerp(va, vb, u.x), lerp(vc, vd, u.x), u.y),
                lerp(lerp(ve, vf, u.x), lerp(vg, vh, u.x), u.y), u.z);
}

// Value noise - simpler, faster
float valueNoise3D(float3 p, float period)
{
    float3 pTiled = wrapCoords(p, period);

    float3 i = floor(pTiled);
    float3 f = frac(pTiled);

    float3 u = quinticInterp(f);

    // Random values at corners
    float n000 = hash31(wrapCoords(i + float3(0, 0, 0), period));
    float n100 = hash31(wrapCoords(i + float3(1, 0, 0), period));
    float n010 = hash31(wrapCoords(i + float3(0, 1, 0), period));
    float n110 = hash31(wrapCoords(i + float3(1, 1, 0), period));
    float n001 = hash31(wrapCoords(i + float3(0, 0, 1), period));
    float n101 = hash31(wrapCoords(i + float3(1, 0, 1), period));
    float n011 = hash31(wrapCoords(i + float3(0, 1, 1), period));
    float n111 = hash31(wrapCoords(i + float3(1, 1, 1), period));

    return lerp(lerp(lerp(n000, n100, u.x), lerp(n010, n110, u.x), u.y),
                lerp(lerp(n001, n101, u.x), lerp(n011, n111, u.x), u.y), u.z);
}

// Worley (cellular) noise - for cloud edge breakup
float worleyNoise3D(float3 p, float period)
{
    float3 pTiled = wrapCoords(p, period);

    float3 i = floor(pTiled);
    float3 f = frac(pTiled);

    float minDist = 1.0;

    // Check 3x3x3 neighborhood
    [unroll]
    for (int x = -1; x <= 1; x++)
    {
        [unroll]
        for (int y = -1; y <= 1; y++)
        {
            [unroll]
            for (int z = -1; z <= 1; z++)
            {
                float3 neighbor = float3(x, y, z);
                float3 cellPos = wrapCoords(i + neighbor, period);

                // Random position within this cell
                float3 cellPoint = hash33(cellPos);
                float3 diff = neighbor + cellPoint - f;
                float dist = dot(diff, diff);
                minDist = min(minDist, dist);
            }
        }
    }

    return sqrt(minDist);
}

// =============================================================================
// FRACTAL BROWNIAN MOTION (FBM) NOISE
// =============================================================================

// Low-frequency FBM for cloud coverage/shape
float fbmLowFreq(float3 p, int octaves, float period)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;
    float totalAmplitude = 0.0;

    [loop]
    for (int i = 0; i < octaves && i < 6; i++)
    {
        value += amplitude * (gradientNoise3D(p * frequency, period * frequency) * 0.5 + 0.5);
        totalAmplitude += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value / max(totalAmplitude, 0.001);
}

// High-frequency FBM with Worley for cloud detail/erosion
float fbmHighFreq(float3 p, int octaves, float period)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;
    float totalAmplitude = 0.0;

    [loop]
    for (int i = 0; i < octaves && i < 4; i++)
    {
        // Mix Worley and value noise for interesting detail
        float worley = worleyNoise3D(p * frequency, period * frequency);
        float perlin = valueNoise3D(p * frequency, period * frequency);
        float mixed = worley * 0.6 + perlin * 0.4;

        value += amplitude * mixed;
        totalAmplitude += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value / max(totalAmplitude, 0.001);
}

// =============================================================================
// HEIGHT-BASED DENSITY PROFILE
// =============================================================================

// Cloud density gradient based on height within cloud layer
// Creates flat dense bottoms and soft tapered tops (cumulus-like shape)
float cloudHeightGradient(float heightFraction)
{
    // Clamp to valid range
    heightFraction = saturate(heightFraction);

    // Bottom gradient: ramps up from 0 to 1 across bottom portion
    float bottomGrad = saturate((heightFraction - 0.0) / CLOUD_HEIGHT_GRADIENT_BOTTOM);

    // Top gradient: ramps down from 1 to 0 across top portion
    float topGrad = saturate((1.0 - heightFraction) / (1.0 - CLOUD_HEIGHT_GRADIENT_TOP));

    // Combine with smooth transitions
    float gradient = bottomGrad * topGrad;

    // Apply anvil effect (clouds spread more at top)
    float anvilFactor = lerp(1.0, 1.0 + heightFraction * 0.5, CLOUD_ANVIL_BIAS * heightFraction);

    return gradient * anvilFactor;
}

// Alternative gradient type selector for variety
float cloudDensityType(float heightFraction, int type)
{
    switch(type)
    {
        case 0: // Cumulus (default) - flat bottom, puffy top
            return cloudHeightGradient(heightFraction);

        case 1: // Stratocumulus - flatter, more uniform
        {
            float bottom = smoothstep(0.0, 0.2, heightFraction);
            float top = smoothstep(1.0, 0.7, heightFraction);
            return bottom * top;
        }

        case 2: // Cumulonimbus - tall towers
        {
            float bottom = smoothstep(0.0, 0.15, heightFraction);
            float top = smoothstep(1.0, 0.6, heightFraction);
            float core = 1.0 - abs(heightFraction - 0.4) * 0.8;
            return bottom * top * saturate(core);
        }

        default:
            return cloudHeightGradient(heightFraction);
    }
}

// =============================================================================
// CLOUD DENSITY SAMPLING
// =============================================================================

// Sample the full cloud density field at a world position
float sampleCloudDensity(float3 worldPos, float time, float lod)
{
    // Check if within cloud layer bounds
    if (worldPos.y < kCloudLayerBottom || worldPos.y > kCloudLayerTop)
        return 0.0;

    // Calculate height fraction within cloud layer (0 at bottom, 1 at top)
    float heightFraction = (worldPos.y - kCloudLayerBottom) / kCloudLayerThickness;

    // Apply wind advection (animate clouds)
    // Scale time appropriately - game time can be large, so use modulo for stability
    // The 0.05 factor converts to reasonable visual speed
    float windTime = fmod(time, 10000.0) * 0.05;
    float2 windDir = normalize(CLOUD_WIND_DIRECTION);
    float2 windOffset = windDir * CLOUD_WIND_SPEED * windTime;
    float3 animatedPos = worldPos;
    animatedPos.xz += windOffset;

    // Also add slight vertical wind shear (higher clouds move faster)
    animatedPos.xz += windOffset * heightFraction * 0.3;

    // Noise tiling period (for seamless wrapping)
    float period = 256.0;

    // ===================
    // Base shape noise (low frequency)
    // ===================
    float3 baseCoord = animatedPos * CLOUD_BASE_NOISE_FREQ;

    // Multi-octave base noise with LOD adjustment
    int baseOctaves = max(2, 4 - int(lod));
    float baseNoise = fbmLowFreq(baseCoord, baseOctaves, period);

    // Apply coverage remapping
    // coverage controls how much of the noise becomes cloud
    float coverageNoise = baseNoise;
    float coverage = CLOUD_COVERAGE + CLOUD_COVERAGE_OFFSET;

    // Remap noise based on coverage (higher coverage = more clouds)
    float baseDensity = saturate((coverageNoise - (1.0 - coverage)) / coverage);

    // Apply height gradient (cumulus shape)
    float heightGrad = cloudHeightGradient(heightFraction);
    baseDensity *= heightGrad;

    // Early out if base density is negligible
    if (baseDensity < 0.01)
        return 0.0;

    // ===================
    // Detail noise (high frequency erosion)
    // ===================
    float3 detailCoord = animatedPos * CLOUD_DETAIL_NOISE_FREQ;

    // Animate detail noise slightly differently for variation
    detailCoord.xz += windDir * CLOUD_WIND_SPEED * windTime * 0.2;

    // Detail noise - erodes cloud edges
    int detailOctaves = max(1, 3 - int(lod * 0.5));
    float detailNoise = fbmHighFreq(detailCoord, detailOctaves, period);

    // Apply detail erosion (subtract from base density)
    // More erosion at edges (lower base density)
    float detailAmount = CLOUD_DETAIL_STRENGTH * (1.0 - baseDensity * 0.5);
    float finalDensity = saturate(baseDensity - detailNoise * detailAmount);

    // Apply overall density multiplier
    finalDensity *= CLOUD_DENSITY * CLOUD_DENSITY_MULTIPLIER;

    return finalDensity;
}

// Cheaper density sample for shadow rays (fewer octaves)
float sampleCloudDensityCheap(float3 worldPos, float time)
{
    if (worldPos.y < kCloudLayerBottom || worldPos.y > kCloudLayerTop)
        return 0.0;

    float heightFraction = (worldPos.y - kCloudLayerBottom) / kCloudLayerThickness;

    // Use same wind time scaling as main density function
    float windTime = fmod(time, 10000.0) * 0.05;
    float2 windOffset = normalize(CLOUD_WIND_DIRECTION) * CLOUD_WIND_SPEED * windTime;
    float3 animatedPos = worldPos;
    animatedPos.xz += windOffset;

    float period = 256.0;
    float3 baseCoord = animatedPos * CLOUD_BASE_NOISE_FREQ;

    // Single octave for speed
    float baseNoise = gradientNoise3D(baseCoord, period) * 0.5 + 0.5;

    float coverage = CLOUD_COVERAGE + CLOUD_COVERAGE_OFFSET;
    float baseDensity = saturate((baseNoise - (1.0 - coverage)) / coverage);
    baseDensity *= cloudHeightGradient(heightFraction);

    return baseDensity * CLOUD_DENSITY * CLOUD_DENSITY_MULTIPLIER * 0.5;
}

// =============================================================================
// RAY-BOX INTERSECTION
// =============================================================================

// Intersect ray with axis-aligned bounding box (cloud layer)
// Returns (tEntry, tExit) where negative values indicate no hit
float2 rayAABBIntersect(float3 rayOrigin, float3 rayDir, float3 boxMin, float3 boxMax)
{
    float3 invDir = 1.0 / (rayDir + sign(rayDir) * 0.0001);

    float3 t0 = (boxMin - rayOrigin) * invDir;
    float3 t1 = (boxMax - rayOrigin) * invDir;

    float3 tmin = min(t0, t1);
    float3 tmax = max(t0, t1);

    float tEntry = max(max(tmin.x, tmin.y), tmin.z);
    float tExit = min(min(tmax.x, tmax.y), tmax.z);

    // Clamp entry to not be behind camera
    tEntry = max(0.0, tEntry);

    // No intersection if exit is before entry
    if (tExit < tEntry)
        return float2(-1.0, -1.0);

    return float2(tEntry, tExit);
}

// Intersect with horizontal cloud layer (infinite in XZ, bounded in Y)
float2 rayCloudLayerIntersect(float3 rayOrigin, float3 rayDir)
{
    // For nearly horizontal rays, use a large XZ extent
    float xzExtent = 100000.0;

    float3 boxMin = float3(-xzExtent, kCloudLayerBottom, -xzExtent);
    float3 boxMax = float3(xzExtent, kCloudLayerTop, xzExtent);

    return rayAABBIntersect(rayOrigin, rayDir, boxMin, boxMax);
}

// =============================================================================
// PHASE FUNCTIONS
// =============================================================================

// Henyey-Greenstein phase function
// g > 0: forward scattering, g < 0: back scattering
float phaseHenyeyGreenstein(float cosTheta, float g)
{
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    denom = max(denom, 0.0001);
    return (1.0 - g2) / (4.0 * kPi * pow(denom, 1.5));
}

// Dual-lobe phase function for realistic cloud scattering
// Combines forward and backward scattering lobes
float cloudPhaseFunction(float cosTheta)
{
    // Primary forward lobe (silver lining effect)
    float forward = phaseHenyeyGreenstein(cosTheta, CLOUD_HG_ASYMMETRY);

    // Secondary backward lobe (subtle backlighting)
    float backward = phaseHenyeyGreenstein(cosTheta, -0.3);

    // Blend the two lobes
    return lerp(backward, forward, CLOUD_PHASE_BLEND);
}

// =============================================================================
// BEER-LAMBERT EXTINCTION
// =============================================================================

// Beer-Lambert law for transmittance
float beerLambert(float opticalDepth)
{
    return exp(-opticalDepth);
}

// Beer-Lambert with powder effect (darkening at cloud base)
float beerLambertPowder(float opticalDepth, float depthInCloud)
{
    float beer = exp(-opticalDepth);

    // Powder effect - darkens where light enters thick cloud
    float powder = 1.0 - exp(-opticalDepth * 2.0);

    return lerp(beer, beer * powder, CLOUD_POWDER_STRENGTH);
}

// =============================================================================
// SELF-SHADOWING (LIGHT MARCHING)
// =============================================================================

// March toward light to estimate self-shadowing
float sampleLightTransmittance(float3 pos, float3 lightDir, float time)
{
    // Calculate how far we need to march toward light
    float maxDist = CLOUD_SHADOW_MARCH_DISTANCE;

    // Limit by cloud layer exit
    float tToTop = (kCloudLayerTop - pos.y) / max(lightDir.y, 0.001);
    if (lightDir.y > 0.0)
        maxDist = min(maxDist, tToTop);

    int numSteps = CLOUD_SHADOW_MARCH_STEPS;
    float stepSize = maxDist / float(numSteps);

    float opticalDepth = 0.0;
    float3 samplePos = pos;

    [loop]
    for (int i = 0; i < numSteps; i++)
    {
        samplePos += lightDir * stepSize;

        // Sample density (use cheap version for performance)
        float density = sampleCloudDensityCheap(samplePos, time);
        opticalDepth += density * stepSize * CLOUD_EXTINCTION_COEFFICIENT;

        // Early out if transmittance is negligible
        if (opticalDepth > 4.0)
            break;
    }

    return beerLambert(opticalDepth);
}

// =============================================================================
// MAIN CLOUD RAYMARCHING
// =============================================================================

struct CloudOutput
{
    float3 color;        // Accumulated in-scattered light
    float transmittance; // Remaining view transmittance
    float depth;         // Depth to first cloud hit
};

CloudOutput renderVolumetricCloudsPhysical(
    float3 rayOrigin,
    float3 rayDir,
    float3 sunDir,
    float3 sunColor,
    float time,
    float maxDistance,
    float2 screenUV)
{
    CloudOutput output;
    output.color = 0.0;
    output.transmittance = 1.0;
    output.depth = maxDistance;

#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    // Intersect ray with cloud layer
    float2 intersection = rayCloudLayerIntersect(rayOrigin, rayDir);
    float tEntry = intersection.x;
    float tExit = min(intersection.y, maxDistance);

    // No intersection with cloud layer
    if (tExit <= tEntry || tEntry < 0.0)
        return output;

    // Limit march distance
    float marchDistance = tExit - tEntry;
    marchDistance = min(marchDistance, 20000.0);

    // Adaptive step count based on distance
    int numSteps = CLOUD_MARCH_STEPS;
    if (marchDistance < 1000.0)
        numSteps = max(16, numSteps / 2);

    float stepSize = marchDistance / float(numSteps);

    // ===================
    // Ray jitter for TAA / anti-banding
    // ===================
    float jitter = 0.0;
    // Blue noise-like jitter using screen position
    float2 jitterSeed = screenUV * 1000.0 + time * 17.0;
    jitter = frac(52.9829189 * frac(dot(jitterSeed, float2(0.06711056, 0.00583715))));
    jitter *= CLOUD_JITTER_STRENGTH;

    // Starting position with jitter
    float3 currentPos = rayOrigin + rayDir * (tEntry + jitter * stepSize);
    float currentT = tEntry + jitter * stepSize;

    // ===================
    // Pre-compute lighting values
    // ===================
    float cosTheta = dot(rayDir, sunDir);
    float phase = cloudPhaseFunction(cosTheta);

    // Time-of-day scaling
    float dayFactor = saturate(sunDir.y * 2.0 + 0.5);
    float timeScale = max(0.05, dayFactor);

    // Ambient light (sky contribution)
    float3 ambientLight = sunColor * CLOUD_AMBIENT_STRENGTH * timeScale;
    ambientLight += float3(0.1, 0.12, 0.15) * (1.0 - dayFactor); // Night ambient

    // LOD based on initial distance
    float lod = saturate(tEntry / 5000.0) * 2.0;

    // ===================
    // Ray marching loop
    // ===================
    float accumulatedOpticalDepth = 0.0;
    bool firstHit = true;

    [loop]
    for (int i = 0; i < numSteps; i++)
    {
        // Early out if nearly opaque
        if (output.transmittance < CLOUD_TRANSMITTANCE_THRESHOLD)
            break;

        // Sample cloud density
        float density = sampleCloudDensity(currentPos, time, lod);

        if (density > 0.001)
        {
            // Record first hit for depth
            if (firstHit)
            {
                output.depth = currentT;
                firstHit = false;
            }

            // ===================
            // Extinction (Beer-Lambert)
            // ===================
            float extinction = density * CLOUD_EXTINCTION_COEFFICIENT * stepSize;
            float transmittance = beerLambertPowder(extinction, accumulatedOpticalDepth);
            accumulatedOpticalDepth += extinction;

            // ===================
            // In-scattering calculation
            // ===================

            // Self-shadowing: march toward sun
            float lightTransmittance = sampleLightTransmittance(currentPos, sunDir, time);

            // Direct sunlight contribution
            float3 directLight = sunColor * phase * lightTransmittance * timeScale;

            // Multi-scatter approximation (brightens cloud interiors)
            float multiScatter = pow(lightTransmittance, 0.25) * CLOUD_MULTISCATTER_STRENGTH;
            float3 scatteredLight = sunColor * multiScatter * 0.3;

            // Silver lining effect (forward scattering at cloud edges)
            float edgeFactor = 1.0 - saturate(density * 10.0);
            float silverLining = pow(saturate(-cosTheta), 3.0) * lightTransmittance * edgeFactor * 0.5;

            // Total lighting
            float3 totalLight = directLight + scatteredLight + ambientLight + sunColor * silverLining;

            // ===================
            // Accumulate color
            // ===================
            float3 scatteredEnergy = totalLight * density * CLOUD_SCATTERING_COEFFICIENT;
            output.color += scatteredEnergy * output.transmittance * stepSize;

            // Update transmittance
            output.transmittance *= transmittance;
        }

        // Advance ray
        currentT += stepSize;
        currentPos += rayDir * stepSize;

        // Exit if past the cloud layer
        if (currentT > tExit)
            break;
    }

    return output;
}

// =============================================================================
// MINECRAFT-STYLE CLOUDS (LEGACY FALLBACK)
// =============================================================================

// Keep the Minecraft-style implementation for CLOUD_RENDERING_MODE == 0

static const float MC_CLOUD_CELL_SIZE = 12.0;
static const float MC_CLOUD_THICKNESS = 4.0;
static const float MC_CLOUD_HEIGHT = 192.0;
static const float MC_CLOUD_SOFT_EDGE = 1.5;

float cellHash(int2 cell)
{
    uint h = hashInt(uint(cell.x) + hashInt(uint(cell.y) * 15731U));
    return float(h) / float(0xFFFFFFFFU);
}

bool hasCloudMC(int2 cell, float time)
{
    float2 windOffset = normalize(CLOUD_WIND_DIRECTION) * CLOUD_WIND_SPEED * time * 0.01;
    int2 shiftedCell = cell - int2(floor(windOffset));
    float h = cellHash(shiftedCell);
    return h < CLOUD_COVERAGE;
}

float getCloudCellDensityMC(float3 localPos, float cellSize, float thickness)
{
    float3 halfSize = float3(cellSize * 0.5, thickness * 0.5, cellSize * 0.5);
    float3 center = halfSize;
    float3 dist = abs(localPos - center);
    float3 edgeDist = halfSize - dist;

    float softX = saturate(edgeDist.x / MC_CLOUD_SOFT_EDGE);
    float softY = saturate(edgeDist.y / MC_CLOUD_SOFT_EDGE);
    float softZ = saturate(edgeDist.z / MC_CLOUD_SOFT_EDGE);

    softX = smoothstep(0.0, 1.0, softX);
    softY = smoothstep(0.0, 1.0, softY);
    softZ = smoothstep(0.0, 1.0, softZ);

    return softX * softY * softZ;
}

float sampleMinecraftCloudDensity(float3 worldPos, float time)
{
    float cloudBottom = MC_CLOUD_HEIGHT;
    float cloudTop = MC_CLOUD_HEIGHT + MC_CLOUD_THICKNESS;

    if (worldPos.y < cloudBottom || worldPos.y > cloudTop)
        return 0.0;

    float2 windOffset = normalize(CLOUD_WIND_DIRECTION) * CLOUD_WIND_SPEED * time * 0.01;
    float3 animatedPos = worldPos;
    animatedPos.xz += windOffset;

    int2 cell = int2(floor(animatedPos.xz / MC_CLOUD_CELL_SIZE));

    if (!hasCloudMC(cell, 0.0))
        return 0.0;

    float3 localPos;
    localPos.xz = fmod(animatedPos.xz, MC_CLOUD_CELL_SIZE);
    if (localPos.x < 0.0) localPos.x += MC_CLOUD_CELL_SIZE;
    if (localPos.z < 0.0) localPos.z += MC_CLOUD_CELL_SIZE;
    localPos.y = worldPos.y - cloudBottom;

    float density = getCloudCellDensityMC(localPos, MC_CLOUD_CELL_SIZE, MC_CLOUD_THICKNESS);
    return density * CLOUD_DENSITY;
}

CloudOutput renderVolumetricCloudsMC(
    float3 rayOrigin,
    float3 rayDir,
    float3 sunDir,
    float3 sunColor,
    float time,
    float maxDistance)
{
    CloudOutput output;
    output.color = 0.0;
    output.transmittance = 1.0;
    output.depth = maxDistance;

#if !ENABLE_VOLUMETRIC_CLOUDS
    return output;
#endif

    float cloudBottom = MC_CLOUD_HEIGHT;
    float cloudTop = MC_CLOUD_HEIGHT + MC_CLOUD_THICKNESS;

    // Ray-layer intersection
    float2 intersection;
    if (abs(rayDir.y) < 0.0001)
    {
        if (rayOrigin.y >= cloudBottom && rayOrigin.y <= cloudTop)
            intersection = float2(0.0, 1000.0);
        else
            return output;
    }
    else
    {
        float tMin = (cloudBottom - rayOrigin.y) / rayDir.y;
        float tMax = (cloudTop - rayOrigin.y) / rayDir.y;
        if (tMin > tMax) { float temp = tMin; tMin = tMax; tMax = temp; }
        tMin = max(0.0, tMin);
        tMax = max(0.0, tMax);
        intersection = float2(tMin, tMax);
    }

    float tEntry = intersection.x;
    float tExit = min(intersection.y, maxDistance);

    if (tExit <= tEntry || tEntry > 20000.0)
        return output;

    float rayLength = tExit - tEntry;
    int numSteps = min(CLOUD_MARCH_STEPS, int(rayLength / 2.0) + 8);
    float stepSize = rayLength / float(numSteps);

    float2 jitterSeed = rayDir.xz * 1000.0 + time * 10.0;
    float jitter = frac(52.9829189 * frac(dot(jitterSeed, float2(0.06711056, 0.00583715))));

    float cosTheta = dot(rayDir, sunDir);
    float phase = cloudPhaseFunction(cosTheta);

    float dayFactor = saturate(sunDir.y * 2.0 + 0.5);
    float timeScale = max(0.05, dayFactor);
    float3 ambientColor = sunColor * 0.2 * timeScale;

    float3 currentPos = rayOrigin + rayDir * (tEntry + jitter * stepSize);
    float currentT = tEntry;

    [loop]
    for (int i = 0; i < numSteps; i++)
    {
        if (output.transmittance < 0.01)
            break;

        float density = sampleMinecraftCloudDensity(currentPos, time);

        if (density > 0.001)
        {
            if (output.depth >= maxDistance)
                output.depth = currentT;

            // Simple light sampling
            float lightTrans = 1.0;
            float3 lightPos = currentPos;
            float lightStepSize = 20.0;
            for (int j = 0; j < 4; j++)
            {
                lightPos += sunDir * lightStepSize;
                float lightDensity = sampleMinecraftCloudDensity(lightPos, time);
                lightTrans *= exp(-lightDensity * lightStepSize * CLOUD_SCATTERING_COEFFICIENT);
            }

            float3 directLight = sunColor * phase * lightTrans;
            float3 scatterLight = sunColor * 0.25 * pow(lightTrans, 0.3) * 0.4;
            float silverLining = pow(saturate(-cosTheta), 4.0) * lightTrans * 0.3;

            float3 totalLight = directLight + scatterLight + ambientColor + sunColor * silverLining;

            float absorption = density * stepSize * CLOUD_SCATTERING_COEFFICIENT;
            float transmittance = exp(-absorption);

            output.color += totalLight * (1.0 - transmittance) * output.transmittance;
            output.transmittance *= transmittance;
        }

        currentT += stepSize;
        currentPos += rayDir * stepSize;

        if (currentT > tExit)
            break;
    }

    return output;
}

// =============================================================================
// UNIFIED CLOUD RENDERING INTERFACE
// =============================================================================

// Main entry point - selects between rendering modes
CloudOutput renderVolumetricClouds(
    float3 rayOrigin,
    float3 rayDir,
    float3 sunDir,
    float3 sunColor,
    float time,
    float maxDistance)
{
#if CLOUD_RENDERING_MODE == 1
    // Physically-inspired volumetric clouds
    return renderVolumetricCloudsPhysical(rayOrigin, rayDir, sunDir, sunColor, time, maxDistance, float2(0, 0));
#else
    // Minecraft-style blocky clouds (legacy)
    return renderVolumetricCloudsMC(rayOrigin, rayDir, sunDir, sunColor, time, maxDistance);
#endif
}

// Version with screen UV for jitter
CloudOutput renderVolumetricCloudsWithUV(
    float3 rayOrigin,
    float3 rayDir,
    float3 sunDir,
    float3 sunColor,
    float time,
    float maxDistance,
    float2 screenUV)
{
#if CLOUD_RENDERING_MODE == 1
    return renderVolumetricCloudsPhysical(rayOrigin, rayDir, sunDir, sunColor, time, maxDistance, screenUV);
#else
    return renderVolumetricCloudsMC(rayOrigin, rayDir, sunDir, sunColor, time, maxDistance);
#endif
}

// =============================================================================
// CLOUD SHADOWS ON TERRAIN
// =============================================================================

// Sample cloud shadow at a world position (for terrain shadowing)
float sampleCloudShadow(float3 worldPos, float3 sunDir, float time)
{
#if !ENABLE_VOLUMETRIC_CLOUDS
    return 1.0;
#endif

    // Project position up to cloud layer along sun direction
    float tToCloudBase = (kCloudLayerBottom - worldPos.y) / max(sunDir.y, 0.001);

    if (tToCloudBase < 0.0 || sunDir.y <= 0.0)
        return 1.0;  // Sun below horizon or position above clouds

    float3 cloudPos = worldPos + sunDir * tToCloudBase;

#if CLOUD_RENDERING_MODE == 1
    // Sample density at projected position
    float density = sampleCloudDensity(cloudPos, time, 1.0);  // Use medium LOD

    // March a bit through cloud layer for better shadow estimate
    float shadowDensity = density;
    float marchDist = min(200.0, kCloudLayerThickness);
    float3 marchPos = cloudPos;

    [unroll]
    for (int i = 0; i < 3; i++)
    {
        marchPos += sunDir * (marchDist / 3.0);
        shadowDensity += sampleCloudDensityCheap(marchPos, time);
    }
    shadowDensity /= 4.0;

    // Convert to shadow factor
    float shadow = exp(-shadowDensity * CLOUD_SHADOW_STRENGTH * 5.0);
    return saturate(shadow);
#else
    // Minecraft-style shadow
    float density = sampleMinecraftCloudDensity(cloudPos, time);
    return 1.0 - saturate(density * CLOUD_SHADOW_STRENGTH * 5.0);
#endif
}

// =============================================================================
// CIRRUS CLOUDS (HIGH-ALTITUDE WISPY LAYER)
// =============================================================================

float cirrusNoise2D(float2 p)
{
    float2 pi = floor(p);
    float2 pf = frac(p);
    float2 w = pf * pf * (3.0 - 2.0 * pf);

    float n00 = hash31(float3(pi, 0.0));
    float n01 = hash31(float3(pi + float2(0, 1), 0.0));
    float n10 = hash31(float3(pi + float2(1, 0), 0.0));
    float n11 = hash31(float3(pi + float2(1, 1), 0.0));

    return lerp(lerp(n00, n10, w.x), lerp(n01, n11, w.x), w.y);
}

float3 renderCirrusClouds(float3 rayDir, float3 sunDir, float3 sunColor, float time)
{
    if (rayDir.y <= 0.01)
        return 0.0;

    float cirrusHeight = CIRRUS_HEIGHT;
    float t = cirrusHeight / rayDir.y;
    float2 hitPos = rayDir.xz * t;

    // Animate with wind (cirrus moves faster)
    float2 windDir = normalize(CLOUD_WIND_DIRECTION);
    float2 animUV = hitPos * 0.0002 + windDir * time * CLOUD_WIND_SPEED * 0.002;

    // Multi-octave noise for wispy appearance
    float density = 0.0;
    density += cirrusNoise2D(animUV * 1.0) * 0.5;
    density += cirrusNoise2D(animUV * 2.0 + 100.0) * 0.25;
    density += cirrusNoise2D(animUV * 4.0 + 200.0) * 0.125;

    // Threshold for wispy pattern
    density = smoothstep(0.35, 0.65, density) * 0.3;

    if (density < 0.001)
        return 0.0;

    // Simple lighting
    float cosTheta = dot(rayDir, sunDir);
    float phase = phaseHenyeyGreenstein(cosTheta, 0.5);

    float dayFactor = saturate(sunDir.y * 2.0 + 0.5);
    float3 cirrusColor = lerp(float3(0.8, 0.85, 1.0), sunColor, 0.3 * dayFactor);

    return cirrusColor * phase * density * CIRRUS_OPACITY;
}

#endif // __OPENRTX_CLOUDS_HLSL__
