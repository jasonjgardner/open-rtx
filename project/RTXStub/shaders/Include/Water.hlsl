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
// OpenRTX Water - Advanced Water Rendering
// Implements parallax wave mapping, volumetric scattering, and caustics
// Based on: "Vibrant Visuals" water customization and ocean rendering techniques
// =============================================================================

#ifndef __OPENRTX_WATER_HLSL__
#define __OPENRTX_WATER_HLSL__

#include "Settings.hlsl"

// =============================================================================
// WAVE FUNCTIONS
// =============================================================================

// Gerstner wave function
float3 gerstnerWave(float2 pos, float time, float2 direction, float wavelength, float steepness, float amplitude)
{
    float k = kTwoPi / wavelength;
    float c = sqrt(9.81 / k); // Phase speed
    float f = k * (dot(direction, pos) - c * time);

    float a = steepness / k;

    return float3(
        direction.x * a * cos(f),
        amplitude * sin(f),
        direction.y * a * cos(f));
}

// Sum of Gerstner waves
float3 sumGerstnerWaves(float2 pos, float time)
{
    float3 offset = 0.0;

    // Wave parameters: direction, wavelength, steepness, amplitude
    // Primary wave
    offset += gerstnerWave(pos, time * WATER_WAVE_SPEED, normalize(float2(1.0, 0.3)), 4.0, 0.15, 0.05 * WATER_WAVE_AMPLITUDE);
    // Secondary wave
    offset += gerstnerWave(pos, time * WATER_WAVE_SPEED * 1.1, normalize(float2(-0.7, 0.6)), 2.5, 0.12, 0.03 * WATER_WAVE_AMPLITUDE);
    // Tertiary wave
    offset += gerstnerWave(pos, time * WATER_WAVE_SPEED * 0.9, normalize(float2(0.3, -0.8)), 1.5, 0.1, 0.02 * WATER_WAVE_AMPLITUDE);
    // Detail waves
    offset += gerstnerWave(pos, time * WATER_WAVE_SPEED * 1.3, normalize(float2(-0.5, -0.5)), 0.8, 0.08, 0.01 * WATER_WAVE_AMPLITUDE);

    return offset;
}

// Hash functions for noise
float hashWater(float2 p)
{
    return frac(sin(dot(p, float2(127.1, 311.7))) * 43758.5453);
}

float hashWater3(float3 p)
{
    p = frac(p * 0.3183099 + 0.1);
    p *= 17.0;
    return frac(p.x * p.y * p.z * (p.x + p.y + p.z));
}

// Simplex-like noise
float waterNoise(float2 p)
{
    float2 i = floor(p);
    float2 f = frac(p);

    float2 u = f * f * (3.0 - 2.0 * f);

    float a = hashWater(i + float2(0.0, 0.0));
    float b = hashWater(i + float2(1.0, 0.0));
    float c = hashWater(i + float2(0.0, 1.0));
    float d = hashWater(i + float2(1.0, 1.0));

    return lerp(lerp(a, b, u.x), lerp(c, d, u.x), u.y);
}

// FBM for water surface
float waterFBM(float2 p, int octaves)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = WATER_WAVE_FREQUENCY;

    [unroll]
    for (int i = 0; i < octaves; i++)
    {
        value += amplitude * waterNoise(p * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value;
}

// =============================================================================
// WAVE NORMAL CALCULATION
// =============================================================================

// Calculate wave height at position
float getWaveHeight(float2 pos, float time)
{
#if ENABLE_WATER_EFFECTS
    float3 gerstner = sumGerstnerWaves(pos, time);
    float noise = waterFBM(pos * 0.5 + time * 0.1 * WATER_WAVE_SPEED, WATER_WAVE_OCTAVES) - 0.5;

    return gerstner.y + noise * WATER_WAVE_AMPLITUDE * 0.5;
#else
    return 0.0;
#endif
}

// Calculate wave normal using finite differences
float3 getWaveNormal(float2 pos, float time)
{
#if ENABLE_WATER_EFFECTS
    float epsilon = 0.02;

    float h0 = getWaveHeight(pos, time);
    float hx = getWaveHeight(pos + float2(epsilon, 0.0), time);
    float hz = getWaveHeight(pos + float2(0.0, epsilon), time);

    float3 tangentX = float3(epsilon, hx - h0, 0.0);
    float3 tangentZ = float3(0.0, hz - h0, epsilon);

    return normalize(cross(tangentZ, tangentX));
#else
    return float3(0.0, 1.0, 0.0);
#endif
}

// =============================================================================
// RAINDROP RIPPLES
// =============================================================================

float raindropRipple(float2 pos, float time, float2 dropPos, float dropTime)
{
    float t = time - dropTime;
    if (t < 0.0 || t > 2.0)
        return 0.0;

    float dist = length(pos - dropPos);
    float waveSpeed = 1.5;
    float wavelength = 0.3;
    float decay = exp(-t * 2.0);

    float phase = dist / wavelength - t * waveSpeed / wavelength;
    float envelope = exp(-dist * 0.5) * decay;

    return sin(phase * kTwoPi) * envelope;
}

// Sum multiple raindrops
float sumRaindrops(float2 pos, float time, float intensity)
{
#if ENABLE_RAIN_RIPPLES
    if (intensity <= 0.0)
        return 0.0;

    float ripples = 0.0;
    int numDrops = int(lerp(4.0, 16.0, intensity));

    for (int i = 0; i < numDrops; i++)
    {
        // Pseudo-random drop position and time
        float2 dropPos = float2(
            hashWater(float2(float(i), 0.0)) * 20.0 - 10.0,
            hashWater(float2(float(i), 1.0)) * 20.0 - 10.0);
        dropPos += floor(pos / 20.0) * 20.0;

        float dropPeriod = lerp(0.5, 1.5, hashWater(float2(float(i), 2.0)));
        float dropTime = floor(time / dropPeriod) * dropPeriod;
        dropTime += hashWater(float2(float(i), 3.0)) * dropPeriod;

        ripples += raindropRipple(pos, time, dropPos, dropTime);
    }

    return ripples * RAIN_RIPPLE_INTENSITY * intensity;
#else
    return 0.0;
#endif
}

// =============================================================================
// WATER ABSORPTION AND SCATTERING
// =============================================================================

// Base water absorption coefficients
float3 getWaterAbsorption()
{
    return float3(WATER_ABSORPTION_R, WATER_ABSORPTION_G, WATER_ABSORPTION_B);
}

// CDOM (Colored Dissolved Organic Matter) absorption
// Yellow-brown tint, absorbs blue light
float3 getCDOMAbsorption()
{
#if ENABLE_WATER_PARTICLES
    return float3(0.01, 0.05, 0.3) * WATER_CDOM_AMOUNT;
#else
    return 0.0;
#endif
}

// Chlorophyll absorption (phytoplankton)
// Green tint from algae, absorbs blue and red
float3 getChlorophyllAbsorption()
{
#if ENABLE_WATER_PARTICLES
    return float3(0.2, 0.02, 0.15) * WATER_CHLOROPHYLL_AMOUNT;
#else
    return 0.0;
#endif
}

// Suspended sediment scattering
// Red-brown turbidity
float3 getSedimentScattering()
{
#if ENABLE_WATER_PARTICLES
    return float3(0.3, 0.15, 0.05) * WATER_SEDIMENT_AMOUNT;
#else
    return 0.0;
#endif
}

// Combined water extinction coefficient
float3 getWaterExtinction()
{
    return getWaterAbsorption() + getCDOMAbsorption() + getChlorophyllAbsorption();
}

// Water transmittance over distance
float3 waterTransmittance(float distance)
{
    return exp(-getWaterExtinction() * distance);
}

// In-scattered light in water
float3 waterInscatter(float distance, float3 lightColor, float lightIntensity)
{
    float3 extinction = getWaterExtinction();
    float3 scattering = getSedimentScattering() + float3(0.003, 0.004, 0.005); // Base water scattering

    // Phase function (Henyey-Greenstein with g~0.9 for water)
    float phase = 0.25; // Simplified isotropic for underwater

    // Inscattered light
    float3 transmittance = waterTransmittance(distance);
    float3 inscatter = scattering * lightColor * lightIntensity * phase * (1.0 - transmittance) / (extinction + 0.001);

    return inscatter;
}

// =============================================================================
// FOAM GENERATION
// =============================================================================

float calculateFoam(float2 pos, float time, float depth, float waveHeight)
{
#if ENABLE_WATER_FOAM
    // Shoreline foam (based on depth)
    float depthFoam = smoothstep(WATER_FOAM_THRESHOLD, 0.0, depth);

    // Wave crest foam (based on wave height)
    float crestFoam = smoothstep(WATER_WAVE_AMPLITUDE * 0.5, WATER_WAVE_AMPLITUDE, waveHeight);

    // Noise for foam texture
    float foamNoise = waterFBM(pos * 2.0 + time * 0.5, 4);
    foamNoise = smoothstep(0.3, 0.7, foamNoise);

    float foam = max(depthFoam, crestFoam) * foamNoise * WATER_FOAM_INTENSITY;

    return saturate(foam);
#else
    return 0.0;
#endif
}

// =============================================================================
// FRESNEL AND REFLECTION
// =============================================================================

// Water Fresnel
float waterFresnel(float NdotV, float eta)
{
    // Schlick approximation
    float f0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
    return f0 + (1.0 - f0) * pow(1.0 - NdotV, 5.0);
}

// =============================================================================
// CAUSTICS
// =============================================================================

// Procedural caustics pattern
float causticPattern(float2 pos, float time)
{
#if ENABLE_CAUSTICS
    float2 uv = pos * CAUSTICS_SCALE;

    // Multiple overlapping caustic layers
    float caustics = 0.0;

    // Layer 1
    float2 uv1 = uv + float2(time * 0.1, time * 0.05);
    float n1 = waterNoise(uv1 * 3.0);
    float n2 = waterNoise(uv1 * 3.0 + 0.5);
    caustics += pow(abs(n1 - n2), 0.5);

    // Layer 2 (different speed and direction)
    float2 uv2 = uv + float2(-time * 0.08, time * 0.12);
    float n3 = waterNoise(uv2 * 2.5);
    float n4 = waterNoise(uv2 * 2.5 + 0.3);
    caustics += pow(abs(n3 - n4), 0.5);

    // Normalize and enhance contrast
    caustics *= 0.5;
    caustics = pow(caustics, 2.0);

    return caustics * CAUSTICS_INTENSITY;
#else
    return 0.0;
#endif
}

// Apply caustics to underwater surface
float3 applyCaustics(float3 worldPos, float3 surfaceColor, float3 lightDir, float time, float waterDepth)
{
#if ENABLE_CAUSTICS
    if (waterDepth <= 0.0)
        return surfaceColor;

    // Project caustic pattern
    float2 causticPos = worldPos.xz;

    // Depth-based falloff
    float falloff = exp(-waterDepth * CAUSTICS_FALLOFF);

    // Sample caustics
    float caustics = causticPattern(causticPos, time * CAUSTICS_SPEED);

    // Light direction affects intensity
    float lightFactor = saturate(lightDir.y);

    return surfaceColor * (1.0 + caustics * falloff * lightFactor);
#else
    return surfaceColor;
#endif
}

// =============================================================================
// MAIN WATER SURFACE EVALUATION
// =============================================================================

struct WaterSurface
{
    float3 position;       // World position with waves
    float3 normal;         // Perturbed normal
    float3 baseColor;      // Water color
    float roughness;       // Surface roughness
    float foam;            // Foam amount
    float fresnel;         // Fresnel term
    float3 refractDir;     // Refraction direction
};

WaterSurface evaluateWaterSurface(float3 basePosition, float3 viewDir, float time, float rainIntensity)
{
    WaterSurface water;

    float2 pos2D = basePosition.xz;

    // Calculate waves
    float waveHeight = getWaveHeight(pos2D, time);

    // Add rain ripples
    float ripples = sumRaindrops(pos2D, time, rainIntensity);
    waveHeight += ripples;

    water.position = basePosition + float3(0.0, waveHeight, 0.0);

    // Calculate normal
    float3 waveNormal = getWaveNormal(pos2D, time);

    // Add ripple normal perturbation
    if (abs(ripples) > 0.001)
    {
        float rippleEpsilon = 0.1;
        float rh0 = sumRaindrops(pos2D, time, rainIntensity);
        float rhx = sumRaindrops(pos2D + float2(rippleEpsilon, 0.0), time, rainIntensity);
        float rhz = sumRaindrops(pos2D + float2(0.0, rippleEpsilon), time, rainIntensity);

        float3 rippleNormal = normalize(float3(rh0 - rhx, 1.0, rh0 - rhz));
        waveNormal = normalize(waveNormal + rippleNormal * 0.3);
    }

    water.normal = waveNormal;

    // Water color (slightly blue-green)
    water.baseColor = float3(0.05, 0.15, 0.2);

    // Roughness
#if WATER_ROUGHNESS_OVERRIDE >= 0.0
    water.roughness = WATER_ROUGHNESS_OVERRIDE;
#else
    water.roughness = 0.02 + rainIntensity * 0.05; // Rougher in rain
#endif

    // Calculate foam
    float depth = 0.0; // Would need actual depth buffer
    water.foam = calculateFoam(pos2D, time, depth, waveHeight);

    // Fresnel
    float NdotV = saturate(dot(water.normal, viewDir));
    water.fresnel = waterFresnel(NdotV, 1.0 / WATER_IOR);

    // Refraction direction
    water.refractDir = refract(-viewDir, water.normal, 1.0 / WATER_IOR);

    return water;
}

// =============================================================================
// UNDERWATER FOG
// =============================================================================

float3 applyUnderwaterFog(float3 color, float distance, float3 lightColor, float lightIntensity)
{
    float3 transmittance = waterTransmittance(distance);
    float3 inscatter = waterInscatter(distance, lightColor, lightIntensity);

    return color * transmittance + inscatter;
}

#endif // __OPENRTX_WATER_HLSL__
