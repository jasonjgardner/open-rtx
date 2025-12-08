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

// Domain warped FBM for organic foam patterns
float foamFBM(float2 p, float time)
{
    // Domain warping for more organic shapes
    float2 q = float2(
        waterFBM(p + float2(0.0, 0.0), 3),
        waterFBM(p + float2(5.2, 1.3), 3)
    );

    float2 r = float2(
        waterFBM(p + 4.0 * q + float2(1.7, 9.2) + time * 0.15, 3),
        waterFBM(p + 4.0 * q + float2(8.3, 2.8) + time * 0.12, 3)
    );

    return waterFBM(p + 4.0 * r, 4);
}

// Voronoi-based foam bubbles
float foamBubbles(float2 p, float time)
{
    float2 i = floor(p);
    float2 f = frac(p);

    float minDist = 1.0;

    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 neighbor = float2(x, y);
            float2 cellCenter = hashWater(i + neighbor) * 0.5 + 0.25;

            // Animate bubble positions slightly
            cellCenter += 0.1 * sin(time * 2.0 + 6.2831 * hashWater(i + neighbor + 100.0));

            float2 diff = neighbor + cellCenter - f;
            float dist = length(diff);
            minDist = min(minDist, dist);
        }
    }

    return 1.0 - smoothstep(0.0, 0.4, minDist);
}

float calculateFoam(float2 pos, float time, float depth, float waveHeight)
{
#if ENABLE_WATER_FOAM
    // Shoreline foam (based on depth) with persistence
    float depthFoam = smoothstep(WATER_FOAM_THRESHOLD, 0.0, depth);
    depthFoam = pow(depthFoam, 0.7);  // Softer falloff

    // Wave crest foam with Jacobian-like fold detection
    float crestFoam = smoothstep(WATER_WAVE_AMPLITUDE * 0.4, WATER_WAVE_AMPLITUDE * 0.9, waveHeight);

    // Calculate wave steepness for breaking wave foam
    float epsilon = 0.05;
    float hLeft = getWaveHeight(pos + float2(-epsilon, 0.0), time);
    float hRight = getWaveHeight(pos + float2(epsilon, 0.0), time);
    float slope = abs(hRight - hLeft) / (2.0 * epsilon);
    float breakingFoam = smoothstep(0.3, 0.8, slope);

    // Domain-warped FBM for organic foam texture
    float foamPattern = foamFBM(pos * 1.5, time);
    foamPattern = smoothstep(0.35, 0.65, foamPattern);

    // Add bubble detail at shorelines
    float bubbles = foamBubbles(pos * 8.0, time) * depthFoam;

    // Combine foam types
    float baseFoam = max(max(depthFoam, crestFoam), breakingFoam * 0.7);
    float foam = baseFoam * foamPattern * WATER_FOAM_INTENSITY;
    foam += bubbles * 0.3 * WATER_FOAM_INTENSITY;

    // Foam lifetime - older foam fades
    float foamAge = frac(time * 0.1 + hashWater(floor(pos)));
    foam *= lerp(0.5, 1.0, 1.0 - foamAge * depthFoam);

    return saturate(foam);
#else
    return 0.0;
#endif
}

// =============================================================================
// FRESNEL AND REFLECTION
// =============================================================================

// Water Fresnel using accurate dielectric Fresnel equations
// For air-to-water: ior = WATER_IOR (1.333)
// Returns reflectance at given angle
float waterFresnel(float NdotV, float ior)
{
    // Use full Fresnel equations for accuracy, especially at grazing angles
    float cosTheta = NdotV;
    float eta = 1.0 / ior;  // Air to water ratio

    float sinThetaTSq = eta * eta * (1.0 - cosTheta * cosTheta);

    // Total internal reflection check (shouldn't happen for air-to-water looking down)
    if (sinThetaTSq > 1.0)
        return 1.0;

    float cosThetaT = sqrt(1.0 - sinThetaTSq);

    // Fresnel equations for s and p polarization
    float rs = (eta * cosTheta - cosThetaT) / (eta * cosTheta + cosThetaT);
    float rp = (cosTheta - eta * cosThetaT) / (cosTheta + eta * cosThetaT);

    // Average of s and p polarization (unpolarized light)
    return (rs * rs + rp * rp) * 0.5;
}

// Schlick approximation for water Fresnel (faster, less accurate at grazing angles)
float waterFresnelSchlick(float NdotV, float ior)
{
    float f0 = pow((ior - 1.0) / (ior + 1.0), 2.0);
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

// Get raw caustics value for a world position (used for underwater surface lighting)
float3 waterCaustics(float3 worldPos, float time, float3 lightDir)
{
#if ENABLE_CAUSTICS
    // Project onto XZ plane for caustic pattern
    float2 causticPos = worldPos.xz;

    // Sample caustics pattern
    float caustics = causticPattern(causticPos, time * CAUSTICS_SPEED);

    // Light direction affects intensity (caustics are stronger with overhead sun)
    float lightFactor = saturate(lightDir.y);

    // Return as light color contribution (white/cyan tinted)
    return float3(0.9, 1.0, 1.0) * caustics * lightFactor;
#else
    return 0.0;
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

    // Roughness - use override if defined, otherwise dynamic based on rain
#ifdef USE_WATER_ROUGHNESS_OVERRIDE
    water.roughness = WATER_ROUGHNESS_OVERRIDE;
#else
    water.roughness = 0.02 + rainIntensity * 0.05; // Rougher in rain
#endif

    // Calculate foam
    float depth = 0.0; // Would need actual depth buffer
    water.foam = calculateFoam(pos2D, time, depth, waveHeight);

    // Fresnel - use accurate dielectric Fresnel with water IOR
    float NdotV = saturate(dot(water.normal, viewDir));
    water.fresnel = waterFresnel(NdotV, WATER_IOR);

    // Refraction direction (HLSL refract expects eta = n1/n2 for air-to-water)
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

// =============================================================================
// UNDERWATER CAMERA DISTORTION
// =============================================================================

// Simple 2D noise for underwater wave distortion
float underwaterNoise(float2 p)
{
    float2 i = floor(p);
    float2 f = frac(p);
    f = f * f * (3.0 - 2.0 * f);  // Smoothstep

    float a = frac(sin(dot(i, float2(12.9898, 78.233))) * 43758.5453);
    float b = frac(sin(dot(i + float2(1, 0), float2(12.9898, 78.233))) * 43758.5453);
    float c = frac(sin(dot(i + float2(0, 1), float2(12.9898, 78.233))) * 43758.5453);
    float d = frac(sin(dot(i + float2(1, 1), float2(12.9898, 78.233))) * 43758.5453);

    return lerp(lerp(a, b, f.x), lerp(c, d, f.x), f.y);
}

// Multi-octave noise for more natural wave pattern
float underwaterFBM(float2 p, float time)
{
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    // Animate UV coordinates
    float2 animatedP = p;

    for (int i = 0; i < 3; i++)
    {
        // Different animation direction per octave
        float2 offset = float2(
            sin(time * 0.7 + i * 1.3) * 0.5,
            cos(time * 0.5 + i * 0.9) * 0.5
        );
        value += amplitude * underwaterNoise((animatedP + offset) * frequency);
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    return value;
}

// Calculate underwater ray direction distortion
// This simulates looking up through the water surface from below
// rayDir: original ray direction (normalized)
// uv: screen UV coordinates (0-1)
// time: animation time
// Returns: distorted ray direction
float3 applyUnderwaterRayDistortion(float3 rayDir, float2 uv, float time)
{
#if !ENABLE_UNDERWATER_DISTORTION
    return rayDir;
#endif

    // IOR-based distortion: rays bend when looking up through water surface
    // Water-to-air refraction (eta = WATER_IOR / 1.0)
    float eta = WATER_IOR;  // Water to air

    // Calculate angle from vertical (how much we're looking up vs sideways)
    float verticalComponent = rayDir.y;

    // Only apply strong IOR distortion when looking upward (toward water surface)
    // Looking down or sideways underwater has minimal distortion
    float upwardFactor = saturate(verticalComponent);

    // Snell's law distortion - rays bend away from normal when going to less dense medium
    // This creates the "dome" effect when looking up from underwater
    float sinTheta = sqrt(1.0 - verticalComponent * verticalComponent);
    float sinThetaRefracted = sinTheta * eta;

    // Check for total internal reflection (critical angle ~48.6° for water)
    bool totalInternalReflection = sinThetaRefracted > 1.0;

    float3 distortedDir = rayDir;

    if (!totalInternalReflection && upwardFactor > 0.01)
    {
        // Calculate refracted angle
        float cosThetaRefracted = sqrt(1.0 - sinThetaRefracted * sinThetaRefracted);

        // Blend between original and refracted based on upward factor and strength setting
        float distortionAmount = upwardFactor * UNDERWATER_DISTORTION_STRENGTH;

        // Create horizontal distortion that increases toward edges (barrel distortion effect)
        float2 centerOffset = uv - 0.5;
        float edgeDist = length(centerOffset);

        // Radial distortion - rays at edges bend more (fisheye-like effect underwater)
        float radialDistortion = edgeDist * distortionAmount * 0.5;

        // Apply radial push outward for underwater dome effect
        float2 radialOffset = normalize(centerOffset + 0.0001) * radialDistortion;

        distortedDir.xz += radialOffset * (1.0 - abs(verticalComponent));
    }

    // Animated wave distortion (always active when underwater)
    float2 waveUV = uv * UNDERWATER_WAVE_SCALE;
    float waveTime = time * UNDERWATER_WAVE_SPEED;

    // Two overlapping wave patterns for more natural look
    float wave1 = underwaterFBM(waveUV, waveTime);
    float wave2 = underwaterFBM(waveUV * 1.3 + 10.0, waveTime * 0.8);
    float combinedWave = (wave1 + wave2) * 0.5 - 0.5;  // Center around 0

    // Apply wave distortion to ray direction
    float waveStrength = UNDERWATER_WAVE_DISTORTION;
    distortedDir.x += combinedWave * waveStrength;
    distortedDir.z += underwaterFBM(waveUV.yx + 5.0, waveTime * 1.1) * waveStrength - waveStrength * 0.5;

    return normalize(distortedDir);
}

// Calculate UV offset for underwater screen-space distortion
// This is for post-process style distortion on the final image
float2 getUnderwaterUVDistortion(float2 uv, float time)
{
#if !ENABLE_UNDERWATER_DISTORTION
    return float2(0, 0);
#endif

    float2 waveUV = uv * UNDERWATER_WAVE_SCALE;
    float waveTime = time * UNDERWATER_WAVE_SPEED;

    // Layered wave distortion
    float2 distortion;
    distortion.x = underwaterFBM(waveUV, waveTime) - 0.5;
    distortion.y = underwaterFBM(waveUV + float2(7.3, 3.7), waveTime * 0.9) - 0.5;

    // Add slower, larger waves
    float2 largeWaveUV = uv * UNDERWATER_WAVE_SCALE * 0.3;
    distortion.x += (underwaterFBM(largeWaveUV, waveTime * 0.5) - 0.5) * 2.0;
    distortion.y += (underwaterFBM(largeWaveUV + float2(13.1, 17.3), waveTime * 0.4) - 0.5) * 2.0;

    return distortion * UNDERWATER_WAVE_DISTORTION;
}

// Get chromatic aberration offsets for underwater RGB separation
float3 getUnderwaterChromaticOffset(float2 uv)
{
    // Early out if chromatic aberration is disabled (check at runtime since it's a float)
    if (UNDERWATER_CHROMATIC_ABERRATION <= 0.0)
        return float3(0, 0, 0);

    // Distance from center affects aberration strength
    float2 centerOffset = uv - 0.5;
    float edgeDist = length(centerOffset);

    // R shifts outward, B shifts inward (like real chromatic aberration)
    return float3(
        edgeDist * UNDERWATER_CHROMATIC_ABERRATION,   // Red channel offset
        0.0,                                           // Green stays centered
        -edgeDist * UNDERWATER_CHROMATIC_ABERRATION   // Blue channel offset
    );
}

#endif // __OPENRTX_WATER_HLSL__
