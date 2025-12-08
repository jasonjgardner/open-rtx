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
// OpenRTX Volumetric Lighting
// Implements volumetric fog, god rays, and atmospheric effects
// =============================================================================

#ifndef __OPENRTX_VOLUMETRIC_LIGHTING_HLSL__
#define __OPENRTX_VOLUMETRIC_LIGHTING_HLSL__

#include "Settings.hlsl"

// =============================================================================
// FOG DENSITY FUNCTIONS
// =============================================================================

// Exponential height fog density
float fogDensityExponential(float height, float baseHeight)
{
    return FOG_DENSITY * exp(-FOG_HEIGHT_FALLOFF * max(0.0, height - baseHeight));
}

// Uniform fog density
float fogDensityUniform()
{
    return FOG_DENSITY;
}

// Density modifier based on height and time-of-day
float calcDensityModifier(float3 worldPos, float3 sunDir, bool isUnderwater)
{
    float densityModifier = 1.0;

    if (!isUnderwater)
    {
        // Reduced height-based falloff for better light transmission
        // Use gentler exponential curve that doesn't cut off light abruptly
        float heightAboveGround = max(0.0, worldPos.y - 64.0);
        float heightFactor = exp(-FOG_HEIGHT_FALLOFF * 0.3 * heightAboveGround); // Reduced falloff
        densityModifier *= heightFactor;

        // Noon fog reduction (less fog at midday)
        // When NOON_FOG_REDUCTION = 1.0, this has no effect
        float dayFactor = saturate(sunDir.y * 2.0 + 0.5);
        float noonFactor = lerp(1.0, NOON_FOG_REDUCTION, dayFactor * dayFactor * 0.5); // Gentler reduction
        densityModifier *= noonFactor;
    }

    return densityModifier;
}

// Combined fog density with all modifiers
float getFogDensity(float3 worldPos, float rainIntensity, float time, float3 sunDir, bool isUnderwater)
{
    float baseDensity = FOG_DENSITY;

    // Apply density modifiers
    baseDensity *= calcDensityModifier(worldPos, sunDir, isUnderwater);

#if ENABLE_RAIN_FOG
    // Increase density during rain (lerp toward STATIC_RAIN_FOG_AMOUNT)
    baseDensity = lerp(baseDensity, baseDensity * (1.0 + RAIN_FOG_AMOUNT), rainIntensity);
#endif

    return baseDensity;
}

// Backwards compatibility - without sun direction
float getFogDensity(float3 worldPos, float rainIntensity, float time)
{
    return getFogDensity(worldPos, rainIntensity, time, float3(0, 1, 0), false);
}

// =============================================================================
// ADVANCED PHASE FUNCTIONS
// =============================================================================

// Henyey-Greenstein phase function with proper normalization
float henyeyGreenstein(float cosTheta, float g)
{
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return (1.0 - g2) / (4.0 * kPi * pow(abs(denom), 1.5));
}

// Cornette-Shanks phase function (improved Mie scattering)
float cornetteShanks(float cosTheta, float g)
{
    float g2 = g * g;
    float num = 3.0 * (1.0 - g2) * (1.0 + cosTheta * cosTheta);
    float denom = (8.0 * kPi) * (2.0 + g2) * pow(abs(1.0 + g2 - 2.0 * g * cosTheta), 1.5);
    return num / max(denom, 1e-7);
}

// Dual-lobe phase function for realistic scattering
float dualLobePhase(float cosTheta, float g1, float g2, float blend)
{
    return lerp(henyeyGreenstein(cosTheta, g2), henyeyGreenstein(cosTheta, g1), blend);
}

// Air fog phase function using AIR_FOG_ASYMMETRY setting
float airFogPhase(float cosTheta)
{
    return henyeyGreenstein(cosTheta, AIR_FOG_ASYMMETRY);
}

// Water fog phase function using WATER_FOG_ASYMMETRY setting  
float waterFogPhase(float cosTheta)
{
    return henyeyGreenstein(cosTheta, WATER_FOG_ASYMMETRY);
}

// Isotropic phase function for GI inscatter (uniform in all directions)
float isotropicPhase()
{
    return 1.0 / (4.0 * kPi);
}

// =============================================================================
// RAINBOW PHYSICS
// =============================================================================
// Rainbows form when sunlight refracts through water droplets.
// Primary rainbow: ~42° from antisolar point (one internal reflection)
// Secondary rainbow: ~51° from antisolar point (two internal reflections)
// Alexander's dark band: 42°-51° region between rainbows (darker sky)

// Physical constants for rainbow angles (in radians)
static const float kPrimaryRainbowAngle = 0.7330;    // 42.0° - center of primary bow
static const float kSecondaryRainbowAngle = 0.8901;  // 51.0° - center of secondary bow

// Angular dispersion causes color separation (red on outside, violet inside for primary)
// These are the deviation angles for each wavelength
static const float kPrimaryRedAngle = 0.7383;     // 42.3° - red (outer edge)
static const float kPrimaryGreenAngle = 0.7243;   // 41.5° - green (middle)
static const float kPrimaryBlueAngle = 0.7086;    // 40.6° - blue/violet (inner edge)

// Secondary rainbow has reversed color order (red inside, violet outside)
static const float kSecondaryRedAngle = 0.8814;   // 50.5° - red (inner edge)
static const float kSecondaryGreenAngle = 0.8901; // 51.0° - green (middle)
static const float kSecondaryBlueAngle = 0.9076;  // 52.0° - blue/violet (outer edge)

// Gaussian intensity profile for rainbow band
// Models the angular distribution of scattered light
float rainbowGaussian(float angle, float centerAngle, float width)
{
    float delta = angle - centerAngle;
    return exp(-delta * delta / (2.0 * width * width));
}

// Supernumerary fringes - interference patterns inside primary rainbow
// Creates subtle colored bands due to wave interference
float supernumeraryPattern(float angle, float baseAngle)
{
    // Supernumeraries appear just inside the primary bow
    float delta = baseAngle - angle;
    if (delta < 0.0 || delta > 0.1)
        return 0.0;

    // Interference creates oscillating intensity
    // Spacing depends on droplet size - we use average spacing
    float fringeSpacing = 0.015; // ~0.9° spacing
    float phase = delta / fringeSpacing * kTwoPi;

    // Damped oscillation - fringes fade toward center
    float envelope = exp(-delta * 8.0);
    return max(0.0, cos(phase) * envelope * 0.3);
}

// Calculate complete rainbow contribution
// Returns additive color to be added to scene (0 = no rainbow visible)
float3 computeRainbow(float3 viewDir, float3 lightDir, float distance)
{
#if !ENABLE_RAINBOW
    return 0.0;
#endif

    // Rainbow appears at antisolar point - opposite the sun
    // Calculate angle from antisolar direction
    float cosAngle = -dot(normalize(viewDir), normalize(lightDir));

    // Clamp to valid range for acos
    cosAngle = clamp(cosAngle, -0.9999, 0.9999);
    float angle = acos(cosAngle);

    // Rainbow only visible when looking away from sun (antisolar hemisphere)
    if (cosAngle > 0.2)
        return 0.0;

    // Angular width of rainbow bands (affected by droplet size variation)
    float bandWidth = 0.025 * RAINBOW_WIDTH;

    float3 rainbowColor = 0.0;

    // ===================
    // PRIMARY RAINBOW (brighter, ~42°)
    // ===================
    float3 primaryIntensity;
    primaryIntensity.r = rainbowGaussian(angle, kPrimaryRedAngle, bandWidth);
    primaryIntensity.g = rainbowGaussian(angle, kPrimaryGreenAngle, bandWidth * 0.9);
    primaryIntensity.b = rainbowGaussian(angle, kPrimaryBlueAngle, bandWidth * 0.8);

    // Add supernumerary fringes (subtle interference patterns)
#if ENABLE_SUPERNUMERARY_FRINGES
    float fringeR = supernumeraryPattern(angle, kPrimaryRedAngle);
    float fringeG = supernumeraryPattern(angle, kPrimaryGreenAngle);
    float fringeB = supernumeraryPattern(angle, kPrimaryBlueAngle);
    primaryIntensity += float3(fringeR, fringeG, fringeB);
#endif

    rainbowColor += primaryIntensity * RAINBOW_PRIMARY_INTENSITY;

    // ===================
    // SECONDARY RAINBOW (dimmer, ~51°, reversed colors)
    // ===================
#if ENABLE_SECONDARY_RAINBOW
    float3 secondaryIntensity;
    // Note: color order is reversed for secondary rainbow
    secondaryIntensity.r = rainbowGaussian(angle, kSecondaryRedAngle, bandWidth * 1.3);
    secondaryIntensity.g = rainbowGaussian(angle, kSecondaryGreenAngle, bandWidth * 1.2);
    secondaryIntensity.b = rainbowGaussian(angle, kSecondaryBlueAngle, bandWidth * 1.1);

    // Secondary is about 43% as bright as primary (due to extra reflection)
    rainbowColor += secondaryIntensity * RAINBOW_SECONDARY_INTENSITY;
#endif

    // ===================
    // ALEXANDER'S DARK BAND
    // ===================
    // Region between primary and secondary rainbows appears darker
    // because no light is scattered into this angular range
#if ENABLE_ALEXANDERS_BAND
    float bandCenter = (kPrimaryRainbowAngle + kSecondaryRainbowAngle) * 0.5;
    float bandHalfWidth = (kSecondaryRainbowAngle - kPrimaryRainbowAngle) * 0.4;
    float darkBand = smoothstep(bandHalfWidth, 0.0, abs(angle - bandCenter));
    // This creates a subtle darkening effect (applied elsewhere in sky rendering)
#endif

    // ===================
    // DISTANCE AND VISIBILITY FACTORS
    // ===================

    // Rainbow needs sufficient distance to be visible (forms in distant rain)
    float distanceFactor = smoothstep(32.0, 128.0, distance);

    // Rainbow visibility depends on viewing geometry
    // Strongest when looking toward horizon (horizontal view)
    float horizonFactor = 1.0 - abs(viewDir.y) * 0.5;

    // Combine factors
    rainbowColor *= distanceFactor * horizonFactor * RAINBOW_INTENSITY;

    // Apply subtle color saturation boost for more vivid appearance
    float luminance = dot(rainbowColor, float3(0.299, 0.587, 0.114));
    rainbowColor = lerp(float3(luminance, luminance, luminance), rainbowColor, RAINBOW_SATURATION);

    return max(rainbowColor, 0.0);
}

// =============================================================================
// VOLUMETRIC RAY MARCHING
// =============================================================================

struct VolumetricOutput
{
    float3 inscatter;     // In-scattered light
    float transmittance;  // Remaining visibility
};

// Enhanced ray marching with improved energy conservation
VolumetricOutput marchVolumetricFog(
    float3 rayOrigin,
    float3 rayDir,
    float maxDistance,
    float3 sunDir,
    float3 sunColor,
    float rainIntensity,
    float time,
    int numSteps)
{
    VolumetricOutput output;
    output.inscatter = 0.0;
    output.transmittance = 1.0;

#if !ENABLE_VOLUMETRIC_LIGHTING
    return output;
#endif

    // Early out if beyond max fog distance
    float effectiveMaxDistance = min(maxDistance, FOG_MAX_DISTANCE);
    if (effectiveMaxDistance <= 0.01)
        return output;

    // Adaptive step sizing: fewer steps for short distances
    int adaptiveSteps = numSteps;
    if (effectiveMaxDistance < 32.0)
    {
        adaptiveSteps = max(4, numSteps / 4);
    }
    else if (effectiveMaxDistance < 64.0)
    {
        adaptiveSteps = max(8, numSteps / 2);
    }

    float stepSize = effectiveMaxDistance / float(adaptiveSteps);
    float3 currentPos = rayOrigin;

    // Calculate phase functions
    float cosTheta = dot(rayDir, sunDir);
    float phase = dualLobePhase(cosTheta, 0.8, -0.3, 0.7);

    // Use proper media properties with reduced coefficients to prevent washout
    // Scale scattering by global density to make it configurable
    float3 scattering = FOG_SCATTERING_COEFFICIENTS * FOG_DENSITY * 10.0;
    float3 absorption = float3(0.0002, 0.0003, 0.0005); // Slight wavelength-dependent absorption
    float3 extinction = scattering + absorption;

    // Track cumulative transmittance for proper inscatter weighting
    float3 cumulativeTransmittance = 1.0;

    for (int i = 0; i < adaptiveSteps; i++)
    {
        // Get density modifier (height falloff, time of day, etc.)
        float densityModifier = calcDensityModifier(currentPos, sunDir, false);

        // Apply rain fog boost
        float effectiveDensity = densityModifier;
#if ENABLE_RAIN_FOG
        effectiveDensity = lerp(effectiveDensity, effectiveDensity * (1.0 + RAIN_FOG_AMOUNT), rainIntensity);
#endif

        if (effectiveDensity > 0.0001)
        {
            // Beer-Lambert law for this step
            float3 stepExtinction = extinction * effectiveDensity * stepSize;
            float3 stepTransmittance = exp(-stepExtinction);

            // In-scattered light calculation
            float3 lightContribution = 0.0;

#if ENABLE_SUN_FOG
            // Sun contribution with terrain shadow approximation
            float3 sunTransmission = 1.0;
            float terrainHeight = 64.0;
            if (sunDir.y < 0.3 && currentPos.y > terrainHeight && rayDir.y < -0.1)
            {
                float shadowDistance = (currentPos.y - terrainHeight) / max(-rayDir.y, 0.01);
                sunTransmission = saturate(exp(-shadowDistance * 0.01));
            }

            // Rainbow effect (only when raining and looking toward antisolar point)
            float3 rainbowFactor = computeRainbow(rayDir, sunDir, effectiveMaxDistance);
            lightContribution += sunColor * phase * sunTransmission * rainbowFactor * SUN_FOG_AMOUNT;
#endif

#if ENABLE_STATIC_GI_FOG
            // Ambient/GI contribution - reduced to prevent washout
            float3 ambientColor = float3(0.08, 0.10, 0.14);
            lightContribution += ambientColor * STATIC_GI_FOG_AMOUNT;
#endif

#if ENABLE_EXPLICIT_LIGHT_SAMPLING
            // Emissive contribution (very subtle)
            float3 emissiveContribution = float3(0.005, 0.003, 0.002) * effectiveDensity;
            lightContribution += emissiveContribution;
#endif

            // Proper inscatter integration using cumulative transmittance
            // This ensures light scattered at position i is attenuated by fog between camera and i
            float3 stepInscatter = lightContribution * scattering * effectiveDensity * stepSize;

            // Weight by cumulative transmittance (light must travel through fog to reach camera)
            output.inscatter += stepInscatter * cumulativeTransmittance;

            // Update cumulative transmittance for next step
            cumulativeTransmittance *= stepTransmittance;
        }

        // Always advance ray
        currentPos += rayDir * stepSize;

        // Early termination when fog is essentially opaque
        float avgTransmittance = dot(cumulativeTransmittance, float3(0.333, 0.333, 0.334));
        if (avgTransmittance < 0.01)
            break;
    }

    // Store final transmittance (average of RGB channels)
    output.transmittance = dot(cumulativeTransmittance, float3(0.333, 0.333, 0.334));

    return output;
}

// =============================================================================
// GOD RAYS / LIGHT SHAFTS
// =============================================================================

// Screen-space god rays (post-process approach)
float3 computeGodRays(
    float2 uv,
    float2 sunScreenPos,
    float3 sceneColor,
    float depth,
    int numSamples)
{
    if (!all(saturate(sunScreenPos) == sunScreenPos))
        return 0.0; // Sun not on screen

    float3 godRays = 0.0;
    float2 deltaUV = (sunScreenPos - uv) / float(numSamples);

    float decay = 0.96;
    float weight = 0.5;
    float illuminationDecay = 1.0;

    float2 sampleUV = uv;

    for (int i = 0; i < numSamples; i++)
    {
        sampleUV += deltaUV;

        // Sample scene brightness (would need actual texture)
        // float3 sampleColor = sceneTexture.Sample(sampleUV);
        float3 sampleColor = sceneColor * (1.0 - float(i) / float(numSamples));

        sampleColor *= illuminationDecay * weight;
        godRays += sampleColor;
        illuminationDecay *= decay;
    }

    return godRays * 0.1;
}

// =============================================================================
// RAINBOW SPECTRAL COLORS
// =============================================================================

// Convert wavelength (380-780nm) to approximate RGB color
// Based on CIE color matching functions approximation
float3 wavelengthToRGB(float wavelength)
{
    float3 color = 0.0;

    if (wavelength >= 380.0 && wavelength < 440.0)
    {
        color.r = -(wavelength - 440.0) / (440.0 - 380.0);
        color.g = 0.0;
        color.b = 1.0;
    }
    else if (wavelength >= 440.0 && wavelength < 490.0)
    {
        color.r = 0.0;
        color.g = (wavelength - 440.0) / (490.0 - 440.0);
        color.b = 1.0;
    }
    else if (wavelength >= 490.0 && wavelength < 510.0)
    {
        color.r = 0.0;
        color.g = 1.0;
        color.b = -(wavelength - 510.0) / (510.0 - 490.0);
    }
    else if (wavelength >= 510.0 && wavelength < 580.0)
    {
        color.r = (wavelength - 510.0) / (580.0 - 510.0);
        color.g = 1.0;
        color.b = 0.0;
    }
    else if (wavelength >= 580.0 && wavelength < 645.0)
    {
        color.r = 1.0;
        color.g = -(wavelength - 645.0) / (645.0 - 580.0);
        color.b = 0.0;
    }
    else if (wavelength >= 645.0 && wavelength <= 780.0)
    {
        color.r = 1.0;
        color.g = 0.0;
        color.b = 0.0;
    }

    // Intensity falloff at spectrum edges
    float intensity = 1.0;
    if (wavelength >= 380.0 && wavelength < 420.0)
        intensity = 0.3 + 0.7 * (wavelength - 380.0) / (420.0 - 380.0);
    else if (wavelength > 700.0 && wavelength <= 780.0)
        intensity = 0.3 + 0.7 * (780.0 - wavelength) / (780.0 - 700.0);

    return color * intensity;
}

// Sample spectral rainbow color at normalized position (0=violet, 1=red)
float3 sampleRainbowSpectrum(float t)
{
    // Map t to wavelength range (violet 400nm to red 700nm)
    float wavelength = lerp(400.0, 700.0, saturate(t));
    return wavelengthToRGB(wavelength);
}

// Rainbow color based on angle from rainbow center
// Used for artistic/legacy rendering
float3 rainbowColor(float angle)
{
    // Primary rainbow angular range
    const float innerAngle = 0.708;   // ~40.6° (violet/blue edge)
    const float outerAngle = 0.738;   // ~42.3° (red edge)
    const float rainbowWidth = outerAngle - innerAngle;

    // Calculate position within rainbow band
    float t = (angle - innerAngle) / rainbowWidth;

    if (t < -0.2 || t > 1.2)
        return 0.0;

    // Smooth edges
    float edgeFade = smoothstep(-0.2, 0.0, t) * smoothstep(1.2, 1.0, t);

    // Sample spectrum
    float3 color = sampleRainbowSpectrum(t);

    return color * edgeFade;
}

// =============================================================================
// FLASHLIGHT VOLUMETRIC
// =============================================================================

#if ENABLE_FLASHLIGHT

float3 computeFlashlightVolumetric(
    float3 rayOrigin,
    float3 rayDir,
    float maxDistance,
    float3 flashlightPos,
    float3 flashlightDir,
    float rainIntensity,
    float time,
    int numSteps)
{
#if !FLASHLIGHT_VOLUMETRIC
    return 0.0;
#endif

    float3 volumetric = 0.0;
    float stepSize = min(maxDistance, FLASHLIGHT_RANGE) / float(numSteps);
    float3 currentPos = rayOrigin;

    for (int i = 0; i < numSteps; i++)
    {
        float density = getFogDensity(currentPos, rainIntensity, time);

        if (density > 0.0001)
        {
            // Direction to flashlight
            float3 toFlashlight = currentPos - flashlightPos;
            float dist = length(toFlashlight);
            toFlashlight /= max(dist, 0.001);

            // Cone attenuation
            float cosAngle = dot(toFlashlight, flashlightDir);
            float innerCone = cos(FLASHLIGHT_INNER_CONE);
            float outerCone = cos(FLASHLIGHT_OUTER_CONE);
            float coneAtten = smoothstep(outerCone, innerCone, cosAngle);

            // Distance attenuation
            float distAtten = 1.0 / (1.0 + dist * dist * 0.01);

            // Range falloff
            float rangeFalloff = 1.0 - saturate(dist / FLASHLIGHT_RANGE);

            // Phase function
            float cosViewLight = dot(rayDir, -toFlashlight);
            float phase = henyeyGreenstein(cosViewLight, 0.5);

            // Accumulate
            float3 lightColor = FLASHLIGHT_COLOR * FLASHLIGHT_INTENSITY;
            volumetric += lightColor * density * coneAtten * distAtten * rangeFalloff * phase * stepSize;
        }

        currentPos += rayDir * stepSize;
    }

    return volumetric;
}

#endif // ENABLE_FLASHLIGHT

// =============================================================================
// DIMENSION-SPECIFIC VOLUMETRICS
// =============================================================================

// Nether fog
float3 computeNetherFog(float3 rayOrigin, float3 rayDir, float distance)
{
    float density = NETHER_FOG_AMOUNT * FOG_DENSITY * 2.0;
    float extinction = density * min(distance, FOG_MAX_DISTANCE * 0.5);

    float3 fogColor = float3(0.3, 0.1, 0.05); // Orange-red nether fog
    float transmittance = exp(-extinction);

    return fogColor * (1.0 - transmittance);
}

// End fog
float3 computeEndFog(float3 rayOrigin, float3 rayDir, float distance)
{
    float density = END_FOG_DENSITY;
    float extinction = density * min(distance, FOG_MAX_DISTANCE);

    float transmittance = exp(-extinction);

    return END_FOG_COLOR * (1.0 - transmittance);
}

// =============================================================================
// MAIN VOLUMETRIC FUNCTION
// =============================================================================

struct VolumetricResult
{
    float3 inscatter;
    float3 extinction;
    float transmittance;
};

VolumetricResult evaluateVolumetrics(
    float3 rayOrigin,
    float3 rayDir,
    float maxDistance,
    float3 sunDir,
    float3 sunColor,
    float rainIntensity,
    float time,
    int dimension) // 0 = overworld, 1 = nether, 2 = end
{
    VolumetricResult result;
    result.inscatter = 0.0;
    result.extinction = 0.0;
    result.transmittance = 1.0;

#if !ENABLE_VOLUMETRIC_LIGHTING
    return result;
#endif

    // Dimension-specific handling
    if (dimension == 1) // Nether
    {
        result.inscatter = computeNetherFog(rayOrigin, rayDir, maxDistance);
        result.transmittance = exp(-NETHER_FOG_AMOUNT * min(maxDistance, FOG_MAX_DISTANCE * 0.5));
        return result;
    }
    else if (dimension == 2) // End
    {
        result.inscatter = computeEndFog(rayOrigin, rayDir, maxDistance);
        result.transmittance = exp(-END_FOG_DENSITY * min(maxDistance, FOG_MAX_DISTANCE));
        return result;
    }

    // Overworld volumetrics
    VolumetricOutput fogOutput = marchVolumetricFog(
        rayOrigin, rayDir, maxDistance,
        sunDir, sunColor,
        rainIntensity, time,
        VOLUMETRIC_STEPS);

    result.inscatter = fogOutput.inscatter;
    result.transmittance = fogOutput.transmittance;

    // Rainbow effect - visible during/after rain when looking toward antisolar point
    // computeRainbow returns additive color (0 when no rainbow visible)
#if ENABLE_RAINBOW
    if (rainIntensity > 0.05)
    {
        // Get rainbow contribution (already handles angle, distance, etc.)
        float3 rainbow = computeRainbow(rayDir, sunDir, maxDistance);

        // Scale by rain intensity - rainbow fades as rain stops
        // Rainbow is most visible just after rain (rainIntensity 0.1-0.5)
        float rainbowVisibility = smoothstep(0.05, 0.2, rainIntensity) * smoothstep(1.0, 0.3, rainIntensity);

        // Add to inscatter, weighted by transmittance and sun brightness
        float sunBrightness = saturate(sunDir.y + 0.1); // Rainbow needs sunlight
        result.inscatter += rainbow * rainbowVisibility * sunBrightness * result.transmittance;
    }
#endif

    return result;
}

// Apply volumetrics to final color with physically-based blending
float3 applyVolumetrics(float3 sceneColor, VolumetricResult volumetrics)
{
    // Physically correct fog blending:
    // finalColor = sceneColor * transmittance + inscatter
    //
    // Where:
    // - transmittance = how much of the original scene is visible through fog
    // - inscatter = light added by the fog itself (sun/ambient scattered toward camera)

    // Apply scene attenuation through fog
    float3 transmittedScene = sceneColor * volumetrics.transmittance;

    // Add in-scattered light from fog
    // The inscatter is already properly weighted by transmittance during marching
    float3 finalColor = transmittedScene + volumetrics.inscatter;

    // Clamp to prevent negative values or excessive brightness
    // Allow inscatter to brighten beyond scene color (for god rays, etc.)
    // but cap at reasonable maximum to prevent HDR blowout
    return max(finalColor, 0.0);
}

#endif // __OPENRTX_VOLUMETRIC_LIGHTING_HLSL__
