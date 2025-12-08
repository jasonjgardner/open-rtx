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
// Following the froxel-based approach from BetterRTX
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
// PHASE FUNCTIONS
// =============================================================================

// Henyey-Greenstein phase function
float henyeyGreenstein(float cosTheta, float g)
{
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return (1.0 - g2) / (4.0 * kPi * pow(abs(denom), 1.5));
}

// Cornette-Shanks phase function (improved Mie)
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
// VOLUMETRIC RAY MARCHING
// =============================================================================

struct VolumetricOutput
{
    float3 inscatter;     // In-scattered light
    float transmittance;  // Remaining visibility
};

// Simple ray marching for volumetric fog
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

    // Phase function for sun
    float cosTheta = dot(rayDir, sunDir);
    float phase = dualLobePhase(cosTheta, 0.8, -0.3, 0.7);

    // Use configured scattering coefficients instead of hardcoded values
    float3 scattering = FOG_SCATTERING_COEFFICIENTS;

    for (int i = 0; i < adaptiveSteps; i++)
    {
        // Use full fog density function with proper sun direction
        float density = getFogDensity(currentPos, rainIntensity, time, sunDir, false);

        if (density > 0.0001)
        {
            // Extinction
            float extinction = density * stepSize;
            float stepTransmittance = exp(-extinction);

            // In-scattered light
            float3 lightContribution = 0.0;

#if ENABLE_SUN_FOG
            // Sun contribution with simple shadow approximation
            // Use height-based shadow for terrain blocking
            float terrainHeight = 64.0; // Approximate terrain height
            float shadowFactor = 1.0;
            
            // Simple terrain shadow: if sun is low and ray goes toward ground, reduce visibility
            if (sunDir.y < 0.3 && currentPos.y > terrainHeight && rayDir.y < -0.1)
            {
                float shadowDistance = (currentPos.y - terrainHeight) / max(-rayDir.y, 0.01);
                shadowFactor = saturate(exp(-shadowDistance * 0.01));
            }
            
            lightContribution += sunColor * phase * shadowFactor * SUN_FOG_AMOUNT;
#endif

#if ENABLE_STATIC_GI_FOG
            // Enhanced ambient/GI contribution with wavelength dependency
            float3 ambientColor = float3(0.15, 0.18, 0.22); // Brighter sky ambient
            lightContribution += ambientColor * STATIC_GI_FOG_AMOUNT;
#endif

            // Improved light accumulation with proper energy conservation
            float3 stepInscatter = lightContribution * scattering * density * stepSize;
            
            // Accumulate inscatter with proper transmittance weighting
            output.inscatter += stepInscatter * output.transmittance;
            output.transmittance *= stepTransmittance;

            // Advance ray only when we have meaningful calculations
            currentPos += rayDir * stepSize;
        }
        else
        {
            // Still advance ray even with negligible density
            currentPos += rayDir * stepSize;
        }

        // Early termination
        if (output.transmittance < 0.01)
            break;
    }

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
// RAINBOW EFFECT
// =============================================================================

// Rainbow color based on angle
float3 rainbowColor(float angle)
{
    // Primary rainbow at 42 degrees
    const float primaryAngle = 0.733; // radians (~42 degrees)
    const float rainbowWidth = 0.035;

    float dist = abs(angle - primaryAngle);
    float rainbow = smoothstep(rainbowWidth, 0.0, dist);

    if (rainbow <= 0.0)
        return 0.0;

    // Spectral colors
    float normalized = saturate((angle - (primaryAngle - rainbowWidth)) / (rainbowWidth * 2.0));

    float3 color;
    if (normalized < 0.166)
        color = lerp(float3(1.0, 0.0, 0.0), float3(1.0, 0.5, 0.0), normalized * 6.0);
    else if (normalized < 0.333)
        color = lerp(float3(1.0, 0.5, 0.0), float3(1.0, 1.0, 0.0), (normalized - 0.166) * 6.0);
    else if (normalized < 0.5)
        color = lerp(float3(1.0, 1.0, 0.0), float3(0.0, 1.0, 0.0), (normalized - 0.333) * 6.0);
    else if (normalized < 0.666)
        color = lerp(float3(0.0, 1.0, 0.0), float3(0.0, 0.0, 1.0), (normalized - 0.5) * 6.0);
    else if (normalized < 0.833)
        color = lerp(float3(0.0, 0.0, 1.0), float3(0.3, 0.0, 0.5), (normalized - 0.666) * 6.0);
    else
        color = lerp(float3(0.3, 0.0, 0.5), float3(0.5, 0.0, 0.5), (normalized - 0.833) * 6.0);

    return color * rainbow;
}

// Calculate rainbow contribution
float3 computeRainbow(float3 viewDir, float3 sunDir, float rainIntensity)
{
#if ENABLE_RAINBOW
    if (rainIntensity <= 0.0)
        return 0.0;

    // Rainbow appears opposite to sun
    float3 antiSunDir = -sunDir;

    // Angle from anti-sun direction
    float cosAngle = dot(viewDir, antiSunDir);
    float angle = acos(saturate(cosAngle));

    float3 rainbow = rainbowColor(angle);

    // Fade based on conditions
    float sunElevation = saturate(sunDir.y * 2.0); // Rainbow visible when sun is low
    float viewElevation = saturate(-viewDir.y + 0.5); // Rainbow in lower sky

    return rainbow * RAINBOW_INTENSITY * rainIntensity * sunElevation * viewElevation;
#else
    return 0.0;
#endif
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

    // Add rainbow
    result.inscatter += computeRainbow(rayDir, sunDir, rainIntensity) * result.transmittance;

    return result;
}

// Apply volumetrics to final color with improved light transmission
float3 applyVolumetrics(float3 sceneColor, VolumetricResult volumetrics)
{
    // Ensure transmittance doesn't go too dark (maintain minimum visibility)
    float minTransmittance = 0.1;
    float effectiveTransmittance = max(volumetrics.transmittance, minTransmittance);
    
    // Apply scene color through fog with proper energy conservation
    float3 transmittedScene = sceneColor * effectiveTransmittance;
    
    // Add inscatter light with wavelength-dependent scattering
    float3 finalColor = transmittedScene + volumetrics.inscatter;
    
    // Prevent oversaturation while maintaining light transmission
    return min(finalColor, sceneColor * 1.5 + volumetrics.inscatter * 1.2);
}

#endif // __OPENRTX_VOLUMETRIC_LIGHTING_HLSL__
