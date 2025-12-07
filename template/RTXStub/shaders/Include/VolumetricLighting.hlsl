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

// Combined fog density
float getFogDensity(float3 worldPos, float rainIntensity, float time)
{
    float baseDensity = fogDensityExponential(worldPos.y, 64.0); // Sea level

#if ENABLE_RAIN_FOG
    // Increase density during rain
    baseDensity *= 1.0 + rainIntensity * RAIN_FOG_AMOUNT;
#endif

    // Optional: time-based density variation
    float timeVariation = sin(time * 0.1) * 0.1 + 1.0;
    baseDensity *= timeVariation;

    return baseDensity;
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

    float stepSize = min(maxDistance, FOG_MAX_DISTANCE) / float(numSteps);
    float3 currentPos = rayOrigin;

    // Phase function for sun
    float cosTheta = dot(rayDir, sunDir);
    float phase = dualLobePhase(cosTheta, 0.8, -0.3, 0.7);

    // Scattering coefficient
    float3 scattering = float3(0.005, 0.006, 0.008); // Slightly blue-shifted

    for (int i = 0; i < numSteps; i++)
    {
        float density = getFogDensity(currentPos, rainIntensity, time);

        if (density > 0.0001)
        {
            // Extinction
            float extinction = density * stepSize;
            float stepTransmittance = exp(-extinction);

            // In-scattered light
            float3 lightContribution = 0.0;

#if ENABLE_SUN_FOG
            // Sun contribution
            float sunVisibility = 1.0; // Would need shadow sampling for accurate results
            lightContribution += sunColor * phase * sunVisibility * SUN_FOG_AMOUNT;
#endif

#if ENABLE_STATIC_GI_FOG
            // Ambient/GI contribution
            float3 ambientColor = float3(0.1, 0.12, 0.15); // Sky ambient
            lightContribution += ambientColor * STATIC_GI_FOG_AMOUNT;
#endif

            // Accumulate
            float3 stepInscatter = lightContribution * scattering * density * stepSize;
            output.inscatter += stepInscatter * output.transmittance;
            output.transmittance *= stepTransmittance;
        }

        currentPos += rayDir * stepSize;

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

// Apply volumetrics to final color
float3 applyVolumetrics(float3 sceneColor, VolumetricResult volumetrics)
{
    return sceneColor * volumetrics.transmittance + volumetrics.inscatter;
}

#endif // __OPENRTX_VOLUMETRIC_LIGHTING_HLSL__
