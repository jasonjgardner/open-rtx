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
// OpenRTX Sky - Physically-Based Atmospheric Scattering
// Implements Rayleigh and Mie scattering with closed-form approximations
// Based on: Hillaire's "A Scalable and Production Ready Sky and Atmosphere Rendering"
// =============================================================================

#ifndef __OPENRTX_SKY_HLSL__
#define __OPENRTX_SKY_HLSL__

#include "Settings.hlsl"

// =============================================================================
// ATMOSPHERIC CONSTANTS
// =============================================================================

// Scattering coefficients
static const float3 kRayleighScattering = float3(RAYLEIGH_SCATTERING_R, RAYLEIGH_SCATTERING_G, RAYLEIGH_SCATTERING_B);
static const float3 kMieScattering = float3(MIE_SCATTERING, MIE_SCATTERING, MIE_SCATTERING);
static const float3 kMieAbsorption = float3(MIE_ABSORPTION, MIE_ABSORPTION, MIE_ABSORPTION);

// Ozone absorption (optional for more realistic atmosphere)
static const float3 kOzoneAbsorption = float3(0.65e-6, 1.881e-6, 0.085e-6);
static const float kOzoneLayerHeight = 25000.0; // meters
static const float kOzoneLayerWidth = 15000.0;  // meters

// =============================================================================
// PHASE FUNCTIONS
// =============================================================================

// Rayleigh phase function
float phaseRayleigh(float cosTheta)
{
    return (3.0 / (16.0 * kPi)) * (1.0 + cosTheta * cosTheta);
}

// Henyey-Greenstein phase function for Mie scattering
float phaseMieHG(float cosTheta, float g)
{
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return (1.0 - g2) / (4.0 * kPi * pow(abs(denom), 1.5));
}

// Cornette-Shanks phase function (improved Mie)
float phaseMieCS(float cosTheta, float g)
{
    float g2 = g * g;
    float num = 3.0 * (1.0 - g2) * (1.0 + cosTheta * cosTheta);
    float denom = (8.0 * kPi) * (2.0 + g2) * pow(abs(1.0 + g2 - 2.0 * g * cosTheta), 1.5);
    return num / max(denom, 1e-7);
}

// Schlick approximation for phase function (faster)
float phaseMieSchlick(float cosTheta, float g)
{
    float k = 1.55 * g - 0.55 * g * g * g;
    float denom = 1.0 - k * cosTheta;
    return (1.0 - k * k) / (4.0 * kPi * denom * denom);
}

// =============================================================================
// RAY-SPHERE INTERSECTION
// =============================================================================

// Returns distances to both intersections, or negative if no intersection
float2 raySphereIntersect(float3 rayOrigin, float3 rayDir, float3 sphereCenter, float sphereRadius)
{
    float3 oc = rayOrigin - sphereCenter;
    float a = dot(rayDir, rayDir);
    float b = 2.0 * dot(oc, rayDir);
    float c = dot(oc, oc) - sphereRadius * sphereRadius;
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0)
        return float2(-1.0, -1.0);

    float sqrtD = sqrt(discriminant);
    return float2((-b - sqrtD) / (2.0 * a), (-b + sqrtD) / (2.0 * a));
}

// =============================================================================
// DENSITY FUNCTIONS
// =============================================================================

// Exponential density falloff
float densityExponential(float height, float scaleHeight)
{
    return exp(-max(0.0, height) / scaleHeight);
}

// Ozone density (tent function)
float densityOzone(float height)
{
    return max(0.0, 1.0 - abs(height - kOzoneLayerHeight) / kOzoneLayerWidth);
}

// =============================================================================
// OPTICAL DEPTH CALCULATION
// =============================================================================

// Numerical integration for optical depth along a ray
float3 computeOpticalDepth(float3 rayOrigin, float3 rayDir, float rayLength, float3 planetCenter, int numSamples)
{
    float stepSize = rayLength / float(numSamples);
    float3 opticalDepth = 0.0;

    for (int i = 0; i < numSamples; i++)
    {
        float3 samplePos = rayOrigin + rayDir * (float(i) + 0.5) * stepSize;
        float height = length(samplePos - planetCenter) - EARTH_RADIUS;

        float densityR = densityExponential(height, RAYLEIGH_SCALE_HEIGHT);
        float densityM = densityExponential(height, MIE_SCALE_HEIGHT);
        float densityO = densityOzone(height);

        opticalDepth += (kRayleighScattering * densityR + kMieScattering * densityM + kOzoneAbsorption * densityO) * stepSize;
    }

    return opticalDepth;
}

// =============================================================================
// CLOSED-FORM APPROXIMATIONS
// =============================================================================

// Chapman function approximation for slant optical depth
float chapmanApprox(float X, float cosZenith)
{
    float c = sqrt(X * kHalfPi);
    if (cosZenith >= 0.0)
    {
        return c * rsqrt(X * cosZenith + 1.0) * exp(-X * (1.0 - cosZenith));
    }
    else
    {
        float sinZenith = sqrt(1.0 - cosZenith * cosZenith);
        return 2.0 * c * exp(X * (sinZenith - 1.0)) - c * rsqrt(1.0 - X * cosZenith) * exp(-X * (1.0 + cosZenith));
    }
}

// Approximate optical depth for exponential atmosphere
float3 approximateOpticalDepth(float height, float cosZenith)
{
    float XR = (EARTH_RADIUS + height) / RAYLEIGH_SCALE_HEIGHT;
    float XM = (EARTH_RADIUS + height) / MIE_SCALE_HEIGHT;

    float chR = chapmanApprox(XR, cosZenith);
    float chM = chapmanApprox(XM, cosZenith);

    float3 depthR = kRayleighScattering * RAYLEIGH_SCALE_HEIGHT * chR;
    float3 depthM = (kMieScattering + kMieAbsorption) * MIE_SCALE_HEIGHT * chM;

    return depthR + depthM;
}

// =============================================================================
// SKY COLOR COMPUTATION
// =============================================================================

// Compute sky color using single scattering
float3 computeSingleScatteringSky(float3 rayDir, float3 sunDir, float3 viewPos, int numSamples, int numLightSamples)
{
    float3 planetCenter = float3(0.0, -EARTH_RADIUS, 0.0);

    // Ray-atmosphere intersection
    float2 atmoIntersect = raySphereIntersect(viewPos, rayDir, planetCenter, ATMOSPHERE_RADIUS);
    if (atmoIntersect.y < 0.0)
        return 0.0; // No intersection

    // Ray-planet intersection (check if we hit the ground)
    float2 planetIntersect = raySphereIntersect(viewPos, rayDir, planetCenter, EARTH_RADIUS);
    float rayLength = (planetIntersect.x > 0.0) ? planetIntersect.x : atmoIntersect.y;
    rayLength = max(0.0, rayLength - max(0.0, atmoIntersect.x));

    float stepSize = rayLength / float(numSamples);
    float3 rayStart = viewPos + rayDir * max(0.0, atmoIntersect.x);

    float3 totalRayleigh = 0.0;
    float3 totalMie = 0.0;
    float3 opticalDepthView = 0.0;

    float cosTheta = dot(rayDir, sunDir);
    float phaseR = phaseRayleigh(cosTheta);
    float phaseM = phaseMieCS(cosTheta, MIE_ASYMMETRY);

    for (int i = 0; i < numSamples; i++)
    {
        float3 samplePos = rayStart + rayDir * (float(i) + 0.5) * stepSize;
        float height = length(samplePos - planetCenter) - EARTH_RADIUS;

        // Local densities
        float densityR = densityExponential(height, RAYLEIGH_SCALE_HEIGHT);
        float densityM = densityExponential(height, MIE_SCALE_HEIGHT);

        // Optical depth from view to sample point
        opticalDepthView += (kRayleighScattering * densityR + (kMieScattering + kMieAbsorption) * densityM) * stepSize;

        // Optical depth from sample to sun
        float3 opticalDepthSun;
        {
            float2 sunAtmoIntersect = raySphereIntersect(samplePos, sunDir, planetCenter, ATMOSPHERE_RADIUS);
            float sunRayLength = sunAtmoIntersect.y;

            // Quick check for shadow
            float2 sunPlanetIntersect = raySphereIntersect(samplePos, sunDir, planetCenter, EARTH_RADIUS);
            if (sunPlanetIntersect.x > 0.0)
            {
                continue; // In shadow
            }

            float sunStepSize = sunRayLength / float(numLightSamples);
            opticalDepthSun = 0.0;

            for (int j = 0; j < numLightSamples; j++)
            {
                float3 sunSamplePos = samplePos + sunDir * (float(j) + 0.5) * sunStepSize;
                float sunHeight = length(sunSamplePos - planetCenter) - EARTH_RADIUS;

                float sunDensityR = densityExponential(sunHeight, RAYLEIGH_SCALE_HEIGHT);
                float sunDensityM = densityExponential(sunHeight, MIE_SCALE_HEIGHT);

                opticalDepthSun += (kRayleighScattering * sunDensityR + (kMieScattering + kMieAbsorption) * sunDensityM) * sunStepSize;
            }
        }

        // Transmittance
        float3 transmittance = exp(-(opticalDepthView + opticalDepthSun));

        // Accumulate scattering
        totalRayleigh += transmittance * densityR * stepSize;
        totalMie += transmittance * densityM * stepSize;
    }

    // Final color
    float3 skyColor = SUN_INTENSITY * (totalRayleigh * kRayleighScattering * phaseR + totalMie * kMieScattering * phaseM);

    return skyColor;
}

// =============================================================================
// FAST SKY APPROXIMATION
// =============================================================================

// Closed-form sky color approximation (O(1) complexity)
float3 computeSkyColorFast(float3 rayDir, float3 sunDir, float height)
{
    float cosTheta = dot(rayDir, sunDir);
    float cosZenith = rayDir.y;

    // Phase functions
    float phaseR = phaseRayleigh(cosTheta);
    float phaseM = phaseMieSchlick(cosTheta, MIE_ASYMMETRY);

    // Approximate optical depth
    float3 opticalDepth = approximateOpticalDepth(height, cosZenith);

    // Transmittance
    float3 transmittance = exp(-opticalDepth);

    // Single scattering approximation
    float3 rayleighColor = kRayleighScattering * phaseR * (1.0 - transmittance.r);
    float3 mieColor = kMieScattering * phaseM * (1.0 - exp(-opticalDepth * 0.1));

    return SUN_INTENSITY * (rayleighColor + mieColor);
}

// =============================================================================
// PRECOMPUTED SKY WITH GRADIENTS
// =============================================================================

// Gradient-based sky approximation (very fast)
float3 computeGradientSky(float3 rayDir, float3 sunDir, float timeOfDay)
{
    float cosTheta = dot(rayDir, sunDir);
    float zenith = rayDir.y;

    // Day colors
    static const float3 dayZenith = float3(0.3, 0.5, 0.85);
    static const float3 dayHorizon = float3(0.6, 0.7, 0.9);

    // Sunset colors
    static const float3 sunsetZenith = float3(0.15, 0.2, 0.4);
    static const float3 sunsetHorizon = float3(1.0, 0.4, 0.1);

    // Night colors
    static const float3 nightZenith = NIGHT_SKY_COLOR;
    static const float3 nightHorizon = float3(0.02, 0.03, 0.05);

    // Determine time of day from sun elevation
    float sunElevation = sunDir.y;
    float dayFactor = saturate(sunElevation * 3.0 + 0.3);
    float sunsetFactor = saturate(1.0 - abs(sunElevation) * 5.0);

    // Blend sky colors based on time
    float3 zenithColor = lerp(nightZenith, dayZenith, dayFactor);
    zenithColor = lerp(zenithColor, sunsetZenith, sunsetFactor * 0.5);

    float3 horizonColor = lerp(nightHorizon, dayHorizon, dayFactor);
    horizonColor = lerp(horizonColor, sunsetHorizon, sunsetFactor);

    // Vertical gradient
    float gradientFactor = saturate(zenith);
    gradientFactor = pow(gradientFactor, 0.5); // Adjust curve

    float3 skyColor = lerp(horizonColor, zenithColor, gradientFactor);

    // Sun glow
    float sunGlow = pow(saturate(cosTheta), 32.0);
    float3 sunColor = lerp(float3(1.0, 0.8, 0.6), float3(1.0, 0.3, 0.1), sunsetFactor);
    skyColor += sunGlow * sunColor * (1.0 - sunsetFactor * 0.5);

    // Mie scattering near sun
    float mieGlow = pow(saturate((cosTheta + 0.1) / 1.1), 8.0);
    float3 mieColor = lerp(float3(1.0, 0.9, 0.7), float3(1.0, 0.5, 0.2), sunsetFactor);
    skyColor += mieGlow * mieColor * 0.3;

    return skyColor;
}

// =============================================================================
// NIGHT SKY WITH STARS
// =============================================================================

// Hash function for procedural stars
float hashStar(float3 p)
{
    p = frac(p * float3(443.897, 441.423, 437.195));
    p += dot(p, p.yzx + 19.19);
    return frac((p.x + p.y) * p.z);
}

// Procedural star field
float3 computeStarField(float3 rayDir, float time)
{
    float3 stars = 0.0;

    // Create star grid
    float3 gridPos = rayDir * 500.0;
    float3 gridId = floor(gridPos);
    float3 gridFrac = frac(gridPos);

    // Check neighboring cells
    for (int x = -1; x <= 1; x++)
    {
        for (int y = -1; y <= 1; y++)
        {
            for (int z = -1; z <= 1; z++)
            {
                float3 offset = float3(x, y, z);
                float3 cellId = gridId + offset;

                // Random star position within cell
                float h = hashStar(cellId);
                float3 starPos = offset + float3(hashStar(cellId + 1.0), hashStar(cellId + 2.0), hashStar(cellId + 3.0)) - gridFrac;

                float dist = length(starPos);

                // Star brightness
                float brightness = smoothstep(0.05, 0.0, dist);
                brightness *= step(0.95, h); // Only some cells have stars

                // Star color variation
                float colorVariation = hashStar(cellId + 4.0);
                float3 starColor = lerp(float3(0.8, 0.85, 1.0), float3(1.0, 0.9, 0.8), colorVariation);

                // Twinkle
                float twinkle = sin(time * 2.0 + h * 100.0) * 0.3 + 0.7;

                stars += brightness * starColor * twinkle;
            }
        }
    }

    return stars * NIGHT_SKY_INTENSITY;
}

// =============================================================================
// COLOR TEMPERATURE (BLACKBODY RADIATION)
// =============================================================================

// Approximate blackbody color from temperature (Kelvin)
float3 blackbodyColor(float temperature)
{
    // Based on Tanner Helland's approximation
    float temp = temperature / 100.0;
    float3 color;

    // Red
    if (temp <= 66.0)
    {
        color.r = 1.0;
    }
    else
    {
        float t = temp - 60.0;
        color.r = 1.29293618606 * pow(t, -0.1332047592);
    }

    // Green
    if (temp <= 66.0)
    {
        color.g = 0.390081579 * log(temp) - 0.631841444;
    }
    else
    {
        float t = temp - 60.0;
        color.g = 1.12989086089 * pow(t, -0.0755148492);
    }

    // Blue
    if (temp >= 66.0)
    {
        color.b = 1.0;
    }
    else if (temp <= 19.0)
    {
        color.b = 0.0;
    }
    else
    {
        float t = temp - 10.0;
        color.b = 0.543206789 * log(t) - 1.19625408914;
    }

    return saturate(color);
}

// =============================================================================
// MAIN SKY FUNCTION
// =============================================================================

struct SkyOutput
{
    float3 color;
    float3 transmittance;
    float3 sunDiskColor;
};

SkyOutput evaluateSky(float3 rayDir, float3 sunDir, float3 moonDir, float time, bool includeSunDisk)
{
    SkyOutput output;
    output.color = 0.0;
    output.transmittance = 1.0;
    output.sunDiskColor = 0.0;

    // Determine day/night
    float sunElevation = sunDir.y;
    float dayFactor = saturate(sunElevation * 4.0 + 0.5);

#if ENABLE_ATMOSPHERIC_SKY
    // Use physically-based sky during day
    if (dayFactor > 0.01)
    {
        // Camera height above sea level (Minecraft Y=64 = sea level)
        float viewHeight = 1000.0; // Assume typical gameplay height

        // Choose quality level based on performance
    #if CLOUD_MARCH_STEPS > 32
        // High quality: full ray marching
        output.color = computeSingleScatteringSky(rayDir, sunDir, float3(0, EARTH_RADIUS + viewHeight, 0), 16, 8);
    #else
        // Lower quality: closed-form approximation
        output.color = computeSkyColorFast(rayDir, sunDir, viewHeight);
    #endif

        // Approximate transmittance for clouds
        float zenith = max(0.0, rayDir.y);
        output.transmittance = exp(-approximateOpticalDepth(viewHeight, zenith) * 0.5);
    }

    // Night sky
    if (dayFactor < 0.99)
    {
        float3 nightColor = computeGradientSky(rayDir, moonDir, time);

        // Add stars
        nightColor += computeStarField(rayDir, time) * (1.0 - dayFactor);

        // Moon glow
        float moonCosAngle = dot(rayDir, moonDir);
        float moonGlow = pow(saturate(moonCosAngle), 64.0);
        nightColor += moonGlow * float3(0.6, 0.65, 0.8) * MOON_INTENSITY;

        output.color = lerp(nightColor, output.color, dayFactor);
    }
#else
    // Fallback to simple gradient sky
    output.color = computeGradientSky(rayDir, sunDir, time);
#endif

    // Sun disk
    if (includeSunDisk && sunElevation > -0.1)
    {
        float sunCosAngle = dot(rayDir, sunDir);
        float sunAngularRadius = kSunRadiusRadians * SUN_RADIUS_MULTIPLIER;

        if (sunCosAngle > cos(sunAngularRadius))
        {
            // Sun limb darkening
            float centerDistance = acos(sunCosAngle) / sunAngularRadius;
            float limbDarkening = sqrt(max(0.0, 1.0 - centerDistance * centerDistance));

            // Temperature-based sun color
            float temp = lerp(SUNRISE_COLOR_TEMPERATURE, SUN_COLOR_TEMPERATURE, saturate(sunElevation * 2.0));
            float3 sunColor = blackbodyColor(temp);

            output.sunDiskColor = sunColor * SUN_INTENSITY * limbDarkening;
        }
    }

    return output;
}

// =============================================================================
// AERIAL PERSPECTIVE
// =============================================================================

// Compute aerial perspective (atmospheric scattering between camera and object)
float3 computeAerialPerspective(float3 surfaceColor, float3 rayDir, float distance, float3 sunDir, float height)
{
    if (distance <= 0.0)
        return surfaceColor;

    float cosZenith = rayDir.y;

    // Optical depth along view ray
    float3 extinction = approximateOpticalDepth(height, cosZenith);
    extinction *= distance / 10000.0; // Scale by distance

    // Transmittance
    float3 transmittance = exp(-extinction);

    // In-scattered light (simplified)
    float cosTheta = dot(rayDir, sunDir);
    float3 inscatter = computeSkyColorFast(rayDir, sunDir, height) * (1.0 - transmittance);

    return surfaceColor * transmittance + inscatter;
}

#endif // __OPENRTX_SKY_HLSL__
