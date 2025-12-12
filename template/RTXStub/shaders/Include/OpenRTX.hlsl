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
// OpenRTX - Main Include File
// This is the primary include that brings together all OpenRTX systems
// =============================================================================

#ifndef __OPENRTX_HLSL__
#define __OPENRTX_HLSL__

// Core settings (must be first)
#include "Settings.hlsl"

// Graphics utilities
#include "GFXHelpers.hlsl"

// Lighting models
#include "BRDF.hlsl"

// Environment rendering
#include "Sky.hlsl"
#include "Clouds.hlsl"
#include "Water.hlsl"
#include "VolumetricLighting.hlsl"

// Post-processing
#include "ToneMapping.hlsl"

// =============================================================================
// OPENRTX RENDERING CONTEXT
// =============================================================================

struct OpenRTXContext
{
    // View info
    float3 viewOrigin;
    float3 viewDir;
    float2 screenUV;

    // Time and animation
    float time;
    float deltaTime;
    uint frameIndex;

    // Environment
    float3 sunDir;
    float3 moonDir;
    float3 sunColor;
    float sunIntensity;
    float moonIntensity;

    // Game-provided sky colors (from resource packs)
    float3 gameSkyColor;
    float3 gameSkyColorUp;
    float3 gameSkyColorDown;
    float skyIntensityAdjustment;
    float3 constantAmbient;

    // Weather
    float rainIntensity;
    float thunderIntensity;

    // Dimension (0=overworld, 1=nether, 2=end)
    int dimension;

    // Exposure
    float exposureEV;
    float avgLuminance;
};

// Initialize context from game state
// This version uses game-provided values from g_view for proper resource pack support
OpenRTXContext initContextFromGame(
    float3 viewOrigin,
    float3 viewDir,
    float2 screenUV,
    float time,
    float3 sunDir,
    float3 gameSunColor,
    float3 gameSkyColor,
    float3 gameSkyColorUp,
    float3 gameSkyColorDown,
    float3 constantAmbient,
    float skyIntensityAdjustment,
    float sunMeshIntensity,
    float moonMeshIntensity,
    float rainIntensity,
    int dimension)
{
    OpenRTXContext ctx;

    ctx.viewOrigin = viewOrigin;
    ctx.viewDir = viewDir;
    ctx.screenUV = screenUV;
    ctx.time = time;
    ctx.deltaTime = 0.016; // Assume 60fps
    ctx.frameIndex = uint(time * 60.0);

    ctx.sunDir = sunDir;
    ctx.moonDir = -sunDir;

    // Use game-provided sun color (from resource pack)
    ctx.sunColor = gameSunColor;
    ctx.sunIntensity = sunMeshIntensity;
    ctx.moonIntensity = moonMeshIntensity;

    // Store game-provided sky colors
    ctx.gameSkyColor = gameSkyColor;
    ctx.gameSkyColorUp = gameSkyColorUp;
    ctx.gameSkyColorDown = gameSkyColorDown;
    ctx.skyIntensityAdjustment = skyIntensityAdjustment;
    ctx.constantAmbient = constantAmbient;

    ctx.rainIntensity = rainIntensity;
    ctx.thunderIntensity = 0.0;

    ctx.dimension = dimension;

    // Default exposure
    ctx.exposureEV = 0.0;
    ctx.avgLuminance = 0.18;

    return ctx;
}

// Legacy initialization (for backwards compatibility)
OpenRTXContext initContext(
    float3 viewOrigin,
    float3 viewDir,
    float2 screenUV,
    float time,
    float3 sunDir,
    float rainIntensity,
    int dimension)
{
    OpenRTXContext ctx;

    ctx.viewOrigin = viewOrigin;
    ctx.viewDir = viewDir;
    ctx.screenUV = screenUV;
    ctx.time = time;
    ctx.deltaTime = 0.016; // Assume 60fps
    ctx.frameIndex = uint(time * 60.0);

    ctx.sunDir = sunDir;
    ctx.moonDir = -sunDir;

    // Determine sun color based on elevation (fallback when game data unavailable)
    float sunElevation = sunDir.y;
    float dayFactor = saturate(sunElevation * 2.0 + 0.5);  // Smoother transition

    // Color temperature based on sun position
    float colorTemp = lerp(SUNRISE_COLOR_TEMPERATURE, SUN_COLOR_TEMPERATURE, saturate(sunElevation * 2.0));
    ctx.sunColor = blackbodyColor(colorTemp);
    ctx.sunIntensity = SUN_INTENSITY * dayFactor;
    ctx.moonIntensity = MOON_INTENSITY * (1.0 - dayFactor);

    // Default sky colors (vanilla-like)
    ctx.gameSkyColor = float3(0.5, 0.7, 1.0);
    ctx.gameSkyColorUp = float3(0.4, 0.6, 0.9);
    ctx.gameSkyColorDown = float3(0.7, 0.8, 1.0);
    ctx.skyIntensityAdjustment = 1.0;
    ctx.constantAmbient = float3(0.1, 0.1, 0.1);

    ctx.rainIntensity = rainIntensity;
    ctx.thunderIntensity = 0.0;

    ctx.dimension = dimension;

    // Default exposure
    ctx.exposureEV = 0.0;
    ctx.avgLuminance = 0.18;

    return ctx;
}

// =============================================================================
// ENHANCED SURFACE SHADING
// =============================================================================

struct EnhancedSurface
{
    float3 position;
    float3 normal;
    float3 geometryNormal;
    float3 viewDir;

    float3 albedo;
    float roughness;
    float metalness;
    float ao;
    float subsurface;
    float emissive;

    bool isWater;
    float waterDepth;
};

// =============================================================================
// TIME-OF-DAY HELPERS
// =============================================================================

// Calculate time-of-day lighting scale from sun position and color
// Returns a value that smoothly transitions between day and night
float getTimeOfDayScale(float3 sunDir, float3 sunColor)
{
    float sunBrightness = dot(sunColor, float3(0.299, 0.587, 0.114));
    float dayFactor = saturate(sunDir.y * 2.0 + 0.5);  // Smoother transition
    return max(0.02, dayFactor * saturate(sunBrightness * 0.5 + 0.5));
}

// Get raw day factor (0 = night, 1 = day) without brightness adjustment
float getDayFactor(float3 sunDir)
{
    return saturate(sunDir.y * 2.0 + 0.5);
}

// =============================================================================
// PBR SURFACE SHADING
// =============================================================================

// Sample cloud shadow attenuation for a surface position
// Returns 0.0 = fully shadowed by clouds, 1.0 = no cloud shadow
float sampleCloudShadowAttenuation(float3 worldPos, float3 sunDir, float time)
{
#if ENABLE_VOLUMETRIC_CLOUDS && CLOUD_SHADOW_STRENGTH > 0.0
    // Sample cloud density above this position toward the sun
    float cloudShadow = sampleCloudShadow(worldPos, sunDir, time);
    // Apply shadow strength setting
    return lerp(1.0, cloudShadow, CLOUD_SHADOW_STRENGTH);
#else
    return 1.0;
#endif
}

// PBR direct lighting only (for shadow application)
// Returns ONLY direct sun/moon lighting, not ambient
float3 shadeSurfaceDirectPBR(EnhancedSurface surface, OpenRTXContext ctx)
{
    float3 result = 0.0;

    float3 N = surface.normal;
    float3 V = surface.viewDir;

    // Sample cloud shadows for this surface position
    float cloudShadow = sampleCloudShadowAttenuation(surface.position, ctx.sunDir, ctx.time);

    // Clouds are fully diffuse - skip all specular calculations
    if (surface.roughness >= 0.99)
    {
        // Pure diffuse lighting for clouds
        float3 L = ctx.sunDir;
        float NdotL = saturate(dot(N, L));

        if (NdotL > 0.0)
        {
            float dayFactor = getDayFactor(ctx.sunDir);
            float timeScale = getTimeOfDayScale(ctx.sunDir, ctx.sunColor);

            if (dayFactor > 0.01)
            {
                // Lambertian diffuse only
                float3 diffuse = surface.albedo * kInvPi * NdotL;
                float3 sunLight = ctx.sunColor * timeScale * 3.0;
                // Cloud shadows don't affect cloud surfaces themselves
                result += diffuse * sunLight * surface.ao;
            }
        }

        return result;
    }

    // Compute F0
    float3 f0 = computeF0(surface.albedo, surface.metalness, 1.5);

    // Get time-of-day factors
    float dayFactor = getDayFactor(ctx.sunDir);
    float timeScale = getTimeOfDayScale(ctx.sunDir, ctx.sunColor);

    // Direct sun lighting
    {
        float3 L = ctx.sunDir;
        float3 H = normalize(V + L);

        float NdotL = saturate(dot(N, L));
        float NdotV = saturate(dot(N, V));
        float NdotH = saturate(dot(N, H));
        float VdotH = saturate(dot(V, H));

        if (NdotL > 0.0 && dayFactor > 0.01)
        {
            // Diffuse
            float3 diffuseAlbedo = surface.albedo * (1.0 - surface.metalness);

#if ENABLE_SUBSURFACE_SCATTERING
            float3 diffuse;
            if (surface.subsurface > 0.0)
            {
                diffuse = subsurfaceDisney(diffuseAlbedo, surface.roughness, NdotV, NdotL, VdotH, surface.subsurface);
            }
            else
            {
                diffuse = evaluateDiffuseBRDF(diffuseAlbedo, surface.roughness, NdotV, NdotL, VdotH, N, V, L);
            }
#else
            float3 diffuse = evaluateDiffuseBRDF(diffuseAlbedo, surface.roughness, NdotV, NdotL, VdotH, N, V, L);
#endif

            // Energy conservation - use metal-aware Fresnel to preserve metal color
            float3 F = fresnelSchlickMetal(VdotH, f0, surface.albedo, surface.metalness);
            diffuse *= (1.0 - F);

            // Specular - use metal-aware BRDF to preserve metal color at grazing angles
            float3 specular = evaluateSpecularBRDFMetal(f0, surface.albedo, surface.metalness,
                                                         surface.roughness, NdotV, NdotL, NdotH, VdotH);

            // Sun light using game-provided color with time-of-day scaling
            float3 sunLight = ctx.sunColor * timeScale * 3.0;  // Scale for PBR response

            // Apply cloud shadow attenuation to direct sunlight
            // This creates soft cloud shadows on terrain and objects
            result += (diffuse + specular) * sunLight * NdotL * surface.ao * cloudShadow;
        }
    }

    // Moon lighting (much dimmer) - also direct lighting
    // Cloud shadows also affect moonlight (use same cloudShadow as sun, reversed direction)
    float nightFactor = 1.0 - dayFactor;
    if (nightFactor > 0.01)
    {
        float3 L = ctx.moonDir;
        float NdotL = saturate(dot(N, L));

        if (NdotL > 0.0)
        {
            float3 moonLight = float3(0.6, 0.65, 0.8) * nightFactor * 0.3;
            float3 diffuse = surface.albedo * (1.0 - surface.metalness) * kInvPi;

            // Sample cloud shadow for moon direction (clouds block moonlight too)
            float moonCloudShadow = sampleCloudShadowAttenuation(surface.position, ctx.moonDir, ctx.time);
            result += diffuse * moonLight * NdotL * surface.ao * moonCloudShadow;
        }
    }

    return result;
}

// Ambient/indirect lighting (never shadowed)
// This approximates what the irradiance cache would provide:
// - Hemisphere-based sky irradiance
// - Ground bounce contribution
// - Horizontal fill light for vertical surfaces
// - Time-of-day scaling
float3 shadeSurfaceAmbientPBR(EnhancedSurface surface, OpenRTXContext ctx)
{
    float3 result = 0.0;

    // Clouds are fully diffuse - simplified ambient lighting
    if (surface.roughness >= 0.99)
    {
        // Simple Lambertian ambient for clouds
        float3 ambient = ctx.constantAmbient;
        ambient = max(ambient, 0.08);  // Increased to match main ambient minimum
        
        float dayFactor = getDayFactor(ctx.sunDir);
        float timeScale = getTimeOfDayScale(ctx.sunDir, ctx.sunColor);
        
        // Basic ambient with time-of-day scaling
        result = surface.albedo * ambient * timeScale * surface.ao;
        return result;
    }

    // Game-provided ambient - ensure reasonable minimum for indoor/shadow areas
    float3 ambient = ctx.constantAmbient;
    ambient = max(ambient, 0.08);  // Further increased minimum to prevent crushed shadows

    // Time-of-day factor for indirect lighting
    float dayFactor = getDayFactor(ctx.sunDir);

    // === INDIRECT DIFFUSE (approximating irradiance cache) ===
    // In the full system, this would come from geometryInfo.incomingIrradiance
    // Here we approximate with hemisphere-based sky/ground/horizontal contribution

    // Normal-based weights for different irradiance sources
    float upDot = dot(surface.normal, float3(0, 1, 0));
    float skyWeight = saturate(upDot * 0.5 + 0.5);       // 0.5-1.0 for upward faces
    float groundWeight = saturate(-upDot * 0.5 + 0.5);   // 0.5-1.0 for downward faces
    float horizontalWeight = 1.0 - abs(upDot);            // 1.0 for vertical, 0.0 for up/down

    // Sky hemisphere contribution
    float3 skyIrradiance = lerp(ctx.gameSkyColorDown, ctx.gameSkyColorUp, skyWeight);
    skyIrradiance = max(skyIrradiance, 0.03);

    // Ground bounce contribution
    float3 groundColor = float3(0.12, 0.10, 0.08);  // Warmer ground albedo
    float3 groundIrradiance = groundColor * ctx.sunColor * dayFactor * 0.6;  // Increased from 0.4

    // Horizontal fill light (for vertical faces like walls)
    // Approximates light bouncing from surrounding environment
    float3 horizonColor = lerp(ctx.gameSkyColorDown, ctx.gameSkyColor, 0.5);
    float3 horizontalIrradiance = horizonColor * 0.4;  // Increased from 0.3

    // Combine all indirect sources
    float3 indirectIrradiance = skyIrradiance * 0.2 * skyWeight
                              + groundIrradiance * groundWeight
                              + horizontalIrradiance * horizontalWeight;

    // Ensure minimum indirect for all surfaces (prevents crushed shadows)
    indirectIrradiance = max(indirectIrradiance, 0.04);  // Increased from 0.02

    // Scale with time of day (but maintain minimum at night)
    indirectIrradiance *= dayFactor * 0.6 + 0.4;

    // Apply to diffuse albedo (following GI pattern: irradiance * diffuseColour)
    result += surface.albedo * (ambient + indirectIrradiance) * surface.ao;

    // === EMISSIVE (never shadowed) ===
    // Emissive surfaces emit light directly
    if (surface.emissive > 0.0)
    {
        result += surface.albedo * surface.emissive * EMISSIVE_INTENSITY;
    }

    return result;
}

// PBR shading with all enhancements (combined for backwards compatibility)
// Uses game-provided lighting values for proper resource pack support
float3 shadeSurfacePBR(EnhancedSurface surface, OpenRTXContext ctx)
{
    return shadeSurfaceDirectPBR(surface, ctx) + shadeSurfaceAmbientPBR(surface, ctx);
}

// =============================================================================
// SKY AND ATMOSPHERE INTEGRATION
// =============================================================================

float3 renderSkyWithClouds(float3 rayDir, OpenRTXContext ctx)
{
    float3 skyColor = 0.0;

    // Use game-provided sky gradient as the base (respects resource packs)
    float gradientT = saturate(rayDir.y * 0.5 + 0.5);
    float3 gameBaseSky = lerp(ctx.gameSkyColorDown, ctx.gameSkyColorUp, gradientT);

    // Check if game sky colors are valid (non-zero)
    float gameSkyLuminance = dot(gameBaseSky, float3(0.299, 0.587, 0.114));

    // If game sky is too dark, use a fallback gradient
    if (gameSkyLuminance < 0.01)
    {
        // Fallback: procedural sky gradient based on sun position
        float dayFactor = getDayFactor(ctx.sunDir);
        float3 dayZenith = float3(0.2, 0.4, 0.8);   // Blue sky
        float3 dayHorizon = float3(0.6, 0.7, 0.9);  // Lighter at horizon
        float3 nightZenith = float3(0.01, 0.01, 0.02);
        float3 nightHorizon = float3(0.02, 0.02, 0.04);

        float3 zenith = lerp(nightZenith, dayZenith, dayFactor);
        float3 horizon = lerp(nightHorizon, dayHorizon, dayFactor);
        gameBaseSky = lerp(horizon, zenith, saturate(rayDir.y));
    }

    // Apply sky intensity adjustment from game, with minimum brightness
    float effectiveSkyIntensity = max(ctx.skyIntensityAdjustment, 0.5);
    gameBaseSky *= effectiveSkyIntensity;

    // Add constant ambient and ensure minimum sky brightness
    gameBaseSky += ctx.constantAmbient;
    gameBaseSky = max(gameBaseSky, 0.02);  // Prevent completely black sky

#if OPENRTX_ENABLED && ENABLE_ATMOSPHERIC_SKY
    // Enhanced atmospheric scattering (additive, not replacing)
    SkyOutput sky = evaluateSky(rayDir, ctx.sunDir, ctx.moonDir, ctx.time, true);

    // Ensure atmospheric sky also has minimum brightness
    sky.color = max(sky.color, 0.01);

    // Blend physical sky with game-provided colors
    float skyDayFactor = getDayFactor(ctx.sunDir);

    // Use game sky as base, blend in atmospheric enhancement
    // Lower blend factor to preserve game sky colors better
    float atmosphericBlend = skyDayFactor * 0.3; // 30% atmospheric during day
    skyColor = lerp(gameBaseSky, gameBaseSky + sky.color * 0.5, atmosphericBlend);

    // Always add sun disk from physical model
    skyColor += sky.sunDiskColor;
#else
    // Use game-provided sky directly (vanilla-like)
    skyColor = gameBaseSky;

    // Simple sun glow
    float sunCosAngle = dot(rayDir, ctx.sunDir);
    float sunGlow = pow(saturate(sunCosAngle), 32.0);
    skyColor += ctx.sunColor * sunGlow * ctx.sunIntensity * 0.1;
#endif

#if ENABLE_VOLUMETRIC_CLOUDS
    // Volumetric clouds
    // Use time-of-day aware intensity for cloud illumination
    float cloudDayFactor = getDayFactor(ctx.sunDir);
    float cloudTimeScale = getTimeOfDayScale(ctx.sunDir, ctx.sunColor);
    float3 cloudSunColor = ctx.sunColor * cloudTimeScale * 2.0;  // Moderate multiplier

    CloudOutput clouds = renderVolumetricClouds(
        ctx.viewOrigin,
        rayDir,
        ctx.sunDir,
        cloudSunColor,
        ctx.time,
        10000.0);

    // Clamp cloud color to prevent blown-out clouds
    clouds.color = min(clouds.color, 2.0);

    // Blend clouds with sky - use proper alpha blending
    // Higher minimum transmittance to keep sky visible
    float minTransmittance = 0.3;
    float effectiveTransmittance = max(clouds.transmittance, minTransmittance);

    // Blend: sky shows through based on transmittance, clouds add on top
    skyColor = skyColor * effectiveTransmittance + clouds.color * (1.0 - clouds.transmittance);

    // Cirrus layer (very subtle)
    float3 cirrus = renderCirrusClouds(rayDir, ctx.sunDir, ctx.sunColor, ctx.time);
    skyColor += cirrus * 0.15;
#endif

    return skyColor;
}

// =============================================================================
// WATER RENDERING INTEGRATION
// =============================================================================

float3 renderWaterSurface(float3 worldPos, float3 viewDir, OpenRTXContext ctx, float3 underwaterColor)
{
#if ENABLE_WATER_EFFECTS
    WaterSurface water = evaluateWaterSurface(worldPos, viewDir, ctx.time, ctx.rainIntensity);

    // Reflection
    float3 reflectDir = reflect(-viewDir, water.normal);
    float3 reflection = renderSkyWithClouds(reflectDir, ctx);

    // Refraction color (from underwater scene)
    float3 refraction = underwaterColor;

    // Apply water absorption to refraction
    float refractionDepth = 2.0; // Approximate depth
    refraction = applyUnderwaterFog(refraction, refractionDepth, ctx.sunColor, ctx.sunIntensity);

    // Blend based on Fresnel
    float3 waterColor = lerp(refraction, reflection, water.fresnel);

    // Add foam
    waterColor = lerp(waterColor, float3(1.0, 1.0, 1.0), water.foam);

    // Specular highlight on water
    float3 H = normalize(viewDir + ctx.sunDir);
    float NdotH = saturate(dot(water.normal, H));
    float specular = pow(NdotH, 256.0) * ctx.sunIntensity;
    waterColor += ctx.sunColor * specular * (1.0 - water.foam);

    return waterColor;
#else
    return underwaterColor;
#endif
}

// =============================================================================
// VOLUMETRIC INTEGRATION
// =============================================================================

float3 applyAtmosphericEffects(float3 sceneColor, float3 rayDir, float hitDistance, OpenRTXContext ctx)
{
    float3 result = sceneColor;

#if ENABLE_VOLUMETRIC_LIGHTING
    VolumetricResult volumetrics = evaluateVolumetrics(
        ctx.viewOrigin,
        rayDir,
        hitDistance,
        ctx.sunDir,
        ctx.sunColor * ctx.sunIntensity,
        ctx.rainIntensity,
        ctx.time,
        ctx.dimension);

    result = applyVolumetrics(result, volumetrics);
#endif

    return result;
}

// =============================================================================
// FINAL COLOR OUTPUT
// =============================================================================

float3 finalizeColor(float3 color, float2 screenUV, OpenRTXContext ctx)
{
#if ENABLE_POST_PROCESSING
    // Apply tonemapping and post-processing
    color = postProcessSimple(color, screenUV, ctx.exposureEV);
#else
    // Just tonemap
    color = tonemap(color);
    color = linearToSRGB(color);
#endif

#if DEBUG_NAN_INF
    // Debug: highlight NaN/Inf
    if (any(isnan(color)) || any(isinf(color)))
    {
        int pattern = int(screenUV.x * 32.0) + int(screenUV.y * 32.0);
        color = (pattern & 1) ? float3(1, 0, 1) : float3(0, 0, 0);
    }
#endif

    return saturate(color);
}

// =============================================================================
// DEBUG VIEWS
// =============================================================================

float3 renderDebugView(EnhancedSurface surface, float depth, float2 motionVector, OpenRTXContext ctx)
{
#if DEBUG_VIEW == 1
    return surface.albedo;
#elif DEBUG_VIEW == 2
    return visualizeNormal(surface.normal);
#elif DEBUG_VIEW == 3
    return surface.roughness;
#elif DEBUG_VIEW == 4
    return surface.metalness;
#elif DEBUG_VIEW == 5
    return surface.emissive;
#elif DEBUG_VIEW == 6
    return surface.ao;
#elif DEBUG_VIEW == 7
    return visualizeDepth(depth, 0.1, 1000.0);
#elif DEBUG_VIEW == 8
    return float3(motionVector * 10.0 + 0.5, 0.0);
#else
    return 0.0;
#endif
}

#endif // __OPENRTX_HLSL__
