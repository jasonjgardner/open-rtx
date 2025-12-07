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
    float dayFactor = saturate(sunElevation * 4.0 + 0.5);

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

// PBR shading with all enhancements
// Uses game-provided lighting values for proper resource pack support
float3 shadeSurfacePBR(EnhancedSurface surface, OpenRTXContext ctx)
{
    float3 result = 0.0;

    float3 N = surface.normal;
    float3 V = surface.viewDir;

    // Compute F0
    float3 f0 = computeF0(surface.albedo, surface.metalness, 1.5);

    // Determine day/night factor from sun elevation
    float sunElevation = ctx.sunDir.y;
    float dayFactor = saturate(sunElevation * 4.0 + 0.5);

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

            // Energy conservation
            float3 F = fresnelSchlick(VdotH, f0);
            diffuse *= (1.0 - F);

            // Specular
            float3 specular = evaluateSpecularBRDF(f0, surface.roughness, NdotV, NdotL, NdotH, VdotH);

            // Sun light using game-provided color with reasonable intensity
            // Scale to match vanilla brightness expectations (vanilla uses 0.45-1.0 range)
            float3 sunLight = ctx.sunColor * dayFactor * 2.0;  // Multiplier ~2 for PBR response
            result += (diffuse + specular) * sunLight * NdotL * surface.ao;
        }
    }

    // Moon lighting (much dimmer)
    float nightFactor = 1.0 - dayFactor;
    if (nightFactor > 0.01)
    {
        float3 L = ctx.moonDir;
        float NdotL = saturate(dot(N, L));

        if (NdotL > 0.0)
        {
            // Use game moon intensity, scaled appropriately
            float3 moonLight = float3(0.6, 0.65, 0.8) * nightFactor * 0.3;

            // Simplified shading for moonlight
            float3 diffuse = surface.albedo * (1.0 - surface.metalness) * kInvPi;
            result += diffuse * moonLight * NdotL * surface.ao;
        }
    }

    // Ambient lighting from game
    result += surface.albedo * ctx.constantAmbient * surface.ao;

    // Emissive (scaled reasonably)
    if (surface.emissive > 0.0)
    {
        result += surface.albedo * surface.emissive * 2.0;
    }

    return result;
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
    gameBaseSky *= ctx.skyIntensityAdjustment;

    // Add constant ambient
    gameBaseSky += ctx.constantAmbient;

#if OPENRTX_ENABLED && ENABLE_ATMOSPHERIC_SKY
    // Enhanced atmospheric scattering (blended with game sky)
    SkyOutput sky = evaluateSky(rayDir, ctx.sunDir, ctx.moonDir, ctx.time, true);

    // Blend physical sky with game-provided colors
    // Use game colors as base, add atmospheric enhancement
    float sunElevation = ctx.sunDir.y;
    float dayFactor = saturate(sunElevation * 4.0 + 0.5);

    // During day: blend atmospheric scattering with game sky
    // During night: use game sky colors more directly
    float atmosphericBlend = dayFactor * 0.5; // 50% atmospheric during day
    skyColor = lerp(gameBaseSky, sky.color, atmosphericBlend);

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
    CloudOutput clouds = renderVolumetricClouds(
        ctx.viewOrigin,
        rayDir,
        ctx.sunDir,
        ctx.sunColor * ctx.sunIntensity,
        ctx.time,
        10000.0);

    // Blend clouds with sky
    skyColor = skyColor * clouds.transmittance + clouds.color;

    // Cirrus layer
    float3 cirrus = renderCirrusClouds(rayDir, ctx.sunDir, ctx.sunColor, ctx.time);
    skyColor += cirrus * 0.3;
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
