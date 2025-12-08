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
// OpenRTX Final Combine Pass
// Combines denoised buffers and applies post-processing
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/OpenRTX.hlsl"

[numthreads(16, 16, 1)]
void FinalCombine(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID)
{
    // Decode buffer indices from root constant
    uint diffuseDenoisingBufferIndex = g_rootConstant0 & 0xff;
    uint specularDenoisingBufferIndex = (g_rootConstant0 >> 8) & 0xff;
    uint shadowDenoisingBufferIndex = (g_rootConstant0 >> 16) & 0xff;

    // Get pixel coordinates
    uint2 pixelCoord = dispatchThreadID.xy;

    // Check bounds
    if (any(pixelCoord >= g_view.renderResolution))
        return;

    // Calculate UV
    float2 uv = (float2(pixelCoord) + 0.5) * g_view.recipRenderResolution;

    // Read the current color (from ray tracing pass)
    float4 colorAlpha = outputBufferFinal[pixelCoord];
    float3 color = colorAlpha.rgb;

    // Read depth for post-processing
    float depth = outputBufferReprojectedPathLength[pixelCoord];

    // Read motion vectors for effects
    float2 motionVector = outputBufferMotionVectors[pixelCoord];

    // Initialize OpenRTX context for post-processing
    OpenRTXContext ctx;
    ctx.time = g_view.time;
    ctx.screenUV = uv;
    ctx.exposureEV = 0.0; // Would come from auto-exposure

#if AUTO_EXPOSURE_ENABLED && !LOCK_EXPOSURE
    // Auto exposure calculation would go here
    // For now, use default
    ctx.exposureEV = 0.0;
#else
    ctx.exposureEV = LOCKED_EXPOSURE_VALUE;
#endif

    // Apply post-processing
#if ENABLE_POST_PROCESSING

    // ---------------------------
    // Bloom (in HDR before tonemapping)
    // ---------------------------
#if ENABLE_BLOOM
    {
        // Extract and add bloom contribution
        float3 bloom = computeBloomSinglePass(color, uv, float2(g_view.renderResolution));

        // Apply anamorphic effect if enabled
        bloom = computeAnamorphicBloom(bloom, uv);

        // Add bloom to color
        color += bloom;
    }
#endif

    // ---------------------------
    // Tonemapping
    // ---------------------------
    color = applyExposure(color, ctx.exposureEV);
    color = tonemap(color);

    // Convert to gamma space for further processing
    color = linearToSRGB(color);

    // ---------------------------
    // Film Grain (applied after tonemapping in LDR gamma space)
    // ---------------------------
#if ENABLE_FILM_GRAIN
    color = applyFilmGrain(color, uv, g_view.time, FILM_GRAIN_INTENSITY);
#endif

    // ---------------------------
    // Vignette
    // ---------------------------
#if ENABLE_VIGNETTE
    {
        float2 center = uv - 0.5;
        float dist = length(center);
        float vignette = smoothstep(VIGNETTE_RADIUS, VIGNETTE_RADIUS - VIGNETTE_SOFTNESS, dist);
        color *= lerp(1.0 - VIGNETTE_INTENSITY, 1.0, vignette);
    }
#endif

    // ---------------------------
    // Underwater Post-Processing
    // ---------------------------
    bool isUnderwater = g_view.cameraIsUnderWater != 0;
#if ENABLE_UNDERWATER_DISTORTION
    if (isUnderwater)
    {
        // Apply underwater chromatic aberration (wavelength-dependent IOR)
        float3 chromaOffset = getUnderwaterChromaticOffset(uv);
        if (any(chromaOffset != 0))
        {
            float2 center = uv - 0.5;
            float2 dir = normalize(center + 0.0001);

            // Sample at offset positions for R and B channels
            // Note: This is approximate since we don't have a separate texture sample
            // The actual offset is applied through the ray distortion in primary pass
            // This adds additional color fringing for enhanced effect
            float2 uvR = uv + dir * chromaOffset.r;
            float2 uvB = uv + dir * chromaOffset.b;

            // Clamp UVs to valid range
            uvR = saturate(uvR);
            uvB = saturate(uvB);

            // Read offset pixels (approximate - uses neighbor pixels)
            int2 pixelR = int2(uvR * float2(g_view.renderResolution));
            int2 pixelB = int2(uvB * float2(g_view.renderResolution));
            pixelR = clamp(pixelR, int2(0, 0), int2(g_view.renderResolution) - 1);
            pixelB = clamp(pixelB, int2(0, 0), int2(g_view.renderResolution) - 1);

            // Blend in chromatic offset
            float3 colorR = outputBufferFinal[pixelR].rgb;
            float3 colorB = outputBufferFinal[pixelB].rgb;
            color.r = lerp(color.r, colorR.r, 0.5);
            color.b = lerp(color.b, colorB.b, 0.5);
        }

        // Apply subtle color tint for underwater atmosphere
        float3 underwaterTint = float3(0.7, 0.85, 1.0);  // Slight blue-green tint
        color *= underwaterTint;
    }
#endif

    // ---------------------------
    // Chromatic Aberration (general)
    // ---------------------------
#if ENABLE_CHROMATIC_ABERRATION
    if (!isUnderwater)  // Skip if underwater CA was already applied
    {
        float2 center = uv - 0.5;
        float dist = length(center);
        float2 dir = center / max(dist, 0.001);

        // This is a placeholder - actual implementation would sample texture
        // at offset positions for R and B channels
        float aberration = dist * dist * CHROMATIC_ABERRATION_STRENGTH;
        // color.r = sampleTexture(uv + dir * aberration).r;
        // color.b = sampleTexture(uv - dir * aberration).b;
    }
#endif

    // ---------------------------
    // Lens Flare
    // ---------------------------
#if ENABLE_LENS_FLARE
    {
        // Calculate sun screen position
        float4 sunClipPos = mul(float4(g_view.directionToSun * 1000.0 + g_view.viewOriginSteveSpace, 1.0), g_view.viewProj);
        float2 sunScreenPos = sunClipPos.xy / sunClipPos.w * 0.5 + 0.5;
        sunScreenPos.y = 1.0 - sunScreenPos.y;

        // Only add lens flare when sun is visible and not underwater
        if (sunClipPos.w > 0.0 && all(saturate(sunScreenPos) == sunScreenPos) && !isUnderwater)
        {
            float3 sunColor = blackbodyColor(SUN_COLOR_TEMPERATURE);
            float3 lensFlare = computeLensFlare(uv, sunScreenPos, sunColor, LENS_FLARE_INTENSITY);
            color += lensFlare;
        }
    }
#endif

    // ---------------------------
    // Dithering (reduce banding)
    // ---------------------------
    color = applyDither(color, float2(pixelCoord), 1.0);

#else
    // Simple path without post-processing
    color = tonemap(color);
    color = linearToSRGB(color);
#endif

    // ---------------------------
    // Debug Views
    // ---------------------------
#if DEBUG_VIEW > 0
    {
        // Read surface data for debug (would need actual gbuffer in production)
        EnhancedSurface debugSurface;
        debugSurface.albedo = colorAlpha.rgb;
        debugSurface.normal = float3(0, 1, 0);
        debugSurface.roughness = 0.5;
        debugSurface.metalness = 0.0;
        debugSurface.emissive = 0.0;
        debugSurface.ao = 1.0;

        color = renderDebugView(debugSurface, depth, motionVector, ctx);
    }
#endif

    // ---------------------------
    // NaN/Inf Check
    // ---------------------------
#if DEBUG_NAN_INF
    if (any(isnan(color)) || any(isinf(color)))
    {
        int pattern = (pixelCoord.x / 16 + pixelCoord.y / 16) & 1;
        color = pattern ? float3(1, 0, 1) : float3(0, 0, 0);
    }
#endif

    // Final output
    outputBufferFinal[pixelCoord] = float4(saturate(color), colorAlpha.a);
}
