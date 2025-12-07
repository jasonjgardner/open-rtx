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
    // Film Grain
    // ---------------------------
#if ENABLE_FILM_GRAIN
    color = applyFilmGrain(color, uv, g_view.time, FILM_GRAIN_INTENSITY);
#endif

    // ---------------------------
    // Tonemapping
    // ---------------------------
    color = applyExposure(color, ctx.exposureEV);
    color = tonemap(color);

    // Convert to gamma space for further processing
    color = linearToSRGB(color);

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
    // Chromatic Aberration
    // ---------------------------
#if ENABLE_CHROMATIC_ABERRATION
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

        // Only add lens flare when sun is visible
        if (sunClipPos.w > 0.0 && all(saturate(sunScreenPos) == sunScreenPos))
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
