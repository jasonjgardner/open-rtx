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
// OpenRTX Tone Mapping
// Implements multiple tonemapping operators with color grading support
// =============================================================================

#ifndef __OPENRTX_TONEMAPPING_HLSL__
#define __OPENRTX_TONEMAPPING_HLSL__

#include "Settings.hlsl"

// =============================================================================
// COLOR SPACE CONVERSIONS
// =============================================================================

// sRGB to Linear
float3 sRGBToLinear(float3 srgb)
{
    return pow(max(srgb, 0.0), 2.2);
}

// Linear to sRGB
float3 linearToSRGB(float3 linear)
{
    return pow(max(linear, 0.0), 1.0 / 2.2);
}

// Luminance calculation (Rec. 709)
float luminance(float3 color)
{
    return dot(color, float3(0.2126, 0.7152, 0.0722));
}

// Luminance calculation (Rec. 2020 for HDR)
float luminanceRec2020(float3 color)
{
    return dot(color, float3(0.2627, 0.6780, 0.0593));
}

// RGB to XYZ (D65 white point)
float3 RGBtoXYZ(float3 rgb)
{
    float3x3 mat = float3x3(
        0.4124564, 0.3575761, 0.1804375,
        0.2126729, 0.7151522, 0.0721750,
        0.0193339, 0.1191920, 0.9503041);
    return mul(mat, rgb);
}

// XYZ to RGB (D65 white point)
float3 XYZtoRGB(float3 xyz)
{
    float3x3 mat = float3x3(
        3.2404542, -1.5371385, -0.4985314,
        -0.9692660, 1.8760108, 0.0415560,
        0.0556434, -0.2040259, 1.0572252);
    return mul(mat, xyz);
}

// RGB to HSV
float3 RGBtoHSV(float3 rgb)
{
    float4 K = float4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    float4 p = lerp(float4(rgb.bg, K.wz), float4(rgb.gb, K.xy), step(rgb.b, rgb.g));
    float4 q = lerp(float4(p.xyw, rgb.r), float4(rgb.r, p.yzx), step(p.x, rgb.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return float3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

// HSV to RGB
float3 HSVtoRGB(float3 hsv)
{
    float4 K = float4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    float3 p = abs(frac(hsv.xxx + K.xyz) * 6.0 - K.www);
    return hsv.z * lerp(K.xxx, saturate(p - K.xxx), hsv.y);
}

// =============================================================================
// EXPOSURE
// =============================================================================

// EV to exposure multiplier
float EVToExposure(float ev)
{
    return pow(2.0, ev);
}

// Apply exposure
float3 applyExposure(float3 color, float exposureEV)
{
    return color * EVToExposure(exposureEV);
}

// Auto exposure (simplified histogram-based)
float computeAutoExposure(float avgLuminance, float minEV, float maxEV)
{
    // Target middle gray
    const float targetLuminance = 0.18;

    // Compute required EV adjustment
    float ev = log2(max(avgLuminance, 0.001) / targetLuminance);
    ev = clamp(-ev, minEV, maxEV);

    return ev;
}

// =============================================================================
// VANILLA FILMIC TONEMAPPING
// =============================================================================

float3 tonemapVanillaFilmic(float3 color)
{
    // Vanilla Minecraft RTX-style filmic curve
    float3 x = max(0.0, color - 0.004);
    return (x * (6.2 * x + 0.5)) / (x * (6.2 * x + 1.7) + 0.06);
}

// =============================================================================
// ACES TONEMAPPING
// =============================================================================

// ACES input transform (approximate)
float3 ACESInputMat(float3 color)
{
    float3x3 mat = float3x3(
        0.59719, 0.35458, 0.04823,
        0.07600, 0.90834, 0.01566,
        0.02840, 0.13383, 0.83777);
    return mul(mat, color);
}

// ACES output transform (approximate)
float3 ACESOutputMat(float3 color)
{
    float3x3 mat = float3x3(
        1.60475, -0.53108, -0.07367,
        -0.10208, 1.10813, -0.00605,
        -0.00327, -0.07276, 1.07602);
    return mul(mat, color);
}

// RRT and ODT fit
float3 RRTAndODTFit(float3 v)
{
    float3 a = v * (v + 0.0245786) - 0.000090537;
    float3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return a / b;
}

float3 tonemapACES(float3 color)
{
    color *= ACES_INPUT_SCALE;
    color = ACESInputMat(color);
    color = RRTAndODTFit(color);
    color = ACESOutputMat(color);
    color *= ACES_OUTPUT_SCALE;
    return saturate(color);
}

// Simplified ACES (faster)
float3 tonemapACESFast(float3 color)
{
    color *= 0.6;
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    return saturate((color * (a * color + b)) / (color * (c * color + d) + e));
}

// =============================================================================
// UNCHARTED 2 (HABLE) TONEMAPPING
// =============================================================================

float3 hablePartial(float3 x)
{
    float A = 0.15;
    float B = 0.50;
    float C = 0.10;
    float D = 0.20;
    float E = 0.02;
    float F = 0.30;
    return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

float3 tonemapUncharted2(float3 color)
{
    float exposureBias = 2.0;
    float3 curr = hablePartial(color * exposureBias);

    float W = 11.2;
    float3 whiteScale = 1.0 / hablePartial(W);
    return curr * whiteScale;
}

// =============================================================================
// REINHARD TONEMAPPING
// =============================================================================

// Simple Reinhard
float3 tonemapReinhardSimple(float3 color)
{
    return color / (1.0 + color);
}

// Extended Reinhard with white point
float3 tonemapReinhardExtended(float3 color, float whitePoint)
{
    float3 numerator = color * (1.0 + color / (whitePoint * whitePoint));
    return numerator / (1.0 + color);
}

// Luminance-based Reinhard
float3 tonemapReinhardLuminance(float3 color)
{
    float lum = luminance(color);
    float toneMappedLum = lum / (1.0 + lum);
    return color * (toneMappedLum / max(lum, 0.0001));
}

// =============================================================================
// UCHIMURA TONEMAPPING
// =============================================================================

float3 tonemapUchimura(float3 x)
{
    // Attempt to match GT (Gran Turismo) tonemapping
    const float P = 1.0;  // Max brightness
    const float a = 1.0;  // Contrast
    const float m = 0.22; // Linear section start
    const float l = 0.4;  // Linear section length
    const float c = 1.33; // Black tightness
    const float b = 0.0;  // Pedestal

    float l0 = ((P - m) * l) / a;
    float L0 = m - m / a;
    float L1 = m + (1.0 - m) / a;
    float S0 = m + l0;
    float S1 = m + a * l0;
    float C2 = (a * P) / (P - S1);
    float CP = -C2 / P;

    float3 w0 = 1.0 - smoothstep(0.0, m, x);
    float3 w2 = step(m + l0, x);
    float3 w1 = 1.0 - w0 - w2;

    float3 T = m * pow(x / m, c) + b;
    float3 S = P - (P - S1) * exp(CP * (x - S0));
    float3 L = m + a * (x - m);

    return T * w0 + L * w1 + S * w2;
}

// =============================================================================
// AGX TONEMAPPING
// =============================================================================

// AgX base contrast
float3 agxDefaultContrastApprox(float3 x)
{
    float3 x2 = x * x;
    float3 x4 = x2 * x2;

    return +15.5 * x4 * x2
           - 40.14 * x4 * x
           + 31.96 * x4
           - 6.868 * x2 * x
           + 0.4298 * x2
           + 0.1191 * x
           - 0.00232;
}

float3 tonemapAgX(float3 color)
{
    // AgX transform matrix
    float3x3 agxMat = float3x3(
        0.842479062253094, 0.0423282422610123, 0.0423756549057051,
        0.0784335999999992, 0.878468636469772, 0.0784336,
        0.0792237451477643, 0.0791661274605434, 0.879142973793104);

    float3x3 agxMatInv = float3x3(
        1.19687900512017, -0.0528968517574562, -0.0529716355144438,
        -0.0980208811401368, 1.15190312990417, -0.0980434501171241,
        -0.0990297440797205, -0.0989611768448433, 1.15107367264116);

    // Apply AgX Log encoding
    float minEv = -12.47393;
    float maxEv = 4.026069;
    color = mul(agxMat, color);
    color = clamp(log2(color), minEv, maxEv);
    color = (color - minEv) / (maxEv - minEv);

    // Apply sigmoid
    color = agxDefaultContrastApprox(color);

    // Convert back
    color = mul(agxMatInv, color);

    return color;
}

// =============================================================================
// MAIN TONEMAPPING FUNCTION
// =============================================================================

float3 tonemap(float3 color)
{
#if TONEMAPPING_TYPE == 0
    return tonemapVanillaFilmic(color);
#elif TONEMAPPING_TYPE == 1
    return tonemapACES(color);
#elif TONEMAPPING_TYPE == 2
    return tonemapUncharted2(color);
#elif TONEMAPPING_TYPE == 3
    return tonemapReinhardLuminance(color);
#elif TONEMAPPING_TYPE == 4
    return tonemapUchimura(color);
#elif TONEMAPPING_TYPE == 5
    return tonemapAgX(color);
#else
    return saturate(color); // No tonemapping
#endif
}

// =============================================================================
// COLOR GRADING
// =============================================================================

// Contrast adjustment
float3 adjustContrast(float3 color, float contrast)
{
    return (color - 0.5) * contrast + 0.5;
}

// Saturation adjustment
float3 adjustSaturation(float3 color, float saturation)
{
    float lum = luminance(color);
    return lerp(lum, color, saturation);
}

// Brightness adjustment
float3 adjustBrightness(float3 color, float brightness)
{
    return color + brightness;
}

// Lift/Gamma/Gain color correction
float3 liftGammaGain(float3 color, float3 lift, float3 gamma, float3 gain)
{
    color = color * gain + lift;
    color = sign(color) * pow(abs(color), 1.0 / gamma);
    return color;
}

// Color temperature adjustment (simplified)
float3 adjustColorTemperature(float3 color, float temperature)
{
    // Shift towards warm (positive) or cool (negative)
    float3 warm = float3(1.0, 0.9, 0.8);
    float3 cool = float3(0.8, 0.9, 1.0);
    float3 tint = lerp(cool, warm, temperature * 0.5 + 0.5);
    return color * tint;
}

// Vignette
float3 applyVignette(float3 color, float2 uv, float intensity, float radius, float softness)
{
#if ENABLE_VIGNETTE
    float2 center = uv - 0.5;
    float dist = length(center);
    float vignette = smoothstep(radius, radius - softness, dist);
    return color * lerp(1.0 - intensity, 1.0, vignette);
#else
    return color;
#endif
}

// =============================================================================
// COMPLETE POST-PROCESSING PIPELINE
// =============================================================================

struct PostProcessParams
{
    float exposureEV;
    float contrast;
    float saturation;
    float brightness;
    float3 lift;
    float3 gamma;
    float3 gain;
    float temperature;
    float vignetteIntensity;
    float vignetteRadius;
    float vignetteSoftness;
};

float3 postProcess(float3 color, float2 uv, PostProcessParams params)
{
    // Apply exposure
    color = applyExposure(color, params.exposureEV);

    // Tonemapping
    color = tonemap(color);

    // Color grading (in gamma space)
    color = linearToSRGB(color);

    color = liftGammaGain(color, params.lift, params.gamma, params.gain);
    color = adjustContrast(color, params.contrast);
    color = adjustSaturation(color, params.saturation);
    color = adjustBrightness(color, params.brightness);
    color = adjustColorTemperature(color, params.temperature);

    // Vignette
    color = applyVignette(color, uv, params.vignetteIntensity, params.vignetteRadius, params.vignetteSoftness);

    // Ensure valid output
    return saturate(color);
}

// Simplified post-processing with defaults
float3 postProcessSimple(float3 color, float2 uv, float exposureEV)
{
    PostProcessParams params;
    params.exposureEV = exposureEV;
    params.contrast = 1.0;
    params.saturation = 1.0;
    params.brightness = 0.0;
    params.lift = 0.0;
    params.gamma = 1.0;
    params.gain = 1.0;
    params.temperature = 0.0;
    params.vignetteIntensity = VIGNETTE_INTENSITY;
    params.vignetteRadius = VIGNETTE_RADIUS;
    params.vignetteSoftness = VIGNETTE_SOFTNESS;

    return postProcess(color, uv, params);
}

#endif // __OPENRTX_TONEMAPPING_HLSL__
