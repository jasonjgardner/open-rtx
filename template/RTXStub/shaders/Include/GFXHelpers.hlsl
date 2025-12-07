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
// OpenRTX Graphics Helpers
// Common graphics utility functions and post-processing effects
// =============================================================================

#ifndef __OPENRTX_GFXHELPERS_HLSL__
#define __OPENRTX_GFXHELPERS_HLSL__

#include "Settings.hlsl"

// =============================================================================
// RANDOM NUMBER GENERATION
// =============================================================================

// PCG random number generator
uint pcgHash(uint input)
{
    uint state = input * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

// Float random [0, 1)
float randomFloat(uint seed)
{
    return float(pcgHash(seed)) / 4294967296.0;
}

// Float2 random [0, 1)^2
float2 randomFloat2(uint seed)
{
    return float2(randomFloat(seed), randomFloat(seed + 1u));
}

// Float3 random [0, 1)^3
float3 randomFloat3(uint seed)
{
    return float3(randomFloat(seed), randomFloat(seed + 1u), randomFloat(seed + 2u));
}

// Blue noise hash
float blueNoise(float2 coord, float time)
{
    float2 p = coord + time * 120.0;
    return frac(sin(dot(p, float2(12.9898, 78.233))) * 43758.5453);
}

// =============================================================================
// SAMPLING FUNCTIONS
// =============================================================================

// Cosine-weighted hemisphere sampling
float3 sampleCosineHemisphere(float2 u, float3 normal)
{
    float phi = kTwoPi * u.x;
    float cosTheta = sqrt(u.y);
    float sinTheta = sqrt(1.0 - u.y);

    float3 tangent, bitangent;
    if (abs(normal.y) < 0.999)
    {
        tangent = normalize(cross(float3(0, 1, 0), normal));
    }
    else
    {
        tangent = normalize(cross(float3(1, 0, 0), normal));
    }
    bitangent = cross(normal, tangent);

    float3 localDir = float3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
    return tangent * localDir.x + bitangent * localDir.y + normal * localDir.z;
}

// Uniform hemisphere sampling
float3 sampleUniformHemisphere(float2 u, float3 normal)
{
    float phi = kTwoPi * u.x;
    float cosTheta = u.y;
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    float3 tangent, bitangent;
    if (abs(normal.y) < 0.999)
    {
        tangent = normalize(cross(float3(0, 1, 0), normal));
    }
    else
    {
        tangent = normalize(cross(float3(1, 0, 0), normal));
    }
    bitangent = cross(normal, tangent);

    float3 localDir = float3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
    return tangent * localDir.x + bitangent * localDir.y + normal * localDir.z;
}

// GGX importance sampling
float3 sampleGGX(float2 u, float3 normal, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;

    float phi = kTwoPi * u.x;
    float cosTheta = sqrt((1.0 - u.y) / (1.0 + (a2 - 1.0) * u.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    float3 H = float3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    float3 up = abs(normal.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, normal));
    float3 bitangent = cross(normal, tangent);

    return normalize(tangent * H.x + bitangent * H.y + normal * H.z);
}

// =============================================================================
// POST-PROCESSING EFFECTS
// =============================================================================

// Film grain
float3 applyFilmGrain(float3 color, float2 uv, float time, float intensity)
{
#if ENABLE_FILM_GRAIN
    float grain = blueNoise(uv * 500.0, time) - 0.5;

    // Luminance-responsive grain (less grain in bright areas)
    float lum = dot(color, float3(0.2126, 0.7152, 0.0722));
    float grainResponse = lerp(1.0, 0.0, saturate(lum) * FILM_GRAIN_LUMINANCE_RESPONSE);

    return color + grain * intensity * grainResponse;
#else
    return color;
#endif
}

// Chromatic aberration
float3 applyChromaticAberration(float3 color, float2 uv, float strength)
{
#if ENABLE_CHROMATIC_ABERRATION
    float2 center = uv - 0.5;
    float dist = length(center);
    float2 dir = center / max(dist, 0.001);

    // Sample each channel at slightly different positions
    // Note: Would need actual texture sampling in real implementation
    float2 redOffset = uv + dir * dist * strength;
    float2 blueOffset = uv - dir * dist * strength;

    // Placeholder - actual implementation would sample texture
    return color; // Return unmodified for now
#else
    return color;
#endif
}

// Lens flare artifacts
float3 computeLensFlare(float2 uv, float2 lightPos, float3 lightColor, float intensity)
{
#if ENABLE_LENS_FLARE
    float3 flare = 0.0;

    // Ghost artifacts
    int numGhosts = 4;
    for (int i = 0; i < numGhosts; i++)
    {
        float scale = 1.0 - float(i) * 0.15;
        float2 ghostPos = lightPos + (uv - lightPos) * (float(i + 1) * 0.5);

        float dist = length(ghostPos - uv);
        float ghost = smoothstep(0.1 * scale, 0.0, dist);

        // Chromatic shift
        float3 chromaShift = float3(
            smoothstep(0.08 * scale, 0.0, dist - 0.01),
            smoothstep(0.08 * scale, 0.0, dist),
            smoothstep(0.08 * scale, 0.0, dist + 0.01));

        flare += ghost * chromaShift * lightColor * (0.5 / float(i + 1));
    }

    // Halo
    float2 toLight = lightPos - uv;
    float dist = length(toLight);
    float halo = smoothstep(0.4, 0.3, dist) * smoothstep(0.2, 0.3, dist);
    flare += halo * lightColor * 0.3;

    return flare * intensity;
#else
    return 0.0;
#endif
}

// =============================================================================
// DEPTH OF FIELD
// =============================================================================

struct DOFParams
{
    float focusDistance;
    float aperture;
    float nearTransition;
    float farTransition;
};

// Calculate Circle of Confusion
float calculateCoC(float depth, DOFParams params)
{
    float coc = params.aperture * abs(depth - params.focusDistance) / depth;
    return coc;
}

// =============================================================================
// MOTION BLUR
// =============================================================================

float3 applyMotionBlur(float3 color, float2 uv, float2 velocity, int numSamples)
{
#if ENABLE_MOTION_BLUR
    if (length(velocity) < 0.001)
        return color;

    float3 result = color;
    float2 step = velocity * MOTION_BLUR_STRENGTH / float(numSamples);

    // Note: Would need actual texture sampling in real implementation
    for (int i = 1; i < numSamples; i++)
    {
        float2 sampleUV = uv + step * float(i);
        // result += textureSample(sampleUV);
    }

    return result / float(numSamples);
#else
    return color;
#endif
}

// =============================================================================
// BLOOM
// =============================================================================

// Extract bright pixels for bloom
float3 extractBloomBrightness(float3 color, float threshold)
{
#if ENABLE_BLOOM
    float brightness = max(max(color.r, color.g), color.b);
    float contribution = max(0.0, brightness - threshold);
    contribution /= max(brightness, 0.001);
    return color * contribution;
#else
    return 0.0;
#endif
}

// =============================================================================
// ANTI-ALIASING HELPERS
// =============================================================================

// Edge detection for FXAA
float2 computeGradient(float3 colorN, float3 colorS, float3 colorE, float3 colorW)
{
    float lumaN = dot(colorN, float3(0.299, 0.587, 0.114));
    float lumaS = dot(colorS, float3(0.299, 0.587, 0.114));
    float lumaE = dot(colorE, float3(0.299, 0.587, 0.114));
    float lumaW = dot(colorW, float3(0.299, 0.587, 0.114));

    return float2(lumaE - lumaW, lumaN - lumaS);
}

// =============================================================================
// DITHERING
// =============================================================================

// Bayer dithering pattern
float bayerDither(float2 coord)
{
    int2 c = int2(coord) % 8;
    int index = c.x + c.y * 8;

    // 8x8 Bayer matrix values (normalized)
    static const float bayer[64] = {
        0.0/64.0, 32.0/64.0,  8.0/64.0, 40.0/64.0,  2.0/64.0, 34.0/64.0, 10.0/64.0, 42.0/64.0,
       48.0/64.0, 16.0/64.0, 56.0/64.0, 24.0/64.0, 50.0/64.0, 18.0/64.0, 58.0/64.0, 26.0/64.0,
       12.0/64.0, 44.0/64.0,  4.0/64.0, 36.0/64.0, 14.0/64.0, 46.0/64.0,  6.0/64.0, 38.0/64.0,
       60.0/64.0, 28.0/64.0, 52.0/64.0, 20.0/64.0, 62.0/64.0, 30.0/64.0, 54.0/64.0, 22.0/64.0,
        3.0/64.0, 35.0/64.0, 11.0/64.0, 43.0/64.0,  1.0/64.0, 33.0/64.0,  9.0/64.0, 41.0/64.0,
       51.0/64.0, 19.0/64.0, 59.0/64.0, 27.0/64.0, 49.0/64.0, 17.0/64.0, 57.0/64.0, 25.0/64.0,
       15.0/64.0, 47.0/64.0,  7.0/64.0, 39.0/64.0, 13.0/64.0, 45.0/64.0,  5.0/64.0, 37.0/64.0,
       63.0/64.0, 31.0/64.0, 55.0/64.0, 23.0/64.0, 61.0/64.0, 29.0/64.0, 53.0/64.0, 21.0/64.0
    };

    return bayer[index];
}

// Apply dithering to reduce banding
float3 applyDither(float3 color, float2 coord, float strength)
{
    float dither = bayerDither(coord) - 0.5;
    return color + dither * strength / 255.0;
}

// =============================================================================
// DEBUG VISUALIZATION
// =============================================================================

// Heat map visualization
float3 heatMap(float value)
{
    float3 color;

    if (value < 0.25)
        color = lerp(float3(0, 0, 0), float3(0, 0, 1), value * 4.0);
    else if (value < 0.5)
        color = lerp(float3(0, 0, 1), float3(0, 1, 0), (value - 0.25) * 4.0);
    else if (value < 0.75)
        color = lerp(float3(0, 1, 0), float3(1, 1, 0), (value - 0.5) * 4.0);
    else
        color = lerp(float3(1, 1, 0), float3(1, 0, 0), (value - 0.75) * 4.0);

    return color;
}

// Normal visualization
float3 visualizeNormal(float3 normal)
{
    return normal * 0.5 + 0.5;
}

// Depth visualization
float3 visualizeDepth(float depth, float nearPlane, float farPlane)
{
    float linearDepth = (depth - nearPlane) / (farPlane - nearPlane);
    return linearDepth;
}

#endif // __OPENRTX_GFXHELPERS_HLSL__
