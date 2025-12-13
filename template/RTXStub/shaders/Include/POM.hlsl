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

#ifndef __POM_HLSL__
#define __POM_HLSL__

#include "Settings.hlsl"

// Small epsilon for comparisons
static const float kPomEpsilon = 1e-5;

// =============================================================================
// POM Helper Structures
// =============================================================================

struct POMResult
{
    float2 uv;           // Final displaced UV coordinates
    float depth;         // Depth at hit point (0 = surface, 1 = deepest)
    float3 normal;       // Computed normal (tangent space)
    float shadow;        // Self-shadow factor (0 = shadowed, 1 = lit)
};

struct POMContext
{
    // Texture and sampler
    Texture2D heightTex;
    SamplerState heightSampler;

    // Atlas tile bounds for UV clamping (prevents bleeding)
    float2 tileMinUV;
    float2 tileMaxUV;
    float2 texelSize;    // 1.0 / textureSize

    // View direction in tangent space (must be normalized)
    float3 viewDirTS;

    // Light direction in tangent space (for shadows)
    float3 lightDirTS;

    // Base UV coordinates
    float2 baseUV;

    // View distance for LOD/fade
    float viewDistance;
};

// =============================================================================
// UV Clamping - Prevents atlas tile bleeding
// =============================================================================

float2 clampUVToTile(float2 uv, float2 tileMin, float2 tileMax, float2 texelSize)
{
    // Inset by half texel to prevent sampling outside tile bounds
    float2 safeMin = tileMin + texelSize * 0.5;
    float2 safeMax = tileMax - texelSize * 0.5;
    return clamp(uv, safeMin, safeMax);
}

// =============================================================================
// Height Sampling with Atlas Clamping
// =============================================================================

float sampleHeight(Texture2D tex, SamplerState samp, float2 uv, float2 tileMin, float2 tileMax, float2 texelSize)
{
    float2 clampedUV = clampUVToTile(uv, tileMin, tileMax, texelSize);
    return tex.SampleLevel(samp, clampedUV, 0).a;
}

// =============================================================================
// Calculate Parallax Offset Direction
// =============================================================================

float2 getParallaxOffset(float3 viewDirTS, float parallaxDepth)
{
    // Project view direction onto the surface plane and scale by depth
    // viewDirTS.z is the component perpendicular to the surface
    float2 parallaxDir = viewDirTS.xy / max(abs(viewDirTS.z), kPomEpsilon);
    return parallaxDir * parallaxDepth;
}

// =============================================================================
// Main POM Ray Marching
// =============================================================================

POMResult computePOM(
    Texture2D heightTex,
    SamplerState heightSampler,
    float2 baseUV,
    float3 viewDirTS,
    float2 tileMinUV,
    float2 tileMaxUV,
    float2 texelSize,
    float viewDistance
)
{
    POMResult result;
    result.uv = baseUV;
    result.depth = 0.0;
    result.normal = float3(0, 0, 1);
    result.shadow = 1.0;

#if ENABLE_POM
    // Calculate distance fade factor (runtime check since POM_FADE_DISTANCE is float)
    float fadeFactor = 1.0;
    if (POM_FADE_DISTANCE > 0.0)
    {
        fadeFactor = saturate(1.0 - (viewDistance - POM_FADE_DISTANCE) / max(POM_FADE_RANGE, 0.001));
        if (fadeFactor <= 0.0)
        {
            return result;
        }
    }

    // Early out for flat surfaces (height near 1.0 means no displacement)
    float initialHeight = sampleHeight(heightTex, heightSampler, baseUV, tileMinUV, tileMaxUV, texelSize);
    if (initialHeight > 1.0 - kPomEpsilon)
    {
        return result;
    }

    // Adjust sample count based on view angle (more samples at grazing angles)
    float NoV = abs(viewDirTS.z);
    int numSamples = (int)lerp((float)POM_SAMPLES_MAX, (float)POM_SAMPLES_MIN, NoV);
    numSamples = max(numSamples, 4);

    // Calculate step size and UV offset per step
    float stepSize = 1.0 / (float)numSamples;
    float2 uvOffset = getParallaxOffset(viewDirTS, POM_DEPTH * fadeFactor);
    float2 uvStep = uvOffset * stepSize;

    // Ray march from surface (depth=0) downward (depth=1)
    float2 currentUV = baseUV;
    float currentDepth = 0.0;
    float currentHeight = initialHeight;

    float2 prevUV = baseUV;
    float prevDepth = 0.0;
    float prevHeight = 1.0;

    // Step through the heightfield
    [loop]
    for (int i = 0; i < numSamples; ++i)
    {
        currentHeight = sampleHeight(heightTex, heightSampler, currentUV, tileMinUV, tileMaxUV, texelSize);

        // Check if ray has gone below the heightfield surface
        // Height of 1.0 = top surface, 0.0 = bottom
        // We compare ray depth against (1.0 - height) to get intersection
        if (currentDepth >= (1.0 - currentHeight))
        {
            break;
        }

        prevUV = currentUV;
        prevDepth = currentDepth;
        prevHeight = currentHeight;

        currentUV -= uvStep;
        currentDepth += stepSize;
    }

    // Linear interpolation for more accurate intersection
#if POM_LINEAR_SEARCH
    float prevRayDepth = prevDepth;
    float prevSurfaceDepth = 1.0 - prevHeight;
    float currentRayDepth = currentDepth;
    float currentSurfaceDepth = 1.0 - currentHeight;

    float delta1 = currentRayDepth - currentSurfaceDepth;
    float delta2 = prevSurfaceDepth - prevRayDepth;
    float t = delta1 / max(delta1 + delta2, kPomEpsilon);

    result.uv = lerp(currentUV, prevUV, t);
    result.depth = lerp(currentDepth, prevDepth, t);
#else
    result.uv = prevUV;
    result.depth = prevDepth;
#endif

    // Clamp final UV to tile bounds
    result.uv = clampUVToTile(result.uv, tileMinUV, tileMaxUV, texelSize);

#endif // ENABLE_POM

    return result;
}

// =============================================================================
// POM Self-Shadowing
// =============================================================================

float computePOMShadow(
    Texture2D heightTex,
    SamplerState heightSampler,
    float2 hitUV,
    float hitDepth,
    float3 lightDirTS,
    float2 tileMinUV,
    float2 tileMaxUV,
    float2 texelSize
)
{
#if ENABLE_POM && ENABLE_POM_SHADOWS
    // If light is below horizon in tangent space, fully shadowed
    if (lightDirTS.z <= 0.0)
    {
        return 0.0;
    }

    // Calculate how many steps we need to trace to the surface
    int numSamples = (int)((1.0 - hitDepth) * POM_SHADOW_SAMPLES);
    if (numSamples < 1)
    {
        return 1.0;
    }

    float stepSize = (1.0 - hitDepth) / (float)numSamples;
    float2 uvOffset = getParallaxOffset(lightDirTS, POM_DEPTH);
    float2 uvStep = uvOffset * stepSize;

    float2 currentUV = hitUV;
    float currentDepth = hitDepth;
    float maxShadow = 0.0;

    [loop]
    for (int i = 0; i < numSamples; ++i)
    {
        currentUV += uvStep;
        currentDepth -= stepSize;

        float height = sampleHeight(heightTex, heightSampler, currentUV, tileMinUV, tileMaxUV, texelSize);
        float surfaceDepth = 1.0 - height;

        // Check if we're below the surface (in shadow)
        float shadowAmount = (surfaceDepth - currentDepth) * (1.0 / POM_SHADOW_SOFTNESS);

        if (shadowAmount > 0.0)
        {
            // Soft shadow based on how much the surface occludes
            maxShadow = max(maxShadow, saturate(shadowAmount));
            if (maxShadow >= 1.0 - kPomEpsilon)
            {
                break;
            }
        }
    }

    return 1.0 - maxShadow;
#else
    return 1.0;
#endif
}

// =============================================================================
// Slope Normal Calculation (Sharp Minecraft-style edges)
// =============================================================================

float3 computeSlopeNormal(
    Texture2D heightTex,
    SamplerState heightSampler,
    float2 hitUV,
    float hitDepth,
    float2 stepDir,
    float2 tileMinUV,
    float2 tileMaxUV,
    float2 texelSize
)
{
#if ENABLE_POM && ENABLE_POM_SLOPE_NORMALS
    // Get pixel-aligned UV
    float2 texSize = 1.0 / texelSize;
    float2 pixelUV = floor(hitUV * texSize) * texelSize;
    float2 subPixel = hitUV - pixelUV - 0.5 * texelSize;

    float2 stepSign = sign(stepDir);

    // Sample neighboring heights
    float2 uvX = pixelUV + float2(texelSize.x * stepSign.x, 0.0);
    float heightX = sampleHeight(heightTex, heightSampler, uvX, tileMinUV, tileMaxUV, texelSize);
    bool hasX = hitDepth > (1.0 - heightX) && sign(subPixel.x) == stepSign.x;

    float2 uvY = pixelUV + float2(0.0, texelSize.y * stepSign.y);
    float heightY = sampleHeight(heightTex, heightSampler, uvY, tileMinUV, tileMaxUV, texelSize);
    bool hasY = hitDepth > (1.0 - heightY) && sign(subPixel.y) == stepSign.y;

    // Aspect ratio correction
    float2 aspectSubPixel = subPixel;
    aspectSubPixel.x *= texSize.x / texSize.y;

    // Choose normal based on which edge we're closer to
    if (abs(aspectSubPixel.x) < abs(aspectSubPixel.y))
    {
        if (hasY) return normalize(float3(0.0, stepSign.y, 0.5));
        if (hasX) return normalize(float3(stepSign.x, 0.0, 0.5));
    }
    else
    {
        if (hasX) return normalize(float3(stepSign.x, 0.0, 0.5));
        if (hasY) return normalize(float3(0.0, stepSign.y, 0.5));
    }

    // Fallback: compute normal from step direction
    float s = step(abs(stepDir.y), abs(stepDir.x));
    return normalize(float3(float2(1.0 - s, s) * stepSign, 0.5));
#else
    return float3(0, 0, 1);
#endif
}

// =============================================================================
// Normal from Height (for areas without slope normals)
// =============================================================================

float3 computeNormalFromHeight(
    Texture2D heightTex,
    SamplerState heightSampler,
    float2 uv,
    float2 tileMinUV,
    float2 tileMaxUV,
    float2 texelSize,
    float heightScale
)
{
    // Sample 4 neighboring heights
    float2 offset = texelSize;

    float hL = sampleHeight(heightTex, heightSampler, uv - float2(offset.x, 0), tileMinUV, tileMaxUV, texelSize);
    float hR = sampleHeight(heightTex, heightSampler, uv + float2(offset.x, 0), tileMinUV, tileMaxUV, texelSize);
    float hD = sampleHeight(heightTex, heightSampler, uv - float2(0, offset.y), tileMinUV, tileMaxUV, texelSize);
    float hU = sampleHeight(heightTex, heightSampler, uv + float2(0, offset.y), tileMinUV, tileMaxUV, texelSize);

    // Compute normal from height differences
    float3 normal;
    normal.x = (hL - hR) * heightScale;
    normal.y = (hD - hU) * heightScale;
    normal.z = 1.0;

    return normalize(normal);
}

// =============================================================================
// Full POM Pipeline (convenience function)
// =============================================================================

POMResult computeFullPOM(
    Texture2D heightTex,
    SamplerState heightSampler,
    float2 baseUV,
    float3 viewDirTS,
    float3 lightDirTS,
    float2 tileMinUV,
    float2 tileMaxUV,
    float2 texelSize,
    float viewDistance
)
{
    // Compute displaced UV and depth
    POMResult result = computePOM(
        heightTex, heightSampler,
        baseUV, viewDirTS,
        tileMinUV, tileMaxUV, texelSize,
        viewDistance
    );

#if ENABLE_POM
    // Compute normal from heightfield
    float2 stepDir = -viewDirTS.xy;

#if ENABLE_POM_SLOPE_NORMALS
    result.normal = computeSlopeNormal(
        heightTex, heightSampler,
        result.uv, result.depth, stepDir,
        tileMinUV, tileMaxUV, texelSize
    );
#else
    result.normal = computeNormalFromHeight(
        heightTex, heightSampler,
        result.uv,
        tileMinUV, tileMaxUV, texelSize,
        POM_DEPTH * 10.0
    );
#endif

    // Compute self-shadowing
    result.shadow = computePOMShadow(
        heightTex, heightSampler,
        result.uv, result.depth, lightDirTS,
        tileMinUV, tileMaxUV, texelSize
    );
#endif

    return result;
}

#endif // __POM_HLSL__
