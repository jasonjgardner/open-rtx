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
// AccumulateGIInscatter - Temporal Accumulation for GI Volumetrics
// Similar to AccumulateInscatter but for the lower-resolution GI grid
// Uses more aggressive temporal filtering since GI changes slowly
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/Settings.hlsl"

// =============================================================================
// GI TEMPORAL ACCUMULATION PARAMETERS
// =============================================================================

// GI froxel grid dimensions (must match CalculateGIInscatterInline.hlsl)
static const uint3 kGIFroxelDimensions = uint3(128, 64, 32);

// GI temporal blend factor (lower than direct light for more stability)
#ifndef GI_INSCATTER_TEMPORAL_BLEND
#define GI_INSCATTER_TEMPORAL_BLEND 0.05
#endif

// Maximum blend (when history is invalid)
#ifndef GI_INSCATTER_TEMPORAL_BLEND_MAX
#define GI_INSCATTER_TEMPORAL_BLEND_MAX 0.5
#endif

// History rejection threshold (more lenient for GI)
#ifndef GI_HISTORY_REJECTION_THRESHOLD
#define GI_HISTORY_REJECTION_THRESHOLD 0.8
#endif

// GI froxel depth parameters (must match CalculateGIInscatterInline)
#ifndef GI_FROXEL_MAX_DEPTH
#define GI_FROXEL_MAX_DEPTH 192.0
#endif

#ifndef GI_FROXEL_DEPTH_EXPONENT
#define GI_FROXEL_DEPTH_EXPONENT 2.5
#endif

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Convert GI froxel slice to depth
float giSliceToDepth(float slice, float maxSlices)
{
    float normalizedSlice = slice / maxSlices;
    return GI_FROXEL_MAX_DEPTH * pow(normalizedSlice, GI_FROXEL_DEPTH_EXPONENT);
}

// Get view-space position from GI froxel coordinate
float3 giFroxelToViewPos(uint3 froxelCoord, float2 fovTan)
{
    float2 uv = (float2(froxelCoord.xy) + 0.5) / float2(kGIFroxelDimensions.xy);
    float2 ndc = uv * 2.0 - 1.0;
    float depth = giSliceToDepth(float(froxelCoord.z) + 0.5, float(kGIFroxelDimensions.z));

    return float3(ndc.x * fovTan.x * depth, ndc.y * fovTan.y * depth, -depth);
}

// Simplified reprojection for GI (more tolerant of approximations)
float3 reprojectGI(float3 currentFroxel)
{
    // For GI, simple temporal reprojection is often sufficient
    // since GI changes slowly and doesn't need per-pixel accuracy
    return currentFroxel;
}

// Manual trilinear sample from 3D texture using Load operations
// This avoids dependency on a linear sampler
float4 sampleGITrilinear(float3 froxelCoord, Texture3D<float4> tex)
{
    // Clamp to valid range
    float3 clampedCoord = clamp(froxelCoord, float3(0.5, 0.5, 0.5), float3(kGIFroxelDimensions) - 0.5);

    // Get integer and fractional parts
    float3 floorCoord = floor(clampedCoord - 0.5);
    float3 frac3 = clampedCoord - 0.5 - floorCoord;

    // Clamp integer coordinates to valid range
    int3 coord0 = clamp(int3(floorCoord), int3(0, 0, 0), int3(kGIFroxelDimensions) - 1);
    int3 coord1 = clamp(coord0 + int3(1, 1, 1), int3(0, 0, 0), int3(kGIFroxelDimensions) - 1);

    // Sample 8 corners of the trilinear cell
    float4 c000 = tex.Load(int4(coord0.x, coord0.y, coord0.z, 0));
    float4 c100 = tex.Load(int4(coord1.x, coord0.y, coord0.z, 0));
    float4 c010 = tex.Load(int4(coord0.x, coord1.y, coord0.z, 0));
    float4 c110 = tex.Load(int4(coord1.x, coord1.y, coord0.z, 0));
    float4 c001 = tex.Load(int4(coord0.x, coord0.y, coord1.z, 0));
    float4 c101 = tex.Load(int4(coord1.x, coord0.y, coord1.z, 0));
    float4 c011 = tex.Load(int4(coord0.x, coord1.y, coord1.z, 0));
    float4 c111 = tex.Load(int4(coord1.x, coord1.y, coord1.z, 0));

    // Trilinear interpolation
    float4 c00 = lerp(c000, c100, frac3.x);
    float4 c10 = lerp(c010, c110, frac3.x);
    float4 c01 = lerp(c001, c101, frac3.x);
    float4 c11 = lerp(c011, c111, frac3.x);

    float4 c0 = lerp(c00, c10, frac3.y);
    float4 c1 = lerp(c01, c11, frac3.y);

    return lerp(c0, c1, frac3.z);
}

// Luminance-based clamping for GI (softer than direct light)
float4 softClamp(float4 value, float4 minVal, float4 maxVal, float softness)
{
    float4 center = (minVal + maxVal) * 0.5;
    float4 range = (maxVal - minVal) * 0.5 * (1.0 + softness);

    float4 clamped = clamp(value, center - range, center + range);
    return lerp(value, clamped, 0.8);  // Soft blend toward clamped
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.AccumulateGIInscatter.group_size}})]
void AccumulateGIInscatter(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID
)
{
    // Note that g_rootConstant0 from BlurGIInscatter pass is accessible here
    uint sourceBufferIndex = g_rootConstant0;

#if !ENABLE_VOLUMETRIC_LIGHTING
    // Pass through if disabled
    if (all(dispatchThreadID < kGIFroxelDimensions))
    {
        volumetricGIResolvedInscatterRW[dispatchThreadID] = volumetricGIInscatterRW[sourceBufferIndex][dispatchThreadID];
    }
    return;
#endif

    // Bounds check
    if (any(dispatchThreadID >= kGIFroxelDimensions))
        return;

    // Get current frame GI inscatter (from blurred buffer)
    float4 currentGI = volumetricGIInscatterRW[sourceBufferIndex][dispatchThreadID];

    // Calculate neighborhood bounds for clamping
    float4 neighborMin = currentGI;
    float4 neighborMax = currentGI;

    // Sample 6-neighborhood
    [unroll]
    for (int axis = 0; axis < 3; axis++)
    {
        int3 offsetNeg = int3(0, 0, 0);
        int3 offsetPos = int3(0, 0, 0);
        offsetNeg[axis] = -1;
        offsetPos[axis] = 1;

        int3 coordNeg = clamp(int3(dispatchThreadID) + offsetNeg, int3(0, 0, 0), int3(kGIFroxelDimensions) - 1);
        int3 coordPos = clamp(int3(dispatchThreadID) + offsetPos, int3(0, 0, 0), int3(kGIFroxelDimensions) - 1);

        float4 sampleNeg = volumetricGIInscatterRW[sourceBufferIndex][uint3(coordNeg)];
        float4 samplePos = volumetricGIInscatterRW[sourceBufferIndex][uint3(coordPos)];

        neighborMin = min(neighborMin, min(sampleNeg, samplePos));
        neighborMax = max(neighborMax, max(sampleNeg, samplePos));
    }

    // Expand neighborhood bounds (GI is smoother, allow more variance)
    float4 neighborRange = neighborMax - neighborMin;
    neighborMin -= neighborRange * 0.5;
    neighborMax += neighborRange * 0.5;

    // Sample history
    float3 prevFroxelCoord = reprojectGI(float3(dispatchThreadID));
    float4 historyGI = float4(0, 0, 0, 1);
    float blendFactor = GI_INSCATTER_TEMPORAL_BLEND;

    // Check validity
    bool historyValid = all(prevFroxelCoord >= 0) && all(prevFroxelCoord < float3(kGIFroxelDimensions));

    if (historyValid)
    {
        // Sample previous frame with trilinear filtering (using manual trilinear)
        historyGI = sampleGITrilinear(prevFroxelCoord, volumetricGIInscatterPrevious);

        // Soft clamp history (GI should be smooth, harsh clamping causes flicker)
        historyGI = softClamp(historyGI, neighborMin, neighborMax, 0.5);

        // Calculate difference for confidence
        float4 diff = abs(currentGI - historyGI);
        float maxDiff = max(max(diff.r, diff.g), max(diff.b, diff.a));

        // GI confidence (more lenient than direct light)
        float confidence = 1.0 - saturate(maxDiff / GI_HISTORY_REJECTION_THRESHOLD);

        // Smooth confidence curve
        confidence = smoothstep(0.0, 1.0, confidence);

        // Adjust blend factor
        blendFactor = lerp(GI_INSCATTER_TEMPORAL_BLEND_MAX, GI_INSCATTER_TEMPORAL_BLEND, confidence);

        // Additional stability: reduce blend when camera is moving fast
        // This helps prevent ghosting during rapid camera motion
        // (Would use motion vector magnitude if available)
    }
    else
    {
        // No valid history - but don't use 100% current for GI
        // Use a softer reset to prevent flicker
        blendFactor = GI_INSCATTER_TEMPORAL_BLEND_MAX;
    }

    // Temporal accumulation
    float4 accumulatedGI = lerp(historyGI, currentGI, blendFactor);

    // Firefly suppression (softer for GI)
    float luminance = dot(accumulatedGI.rgb, float3(0.299, 0.587, 0.114));
    float maxLuminance = 5.0;  // Lower threshold for GI
    if (luminance > maxLuminance)
    {
        accumulatedGI.rgb *= maxLuminance / luminance;
    }

    // Ensure non-negative
    accumulatedGI = max(accumulatedGI, 0.0);

    // Preserve transmittance
    accumulatedGI.a = lerp(historyGI.a, currentGI.a, max(blendFactor, 0.2));

    // Write result
    volumetricGIResolvedInscatterRW[dispatchThreadID] = accumulatedGI;
}
