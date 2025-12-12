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
// AccumulateInscatter - Temporal Accumulation for Direct Light Volumetrics
// Blends current frame inscatter with previous frame for temporal stability
// Includes motion-aware reprojection for reduced ghosting
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/Settings.hlsl"

// =============================================================================
// TEMPORAL ACCUMULATION PARAMETERS
// =============================================================================

// Froxel grid dimensions (must match CalculateInscatterInline.hlsl)
static const uint3 kFroxelDimensions = uint3(256, 128, 64);

// Temporal blend factor (lower = more accumulation, smoother but more ghosting)
#ifndef INSCATTER_TEMPORAL_BLEND
#define INSCATTER_TEMPORAL_BLEND 0.1
#endif

// Maximum temporal blend (used when history is invalid)
#ifndef INSCATTER_TEMPORAL_BLEND_MAX
#define INSCATTER_TEMPORAL_BLEND_MAX 1.0
#endif

// History rejection threshold (reject history if too different)
#ifndef INSCATTER_HISTORY_REJECTION_THRESHOLD
#define INSCATTER_HISTORY_REJECTION_THRESHOLD 0.5
#endif

// Froxel max depth (must match CalculateInscatterInline)
#ifndef FROXEL_MAX_DEPTH
#define FROXEL_MAX_DEPTH 256.0
#endif

#ifndef FROXEL_DEPTH_EXPONENT
#define FROXEL_DEPTH_EXPONENT 2.0
#endif

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Convert froxel coordinate to world position (for reprojection)
float froxelSliceToDepth(float slice, float maxSlices)
{
    float normalizedSlice = slice / maxSlices;
    return FROXEL_MAX_DEPTH * pow(normalizedSlice, FROXEL_DEPTH_EXPONENT);
}

// Get view-space position from froxel coordinate
float3 froxelToViewPos(uint3 froxelCoord, float2 fovTan)
{
    // Normalized froxel coordinates [0, 1]
    float2 uv = (float2(froxelCoord.xy) + 0.5) / float2(kFroxelDimensions.xy);

    // Convert to NDC [-1, 1]
    float2 ndc = uv * 2.0 - 1.0;

    // Get depth
    float depth = froxelSliceToDepth(float(froxelCoord.z) + 0.5, float(kFroxelDimensions.z));

    // View-space position (looking down -Z)
    return float3(ndc.x * fovTan.x * depth, ndc.y * fovTan.y * depth, -depth);
}

// Reproject view position to previous frame froxel coordinate
float3 reprojectToPreviousFroxel(float3 viewPos, float4x4 currentViewProj, float4x4 prevViewProjInv, float2 fovTan)
{
    // Transform to world space (current frame)
    // This is simplified - in production, use proper inverse matrices
    float4 clipPos = mul(float4(viewPos, 1.0), currentViewProj);
    float3 ndcPos = clipPos.xyz / clipPos.w;

    // For now, return current position (proper reprojection needs prev frame matrices)
    // This is a placeholder - real implementation would use motion vectors or prev matrices
    float2 prevUV = ndcPos.xy * 0.5 + 0.5;
    float prevDepth = -viewPos.z;

    // Convert back to froxel coordinates
    float3 prevFroxel;
    prevFroxel.xy = prevUV * float2(kFroxelDimensions.xy);
    prevFroxel.z = float(kFroxelDimensions.z) * pow(saturate(prevDepth / FROXEL_MAX_DEPTH), 1.0 / FROXEL_DEPTH_EXPONENT);

    return prevFroxel;
}

// Trilinear sample from 3D texture with bounds checking
float4 sampleInscatterTrilinear(float3 froxelCoord, Texture3D<float4> tex, SamplerState samp)
{
    float3 normalizedCoord = froxelCoord / float3(kFroxelDimensions);

    // Clamp to valid range
    normalizedCoord = saturate(normalizedCoord);

    return tex.SampleLevel(samp, normalizedCoord, 0);
}

// Color clamping for anti-ghosting (neighborhood clamp)
float4 clampToNeighborhood(float4 currentSample, float4 historySample, float4 neighborMin, float4 neighborMax)
{
    // Clamp history to current frame's neighborhood
    return clamp(historySample, neighborMin, neighborMax);
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.AccumulateInscatter.group_size}})]
void AccumulateInscatter(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID
)
{
    // Note that g_rootConstant0 from BlurGIInscatter pass is accessible here

#if !ENABLE_VOLUMETRIC_LIGHTING
    // Pass through if disabled
    if (all(dispatchThreadID < kFroxelDimensions))
    {
        volumetricResolvedInscatterRW[dispatchThreadID] = volumetricInscatterRW[dispatchThreadID];
    }
    return;
#endif

    // Bounds check
    if (any(dispatchThreadID >= kFroxelDimensions))
        return;

    // Get current frame inscatter
    float4 currentInscatter = volumetricInscatterRW[dispatchThreadID];

    // Calculate neighbor min/max for clamping (anti-ghosting)
    float4 neighborMin = currentInscatter;
    float4 neighborMax = currentInscatter;

    // Sample 6-neighborhood for variance estimation
    [unroll]
    for (int axis = 0; axis < 3; axis++)
    {
        int3 offsetNeg = int3(0, 0, 0);
        int3 offsetPos = int3(0, 0, 0);
        offsetNeg[axis] = -1;
        offsetPos[axis] = 1;

        int3 coordNeg = int3(dispatchThreadID) + offsetNeg;
        int3 coordPos = int3(dispatchThreadID) + offsetPos;

        // Clamp to bounds
        coordNeg = clamp(coordNeg, int3(0, 0, 0), int3(kFroxelDimensions) - 1);
        coordPos = clamp(coordPos, int3(0, 0, 0), int3(kFroxelDimensions) - 1);

        float4 sampleNeg = volumetricInscatterRW[uint3(coordNeg)];
        float4 samplePos = volumetricInscatterRW[uint3(coordPos)];

        neighborMin = min(neighborMin, min(sampleNeg, samplePos));
        neighborMax = max(neighborMax, max(sampleNeg, samplePos));
    }

    // Add some margin to neighborhood bounds
    float4 neighborRange = neighborMax - neighborMin;
    neighborMin -= neighborRange * 0.25;
    neighborMax += neighborRange * 0.25;

    // Get FOV for reprojection
    float2 fovTan = float2(1.0 / g_view.proj[0][0], 1.0 / g_view.proj[1][1]);

    // Get view-space position
    float3 viewPos = froxelToViewPos(dispatchThreadID, fovTan);

    // Reproject to previous frame (simplified - uses current frame as approximation)
    float3 prevFroxelCoord = reprojectToPreviousFroxel(viewPos, g_view.viewProj, g_view.viewProj, fovTan);

    // Sample history buffer
    float4 historyInscatter = float4(0, 0, 0, 1);
    float blendFactor = INSCATTER_TEMPORAL_BLEND;

    // Check if reprojected position is valid
    bool historyValid = all(prevFroxelCoord >= 0) && all(prevFroxelCoord < float3(kFroxelDimensions));

    if (historyValid)
    {
        // Sample previous frame with trilinear filtering
        float3 normalizedCoord = prevFroxelCoord / float3(kFroxelDimensions);
        historyInscatter = volumetricInscatterPrevious.SampleLevel(linearClampSampler, normalizedCoord, 0);

        // Clamp history to neighborhood (anti-ghosting)
        historyInscatter = clampToNeighborhood(currentInscatter, historyInscatter, neighborMin, neighborMax);

        // Calculate confidence based on how similar history is
        float4 diff = abs(currentInscatter - historyInscatter);
        float maxDiff = max(max(diff.r, diff.g), max(diff.b, diff.a));
        float confidence = 1.0 - saturate(maxDiff / INSCATTER_HISTORY_REJECTION_THRESHOLD);

        // Adjust blend factor based on confidence
        blendFactor = lerp(INSCATTER_TEMPORAL_BLEND_MAX, INSCATTER_TEMPORAL_BLEND, confidence);
    }
    else
    {
        // No valid history - use current frame only
        blendFactor = 1.0;
    }

    // Temporal blend
    float4 accumulatedInscatter = lerp(historyInscatter, currentInscatter, blendFactor);

    // Firefly reduction - clamp extreme values
    float luminance = dot(accumulatedInscatter.rgb, float3(0.299, 0.587, 0.114));
    float maxLuminance = 10.0;
    if (luminance > maxLuminance)
    {
        accumulatedInscatter.rgb *= maxLuminance / luminance;
    }

    // Write result
    volumetricResolvedInscatterRW[dispatchThreadID] = accumulatedInscatter;

    // Also copy transmittance
    volumetricResolvedTransmissionRW[dispatchThreadID] = float4(accumulatedInscatter.a, accumulatedInscatter.a, accumulatedInscatter.a, 1.0);
}
