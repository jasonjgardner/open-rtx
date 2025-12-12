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
// BlurGIInscatter - Spatial Blur Pass for GI Volumetrics
// Applies a 3D separable Gaussian blur to reduce noise in GI volumetric data
// Uses ping-pong buffer approach for multiple passes
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/Settings.hlsl"

// =============================================================================
// BLUR PARAMETERS
// =============================================================================

// GI froxel grid dimensions (must match CalculateGIInscatterInline.hlsl)
static const uint3 kGIFroxelDimensions = uint3(128, 64, 32);

// Blur kernel radius (in froxels)
#ifndef GI_BLUR_RADIUS
#define GI_BLUR_RADIUS 2
#endif

// Blur sigma for Gaussian kernel
#ifndef GI_BLUR_SIGMA
#define GI_BLUR_SIGMA 1.5
#endif

// Depth-aware blur threshold (prevents blurring across depth discontinuities)
#ifndef GI_BLUR_DEPTH_THRESHOLD
#define GI_BLUR_DEPTH_THRESHOLD 0.2
#endif

// =============================================================================
// GAUSSIAN KERNEL
// =============================================================================

// Pre-computed Gaussian weights for radius 2
// sigma = 1.5, normalized
static const float kGaussianWeights[5] = {
    0.153388,   // center
    0.221461,   // offset 1
    0.221461,
    0.100970,   // offset 2
    0.100970
};

// Compute Gaussian weight at given distance
float gaussianWeight(float distance, float sigma)
{
    return exp(-(distance * distance) / (2.0 * sigma * sigma));
}

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Safe texture fetch with bounds clamping
float4 sampleGIFroxel(uint3 coord, uint bufferIndex)
{
    // Clamp to valid bounds
    coord = clamp(coord, uint3(0, 0, 0), kGIFroxelDimensions - 1);
    return volumetricGIInscatterRW[bufferIndex][coord];
}

// Convert froxel Z to linear depth for depth-aware blur
float froxelZToDepth(uint z)
{
    float normalizedZ = (float(z) + 0.5) / float(kGIFroxelDimensions.z);
    // Match the exponential distribution from CalculateGIInscatterInline
    return 192.0 * pow(normalizedZ, 2.5);
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.BlurGIInscatter.group_size}})]
void BlurGIInscatter(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID
)
{
    // Buffer ping-pong indices
    uint volumetricGIInscatterBufferIndexFrom = 1 - g_rootConstant0;
    uint volumetricGIInscatterBufferIndexTo = g_rootConstant0;

#if !ENABLE_VOLUMETRIC_LIGHTING
    // Pass through if volumetrics disabled
    if (all(dispatchThreadID < kGIFroxelDimensions))
    {
        volumetricGIInscatterRW[volumetricGIInscatterBufferIndexTo][dispatchThreadID] =
            volumetricGIInscatterRW[volumetricGIInscatterBufferIndexFrom][dispatchThreadID];
    }
    return;
#endif

    // Bounds check
    if (any(dispatchThreadID >= kGIFroxelDimensions))
        return;

    // Get center sample
    float4 centerSample = sampleGIFroxel(dispatchThreadID, volumetricGIInscatterBufferIndexFrom);
    float centerDepth = froxelZToDepth(dispatchThreadID.z);

    // Accumulated blur result
    float4 blurResult = float4(0, 0, 0, 0);
    float totalWeight = 0.0;

    // 3D separable blur
    // For performance, we do a single pass with a small 3D kernel
    // instead of 3 separate 1D passes

    [unroll]
    for (int dz = -GI_BLUR_RADIUS; dz <= GI_BLUR_RADIUS; dz++)
    {
        [unroll]
        for (int dy = -GI_BLUR_RADIUS; dy <= GI_BLUR_RADIUS; dy++)
        {
            [unroll]
            for (int dx = -GI_BLUR_RADIUS; dx <= GI_BLUR_RADIUS; dx++)
            {
                // Calculate offset coordinate
                int3 offset = int3(dx, dy, dz);
                int3 sampleCoord = int3(dispatchThreadID) + offset;

                // Skip out of bounds (clamping handled in fetch)
                if (any(sampleCoord < 0) || any(sampleCoord >= int3(kGIFroxelDimensions)))
                    continue;

                // Sample neighbor
                float4 neighborSample = sampleGIFroxel(uint3(sampleCoord), volumetricGIInscatterBufferIndexFrom);
                float neighborDepth = froxelZToDepth(uint(sampleCoord.z));

                // Calculate spatial Gaussian weight
                float dist = length(float3(offset));
                float spatialWeight = gaussianWeight(dist, GI_BLUR_SIGMA);

                // Depth-aware weight (reduce blur across depth discontinuities)
                float depthDiff = abs(centerDepth - neighborDepth) / max(centerDepth, 1.0);
                float depthWeight = exp(-depthDiff * depthDiff / (2.0 * GI_BLUR_DEPTH_THRESHOLD * GI_BLUR_DEPTH_THRESHOLD));

                // Transmittance similarity weight
                // Prevents blurring across different fog densities
                float transmittanceDiff = abs(centerSample.a - neighborSample.a);
                float transmittanceWeight = exp(-transmittanceDiff * transmittanceDiff * 4.0);

                // Combined weight
                float weight = spatialWeight * depthWeight * transmittanceWeight;

                // Accumulate
                blurResult += neighborSample * weight;
                totalWeight += weight;
            }
        }
    }

    // Normalize
    if (totalWeight > 0.0001)
    {
        blurResult /= totalWeight;
    }
    else
    {
        blurResult = centerSample;
    }

    // Preserve transmittance channel more accurately
    // Blend between blurred and original transmittance
    blurResult.a = lerp(centerSample.a, blurResult.a, 0.5);

    // Output to destination buffer
    volumetricGIInscatterRW[volumetricGIInscatterBufferIndexTo][dispatchThreadID] = blurResult;
}
