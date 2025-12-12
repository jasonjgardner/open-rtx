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
// CalculateGIInscatterInline - Global Illumination Volumetric Inscatter Pass
// Calculates indirect light (bounced light, emissives) scattered through fog
// Uses a lower resolution 3D grid (128x64x32) for performance
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/Settings.hlsl"
#include "Include/Sky.hlsl"
#include "Include/VolumetricLighting.hlsl"
#include "Include/Util.hlsl"

// =============================================================================
// GI FROXEL GRID PARAMETERS
// =============================================================================

// GI froxel grid is lower resolution for performance
// volumetricGIInscatterRW is 128x64x32 (R16G16B16A16_FLOAT)
static const uint3 kGIFroxelDimensions = uint3(128, 64, 32);

// GI-specific depth distribution
#ifndef GI_FROXEL_DEPTH_EXPONENT
#define GI_FROXEL_DEPTH_EXPONENT 2.5
#endif

// GI maximum depth (can be larger since resolution is lower)
#ifndef GI_FROXEL_MAX_DEPTH
#define GI_FROXEL_MAX_DEPTH 192.0
#endif

// Number of GI probe directions for spherical sampling
#ifndef GI_INSCATTER_DIRECTIONS
#define GI_INSCATTER_DIRECTIONS 6
#endif

// GI inscatter intensity multiplier
#ifndef GI_INSCATTER_INTENSITY
#define GI_INSCATTER_INTENSITY 1.0
#endif

// Enable emissive light sampling in fog
#ifndef ENABLE_EMISSIVE_FOG
#define ENABLE_EMISSIVE_FOG 1
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

// Convert depth to GI froxel slice
float depthToGISlice(float depth, float maxSlices)
{
    float normalizedDepth = saturate(depth / GI_FROXEL_MAX_DEPTH);
    return maxSlices * pow(normalizedDepth, 1.0 / GI_FROXEL_DEPTH_EXPONENT);
}

// Get world position from GI froxel coordinate
float3 giFroxelToWorldPos(uint3 froxelCoord, float3 cameraPos, float3 cameraForward, float3 cameraRight, float3 cameraUp, float2 fovTan)
{
    // Normalized froxel coordinates [0, 1]
    float2 uv = (float2(froxelCoord.xy) + 0.5) / float2(kGIFroxelDimensions.xy);

    // Convert to view-space direction [-1, 1]
    float2 ndc = uv * 2.0 - 1.0;

    // Calculate view ray direction
    float3 viewDir = normalize(
        cameraForward +
        cameraRight * (ndc.x * fovTan.x) +
        cameraUp * (ndc.y * fovTan.y)
    );

    // Get depth for this slice
    float depth = giSliceToDepth(float(froxelCoord.z) + 0.5, float(kGIFroxelDimensions.z));

    return cameraPos + viewDir * depth;
}

// Spherical Fibonacci distribution for uniform sphere sampling
float3 sphericalFibonacci(uint index, uint numSamples)
{
    const float goldenRatio = 1.618033988749895;
    float i = float(index) + 0.5;

    float phi = 2.0 * kPi * frac(i / goldenRatio);
    float cosTheta = 1.0 - 2.0 * i / float(numSamples);
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));

    return float3(
        sinTheta * cos(phi),
        cosTheta,
        sinTheta * sin(phi)
    );
}

// Hash function for random rotation
float3x3 randomRotation(uint3 coord, uint frameIndex)
{
    uint seed = coord.x * 73856093 ^ coord.y * 19349663 ^ coord.z * 83492791 ^ frameIndex * 12345;
    seed = seed * 0x27d4eb2d;

    float angle = float(seed) / 4294967295.0 * kTwoPi;
    float c = cos(angle);
    float s = sin(angle);

    // Simple rotation around Y axis
    return float3x3(
        c, 0, s,
        0, 1, 0,
        -s, 0, c
    );
}

// Estimate ambient occlusion for a position (simplified)
float estimateAO(float3 worldPos)
{
    // Simple height-based AO approximation
    float heightAboveGround = max(0.0, worldPos.y - 64.0);
    float heightAO = saturate(heightAboveGround / 32.0);

    return lerp(0.5, 1.0, heightAO);
}

// Sample sky radiance for a direction
float3 sampleSkyRadiance(float3 direction, float3 sunDir, float3 sunColor)
{
    // Simplified sky model
    float sunDot = saturate(dot(direction, sunDir));
    float horizonFactor = saturate(direction.y + 0.1);

    // Sky gradient
    float3 zenithColor = float3(0.2, 0.4, 0.8);
    float3 horizonColor = float3(0.6, 0.7, 0.9);
    float3 skyColor = lerp(horizonColor, zenithColor, saturate(direction.y));

    // Sun contribution
    float sunGlow = pow(sunDot, 32.0) * 0.5;
    skyColor += sunColor * sunGlow * 0.1;

    // Night sky
    float nightFactor = saturate(-sunDir.y);
    skyColor = lerp(skyColor, NIGHT_SKY_COLOR, nightFactor);

    return skyColor * horizonFactor;
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.CalculateGIInscatterInline.group_size}})]
void CalculateGIInscatterInline(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID
)
{
    // Note that g_rootConstant0 from BlurGradients pass is accessible here
    uint outputBufferIndex = g_rootConstant0;

#if !ENABLE_VOLUMETRIC_LIGHTING
    volumetricGIInscatterRW[outputBufferIndex][dispatchThreadID] = float4(0, 0, 0, 1);
    return;
#endif

    // Bounds check
    if (any(dispatchThreadID >= kGIFroxelDimensions))
        return;

    // Get camera parameters
    float3 cameraPos = g_view.viewOriginSteveSpace;

    // Derive view vectors from getViewDirection (using NDC corners)
    float3 cameraForward = getViewDirection(float2(0.0, 0.0));
    float3 cameraRight = normalize(getViewDirection(float2(0.1, 0.0)) - getViewDirection(float2(-0.1, 0.0)));
    float3 cameraUp = normalize(getViewDirection(float2(0.0, -0.1)) - getViewDirection(float2(0.0, 0.1)));

    float2 fovTan = float2(1.0 / g_view.proj[0][0], 1.0 / g_view.proj[1][1]);

    // Get world position for this froxel
    float3 worldPos = giFroxelToWorldPos(
        dispatchThreadID,
        cameraPos,
        cameraForward,
        cameraRight,
        cameraUp,
        fovTan
    );

    float3 rayDir = normalize(worldPos - cameraPos);
    float depth = length(worldPos - cameraPos);

    // Get lighting parameters
    float3 sunDir = g_view.directionToSun;
    float3 sunColor = blackbodyColor(SUN_COLOR_TEMPERATURE) * SUN_INTENSITY;

    // Calculate fog properties
    float baseDensity = getFogDensity(worldPos, g_view.rainLevel, g_view.time, sunDir, g_view.cameraIsUnderWater != 0);
    float3 scattering = FOG_SCATTERING_COEFFICIENTS * baseDensity * 10.0;
    float3 extinction = scattering + float3(0.0002, 0.0003, 0.0005);

    // Accumulated GI inscatter
    float3 giInscatter = float3(0, 0, 0);

    // Random rotation for this froxel (for temporal stability)
    // Use time to derive a pseudo-frame index for temporal variation
    uint pseudoFrameIndex = uint(g_view.time * 60.0) % 256;
    float3x3 rotation = randomRotation(dispatchThreadID, pseudoFrameIndex);

    // Sample incoming radiance from multiple directions
    [unroll]
    for (uint i = 0; i < GI_INSCATTER_DIRECTIONS; i++)
    {
        // Get sample direction using spherical Fibonacci
        float3 sampleDir = sphericalFibonacci(i, GI_INSCATTER_DIRECTIONS);
        sampleDir = mul(rotation, sampleDir);

        // Sample sky radiance for this direction
        float3 incomingRadiance = sampleSkyRadiance(sampleDir, sunDir, sunColor);

        // Phase function for this direction pair (isotropic for GI)
        float phase = isotropicPhase();

        // Accumulate contribution
        giInscatter += incomingRadiance * phase;
    }

    // Normalize by number of samples
    giInscatter /= float(GI_INSCATTER_DIRECTIONS);

    // Apply scattering coefficient
    giInscatter *= scattering;

    // Apply ambient occlusion estimate
    float ao = estimateAO(worldPos);
    giInscatter *= ao;

    // Apply GI intensity
    giInscatter *= GI_INSCATTER_INTENSITY;

#if ENABLE_EMISSIVE_FOG
    // Add contribution from nearby emissive blocks (simplified)
    // In a full implementation, this would sample the irradiance cache
    float emissiveProbe = 0.0;

    // Torch/emissive contribution based on distance from likely light positions
    // This is a placeholder - real implementation would query actual light positions
    float3 emissiveColor = float3(1.0, 0.7, 0.4) * EMISSIVE_GI_STRENGTH * 0.1;
    giInscatter += emissiveColor * scattering * emissiveProbe;
#endif

    // Calculate slice thickness for integration
    float sliceNear = giSliceToDepth(float(dispatchThreadID.z), float(kGIFroxelDimensions.z));
    float sliceFar = giSliceToDepth(float(dispatchThreadID.z + 1), float(kGIFroxelDimensions.z));
    float sliceThickness = sliceFar - sliceNear;

    // Calculate transmittance
    float3 sliceExtinction = extinction * sliceThickness;
    float transmittance = exp(-dot(sliceExtinction, float3(0.333, 0.333, 0.334)));

    // Scale by slice thickness
    giInscatter *= sliceThickness;

    // Store result
    volumetricGIInscatterRW[outputBufferIndex][dispatchThreadID] = float4(giInscatter, transmittance);
}
