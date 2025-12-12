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
// CalculateInscatterInline - Direct Light Volumetric Inscatter Pass
// Calculates sun/moon light scattered through atmospheric fog into a 3D froxel grid
// =============================================================================

#include "Include/Generated/Signature.hlsl"
#include "Include/Settings.hlsl"
#include "Include/VolumetricLighting.hlsl"

// =============================================================================
// FROXEL GRID PARAMETERS
// =============================================================================

// Froxel grid dimensions (matches buffer layout from Signature.hlsl)
// volumetricInscatterRW is 256x128x64 (R16G16B16A16_FLOAT)
static const uint3 kFroxelDimensions = uint3(256, 128, 64);

// Depth distribution exponent for non-linear z-slicing
// Higher values place more slices near the camera
#ifndef FROXEL_DEPTH_EXPONENT
#define FROXEL_DEPTH_EXPONENT 2.0
#endif

// Maximum depth for froxel grid (in blocks/meters)
#ifndef FROXEL_MAX_DEPTH
#define FROXEL_MAX_DEPTH 256.0
#endif

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Convert froxel grid coordinate to world depth using exponential distribution
float froxelSliceToDepth(float slice, float maxSlices)
{
    float normalizedSlice = slice / maxSlices;
    // Exponential distribution places more detail near camera
    return FROXEL_MAX_DEPTH * pow(normalizedSlice, FROXEL_DEPTH_EXPONENT);
}

// Convert depth to froxel slice (inverse of above)
float depthToFroxelSlice(float depth, float maxSlices)
{
    float normalizedDepth = saturate(depth / FROXEL_MAX_DEPTH);
    return maxSlices * pow(normalizedDepth, 1.0 / FROXEL_DEPTH_EXPONENT);
}

// Get world position from froxel grid coordinate
float3 froxelToWorldPos(uint3 froxelCoord, float3 cameraPos, float3 cameraForward, float3 cameraRight, float3 cameraUp, float2 fovTan)
{
    // Normalized froxel coordinates [0, 1]
    float2 uv = (float2(froxelCoord.xy) + 0.5) / float2(kFroxelDimensions.xy);

    // Convert to view-space direction [-1, 1]
    float2 ndc = uv * 2.0 - 1.0;

    // Calculate view ray direction using FOV tangent
    float3 viewDir = normalize(
        cameraForward +
        cameraRight * (ndc.x * fovTan.x) +
        cameraUp * (ndc.y * fovTan.y)
    );

    // Get depth for this slice
    float depth = froxelSliceToDepth(float(froxelCoord.z) + 0.5, float(kFroxelDimensions.z));

    // Calculate world position
    return cameraPos + viewDir * depth;
}

// Blue noise dithering for temporal stability
float blueNoiseDither(uint3 coord, uint frameIndex)
{
    // Simple hash-based blue noise approximation
    uint seed = coord.x + coord.y * 256 + coord.z * 65536 + frameIndex * 16777216;
    seed = (seed ^ 61) ^ (seed >> 16);
    seed = seed + (seed << 3);
    seed = seed ^ (seed >> 4);
    seed = seed * 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return float(seed) / 4294967295.0;
}

// =============================================================================
// MAIN COMPUTE SHADER
// =============================================================================

[numthreads({{RTXStub.passes.CalculateInscatterInline.group_size}})]
void CalculateInscatterInline(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID
)
{
#if !ENABLE_VOLUMETRIC_LIGHTING
    // Early out if volumetric lighting is disabled
    volumetricInscatterRW[dispatchThreadID] = float4(0, 0, 0, 1);
    return;
#endif

    // Bounds check
    if (any(dispatchThreadID >= kFroxelDimensions))
        return;

    // Get camera parameters from view constants
    float3 cameraPos = g_view.viewOriginSteveSpace;
    float3 cameraForward = -g_view.viewDir;  // View matrix looks down -Z
    float3 cameraRight = g_view.viewRight;
    float3 cameraUp = g_view.viewUp;

    // Calculate FOV tangent from projection matrix
    // Using aspect ratio and vertical FOV from projection
    float2 fovTan = float2(1.0 / g_view.proj[0][0], 1.0 / g_view.proj[1][1]);

    // Get world position for this froxel
    float3 worldPos = froxelToWorldPos(
        dispatchThreadID,
        cameraPos,
        cameraForward,
        cameraRight,
        cameraUp,
        fovTan
    );

    // Calculate ray direction from camera to froxel
    float3 rayDir = normalize(worldPos - cameraPos);
    float depth = length(worldPos - cameraPos);

    // Get sun/light parameters
    float3 sunDir = g_view.directionToSun;
    float3 sunColor = blackbodyColor(SUN_COLOR_TEMPERATURE) * SUN_INTENSITY;

    // Apply time-of-day color temperature
    float sunHeight = sunDir.y;
    if (sunHeight < 0.3)
    {
        // Sunrise/sunset warm tint
        float warmth = smoothstep(0.3, -0.1, sunHeight);
        float3 warmColor = blackbodyColor(lerp(SUN_COLOR_TEMPERATURE, SUNRISE_COLOR_TEMPERATURE, warmth));
        sunColor = warmColor * SUN_INTENSITY * saturate(sunHeight + 0.1) * 5.0;
    }

    // Blue noise jitter for temporal stability
    float jitter = blueNoiseDither(dispatchThreadID, g_view.frameIndex);
    float3 jitteredPos = worldPos + rayDir * (jitter - 0.5) * (depth / float(kFroxelDimensions.z));

    // Calculate fog density at this position
    float baseDensity = getFogDensity(jitteredPos, g_view.rainLevel, g_view.time, sunDir, g_view.cameraIsUnderWater != 0);

    // Phase function - how light scatters toward the viewer
    float cosTheta = dot(rayDir, sunDir);
    float phase = dualLobePhase(cosTheta, 0.8, -0.3, 0.7);

    // Calculate scattering and extinction
    float3 scattering = FOG_SCATTERING_COEFFICIENTS * baseDensity * 10.0;
    float3 extinction = scattering + float3(0.0002, 0.0003, 0.0005);

    // In-scattered light from sun
    float3 inscatter = float3(0, 0, 0);

#if ENABLE_SUN_FOG
    // Sun visibility (simplified shadow test using height)
    float sunVisibility = 1.0;

    // Basic terrain shadow approximation
    float terrainHeight = 64.0;  // Sea level
    if (sunDir.y < 0.2 && jitteredPos.y < terrainHeight + 32.0)
    {
        // Approximate shadow from terrain when sun is low
        float shadowFactor = smoothstep(terrainHeight, terrainHeight + 32.0, jitteredPos.y);
        sunVisibility *= shadowFactor;
    }

    // Cloud shadow approximation
    if (sunDir.y > 0.0)
    {
        float cloudShadow = 1.0;
        #if ENABLE_VOLUMETRIC_CLOUDS
        // Sample cloud density above this point (simplified)
        float3 cloudSamplePos = jitteredPos + sunDir * 100.0;
        if (cloudSamplePos.y > CLOUD_MIN_HEIGHT && cloudSamplePos.y < CLOUD_MAX_HEIGHT)
        {
            cloudShadow = 1.0 - CLOUD_SHADOW_STRENGTH * 0.5;
        }
        #endif
        sunVisibility *= cloudShadow;
    }

    inscatter += sunColor * phase * sunVisibility * scattering * SUN_FOG_AMOUNT;
#endif

#if ENABLE_STATIC_GI_FOG
    // Ambient/sky contribution (isotropic)
    float3 ambientColor = float3(0.08, 0.10, 0.14);

    // Modulate ambient by sun visibility (darker in shadows)
    float ambientOcclusion = lerp(0.3, 1.0, saturate(sunDir.y + 0.2));
    inscatter += ambientColor * scattering * STATIC_GI_FOG_AMOUNT * ambientOcclusion;
#endif

    // Calculate transmittance for this slice
    // Use slice thickness as the integration distance
    float sliceNear = froxelSliceToDepth(float(dispatchThreadID.z), float(kFroxelDimensions.z));
    float sliceFar = froxelSliceToDepth(float(dispatchThreadID.z + 1), float(kFroxelDimensions.z));
    float sliceThickness = sliceFar - sliceNear;

    float3 sliceExtinction = extinction * sliceThickness;
    float transmittance = exp(-dot(sliceExtinction, float3(0.333, 0.333, 0.334)));

    // Scale inscatter by slice thickness for proper integration
    inscatter *= sliceThickness;

    // Store result: RGB = inscattered light, A = transmittance
    volumetricInscatterRW[dispatchThreadID] = float4(inscatter, transmittance);
}
