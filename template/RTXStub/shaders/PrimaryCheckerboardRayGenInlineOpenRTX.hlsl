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
// OpenRTX Primary Ray Generation
// Enhanced ray tracing with physically-based rendering
// =============================================================================

// Enable glass backface culling like vanilla
#define CULL_GLASS_BACK_FACES 1

#include "Include/Renderer.hlsl"
#include "Include/Util.hlsl"
#include "Include/OpenRTX.hlsl"

// =============================================================================
// ENHANCED RAY STATE
// =============================================================================

struct OpenRTXRayState
{
    RayDesc rayDesc;

    float3 color;
    float3 throughput;

    float distance;
    float3 motion;

    uint instanceMask;

    // OpenRTX additions
    float3 normal;
    float roughness;
    float metalness;
    bool hitWater;
    float waterDepth;

    void Init()
    {
        color = 0;
        throughput = 1;
        distance = 0;
        motion = 0;
        instanceMask = 0xff;

        normal = float3(0, 1, 0);
        roughness = 0.5;
        metalness = 0;
        hitWater = false;
        waterDepth = 0;
    }
};

// =============================================================================
// ENHANCED RENDERING FUNCTIONS
// =============================================================================

void RenderSkyOpenRTX(inout OpenRTXRayState rayState, OpenRTXContext ctx)
{
    if (all(rayState.throughput == 0))
        return;

#if OPENRTX_ENABLED && ENABLE_ATMOSPHERIC_SKY
    // Use OpenRTX atmospheric sky
    float3 skyColor = renderSkyWithClouds(rayState.rayDesc.Direction, ctx);
    rayState.color += rayState.throughput * skyColor;
#else
    // Fallback to vanilla-style sky
    const float3 skyColor = float3(170, 209, 254) / 255;
    const float3 gradientColor = float3(121, 167, 255) / 255;

    const float3 nightSkyColor = float3(10, 12, 22) / 255;
    const float3 nightGradientColor = float3(1, 1, 2) / 255;

    float gradientLerp = max(0.0, lerp(-0.15, 1.0, rayState.rayDesc.Direction.y));
    gradientLerp = pow(gradientLerp, 0.5);

    const float nightThreshold = -0.3;
    const float dayThreshold = 0.2;
    float timeOfDayLerp = saturate((ctx.sunDir.y - nightThreshold) / (dayThreshold - nightThreshold));

    float3 dayColor = lerp(skyColor, gradientColor, gradientLerp);
    float3 nightColor = lerp(nightSkyColor, nightGradientColor, gradientLerp);

    float3 finalColor = lerp(nightColor, dayColor, timeOfDayLerp);
    rayState.color += rayState.throughput * finalColor;
#endif
}

void RenderVanillaOpenRTX(HitInfo hitInfo, inout OpenRTXRayState rayState, OpenRTXContext ctx)
{
    ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
    GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
    SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

    float3 worldPos = surfaceInfo.position - g_view.waveWorksOriginInSteveSpace;
    worldPos = worldPos - floor(worldPos / 1024) * 1024;

    // Store surface properties
    rayState.normal = surfaceInfo.normal;
    rayState.roughness = surfaceInfo.roughness;
    rayState.metalness = surfaceInfo.metalness;

    // Check for water
    bool isWater = hitInfo.materialType == MATERIAL_TYPE_WATER;
    if (isWater)
    {
        rayState.hitWater = true;
    }

#if OPENRTX_ENABLED && ENABLE_PBR_LIGHTING
    // Enhanced PBR lighting
    EnhancedSurface surface;
    surface.position = surfaceInfo.position;
    surface.normal = surfaceInfo.normal;
    surface.geometryNormal = geometryInfo.geometryNormal;
    surface.viewDir = -rayState.rayDesc.Direction;
    surface.albedo = surfaceInfo.color;
    surface.roughness = max(surfaceInfo.roughness, MIN_ROUGHNESS);
    surface.metalness = surfaceInfo.metalness;
    surface.ao = 1.0; // Would need AO buffer
    surface.subsurface = surfaceInfo.subsurface;
    surface.emissive = surfaceInfo.emissive;
    surface.isWater = isWater;
    surface.waterDepth = 0.0;

    // Apply default roughness fix
#if FIX_DEFAULT_MATERIAL
    if (geometryInfo.pbrTextureDataIndex == kInvalidPBRTextureHandle)
    {
        surface.roughness = DEFAULT_ROUGHNESS;
    }
#endif

    float3 light = shadeSurfacePBR(surface, ctx);
#else
    // Vanilla-like shading
    float3 light = lerp(
        lerp(0.6, 0.8, abs(dot(surfaceInfo.normal, float3(0, 0, 1)))),
        lerp(0.45, 1, mad(dot(surfaceInfo.normal, float3(0, 1, 0)), 0.5, 0.5)),
        abs(dot(surfaceInfo.normal, float3(0, 1, 0))));

    // Apply emissive
    light = lerp(light, 1, surfaceInfo.emissive);
#endif

    // Force full alpha for opaque/alphatest
    if (hitInfo.materialType == MATERIAL_TYPE_OPAQUE || hitInfo.materialType == MATERIAL_TYPE_ALPHA_TEST)
        surfaceInfo.alpha = 1;

    // Cloud special handling
    if (objectInstance.flags & kObjectInstanceFlagClouds)
    {
#if OPENRTX_ENABLED
        light = geometryInfo.color.rgb;
#else
        light = geometryInfo.color.rgb;
#endif
        surfaceInfo.alpha = 0.7;
    }

    // Point lights
    for (int i = 0; i < min(10, g_view.cpuLightsCount); i++)
    {
        LightInfo lightInfo = inputLightsBuffer[i];
        LightData lightData = UnpackLight(lightInfo.packedData);

        float3 lDir = lightInfo.position - surfaceInfo.position;
        float lDist = length(lDir);
        lDir /= lDist;

        float attenuation = max(0, dot(surfaceInfo.normal, lDir)) / (lDist * lDist);

#if OPENRTX_ENABLED && ENABLE_PBR_LIGHTING
        // PBR point light contribution
        float3 H = normalize(-rayState.rayDesc.Direction + lDir);
        float NdotL = saturate(dot(surfaceInfo.normal, lDir));
        float NdotH = saturate(dot(surfaceInfo.normal, H));

        light += 100 * attenuation * lightData.intensity * lightData.color;
#else
        light += 100 * attenuation * lightData.intensity * lightData.color;
#endif
    }

    // Determine blending mode
    const bool isBlockBreakingOverlay = objectInstance.flags == (kObjectInstanceFlagAlphaTestThresholdHalf | kObjectInstanceFlagTextureAlphaControlsVertexColor);

    float3 throughput;
    float3 emission;

    if (objectInstance.flags & (kObjectInstanceFlagSun | kObjectInstanceFlagMoon))
    {
        // Additive blending for sun/moon
        throughput = 1;
        float meshIntensity = (objectInstance.flags & kObjectInstanceFlagSun) ? g_view.sunMeshIntensity : g_view.moonMeshIntensity;
        emission = surfaceInfo.color * meshIntensity * surfaceInfo.alpha;
    }
    else if (isBlockBreakingOverlay)
    {
#if FIX_BLOCK_BREAKING_OVERLAY
        // Multiplicative blending
        throughput = surfaceInfo.color;
        emission = 0;
#else
        throughput = surfaceInfo.color;
        emission = 0;
#endif
    }
    else
    {
        // Standard alpha blend
        throughput = 1 - surfaceInfo.alpha;
        emission = surfaceInfo.color * surfaceInfo.alpha * light;
    }

    // Glint effect
    if (objectInstance.flags & kObjectInstanceFlagGlint)
        emission += (sin(3.0 * g_view.time) * 0.5 + 0.5) * (float3(077, 23, 255) / 255.0);

    // Advance ray
    rayState.rayDesc.TMin = hitInfo.rayT;

    // Accumulate
    rayState.color += emission * rayState.throughput;
    rayState.throughput *= throughput;

    // Update state
    rayState.distance = hitInfo.rayT;
    rayState.motion = surfaceInfo.position - surfaceInfo.prevPosition;
}

float3 RenderRayOpenRTX(RayDesc rayDesc, out float outputDistance, out float3 outputMotion, OpenRTXContext ctx)
{
    RayQuery<RAY_FLAG_NONE> q;

    OpenRTXRayState rayState;
    rayState.Init();
    rayState.rayDesc = rayDesc;

    // Limit translucent surfaces
    for (int i = 0; i < 100; i++)
    {
        q.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, rayState.instanceMask, rayState.rayDesc);

        while (q.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(q);
            if (AlphaTestHitLogic(hitInfo))
            {
                q.CommitNonOpaqueTriangleHit();
            }
        }

        if (q.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(q);
            RenderVanillaOpenRTX(hitInfo, rayState, ctx);
        }
        else
        {
            break;
        }

        if (all(rayState.throughput == 0))
            break;
    }

    const float maxDistance = 65504;

    if (all(rayState.throughput == 0))
    {
        outputDistance = min(rayState.distance, maxDistance);
        outputMotion = rayState.motion;
    }
    else
    {
        outputDistance = maxDistance;
        outputMotion = 0;
    }

    // Render sky
    RenderSkyOpenRTX(rayState, ctx);

    // Apply volumetric effects
#if OPENRTX_ENABLED && ENABLE_VOLUMETRIC_LIGHTING
    float3 finalColor = applyAtmosphericEffects(rayState.color, rayDesc.Direction, rayState.distance, ctx);
#else
    float3 finalColor = rayState.color;
#endif

    return finalColor;
}

// =============================================================================
// MAIN ENTRY POINT
// =============================================================================

[numthreads(4, 8, 1)]
void PrimaryCheckerboardRayGenInline(
    uint3 dispatchThreadID : SV_DispatchThreadID,
    uint3 groupThreadID : SV_GroupThreadID,
    uint groupIndex : SV_GroupIndex,
    uint3 groupID : SV_GroupID)
{
    if (any(dispatchThreadID.xy >= g_view.renderResolution))
        return;

    // Set up ray
    RayDesc rayDesc;
    rayDesc.Direction = rayDirFromNDC(getNDCjittered(dispatchThreadID.xy));
    rayDesc.Origin = g_view.viewOriginSteveSpace;
    rayDesc.TMin = 0;
    rayDesc.TMax = 10000;

    // Initialize OpenRTX context with game-provided values
    float3 sunDir = getTrueDirectionToSun();

    OpenRTXContext ctx = initContextFromGame(
        rayDesc.Origin,
        rayDesc.Direction,
        (float2(dispatchThreadID.xy) + 0.5) * g_view.recipRenderResolution,
        g_view.time,
        sunDir,
        g_view.sunColour,           // Game-provided sun color
        g_view.skyColor,            // Game-provided sky color
        g_view.skyColorUp,          // Upper sky gradient
        g_view.skyColorDown,        // Lower sky gradient
        g_view.constantAmbient,     // Ambient light
        g_view.skyIntensityAdjustment,
        g_view.sunMeshIntensity,
        g_view.moonMeshIntensity,
        g_view.rainLevel,           // Rain intensity from game
        0); // Overworld dimension

    // Render
    float hitDist;
    float3 objMotion;
    float3 color = RenderRayOpenRTX(rayDesc, hitDist, objMotion, ctx);

    // Calculate motion vector
    float2 motionVector = computeMotionVector(rayDesc.Origin + rayDesc.Direction * hitDist, objMotion);

    // Debug NaN/Inf visualization
#if DEBUG_NAN_INF
    if (any(isinf(color)) || any(isinf(motionVector)) || isinf(hitDist))
        color = (dispatchThreadID.x / 32 + dispatchThreadID.y / 32) & 1 ? float3(1, 1, 0) : 0;
    if (any(isnan(color)) || any(isnan(motionVector)) || isnan(hitDist))
        color = (dispatchThreadID.x / 32 + dispatchThreadID.y / 32) & 1 ? float3(1, 0, 1) : 0;
#endif

    // Output
    outputBufferFinal[dispatchThreadID.xy] = float4(color, 1);
    outputBufferMotionVectors[dispatchThreadID.xy] = motionVector;
    outputBufferReprojectedPathLength[dispatchThreadID.xy] = hitDist;
}
