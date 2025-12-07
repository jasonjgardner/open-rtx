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
    float3 position;
    float3 albedo;
    float roughness;
    float metalness;
    bool hitWater;
    float waterDepth;
    bool hitOpaque;

    void Init()
    {
        color = 0;
        throughput = 1;
        distance = 0;
        motion = 0;
        instanceMask = 0xff;

        normal = float3(0, 1, 0);
        position = 0;
        albedo = 0.5;
        roughness = 0.5;
        metalness = 0;
        hitWater = false;
        waterDepth = 0;
        hitOpaque = false;
    }
};

// =============================================================================
// SHADOW RAY TRACING WITH COLOR TRANSMISSION
// =============================================================================

struct ShadowResult
{
    float visibility;       // 0 = fully shadowed, 1 = fully lit
    float3 transmission;    // Color transmitted through translucent surfaces
};

// Trace shadow ray with color transmission through translucent surfaces
ShadowResult traceShadowRay(float3 origin, float3 direction, float maxDist, OpenRTXContext ctx)
{
    ShadowResult result;
    result.visibility = 1.0;
    result.transmission = 1.0;

#if ENABLE_RAYTRACED_SHADOWS
    RayQuery<RAY_FLAG_NONE> shadowQuery;
    RayDesc shadowRay;
    shadowRay.Origin = origin + direction * 0.001;  // Bias to avoid self-intersection
    shadowRay.Direction = direction;
    shadowRay.TMin = 0.0;
    shadowRay.TMax = maxDist;

    float3 colorAccum = 1.0;

    // Trace through multiple translucent surfaces
    for (int bounce = 0; bounce < 8; bounce++)
    {
        shadowQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, shadowRay);

        while (shadowQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(shadowQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                shadowQuery.CommitNonOpaqueTriangleHit();
            }
        }

        if (shadowQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(shadowQuery);
            ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];

            // Skip cloud geometry for shadow tests (vanilla mesh clouds)
            if (objectInstance.flags & kObjectInstanceFlagClouds)
            {
                // Continue past clouds
                shadowRay.Origin = shadowRay.Origin + shadowRay.Direction * (hitInfo.rayT + 0.001);
                shadowRay.TMax -= hitInfo.rayT;
                if (shadowRay.TMax <= 0.0)
                    break;
                continue;
            }

            GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
            SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

            // Check if opaque hit (but not clouds)
            if (hitInfo.materialType == MATERIAL_TYPE_OPAQUE ||
                hitInfo.materialType == MATERIAL_TYPE_ALPHA_TEST)
            {
                // Opaque surface blocks light completely
                result.visibility = 0.0;
                result.transmission = 0.0;
                return result;
            }

            // Translucent surface - accumulate color transmission
#if ENABLE_COLORED_SHADOWS
            // Apply Beer's law for color absorption
            float3 surfaceColor = surfaceInfo.color;
            float alpha = surfaceInfo.alpha;

            // Color transmission through translucent material
            float3 transmittedColor = lerp(1.0, surfaceColor, alpha * COLOR_TRANSMISSION_STRENGTH);
            colorAccum *= transmittedColor;

            // Reduce visibility based on alpha
            result.visibility *= (1.0 - alpha * 0.5);
#else
            result.visibility *= (1.0 - surfaceInfo.alpha);
#endif

            // Continue ray past this surface
            shadowRay.Origin = shadowRay.Origin + shadowRay.Direction * (hitInfo.rayT + 0.001);
            shadowRay.TMax -= hitInfo.rayT;

            if (shadowRay.TMax <= 0.0 || result.visibility < 0.01)
                break;
        }
        else
        {
            // No more hits - reached the light
            break;
        }
    }

    result.transmission = colorAccum;
#endif

    return result;
}

// Generate soft shadow sample direction
float3 getSoftShadowDirection(float3 lightDir, float2 xi, float softness)
{
    // Create orthonormal basis around light direction
    float3 tangent = normalize(cross(lightDir, abs(lightDir.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0)));
    float3 bitangent = cross(lightDir, tangent);

    // Disk sampling for soft shadows
    float r = sqrt(xi.x) * softness;
    float theta = xi.y * 2.0 * kPi;

    float3 offset = tangent * (r * cos(theta)) + bitangent * (r * sin(theta));
    return normalize(lightDir + offset);
}

// =============================================================================
// REFLECTION RAY TRACING
// =============================================================================

// Simple hash for noise generation
float hashReflection(float2 p)
{
    float3 p3 = frac(float3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.x + p3.y) * p3.z);
}

// GGX importance sampling for rough reflections
float3 importanceSampleGGX(float2 xi, float3 N, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;

    float phi = 2.0 * kPi * xi.x;
    float cosTheta = sqrt((1.0 - xi.y) / (1.0 + (a2 - 1.0) * xi.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    // Spherical to cartesian
    float3 H;
    H.x = cos(phi) * sinTheta;
    H.y = sin(phi) * sinTheta;
    H.z = cosTheta;

    // Tangent space to world space
    float3 up = abs(N.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, N));
    float3 bitangent = cross(N, tangent);

    return normalize(tangent * H.x + bitangent * H.y + N * H.z);
}

// Trace a single reflection sample
float3 traceSingleReflection(float3 origin, float3 reflectDir, float3 normal, OpenRTXContext ctx)
{
    float3 reflectionColor = 0.0;

    // Trace reflection ray
    RayQuery<RAY_FLAG_NONE> reflectQuery;
    RayDesc reflectRay;
    reflectRay.Origin = origin + normal * 0.01;  // Bias along normal
    reflectRay.Direction = reflectDir;
    reflectRay.TMin = 0.0;
    reflectRay.TMax = 1000.0;

    reflectQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, reflectRay);

    while (reflectQuery.Proceed())
    {
        HitInfo hitInfo = GetCandidateHitInfo(reflectQuery);
        if (AlphaTestHitLogic(hitInfo))
        {
            reflectQuery.CommitNonOpaqueTriangleHit();
        }
    }

    if (reflectQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        HitInfo hitInfo = GetCommittedHitInfo(reflectQuery);
        ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
        GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
        SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

        // Simple shading for reflected surface
        float3 N = surfaceInfo.normal;
        float NdotL = saturate(dot(N, ctx.sunDir));

        // Direct lighting on reflected surface
        float3 directLight = surfaceInfo.color * ctx.sunColor * NdotL * 0.5;

        // Add ambient
        directLight += surfaceInfo.color * ctx.constantAmbient;

        // Add emissive
        directLight += surfaceInfo.color * surfaceInfo.emissive * 2.0;

        reflectionColor = directLight;
    }
    else
    {
        // Hit sky
        reflectionColor = renderSkyWithClouds(reflectDir, ctx);
    }

    return reflectionColor;
}

// Trace reflection ray with inline denoising (multiple samples + filtering)
float3 traceReflection(float3 origin, float3 viewDir, float3 normal, float roughness,
                       float metalness, float3 albedo, float2 pixelCoord, OpenRTXContext ctx, int bounceDepth)
{
    float3 reflectionColor = 0.0;

#if ENABLE_RAYTRACED_REFLECTIONS
    // Skip if too rough
    float effectiveRoughness = max(roughness - REFLECTION_ROUGHNESS_BIAS, MIN_ROUGHNESS);
    if (effectiveRoughness > REFLECTION_MAX_ROUGHNESS || bounceDepth >= REFLECTION_MAX_BOUNCES)
    {
        // Fall back to sky reflection for rough surfaces
        float3 reflectDir = reflect(viewDir, normal);
        return renderSkyWithClouds(reflectDir, ctx) * 0.1;
    }

    // Determine number of samples based on roughness (more samples for rougher surfaces)
    int numSamples = REFLECTION_SAMPLES;
    if (effectiveRoughness < 0.05)
    {
        numSamples = 1;  // Mirror-like surfaces need only one sample
    }

    float3 accumulatedColor = 0.0;
    float totalWeight = 0.0;

    // Use temporal jitter based on frame for better convergence
    float frameJitter = frac(g_view.time * 7.31);

    [unroll]
    for (int s = 0; s < REFLECTION_SAMPLES; s++)
    {
        if (s >= numSamples) break;

        // Generate sample direction with temporal variation
        float2 xi = float2(
            hashReflection(pixelCoord + float2(s * 17.0 + frameJitter * 100.0, 0.0)),
            hashReflection(pixelCoord + float2(0.0, s * 31.0 + frameJitter * 100.0 + 1000.0))
        );

        float3 H = importanceSampleGGX(xi, normal, effectiveRoughness);
        float3 reflectDir = reflect(viewDir, H);

        // Ensure reflection direction points away from surface
        if (dot(reflectDir, normal) < 0.0)
            reflectDir = reflect(viewDir, normal);

        // Calculate sample weight based on PDF
        float NdotH = saturate(dot(normal, H));
        float VdotH = saturate(dot(-viewDir, H));
        float weight = 1.0;  // GGX already importance sampled

        // Trace and accumulate
        float3 sampleColor = traceSingleReflection(origin, reflectDir, normal, ctx);
        accumulatedColor += sampleColor * weight;
        totalWeight += weight;
    }

    // Average samples
    if (totalWeight > 0.0)
    {
        reflectionColor = accumulatedColor / totalWeight;
    }

    // Compute Fresnel
    float3 f0 = lerp(0.04, albedo, metalness);
    float NdotV = saturate(dot(normal, -viewDir));
    float3 fresnel = f0 + (1.0 - f0) * pow(1.0 - NdotV, 5.0);

    // Apply fresnel to reflection
    reflectionColor *= fresnel;

    // Fade based on roughness
    float roughnessFade = 1.0 - saturate(effectiveRoughness / REFLECTION_MAX_ROUGHNESS);
    reflectionColor *= roughnessFade;

    // Apply simple edge-aware filtering to reduce noise
    // Using a subtle contrast-based smoothing factor
    float luminance = dot(reflectionColor, float3(0.299, 0.587, 0.114));
    float smoothFactor = saturate(1.0 - effectiveRoughness * 2.0);
    reflectionColor = lerp(reflectionColor, reflectionColor * smoothFactor + luminance * (1.0 - smoothFactor) * fresnel, effectiveRoughness);
#endif

    return reflectionColor;
}

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

void RenderVanillaOpenRTX(HitInfo hitInfo, inout OpenRTXRayState rayState, OpenRTXContext ctx, float2 pixelCoord)
{
    ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
    GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
    SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

    float3 worldPos = surfaceInfo.position - g_view.waveWorksOriginInSteveSpace;
    worldPos = worldPos - floor(worldPos / 1024) * 1024;

    // Store surface properties
    rayState.normal = surfaceInfo.normal;
    rayState.position = surfaceInfo.position;
    rayState.albedo = surfaceInfo.color;
    rayState.roughness = surfaceInfo.roughness;
    rayState.metalness = surfaceInfo.metalness;

    // Check for water
    bool isWater = hitInfo.materialType == MATERIAL_TYPE_WATER;
    if (isWater)
    {
        rayState.hitWater = true;
    }

    // Check if this is an opaque hit (for reflections later)
    bool isOpaque = (hitInfo.materialType == MATERIAL_TYPE_OPAQUE ||
                     hitInfo.materialType == MATERIAL_TYPE_ALPHA_TEST);
    if (isOpaque && !rayState.hitOpaque)
    {
        rayState.hitOpaque = true;
    }

    // =================================================================
    // RAYTRACED SHADOWS WITH COLOR TRANSMISSION
    // =================================================================
    ShadowResult shadow;
    shadow.visibility = 1.0;
    shadow.transmission = 1.0;

#if OPENRTX_ENABLED && ENABLE_RAYTRACED_SHADOWS
    // Only trace shadows for surfaces facing the sun
    float NdotL = dot(surfaceInfo.normal, ctx.sunDir);
    if (NdotL > 0.0 && isOpaque)
    {
#if ENABLE_SOFT_SHADOWS && SHADOW_SAMPLES > 1
        // Soft shadows with multiple samples
        float totalVisibility = 0.0;
        float3 totalTransmission = 0.0;

        [unroll]
        for (int s = 0; s < SHADOW_SAMPLES; s++)
        {
            float2 xi = float2(
                hashReflection(pixelCoord + float2(s * 17.0, g_view.time * 50.0)),
                hashReflection(pixelCoord + float2(g_view.time * 50.0, s * 31.0))
            );
            float3 shadowDir = getSoftShadowDirection(ctx.sunDir, xi, SHADOW_SOFTNESS);
            ShadowResult sampleShadow = traceShadowRay(surfaceInfo.position, shadowDir, 1000.0, ctx);
            totalVisibility += sampleShadow.visibility;
            totalTransmission += sampleShadow.transmission;
        }
        shadow.visibility = totalVisibility / float(SHADOW_SAMPLES);
        shadow.transmission = totalTransmission / float(SHADOW_SAMPLES);
#else
        // Single hard shadow sample
        shadow = traceShadowRay(surfaceInfo.position, ctx.sunDir, 1000.0, ctx);
#endif
    }
#endif

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
    surface.ao = 1.0;
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

    // Get direct and ambient lighting separately
    float3 directLight = shadeSurfaceDirectPBR(surface, ctx);
    float3 ambientLight = shadeSurfaceAmbientPBR(surface, ctx);

    // Apply shadow only to direct lighting, with color transmission
    float3 shadowedDirect = directLight * shadow.visibility * shadow.transmission;

    // Combine: shadowed direct + unshadowed ambient
    float3 light = shadowedDirect + ambientLight;
#else
    // Vanilla-like shading
    float3 light = lerp(
        lerp(0.6, 0.8, abs(dot(surfaceInfo.normal, float3(0, 0, 1)))),
        lerp(0.45, 1, mad(dot(surfaceInfo.normal, float3(0, 1, 0)), 0.5, 0.5)),
        abs(dot(surfaceInfo.normal, float3(0, 1, 0))));

    // Apply shadow
    light *= lerp(0.3, 1.0, shadow.visibility);

    // Apply color transmission from stained glass etc.
    light *= shadow.transmission;

    // Apply emissive
    light = lerp(light, 1, surfaceInfo.emissive);
#endif

    // =================================================================
    // RAYTRACED REFLECTIONS
    // =================================================================
    float3 reflectionContrib = 0.0;

#if OPENRTX_ENABLED && ENABLE_RAYTRACED_REFLECTIONS
    // Add reflections for smooth/metallic surfaces
    float reflectivity = surfaceInfo.metalness + (1.0 - surface.roughness) * 0.3;
    if (reflectivity > 0.1 && isOpaque)
    {
        reflectionContrib = traceReflection(
            surfaceInfo.position,
            rayState.rayDesc.Direction,
            surfaceInfo.normal,
            surface.roughness,
            surfaceInfo.metalness,
            surfaceInfo.color,
            pixelCoord,
            ctx,
            0  // Bounce depth
        );

        // For metals, reflection replaces diffuse; for dielectrics, it adds on top
        if (surfaceInfo.metalness > 0.5)
        {
            light = lerp(light, reflectionContrib, surfaceInfo.metalness);
        }
        else
        {
            light += reflectionContrib * (1.0 - surface.roughness);
        }
    }
#endif

    // Force full alpha for opaque/alphatest
    if (isOpaque)
        surfaceInfo.alpha = 1;

    // Cloud special handling
    if (objectInstance.flags & kObjectInstanceFlagClouds)
    {
        light = geometryInfo.color.rgb;
        surfaceInfo.alpha = 0.7;
    }

    // Point lights with shadow
    for (int i = 0; i < min(10, g_view.cpuLightsCount); i++)
    {
        LightInfo lightInfo = inputLightsBuffer[i];
        LightData lightData = UnpackLight(lightInfo.packedData);

        float3 lDir = lightInfo.position - surfaceInfo.position;
        float lDist = length(lDir);
        lDir /= max(lDist, 0.001);

        float attenuation = max(0, dot(surfaceInfo.normal, lDir)) / max(lDist * lDist, 0.001);
        light += 100 * attenuation * lightData.intensity * lightData.color;
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
        // Multiplicative blending
        throughput = surfaceInfo.color;
        emission = 0;
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

float3 RenderRayOpenRTX(RayDesc rayDesc, out float outputDistance, out float3 outputMotion, OpenRTXContext ctx, float2 pixelCoord)
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
            RenderVanillaOpenRTX(hitInfo, rayState, ctx, pixelCoord);
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
    float3 color = RenderRayOpenRTX(rayDesc, hitDist, objMotion, ctx, float2(dispatchThreadID.xy));

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
