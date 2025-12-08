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

// =============================================================================
// GLASS REFRACTION HELPERS
// =============================================================================

// Wavelength-dependent IOR using Cauchy's equation for dispersion
// Returns IOR for red, green, blue wavelengths
float3 getCauchyIOR(float baseIOR, float dispersionStrength)
{
    // Cauchy coefficients approximation
    // n(λ) = A + B/λ² where λ in micrometers
    // Red ~0.65μm, Green ~0.55μm, Blue ~0.45μm
    float A = baseIOR - dispersionStrength * 0.5;
    float B = dispersionStrength * 0.01;

    float3 wavelengths = float3(0.65, 0.55, 0.45);  // RGB wavelengths in μm
    float3 wavelengthsSq = wavelengths * wavelengths;

    return A + B / wavelengthsSq;
}

// Refract with specific IOR ratio (eta = n1/n2)
// For air to glass: eta = 1.0 / glassIOR
// For glass to air: eta = glassIOR
float3 refractWithEta(float3 incident, float3 normal, float eta, out bool tir)
{
    float cosI = -dot(incident, normal);

    // Ensure normal faces against incident direction
    if (cosI < 0.0)
    {
        normal = -normal;
        cosI = -cosI;
    }

    float sinT2 = eta * eta * (1.0 - cosI * cosI);

    if (sinT2 > 1.0)
    {
        tir = true;
        return reflect(incident, normal);
    }

    tir = false;
    float cosT = sqrt(1.0 - sinT2);
    return eta * incident + (eta * cosI - cosT) * normal;
}

// Convenience wrapper: refract entering glass from air
float3 refractEnterGlass(float3 incident, float3 normal, float glassIOR, out bool tir)
{
    return refractWithEta(incident, normal, 1.0 / glassIOR, tir);
}

// Convenience wrapper: refract exiting glass to air
float3 refractExitGlass(float3 incident, float3 normal, float glassIOR, out bool tir)
{
    return refractWithEta(incident, normal, glassIOR, tir);
}

// Legacy wrapper for compatibility
float3 refractWithIOR(float3 incident, float3 normal, float ior, out bool tir)
{
    return refractEnterGlass(incident, normal, ior, tir);
}

// Calculate Fresnel reflectance for dielectric (glass) at given angle
float glassFresnelReflectance(float cosTheta, float ior)
{
    // Use full Fresnel equations for accurate glass
    float eta = 1.0 / ior;  // Air to glass
    float sinThetaTSq = eta * eta * (1.0 - cosTheta * cosTheta);

    // Total internal reflection check
    if (sinThetaTSq > 1.0)
        return 1.0;

    float cosThetaT = sqrt(max(0.0, 1.0 - sinThetaTSq));

    // Fresnel equations for s and p polarization
    float rs = (eta * cosTheta - cosThetaT) / (eta * cosTheta + cosThetaT);
    float rp = (cosTheta - eta * cosThetaT) / (cosTheta + eta * cosThetaT);

    return (rs * rs + rp * rp) * 0.5;
}

// Calculate refracted direction through glass
float3 glassRefract(float3 incident, float3 normal, float ior, out bool totalInternalReflection)
{
    float eta = 1.0 / ior;  // Air to glass ratio
    float cosI = -dot(incident, normal);

    // Check if we're inside the glass (ray exiting)
    if (cosI < 0.0)
    {
        // Flip normal and use inverse ratio
        normal = -normal;
        cosI = -cosI;
        eta = ior;  // Glass to air
    }

    float sinT2 = eta * eta * (1.0 - cosI * cosI);

    // Total internal reflection
    if (sinT2 > 1.0)
    {
        totalInternalReflection = true;
        return reflect(incident, normal);
    }

    totalInternalReflection = false;
    float cosT = sqrt(1.0 - sinT2);
    return eta * incident + (eta * cosI - cosT) * normal;
}

// Physically accurate glass refraction with proper entry/exit handling
// Returns: color of scene behind glass, attenuated by Beer-Lambert absorption
float3 traceGlassRefractionWithIOR(float3 origin, float3 incidentDir, float3 entryNormal,
                                    float3 glassColor, float glassAlpha, OpenRTXContext ctx,
                                    int maxBounces, float ior)
{
#if !ENABLE_GLASS_REFRACTION
    return 0.0;
#endif

    float3 throughput = 1.0;
    float3 currentOrigin = origin;
    float3 currentDir = incidentDir;
    bool insideGlass = true;  // We start by entering glass
    float currentIOR = ior;
    float totalGlassDistance = 0.0;  // Track distance traveled inside glass

    // Convert glass color to absorption coefficients using Beer-Lambert law
    // For colored glass: absorption = -ln(transmittance) / reference_distance
    // We use the glass color as the transmittance at 1 unit distance
    // Darker colors = higher absorption, brighter colors = lower absorption
    float3 glassTransmittance = saturate(glassColor);
    // Prevent log(0) - minimum transmittance of 1% per unit
    glassTransmittance = max(glassTransmittance, 0.01);
    // Calculate absorption coefficient (higher alpha = more colored/opaque glass)
    float3 absorptionCoeff = -log(glassTransmittance) * lerp(0.5, 2.0, glassAlpha);

    // First, refract the incident ray as it enters the glass (air to glass)
    bool tir;
    currentDir = refractEnterGlass(incidentDir, entryNormal, currentIOR, tir);
    if (tir)
    {
        // Shouldn't happen when entering glass from air, but handle it
        return 0.0;
    }

    // Trace through glass with multiple bounces
    [loop]
    for (int bounce = 0; bounce < maxBounces; bounce++)
    {
        RayQuery<RAY_FLAG_NONE> refractQuery;
        RayDesc refractRay;
        refractRay.Origin = currentOrigin + currentDir * 0.001;
        refractRay.Direction = currentDir;
        refractRay.TMin = 0.0;
        refractRay.TMax = 1000.0;

        refractQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, refractRay);

        while (refractQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(refractQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                refractQuery.CommitNonOpaqueTriangleHit();
            }
        }

        if (refractQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(refractQuery);
            ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
            GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
            SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

            float distance = hitInfo.rayT;

            // Apply Beer-Lambert absorption while traveling through glass
            // T = exp(-absorption_coefficient * distance)
            if (insideGlass)
            {
                float3 absorption = exp(-absorptionCoeff * distance * GLASS_ABSORPTION_SCALE);
                throughput *= absorption;
                totalGlassDistance += distance;
            }

            // Determine if we hit glass or something else
            bool hitGlass = (hitInfo.materialType == MATERIAL_TYPE_ALPHA_BLEND);
            bool hitWater = (hitInfo.materialType == MATERIAL_TYPE_WATER);

            // Only treat as glass if it's alpha blend and not water
            if (hitGlass && !hitWater && insideGlass)
            {
                float3 hitNormal = surfaceInfo.normal;

                // Determine if we're exiting or entering glass based on normal direction
                bool exitingGlass = dot(currentDir, hitNormal) > 0;

                if (exitingGlass)
                {
                    // For very thin glass (like panes), skip exit refraction
                    // This maintains the visual distortion effect
                    if (totalGlassDistance < THIN_GLASS_THRESHOLD)
                    {
                        // Apply minimum absorption for thin glass visibility
                        // Use a virtual minimum thickness for color tinting
                        float virtualThickness = max(totalGlassDistance, GLASS_MIN_THICKNESS);
                        float additionalDist = virtualThickness - totalGlassDistance;
                        if (additionalDist > 0.0)
                        {
                            float3 additionalAbsorption = exp(-absorptionCoeff * additionalDist * GLASS_ABSORPTION_SCALE);
                            throughput *= additionalAbsorption;
                        }

                        // Skip exit refraction for thin glass - maintain distorted direction
                        currentOrigin = surfaceInfo.position;
                        insideGlass = false;
                        continue;
                    }

                    // Exiting glass to air - use proper glass-to-air refraction
                    // Normal should point into the glass (opposite ray direction)
                    float3 exitNormal = -hitNormal;

                    bool exitTir;
                    float3 exitDir = refractExitGlass(currentDir, exitNormal, currentIOR, exitTir);

                    if (exitTir)
                    {
                        // Total internal reflection - reflect and stay inside
                        currentDir = reflect(currentDir, exitNormal);
                        currentOrigin = surfaceInfo.position;
                        continue;
                    }

                    // Successfully exited glass
                    currentDir = exitDir;
                    currentOrigin = surfaceInfo.position;
                    insideGlass = false;
                }
                else
                {
                    // Entering another glass surface (stacked glass)
                    // Update absorption coefficients for the new glass color
                    float3 newGlassTransmittance = max(saturate(surfaceInfo.color), 0.01);
                    absorptionCoeff = -log(newGlassTransmittance) * lerp(0.5, 2.0, surfaceInfo.alpha);

                    bool enterTir;
                    currentDir = refractEnterGlass(currentDir, hitNormal, currentIOR, enterTir);
                    currentOrigin = surfaceInfo.position;
                    insideGlass = true;
                }
                continue;
            }
            else
            {
                // Hit a non-glass surface - shade it and return
                float3 N = surfaceInfo.normal;
                float NdotL = saturate(dot(N, ctx.sunDir));

                // Full shading for the surface behind glass
                float3 surfaceLight = surfaceInfo.color * ctx.sunColor * NdotL;
                surfaceLight += surfaceInfo.color * ctx.constantAmbient;
                surfaceLight += surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;

                return throughput * surfaceLight;
            }
        }
        else
        {
            // Hit sky
            return throughput * renderSkyWithClouds(currentDir, ctx);
        }
    }

    // Ran out of bounces - return what we have
    return throughput * renderSkyWithClouds(currentDir, ctx);
}

// =============================================================================
// CHROMATIC DISPERSION AND FROSTED GLASS
// =============================================================================

// Generate importance-sampled direction for rough glass
float3 sampleFrostedDirection(float3 direction, float3 normal, float roughness, float2 randomSeed)
{
    if (roughness < 0.001)
        return direction;

    // GGX-like sampling for frosted glass
    float2 xi = float2(
        frac(sin(dot(randomSeed, float2(12.9898, 78.233))) * 43758.5453),
        frac(sin(dot(randomSeed * 1.5, float2(127.1, 311.7))) * 43758.5453)
    );

    float a = roughness * roughness;
    float theta = atan(a * sqrt(xi.x) / sqrt(1.0 - xi.x));
    float phi = 2.0 * 3.14159 * xi.y;

    // Microfacet normal in tangent space
    float3 H_tangent = float3(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    );

    // Build tangent space around refracted direction
    float3 up = abs(direction.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
    float3 tangent = normalize(cross(up, direction));
    float3 bitangent = cross(direction, tangent);

    // Transform to world space and blend with original direction
    float3 scatteredDir = normalize(
        H_tangent.x * tangent +
        H_tangent.y * bitangent +
        H_tangent.z * direction
    );

    return scatteredDir;
}

// Trace glass with chromatic dispersion (traces R, G, B separately)
float3 traceGlassWithDispersion(float3 origin, float3 incidentDir, float3 normal,
                                 float3 glassColor, float glassAlpha, OpenRTXContext ctx,
                                 int maxBounces, float roughness, float2 pixelCoord)
{
#if !ENABLE_GLASS_REFRACTION
    return 0.0;
#endif

    // Apply roughness to incident direction if needed (frosted glass)
    float3 direction = incidentDir;
    if (roughness > 0.001)
    {
        direction = sampleFrostedDirection(incidentDir, normal, roughness, pixelCoord);
    }

#if ENABLE_GLASS_DISPERSION
    // Get wavelength-dependent IOR for chromatic dispersion
    float3 iorRGB = getCauchyIOR(GLASS_IOR, GLASS_DISPERSION_STRENGTH);

    float3 resultRGB = 0.0;

    // Trace each wavelength with its own IOR
    float3 colorR = traceGlassRefractionWithIOR(origin, direction, normal, glassColor, glassAlpha, ctx, maxBounces, iorRGB.r);
    float3 colorG = traceGlassRefractionWithIOR(origin, direction, normal, glassColor, glassAlpha, ctx, maxBounces, iorRGB.g);
    float3 colorB = traceGlassRefractionWithIOR(origin, direction, normal, glassColor, glassAlpha, ctx, maxBounces, iorRGB.b);

    resultRGB.r = colorR.r;
    resultRGB.g = colorG.g;
    resultRGB.b = colorB.b;

    return resultRGB;
#else
    // No dispersion - use standard IOR
    return traceGlassRefractionWithIOR(origin, direction, normal, glassColor, glassAlpha, ctx, maxBounces, GLASS_IOR);
#endif
}

// =============================================================================
// MEDIUM TRACKING FOR RAY TRANSITIONS
// =============================================================================

// Medium types for IOR tracking
#define MEDIUM_AIR   0
#define MEDIUM_WATER 1
#define MEDIUM_GLASS 2

// Get IOR for medium transition
float getMediumIOR(int fromMedium, int toMedium)
{
    // IOR values
    static const float IOR_AIR = 1.0;
    static const float IOR_WATER = WATER_IOR;
    static const float IOR_GLASS = GLASS_IOR;

    float fromIOR = IOR_AIR;
    float toIOR = IOR_AIR;

    if (fromMedium == MEDIUM_WATER) fromIOR = IOR_WATER;
    else if (fromMedium == MEDIUM_GLASS) fromIOR = IOR_GLASS;

    if (toMedium == MEDIUM_WATER) toIOR = IOR_WATER;
    else if (toMedium == MEDIUM_GLASS) toIOR = IOR_GLASS;

    return toIOR / fromIOR;
}

// Refract ray through medium boundary
float3 refractMediumBoundary(float3 incident, float3 normal, int fromMedium, int toMedium, out bool tir)
{
    float relativeIOR = getMediumIOR(fromMedium, toMedium);

    float eta = 1.0 / relativeIOR;
    float cosI = -dot(incident, normal);

    // Handle backface
    if (cosI < 0.0)
    {
        normal = -normal;
        cosI = -cosI;
        eta = relativeIOR;
    }

    float sinT2 = eta * eta * (1.0 - cosI * cosI);

    if (sinT2 > 1.0)
    {
        tir = true;
        return reflect(incident, normal);
    }

    tir = false;
    float cosT = sqrt(1.0 - sinT2);
    return eta * incident + (eta * cosI - cosT) * normal;
}

// Fresnel for medium transition
float fresnelMediumBoundary(float cosTheta, int fromMedium, int toMedium)
{
    float relativeIOR = getMediumIOR(fromMedium, toMedium);
    float eta = 1.0 / relativeIOR;
    float sinThetaTSq = eta * eta * (1.0 - cosTheta * cosTheta);

    if (sinThetaTSq > 1.0)
        return 1.0;

    float cosThetaT = sqrt(max(0.0, 1.0 - sinThetaTSq));
    float rs = (eta * cosTheta - cosThetaT) / (eta * cosTheta + cosThetaT);
    float rp = (cosTheta - eta * cosThetaT) / (cosTheta + eta * cosThetaT);

    return (rs * rs + rp * rp) * 0.5;
}

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
            if (hitInfo.materialType == MATERIAL_TYPE_OPAQUE)
            {
                // Opaque surface blocks light completely
                result.visibility = 0.0;
                result.transmission = 0.0;
                return result;
            }
            else if (hitInfo.materialType == MATERIAL_TYPE_ALPHA_TEST)
            {
                // Alpha test materials (leaves, fences, etc.) should block most light
                // but allow minimal transmission to prevent complete darkness
                result.visibility *= 0.1; // Allow 10% light through alpha test
                result.transmission *= 0.05; // Minimal color transmission
                return result;
            }

            // Translucent surface - accumulate color transmission (only for alpha blend materials)
            if (hitInfo.materialType == MATERIAL_TYPE_ALPHA_BLEND)
            {
#if ENABLE_COLORED_SHADOWS
                // Apply Beer's law for color absorption
                float3 surfaceColor = surfaceInfo.color;
                float alpha = surfaceInfo.alpha;

                // Enhanced color transmission through translucent material
                // Use stronger color mixing for better colored light transmission
                float transmissionFactor = alpha * COLOR_TRANSMISSION_STRENGTH;
                float3 transmittedColor = lerp(1.0, surfaceColor, transmissionFactor * 1.5); // Boosted transmission
                colorAccum *= transmittedColor;

                // Reduce visibility based on alpha, but allow more colored light through
                result.visibility *= (1.0 - alpha * 0.3); // Reduced visibility penalty for colored glass
#else
                result.visibility *= (1.0 - surfaceInfo.alpha);
#endif
            }

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
// EMISSIVE LIGHT SAMPLING (Global Illumination from emissive blocks)
// =============================================================================

// Sample indirect light from emissive surfaces
float3 sampleEmissiveLight(float3 origin, float3 normal, float2 pixelCoord, OpenRTXContext ctx, int numSamples)
{
#if !ENABLE_EMISSIVE_GI
    return 0.0;
#endif

    float3 emissiveLight = 0.0;

    // Build tangent space for hemisphere sampling
    float3 tangent = normalize(cross(normal, abs(normal.y) < 0.999 ? float3(0, 1, 0) : float3(1, 0, 0)));
    float3 bitangent = cross(normal, tangent);

    float totalWeight = 0.0;

    [loop]
    for (int i = 0; i < numSamples; i++)
    {
        // Generate random direction in hemisphere (cosine weighted)
        float2 xi = float2(
            frac(sin(dot(pixelCoord + float2(i * 17.0, g_view.time * 50.0), float2(12.9898, 78.233))) * 43758.5453),
            frac(sin(dot(pixelCoord + float2(g_view.time * 50.0, i * 31.0), float2(127.1, 311.7))) * 43758.5453)
        );

        // Cosine-weighted hemisphere sampling
        float r = sqrt(xi.x);
        float theta = 2.0 * kPi * xi.y;
        float x = r * cos(theta);
        float y = r * sin(theta);
        float z = sqrt(max(0.0, 1.0 - xi.x));

        float3 sampleDir = normalize(x * tangent + y * bitangent + z * normal);

        // Trace ray
        RayQuery<RAY_FLAG_NONE> emissiveQuery;
        RayDesc emissiveRay;
        emissiveRay.Origin = origin + normal * 0.01;
        emissiveRay.Direction = sampleDir;
        emissiveRay.TMin = 0.0;
        emissiveRay.TMax = EMISSIVE_GI_RANGE;

        emissiveQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, emissiveRay);

        while (emissiveQuery.Proceed())
        {
            HitInfo hitInfo = GetCandidateHitInfo(emissiveQuery);
            if (AlphaTestHitLogic(hitInfo))
            {
                emissiveQuery.CommitNonOpaqueTriangleHit();
            }
        }

        if (emissiveQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
        {
            HitInfo hitInfo = GetCommittedHitInfo(emissiveQuery);
            ObjectInstance objectInstance = objectInstances[hitInfo.objectInstanceIndex];
            GeometryInfo geometryInfo = GetGeometryInfo(hitInfo, objectInstance);
            SurfaceInfo surfaceInfo = MaterialVanilla(hitInfo, geometryInfo, objectInstance);

            // If we hit an emissive surface, accumulate its light
            if (surfaceInfo.emissive > 0.01)
            {
                float distance = hitInfo.rayT;
                float falloff = 1.0 / (1.0 + distance * distance * EMISSIVE_GI_FALLOFF);

                // Weight by emissive strength and distance
                float3 emitted = surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY;
                emissiveLight += emitted * falloff * EMISSIVE_GI_STRENGTH;
                totalWeight += 1.0;
            }
            else
            {
                // Non-emissive surface - add small ambient bounce
                float3 bounceLight = surfaceInfo.color * ctx.constantAmbient * 0.2;
                emissiveLight += bounceLight;
                totalWeight += 1.0;
            }
        }
    }

    if (totalWeight > 0.0)
    {
        emissiveLight /= totalWeight;
    }

    return emissiveLight;
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

        // Add emissive with indirect boost (allows emissive to cast more GI light)
        // Uses INDIRECT_EMISSIVE_BOOST to make torches/glowstone light up surroundings
        directLight += surfaceInfo.color * surfaceInfo.emissive * EMISSIVE_INTENSITY * INDIRECT_EMISSIVE_BOOST;

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
    // Calculate distance for LOD
    float distanceFromCamera = length(origin - ctx.viewOrigin);
    float distanceLOD = saturate(distanceFromCamera / 256.0);  // LOD starts at ~256 blocks

    // Skip if too rough
    float effectiveRoughness = max(roughness - REFLECTION_ROUGHNESS_BIAS, MIN_ROUGHNESS);

    // Increase roughness threshold at distance (skip more reflections far away)
    float adjustedMaxRoughness = lerp(REFLECTION_MAX_ROUGHNESS, REFLECTION_MAX_ROUGHNESS * 0.5, distanceLOD);

    if (effectiveRoughness > adjustedMaxRoughness || bounceDepth >= REFLECTION_MAX_BOUNCES)
    {
        // Fall back to sky reflection for rough/distant surfaces
        float3 reflectDir = reflect(viewDir, normal);
        return renderSkyWithClouds(reflectDir, ctx) * 0.1;
    }

    // Determine number of samples - reduce at distance for performance
    int maxSamples = int(lerp(float(REFLECTION_SAMPLES), 1.0, distanceLOD));
    int numSamples = maxSamples;
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

    // Compute Fresnel with metal-aware F90
    // For metals, F90 should be tinted rather than pure white to preserve color at grazing angles
    float3 f0 = lerp(0.04, albedo, metalness);
    // F90: For dielectrics = white, for metals = slightly brighter version of F0
    // This preserves metal color at grazing angles instead of going pure white
    float3 f90 = lerp(1.0, saturate(albedo * 1.2 + 0.2), metalness);
    float NdotV = saturate(dot(normal, -viewDir));
    float fresnelPow = pow(1.0 - NdotV, 5.0);
    float3 fresnel = f0 + (f90 - f0) * fresnelPow;

    // Apply fresnel to reflection
    reflectionColor *= fresnel;

    // Fade based on roughness
    float roughnessFade = 1.0 - saturate(effectiveRoughness / REFLECTION_MAX_ROUGHNESS);
    reflectionColor *= roughnessFade;

    // Apply edge-aware filtering for noise reduction (preserving color for metals)
    float smoothFactor = saturate(1.0 - effectiveRoughness * 2.0);
    // For metals, preserve color; for dielectrics, allow some desaturation
    float3 smoothedColor = lerp(reflectionColor * smoothFactor, reflectionColor, metalness);
    reflectionColor = lerp(reflectionColor, smoothedColor, effectiveRoughness);
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

    // Check for cloud geometry
    bool isCloud = (objectInstance.flags & kObjectInstanceFlagClouds) != 0;

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
    // Kill diffuse for emissive surfaces (they emit, not reflect)
    surface.albedo = lerp(surfaceInfo.color, 0.0, surfaceInfo.emissive);
    surface.roughness = max(surfaceInfo.roughness, MIN_ROUGHNESS);
    surface.metalness = surfaceInfo.metalness;
    surface.ao = 1.0;
    surface.subsurface = surfaceInfo.subsurface;
    surface.emissive = surfaceInfo.emissive;
    surface.isWater = isWater;
    surface.waterDepth = 0.0;

    // Clouds should be fully diffuse (no specular reflections)
    // They're volumetric scattering surfaces, not smooth reflective materials
    // This MUST come after the MIN_ROUGHNESS clamp to override it
    if (isCloud)
    {
        surface.roughness = 1.0;  // Fully rough - no specular (overrides MIN_ROUGHNESS)
        surface.metalness = 0.0;  // Non-metallic
    }

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

    // Apply shadow with proper linear falloff to prevent crushing
    // Use direct multiplication instead of lerp for accurate shadow values
    light *= shadow.visibility;

    // Apply color transmission from stained glass etc.
    light *= shadow.transmission;

    // Apply emissive - boost light for emissive surfaces
    light = lerp(light, EMISSIVE_INTENSITY, surfaceInfo.emissive);
#endif

    // =================================================================
    // RAYTRACED REFLECTIONS
    // =================================================================
    float3 reflectionContrib = 0.0;

#if OPENRTX_ENABLED && ENABLE_RAYTRACED_REFLECTIONS
    // Add reflections for smooth/metallic surfaces (not clouds - they're diffuse scattering)
    float reflectivity = surface.metalness + (1.0 - surface.roughness) * 0.3;
    if (reflectivity > 0.1 && isOpaque && !isCloud)
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

    // =================================================================
    // EMISSIVE GI (Light from glowing blocks illuminating surroundings)
    // =================================================================
#if OPENRTX_ENABLED && ENABLE_EMISSIVE_GI
    // Only apply emissive GI to non-emissive opaque surfaces
    if (isOpaque && surfaceInfo.emissive < 0.1)
    {
        float3 emissiveGI = sampleEmissiveLight(
            surfaceInfo.position,
            surfaceInfo.normal,
            pixelCoord,
            ctx,
            EMISSIVE_GI_SAMPLES
        );

        // Add emissive light contribution to diffuse surfaces
        light += surfaceInfo.color * emissiveGI;
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
        // Glass/transparent surface with IOR-based Fresnel
#if ENABLE_GLASS_REFRACTION
        // Check if this is actual glass: must be ALPHA_BLEND material type, not water
        // Opaque and alpha-test materials should never get glass treatment
        bool isTransparentMaterial = (hitInfo.materialType == MATERIAL_TYPE_ALPHA_BLEND);
        bool isGlassLike = isTransparentMaterial && !isWater && surfaceInfo.alpha > 0.1;

        if (isGlassLike)
        {
            float3 viewDir = -rayState.rayDesc.Direction;
            float3 N = surfaceInfo.normal;
            float NdotV = abs(dot(N, viewDir));

            // Calculate Fresnel reflectance using glass IOR
            // Apply strength parameter and cap to prevent excessive washout at grazing angles
            float fresnelRaw = glassFresnelReflectance(NdotV, GLASS_IOR);
            // Calculate proper F0 from glass IOR instead of hardcoded value
            float3 glassF0 = f0FromIOR(GLASS_IOR);
            float fresnelBase = glassF0.r;  // Base reflectance at normal incidence from IOR
            // Interpolate between base and full Fresnel based on strength parameter
            float fresnel = lerp(fresnelBase, fresnelRaw, GLASS_FRESNEL_STRENGTH);
            // Cap maximum reflectance to preserve color transmission
            fresnel = min(fresnel, GLASS_FRESNEL_MAX);

            // For clear glass (low alpha), most light transmits with refraction
            // For tinted glass (high alpha), more absorption/tinting
            float glassOpacity = surfaceInfo.alpha;

            // Check for total internal reflection
            bool tir;
            float3 testDir = glassRefract(rayState.rayDesc.Direction, N, GLASS_IOR, tir);

            // Trace refracted ray to get what's behind the glass
            // Beer-Lambert absorption is applied inside the trace function
            float3 refractedColor = 0.0;
            if (!tir)
            {
                // Use dispersion-aware function (handles frosted glass too)
                // Color absorption is handled by Beer-Lambert law in the trace
                refractedColor = traceGlassWithDispersion(
                    surfaceInfo.position, rayState.rayDesc.Direction, N,
                    surfaceInfo.color, glassOpacity, ctx,
                    GLASS_MAX_BOUNCES, GLASS_ROUGHNESS,
                    pixelCoord);
            }
            else
            {
                // Total internal reflection - use reflection instead
                refractedColor = renderSkyWithClouds(reflect(rayState.rayDesc.Direction, N), ctx);
            }

            // Combine: Fresnel reflection + transmitted/refracted light
            // Apply additional glass color tinting for stronger colored transmission
            float transmissionFactor = 1.0 - fresnel;
            float3 glassTint = lerp(1.0, surfaceInfo.color, surfaceInfo.alpha * 0.8); // Stronger color tint
            float3 transmittedLight = refractedColor * transmissionFactor * glassTint;

            // Trace actual reflection for glass specular (not just sky)
            // This prevents washed-out appearance in enclosed spaces
            float3 reflectDir = reflect(rayState.rayDesc.Direction, N);
            float3 reflectedLight = traceSingleReflection(
                surfaceInfo.position, reflectDir, N, ctx);
            float3 glassSpecular = reflectedLight * fresnel;

            emission = transmittedLight + glassSpecular;

            // Throughput for any additional layers
            throughput = 0.0;  // We've fully handled this with refraction trace
        }
        else
#endif
        // Water surface refraction (looking down into water from above)
#if ENABLE_WATER_REFRACTION
        if (isWater && !g_view.cameraIsUnderWater)
        {
            float3 viewDir = -rayState.rayDesc.Direction;

            // Get enhanced water surface properties with waves
            WaterSurface water = evaluateWaterSurface(
                surfaceInfo.position,
                viewDir,
                ctx.time,
                ctx.rainIntensity
            );

            // Use wave-displaced normal for refraction
            float3 N = water.normal;
            float NdotV = saturate(dot(N, viewDir));

            // Get Fresnel reflectance
            float fresnel = water.fresnel;

            // Calculate refracted ray direction (air to water)
            float3 refractDir = water.refractDir;

            // Apply wave-based distortion strength
            float distortionStrength = WATER_REFRACTION_STRENGTH;

            // Trace refracted ray to see what's underwater
            float3 refractedColor = 0.0;
            float underwaterDepth = 0.0;

            // Create refracted ray
            RayDesc refractRay;
            refractRay.Origin = surfaceInfo.position + refractDir * 0.01;  // Small offset to avoid self-intersection
            refractRay.Direction = refractDir;
            refractRay.TMin = 0.0;
            refractRay.TMax = WATER_REFRACTION_DEPTH_FADE * 2.0;  // Limit underwater trace distance

            // Trace the refracted ray
            RayQuery<RAY_FLAG_NONE> refractQuery;
            refractQuery.TraceRayInline(SceneBVH, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, 0xff, refractRay);

            while (refractQuery.Proceed())
            {
                HitInfo refractHit = GetCandidateHitInfo(refractQuery);
                // Skip water surfaces in refraction trace
                if (refractHit.materialType != MATERIAL_TYPE_WATER && AlphaTestHitLogic(refractHit))
                {
                    refractQuery.CommitNonOpaqueTriangleHit();
                }
            }

            if (refractQuery.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
            {
                HitInfo underwaterHit = GetCommittedHitInfo(refractQuery);
                underwaterDepth = underwaterHit.rayT;

                // Get surface info for underwater object
                ObjectInstance underwaterObj = objectInstances[underwaterHit.objectInstanceIndex];
                GeometryInfo underwaterGeom = GetGeometryInfo(underwaterHit, underwaterObj);
                SurfaceInfo underwaterSurf = MaterialVanilla(underwaterHit, underwaterGeom, underwaterObj);

                if (!underwaterSurf.shouldDiscard)
                {
                    // Simple lighting for underwater surface
                    float3 underwaterNormal = underwaterSurf.normal;
                    float underwaterNdotL = saturate(dot(underwaterNormal, ctx.sunDir));

                    // Underwater receives attenuated sunlight
                    float3 underwaterLight = ctx.sunColor * ctx.sunIntensity * underwaterNdotL * 0.5;
                    underwaterLight += ctx.constantAmbient * 0.3;  // Ambient

                    // Add caustics effect on underwater surfaces
                    float3 caustics = waterCaustics(underwaterSurf.position, ctx.time, ctx.sunDir);
                    underwaterLight += caustics * ctx.sunIntensity * 0.5;

                    refractedColor = underwaterSurf.color * underwaterLight;

                    // Apply water absorption based on depth traveled
                    float3 waterAbsorption = waterTransmittance(underwaterDepth);
                    refractedColor *= waterAbsorption;

                    // Add underwater inscatter (blue tint from water)
                    float3 inscatter = waterInscatter(underwaterDepth, ctx.sunColor, ctx.sunIntensity * 0.3);
                    refractedColor += inscatter;
                }
            }
            else
            {
                // No underwater hit - show deep water color
                underwaterDepth = WATER_REFRACTION_DEPTH_FADE;
                refractedColor = waterInscatter(underwaterDepth, ctx.sunColor, ctx.sunIntensity * 0.5);
            }

            // Fade refraction effect with depth (distant underwater objects less distorted)
            float depthFade = saturate(underwaterDepth / WATER_REFRACTION_DEPTH_FADE);

            // Trace reflection on water surface
            float3 reflectDir = reflect(rayState.rayDesc.Direction, N);
            float3 reflectedColor = renderSkyWithClouds(reflectDir, ctx);

            // Also trace actual scene reflection for nearby objects
            float3 sceneReflection = traceSingleReflection(surfaceInfo.position, reflectDir, N, ctx);
            reflectedColor = lerp(sceneReflection, reflectedColor, 0.3);  // Blend scene and sky reflection

            // Combine reflection and refraction using Fresnel
            emission = lerp(refractedColor, reflectedColor, fresnel);

            // Water is fully handled
            throughput = 0.0;
        }
        else
#endif
        {
            // Standard alpha blend for very transparent surfaces
            throughput = 1 - surfaceInfo.alpha;
            emission = surfaceInfo.color * surfaceInfo.alpha * light;
        }
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

    // Determine distance for outputs and volumetric effects
    float effectiveDistance;
    if (all(rayState.throughput == 0))
    {
        // Hit opaque surface - use actual hit distance
        outputDistance = min(rayState.distance, maxDistance);
        outputMotion = rayState.motion;
        effectiveDistance = outputDistance;
    }
    else
    {
        // Sky ray or partially transparent - use far distance for proper fog
        outputDistance = maxDistance;
        outputMotion = 0;
        effectiveDistance = maxDistance;
    }

    // Render sky
    RenderSkyOpenRTX(rayState, ctx);

    // Apply volumetric effects using effective distance (not rayState.distance which may be 0)
#if OPENRTX_ENABLED && ENABLE_VOLUMETRIC_LIGHTING
    float3 finalColor = applyAtmosphericEffects(rayState.color, rayDesc.Direction, effectiveDistance, ctx);
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

    // Apply underwater camera distortion when submerged
    float2 screenUV = (float2(dispatchThreadID.xy) + 0.5) * g_view.recipRenderResolution;
    bool isUnderwater = g_view.cameraIsUnderWater != 0;

#if ENABLE_UNDERWATER_DISTORTION
    if (isUnderwater)
    {
        // Apply IOR-based ray distortion for underwater view
        rayDesc.Direction = applyUnderwaterRayDistortion(rayDesc.Direction, screenUV, g_view.time);
    }
#endif

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
