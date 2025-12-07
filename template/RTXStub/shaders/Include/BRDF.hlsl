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
// OpenRTX BRDF - Physically-Based Bidirectional Reflectance Distribution Functions
// Implements multiple diffuse and specular models with energy conservation
// =============================================================================

#ifndef __OPENRTX_BRDF_HLSL__
#define __OPENRTX_BRDF_HLSL__

#include "Settings.hlsl"

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

// Clamp dot product to avoid negative values and numerical issues
float saturateDot(float3 a, float3 b)
{
    return saturate(dot(a, b));
}

// Safe normalization
float3 safeNormalize(float3 v)
{
    float len = length(v);
    return len > 0.0 ? v / len : float3(0, 0, 1);
}

// Convert roughness to alpha (squared roughness)
float roughnessToAlpha(float roughness)
{
    float r = max(roughness, MIN_ROUGHNESS);
    return r * r;
}

// Remapping function for Disney roughness
float disneyRoughnessToAlpha(float roughness)
{
    float r = max(roughness, MIN_ROUGHNESS);
    return r * r;
}

// =============================================================================
// FRESNEL FUNCTIONS
// =============================================================================

// Schlick's approximation for Fresnel
float3 fresnelSchlick(float cosTheta, float3 f0)
{
    float t = 1.0 - cosTheta;
    float t2 = t * t;
    float t5 = t2 * t2 * t;
    return f0 + (1.0 - f0) * t5;
}

// Schlick's approximation with roughness for IBL
float3 fresnelSchlickRoughness(float cosTheta, float3 f0, float roughness)
{
    float t = 1.0 - cosTheta;
    float t2 = t * t;
    float t5 = t2 * t2 * t;
    float3 fr = max(float3(1.0 - roughness, 1.0 - roughness, 1.0 - roughness), f0);
    return f0 + (fr - f0) * t5;
}

// Full Fresnel for dielectrics
float fresnelDielectric(float cosThetaI, float eta)
{
    float sinThetaTSq = eta * eta * (1.0 - cosThetaI * cosThetaI);

    if (sinThetaTSq > 1.0)
        return 1.0; // Total internal reflection

    float cosThetaT = sqrt(max(0.0, 1.0 - sinThetaTSq));

    float rs = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);
    float rp = (cosThetaI - eta * cosThetaT) / (cosThetaI + eta * cosThetaT);

    return (rs * rs + rp * rp) * 0.5;
}

// Get F0 for a given IOR
float3 f0FromIOR(float ior)
{
    float f = (ior - 1.0) / (ior + 1.0);
    return float3(f * f, f * f, f * f);
}

// Get F0 for metals using base color
float3 f0ForMetal(float3 baseColor)
{
    return baseColor;
}

// Blend F0 between dielectric and metal
float3 computeF0(float3 baseColor, float metalness, float ior)
{
    float3 dielectricF0 = f0FromIOR(ior);
    return lerp(dielectricF0, baseColor, metalness);
}

// =============================================================================
// NORMAL DISTRIBUTION FUNCTIONS (NDF)
// =============================================================================

// GGX/Trowbridge-Reitz NDF
float D_GGX(float NdotH, float alpha)
{
    float a2 = alpha * alpha;
    float NdotH2 = NdotH * NdotH;

    float denom = NdotH2 * (a2 - 1.0) + 1.0;
    denom = kPi * denom * denom;

    return a2 / max(denom, 1e-7);
}

// Beckmann NDF
float D_Beckmann(float NdotH, float alpha)
{
    float a2 = alpha * alpha;
    float NdotH2 = NdotH * NdotH;

    float exponent = (NdotH2 - 1.0) / (a2 * NdotH2);
    return exp(exponent) / (kPi * a2 * NdotH2 * NdotH2);
}

// Blinn-Phong NDF (legacy)
float D_BlinnPhong(float NdotH, float alpha)
{
    float shininess = 2.0 / (alpha * alpha) - 2.0;
    return (shininess + 2.0) * kInvTwoPi * pow(max(NdotH, 0.0), shininess);
}

// Select NDF based on settings
float D_Specular(float NdotH, float alpha)
{
#if SPECULAR_DISTRIBUTION == 0
    return D_GGX(NdotH, alpha);
#elif SPECULAR_DISTRIBUTION == 1
    return D_Beckmann(NdotH, alpha);
#else
    return D_BlinnPhong(NdotH, alpha);
#endif
}

// =============================================================================
// GEOMETRY/VISIBILITY FUNCTIONS
// =============================================================================

// Smith GGX single-directional shadowing
float G1_GGX(float NdotV, float alpha)
{
    float a2 = alpha * alpha;
    float NdotV2 = NdotV * NdotV;

    return 2.0 * NdotV / (NdotV + sqrt(a2 + (1.0 - a2) * NdotV2));
}

// Smith GGX separable (uncorrelated)
float G_SmithGGX(float NdotV, float NdotL, float alpha)
{
    return G1_GGX(NdotV, alpha) * G1_GGX(NdotL, alpha);
}

// Smith GGX Height-Correlated (more accurate)
float G_SmithGGXCorrelated(float NdotV, float NdotL, float alpha)
{
    float a2 = alpha * alpha;

    float GGXV = NdotL * sqrt(NdotV * NdotV * (1.0 - a2) + a2);
    float GGXL = NdotV * sqrt(NdotL * NdotL * (1.0 - a2) + a2);

    float GGX = GGXV + GGXL;
    return GGX > 0.0 ? 0.5 / GGX : 0.0;
}

// Fast approximation of height-correlated Smith
float G_SmithGGXCorrelatedFast(float NdotV, float NdotL, float alpha)
{
    float a = alpha;
    float GGXV = NdotL * (NdotV * (1.0 - a) + a);
    float GGXL = NdotV * (NdotL * (1.0 - a) + a);
    return 0.5 / max(GGXV + GGXL, 1e-7);
}

// Kelemen geometry term
float G_Kelemen(float VdotH)
{
    return 1.0 / (VdotH * VdotH);
}

// Select geometry function based on settings
float G_Specular(float NdotV, float NdotL, float alpha, float VdotH)
{
#if GEOMETRY_TERM == 0
    return G_SmithGGX(NdotV, NdotL, alpha);
#elif GEOMETRY_TERM == 1
    return G_SmithGGXCorrelated(NdotV, NdotL, alpha);
#else
    return G_Kelemen(VdotH) * NdotL * NdotV;
#endif
}

// Visibility function (G / (4 * NdotV * NdotL))
float V_SmithGGXCorrelated(float NdotV, float NdotL, float alpha)
{
    float a2 = alpha * alpha;

    float GGXV = NdotL * sqrt(NdotV * NdotV * (1.0 - a2) + a2);
    float GGXL = NdotV * sqrt(NdotL * NdotL * (1.0 - a2) + a2);

    return 0.5 / max(GGXV + GGXL, 1e-7);
}

// =============================================================================
// DIFFUSE BRDF MODELS
// =============================================================================

// Lambertian diffuse (simplest)
float3 diffuseLambertian(float3 albedo)
{
    return albedo * kInvPi;
}

// Disney/Burley diffuse
float3 diffuseBurley(float3 albedo, float roughness, float NdotV, float NdotL, float VdotH)
{
    float FD90 = 0.5 + 2.0 * roughness * VdotH * VdotH;

    float FdV = 1.0 + (FD90 - 1.0) * pow(1.0 - NdotV, 5.0);
    float FdL = 1.0 + (FD90 - 1.0) * pow(1.0 - NdotL, 5.0);

    return albedo * kInvPi * FdV * FdL;
}

// Frostbite energy-conserving Burley diffuse
float3 diffuseFrostbiteBurley(float3 albedo, float roughness, float NdotV, float NdotL, float VdotH)
{
    float energyBias = lerp(0.0, 0.5, roughness);
    float energyFactor = lerp(1.0, 1.0 / 1.51, roughness);

    float FD90 = energyBias + 2.0 * roughness * VdotH * VdotH;

    float FdV = 1.0 + (FD90 - 1.0) * pow(1.0 - NdotV, 5.0);
    float FdL = 1.0 + (FD90 - 1.0) * pow(1.0 - NdotL, 5.0);

    return albedo * kInvPi * FdV * FdL * energyFactor;
}

// Oren-Nayar diffuse (good for rough surfaces)
float3 diffuseOrenNayar(float3 albedo, float roughness, float NdotV, float NdotL, float3 N, float3 V, float3 L)
{
    float sigma2 = roughness * roughness;

    float A = 1.0 - 0.5 * sigma2 / (sigma2 + 0.33);
    float B = 0.45 * sigma2 / (sigma2 + 0.09);

    // Compute cos(phi_i - phi_o)
    // Guard against zero-length vectors when V or L is aligned with N
    float3 VtRaw = V - N * NdotV;
    float3 LtRaw = L - N * NdotL;
    float VtLen = length(VtRaw);
    float LtLen = length(LtRaw);

    // When vectors are too short, cosPhi contribution is negligible
    float cosPhi = 0.0;
    if (VtLen > 0.001 && LtLen > 0.001)
    {
        cosPhi = max(0.0, dot(VtRaw / VtLen, LtRaw / LtLen));
    }

    float thetaI = acos(saturate(NdotL));
    float thetaR = acos(saturate(NdotV));
    float alpha = max(thetaI, thetaR);
    float beta = min(thetaI, thetaR);

    return albedo * kInvPi * (A + B * cosPhi * sin(alpha) * tan(beta));
}

// Select diffuse BRDF based on settings
float3 evaluateDiffuseBRDF(float3 albedo, float roughness, float NdotV, float NdotL, float VdotH, float3 N, float3 V, float3 L)
{
#if DIFFUSE_BRDF == 0
    return diffuseLambertian(albedo);
#elif DIFFUSE_BRDF == 1
    return diffuseBurley(albedo, roughness, NdotV, NdotL, VdotH);
#elif DIFFUSE_BRDF == 2
    return diffuseFrostbiteBurley(albedo, roughness, NdotV, NdotL, VdotH);
#else
    return diffuseOrenNayar(albedo, roughness, NdotV, NdotL, N, V, L);
#endif
}

// =============================================================================
// SPECULAR BRDF
// =============================================================================

// Full specular BRDF evaluation
float3 evaluateSpecularBRDF(float3 f0, float roughness, float NdotV, float NdotL, float NdotH, float VdotH)
{
    float alpha = roughnessToAlpha(roughness);

    float D = D_Specular(NdotH, alpha);
    float3 F = fresnelSchlick(VdotH, f0);
    float V = V_SmithGGXCorrelated(NdotV, NdotL, alpha);

    return D * F * V;
}

// Specular BRDF with separate geometry function
float3 evaluateSpecularBRDFSeparate(float3 f0, float roughness, float NdotV, float NdotL, float NdotH, float VdotH)
{
    float alpha = roughnessToAlpha(roughness);

    float D = D_Specular(NdotH, alpha);
    float3 F = fresnelSchlick(VdotH, f0);
    float G = G_Specular(NdotV, NdotL, alpha, VdotH);

    float denom = 4.0 * NdotV * NdotL;
    return denom > 0.0 ? (D * F * G) / denom : 0.0;
}

// =============================================================================
// MULTI-SCATTER ENERGY COMPENSATION
// =============================================================================

// Approximation for multi-scatter energy
float3 computeMultiScatterEnergy(float3 f0, float roughness, float NdotV)
{
#if ENABLE_MULTISCATTER_COMPENSATION
    // Energy loss due to multiple bounces
    float3 FssEss = fresnelSchlickRoughness(NdotV, f0, roughness);

    // Average Fresnel
    float Ess = 1.0 - roughness; // Simplified approximation
    float3 Fav = f0 + (1.0 - f0) / 21.0;

    // Multi-scatter contribution
    float3 FmsEms = (1.0 - FssEss) * FssEss * Fav / (1.0 - Fav * (1.0 - Ess));

    return FssEss + FmsEms;
#else
    return fresnelSchlickRoughness(NdotV, f0, roughness);
#endif
}

// =============================================================================
// SUBSURFACE SCATTERING
// =============================================================================

// Cheap subsurface approximation (wrap lighting)
float3 subsurfaceWrap(float3 albedo, float NdotL, float subsurface, float3 subsurfaceColor)
{
    float wrap = 0.5;
    float wrappedNdotL = saturate((NdotL + wrap) / (1.0 + wrap));
    float scatter = wrappedNdotL * subsurface;

    return albedo * subsurfaceColor * scatter;
}

// Disney subsurface approximation
float3 subsurfaceDisney(float3 albedo, float roughness, float NdotV, float NdotL, float VdotH, float subsurface)
{
    float FL = pow(1.0 - NdotL, 5.0);
    float FV = pow(1.0 - NdotV, 5.0);

    float Fss90 = VdotH * VdotH * roughness;
    float Fss = lerp(1.0, Fss90, FL) * lerp(1.0, Fss90, FV);
    // Guard against division by zero at grazing angles
    float ss = 1.25 * (Fss * (1.0 / max(NdotL + NdotV, 0.001) - 0.5) + 0.5);

    return albedo * kInvPi * lerp(1.0, ss, subsurface);
}

// =============================================================================
// COMPLETE BRDF EVALUATION
// =============================================================================

struct BRDFInput
{
    float3 albedo;
    float3 normal;
    float roughness;
    float metalness;
    float subsurface;
    float ao;
    float3 viewDir;
    float3 lightDir;
    float ior;
};

struct BRDFOutput
{
    float3 diffuse;
    float3 specular;
    float3 total;
};

BRDFOutput evaluateBRDF(BRDFInput input)
{
    BRDFOutput output;
    output.diffuse = 0;
    output.specular = 0;
    output.total = 0;

    float3 N = input.normal;
    float3 V = input.viewDir;
    float3 L = input.lightDir;
    float3 H = safeNormalize(V + L);

    float NdotV = saturateDot(N, V);
    float NdotL = saturateDot(N, L);
    float NdotH = saturateDot(N, H);
    float VdotH = saturateDot(V, H);

    // Early out for back-facing lights
    if (NdotL <= 0.0)
        return output;

    // Compute F0
    float3 f0 = computeF0(input.albedo, input.metalness, input.ior);

    // Fresnel at normal incidence for diffuse energy conservation
    float3 F = fresnelSchlick(VdotH, f0);

    // Diffuse contribution (metals have no diffuse)
    float3 diffuseAlbedo = input.albedo * (1.0 - input.metalness);

#if ENABLE_SUBSURFACE_SCATTERING
    if (input.subsurface > 0.0)
    {
        output.diffuse = subsurfaceDisney(diffuseAlbedo, input.roughness, NdotV, NdotL, VdotH, input.subsurface);
    }
    else
#endif
    {
        output.diffuse = evaluateDiffuseBRDF(diffuseAlbedo, input.roughness, NdotV, NdotL, VdotH, N, V, L);
    }

    // Energy conservation: diffuse is reduced by specular reflection
    output.diffuse *= (1.0 - F);

    // Specular contribution
    output.specular = evaluateSpecularBRDF(f0, input.roughness, NdotV, NdotL, NdotH, VdotH);

    // Apply NdotL
    output.diffuse *= NdotL;
    output.specular *= NdotL;

    // Ambient occlusion affects diffuse
    output.diffuse *= input.ao;

#if ENABLE_REFLECTION_OCCLUSION
    // Specular occlusion (approximate)
    float specOcclusion = saturate(pow(NdotV + input.ao, input.roughness * input.roughness) - 1.0 + input.ao);
    output.specular *= specOcclusion;
#endif

    output.total = output.diffuse + output.specular;

    return output;
}

// =============================================================================
// ENVIRONMENT/IBL BRDF
// =============================================================================

// Pre-integrated BRDF for split-sum approximation
float2 integrateBRDF(float NdotV, float roughness)
{
    // Analytical approximation of the DFG LUT
    // Based on: https://www.unrealengine.com/blog/physically-based-shading-on-mobile

    float4 c0 = float4(-1.0, -0.0275, -0.572, 0.022);
    float4 c1 = float4(1.0, 0.0425, 1.04, -0.04);

    float4 r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;

    return float2(-1.04, 1.04) * a004 + r.zw;
}

// Environment BRDF evaluation
float3 evaluateEnvironmentBRDF(float3 f0, float roughness, float NdotV, float3 irradiance, float3 prefilteredColor, float ao)
{
    float2 brdf = integrateBRDF(NdotV, roughness);

    float3 Fr = max(float3(1.0 - roughness, 1.0 - roughness, 1.0 - roughness), f0) - f0;
    float3 kS = f0 + Fr * pow(1.0 - NdotV, 5.0);

#if ENABLE_MULTISCATTER_COMPENSATION
    float3 FssEss = kS * brdf.x + brdf.y;
    float Ess = brdf.x + brdf.y;
    float Ems = 1.0 - Ess;
    float3 Fav = f0 + (1.0 - f0) / 21.0;
    float3 FmsEms = Ems * FssEss * Fav / (1.0 - Fav * Ems);
    float3 specular = prefilteredColor * (FssEss + FmsEms);
#else
    float3 specular = prefilteredColor * (kS * brdf.x + brdf.y);
#endif

    float3 kD = (1.0 - kS) * (1.0 - 0.0); // Assuming non-metal for diffuse
    float3 diffuse = kD * irradiance;

#if ENABLE_REFLECTION_OCCLUSION
    float specOcclusion = saturate(pow(NdotV + ao, roughness * roughness) - 1.0 + ao);
    specular *= specOcclusion;
#endif

    return (diffuse + specular) * ao;
}

// =============================================================================
// UTILITY FUNCTIONS FOR MATERIALS
// =============================================================================

// Remap perceptual roughness to linear
float perceptualToLinearRoughness(float perceptualRoughness)
{
    return perceptualRoughness * perceptualRoughness;
}

// Clamp roughness to avoid numerical issues
float clampRoughness(float roughness)
{
    return clamp(roughness, MIN_ROUGHNESS, MAX_ROUGHNESS);
}

// Compute dominant reflection direction for anisotropic materials
float3 getDominantReflectionDir(float3 N, float3 R, float roughness)
{
    float s = 1.0 - roughness;
    return lerp(N, R, s * (sqrt(s) + roughness));
}

// Compute reflection vector
float3 computeReflection(float3 V, float3 N)
{
    return reflect(-V, N);
}

// Compute refraction vector
float3 computeRefraction(float3 V, float3 N, float eta)
{
    return refract(-V, N, eta);
}

#endif // __OPENRTX_BRDF_HLSL__
