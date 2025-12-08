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
// OpenRTX Settings - Comprehensive shader configuration
// All settings use #ifndef guards for external overrides
// =============================================================================

#ifndef __OPENRTX_SETTINGS_HLSL__
#define __OPENRTX_SETTINGS_HLSL__

// =============================================================================
// FEATURE TOGGLES
// =============================================================================

// Master toggle for advanced features (set to 0 for vanilla-like rendering)
#ifndef OPENRTX_ENABLED
#define OPENRTX_ENABLED 1
#endif

// -----------------------------------------------------------------------------
// Rendering Features
// -----------------------------------------------------------------------------

#ifndef ENABLE_PBR_LIGHTING
#define ENABLE_PBR_LIGHTING 1  // Enable physically-based BRDF
#endif

#ifndef ENABLE_ATMOSPHERIC_SKY
#define ENABLE_ATMOSPHERIC_SKY 1  // Rayleigh/Mie atmospheric scattering
#endif

#ifndef ENABLE_VOLUMETRIC_CLOUDS
#define ENABLE_VOLUMETRIC_CLOUDS 0  // Disabled: use vanilla mesh with puffy enhancement instead
#endif

#ifndef ENABLE_VOLUMETRIC_LIGHTING
#define ENABLE_VOLUMETRIC_LIGHTING 1  // Volumetric fog and light shafts
#endif

#ifndef ENABLE_WATER_EFFECTS
#define ENABLE_WATER_EFFECTS 1  // Advanced water rendering with waves
#endif

#ifndef ENABLE_SUBSURFACE_SCATTERING
#define ENABLE_SUBSURFACE_SCATTERING 1  // SSS for translucent materials
#endif

#ifndef ENABLE_CAUSTICS
#define ENABLE_CAUSTICS 1  // Water caustics
#endif

#ifndef ENABLE_RAIN_WETNESS
#define ENABLE_RAIN_WETNESS 1  // Rain wetness and puddles
#endif

#ifndef ENABLE_SNOW_COVER
#define ENABLE_SNOW_COVER 0  // Dynamic snow accumulation (experimental)
#endif

#ifndef ENABLE_POST_PROCESSING
#define ENABLE_POST_PROCESSING 1  // Post-processing effects
#endif

#ifndef ENABLE_FLASHLIGHT
#define ENABLE_FLASHLIGHT 0  // Player-attached flashlight
#endif

// -----------------------------------------------------------------------------
// Bug Fixes (all enabled by default)
// -----------------------------------------------------------------------------

#ifndef FIX_DEFAULT_MATERIAL
#define FIX_DEFAULT_MATERIAL 1  // 100% rough default instead of 0% (vanilla: 0)
#endif

#ifndef FIX_BACKFACE_CULLING
#define FIX_BACKFACE_CULLING 1  // Skip back faces of opaque blocks
#endif

#ifndef FIX_BLEND_SURFACES
#define FIX_BLEND_SURFACES 1  // Skip fully transparent blend objects
#endif

#ifndef FIX_BANNER_UVS
#define FIX_BANNER_UVS 1  // Correct 1px banner UV offset (vanilla: 1)
#endif

#ifndef FIX_DENOISER_PARAMETERS
#define FIX_DENOISER_PARAMETERS 1  // Bounds check denoiser indices
#endif

#ifndef FIX_ITEM_ALBEDO_DARKENING
#define FIX_ITEM_ALBEDO_DARKENING 1  // Bow/crossbow fix
#endif

#ifndef FIX_TINT_SHADED_SURFACES
#define FIX_TINT_SHADED_SURFACES 1  // Remove baked vanilla shading (vanilla: 1)
#endif

#ifndef FIX_INDIRECT_BLEND_EMISSIVES
#define FIX_INDIRECT_BLEND_EMISSIVES 1  // Alpha blend emissive contribution
#endif

#ifndef FIX_SHADOW_ALIASING
#define FIX_SHADOW_ALIASING 1  // Disable pre-interleave shadow blend
#endif

#ifndef FIX_BLOCK_BREAKING_OVERLAY
#define FIX_BLOCK_BREAKING_OVERLAY 1  // Multiplicative blend
#endif

#ifndef FIX_TRANSLUCENT_EXPOSURE
#define FIX_TRANSLUCENT_EXPOSURE 1  // Exposure through translucents
#endif

#ifndef FIX_ENTITY_EMISSIVE
#define FIX_ENTITY_EMISSIVE 1  // Detect emissive entities via alpha
#endif

// Emissive intensity multiplier (kEmissiveBrightnessDirect equivalent)
// Controls how bright emissive surfaces appear in direct rendering
#ifndef EMISSIVE_INTENSITY
#define EMISSIVE_INTENSITY 5.0  // Default: 5.0 (matches BetterRTX)
#endif

// Indirect emissive boost for GI contribution
// Allows emissive surfaces to cast more light than they visually appear
#ifndef INDIRECT_EMISSIVE_BOOST
#define INDIRECT_EMISSIVE_BOOST 2.0  // GI boost multiplier
#endif

// =============================================================================
// EMISSIVE GI SETTINGS (Light bouncing from glowing blocks)
// =============================================================================

// Enable raytraced emissive global illumination
#ifndef ENABLE_EMISSIVE_GI
#define ENABLE_EMISSIVE_GI 1
#endif

// Number of samples for emissive GI (higher = better quality, slower)
#ifndef EMISSIVE_GI_SAMPLES
#define EMISSIVE_GI_SAMPLES 4
#endif

// Maximum range for emissive light sampling (in blocks)
#ifndef EMISSIVE_GI_RANGE
#define EMISSIVE_GI_RANGE 32.0
#endif

// Falloff rate for emissive light (higher = faster falloff)
#ifndef EMISSIVE_GI_FALLOFF
#define EMISSIVE_GI_FALLOFF 0.05
#endif

// Strength multiplier for emissive GI
#ifndef EMISSIVE_GI_STRENGTH
#define EMISSIVE_GI_STRENGTH 1.5
#endif

// =============================================================================
// BRDF SETTINGS
// =============================================================================

// Diffuse BRDF selection
// 0: Lambertian (fastest, least realistic)
// 1: Burley (Disney, good balance)
// 2: Frostbite Burley (energy conserving)
// 3: Oren-Nayar (good for rough surfaces)
#ifndef DIFFUSE_BRDF
#define DIFFUSE_BRDF 2  // Frostbite Burley default
#endif

// Specular distribution function
// 0: GGX (industry standard)
// 1: Beckmann
// 2: Blinn-Phong (legacy)
#ifndef SPECULAR_DISTRIBUTION
#define SPECULAR_DISTRIBUTION 0  // GGX default
#endif

// Fresnel approximation
// 0: Schlick (standard)
// 1: Schlick with roughness adjustment
#ifndef FRESNEL_MODEL
#define FRESNEL_MODEL 1
#endif

// Geometry/visibility term
// 0: Smith GGX
// 1: Smith GGX Height-Correlated
// 2: Kelemen
#ifndef GEOMETRY_TERM
#define GEOMETRY_TERM 1  // Smith GGX Height-Correlated
#endif

// Multi-scatter energy compensation
#ifndef ENABLE_MULTISCATTER_COMPENSATION
#define ENABLE_MULTISCATTER_COMPENSATION 1
#endif

// =============================================================================
// MATERIAL SETTINGS
// =============================================================================

// Default roughness for surfaces without PBR data (vanilla uses 0.0)
#ifndef DEFAULT_ROUGHNESS
#define DEFAULT_ROUGHNESS 0.5  // More realistic default
#endif

// Minimum roughness to avoid specular aliasing
#ifndef MIN_ROUGHNESS
#define MIN_ROUGHNESS 0.045
#endif

// Maximum roughness
#ifndef MAX_ROUGHNESS
#define MAX_ROUGHNESS 1.0
#endif

// Ambient occlusion settings
#ifndef AO_STRENGTH
#define AO_STRENGTH 1.0
#endif

#ifndef ENABLE_BENT_NORMALS
#define ENABLE_BENT_NORMALS 1  // Improved indirect lighting
#endif

#ifndef ENABLE_MULTI_BOUNCE_AO
#define ENABLE_MULTI_BOUNCE_AO 1  // AO approximation for colored surfaces
#endif

#ifndef ENABLE_REFLECTION_OCCLUSION
#define ENABLE_REFLECTION_OCCLUSION 1  // Occlusion for specular reflections
#endif

// =============================================================================
// SUN/MOON SETTINGS
// =============================================================================

#ifndef SUN_RADIUS_MULTIPLIER
#define SUN_RADIUS_MULTIPLIER 1.0  // Sun angular size (vanilla: 1.0)
#endif

#ifndef MOON_RADIUS_MULTIPLIER
#define MOON_RADIUS_MULTIPLIER 1.0  // Moon angular size
#endif

#ifndef SUN_INTENSITY
#define SUN_INTENSITY 100.0  // Sun luminance (vanilla: ~100)
#endif

#ifndef MOON_INTENSITY
#define MOON_INTENSITY 0.5  // Moon luminance (vanilla: ~0.5)
#endif

// Sun/Moon color temperature (Kelvin)
#ifndef SUN_COLOR_TEMPERATURE
#define SUN_COLOR_TEMPERATURE 5778.0  // Actual sun temperature
#endif

#ifndef SUNRISE_COLOR_TEMPERATURE
#define SUNRISE_COLOR_TEMPERATURE 3000.0  // Warm sunrise
#endif

#ifndef SUNSET_COLOR_TEMPERATURE
#define SUNSET_COLOR_TEMPERATURE 2500.0  // Warm sunset
#endif

#ifndef MOON_COLOR_TEMPERATURE
#define MOON_COLOR_TEMPERATURE 4100.0  // Slightly cool moonlight
#endif

// Custom sun angles (set to -1 to use vanilla)
#ifndef CUSTOM_SUN_ZENITH
#define CUSTOM_SUN_ZENITH -1.0  // -1 = use vanilla
#endif

#ifndef CUSTOM_SUN_AZIMUTH
#define CUSTOM_SUN_AZIMUTH -1.0  // -1 = use vanilla
#endif

// =============================================================================
// ATMOSPHERIC SKY SETTINGS
// =============================================================================

// Planet parameters
#ifndef EARTH_RADIUS
#define EARTH_RADIUS 6360000.0  // meters (6360 km)
#endif

#ifndef ATMOSPHERE_RADIUS
#define ATMOSPHERE_RADIUS 6420000.0  // meters (6420 km)
#endif

// Scale heights (in meters)
#ifndef RAYLEIGH_SCALE_HEIGHT
#define RAYLEIGH_SCALE_HEIGHT 7994.0
#endif

#ifndef MIE_SCALE_HEIGHT
#define MIE_SCALE_HEIGHT 1200.0
#endif

// Scattering coefficients at sea level
#ifndef RAYLEIGH_SCATTERING_R
#define RAYLEIGH_SCATTERING_R 5.5e-6
#endif

#ifndef RAYLEIGH_SCATTERING_G
#define RAYLEIGH_SCATTERING_G 13.0e-6
#endif

#ifndef RAYLEIGH_SCATTERING_B
#define RAYLEIGH_SCATTERING_B 22.4e-6
#endif

#ifndef MIE_SCATTERING
#define MIE_SCATTERING 21.0e-6
#endif

#ifndef MIE_ABSORPTION
#define MIE_ABSORPTION 4.4e-6
#endif

// Mie phase function asymmetry (g parameter)
#ifndef MIE_ASYMMETRY
#define MIE_ASYMMETRY 0.76  // Forward scattering bias (vanilla: ~0.8)
#endif

// Ray marching quality settings
#ifndef RAYLEIGH_PRIMARY_INTEGRAL_STEPS
#define RAYLEIGH_PRIMARY_INTEGRAL_STEPS 16  // Primary ray samples (8-32, vanilla: 16)
#endif

#ifndef RAYLEIGH_LIGHT_INTEGRAL_STEPS
#define RAYLEIGH_LIGHT_INTEGRAL_STEPS 0  // Light ray samples (0 = fast approximation, 4-8 for quality)
#endif

// Night sky
#ifndef NIGHT_SKY_INTENSITY
#define NIGHT_SKY_INTENSITY 0.002  // Star visibility
#endif

#ifndef NIGHT_SKY_COLOR
#define NIGHT_SKY_COLOR float3(0.04, 0.05, 0.09)
#endif

// =============================================================================
// VOLUMETRIC CLOUD SETTINGS
// =============================================================================

// Layer altitudes (meters relative to sea level, Minecraft Y=0)
#ifndef CLOUD_MIN_HEIGHT
#define CLOUD_MIN_HEIGHT 1500.0  // Cumulus base (vanilla clouds at Y=192)
#endif

#ifndef CLOUD_MAX_HEIGHT
#define CLOUD_MAX_HEIGHT 4000.0  // Cumulus top
#endif

#ifndef CIRRUS_HEIGHT
#define CIRRUS_HEIGHT 8000.0  // Cirrus layer
#endif

// Cloud density and coverage
#ifndef CLOUD_COVERAGE
#define CLOUD_COVERAGE 0.5  // 0-1, overall cloud coverage
#endif

#ifndef CLOUD_DENSITY
#define CLOUD_DENSITY 0.08  // Optical density multiplier
#endif

#ifndef CLOUD_DENSITY_MULTIPLIER
#define CLOUD_DENSITY_MULTIPLIER 2.5  // Cloud thickness multiplier
#endif

#ifndef CLOUD_DETAIL_STRENGTH
#define CLOUD_DETAIL_STRENGTH 0.35  // Worley noise erosion
#endif

#ifndef CLOUD_SHADOW_STRENGTH
#define CLOUD_SHADOW_STRENGTH 0.7  // Terrain shadow intensity from clouds
#endif

#ifndef CIRRUS_OPACITY
#define CIRRUS_OPACITY 0.25  // High cloud layer visibility
#endif

// Ray marching settings
#ifndef CLOUD_MARCH_STEPS
#define CLOUD_MARCH_STEPS 64  // Quality vs performance
#endif

#ifndef CLOUD_LIGHT_MARCH_STEPS
#define CLOUD_LIGHT_MARCH_STEPS 6  // Light sampling steps
#endif

// Cloud animation
#ifndef CLOUD_WIND_SPEED
#define CLOUD_WIND_SPEED 20.0  // m/s
#endif

#ifndef CLOUD_WIND_DIRECTION
#define CLOUD_WIND_DIRECTION float2(1.0, 0.3)  // Normalized direction
#endif

// Scattering
#ifndef CLOUD_SCATTERING_COEFFICIENT
#define CLOUD_SCATTERING_COEFFICIENT 0.5  // Reduced for more translucent clouds
#endif

#ifndef CLOUD_PHASE_G1
#define CLOUD_PHASE_G1 0.8  // Forward scattering
#endif

#ifndef CLOUD_PHASE_G2
#define CLOUD_PHASE_G2 -0.3  // Back scattering
#endif

#ifndef CLOUD_PHASE_BLEND
#define CLOUD_PHASE_BLEND 0.7  // Blend between G1 and G2
#endif

// Powder effect
#ifndef CLOUD_POWDER_STRENGTH
#define CLOUD_POWDER_STRENGTH 0.5
#endif

// =============================================================================
// MESH CLOUD SETTINGS (Vanilla Cloud Geometry Enhancement)
// =============================================================================

// Enable puffy cloud mesh geometry (normal perturbation for rounded look)
#ifndef ENABLE_PUFFY_CLOUD_GEOMETRY
#define ENABLE_PUFFY_CLOUD_GEOMETRY 1
#endif

// Noise scale for displacement pattern (smaller = larger features)
#ifndef PUFFY_CLOUD_NOISE_SCALE
#define PUFFY_CLOUD_NOISE_SCALE 0.02
#endif

// Billow frequency - how many spherical bulges per cloud
#ifndef PUFFY_CLOUD_BILLOW_FREQ
#define PUFFY_CLOUD_BILLOW_FREQ 2.0
#endif

// Maximum displacement distance (visual strength)
#ifndef PUFFY_CLOUD_DISPLACEMENT
#define PUFFY_CLOUD_DISPLACEMENT 8.0
#endif

// Vertical bias - how much displacement favors upward direction
#ifndef PUFFY_CLOUD_VERTICAL_BIAS
#define PUFFY_CLOUD_VERTICAL_BIAS 2.0
#endif

// Cloud mesh texture scale
#ifndef CLOUD_MESH_TEXTURE_SCALE
#define CLOUD_MESH_TEXTURE_SCALE 1.0
#endif

// Cloud shadow opacity (for volumetric cloud shadows on vanilla mesh)
#ifndef CLOUD_SHADOW_OPACITY
#define CLOUD_SHADOW_OPACITY 0.3
#endif

// =============================================================================
// VOLUMETRIC LIGHTING SETTINGS
// =============================================================================

// Fog parameters
#ifndef FOG_DENSITY
#define FOG_DENSITY 0.001
#endif

#ifndef FOG_HEIGHT_FALLOFF
#define FOG_HEIGHT_FALLOFF 0.1
#endif

#ifndef FOG_MAX_DISTANCE
#define FOG_MAX_DISTANCE 256.0  // Blocks (VOLUMETRIC_FOG_RANGE equivalent)
#endif

// Phase function asymmetry (Henyey-Greenstein g parameter)
// g > 0: Forward scattering (sun glow when looking toward light)
// g = 0: Isotropic (uniform in all directions)
// g < 0: Back scattering
#ifndef AIR_FOG_ASYMMETRY
#define AIR_FOG_ASYMMETRY 0.075  // Subtle forward scattering in air
#endif

#ifndef WATER_FOG_ASYMMETRY
#define WATER_FOG_ASYMMETRY 0.23  // Stronger scattering underwater
#endif

// Fog scattering coefficients (RGB wavelength-dependent)
#ifndef FOG_SCATTERING_COEFFICIENTS
#define FOG_SCATTERING_COEFFICIENTS float3(0.002, 0.00184, 0.0015)  // Rayleigh-like
#endif

// Noon fog reduction (1.0 = no reduction at midday)
#ifndef NOON_FOG_REDUCTION
#define NOON_FOG_REDUCTION 0.7  // 30% less fog at noon
#endif

// Volumetric ray marching
#ifndef VOLUMETRIC_STEPS
#define VOLUMETRIC_STEPS 32
#endif

// Static GI fog
#ifndef ENABLE_STATIC_GI_FOG
#define ENABLE_STATIC_GI_FOG 1
#endif

#ifndef STATIC_GI_FOG_AMOUNT
#define STATIC_GI_FOG_AMOUNT 0.3
#endif

// Sun fog
#ifndef ENABLE_SUN_FOG
#define ENABLE_SUN_FOG 1
#endif

#ifndef SUN_FOG_AMOUNT
#define SUN_FOG_AMOUNT 0.5
#endif

// Rain fog
#ifndef ENABLE_RAIN_FOG
#define ENABLE_RAIN_FOG 1
#endif

#ifndef RAIN_FOG_AMOUNT
#define RAIN_FOG_AMOUNT 0.8
#endif

// Rainbow
#ifndef ENABLE_RAINBOW
#define ENABLE_RAINBOW 1
#endif

#ifndef RAINBOW_INTENSITY
#define RAINBOW_INTENSITY 0.3
#endif

// =============================================================================
// WATER SETTINGS
// =============================================================================

// Physical constants
#ifndef WATER_IOR
#define WATER_IOR 1.333  // Index of refraction
#endif

#ifndef WATER_ABSORPTION_R
#define WATER_ABSORPTION_R 0.45  // Red absorption (per meter)
#endif

#ifndef WATER_ABSORPTION_G
#define WATER_ABSORPTION_G 0.029  // Green absorption
#endif

#ifndef WATER_ABSORPTION_B
#define WATER_ABSORPTION_B 0.018  // Blue absorption
#endif

// Wave parameters
#ifndef WATER_WAVE_AMPLITUDE
#define WATER_WAVE_AMPLITUDE 0.05  // Wave height multiplier
#endif

#ifndef WATER_WAVE_FREQUENCY
#define WATER_WAVE_FREQUENCY 1.0  // Wave frequency multiplier
#endif

#ifndef WATER_WAVE_SPEED
#define WATER_WAVE_SPEED 1.0  // Wave animation speed
#endif

#ifndef WATER_WAVE_OCTAVES
#define WATER_WAVE_OCTAVES 5  // Detail level
#endif

// Roughness override (-1 to use PBR)
#ifndef WATER_ROUGHNESS_OVERRIDE
#define WATER_ROUGHNESS_OVERRIDE 0.02  // Smooth water surface
#endif

// Underwater camera distortion (IOR-based refraction when camera is submerged)
#ifndef ENABLE_UNDERWATER_DISTORTION
#define ENABLE_UNDERWATER_DISTORTION 1
#endif

#ifndef UNDERWATER_DISTORTION_STRENGTH
#define UNDERWATER_DISTORTION_STRENGTH 0.15  // How much the image warps (0 = none, 1 = full IOR)
#endif

#ifndef UNDERWATER_WAVE_DISTORTION
#define UNDERWATER_WAVE_DISTORTION 0.02  // Animated wave distortion intensity
#endif

#ifndef UNDERWATER_WAVE_SPEED
#define UNDERWATER_WAVE_SPEED 2.0  // Speed of underwater wave animation
#endif

#ifndef UNDERWATER_WAVE_SCALE
#define UNDERWATER_WAVE_SCALE 3.0  // Scale of underwater wave pattern
#endif

#ifndef UNDERWATER_CHROMATIC_ABERRATION
#define UNDERWATER_CHROMATIC_ABERRATION 0.003  // RGB channel separation underwater
#endif

// Water surface refraction (looking down into water from above)
#ifndef ENABLE_WATER_REFRACTION
#define ENABLE_WATER_REFRACTION 1
#endif

#ifndef WATER_REFRACTION_STRENGTH
#define WATER_REFRACTION_STRENGTH 1.0  // Multiplier for IOR-based refraction (0 = none, 1 = physically accurate)
#endif

#ifndef WATER_REFRACTION_DEPTH_FADE
#define WATER_REFRACTION_DEPTH_FADE 16.0  // Distance at which refraction effect fades (blocks)
#endif

// Foam settings
#ifndef ENABLE_WATER_FOAM
#define ENABLE_WATER_FOAM 1
#endif

#ifndef WATER_FOAM_THRESHOLD
#define WATER_FOAM_THRESHOLD 0.8  // Edge detection threshold
#endif

#ifndef WATER_FOAM_INTENSITY
#define WATER_FOAM_INTENSITY 0.5
#endif

// Vibrant Visuals style water particles (CDOM, chlorophyll, sediment)
#ifndef ENABLE_WATER_PARTICLES
#define ENABLE_WATER_PARTICLES 1
#endif

// CDOM (Colored Dissolved Organic Matter) - yellow-brown tint
#ifndef WATER_CDOM_AMOUNT
#define WATER_CDOM_AMOUNT 0.1
#endif

// Chlorophyll - green tint from algae
#ifndef WATER_CHLOROPHYLL_AMOUNT
#define WATER_CHLOROPHYLL_AMOUNT 0.05
#endif

// Suspended sediment - red-brown turbidity
#ifndef WATER_SEDIMENT_AMOUNT
#define WATER_SEDIMENT_AMOUNT 0.02
#endif

// =============================================================================
// WATER BIOME PRESETS
// =============================================================================
// Set WATER_PRESET to use predefined water types. Set to 0 for custom values.
// 0 = Custom (uses individual settings above)
// 1 = Clear Ocean (tropical, very clear)
// 2 = Coastal Ocean (moderate visibility)
// 3 = River (greenish, moderate particles)
// 4 = Swamp (murky, high organic matter)
// 5 = Lake (clear, cold water)

#ifndef WATER_PRESET
#define WATER_PRESET 0
#endif

#if WATER_PRESET == 1  // Clear Ocean
    #undef WATER_CDOM_AMOUNT
    #undef WATER_CHLOROPHYLL_AMOUNT
    #undef WATER_SEDIMENT_AMOUNT
    #define WATER_CDOM_AMOUNT 0.02
    #define WATER_CHLOROPHYLL_AMOUNT 0.01
    #define WATER_SEDIMENT_AMOUNT 0.005
#elif WATER_PRESET == 2  // Coastal Ocean
    #undef WATER_CDOM_AMOUNT
    #undef WATER_CHLOROPHYLL_AMOUNT
    #undef WATER_SEDIMENT_AMOUNT
    #define WATER_CDOM_AMOUNT 0.08
    #define WATER_CHLOROPHYLL_AMOUNT 0.06
    #define WATER_SEDIMENT_AMOUNT 0.03
#elif WATER_PRESET == 3  // River
    #undef WATER_CDOM_AMOUNT
    #undef WATER_CHLOROPHYLL_AMOUNT
    #undef WATER_SEDIMENT_AMOUNT
    #define WATER_CDOM_AMOUNT 0.15
    #define WATER_CHLOROPHYLL_AMOUNT 0.08
    #define WATER_SEDIMENT_AMOUNT 0.06
#elif WATER_PRESET == 4  // Swamp
    #undef WATER_CDOM_AMOUNT
    #undef WATER_CHLOROPHYLL_AMOUNT
    #undef WATER_SEDIMENT_AMOUNT
    #define WATER_CDOM_AMOUNT 0.4
    #define WATER_CHLOROPHYLL_AMOUNT 0.2
    #define WATER_SEDIMENT_AMOUNT 0.1
#elif WATER_PRESET == 5  // Lake
    #undef WATER_CDOM_AMOUNT
    #undef WATER_CHLOROPHYLL_AMOUNT
    #undef WATER_SEDIMENT_AMOUNT
    #define WATER_CDOM_AMOUNT 0.05
    #define WATER_CHLOROPHYLL_AMOUNT 0.03
    #define WATER_SEDIMENT_AMOUNT 0.01
#endif

// Raindrop ripples
#ifndef ENABLE_RAIN_RIPPLES
#define ENABLE_RAIN_RIPPLES 1
#endif

#ifndef RAIN_RIPPLE_INTENSITY
#define RAIN_RIPPLE_INTENSITY 0.3
#endif

// =============================================================================
// GLASS SETTINGS
// =============================================================================

// Index of Refraction for glass (real glass is 1.5-1.52)
#ifndef GLASS_IOR
#define GLASS_IOR 1.52
#endif

// Enable glass refraction ray tracing
#ifndef ENABLE_GLASS_REFRACTION
#define ENABLE_GLASS_REFRACTION 1
#endif

// Maximum refraction bounces through glass
#ifndef GLASS_MAX_BOUNCES
#define GLASS_MAX_BOUNCES 4
#endif

// Glass absorption scale multiplier for Beer-Lambert absorption
// Higher values = stronger color tinting over distance
// 1.0 = physically-based default, higher for more dramatic tinting
#ifndef GLASS_ABSORPTION_SCALE
#define GLASS_ABSORPTION_SCALE 1.5
#endif

// Roughness for frosted glass effect (0 = clear, higher = frosted)
#ifndef GLASS_ROUGHNESS
#define GLASS_ROUGHNESS 0.0
#endif

// Enable chromatic dispersion (rainbow effects at edges)
#ifndef ENABLE_GLASS_DISPERSION
#define ENABLE_GLASS_DISPERSION 0  // Expensive, disabled by default
#endif

// Dispersion strength (IOR variation across wavelengths)
#ifndef GLASS_DISPERSION_STRENGTH
#define GLASS_DISPERSION_STRENGTH 0.02
#endif

// Thin glass threshold - glass thinner than this won't apply exit refraction
// This maintains visible distortion for thin panes. Set higher for more distortion.
#ifndef THIN_GLASS_THRESHOLD
#define THIN_GLASS_THRESHOLD 0.25
#endif

// Minimum virtual thickness for glass color tinting (blocks)
// Even paper-thin glass will show at least this much color absorption
#ifndef GLASS_MIN_THICKNESS
#define GLASS_MIN_THICKNESS 0.1
#endif

// Fresnel effect strength (0 = no Fresnel, 1 = full physical Fresnel)
// Lower values preserve more color transmission at grazing angles
#ifndef GLASS_FRESNEL_STRENGTH
#define GLASS_FRESNEL_STRENGTH 0.7
#endif

// Maximum Fresnel reflectance cap (prevents total washout at extreme angles)
// Physical max is 1.0, but capping lower preserves color visibility
#ifndef GLASS_FRESNEL_MAX
#define GLASS_FRESNEL_MAX 0.6
#endif

// =============================================================================
// CAUSTICS SETTINGS
// =============================================================================

#ifndef CAUSTICS_INTENSITY
#define CAUSTICS_INTENSITY 0.5
#endif

#ifndef CAUSTICS_FALLOFF
#define CAUSTICS_FALLOFF 0.1  // Distance falloff
#endif

#ifndef CAUSTICS_SCALE
#define CAUSTICS_SCALE 0.5  // Pattern scale
#endif

#ifndef CAUSTICS_SPEED
#define CAUSTICS_SPEED 1.0  // Animation speed
#endif

#ifndef ENABLE_REFLECTED_CAUSTICS
#define ENABLE_REFLECTED_CAUSTICS 1  // Caustics on surfaces above water
#endif

// =============================================================================
// RAYTRACED REFLECTIONS SETTINGS
// =============================================================================

#ifndef ENABLE_RAYTRACED_REFLECTIONS
#define ENABLE_RAYTRACED_REFLECTIONS 1  // Master toggle for RT reflections
#endif

#ifndef REFLECTION_MAX_ROUGHNESS
#define REFLECTION_MAX_ROUGHNESS 0.5  // Skip reflections above this roughness
#endif

#ifndef REFLECTION_MAX_BOUNCES
#define REFLECTION_MAX_BOUNCES 2  // Maximum reflection bounces
#endif

#ifndef REFLECTION_ROUGHNESS_BIAS
#define REFLECTION_ROUGHNESS_BIAS 0.0  // Bias to make reflections sharper
#endif

#ifndef ENABLE_REFLECTION_DENOISING
#define ENABLE_REFLECTION_DENOISING 1  // Temporal denoising for reflections
#endif

#ifndef REFLECTION_TEMPORAL_ALPHA
#define REFLECTION_TEMPORAL_ALPHA 0.1  // Lower = more temporal accumulation
#endif

#ifndef REFLECTION_SAMPLES
#define REFLECTION_SAMPLES 2  // Number of reflection samples per pixel (1 = sharp, 2-4 = denoised)
#endif

// =============================================================================
// SHADOW SETTINGS
// =============================================================================

#ifndef ENABLE_RAYTRACED_SHADOWS
#define ENABLE_RAYTRACED_SHADOWS 1  // Master toggle for RT shadows
#endif

#ifndef ENABLE_SOFT_SHADOWS
#define ENABLE_SOFT_SHADOWS 1  // Soft shadow penumbras
#endif

#ifndef SHADOW_SOFTNESS
#define SHADOW_SOFTNESS 0.02  // Sun angular radius for soft shadows
#endif

#ifndef SHADOW_SAMPLES
#define SHADOW_SAMPLES 1  // Shadow ray samples (1 = hard, 4+ = soft)
#endif

#ifndef ENABLE_TEMPORAL_SHADOW_FILTER
#define ENABLE_TEMPORAL_SHADOW_FILTER 1
#endif

#ifndef SHADOW_BLUR_STRENGTH
#define SHADOW_BLUR_STRENGTH 0.5
#endif

#ifndef ENABLE_COLORED_SHADOWS
#define ENABLE_COLORED_SHADOWS 1  // Through stained glass
#endif

#ifndef ENABLE_SHADOW_REFRACTION
#define ENABLE_SHADOW_REFRACTION 1  // Underwater shadow rays
#endif

#ifndef COLOR_TRANSMISSION_STRENGTH
#define COLOR_TRANSMISSION_STRENGTH 1.5  // Strength of colored light transmission
#endif

// =============================================================================
// RAIN/WETNESS SETTINGS
// =============================================================================

#ifndef RAIN_WETNESS_STRENGTH
#define RAIN_WETNESS_STRENGTH 0.8
#endif

#ifndef RAIN_WETNESS_ROUGHNESS
#define RAIN_WETNESS_ROUGHNESS 0.1  // Wet surface roughness
#endif

#ifndef RAIN_WETNESS_DARKENING
#define RAIN_WETNESS_DARKENING 0.3  // Albedo reduction
#endif

// Puddles
#ifndef ENABLE_PUDDLES
#define ENABLE_PUDDLES 1
#endif

#ifndef PUDDLE_FREQUENCY
#define PUDDLE_FREQUENCY 1.0  // Noise frequency
#endif

#ifndef PUDDLE_COVERAGE
#define PUDDLE_COVERAGE 0.3  // Maximum coverage
#endif

#ifndef PUDDLE_ROUGHNESS
#define PUDDLE_ROUGHNESS 0.05
#endif

// =============================================================================
// SNOW SETTINGS (experimental)
// =============================================================================

#ifndef SNOW_INTENSITY
#define SNOW_INTENSITY 0.5  // Based on weather
#endif

#ifndef SNOW_HEIGHT_THRESHOLD
#define SNOW_HEIGHT_THRESHOLD 100.0  // Y level for snow line
#endif

#ifndef SNOW_NORMAL_THRESHOLD
#define SNOW_NORMAL_THRESHOLD 0.7  // Upward-facing requirement
#endif

#ifndef SNOW_EDGE_ONLY
#define SNOW_EDGE_ONLY 0  // Only snow on edges
#endif

#ifndef SNOW_COLOR
#define SNOW_COLOR float3(0.95, 0.97, 1.0)
#endif

#ifndef SNOW_ROUGHNESS
#define SNOW_ROUGHNESS 0.8
#endif

// =============================================================================
// TONEMAPPING SETTINGS
// =============================================================================

// Tonemapping operator
// 0: Vanilla (filmic)
// 1: ACES
// 2: Uncharted 2 (Hable)
// 3: Reinhard
// 4: Uchimura
// 5: AgX
// 6: GTA V (cinematic filmic)
// 7: Watch Dogs (orange/teal cinematic)
// 8: None (linear)
#ifndef TONEMAPPING_TYPE
#define TONEMAPPING_TYPE 1  // ACES default
#endif

// ACES parameters
#ifndef ACES_INPUT_SCALE
#define ACES_INPUT_SCALE 0.6
#endif

#ifndef ACES_OUTPUT_SCALE
#define ACES_OUTPUT_SCALE 1.0
#endif

// GTA V tonemapper settings
#ifndef GTAV_HIGHLIGHT_DESAT
#define GTAV_HIGHLIGHT_DESAT 0.3  // Highlight desaturation (film look)
#endif

// Watch Dogs tonemapper settings
#ifndef WATCHDOGS_HIGHLIGHT_DESAT
#define WATCHDOGS_HIGHLIGHT_DESAT 0.4  // Stronger highlight desaturation
#endif

#ifndef WATCHDOGS_CONTRAST
#define WATCHDOGS_CONTRAST 1.05  // Extra contrast punch
#endif

// =============================================================================
// BLOOM SETTINGS
// =============================================================================

// Enable bloom effect
#ifndef ENABLE_BLOOM
#define ENABLE_BLOOM 1
#endif

// Brightness threshold for bloom
#ifndef BLOOM_THRESHOLD
#define BLOOM_THRESHOLD 0.8
#endif

// Soft knee for smooth threshold transition
#ifndef BLOOM_KNEE
#define BLOOM_KNEE 0.5
#endif

// Bloom intensity
#ifndef BLOOM_INTENSITY
#define BLOOM_INTENSITY 0.3
#endif

// Bloom color tint (warm for cinematic look)
#ifndef BLOOM_TINT
#define BLOOM_TINT float3(1.1, 1.0, 0.9)
#endif

// Enable anamorphic bloom (horizontal streaks)
#ifndef ENABLE_ANAMORPHIC_BLOOM
#define ENABLE_ANAMORPHIC_BLOOM 0
#endif

// Anamorphic streak intensity
#ifndef ANAMORPHIC_INTENSITY
#define ANAMORPHIC_INTENSITY 0.5
#endif

// Anamorphic spread (higher = tighter streaks)
#ifndef ANAMORPHIC_SPREAD
#define ANAMORPHIC_SPREAD 3.0
#endif

// Exposure
#ifndef AUTO_EXPOSURE_ENABLED
#define AUTO_EXPOSURE_ENABLED 1
#endif

#ifndef EXPOSURE_MIN_EV
#define EXPOSURE_MIN_EV -4.0
#endif

#ifndef EXPOSURE_MAX_EV
#define EXPOSURE_MAX_EV 12.0
#endif

#ifndef EXPOSURE_ADAPTATION_SPEED
#define EXPOSURE_ADAPTATION_SPEED 1.0
#endif

#ifndef EXPOSURE_SKY_CONTRIBUTION
#define EXPOSURE_SKY_CONTRIBUTION 1  // Include sky in metering
#endif

#ifndef EXPOSURE_EMISSIVE_CONTRIBUTION
#define EXPOSURE_EMISSIVE_CONTRIBUTION 1  // Include emissives in metering
#endif

#ifndef LOCK_EXPOSURE
#define LOCK_EXPOSURE 0  // For content creators
#endif

#ifndef LOCKED_EXPOSURE_VALUE
#define LOCKED_EXPOSURE_VALUE 0.0  // EV when locked
#endif

// =============================================================================
// POST-PROCESSING SETTINGS
// =============================================================================

// Vignette
#ifndef ENABLE_VIGNETTE
#define ENABLE_VIGNETTE 1
#endif

#ifndef VIGNETTE_INTENSITY
#define VIGNETTE_INTENSITY 0.3
#endif

#ifndef VIGNETTE_RADIUS
#define VIGNETTE_RADIUS 0.8
#endif

#ifndef VIGNETTE_SOFTNESS
#define VIGNETTE_SOFTNESS 0.5
#endif

// Lens flare
#ifndef ENABLE_LENS_FLARE
#define ENABLE_LENS_FLARE 1
#endif

#ifndef LENS_FLARE_INTENSITY
#define LENS_FLARE_INTENSITY 0.2
#endif

#ifndef LENS_FLARE_THRESHOLD
#define LENS_FLARE_THRESHOLD 0.8
#endif

#ifndef LENS_FLARE_CHROMATIC_DISTORTION
#define LENS_FLARE_CHROMATIC_DISTORTION 0.02
#endif

// Film grain
#ifndef ENABLE_FILM_GRAIN
#define ENABLE_FILM_GRAIN 0  // Disabled by default
#endif

#ifndef FILM_GRAIN_INTENSITY
#define FILM_GRAIN_INTENSITY 0.05
#endif

#ifndef FILM_GRAIN_LUMINANCE_RESPONSE
#define FILM_GRAIN_LUMINANCE_RESPONSE 0.5  // Less grain in bright areas
#endif

// Chromatic aberration
#ifndef ENABLE_CHROMATIC_ABERRATION
#define ENABLE_CHROMATIC_ABERRATION 0  // Disabled by default
#endif

#ifndef CHROMATIC_ABERRATION_STRENGTH
#define CHROMATIC_ABERRATION_STRENGTH 0.005
#endif

// Motion blur
#ifndef ENABLE_MOTION_BLUR
#define ENABLE_MOTION_BLUR 0  // Disabled by default
#endif

#ifndef MOTION_BLUR_STRENGTH
#define MOTION_BLUR_STRENGTH 0.5
#endif

#ifndef MOTION_BLUR_SAMPLES
#define MOTION_BLUR_SAMPLES 8
#endif

// Bloom
#ifndef ENABLE_BLOOM
#define ENABLE_BLOOM 1
#endif

#ifndef BLOOM_INTENSITY
#define BLOOM_INTENSITY 0.3
#endif

#ifndef BLOOM_THRESHOLD
#define BLOOM_THRESHOLD 1.0
#endif

#ifndef BLOOM_RADIUS
#define BLOOM_RADIUS 1.0
#endif

// Depth of field
#ifndef ENABLE_DOF
#define ENABLE_DOF 0  // Disabled by default, expensive
#endif

#ifndef DOF_APERTURE
#define DOF_APERTURE 0.05  // f-stop simulation
#endif

#ifndef DOF_FOCUS_DISTANCE
#define DOF_FOCUS_DISTANCE -1.0  // -1 = auto focus
#endif

#ifndef DOF_AUTO_FOCUS_SPEED
#define DOF_AUTO_FOCUS_SPEED 2.0
#endif

#ifndef DOF_NEAR_TRANSITION
#define DOF_NEAR_TRANSITION 0.5
#endif

#ifndef DOF_FAR_TRANSITION
#define DOF_FAR_TRANSITION 2.0
#endif

// =============================================================================
// DIMENSION-SPECIFIC SETTINGS
// =============================================================================

// Nether
#ifndef NETHER_EXPOSURE_MIN
#define NETHER_EXPOSURE_MIN -2.0
#endif

#ifndef NETHER_EXPOSURE_MAX
#define NETHER_EXPOSURE_MAX 6.0
#endif

#ifndef NETHER_EXPLICIT_LIGHT_INTENSITY
#define NETHER_EXPLICIT_LIGHT_INTENSITY 2.0
#endif

#ifndef NETHER_AMBIENT_LIGHT
#define NETHER_AMBIENT_LIGHT float3(0.15, 0.05, 0.02)
#endif

#ifndef NETHER_FOG_AMOUNT
#define NETHER_FOG_AMOUNT 0.5
#endif

// The End
#ifndef END_SKY_ENABLED
#define END_SKY_ENABLED 1
#endif

#ifndef END_SUN_DIRECTION
#define END_SUN_DIRECTION float3(0.0, 0.8, 0.6)
#endif

#ifndef END_SUN_INTENSITY
#define END_SUN_INTENSITY 20.0
#endif

#ifndef END_SUN_RADIUS
#define END_SUN_RADIUS 0.02  // Radians
#endif

#ifndef END_MIE_SCATTERING
#define END_MIE_SCATTERING 1
#endif

#ifndef END_LIGHTNING_FLASHES
#define END_LIGHTNING_FLASHES 0  // Experimental
#endif

#ifndef END_FOG_DENSITY
#define END_FOG_DENSITY 0.002
#endif

#ifndef END_FOG_COLOR
#define END_FOG_COLOR float3(0.05, 0.02, 0.1)
#endif

// =============================================================================
// STATUS EFFECT SETTINGS
// =============================================================================

#ifndef NIGHT_VISION_EXPOSURE_BOOST
#define NIGHT_VISION_EXPOSURE_BOOST 4.0  // EV boost
#endif

#ifndef NIGHT_VISION_AMBIENT_BOOST
#define NIGHT_VISION_AMBIENT_BOOST 0.3
#endif

#ifndef DARKNESS_PULSE_ENABLED
#define DARKNESS_PULSE_ENABLED 1
#endif

#ifndef BLINDNESS_SKY_FIX
#define BLINDNESS_SKY_FIX 1
#endif

// =============================================================================
// FLASHLIGHT SETTINGS
// =============================================================================

#ifndef FLASHLIGHT_INTENSITY
#define FLASHLIGHT_INTENSITY 50.0
#endif

#ifndef FLASHLIGHT_INNER_CONE
#define FLASHLIGHT_INNER_CONE 0.3  // Radians
#endif

#ifndef FLASHLIGHT_OUTER_CONE
#define FLASHLIGHT_OUTER_CONE 0.5  // Radians
#endif

#ifndef FLASHLIGHT_COLOR
#define FLASHLIGHT_COLOR float3(1.0, 0.95, 0.9)
#endif

#ifndef FLASHLIGHT_RANGE
#define FLASHLIGHT_RANGE 50.0  // Blocks
#endif

#ifndef FLASHLIGHT_VOLUMETRIC
#define FLASHLIGHT_VOLUMETRIC 1  // Beam through fog
#endif

#ifndef FLASHLIGHT_BOUNCE
#define FLASHLIGHT_BOUNCE 1  // Sway when walking
#endif

// =============================================================================
// PERFORMANCE SETTINGS
// =============================================================================

// Distance-based quality reduction
#ifndef ENABLE_DISTANCE_LOD
#define ENABLE_DISTANCE_LOD 1
#endif

#ifndef LOD_DISTANCE_THRESHOLD
#define LOD_DISTANCE_THRESHOLD 64.0  // Blocks
#endif

#ifndef LOD_SKIP_SPECULAR_ROUGH
#define LOD_SKIP_SPECULAR_ROUGH 1  // Skip specular for rough distant surfaces
#endif

#ifndef LOD_ROUGH_THRESHOLD
#define LOD_ROUGH_THRESHOLD 0.7  // Roughness threshold for skipping
#endif

#ifndef LOD_GI_EARLY_TERMINATION
#define LOD_GI_EARLY_TERMINATION 1  // Early exit for distant GI
#endif

// Pixelated mode (retro)
#ifndef ENABLE_PIXELATED_MODE
#define ENABLE_PIXELATED_MODE 0
#endif

#ifndef PIXELATED_SHADOWS
#define PIXELATED_SHADOWS 0
#endif

#ifndef PIXELATED_REFLECTIONS
#define PIXELATED_REFLECTIONS 0
#endif

#ifndef PIXELATED_GI
#define PIXELATED_GI 0
#endif

// =============================================================================
// DEBUG SETTINGS
// =============================================================================

#ifndef DEBUG_VIEW
// 0: None
// 1: Albedo
// 2: Normal
// 3: Roughness
// 4: Metalness
// 5: Emissive
// 6: AO
// 7: Depth
// 8: Motion vectors
// 9: Overdraw
// 10: Performance heatmap
#define DEBUG_VIEW 0
#endif

#ifndef DEBUG_NAN_INF
#define DEBUG_NAN_INF 1  // Highlight NaN/Inf with checker pattern
#endif

// =============================================================================
// PHYSICS CONSTANTS (do not modify unless you know what you're doing)
// =============================================================================

static const float kPi = 3.14159265358979323846;
static const float kTwoPi = 6.28318530717958647692;
static const float kHalfPi = 1.57079632679489661923;
static const float kInvPi = 0.31830988618379067154;
static const float kInvTwoPi = 0.15915494309189533577;
static const float kSqrtTwo = 1.41421356237309504880;
static const float kInvSqrtTwo = 0.70710678118654752440;

// Speed of light in vacuum (m/s)
static const float cSpeedOfLight = 299792458.0;

// Planck's constant (J*s)
static const float cPlanckConstant = 6.62607015e-34;

// Boltzmann constant (J/K)
static const float cBoltzmannConstant = 1.380649e-23;

// Stefan-Boltzmann constant (W/(m^2*K^4))
static const float cStefanBoltzmann = 5.670374419e-8;

// Standard gravity (m/s^2)
static const float cGravity = 9.80665;

// Water refractive index
static const float cWaterRefractiveIndex = 1.333;

// Glass refractive index (uses GLASS_IOR setting)
static const float cGlassRefractiveIndex = GLASS_IOR;

// Air refractive index at sea level
static const float cAirRefractiveIndex = 1.000293;

// Sun angular radius (radians)
static const float kSunRadiusRadians = 0.00465; // ~0.533 degrees

// Moon angular radius (radians)
static const float kMoonRadiusRadians = 0.00452; // ~0.517 degrees

#endif // __OPENRTX_SETTINGS_HLSL__
