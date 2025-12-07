# OpenRTX - Open-Source BetterRTX Parity Shader

OpenRTX is an open-source Minecraft Bedrock RTX shader that aims for feature parity with BetterRTX. Built on top of veka0's mcrtx-shader-template, it provides a comprehensive, community-maintained alternative with physically-based rendering.

## Features

### Core Rendering Enhancements

#### Settings System (`Settings.hlsl`)
- Comprehensive configuration with 200+ settings
- All settings use `#ifndef` guards for external overrides
- Organized by category with vanilla comparison values
- Feature toggles for all major systems

#### BRDF System (`BRDF.hlsl`)
Multiple diffuse BRDF options:
- **Lambertian** (fastest)
- **Burley/Disney** (good balance)
- **Frostbite Burley** (energy conserving, default)
- **Oren-Nayar** (good for rough surfaces)

Specular features:
- GGX/Trowbridge-Reitz distribution
- Beckmann and Blinn-Phong options
- Smith GGX height-correlated geometry
- Schlick Fresnel with roughness adjustment
- Multi-scatter energy compensation

#### Enhanced Material System
- Subsurface scattering (SSS) support
- Multi-bounce AO approximation
- Reflection occlusion
- Default material roughness fix (100% rough instead of 0%)
- Entity emissive support

### Atmospheric Rendering

#### Sky System (`Sky.hlsl`)
Physically-based atmospheric scattering:
- Rayleigh scattering for blue sky
- Mie scattering for sun glow and haze
- Henyey-Greenstein and Cornette-Shanks phase functions
- Closed-form approximations for performance
- Day/night transitions with color temperature
- Procedural star field
- Blackbody color temperature for sun/moon

#### Volumetric Clouds (`Clouds.hlsl`)
Real-time ray-marched volumetric clouds:
- **Cumulus layer** (1500m-4000m):
  - Multi-octave FBM noise (5 octaves)
  - Worley noise erosion for detail
  - Dual-lobe Henyey-Greenstein phase function
  - Beer-Lambert attenuation
  - Powder effect for dark edges
  - Adaptive step sizing
- **Cirrus layer** (~8000m):
  - 2D texture with parallax
  - Forward scattering effects
- Cloud shadows on terrain
- Wind animation with configurable speed/direction

### Water Rendering (`Water.hlsl`)

Advanced water effects:
- **Gerstner wave** displacement with parallax mapping
- **FBM noise** for fine detail
- **Vibrant Visuals parity** water particles:
  - CDOM (yellow-brown tint, absorbs blue)
  - Chlorophyll (green tint from algae)
  - Suspended sediment (red-brown turbidity)
- **Seafoam** at shorelines
- **Raindrop ripples** during weather
- Configurable roughness override
- Physically-based Fresnel and refraction

### Volumetric Lighting (`VolumetricLighting.hlsl`)

- Height-based exponential fog
- Static GI fog contribution
- Sun/moon fog with phase functions
- Rain fog with weather scaling
- Rainbow effects during weather transitions
- Dimension-specific fog (Nether, End)
- Flashlight volumetric beam (optional)

### Post-Processing (`ToneMapping.hlsl`)

#### Tonemapping Operators
1. **Vanilla** (filmic) - Original MC RTX style
2. **ACES** - Academy Color Encoding System
3. **Uncharted 2** (Hable) - Cinematic look
4. **Reinhard** - Classic operator
5. **Uchimura** - Gran Turismo style
6. **AgX** - Modern filmic look
7. **None** - Linear output

#### Color Grading
- Lift/Gamma/Gain adjustment
- Contrast and saturation
- Color temperature shift
- Vignette effect

#### Post-Processing Effects
- **Vignette** with configurable intensity/radius/softness
- **Lens flare** with ghost artifacts and chromatic distortion
- **Film grain** (animated, luminance-responsive)
- **Chromatic aberration** (edge distortion)
- **Motion blur** (velocity-based)
- **Bloom** (multi-pass with threshold)
- **Depth of Field** (ray-traced with aperture control)
- **Dithering** (Bayer matrix for banding reduction)

### Bug Fixes

All vanilla bug fixes included:
- `FIX_DEFAULT_MATERIAL` - 100% rough default instead of 0%
- `FIX_BACKFACE_CULLING` - Skip back faces of opaque blocks
- `FIX_BLEND_SURFACES` - Skip fully transparent blend objects
- `FIX_BANNER_UVS` - Correct 1px banner UV offset
- `FIX_ITEM_ALBEDO_DARKENING` - Bow/crossbow fix
- `FIX_TINT_SHADED_SURFACES` - Remove baked vanilla shading
- `FIX_BLOCK_BREAKING_OVERLAY` - Multiplicative blend
- `FIX_ENTITY_EMISSIVE` - Entity emissive detection

### Dimension Support

#### Nether
- Custom exposure limits
- Enhanced explicit light intensity
- Ambient light for visibility
- Orange-red fog

#### The End
- Custom procedural sky
- End sun with Mie scattering
- Purple fog atmosphere

### Performance Features

- Distance-based quality reduction (LOD)
- Skip specular for rough distant surfaces
- GI ray early termination
- Closed-form atmospheric approximations
- Adaptive cloud step sizing

## File Structure

```
RTXStub/shaders/Include/
├── Settings.hlsl           # All configuration parameters
├── BRDF.hlsl               # Lighting models
├── Sky.hlsl                # Atmospheric scattering
├── Clouds.hlsl             # Volumetric clouds
├── Water.hlsl              # Water rendering
├── VolumetricLighting.hlsl # Fog and god rays
├── ToneMapping.hlsl        # Tonemapping operators
├── GFXHelpers.hlsl         # Graphics utilities
├── OpenRTX.hlsl            # Main integration file
└── Generated/              # Auto-generated signatures
```

## Configuration

### Quick Start

To enable/disable major features, edit `Settings.hlsl`:

```hlsl
#define OPENRTX_ENABLED 1           // Master toggle
#define ENABLE_PBR_LIGHTING 1       // Physically-based BRDF
#define ENABLE_ATMOSPHERIC_SKY 1    // Rayleigh/Mie sky
#define ENABLE_VOLUMETRIC_CLOUDS 1  // Ray-marched clouds
#define ENABLE_VOLUMETRIC_LIGHTING 1 // Fog effects
#define ENABLE_WATER_EFFECTS 1      // Wave animation
#define ENABLE_POST_PROCESSING 1    // All post-FX
```

### Presets

**Performance Preset**:
```hlsl
#define CLOUD_MARCH_STEPS 32
#define VOLUMETRIC_STEPS 16
#define WATER_WAVE_OCTAVES 3
#define ENABLE_VOLUMETRIC_CLOUDS 0
```

**Quality Preset** (default):
```hlsl
#define CLOUD_MARCH_STEPS 64
#define VOLUMETRIC_STEPS 32
#define WATER_WAVE_OCTAVES 5
```

**Cinematic Preset**:
```hlsl
#define CLOUD_MARCH_STEPS 128
#define VOLUMETRIC_STEPS 64
#define ENABLE_DOF 1
#define ENABLE_MOTION_BLUR 1
#define TONEMAPPING_TYPE 1  // ACES
```

### Tonemapping Selection

```hlsl
#define TONEMAPPING_TYPE 0  // Vanilla
#define TONEMAPPING_TYPE 1  // ACES (default)
#define TONEMAPPING_TYPE 2  // Uncharted 2
#define TONEMAPPING_TYPE 3  // Reinhard
#define TONEMAPPING_TYPE 4  // Uchimura
#define TONEMAPPING_TYPE 5  // AgX
```

### BRDF Selection

```hlsl
#define DIFFUSE_BRDF 0      // Lambertian
#define DIFFUSE_BRDF 1      // Burley
#define DIFFUSE_BRDF 2      // Frostbite Burley (default)
#define DIFFUSE_BRDF 3      // Oren-Nayar
```

## Build Instructions

1. Install requirements:
   - Python 3.10+
   - [Lazurite](https://github.com/veka0/lazurite)
   - [DXC compiler](https://github.com/microsoft/DirectXShaderCompiler/releases)
   - [Shaderc compiler](https://github.com/veka0/bgfx-mcbe/releases/tag/binaries)

2. Place vanilla material.bin files in `vanilla/`:
   - `RTXStub.material.bin`
   - `RTXPostFX.Bloom.material.bin`
   - `RTXPostFX.Tonemapping.material.bin`

3. Build:
```bash
lazurite build project/ -o ./
```

## Docker Build

```dockerfile
FROM jasongardner/dxc:latest AS dxc
# Add your build steps here
```

## Compatibility

- Minecraft Bedrock 1.21+
- NVIDIA RTX (all generations)
- AMD RDNA2+ (RX 6000+)
- Intel ARC (limited testing)
- DLSS/FSR upscaling support

## License

MIT License - See LICENSE file for details.

## Credits

- **veka0** - mcrtx-shader-template foundation
- **OpenRTX Contributors** - Shader implementation
- **BetterRTX Team** - Feature inspiration
- **Guerrilla Games** - Horizon cloud rendering technique
- **Epic/Unreal** - BRDF implementations
- **Academy** - ACES tonemapping

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description

## Known Issues

- Cirrus clouds are simplified (2D plane vs. 3D)
- Caustics are procedural approximation
- Some Intel ARC features may need adjustment
- Motion blur requires velocity buffer (placeholder)

## Roadmap

- [ ] Full GI integration with denoiser
- [ ] Screen-space reflections fallback
- [ ] Improved caustics with ray tracing
- [ ] Weather system integration
- [ ] Performance profiling tools
