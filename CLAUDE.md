# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OpenRTX is an open-source Minecraft Bedrock RTX shader built on veka0's mcrtx-shader-template. It uses a templating engine to extract data from vanilla material.bin files and generate compilable HLSL shader code.

## Build Commands

### Generate project from template (required before first build)
```bash
python src/script.py
```

### Compile shaders
```bash
python -m lazurite build ./project -o ./output
```

## Repository Structure

- `template/` - Source shader code with templating tokens (e.g., `{{RTXStub.passes.TAA.group_size}}`)
- `project/` - Auto-generated compilable project (do not edit directly - changes will be overwritten)
- `src/` - Python scripts for template processing and material analysis
- `vanilla/` - Place vanilla material.bin files here (RTXStub, RTXPostFX.Bloom, RTXPostFX.Tonemapping)

**Key distinction**: Edit files in `template/`, not `project/`. Run `python src/script.py` to regenerate `project/` after template changes.

## Key Shader Files

| File | Purpose |
|------|---------|
| `Settings.hlsl` | All configurable shader parameters (~200+ settings) |
| `OpenRTX.hlsl` | Main integration file |
| `BRDF.hlsl` | Physically-based lighting models |
| `Clouds.hlsl` | Volumetric cloud rendering |
| `Sky.hlsl` | Atmospheric scattering |
| `Water.hlsl` | Water effects with Gerstner waves |
| `VolumetricLighting.hlsl` | Fog and god rays |
| `ToneMapping.hlsl` | Tonemapping operators and color grading |
| `Generated/Signature.hlsl` | Auto-generated resource signatures |
| `Generated/Structs.hlsl` | Auto-generated shader structs |

## HLSL Coding Conventions

### Naming
- Variables: `camelCase`
- Constants/Defines: `UPPER_SNAKE_CASE`
- Functions: `camelCase`
- Static constants: `kPascalCase` (e.g., `kEarthRadius`)

### Settings Pattern
Always wrap new settings with `#ifndef` guards for external overrides:
```hlsl
#ifndef MY_NEW_SETTING
#define MY_NEW_SETTING (default_value)
#endif
```

### Performance
- Use `saturate()` over `clamp(x, 0, 1)`
- Use `isLightSample` parameter to skip expensive calculations for shadow rays
- Implement early-out conditions and adaptive step sizing for ray marching

## Material Format

- **MERS texture format**: Metallic / Emissive / Roughness / Subsurface (RGBA channels)
- **Normal maps**: DirectX format with optional AO in blue channel and heightmap for POM in alpha

## Dependencies

- Python 3.10+ (3.12 recommended)
- [Lazurite](https://github.com/veka0/lazurite) - Python library for Bedrock shaders
- [DXC](https://github.com/microsoft/DirectXShaderCompiler/releases) - DirectX Shader Compiler (for RTXStub)
- [Shaderc](https://github.com/veka0/bgfx-mcbe/releases/tag/binaries) - BGFX shader compiler (for PostFX)

Place compiler executables in the repository root directory.

## Pipeline Overview

The RTXStub pipeline runs ~30+ compute passes in sequence. Key passes for development:
- `PreBlasSkinning` - Entity animation and motion vectors
- `PrimaryCheckerboardRayGenInline` - Main ray tracing logic
- `FinalCombine` - Composites all lighting contributions
- `CopyToFinal` - Transfers upscaled image to output

PostFX renders after RTXStub: BloomDownscale → BloomUpscale → Tonemap

## Best Practices

See [HLSL Best Practices](docs/datasets/HLSL-BEST-PRACTICES.md) for detailed optimization guidelines.

### Template System Usage
- **Hierarchical Token Naming**: Use `{{MaterialName.category.specific_item}}` pattern (e.g., `{{RTXStub.passes.TAA.group_size}}`)
- **Centralized Resources**: All resource bindings declared in `Generated/Signature.hlsl` using tokens
- **Dynamic Configuration**: Compute shader `[numthreads]` attributes can be set via tokens
- **Token Validation**: Ensure all tokens exist in `template_token_mapping` or they'll be skipped during generation

### Pipeline Architecture
The complete RTX pipeline executes ~55-60 compute shader dispatches per frame:

1. **Geometry Preprocessing** (6 passes) - Skin animation, face data, irradiance cache
2. **BVH Construction** (API call)
3. **Light Metering** (2 passes) - Histogram and tone curve calculation
4. **Ray Tracing** (8 passes) - Primary, shadows, diffuse, specular, refraction
5. **Volumetric Effects** (6 passes) - Inscatter, GI, blur, accumulation
6. **Denoising** (21+ passes) - Atrous, temporal, shadow, firefly filtering
7. **Final Composition** (5-6 passes) - Combine, TAA/DLSS, copy, histogram
8. **Post-Processing** (9+ passes) - Bloom pyramid, tonemapping

### Critical Pipeline Integration Points
- **Upscaling Support**: Must output three buffers from primary pass:
  - `outputBufferFinal` (HDR color)
  - `outputBufferMotionVectors` (screen-space velocity)
  - `outputBufferReprojectedPathLength` (ray hit distance)
- **Root Constants**: Use `g_rootConstant0` etc. for per-dispatch configuration without constant buffer updates
- **Buffer Dependencies**: Track data flow between passes - many passes depend on previous pass outputs

### Common Pitfalls
- **Missing Tokens**: Unmatched `{{token}}` placeholders cause compilation failures
- **Binary Files**: Tokens only work in text files, not binary resources
- **Upscaling Padding**: `CopyToFinal` dispatch needs ~1/7 extra resolution when upscaling enabled
- **TAA/DLSS Branch**: Pipeline behavior changes significantly based on `isUpscalingEnabled()`

Utilize skills and MCP tools!