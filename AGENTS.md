# OpenRTX

An open-source Minecraft Bedrock RTX shader.

## Development

- Reference `docs/HLSL-BEST-PRACTICES.md` for best practices.
- Always utilize @web search and MCP tools.

## Building

`python -m lazurite build ./project -o ./output`

## Materials

MERS: Metallic / Emissive / Roughness / Subsurface
Assigned to RGBA, respective channels.

Normal map is DirectX. It includes an optional AO channel in the blue channel and an optional heightmap for POM in the alpha channel.

---

## **Authoritative Guidelines for AI Agents Writing HLSL Shader Code**

This document defines the **non-negotiable standards**, **best practices**, and **mandatory conventions** for **AI agents** (or agentic systems) involved in the generation of HLSL shader code. The objective is to ensure correctness, performance, readability, and compatibility with real-time graphics pipelines (i.e. DirectX 12 / RenderDragon)

## ✅ Absolute Rules (MUST-FOLLOW)

### 1. **Shader Compilation Must Succeed**

* Never emit incomplete or un-compilable shader code.
* Always include valid `float`, `float2/3/4`, `sampler`, `cbuffer`, `struct`, and return types.
* Entry functions must match shader type requirements (`VS`, `PS`, `CS`, etc.).

  * E.g., vertex shaders must return a struct with at least a `SV_POSITION` semantic.

### 2. **Respect HLSL Syntax & Grammar**

* Terminate every statement with `;`.
* Avoid implicit casting (e.g., don’t assume `int -> float` will work safely).
* Use correct types (`float4`, not `vec4`), correct semantics (`SV_Target`, not `gl_FragColor`).
* Use lowercase keywords (`float`, `return`, `if`, etc.), PascalCase for types.

### 3. **Define All Inputs/Outputs Explicitly**

* Always define vertex input/output structs clearly:

  * Use proper **semantics** like `POSITION`, `TEXCOORD`, `NORMAL`, `SV_Position`, etc.
* Fragment/Pixel Shader outputs must be bound to `SV_Target`.

### 4. **Maintain Cross-Hardware Compatibility**

* Avoid use of vendor-specific extensions or types (e.g., NVIDIA-only intrinsics).
* Avoid undefined behavior (e.g., out-of-bounds texture reads, `NaN` propagation).
* Use standard precision (`float`) unless `min16float` is intentional for mobile.

### 5. **Comment Critical Sections Clearly**

* Include comments before complex math, especially lighting, BRDFs, UV distortions, etc.
* Use `//` for inline or block context, never leave “magic numbers” unexplained.

---

## 🔧 Best Practices (HIGHLY RECOMMENDED)

### 1. **Modular Code Organization**

* Separate lighting, material, and utility functions.
* Favor `#include` or well-defined function blocks to isolate responsibility.

### 2. **Use Constant Buffers for Uniforms**

* Group related uniforms in `cbuffer` blocks.

```hlsl
cbuffer PerFrame : register(b0)
{
    float4x4 ViewProjection;
    float3 CameraPosition;
    float Time;
};
```

### 3. **Use Linear Space for Lighting**

* Convert inputs from gamma if necessary (`pow(color, 2.2)`).
* Do lighting calculations in linear space and convert final result to gamma if required.

### 4. **Profile and Optimize**

* Minimize expensive operations (`normalize`, `pow`, `exp`, `dot` if repeated).
* Reduce dependent texture reads in pixel shaders.
* Leverage shader model capabilities (e.g., SM 5.0, SM 6.0 features if target allows).

### 5. **Use Precise Interpolation**

* For interpolators, specify interpolation modifiers where needed:

  * `nointerpolation`, `linear`, `centroid`, `noperspective`

---

## ⚠️ Common Mistakes to Avoid

| Mistake                                 | Consequence                                       |
| --------------------------------------- | ------------------------------------------------- |
| Writing `vec3`, `mix`, `fract`          | GLSL syntax — not valid in HLSL                   |
| Forgetting `SV_Target` or `SV_Position` | Shader won’t compile or render                    |
| Using `float3 pos : POSITION;` in PS    | POSITION semantic is vertex-only                  |
| Ignoring shader stage roles             | E.g., lighting in vertex shader = visual artifact |
| Hardcoded magic numbers                 | Breaks modularity and tuning flexibility          |

---

## 🧠 Critical Shader Knowledge AI Agents MUST Understand

1. **Shader Pipeline Stages**

   * Vertex Shader: transforms mesh vertex data.
   * Pixel Shader: computes color per pixel.
   * Optional: Geometry, Hull, Domain, Compute.

2. **Semantics and Binding**

   * Semantics bind data between stages. E.g., `POSITION`, `TEXCOORD0`, `NORMAL`, `SV_Target`, etc.
   * Registers and resource binding (`t0`, `s0`, `b0`) map shader variables to engine-side data.

3. **Math Functions**

   * Prefer HLSL-native intrinsics (`lerp`, `saturate`, `dot`, `normalize`).
   * Understand coordinate space transforms (world, view, projection).
   * Know lighting models: Lambert, Blinn-Phong, Cook-Torrance (PBR).

4. **Precision & GPU Cost**

   * Each `normalize` or `pow` carries ALU cost.
   * Dependent texture reads stall texture units.
   * Use `min16float` and loop unrolling wisely in performance-critical shaders.

5. **Debugging Shader Issues**

   * Visual artifacts may arise from:

     * Incorrect interpolators
     * Uninitialized variables
     * NaN propagation (e.g., normalize(0))
     * Mismatched input/output structures

---

## 🧩 Sample Template (Pixel Shader)

```hlsl
Texture2D DiffuseMap : register(t0);
SamplerState Sampler : register(s0);

struct VS_OUTPUT {
    float4 Position : SV_POSITION;
    float2 UV       : TEXCOORD0;
};

float4 main(VS_OUTPUT input) : SV_Target
{
    float4 baseColor = DiffuseMap.Sample(Sampler, input.UV);
    return baseColor;
}
```

---

## 📌 Final Reminders for Agents

* **Be deterministic**: Output the same shader given the same intent.
* **Avoid randomness** unless procedural noise is specifically requested.
* **Use concise, human-readable naming**.
* Always verify code logic **before emitting**.

---

## 🧪 Validation Checklist

Before emitting HLSL shader code, ensure:

* [ ] Shader compiles with `dxc`.
* [ ] Entry point (`main`) is defined and correctly bound.
* [ ] All semantics are valid for the target stage.
* [ ] All types are HLSL-specific (no GLSL/WebGL syntax).
* [ ] No unresolved functions or parameters exist.
* [ ] No undefined behavior (e.g., div by 0, normalize(0)).

---

## 🔄 Continuous Evolution

This `AGENTS.md` is a living document. Feedback and contributions from developers and agent performance reviews will continuously refine its contents to ensure shaders produced are both **state-of-the-art** and **production-safe**.
