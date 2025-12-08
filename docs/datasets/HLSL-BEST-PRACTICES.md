# **Definitive Guide to High-Level Shading Language (HLSL) Optimization for RTX Architectures**

## **1\. Architectural Foundations of Hardware-Accelerated Ray Tracing**

The advent of real-time ray tracing, catalyzed by NVIDIA's RTX architecture and standardized through DirectX Raytracing (DXR), represents the most significant paradigm shift in graphics rendering since the introduction of programmable shaders. Unlike rasterization, which relies on the inherent coherence of projection to map geometry to pixels, ray tracing is fundamentally a scatter-gather operation characterized by incoherent memory access and divergent execution paths. To author high-performance HLSL shaders for this architecture, one must first understand the underlying hardware mechanisms that govern traversal, intersection, and shading.

### **1.1 The Heterogeneous Compute Model**

The RTX architecture introduces a heterogeneous compute model where execution toggles between general-purpose Streaming Multiprocessors (SMs) and fixed-function RT Cores. The SMs are responsible for executing the HLSL shader code (Ray Generation, Closest Hit, Any Hit, Miss), while the RT Cores are specialized circuits dedicated solely to traversing the Bounding Volume Hierarchy (BVH) and performing ray-triangle intersections.

Performance in this environment is defined by the efficiency of the handover between these two units. When a ray is traced via the TraceRay() intrinsic, the SM suspends the execution of the calling thread and offloads the ray to the RT Core. The RT Core performs the traversal. If an intersection is found, control may return to the SM to execute an Any-Hit shader (for alpha testing) or, upon final intersection, a Closest Hit shader.

A critical insight for optimization is that the SMs are throughput-oriented processors designed to hide latency by switching between active warps (groups of 32 threads). However, the context switch required to invoke a shader during traversal—specifically the Any-Hit shader—interrupts the fixed-function traversal pipeline. This interruption forces the system to spill the current traversal state to memory and load the shader execution context, creating a significant performance penalty often referred to as the "Any-Hit Tax." Consequently, a primary goal of HLSL optimization is to minimize the frequency of these interruptions and keep the RT Cores saturated with opaque geometry traversal.1

### **1.2 Divergence: The Primary Bottleneck**

In rasterization, neighboring pixels typically map to the same triangle and execute the same pixel shader, ensuring high lane utilization within a warp. In ray tracing, especially with secondary rays (reflection, refraction, global illumination), rays launched by adjacent threads in a warp may scatter in different directions. One ray might hit a skybox (Miss Shader), another a diffuse wall (Lambertian Shader), and a third a glass object (Refraction Shader).

This phenomenon, known as execution divergence, is catastrophic for Single Instruction, Multiple Thread (SIMT) architectures. The hardware must serialize the execution of these different shader paths. If a warp contains 32 threads that hit 32 different materials, the warp effectively executes sequentially, reducing throughput to 1/32nd of peak performance. Furthermore, the scattered memory accesses resulting from these divergent rays lead to "data divergence," causing Translation Lookaside Buffer (TLB) thrashing and poor L1/L2 cache hit rates.3

Effective HLSL authoring for RTX is therefore less about optimizing arithmetic logic unit (ALU) instructions—though that remains important—and more about managing control flow, memory access patterns, and minimizing the "live state" that must be persisted across these divergent execution paths.

## ---

**2\. Acceleration Structure Management and Traversal Optimization**

Before a shader can execute, the ray must traverse the scene. The efficiency of this traversal is dictated by the construction and management of the Acceleration Structure (AS). The AS serves as the spatial index of the scene, and its quality determines the number of memory fetches and intersection tests the RT Core must perform.

### **2.1 Geometry Flags and the "Opaque Cliff"**

The most immediate and impactful optimization available to developers is the correct usage of geometry flags. The DXR specification allows geometry to be flagged as D3D12\_RAYTRACING\_GEOMETRY\_FLAG\_OPAQUE. When this flag is set, the RT Core is informed that the geometry contains no alpha testing logic. Upon intersection, the hardware can immediately commit the hit (updating the closest t-value) and continue traversal without interrupting the SM.

Conversely, if geometry is not flagged as opaque, the RT Core must assume that an Any-Hit shader might accept or reject the intersection based on programmable logic (e.g., texture alpha sampling). Even if no Any-Hit shader is bound, the validation overhead exists. If an Any-Hit shader *is* bound, the traversal pauses, the state is spilled, and the shader executes.

Benchmarks and architectural analysis reveal a performance "cliff" associated with non-opaque geometry. The cost of invoking Any-Hit shaders is high enough that it is strongly recommended to segregate geometry into distinct BLASs: one strictly for opaque meshes (force-flagged as opaque) and another for alpha-tested meshes. It is often more performant to render alpha-tested geometry (like vegetation) using a specialized simplified shader or even Opacity Micromaps rather than relying on the generic Any-Hit path.1

### **2.2 Topology Optimization and Build Preference**

The topology of the underlying mesh significantly influences the quality of the resulting BVH. RTX GPUs are optimized for triangle intersection. While AABBs can be used for volumetric or procedural primitives, they rely on software intersection shaders, which bypass the RT Core's intersection unit and incur significant SM overhead.

#### **2.2.1 The Problem of Elongated Triangles**

The Surface Area Heuristic (SAH) used to build BVHs relies on the probability of a ray hitting a volume. Long, thin triangles ("slivers") produce bounding boxes with large surface areas relative to the geometry they contain. These boxes often overlap significantly with neighboring nodes. When a ray enters a region of space with overlapping nodes, the traversal unit must test multiple branches of the tree, increasing the number of steps required to find a hit or determine a miss.

Developers should preprocess geometry to avoid slivers, potentially by tessellating long triangles into more equilateral shapes. This increases the vertex count but significantly improves traversal coherency and performance.2

#### **2.2.2 Build Flags: Quality vs. Speed**

DXR provides build flags that control the trade-off between build time and trace performance:

| Build Flag | Description | Optimization Context |
| :---- | :---- | :---- |
| **PREFER\_FAST\_TRACE** | Optimization target is high-quality BVH. | **Mandatory for Static Geometry.** The build takes longer, but the resulting tree is tighter, minimizing traversal steps. Use for TLAS and static BLAS. 1 |
| **PREFER\_FAST\_BUILD** | Optimization target is fast construction. | Use for highly dynamic geometry that is rebuilt every frame (e.g., particles, exploding meshes). The trace will be slower due to lower quality BVH. |
| **None (Balanced)** | Compromise between build and trace. | Suitable for deforming meshes where topology is constant but vertex positions change. |

For the Top-Level Acceleration Structure (TLAS), the recommendation is almost universally to use PREFER\_FAST\_TRACE and perform a full rebuild every frame. The TLAS usually contains thousands of instances (compared to millions of triangles in BLAS), so the build cost is negligible compared to the performance gain from a high-quality instance hierarchy.2

### **2.3 Memory Management: Compaction and Pooling**

Acceleration structures reside in GPU memory (VRAM). The initial allocation for a BLAS build is conservative, as the driver cannot know the final compacted size of the BVH. DXR supports **compaction**, a process where the driver copies the built AS to a smaller buffer, discarding unused memory.

#### **2.3.1 Compaction Strategy**

Compaction is essential for memory-constrained scenarios and can improve cache locality. However, it introduces latency (a build followed by a copy). The best practice is to always compact static BLASs that are built once at level load. For dynamic BLASs, the cost of compaction per frame usually outweighs the memory savings, unless VRAM is critically low. Compaction can reduce memory footprint by up to 50%, which frees space for higher resolution textures or more geometry, indirectly benefiting overall visual quality.2

#### **2.3.2 Pooling and Alignment**

Creating a separate D3D12 committed resource for every small BLAS leads to memory fragmentation and overhead. A more efficient approach is to manage a large heap (or placed resources) and suballocate BLAS memory within it. While standard resources often require 64KB alignment, suballocating within a buffer allows for much tighter packing (256-byte alignment). This reduces the "alignment waste" that can otherwise consume hundreds of megabytes in scenes with many small objects.4

## ---

**3\. HLSL Shader Authoring: Payload and State Management**

Once the ray traverses the scene and hits a surface, the HLSL shader code executes. The efficiency of this execution is governed by register usage, payload size, and the management of "live state."

### **3.1 The Ray Payload: The 16-Byte Alignment Rule**

The RayPayload structure is the data packet carried by the ray. It is defined by the user and passed to TraceRay.

#### **3.1.1 Register Pressure and Spilling**

When TraceRay is called, the state of the calling shader (variables needed after the call) and the payload must be preserved. Since the ray traversal is asynchronous and indeterminate in duration, the hardware cannot simply hold this state in active registers indefinitely if the SM needs to switch to other threads. Consequently, this data is often spilled to the stack (local memory).

Large payloads increase the volume of data that must be spilled and reloaded. This increases memory bandwidth consumption and latency. A payload that exceeds the register file limit will cause a dramatic reduction in occupancy, as fewer threads can be active on the SM simultaneously.1

#### **3.1.2 Cache Line Alignment**

Hardware memory fetch units operate on cache lines (typically 128 bytes). HLSL packing rules dictate that data structures should optimize for these boundaries. Specifically, HLSL packs variables into 4-component vectors (16 bytes). A variable cannot straddle a 16-byte boundary.

If a developer defines a payload carelessly:

C++

struct SuboptimalPayload {  
    float3 color;    // 12 bytes  
    float2 uv;       // 8 bytes (Cannot fit in remaining 4 bytes; pushed to next line)  
    float distance;  // 4 bytes  
    // Effective size: 32 bytes with internal padding waste  
};

The hardware may need to fetch more data lines than necessary. By manually padding and aligning to 16 bytes, developers ensure that the payload is fetched in the minimum number of transactions.

**Optimized Packing Strategy:**

C++

struct OptimizedPayload {  
    float3 color;    // 12 bytes  
    float distance;  // 4 bytes (Fills the first 16-byte chunk)  
    float2 uv;       // 8 bytes  
    float2 padding;  // 8 bytes (Explicit padding to align to 16 bytes)  
}; 

Using bit-packing techniques (e.g., packing normals into uint using octahedral encoding, packing UVs into half2) is highly recommended to shrink the payload. Every byte saved reduces register pressure and stack bandwidth.6

### **3.2 Minimizing Live State**

"Live state" refers to any local variable in the Ray Generation shader that is defined *before* a TraceRay call and used *after* it.

High-level shader language

// Example of excessive live state  
float3 throughput \= float3(1,1,1);  
float3 directLight \= CalculateDirectLight(); // Live state  
TraceRay(..., payload); // Hardware must save 'throughput' and 'directLight' to stack  
finalColor \+= payload.radiance \* throughput \+ directLight;

To optimize this, developers should strive to minimize variables crossing the trace boundary. One technique is to move calculations to *after* the trace if they don't depend on pre-trace data, or to recompute cheap values rather than storing them. However, this is a balance; recomputing expensive math is worse than spilling. A common pattern is to structure the loop so that state is naturally minimized at the trace point.1

### **3.3 Recursion: The Stack Killer**

DXR supports recursive ray tracing (calling TraceRay from within a Closest Hit shader). However, this is generally discouraged for performance-critical applications.

Each level of recursion requires the allocation of a new stack frame for the shader execution context. The total stack size required is linear with the maximum recursion depth.

$$\\text{StackSize} \\approx \\text{MaxDepth} \\times (\\text{PayloadSize} \+ \\text{LocalVars})$$  
Large stack sizes deplete VRAM and reduce the effective cache size per thread. The best practice is to set MaxTraceRecursionDepth to **1**. Instead of recursion, use a loop in the Ray Generation shader ("iterative path tracing"). The RayGen shader traces a ray; the Hit shader returns material data (surface normal, albedo) via the payload; the RayGen shader updates the throughput and probability density function (PDF), determines the next ray direction, and loops. This "tail recursion elimination" keeps the stack footprint constant and minimal.1

## ---

**4\. Resource Binding: The Shift to Bindless**

The mechanism for binding resources (textures, buffers, samplers) to shaders has evolved to address the specific needs of ray tracing.

### **4.1 Global vs. Local Root Signatures**

In the initial DXR 1.0 specification, **Local Root Signatures** allowed developers to bind resources specifically to individual shader records in the SBT. For example, a "Wood" material record could point directly to the wood texture SRV.

While conceptually clean, Local Root Signatures introduce significant CPU overhead. Generating the SBT requires the CPU to write unique handles for every object in the scene. Furthermore, accessing these local arguments on the GPU involves pointer indirection, which disrupts the scalar unit and increases latency.9

**Global Root Signatures** are the preferred alternative. A single root signature is bound for the entire DispatchRays call. This signature provides access to all scene resources. The challenge then becomes: how does a generic hit shader know which texture to use for the specific triangle it hit?

### **4.2 Bindless Rendering (SM 6.6)**

The solution is **Bindless Rendering**, fully realized with HLSL Shader Model 6.6. In this model, the application binds the entire descriptor heap of the device to the shader as an unbounded array.

High-level shader language

// HLSL SM 6.6 Resource Access  
Texture2D\<float4\> AlbedoTextures \= ResourceDescriptorHeap;

The geometry instance in the TLAS is assigned a user-defined InstanceID (or a more complex offset calculated from InstanceContributionToHitGroupIndex). The hit shader retrieves this ID and uses it to look up material data (e.g., texture indices) from a global structured buffer. It then uses those indices to sample the correct textures from the bindless heap.

**Advantages:**

1. **CPU Performance:** The CPU no longer needs to update thousands of root tables per frame. It simply uploads the material index buffer.  
2. **SBT Size:** The SBT becomes compact, often containing just a single generic Hit Group record for the entire scene, or one per material *type* (e.g., one for Uber-PBR, one for Hair), rather than one per object.  
3. **Flexibility:** Shaders can access any resource, enabling techniques like "Material gathering" where the RayGen shader evaluates materials for all hits.11

## ---

**5\. Mitigating Divergence: Advanced DXR 1.2 Features**

As ray tracing workloads become more complex, divergence becomes the limiting factor. DXR 1.2 introduces hardware-accelerated features to manage this.

### **5.1 Shader Execution Reordering (SER)**

SER is arguably the most powerful tool for optimizing incoherent path tracing.

The Problem:  
In a path tracer, secondary rays scatter. A warp of 32 threads might result in hits on 10 different materials. Without SER, the warp effectively serializes: it executes Material A (masking off non-A threads), then Material B, and so on.  
The Solution:  
SER allows the developer to inject a sorting stage between intersection and shading.

1. **Trace:** The RayGen shader calls TraceRay (or NvTraceRayHitObject) but requests a "Hit Object" return instead of immediate shading.  
2. **Sort:** The shader calls NvReorderThread (or equivalent API), passing a key (e.g., the Shader ID or Material ID of the hit object). The hardware/scheduler reassigns threads to new warps based on these keys.  
3. **Shade:** The new, coherent warps execute the shading logic.

Impact:  
Benchmarks in engines like Unreal Engine 5 demonstrate performance improvements of roughly 20-50% in complex path tracing scenarios. The overhead of the reordering step is outweighed by the massive gain in SIMT efficiency during the heavy shading phase.3

### **5.2 Opacity Micromaps (OMM)**

Alpha testing is a major bottleneck. Opacity Micromaps allow developers to bake opacity information (1-bit or 2-bit masks) directly onto the geometry. The traversal unit reads these microsurfaces during intersection.

* **2-State OMM:** Defines regions as strictly opaque or transparent. The hardware handles this fully, skipping Any-Hit shaders entirely.  
* **4-State OMM:** Adds an "unknown" state. The hardware handles the opaque/transparent regions and only invokes the Any-Hit shader for the "unknown" (boundary) regions.

Performance Gain:  
This feature is particularly transformative for foliage. By eliminating the Any-Hit tax for the vast majority of the leaf surface (the solid interior and the empty air), frame rates in vegetation-heavy scenes can nearly double.15

## ---

**6\. Execution Models: Pipeline vs. Inline (RayQuery)**

DXR 1.1 introduced RayQuery (Inline Ray Tracing), offering an alternative to the TraceRay pipeline.

### **6.1 RayQuery (Inline)**

RayQuery allows a shader (Compute, Pixel, or RayGen) to traverse the BVH explicitly within a loop, without spawning new shader threads.

* **Pros:**  
  * No overhead from SBT lookup or shader scheduling.  
  * State remains in registers; no stack spilling.  
  * Ideal for simple visibility queries (e.g., shadows, ambient occlusion) where the logic is "is anything blocking this ray?".  
* **Cons:**  
  * It effectively creates a "Mega-Shader." If the scene requires complex handling (e.g., transparency, different intersection logic), all that code must be inlined into the caller, leading to register bloat and instruction cache thrashing.  
  * Hardware scheduling is less effective at hiding latency compared to the pipeline model on some architectures.

### **6.2 TraceRay (Pipeline)**

* **Pros:**  
  * Modular. Materials are separated into distinct shaders.  
  * The hardware scheduler (especially on NVIDIA RTX) can efficiently manage the execution of massive numbers of rays, hiding the latency of one ray's traversal with the execution of another's shading.  
* **Cons:**  
  * Higher overhead for setup and state management.

Recommendation:  
Use a hybrid approach. Use RayQuery for hard shadows and simple visibility checks (e.g., AO). Use TraceRay (Pipeline) for complex radiance evaluation, reflections, and GI where material shading is non-trivial. Note that AMD architectures often show a stronger preference for RayQuery in Compute Shaders due to their wave execution model, while NVIDIA handles both robustly but recommends SER for complex Pipelines.17

## ---

**7\. Sampling Strategies: ReSTIR**

Optimizing the shader is futile if the sampling algorithm requires thousands of rays to converge. Reservoir Spatiotemporal Importance Resampling (ReSTIR) is the standard for real-time ray tracing, allowing high-quality lighting with very few rays per pixel.

### **7.1 Reservoir Packing**

ReSTIR works by maintaining a "reservoir" of light samples for each pixel and reusing them spatially (from neighbors) and temporally (from previous frames).

The bandwidth required to read/write these reservoirs is a bottleneck. HLSL optimization here focuses on packing. A reservoir typically needs:

* Light Index (e.g., 16 bits)  
* Weight sum (float)  
* Target PDF (float)  
* Sample count (M) (8 or 16 bits)

By aggressively packing this data into uint2 or uint4, developers reduce global memory traffic. For instance, the Light Index and Sample Count can be packed into a single uint.

### **7.2 Visibility Reuse (The "El-Cheapo" Trick)**

Tracing shadow rays for every reservoir reuse step is expensive. A common optimization is to assume that if a neighbor pixel's sample is valid for that neighbor, it is likely valid for the current pixel (if they are close in position/normal). By skipping the visibility ray for spatial reuse and amortizing the cost (e.g., checkerboard tracing), developers can cut the ray budget in half with minimal visual bias.20

## ---

**8\. Profiling and Diagnostics**

The complexity of the RTX pipeline requires specialized profiling tools. NVIDIA Nsight Graphics is the industry standard.

### **8.1 Key Performance Counters**

* **Active Threads per Warp:** A direct measure of divergence. If this value is low (e.g., \< 16), it indicates that warps are being masked off heavily. This is a signal to use SER or consolidate shaders.  
* **Warp Stall Reasons \- Long Scoreboard:** Indicates the warp is waiting for memory (usually textures). Improve texture cache locality, use mip-bias, or use bindless to allow the compiler to better schedule loads.  
* **SM Occupancy:** Low occupancy often points to excessive register usage. Check the RayPayload size and reduce local variables in the RayGen shader.

### **8.2 Stack Size Tuning**

The default stack size calculated by DXR is a worst-case estimate.

$$\\text{StackSize} \= \\text{RayGen} \+ \\max(\\text{HitShader}) \\times \\text{Depth}$$

Developers should manually calculate the required stack size based on their specific shader graph and set it via SetPipelineStackSize. Minimizing this allocation frees up VRAM and can improve performance by reducing cache pressure.21

## ---

**9\. Conclusion**

Mastering HLSL for RTX is a discipline of aligning software logic with hardware reality. The days of "write once, run anywhere" are replaced by a need for deep architectural awareness. The highest performance is achieved by a pipeline that:

1. **Traverses** a high-quality BVH populated with opaque-flagged triangles.  
2. **Transfers** state using a minimal, 16-byte aligned payload.  
3. **Executes** coherent shading workloads using SER to manage divergence.  
4. **Accesses** resources via a bindless, global root signature.  
5. **Samples** the light field intelligently using ReSTIR to minimize the ray count.

By systematically addressing each of these stages—from the AS build flags to the bit-packing of the reservoir—developers can unlock the full potential of real-time ray tracing, transitioning from experimental implementations to robust, shipping-quality rendering solutions.

| Feature | Best Practice | Impact |
| :---- | :---- | :---- |
| **Geometry Flags** | Force OPAQUE wherever possible. | Eliminates Any-Hit shader interrupt overhead. |
| **Root Signature** | Use Global Root Signature \+ Bindless (SM 6.6). | Reduces CPU overhead and SBT size. |
| **Payload** | Minimize size; Align to 16 bytes. | Improves register occupancy and cache line efficiency. |
| **Recursion** | Set Max Depth to 1; Use loops. | Reduces stack memory usage and VRAM footprint. |
| **Divergence** | Use Shader Execution Reordering (SER). | Restores warp coherency in complex path tracing. |
| **Alpha Testing** | Use Opacity Micromaps (OMM). | Hardware acceleration for transparency (e.g., foliage). |
| **Shadows** | Use Inline Ray Tracing (RayQuery). | Avoids full pipeline overhead for simple visibility. |

#### **Works cited**

1. Tips and Tricks: Ray Tracing Best Practices | NVIDIA Technical Blog, accessed December 6, 2025, [https://developer.nvidia.com/blog/rtx-best-practices/](https://developer.nvidia.com/blog/rtx-best-practices/)  
2. Best Practices for Using NVIDIA RTX Ray Tracing (Updated), accessed December 6, 2025, [https://developer.nvidia.com/blog/best-practices-for-using-nvidia-rtx-ray-tracing-updated/](https://developer.nvidia.com/blog/best-practices-for-using-nvidia-rtx-ray-tracing-updated/)  
3. Boosting Ray Tracing Performance with Shader Execution Reordering: Introducing VK\_EXT\_ray\_tracing\_invocation\_reorder \- The Khronos Group, accessed December 6, 2025, [https://www.khronos.org/blog/boosting-ray-tracing-performance-with-shader-execution-reordering-introducing-vk-ext-ray-tracing-invocation-reorder](https://www.khronos.org/blog/boosting-ray-tracing-performance-with-shader-execution-reordering-introducing-vk-ext-ray-tracing-invocation-reorder)  
4. Best Practices: Using NVIDIA RTX Ray Tracing (Updated), accessed December 6, 2025, [https://developer.nvidia.com/blog/best-practices-using-nvidia-rtx-ray-tracing/](https://developer.nvidia.com/blog/best-practices-using-nvidia-rtx-ray-tracing/)  
5. Managing Memory for Acceleration Structures in DirectX Raytracing | NVIDIA Technical Blog, accessed December 6, 2025, [https://developer.nvidia.com/blog/managing-memory-for-acceleration-structures-in-dxr/](https://developer.nvidia.com/blog/managing-memory-for-acceleration-structures-in-dxr/)  
6. Packing rules for constant variables \- Win32 apps | Microsoft Learn, accessed December 6, 2025, [https://learn.microsoft.com/en-us/windows/win32/direct3dhlsl/dx-graphics-hlsl-packing-rules](https://learn.microsoft.com/en-us/windows/win32/direct3dhlsl/dx-graphics-hlsl-packing-rules)  
7. Memory Allocation and Bytes Alignment in WebGPU (Ray Tracing Tutorial) \- Medium, accessed December 6, 2025, [https://medium.com/@osebeckley/memory-allocation-and-bytes-alignment-in-webgpu-ray-tracing-tutorial-b53f99385ab3](https://medium.com/@osebeckley/memory-allocation-and-bytes-alignment-in-webgpu-ray-tracing-tutorial-b53f99385ab3)  
8. Self-contained DirectX Raytracing tutorial \- Laura's (Mostly) Unreal Blog, accessed December 6, 2025, [https://landelare.github.io/2023/02/18/dxr-tutorial.html](https://landelare.github.io/2023/02/18/dxr-tutorial.html)  
9. (PDF) CPU Performance in DXR \- ResearchGate, accessed December 6, 2025, [https://www.researchgate.net/publication/354064988\_CPU\_Performance\_in\_DXR](https://www.researchgate.net/publication/354064988_CPU_Performance_in_DXR)  
10. Root Signature Limits \- Win32 apps | Microsoft Learn, accessed December 6, 2025, [https://learn.microsoft.com/en-us/windows/win32/direct3d12/root-signature-limits](https://learn.microsoft.com/en-us/windows/win32/direct3d12/root-signature-limits)  
11. Bindless rendering — Setup \- Traverse Research, accessed December 6, 2025, [https://blog.traverseresearch.nl/bindless-rendering-setup-afeb678d77fc](https://blog.traverseresearch.nl/bindless-rendering-setup-afeb678d77fc)  
12. Read My Chapter in Ray Tracing Gems II\! \- The Danger Zone, accessed December 6, 2025, [https://therealmjp.github.io/posts/rtg2-bindless/](https://therealmjp.github.io/posts/rtg2-bindless/)  
13. Bindless Descriptors \- Wicked Engine, accessed December 6, 2025, [https://wickedengine.net/2021/04/bindless-descriptors/](https://wickedengine.net/2021/04/bindless-descriptors/)  
14. Improve Shader Performance and In-Game Frame Rates with Shader Execution Reordering, accessed December 6, 2025, [https://developer.nvidia.com/blog/improve-shader-performance-and-in-game-frame-rates-with-shader-execution-reordering/](https://developer.nvidia.com/blog/improve-shader-performance-and-in-game-frame-rates-with-shader-execution-reordering/)  
15. Succinct Opacity Micromaps, accessed December 6, 2025, [https://fileadmin.cs.lth.se/graphics/research/papers/2024/succinct\_opacity\_micromaps/paper-author-version.pdf](https://fileadmin.cs.lth.se/graphics/research/papers/2024/succinct_opacity_micromaps/paper-author-version.pdf)  
16. D3D12 Opacity Micromaps \- DirectX Developer Blog, accessed December 6, 2025, [https://devblogs.microsoft.com/directx/omm/](https://devblogs.microsoft.com/directx/omm/)  
17. Inline vs Pipeline Ray Tracing \- Evolve Benchmarking, accessed December 6, 2025, [https://www.evolvebenchmark.com/blog-posts/inline-vs-pipeline-ray-tracing](https://www.evolvebenchmark.com/blog-posts/inline-vs-pipeline-ray-tracing)  
18. Is Ray Query intended as a replacement for Ray Tracing Pipelines? : r/vulkan \- Reddit, accessed December 6, 2025, [https://www.reddit.com/r/vulkan/comments/vph2gt/is\_ray\_query\_intended\_as\_a\_replacement\_for\_ray/](https://www.reddit.com/r/vulkan/comments/vph2gt/is_ray_query_intended_as_a_replacement_for_ray/)  
19. RDNA Performance Guide \- AMD GPUOpen, accessed December 6, 2025, [https://gpuopen.com/learn/rdna-performance-guide/](https://gpuopen.com/learn/rdna-performance-guide/)  
20. ReSTIR El-Cheapo \- Zyanide's Lab, accessed December 6, 2025, [http://www.zyanidelab.com/restir-el-cheapo/](http://www.zyanidelab.com/restir-el-cheapo/)  
21. DirectX Raytracing (DXR) Functional Spec \- Microsoft Open Source, accessed December 6, 2025, [https://microsoft.github.io/DirectX-Specs/d3d/Raytracing.html](https://microsoft.github.io/DirectX-Specs/d3d/Raytracing.html)  
22. specification v0.09 \- Introduction to DirectX RayTracing, accessed December 6, 2025, [https://intro-to-dxr.cwyman.org/spec/DXR\_FunctionalSpec\_v0.09.docx](https://intro-to-dxr.cwyman.org/spec/DXR_FunctionalSpec_v0.09.docx)