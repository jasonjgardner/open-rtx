# OpenRTX Shader Build Environment
# Based on jasongardner/dxc for DirectX Shader Compiler

FROM jasongardner/dxc:latest AS dxc

# Install Python and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    curl \
    zip \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment and install Lazurite
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir lazurite

# Set working directory
WORKDIR /shader

# Copy project files
COPY . .

# Download vanilla materials (can be overridden at build time)
ARG VANILLA_BASE_URL=https://cdn.bedrock.graphics/vanilla/v23
RUN mkdir -p vanilla && \
    curl -sL "${VANILLA_BASE_URL}/RTXStub.material.bin" -o vanilla/RTXStub.material.bin && \
    curl -sL "${VANILLA_BASE_URL}/RTXPostFX.material.bin" -o vanilla/RTXPostFX.material.bin && \
    curl -sL "${VANILLA_BASE_URL}/RTXPostFX.Tonemapping.material.bin" -o vanilla/RTXPostFX.Tonemapping.material.bin && \
    curl -sL "${VANILLA_BASE_URL}/RTXPostFX.Bloom.material.bin" -o vanilla/RTXPostFX.Bloom.material.bin

# Build shaders
RUN mkdir -p build && \
    lazurite build project/ -o build/

# Final stage - minimal image with just the built shaders
FROM alpine:latest AS output

WORKDIR /output
COPY --from=dxc /shader/build/*.material.bin ./

# Default command shows built files
CMD ["ls", "-la", "/output"]
