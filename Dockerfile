# OpenRTX Shader Build Environment
# Version pins match .github/workflows/build-shader.yml

FROM ubuntu:22.04 AS builder

# Version pins (sync with CI/CD)
ARG DXC_VERSION=1.8.2505.1
ARG LAZURITE_VERSION=0.7.0
ARG VANILLA_VERSION=v23

# Build arg for vanilla materials URL (override at build time)
ARG VANILLA_BASE_URL

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    curl \
    wget \
    unzip \
    zip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment and install Lazurite
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir lazurite==${LAZURITE_VERSION}

# Download and install DXC (DirectX Shader Compiler)
RUN mkdir -p /opt/dxc && \
    curl -sL "https://github.com/microsoft/DirectXShaderCompiler/releases/download/v${DXC_VERSION}/dxc_2025_07_14.zip" -o /tmp/dxc.zip && \
    unzip -q /tmp/dxc.zip -d /opt/dxc && \
    rm /tmp/dxc.zip
ENV DXC_PATH=/opt/dxc/bin/x64/dxc

# Download shaderc (veka0's bgfx-mcbe version for Minecraft PostFX shaders)
RUN curl -sL "https://github.com/veka0/bgfx-mcbe/releases/download/binaries/shaderc" -o /usr/local/bin/shaderc && \
    chmod +x /usr/local/bin/shaderc
ENV SHADERC_PATH=/usr/local/bin/shaderc

# Set working directory
WORKDIR /shader

# Copy project files
COPY . .

# Download vanilla materials if URL provided
RUN if [ -n "$VANILLA_BASE_URL" ]; then \
        mkdir -p vanilla && \
        echo "Downloading vanilla materials from provided URL..." && \
        curl -sL "${VANILLA_BASE_URL}/RTXStub.material.bin" -o vanilla/RTXStub.material.bin && \
        curl -sL "${VANILLA_BASE_URL}/RTXPostFX.Tonemapping.material.bin" -o vanilla/RTXPostFX.Tonemapping.material.bin && \
        curl -sL "${VANILLA_BASE_URL}/RTXPostFX.Bloom.material.bin" -o vanilla/RTXPostFX.Bloom.material.bin && \
        echo "Downloaded vanilla materials:" && ls -la vanilla/; \
    else \
        echo "VANILLA_BASE_URL not provided - skipping vanilla materials download"; \
        echo "Mount vanilla/ directory or provide VANILLA_BASE_URL at build time"; \
    fi

# Run template generation script if it exists
RUN if [ -f "src/script.py" ]; then \
        cd src && python script.py && cd ..; \
    fi

# Build shaders
RUN mkdir -p build && \
    if [ -d "vanilla" ] && [ "$(ls -A vanilla 2>/dev/null)" ]; then \
        lazurite build project/ -o build/ --shaderc $SHADERC_PATH --dxc $DXC_PATH && \
        echo "Build complete:" && ls -la build/; \
    else \
        echo "Vanilla materials not found - cannot build"; \
        exit 1; \
    fi

# Final stage - minimal image with just the built shaders
FROM alpine:latest AS output

WORKDIR /output
COPY --from=builder /shader/build/*.material.bin ./

# Verify output
RUN echo "Built shader files:" && ls -la /output/

# Default command shows built files
CMD ["ls", "-la", "/output"]
