# OpenRTX Shader Development Container
# Mirrors the CI/CD workflow environment for local development

FROM ubuntu:22.04

# Version pins (match .github/workflows/build-shader.yml)
ARG DXC_VERSION=1.8.2505.1
ARG LAZURITE_VERSION=0.7.0
ARG PYTHON_VERSION=3.12

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    curl \
    wget \
    unzip \
    zip \
    git \
    ca-certificates \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    python${PYTHON_VERSION} \
    python${PYTHON_VERSION}-venv \
    python${PYTHON_VERSION}-dev \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.12 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python${PYTHON_VERSION} 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python${PYTHON_VERSION} 1

# Create virtual environment and install Lazurite
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir lazurite==${LAZURITE_VERSION}

# Download and install DXC (DirectX Shader Compiler)
RUN mkdir -p /opt/dxc && \
    curl -sL "https://github.com/microsoft/DirectXShaderCompiler/releases/download/v${DXC_VERSION}/dxc_2025_07_14.zip" -o /tmp/dxc.zip && \
    unzip -q /tmp/dxc.zip -d /opt/dxc && \
    rm /tmp/dxc.zip && \
    chmod +x /opt/dxc/bin/x64/dxc || true
ENV DXC_PATH=/opt/dxc/bin/x64/dxc

# Download shaderc (veka0's bgfx-mcbe version for Minecraft shaders)
RUN curl -sL "https://github.com/veka0/bgfx-mcbe/releases/download/binaries/shaderc" -o /usr/local/bin/shaderc && \
    chmod +x /usr/local/bin/shaderc
ENV SHADERC_PATH=/usr/local/bin/shaderc

# Create workspace directory
WORKDIR /workspaces/mcrtx-shader-template

# Create vanilla materials directory (to be populated via secrets or manual download)
RUN mkdir -p /workspaces/mcrtx-shader-template/vanilla

# Environment variables for build scripts
ENV VANILLA_VERSION=v23
ENV DXC_VERSION=${DXC_VERSION}
ENV LAZURITE_VERSION=${LAZURITE_VERSION}

# Verify installations
RUN echo "=== Environment ===" && \
    echo "Python: $(python3 --version)" && \
    echo "Lazurite: $(pip show lazurite | grep Version)" && \
    echo "Shaderc: $(shaderc --help 2>&1 | head -1 || echo 'installed')" && \
    echo "DXC: $(ls -la /opt/dxc/bin/x64/ 2>/dev/null || echo 'check platform')"

# Default shell
CMD ["/bin/bash"]
