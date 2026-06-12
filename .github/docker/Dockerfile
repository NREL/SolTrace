# Start from NVIDIA's official CUDA development image matching your host GPU runner (Ubuntu 24.04 & CUDA 13.1)
FROM nvidia/cuda:13.1.0-devel-ubuntu24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install compiler tools, CMake, and git 
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

# 1. Copy the raw .sh installer file into a temporary container directory
COPY NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64.sh /tmp/NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64.sh

# 2. Make it executable, run it to install to /opt, and delete the installer file to save space
RUN chmod +x /tmp/NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64.sh && \
    /tmp/NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64.sh --skip-license --prefix=/opt && \
    rm /tmp/NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64.sh

# Set environment variable so CMake can automatically locate the SDK
# Note: The installer extracts into a subdirectory named "NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64" inside /opt
ENV OptiX_INSTALL_DIR=/opt/NVIDIA-OptiX-SDK-9.0.0-linux64-x86_64
