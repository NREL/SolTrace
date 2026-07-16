# OptixCSP


**OptixCSP** is a modern, cross-platform (Windows and Linux) GPU‑accelerated ray‑tracing middleware library for modeling and simulating optical systems in engineering applications, such as **Concentrated Solar Power (CSP)**. It leverages [NVIDIA OptiX](https://developer.nvidia.com/optix) and can be used as a stand‑alone library or integrated into CSP workflows (e.g., [SolTrace](https://github.com/NREL/SolTrace)) to significantly accelerate ray‑tracing tasks.

OptixCSP focuses on accuracy, performance, and extensibility so you can:

* Rapidly simulate optical behavior in CSP systems: heliostats, receivers, sun modeling, optical properties, etc.
* Leverage GPU acceleration for complex scenes, large‑scale simulations, and dynamic/annual analyses.
* Integrate OptixCSP into existing CSP software as a high‑performance backend.
* Build standalone ray‑tracing applications using the OptixCSP core.


## Prerequisites

### Hardware & Drivers

* NVIDIA Turing-class or newer GPU with compute capability 7.5 or higher.
  RTX3050, RTX4070Ti, and H100 have been tested.
* NVIDIA driver, compatible with your target OptiX & CUDA toolkit.

### Software

* **CUDA Toolkit**: 12.0 or newer.
* **NVIDIA OptiX SDK**: 8.x or newer.
* **CMake**: 3.18 or higher.
* **C++ Compiler**: C++17 or newer.

> **Note:** The **runtime** (`libnvoptix.so` on Linux, `nvoptix.dll` on Windows) is provided by the NVIDIA driver. The OptiX SDK provides headers only. Make sure the optix version is supported by the driver on your machine.

### Verified build/test configurations

* **Windows**: MS Visual Studio 2022, CUDA 12.8, OptiX 9.0
* **Windows**: MS Visual Studio 2022, CUDA 12.8, OptiX 8.1
* **Windows**: MS Visual Studio 2022, CUDA 12.3, OptiX 8.1
* **Linux (redhat 8.0)**: gcc 11.2, CUDA 12.3, Optix 8.0
* **Linux (redhat 8.0)**: gcc 12.1, CUDA 12.3, Optix 8.0
* **Linux (ubuntu 22.04)**: gcc 11.4, CUDA 12.8, Optix 9.0

### Quick checks (Linux)

Use these to confirm the system is ready before building:

```bash
# 1) GPU & driver
nvidia-smi
# Look for your GPU and the CUDA Version (max toolkit version supported by the driver)

# 2) OptiX runtime library from the driver, you should see libnvoptix.so.1
ldconfig -p | grep nvoptix

# 3) CUDA compiler (if installed on the node)
nvcc --version
```
### Installing the OptiX SDK

Download the OptiX SDK from NVIDIA ([link here](https://developer.nvidia.com/designworks/optix/download)). If your NVIDIA driver does not support the latest OptiX, please download the compatible version from the [legacy download page](https://developer.nvidia.com/designworks/optix/downloads/legacy). For Windows, follow the instructions; on Linux, run the provided script.


```bash
chmod +x OptiX-SDK-8.0.x-linux64.sh
./OptiX-SDK-8.0.x-linux64.sh --skip-license --prefix="$HOME/Optix8.0"
```

This creates a tree like:

```
$HOME/Optix8.0/
  include/
  SDK/
```

## Building OptixCSP
Clone the repo, 
```bash
git clone git@github.com:uwsbel/optixCSP.git
cd optixCSP
```

### Linux

```bash
# configure cmake
mkdir -p build
cd build
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DOptiX_INCLUDE="$HOME/Optix8.0/include" \
  -DOptiX_INSTALL_DIR="$HOME/Optix8.0/"
  ..

# build
cmake --build . -j
```
>**Notes**: On kestrel, consider loading `CUDA Toolkit` and `gcc` first before running cmake: `module load cuda/12.4 gcc/12`. 

### Windows
* In CMake GUI, set `OptiX_INCLUDE` to the Optix SDK's `include/` folder (e.g., `C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.1.0/include`).
* Generate for Visual Studio 2022, then build the `Release` and `Debug` configuration.

The build targets CUDA compute capability 7.5 by default, independent of the
GPU installed on the build host. Set
`SOLTRACE_OPTIX_COMPUTE_CAPABILITY` to a numeric architecture such as `86` only
when intentionally producing a specialized build with a higher minimum.


## Running the demos

After building, demos/executables are placed under your build output (e.g., `build/bin/`). Run them from the command line or your IDE. 

## Using OptixCSP as a library

You can link `OptixCSP` to third‑party software such as **SolTrace**.

1. Build and install (or just build) OptixCSP per the instructions above. The build exports a package config file `OptixCSPConfig.cmake` in its build tree (e.g., `build/cmake/`).

2. In your project’s `CMakeLists.txt`, configure dependencies and find OptixCSP:

```cmake
cmake_minimum_required(VERSION 3.18)
project(MyApp CXX)

find_package(CUDAToolkit REQUIRED)  # Provides CUDA::cuda_driver
find_package(OptiX REQUIRED) # Provides optix headers

# Tell CMake where to find OptixCSP's package config
set(OptixCSP_DIR "/path/to/OptixCSP/build/cmake")
find_package(OptixCSP REQUIRED)

add_executable(MyApp main.cpp)

target_link_libraries(MyApp
  PUBLIC
    OptixCSP::OptixCSP_core
    CUDA::cuda_driver
)

# If you need headers from the OptiX SDK directly:
# target_include_directories(MyApp PRIVATE "$ENV{HOME}/optix-8.0/include")
```
