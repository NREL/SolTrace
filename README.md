# SolTrace<sup>TM</sup>
![Build](https://github.com/NLR-SolTrace/SolTrace/actions/workflows/CI.yml/badge.svg)[![Coverage Status](https://coveralls.io/repos/github/NLR-SolTrace/SolTrace/badge.svg)](https://coveralls.io/github/NLR-SolTrace/SolTrace?branch=develop)


The SolTrace<sup>TM</sup> Open Source Project repository contains the source code, tools, and instructions to build a desktop version of the National Renewable Energy Laboratory's SolTrace. SolTrace is a software tool developed at NREL to model concentrating solar power (CSP) systems and analyze their optical performance. Although ideally suited for solar applications, the code can also be used to model and characterize many general optical systems. The creation of the code evolved out of a need to model more complex solar optical systems than could be modeled with existing tools. For more details about SolTrace's capabilities, see the [SolTrace website](https://www.nrel.gov/csp/soltrace.html). For details on integration with SAM, see the [SAM website](https://sam.nrel.gov).

## SolTrace Qt GUI

SolTrace includes a Qt-based graphical interface in [`gui/`](gui/). The legacy wxWidgets/LK/WEX user-interface build instructions have been retired from this README. For current GUI requirements, Qt Creator workflow, command-line build instructions, deployment notes, and translation guidance, see the [SolTrace Qt GUI README](gui/README.md).

## SolTrace Python API

Users can also run simulations via the [`pysoltrace`](https://github.com/NREL/SolTrace/blob/develop/app/deploy/api/pysoltrace.py) SolTrace Python API found in the folder `app/deploy/api`. Example files for running the API are found in the `app/deploy/api/examples` subfolder. Documentation is available in HTML or PDF format in the corresponding API subfolder. 

The `pysoltrace` API is capable of running multi-threaded simulations, generating flux maps, creating 3D interactive trace plots, and provides other capabilities that are found in the SolTrace graphical interface. The functionality and flexibility of the API generally exceeds that of the graphical interface. 

The API requires the compiled coretrace library. Project files for building this library are generated using CMake as outlined in the core build steps below. To work without the graphical interface, leave the GUI build option disabled and build the coretrace API target from the generated project.

## Steps for Building SolTrace Core

Building SolTrace core libraries and command-line targets requires a C++17-capable compiler and CMake 3.19 or greater. Once these are available, building can be done in the normal pattern of configure and build:

```sh
git clone https://github.com/NatLabRockies/SolTrace.git
cd SolTrace
mkdir build
cd build
cmake ..
cmake --build . -j4
```

Note the `-j4` instructs cmake to use 4 processes to compile the source code. This can be adjusted if the number of processors available is greater (or less) than 4.

### Building with Intel's Embree Ray Tracing Library

Information about Embree can be found on the [Embree webpage](https://www.embree.org/) or on the [Embree github page](https://github.com/RenderKit/embree).

Prior to building SolTrace, you will need to install Embree v4.x.x. On Linux this is best done with a package manager. E.g.,

```sh
sudo apt-get update
sudo apt-get install libembree-dev
```
On Mac, you can use Homebrew

```sh
brew install embree
```

On Windows (this works for Linux and MacOS as well), you need to download the binaries from the [github page](https://github.com/RenderKit/embree) (under the appropriate installation header). Follow the corresponding install instructions found there making sure to add the location of the embree DLL's to your system path.

Once Embree is installed, clone the SolTrace repo, configure with embree enabled, and build:

```sh
git clone https://github.com/NatLabRockies/SolTrace.git
cd SolTrace
mkdir build
cd build
cmake .. -DSOLTRACE_BUILD_EMBREE_SUPPORT=ON
cmake --build . -j4
```

If cmake is having trouble locating the Embree install, you specify Embree's install location passing the `embree_DIR` variable to cmake. In this case the configure command would be

```sh
cmake .. -DSOLTRACE_BUILD_EMBREE_SUPPORT=ON -Dembree_DIR=<EMBREE_INSTALL_DIR>
```

### Building with Nvidia's Optix Ray Tracing Library

SolTrace includes an OptiX-based runner built for GPU-accelerated ray tracing.

#### Prerequisites

* NVIDIA Turing-class or newer GPU with compute capability 7.5 or higher.
* NVIDIA driver compatible with your target CUDA toolkit and OptiX SDK.
* CUDA Toolkit 12.0 or newer.
* NVIDIA OptiX SDK 8.x or newer.
* CMake 3.19 or newer.
* C++17-capable compiler.

Verified build and test configurations from GPU runner development include:

* Windows: Visual Studio 2022, CUDA 12.8, OptiX 9.0
* Windows: Visual Studio 2022, CUDA 12.8, OptiX 8.1
* Windows: Visual Studio 2022, CUDA 12.3, OptiX 8.1
* Linux (Red Hat 8.0): gcc 11.2, CUDA 12.3, OptiX 8.0
* Linux (Red Hat 8.0): gcc 12.1, CUDA 12.3, OptiX 8.0
* Linux (Ubuntu 22.04): gcc 11.4, CUDA 12.8, OptiX 9.0

> Note: The OptiX runtime library is provided by the NVIDIA driver. The OptiX SDK provides the headers used at build time.

#### Quick checks

Before building, confirm the GPU, driver, and CUDA toolchain are visible:

```sh
nvidia-smi
nvcc --version
```

On Linux, you can also confirm the OptiX runtime library is available from the driver:

```sh
ldconfig -p | grep nvoptix
```

#### Install the CUDA Toolkit

Install the CUDA Toolkit from NVIDIA before configuring SolTrace: <https://developer.nvidia.com/cuda-downloads>.

Use `nvidia-smi` to check the maximum CUDA version supported by your installed NVIDIA driver, then choose a compatible CUDA Toolkit release.

#### Install the OptiX SDK

Download the OptiX SDK from NVIDIA: <https://developer.nvidia.com/designworks/optix/download>.

If your NVIDIA driver does not support the latest OptiX release, use the legacy downloads page: <https://developer.nvidia.com/designworks/optix/downloads/legacy>.

#### Configure and build SolTrace with OptiX support

Linux example:

```sh
git clone https://github.com/NREL/SolTrace.git
cd SolTrace
mkdir build
cd build
cmake .. \
    -DSOLTRACE_BUILD_OPTIX_SUPPORT=ON \
    -DOptiX_INSTALL_DIR=/path/to/NVIDIA-OptiX-SDK
cmake --build . -j4
```

Windows example:

```bat
mkdir build
cd build
cmake .. ^
    -G "Visual Studio 17 2022" ^
    -DSOLTRACE_BUILD_OPTIX_SUPPORT=ON ^
    -DOptiX_INSTALL_DIR="C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.1.0"
cmake --build . --config Release -j
```

If CMake cannot find the OptiX SDK automatically, set `OptiX_INSTALL_DIR` to the SDK root containing `include/optix.h`, or add that location to `CMAKE_PREFIX_PATH`.

OptiX builds target compute capability 7.5 by default. This produces native
code for Turing GPUs and virtual PTX that newer NVIDIA GPUs can compile at run
time. Build machines do not need an NVIDIA GPU, and their installed hardware
does not change the resulting artifact. To intentionally raise the minimum for
a specialized build, configure with, for example,
`-DSOLTRACE_OPTIX_COMPUTE_CAPABILITY=86`.

## Contributing

If you would like to report an issue with SolTrace or make a feature request, please let us know by adding a new issue on the [issues page](https://github.com/NREL/SolTrace/issues).

If you would like to submit code to fix an issue or add a feature, you can use GitHub to do so. Please see [Contributing](CONTRIBUTING.md) for instructions.

## License

SolTrace's open source code is copyrighted by the Alliance for Sustainable Energy and licensed under a [mixed MIT and GPLv3 license](LICENSE.md). It allows for-profit and not-for-profit organizations to develop and redistribute software based on SolTrace under terms of an MIT license and requires that research entities including national laboratories, colleges and universities, and non-profit organizations make the source code of any redistribution publicly available under terms of a GPLv3 license.

## Citing SolTrace

We appreciate your use of SolTrace, and ask that you appropriately cite the software in exchange for its open-source publication. Please use one of the following references in documentation that you provide on your work. For general usage citations, the preferred option is:

> Wendelin, T. (2003). "SolTRACE: A New Optical Modeling Tool for Concentrating Solar Optics." Proceedings of the ISEC 2003: International Solar Energy Conference, 15-18 March 2003, Kohala Coast, Hawaii. New York: American Society of Mechanical Engineers, pp. 253-260; NREL Report No. CP-550-32866.

For citations in work that involves substantial development or extension of the existing code, the preferred option is:

> Wendelin, T., Wagner, M.J. (2018). "SolTrace Open-Source Software Project: [github.com/NREL/SolTrace](https://github.com/NREL/SolTrace)". National Renewable Energy Laboratory. Golden, Colorado.
