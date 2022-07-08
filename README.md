# SolTrace

The SolTrace Open Source Project repository contains the source code, tools, and instructions to build a desktop version of the National Renewable Energy Laboratory's SolTrace. SolTrace is a software tool developed at NREL to model concentrating solar power (CSP) systems and analyze their optical performance. Although ideally suited for solar applications, the code can also be used to model and characterize many general optical systems. The creation of the code evolved out of a need to model more complex solar optical systems than could be modeled with existing tools. For more details about SolTrace's capabilities, see the [SolTrace website](https://www.nrel.gov/csp/soltrace.html). For details on integration with SAM, see the [SAM website](https://sam.nrel.gov).

The desktop version of SolTrace for Windows or Linux builds from the following open source projects:

* [LK](https://github.com/nrel/lk) is a scripting language that is integrated into SAM and allows users to add functionality to the program.

* [wxWidgets](https://www.wxwidgets.org/) is a cross-platform graphical user interface platform used for SAM's user interface, and for the development tools included with SSC (SDKtool) and LK (LKscript). The current version of SAM uses wxWidgets 3.1.0.

* [WEX](https://github.com/nrel/wex) is a set of extensions to wxWidgets for custom user-interface elements used by SAM, and by LKscript and DView, which are integrated into SAM.

* This repository, **SolTrace**, provides the user interface to assign values to inputs of the computational modules, run the modules in the correct order, and display calculation results. It also includes tools for editing LK scripts and viewing ray intersection and flux map data.

## Quick Steps for Building SolTrace

For detailed build instructions see the [wiki](https://github.com/NREL/SolTrace/wiki), with specific instructions for:

* [Windows](https://github.com/NREL/SolTrace/wiki/build-windows)
* [OSX](https://github.com/NREL/SolTrace/wiki/build-osx)
* [Linux](https://github.com/NREL/SolTrace/wiki/build-linux)

These are the general quick steps you need to follow to set up your computer for developing SolTrace:

1. Set up your development tools:

    * Windows: Visual Studio 2019 Community or other editions available at [https://www.visualstudio.com/](https://www.visualstudio.com/).
    * Linux: g++ compiler available at [http://www.cprogramming.com/g++.html](http://www.cprogramming.com/g++.html) or as part of the Linux distribution.

2. Download and install CMake 3.19 or higher from [https://cmake.org/download/](https://cmake.org/download/) with the ```Add CMake to the System Path for ...``` option selected.

3. Download the wxWidgets 3.1.5 source code for your operating system from [https://www.wxwidgets.org/downloads/](https://www.wxwidgets.org/downloads/).

4. Build wxWidgets.

5. In Windows, create the WXMSW3 environment variable on your computer to point to the wxWidgets installation folder, or Linux, create the dynamic link `/usr/<USERNAME>/local/bin/wx-config-3` to point to `/path/to/wxWidgets/bin/wx-config`.

6. As you did for wxWidgets, for each of the following projects, clone (download) the repository and then (Windows only) create an environment variable pointing to the project folder. 

<table>
<tr><th>Project</th><th>Repository URL</th><th>Windows Environment Variable</th></tr>
<tr><td>LK</td><td>https://github.com/NREL/lk</td><td>LKDIR</td></tr>
<tr><td>WEX</td><td>https://github.com/NREL/wex</td><td>WEXDIR</td></tr>
</table>

7. Run CMake to create the project build files
    1. Copy the file ```parent-dir-CMakeLists.txt``` into the parent directory also containing ```soltrace/ lk/ wex/``` and ```wxwidgets-3.x.x/``` folders.
    
    2. Rename this file to ```CMakeLists.txt``` before running cmake. You may need to temporarily rename any other file in this directory with the same name. 
    
        E.g., the file should be at ```C:/stdev/CMakeLists.txt```

    3. Create a directory in the main parent folder to store the build files. 
    E.g., ```C:/stdev/build-soltrace/```
    
    4. Open a shell or command window in the build folder from step 3

    5. Copy the following cmake command to the shell and run. Replace the cmake target with a [supported generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#manual:cmake-generators(7))
    
        ```> cmake -G "Visual Studio 16 2019" -DCMAKE_CONFIGURATION_TYPES="Debug;Release" -DCMAKE_SYSTEM_VERSION=10.0 -DSAM_SKIP_TOOLS=1 .. ```

    6. Confirm the project files built. If running visual studio, you should see a ```soltrace_ui.sln``` file in the build-soltrace/ directory.
    
    7. Build all files. The output is stored in the soltrace repository folder, e.g., ```C:/stdev/soltrace/app/deploy/soltrace.exe```. 

        Note that output is NOT stored in the ```build-soltrace/``` directory!

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