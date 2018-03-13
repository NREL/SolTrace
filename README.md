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

    * Windows: Visual Studio 2017 Community or other editions available at [https://www.visualstudio.com/](https://www.visualstudio.com/).
    * Linux: g++ compiler available at [http://www.cprogramming.com/g++.html](http://www.cprogramming.com/g++.html) or as part of the Linux distribution.

2. Download the wxWidgets 3.1.0 source code for your operating system from [https://www.wxwidgets.org/downloads/](https://www.wxwidgets.org/downloads/).

3. Build wxWidgets.

4. In Windows, create the WXMSW3 environment variable on your computer to point to the wxWidgets installation folder, or Linux, create the dynamic link `/usr/<USERNAME>/local/bin/wx-config-3` to point to `/path/to/wxWidgets/bin/wx-config`.

5. As you did for wxWidgets, for each of the following projects, clone (download) the repository, build the project, and then (Windows only) create an environment variable pointing to the project folder. Build the projects in the following order, and assign the environment variable for each project before you build the next one:

<table>
<tr><th>Project</th><th>Repository URL</th><th>Windows Environment Variable</th></tr>
<tr><td>LK</td><td>https://github.com/NREL/lk</td><td>LKDIR</td></tr>
<tr><td>WEX</td><td>https://github.com/NREL/wex</td><td>WEXDIR</td></tr>
</table>

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