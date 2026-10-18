# Directory for building second version of SolTrace API, stapi_v2.

### Build

Expects SolTrace with `SOLTRACE_BUILD_API="ON"` to be built and compiled before building the wheel.
The bindings are not relient on the python interpreter, but the wheel is platform dependent because the shared libraries are shipped. 
Need build dependancies, [build](https://build.pypa.io/en/stable/) and [hatchling](https://hatch.pypa.io/1.13/), run `pip install build hatchling`.
In `SolTrace/api/pysoltrace`, run `py -m build --wheel --outdir ..\dist\`.
For `nvsoltrace`, expects SolTrace to be built with `SOLTRACE_BUILD_OPTIX_SUPPORT="ON`. In `SolTrace/api/nvsoltrace`, run `py -m build --wheel --outdir ..\dist\`.

Uses `pyproject.toml` build the wheel. "A wheel is a ZIP-format archive with a specially formatted file name and the .whl extension." ([see here for more](https://packaging.python.org/en/latest/specifications/binary-distribution-format/)). Once built, a path is set for searching for ptx files to a temporary directory. Rebuilding the SolTrace locally will set that to the typical build path again.

When using the python wrapper, if you see `OSError: [WinError -529697949] Windows Error 0xe06d7363` means either a bad memory operation, i.e. null pointer dereference, or runtime linking error, i.e. can't locate PTX. The error number is some MSVC code, don't know   equivalent error code for gcc/clang on Linux/macOS.

The `Microsoft.CppBuild.targets(548,5): warning MSB8029` is about building in a temporary directory and seems to be harmless.

Note:
View wheel contents on Windows: `Add-Type -A "System.IO.Compression.FileSystem"; [IO.Compression.ZipFile]::OpenRead("C:\abs\path\to\SolTrace\api\dist\pysoltrace-0.0.0-py3-none-win_amd64.whl").Entries.FullName`
Linux: unzip -l pysoltrace-0.0.0-py3-none-manylinux.whl

View dll contents on Windows in Developer PowerShell: `dumpbin /exports stapi_v2.dll`
