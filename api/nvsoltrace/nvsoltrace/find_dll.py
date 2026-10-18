import pathlib, sys

def find_dll() -> pathlib.Path:
    _path = pathlib.Path(__file__).parent.resolve() / 'bin'
    assert _path.exists(), f"Could not find SolTrace DLL at {_path}"
                
    if sys.platform == "win32":
        _lib_name = "stapi_v2.dll"
    elif sys.platform == "darwin":
        _lib_name = "libstapi_v2.dylib"
    else:
        _lib_name = "libstapi_v2.so"

    return _path / _lib_name