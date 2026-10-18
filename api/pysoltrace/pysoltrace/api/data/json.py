from pathlib import Path, WindowsPath
import orjson

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

##############################################################
# functions for simulation data management thru json strings #
##############################################################
class json(context):
    @st_function
    def load(self, input_json: bytes | str | dict | Path) -> None:
        assert isinstance(input_json, (bytes | str, dict, Path)), f'input_json must be a str, dict, or Path, got {type(input_json)}'

        _json = input_json
        loader = self._pdll.st_read_input_json
        if isinstance(input_json, str):    _json = input_json.encode('utf-8')
        elif isinstance(input_json, dict): _json = orjson.dumps(input_json)
        elif isinstance(input_json, Path):
            _json = bytes(input_json.resolve())
            loader = self._pdll.st_read_input_json_file

        return loader(self._pcxt, _json)

    @st_function
    def dump(self, filename: str):
        assert isinstance(filename, (str, Path, bytes)), f'filename must be a str, bytes, or Path, got {type(filename)}'

        _filename = filename
        if isinstance(_filename, str): _filename = _filename.encode()
        elif isinstance(_filename, Path): _filename = str(_filename).encode()

        return self._pdll.st_export_json_file(self._pcxt, _filename)