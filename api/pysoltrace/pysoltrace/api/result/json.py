import ctypes
from pathlib import Path
from typing import Literal
import orjson

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

# TODO: for when group results are implemented

##############################################################
# functions for simulation data management thru json strings #
##############################################################
class json(context):
    @st_function
    def load(self, input_json: str | dict) -> None:
        # TODO: add Path arg for filename
        # if isinstance(input_json, str):
        #     f = open(input_json, mode='rb')
        #     code = self.__pdll.st_read_input_json(self.__pcxt, f.read())
        #     f.close()
        # else:
        #     code = self.__pdll.st_read_input_json(self.__pcxt, orjson.dumps(input_json))
        # self.__check_return_code(code)

        assert isinstance(input_json, (str, dict, Path)), f'input_json must be a str, dict, or Path, got {type(input_json)}'

        loader = self._pdll.st_read_input_json
        if isinstance(input_json, str):
            f = open(input_json, mode='rb')
            _json = f.read()
            f.close()
        elif isinstance(input_json, dict):
            _json = orjson.dumps(input_json)
        # TODO: implement below
        # else:
        #     _json = input_json.name.encode('utf-8')
        #     loader = self._pdll.st_read_input_json_by_name

        return loader(self._pcxt, _json)

    @st_function
    def dump(self, filename: str):
        return self._pdll.st_export_json_file(self._pcxt, filename.encode())