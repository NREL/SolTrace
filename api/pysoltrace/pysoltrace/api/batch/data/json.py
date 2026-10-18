from pathlib import Path
import orjson

from pysoltrace import dot_h
from pysoltrace.api.batch.utils import batcher, generate_api_call

##############################################################
# functions for simulation data management thru json strings #
##############################################################
class json(batcher):
    def load(self, input_json: bytes | str | dict | Path) -> None:
        assert isinstance(input_json, (bytes | str, dict, Path)), \
            f'input_json must be a str, dict, or Path, got {type(input_json)}'

        _json = input_json
        loader = dot_h.st_api_call.CALL_ST_READ_INPUT_JSON
        if isinstance(input_json, str):    _json = input_json.encode('utf-8')
        elif isinstance(input_json, dict): _json = orjson.dumps(input_json)
        elif isinstance(input_json, Path):
            _json = bytes(input_json.resolve())
            loader = dot_h.st_api_call.CALL_ST_READ_INPUT_JSON_FILE

        return self.adder(generate_api_call(loader, _json))

    def dump(self, filename: str):
        return self.adder(generate_api_call(dot_h.st_api_call.CALL_ST_EXPORT_JSON_FILE,
                                            filename.encode()))