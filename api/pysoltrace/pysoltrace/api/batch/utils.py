from collections.abc import Callable
from pysoltrace import dot_h

def generate_api_call(call_type: int, *args):
    assert call_type < dot_h.st_api_call.API_CALL_COUNT, \
        f'Invalid st_api_v2 batch call ({call_type}).'

    rt = dot_h.st_api_call_args()
    rt.type = call_type
    args_name, args_cls = rt.payload._fields_[call_type]
    rt.payload = dot_h.st_api_call_args.payload.type(**{
        args_name: getattr(dot_h, args_cls.__name__)(*args)
    })

    return rt

class batcher:
    def __init__(self, adder: Callable):
        self.adder: Callable = adder