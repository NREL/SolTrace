import ctypes

from pysoltrace import dot_h
from pysoltrace.api.dll import context
from pysoltrace.api.utils import check_return_code, st_function
from pysoltrace.api.batch.parameters import parameters
from pysoltrace.api.batch.data import data
from pysoltrace.api.batch.runner import runner
from pysoltrace.api.batch.result import result
from pysoltrace.api.batch.legacy import legacy

class batch_record:
    __slots__ = ('ready', 'batch_call', '__get_res')

    def __init__(self, batch_call, get_res = None):
        self.ready = False
        self.batch_call = batch_call
        self.__get_res = get_res if get_res != None else lambda: None

    @property
    def value(self): return self.__get_res() if self.ready else None

    def __repr__(self):
        call = dot_h.st_api_call(self.batch_call.type)
        return f'{call.name} ({self.ready}): {self.value}'

class batch(context):
    def __init__(self, pdll, pcxt):
        super().__init__(pdll, pcxt)

        self.parameters = parameters(self.add)
        self.data       = data(self.add)
        self.runner     = runner(self.add)
        self.result     = result(self.add)
        self.legacy     = legacy(self.add)

        self.__called = 0
        # TODO: speed test without the internal list, ie make __call__ take a list
        self.__calls: list[batch_record] = []

    def add(self, api_call, get_res = None):
        br = batch_record(api_call, get_res)
        self.__calls.append(br)
        return br

    def clear(self):
        self.__called = 0
        self.__calls: list[batch_record] = []

    def __call__(self, verbose = False):
        num_calls = len(self.__calls) - self.__called
        assert num_calls > 0, "No batch calls recorded to call"

        callable_brs = self.__calls[self.__called:]
        void_cast = lambda c: ctypes.cast(ctypes.byref(c), ctypes.c_void_p)
        arr = lambda calls, num: (ctypes.c_void_p * num)(*[
            ctypes.cast(ctypes.byref(br.batch_call), ctypes.c_void_p) for br in calls
        ])
        fail_iteration = ctypes.c_uint(0)

        code = self._pdll.st_batch(self._pcxt,
                                   arr(callable_brs, num_calls),
                                   num_calls,
                                   ctypes.byref(fail_iteration),
                                   verbose)
        check_return_code(code, fail_iteration.value)
        for br in callable_brs: br.ready = True
        self.__called += num_calls

__all__ = ['batch', 'batch_record']