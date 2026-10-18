from pysoltrace.api.batch.data.element import element
from pysoltrace.api.batch.data.optic import optic
from pysoltrace.api.batch.data.sun import sun
from pysoltrace.api.batch.data.json import json

class data:
    __slots__ = ('element', 'optic', 'sun', 'json')
    
    def __init__(self, adder):
        self.element = element(adder)
        self.optic   = optic(adder)
        self.sun     = sun(adder)
        self.json    = json(adder)

__all__ = ['data']