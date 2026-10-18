from pysoltrace.api.data.element import element
from pysoltrace.api.data.optic import optic
from pysoltrace.api.data.sun import sun
from pysoltrace.api.data.json import json

class data:
    __slots__ = ('element', 'optic', 'sun', 'json')
    
    def __init__(self, pdll, pcxt):
        self.element = element(pdll, pcxt)
        self.optic   = optic(pdll, pcxt)
        self.sun     = sun(pdll, pcxt)
        self.json    = json(pdll, pcxt)

__all__ = ['data']