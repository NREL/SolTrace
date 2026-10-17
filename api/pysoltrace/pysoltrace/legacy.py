# GOAL: recreate the python side of this with STAPIv2, use the link below to test implementation
# https://github.com/uwsbel/SolTrace/blob/5e08741c1da9e6e68cef5232d9da401586674618/app/deploy/api/examples/tower-demo.py

# from datetime import datetime
# import sys, os, copy
# from typing import List
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy
# from ctypes import *
# c_number = c_double   #must be either c_double or c_float depending on coretrace definition
# import multiprocessing
# import time
# import math
# import random

# # Callback to print command line progress messages
# @CFUNCTYPE(c_int, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32, c_void_p)
# def api_callback(ntracedtotal, ntraced, ntotrace, curstage, nstages, thread_id):
#     if thread_id:
#         if thread_id > 1:
#             return 1

#     w = 50
#     prog = 100.*float(ntraced)/float(ntotrace)
#     pbprog = int(prog * w / 100.)
#     pbar = pbprog *"▮" + (w-pbprog)*"▯"

#     print("{:s}  | Stage ({:d}/{:d}) - Complete {:.2f}%".format(pbar, curstage, nstages, prog), end='\r')
#     return 1

# @CFUNCTYPE(c_int, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32, c_void_p)
# def no_api_callback(ntracedtotal, ntraced, ntotrace, curstage, nstages, thread_id):
#     return 1


# def _thread_func(pobj, as_pt, seed, id, no_callback):
#     pobj.run(seed, as_pt, 0, id, no_callback)
#     return copy.deepcopy(pobj.raydata), copy.copy(pobj.sunstats)
import ctypes, random
from collections.abc import Callable
from typing import Literal
import numpy as np
import pandas as pd
from datetime import datetime

from pysoltrace import api, dot_h, soltrace_json as st_json, math_utils, Point, soltrace_structs as _STS

_default_cls_arg = lambda arg, default_cls, *args: arg if arg != None else default_cls(*args)

def _format_class(prefix, 
                  pairs, 
                  underline: bool = True, 
                  buffer:    int  = 1,
                  force_max: int  = 0):
    if underline: prefix += ':\n' + '-'*len(prefix) + '\n'
    else: prefix += ':\n'

    max_key_len = force_max if force_max > 0 else max(len(str(p[0])) for p in pairs)
    fmt = lambda k, v: f'{str(k + ':'):<{max_key_len}} {v}'

    buf_str = ' ' * buffer

    return f'{prefix}{buf_str}' + f'\n{buf_str}'.join([fmt(k, v) for k, v in pairs])

# i think these underscore classes can go one layer up
# interfacing with SolTrace thru these formats: 
#  - stinput
#  - json
#  - an individial member
#  - all members -> look at what happens after sscanf portions of simdata_io.cpp and expose

######################
# Optical Face class #
######################
class _Face:
    CHAR_TO_NAME = {
        'g': st_json.ERROR_GAUSSIAN, 
        'p': st_json.ERROR_PILLBOX, 
        'f': st_json.ERROR_DIFFUSE, 
        'd': st_json.ERROR_USER_DEFINED
    }

    def __init__(self,
                    dist_type:       Literal['g', 'p', 'f', 'd'] = 'g',     #One of 'g'->Gauss 'p'->Pillbox 'd'->Diffuse # SolTrace::Data::char_to_distribution
                    refraction_real: float = 1.1,         #real component of the refraction index
                    reflectivity:    float = 0.96,         #reflectivity
                    transmissivity:  float = 0.,         #transmissivity
                    slope_error:     float = 0.95,          #RMS slope error [mrad]
                    spec_error:      float = 0.2,            #RMS specularity error [mrad]
                    userefltable:    bool  = False,             #Flag [bool] use reflectivity table
                    refltable:       list  = [],  #[[angle1,refl1],[...]]
                    usetranstable:   bool  = False,             #Flag [bool] use transmissivity table
                    transtable:      list  = [],  #[[angle1,trans1],[...]]
                ):
        self.dist_type       = dist_type
        self.refraction_real = refraction_real
        self.reflectivity    = reflectivity
        self.transmissivity  = transmissivity
        self.slope_error     = slope_error
        self.spec_error      = spec_error
        self.userefltable    = userefltable
        self.refltable       = refltable
        self.usetranstable   = usetranstable
        self.transtable      = transtable
        ## Distribution type for surface interactions. One of:
        # {'g':Gaussian, 'p':Pillbox, 'd':Diffuse }
        ## Real component of the refraction index
        ## [0..1] Surface reflectivity
        ## [0..1] Surface transmissivity
        ## [mrad] Surface RMS slope error, half-angle
        ## [mrad] Surface specularity error, half-angle
        ## Flag specifying use of user reflectivity table to modify reflectivity as a function of incidence angle
        ## [mrad,0..1] 2D list containing pairs of [angle,reflectivity] values.
        ## Flag specifying use of user transmissivity table to modify transmissivity as a function of incidence angle
        ## [mrad,0..1] 2D list containing pairs of [angle,transmissivity] values.
    
    def __repr__(self) -> str:
        pairs = [
            ('Distribution Type',        _Face.CHAR_TO_NAME[self.dist_type]),
            ('Refraction Real',          self.refraction_real),
            ('Reflectivity',             self.reflectivity),
            ('Transmissivity',           self.transmissivity),
            ('Slope Error',              self.slope_error),
            ('Specularity Error',        self.spec_error),
            ('Use Reflectivity Table',   self.userefltable), 
            ('Use Transmissivity Table', self.usetranstable),
        ]
        if self.usetranstable: pairs.append(('Transmissivity Table', self.transtable))
        if self.userefltable: pairs.append(('Reflectivity Table', self.refltable))
        return _format_class('Optical Face', pairs, False, 2, len('Use Transmissivity Table') + 1)

    def copy(self, fnew: _Face = None) -> _Face:
        new_face = _Face(self.dist_type,
                         self.refraction_real,
                         self.reflectivity,
                         self.transmissivity,
                         self.slope_error,
                         self.spec_error,
                         self.userefltable,
                         self.refltable.copy(),
                         self.usetranstable,
                         self.transtable.copy())
        if fnew != None: fnew = new_face
        else:            return new_face

################
# Optics class #
################
class _Optics:
    Face = _Face

    ########################
    # Optics class methods #
    ########################
    def __init__(self, 
                 id:    int, 
                 name:  str   = 'New Optic',
                 type: Literal[1, 2] = 2,
                 front: _Face = None,
                 back:  _Face = None):
        self.id    = id
        self.name  = name
        self.type  = type
        self.front = _default_cls_arg(front, _Face)
        self.back  = _default_cls_arg(back, _Face)

    def __repr__(self) -> str:
        prefix = f'{self.name} - Optical Properties'
        prefix += '\n' + '-'*len(prefix) + '\n'
        return prefix + f' id: {self.id}\n Front {self.front}\n Back {self.back}'

    def copy(self, onew: _Optics = None) -> _Optics:
        new_optics = _Optics(self.id,
                             self.name,
                             self.type,
                             self.front.copy(),
                             self.back.copy())
        if onew != None: onew = new_optics
        else:            return new_optics

    def Create(self, stapi: api, _do: bool = False):
        opt_set = _STS.optical_properties_set(self.name.encode(), 
                                              self.front.refraction_real,
                                              self.back.refraction_real,
                                              self.type)
        front   = _STS.optical_properties_face(self.front.transmissivity,
                                               self.front.reflectivity,
                                               self.front.slope_error,
                                               self.front.spec_error,
                                               self.front.dist_type[0].encode())
        back    = _STS.optical_properties_face(self.back.transmissivity,
                                               self.back.reflectivity,
                                               self.back.slope_error,
                                               self.back.spec_error,
                                               self.back.dist_type[0].encode())
        if _do: return stapi.data.optic.add(opt_set, front, back)
        return stapi.batch.data.optic.add(opt_set, front, back)
    
        # call = stapi.generate_api_call(dot_h.st_api_call.CALL_ST_ADD_OPTICAL_PROPERIES_SET,
        #                                ctypes.pointer(opt_set),
        #                                ctypes.pointer(front),
        #                                ctypes.pointer(back),
        #                                ctypes.pointer(ctypes.c_uint64()))
        # if _do: return stapi.batch([call])
        # return call

#############
# Sun class #
#############
class _Sun:
    CHAR_TO_NAME = {
        'g': st_json.SUN_GAUSSIAN, 
        'p': st_json.SUN_PILLBOX, 
        'd': st_json.SUN_USER_DEFINED
    }

    def __init__(self,
                 point_source = False,
                 shape: Literal['g', 'p', 'b', 'd'] = 'g',
                 sigma: float = 4.65,
                 position: Point = None,
                 user_intensity_table: list = []):
        self.point_source = point_source
        self.shape        = shape
        self.sigma        = sigma
        self.position     = _default_cls_arg(position, Point, 0, 0, 100)
        self.user_intensity_table = user_intensity_table
        #  ## Flag indicating whether the sun is modeled as a point source at a finite distance.
        # ## Sun shape model. One of: {'p':Pillbox, 'g':Gaussian, 'd':data table, 'f':gray diffuse}
        # ## [mrad] Half-width or std. dev. of the error distribution
        # ## Location of the sun/sun vector in global coordinates
        # ## [mrad, 0..1] 2D list containing pairs of
        # # angle deviation from sun vector and irradiation intensity.
        # # A typical table will have angles spanning 0->~5mrad, and inten-
        # # sities starting at 1 and decreasing to zero. The table must
        # # contain at least 2 entries.

    def __repr__(self) -> str:
        pairs = [
            ('Point Source', self.point_source),
            ('Shape', _Sun.CHAR_TO_NAME[self.shape]),
            ('Sigma', self.sigma),
            ('Position', self.position),
        ]
        if len(self.user_intensity_table): pairs.append(('Intensity Table', self.user_intensity_table))
        return _format_class('Sun', pairs)

    def copy(self, snew: _Sun = None) -> _Sun:
        new_sun = _Sun(self.point_source,
                       self.shape,
                       self.sigma,
                       self.position.copy(),
                       self.user_intensity_table.copy())
        if snew != None: snew = new_sun
        else:            return new_sun

    def Create(self, stapi: api, _do: bool = False):
        npoints = len(self.user_intensity_table)
        args = _STS.sun(npoints,
                        self.position.x,
                        self.position.y,
                        self.position.z,
                        self.sigma,
                        self.shape[0].encode())

        _angle     = []
        _intensity = []
        if npoints > 2 and self.shape.lower()[0] == 'd':
            cleaned = list(zip(*self.user_intensity_table))
            _angle     = list(cleaned[0])
            _intensity = list(cleaned[1])

        if _do: return stapi.data.sun.add(args, _angle, _intensity)
        return stapi.batch.data.sun.add(args, _angle, _intensity)

        # call = stapi.generate_api_call(dot_h.st_api_call.CALL_ST_ADD_SUN,
        #                                ctypes.pointer(args),
        #                                _angle,
        #                                _intensity)
        # if _do: return stapi.batch([call])
        # return call

    def calc_sun_vector(self): pass # TODO: expose SolTrace::Data::SolarPositionCalculator

#################
# Element class #
#################
class _Element:
    PARAM_RANGE = range(8)
    APERTURE_CHAR_TO_NAME = {
        'c': st_json.CIRCLE,
        'h': st_json.HEXAGON,
        't': st_json.EQUILATERAL_TRIANGLE,
        'r': st_json.RECTANGLE,
        'a': st_json.ANNULUS,
        'l': st_json.SINGLE_AXIS_CURVATURE_SECTION,
        'i': st_json.IRREGULAR_TRIANGLE,
        'q': st_json.IRREGULAR_QUADRILATERAL,
    }

    SURFACE_CHAR_TO_NAME = {
        's': st_json.SPHERE,
        'p': st_json.PARABOLA,
        'o': st_json.HYPER,
        'g': st_json.GENERAL_SPENCER_MURTY,
        'f': st_json.FLAT,
        'c': st_json.CONE,
        't': st_json.CYLINDER,
        'd': st_json.TORUS
    }

    def __init__(self, 
                 parent_stage:    _Stage, 
                 element_id:      int,
                 enabled:         bool = True,
                 virtual:         bool = False,
                 position:        Point = None,
                 aim:             Point = None,
                 zrot:            float = 0.,
                 aperture:        Literal['c', 'h', 't', 'r', 'a', 'l', 'i', 'q'] = 'r',
                 aperture_params: list = None,
                 surface:         Literal['s', 'p', 'o', 'g', 'f', 'c', 't', 'd'] = 'f',
                 surface_params:  list = None,
                 surface_file:    str = None,
                 interaction:     Literal[1, 2] = 1, #1=refract, 2=reflect
                 optic:           _Optics = None):
        self.stage_id        = parent_stage.id
        self.id              = element_id
        self.enabled         = enabled
        self.virtual         = virtual
        self.position        = _default_cls_arg(position, Point)
        self.aim             = _default_cls_arg(aim, Point, 0, 0, 1)
        self.zrot            = zrot
        self.aperture        = aperture
        self.aperture_params = aperture_params if aperture_params != None else [0. for _ in _Element.PARAM_RANGE]
        self.surface         = surface
        self.surface_params  = surface_params if surface_params != None else [0. for _ in _Element.PARAM_RANGE]
        self.surface_file    = surface_file
        self.interaction     = interaction #1=refract, 2=reflect
        self.optic           = optic
        # ## Identifying integer associated with the containing stage
        # ## Identifying integer associated with element
        # ## Flag indicating whether the element is included in the model
        # ## Element location in stage coordinates
        # ## Element coordinate system aim point in stage coordinates
        # self.aim.z = 1
        # ## [deg] Rotation of coordinate system around z-axis
        # ## Charater indicating aperture type.
        # ## Up to 8 coefficients defining aperture -- values depend on selection for 'aperture'
        # ## Character indicating surface type.
        # ## Up to 8 coefficients defining surface -- values depend on selection for 'surface'
        # ## Name for surface file, if using compatible type.
        # ## Flag indicating optical interaction type. {1:refraction, 2:reflection}
        # ## Reference to *Optics* instance associated with this element

    def __repr__(self) -> str:
        pairs = [
            ('Stage ID', self.stage_id),
            ('ID', self.id),
            ('Enabled', self.enabled),
            ('Position', self.position),
            ('Aim Point', self.aim),
            ('Z Rotation', self.zrot),
            ('Interaction Type', st_json.REFRACTION if self.interaction == 1 else st_json.REFLECTION),
            ('Aperture Type', _Element.APERTURE_CHAR_TO_NAME[self.aperture]),
            ('Aperture Parameters', self.surface_params),
            ('Surface Type', _Element.SURFACE_CHAR_TO_NAME[self.surface]),
            ('Surface Parameters', self.surface_params),
        ]
        if self.surface_file != None: pairs.append(('Surface File', self.surface_file))
        return '\n'.join([_format_class('Element', pairs, False, 2, len('Aperture Parameters') + 1), str(self.optic)])

    def copy(self, enew: _Element = None) -> _Element:
        new_el = _Element(self.stage_id,
                          self.id,
                          self.enabled,
                          self.virtual,
                          self.position.copy(),
                          self.aim.copy(),
                          self.zrot,
                          self.aperture,
                          self.aperture_params.copy(),
                          self.surface,
                          self.surface_params.copy(),
                          self.surface_file,
                          self.interaction,
                          self.optic)
        if enew != None: enew = new_el
        else:            return new_el

    def Create(self, stapi: api, unstager: Callable | None = None, _do: bool = False):
        # TODO: unstage element position/aim
        el_args = _STS.element(unstager(self.position) if unstager else self.position,
                               unstager(self.aim) if unstager else self.aim,
                               self.zrot,
                               self.enabled,
                               self.virtual)
        geom_args = _STS.geometry(self.aperture_params,
                                  self.surface_params,
                                  self.aperture.encode(),
                                  self.surface.encode())

        if _do: return stapi.data.element.add(el_args, self.optic.id, geom_args)
        return stapi.batch.data.element.add(el_args, self.optic.id, geom_args)
        # _a_params = (ctypes.c_double * 8)(*self.aperture_params)
        # _s_params = (ctypes.c_double * 8)(*self.surface_params)
        # call = stapi.generate_api_call(dot_h.st_api_call.CALL_ST_ADD_ELEMENT,
        #                                ctypes.pointer(el_args),
        #                                self.optic.id,
        #                                _a_params,
        #                                _s_params,
        #                                ctypes.pointer(ctypes.c_uint64()))
        # if _do: return stapi.batch([call])
        # return call

    #####################
    # Surface Functions #
    #####################
    
    def surface_spherical(self, radius: float):
        self.surface = 's'
        self.surface_params[0] = 1. / radius
        self.surface_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def surface_parabolic(self, focal_len_x: float, focal_len_y: float):
        self.surface = 'p'
        # skipping conversion, 1. / (2. * focal_len), as it's 
        # calculated in SolTrace::Data::make_surface_from_type
        self.surface_params[0] = focal_len_x
        self.surface_params[1] = focal_len_y
        self.aperture_params[2:8] = [0., 0., 0., 0., 0., 0.]
        return True
    
    def surface_flat(self):
        self.surface = 'f'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        return True
    
    def surface_hypellip(self, vertex_curv: float, kappa: float):
        self.surface = 'o'
        self.surface_params[0] = vertex_curv
        self.surface_params[1] = kappa
        self.aperture_params[2:8] = [0., 0., 0., 0., 0., 0.]
        return True

    def surface_conical(self, theta: float):
        self.surface = 'c'
        self.surface_params[0] = theta
        self.aperture_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def surface_cylindrical(self, radius: float):
        self.surface = 't'
        self.surface_params[0] = 1./radius
        self.aperture_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def surface_toroid(self, rad_annulus: float, rad_ring: float):
        self.surface = 'd'
        self.surface_params[0] = rad_annulus
        self.surface_params[1] = rad_ring
        self.aperture_params[2:8] = [0., 0., 0., 0., 0., 0.]
        return True

    def surface_zernicke(self, file_path: str):
        self.surface = 'm'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        self.surface_file = file_path
        return True

    def surface_polynomialrev(self, file_path: str):
        self.surface = 'r'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        self.surface_file = file_path
        return True

    def surface_cubicspline(self, file_path: str):
        self.surface = 'i'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        self.surface_file = file_path
        return True

    def surface_finiteelement(self, file_path: str):
        self.surface = 'e'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        self.surface_file = file_path
        return True

    def surface_vshot(self, file_path: str):
        self.surface = 'v'
        self.surface_params = [0. for _ in _Element.PARAM_RANGE]
        self.surface_file = file_path
        return True
    
    ######################
    # Aperture Functions #
    ######################
    
    def aperture_circle(self, diameter: float):
        self.aperture = 'c'
        self.aperture_params[0] = diameter
        self.aperture_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def aperture_hexagon(self, diameter: float):
        self.aperture = 'h'
        self.aperture_params[0] = diameter
        self.aperture_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def aperture_triangle(self, diameter: float):
        self.aperture = 't'
        self.aperture_params[0] = diameter
        self.aperture_params[1:8] = [0., 0., 0., 0., 0., 0., 0.]
        return True

    def aperture_rectangle(self, length_x: float, length_y: float):
        self.aperture = 'r'
        self.aperture_params[0] = length_x
        self.aperture_params[1] = length_y
        self.aperture_params[2:8] = [0., 0., 0., 0., 0., 0.]
        return True
    
    def aperture_annulus(self, r_inner: float, r_outer: float, theta: float):
        self.aperture = 'a'
        self.aperture_params[0] = r_inner
        self.aperture_params[1] = r_outer
        self.aperture_params[2] = theta
        self.aperture_params[3:8] = [0., 0., 0., 0., 0.]
        return True

    def aperture_singleax_curve(self, x1: float, x2: float, L: float):
        self.aperture = 'l'
        self.aperture_params[0] = x1
        self.aperture_params[1] = x2
        self.aperture_params[2] = L
        self.aperture_params[3:8] = [0., 0., 0., 0., 0.]
        return True

    def aperture_irr_triangle(self, x1: float, y1: float, x2: float, y2: float, x3: float, y3: float):
        self.aperture = 'i'
        self.aperture_params[0] = x1
        self.aperture_params[1] = y1
        self.aperture_params[2] = x2
        self.aperture_params[3] = y2
        self.aperture_params[4] = x3
        self.aperture_params[5] = y3
        self.aperture_params[6:8] = [0., 0.]
        return True

    def aperture_quadrilateral(self, x1: float, y1: float, x2: float, y2: float, x3: float, y3: float, x4: float, y4: float):
        self.aperture = 'q'
        self.aperture_params[0] = x1
        self.aperture_params[1] = y1
        self.aperture_params[2] = x2
        self.aperture_params[3] = y2
        self.aperture_params[4] = x3
        self.aperture_params[5] = y3
        self.aperture_params[6] = x4
        self.aperture_params[7] = y4
        return True

###############
# Stage class #
###############
class _Stage:
    """Depreceating stages from SolTrace.
    
    This class exists to maintain structure of legacy 
    Python wrapper and provide unstaging utilities.
    """
    Element = _Element

    def __init__(self, 
                 id:              int,
                 position:        Point = None,
                 aim:             Point = None,
                 zrot:            float = 0.,
                 is_virtual:      bool  = False,
                 is_multihit:     bool  = True,
                 is_tracethrough: bool  = False,
                 name:            str   = None,
                 elements:        list[_Element] = None):
        self.id              = id
        self.position        = _default_cls_arg(position, Point)
        self.aim             = _default_cls_arg(aim, Point, 0, 0, 1)
        self.zrot            = zrot
        self.is_virtual      = is_virtual
        self.is_multihit     = is_multihit
        self.is_tracethrough = is_tracethrough
        self.name            = name if name != None else "stage_{:d}".format(id)
        self.elements        = elements if elements != None else []
        # ## Identifying integer associated with the stage
        # ## Stage location in global coordinates
        # ## Coordinate system aim point in global coordinates
        # self.aim.z = 1
        # ## [deg] Rotation of coordinate system around z-axis
        # ## Flag indicating virtual stage
        # ## Flag indicating that rays can have multiple interactions within a single stage.
        # ## Flag indicating the stage is in trace-through mode
        # ## Descriptive name for this stage
        # ## list of all elements in the stage

    def __repr__(self):
        return '\n'.join([str(e) for e in self.elements])

    def copy(self, snew: _Stage = None) -> _Stage:
        new_stage = _Stage(self.id,
                           self.position.copy(),
                           self.aim.copy(),
                           self.zrot,
                           self.is_virtual,
                           self.is_multihit,
                           self.is_tracethrough,
                           self.name,
                           [e.copy() for e in self.elements])
        if snew != None: snew = new_stage
        else:            return new_stage

    def __get_unstager(self) -> Callable | None:
        return math_utils.get_unstager(self.position, self.aim, 0)

    def Create(self, stapi: api, _do: bool = False):
        unstager = self.__get_unstager()
        # calls = [el.Create(stapi, unstager) for el in self.elements]
        # if _do: return stapi.batch(calls)
        # return calls
        return [el.Create(stapi, unstager, _do) for el in self.elements]

    def add_element(self, enew: _Element = None) -> _Element:
        enew    = _default_cls_arg(enew, _Element, self, None)
        enew.id = len(self.elements)
        self.elements.append(enew)
        return enew

class legacy:
    Optics = _Optics
    Sun    = _Sun
    Stage  = _Stage

    ###################################
    # Legacy PySolTrace class methods #
    ###################################
    def __init__(self,
                 sun:               _Sun = None,
                 optics:            list[_Optics] = [],
                 stages:            list[_Stage]  = [],
                 num_ray_hits:      int   = int(1e5),
                 max_rays_traced:   int   = int(1e7),
                 is_sunshape:       bool  = True,
                 is_surface_errors: bool  = True,
                 dni:               float = 1000.,
                 testing:           bool  = False):
        self.sun = sun
        self.optics = optics
        self.stages = stages
        self.num_ray_hits = num_ray_hits
        self.max_rays_traced = max_rays_traced
        self.is_sunshape = is_sunshape
        self.is_surface_errors = is_surface_errors
        self.dni = dni  #w/m^2
        self.testing = testing
        self.raydata = None
        self.sunstats = None
        self.powerperray = None
        self.__op_registry = st_json.OpticalPropertyRegistry()

    def copy(self, stnew: legacy = None) -> legacy:
        new_st = legacy(self.sun,
                        [o.copy() for o in self.optics],
                        [s.copy() for s in self.stages],
                        self.num_ray_hits,
                        self.max_rays_traced,
                        self.is_sunshape,
                        self.is_surface_errors,
                        self.dni)
        if self.raydata != None:     new_st.raydata     = self.raydata.copy()
        if self.sunstats != None:    new_st.sunstats    = self.sunstats.copy()
        if self.powerperray != None: new_st.powerperray = self.powerperray
        if stnew != None: stnew = new_st
        else:             return new_st

    def Create(self, as_power_tower: bool, _do: bool = False):
        calls = [
            self.stapi.batch.legacy.sim_params(int(self.num_ray_hits),
                                               int(self.max_rays_traced),
                                               as_power_tower),
            self.stapi.batch.parameters.errors(self.is_sunshape,
                                               self.is_surface_errors),
            self.sun.Create(self.stapi),
            *[op.Create(self.stapi) for op in self.optics]
        ]
        for s in self.stages: calls.extend(s.Create(self.stapi))

        if _do: self.stapi.batch()
        return calls
    
    def add_optic(self, 
                  optic_name: str,
                  type: Literal[1, 2] = 2,
                  front: _Optics = None,
                  back: _Optics = None) -> _Optics:
        new_optics = _Optics(len(self.optics),
                             optic_name,
                             type,
                             front,
                             back)
        self.optics.append(new_optics)
        return new_optics

    def clear_optics(self):
        self.optics = []

    def add_sun(self):
        self.sun = _Sun()
        return self.sun

    def add_stage(self, snew: _Stage = None) -> _Stage:
        snew    = _default_cls_arg(snew, _Stage, len(self.stages))
        snew.id = len(self.stages)
        self.stages.append(snew)
        return snew

    def clear_stages(self):
        self.stages = []
    

    def get_ray_dataframe(self):
        """
        Get a pandas dataframe with all of the ray data from the simulation.

        Returns
        ----------
        Pandas.DataFrame
            with columns:
            loc_x   | Ray hit location, x-coordinate
            loc_y   | Ray hit location, y-coordinate
            loc_z   | Ray hit location, z-coordinate
            cos_x   | Ray directional vector, x-component
            cos_y   | Ray directional vector, y-component
            cos_z   | Ray directional vector, z-component
            element | Element associated with ray hit
            stage   | Stage associated with ray hit
            number  | Ray number
        """
        if not self.stapi: raise ValueError("No SolTrace API object (STAPIv2) created. Call PySolTrace.run to initialize the object.")

        if self.raydata: return self.raydata

        results_data = self.stapi.result.get()

        self.raydata = pd.DataFrame({
            "loc_x":   results_data['loc_x'],
            "loc_y":   results_data['loc_y'],
            "loc_z":   results_data['loc_z'],
            "cos_x":   results_data['cos_x'],
            "cos_y":   results_data['cos_y'],
            "cos_z":   results_data['cos_z'],
            "element": results_data['element_map'],
            "stage":   results_data['stage_map'],
            "number":  results_data['ray_numbers']
        })

        return self.raydata

    def get_sun_stats(self):
        """
        Get information on the sun box.

        Returns
        ----------
        dict
            Keys in the return dictionary are:
            'width' --> width of sun box (not set if using OptixRunner)
            'height' --> height of sun box (not set if using OptixRunner)
            'area' --> area of sun box
            'nsunrays' --> Number of sun rays simulated
        """
        if not self.stapi: raise ValueError("No SolTrace API object (STAPIv2) created. Call PySolTrace.run to initialize the object.")
        
        if self.sunstats: return self.sunstats

        self.sunstats = {
            item[0]: item[1] 
            for item in 
            zip(['width', 'height', 'area', 'nsunrays'],
                self.stapi.result.sun_stats())
        }
        self.sunstats['position'] = self.sun.position

        return self.sunstats

    def run(self,
            seed:           int  = -1,
            as_power_tower: bool = False,
            nthread:        int  = 1,
            thread_id:      int  = 0,
            no_callback:    bool = False,
            runner_type:    Literal[0, 1, 2] = 0):
        # Parameters
        # ----------
        # seed : int
        #     Seed for random number generator. [-1] for random seed. Seeding happens
        #     differently for single vs multi-thread modes.
        #         * If nthreads == 1 and seed < 0: a random int is chosen as the seed value.
        #         * If nthreads > 1 and seed < 0: a random int is chosen for the first
        #           thread seed value. Other threads i=1..(nthreads-1) are assigned
        #           (first value) + i*123.
        # as_power_tower : bool
        #     Flag indicating simulation should be processed as power
        #     tower / central receiver type, with corresponding efficiency adjustments.
        # nthread : int
        #     Number of threads to execute. Will be limited by the method to the number
        #     available on the machine.
        #         * If nthreads > 1, the function will call recursively while setting
        #           nthreads=0 for each thread spawned.
        #         * If nthreads == 1, the function will run in single-thread mode. Seed
        #           values are checked.
        #         * If nthreads == 0, the function will run in single-thread mode. Seed
        #           values are not checked and should be handled prior to calling in
        #           this mode.
        # thread_id : int
        #     Argument used by the multi-threading call. Do not manually specify this value.

        # pdll = self.__load_dll()
        self.stapi = api(testing = self.testing)

        if seed<0:
            runseed = random.randint(1,int(1e9))
        else:
            runseed = seed

        if nthread in [0,1]:
            sim_data_calls = self.Create(as_power_tower)
            self.stapi.batch.runner.setup(runner_type, 1, [runseed])
            self.stapi.batch.runner.run()
            self.stapi.batch.runner.report()
            # _seeds = (ctypes.c_uint * 1)(*[runseed])
            # set_up_call = self.stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_SETUP, runner_type, 1, _seeds, 1)
            # run_call = self.stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_RUN_V2)
            # report_call = self.stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_REPORT, 0)

            self.stapi.batch(True)

            self.get_ray_dataframe()
            self.get_sun_stats()
            self.powerperray = self.sunstats['area'] * self.dni / self.sunstats['nsunrays']

            # print(raydata)

        #     """self.Create(pdll, p_data)"""
        #     # self.sun.Create(pdll, p_data)
        #     # for opt in self.optics:
        #     #     opt.Create(pdll, p_data)
        #     # for stage in self.stages:
        #     #     stage.Create(pdll, p_data)

        #     pdll.st_sim_errors.restype = c_int
        #     pdll.st_sim_errors(c_void_p(p_data), c_int(1 if self.is_sunshape else 0), c_int(1 if self.is_surface_errors else 0))

        #     pdll.st_sim_params.restype = c_int
        #     pdll.st_sim_params(c_void_p(p_data), c_int(int(self.num_ray_hits)), c_int(int(self.max_rays_traced)), c_int(as_power_tower))

        #     if thread_id == 0:
        #         tstart = time.time()

        #     pdll.st_sim_run.restype = c_int
        #     if no_callback:
        #         res = pdll.st_sim_run( c_void_p(p_data), c_uint16(runseed), no_api_callback, thread_id)
        #     else:
        #         res = pdll.st_sim_run( c_void_p(p_data), c_uint16(runseed), api_callback, thread_id)
        #         if thread_id == 0:
        #             print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

        #     # Collect simulation output, including raw ray data and sunbox stats
        #     self.raydata = self.__get_ray_dataframe(pdll,p_data)
        #     self.sunstats = self.__get_sun_stats(pdll, p_data)
        #     # Compute and save power per ray
        #     self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni

        #     pdll.st_free_context.restype = c_bool
        #     pdll.st_free_context(c_void_p(p_data))

        #     return res
        # else:
        #     seeds = [seed + i*123 for i in range(nthread)]

        #     P = [[self.copy(), as_power_tower, seeds[i], i+1, no_callback] for i in range(nthread)]

        #     # modify the number of rays to match the required totals
        #     nrpt = int(float(self.num_ray_hits)/float(nthread))
        #     mrpt = int(float(self.max_rays_traced)/float(nthread))

        #     for p in P:
        #         p[0].num_ray_hits = nrpt
        #         p[0].max_rays_traced = mrpt
        #         if p == P[0]:
        #             p[0].num_ray_hits += int(float(self.num_ray_hits) % float(nthread))
        #             p[0].max_rays_traced += int(float(self.max_rays_traced) % float(nthread))

        #     pool = multiprocessing.Pool(nthread)
        #     if not no_callback:
        #         print("Launching {:d} threads...".format(nthread))
        #     tstart = time.time()
        #     res = pool.starmap_async(_thread_func, P)
        #     pool.close()
        #     pool.join()
        #     if not no_callback:
        #         print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

        #     # Modify the ray number for threads 2+ to avoid duplication
        #     try:
        #         dfs = [r[0] for r in res.get()]
        #         rstart = int(dfs[0].iloc[-1].number)
        #     except:
        #         print("Unknown error caused the simulation to fail. Try re-running.")
        #         return

        #     if len(dfs)>1:
        #         for d in dfs[1:]:
        #             d.number = d.number+rstart
        #             rstart = d.number.iloc[-1]

        #     self.raydata = pd.concat(dfs)
        #     self.raydata.reset_index(inplace=True)
        #     self.sunstats = res.get()[0][1]  #take the first thread result
        #     # add up all the sunrays from all threads
        #     srct = 0
        #     for r in res.get():
        #         srct += r[1]['nsunrays']
        #     self.sunstats['nsunrays'] = srct

        #     # Compute and save power per ray [W]
        #     self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni

        #     return 1
        pass

    def write_soltrace_input_file(self, path: str):
        """
        Write a SolTrace input file (.stinput) based on the currently created API objects. This file is written
        using the objects and data in the PySolTrace instance, not necessarily on what has been created in the
        coretrace 'context' data space. The 'context' may not match the PySolTrace instance if not all 'Create()'
        methods have been called.

        Parameters
        ==========
        path : str
            Path and file name to be used to write the resulting .stinput file

        Returns
        =======
        None
        """
        with open(path, 'w') as fout:

            # Header
            dt = datetime.now()
            fout.write( "# SOLTRACE VERSION 2026.07.15 INPUT FILE -- GENERATED BY stapi_v2 v{:s} on {:02d}/{:02d}/{:04d} at {:02d}:{:02d}:{:02d}\n".format(
                "v000",
                dt.day,
                dt.month,
                dt.year,
                dt.hour,
                dt.minute,
                dt.second))

            #------------------- Sun shape
            fout.write( "SUN\tPTSRC\t{:d}\tSHAPE\t{:1s}\tSIGMA\t{:f}\tHALFWIDTH\t{:f}\n".format( 0, self.sun.shape, self.sun.sigma, self.sun.sigma) )
            fout.write( "XYZ\t{:f}\t{:f}\t{:f}\tUSELDH\t{:d}\tLDH\t{:f}\t{:f}\t{:f}\n".format(
                self.sun.position.x,
                self.sun.position.y,
                self.sun.position.z,
                0,
                0.,
                0.,
                0.))
            if( self.sun.shape == 'd' ):
                np = len(self.sun.user_intensity_table)
                fout.write( "USER SHAPE DATA\t{:d}\n".format( np) )
                for i in range(np):
                    fout.write( "{:f}\t{:f}\n".format( *self.sun.user_intensity_table[i] ) )
            else:
                fout.write( "USER SHAPE DATA\t{:d}\n".format( 0) )

            #------------------- Optics list
            fout.write( "OPTICS LIST COUNT\t{:d}\n".format( len(self.optics) ) )
            for optics in self.optics:
                fout.write( "OPTICAL PAIR\t{:s}\n".format( optics.name ) )
                for opt in [optics.front, optics.back]:
                    fout.write( "OPTICAL\t{:1s}\t{:d}\t{:d}\t{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\n".format(
                        opt.dist_type,
                        3,
                        1,
                        4,
                        opt.reflectivity,
                        opt.transmissivity,
                        opt.slope_error,
                        opt.spec_error ,
                        opt.refraction_real,
                        0.,
                        0.,
                        0.,
                        0.,
                        0. ))

            #------------------- loop through Stage list
            fout.write( "STAGE LIST COUNT\t{:d}\n".format( len(self.stages)) )
            for stage in self.stages:
                fout.write( "STAGE\tXYZ\t{:f}\t{:f}\t{:f}\tAIM\t{:f}\t{:f}\t{:f}\tZROT\t{:f}\tVIRTUAL\t{:d}\tMULTIHIT\t{:d}\tELEMENTS\t{:d}\tTRACETHROUGH\t{:d}\n".format(
                    stage.position.x,
                    stage.position.y,
                    stage.position.z,
                    stage.aim.x,
                    stage.aim.y,
                    stage.aim.z,
                    stage.zrot,
                    int(stage.is_virtual),
                    int(stage.is_multihit),
                    len(stage.elements),
                    int(stage.is_tracethrough)))

                fout.write( "{:s}\n".format( stage.name ) )

                #------------------- loop through element list
                for el in stage.elements:
                    # format string
                    fmt = '\t'.join(['{:d}'] + ['{:f}']*7 + (['{:1s}'] + ['{:f}']*8)*2 + ['{:s}', '{:s}', '{:d}', '{:s}\n'])
                    fout.write(fmt.format(
                        int(el.enabled),
                        el.position.x, el.position.y, el.position.z,  #origin
                        el.aim.x, el.aim.y, el.aim.z,  #aim
                        el.zrot,
                        el.aperture,
                        # Ap_A, Ap_B, Ap_C, Ap_D, Ap_E, Ap_F, Ap_G, Ap_H,
                        *el.aperture_params,
                        el.surface,
                        # Su_A, Su_B, Su_C, Su_D, Su_E, Su_F, Su_G, Su_H,
                        *el.surface_params,
                        "" if not el.surface_file else el.surface_file, #Surface geometry file
                        el.optic.name, el.interaction,
                        "" ) )
        return


    ####################################
    # utility transform/math functions #
    ####################################
    @staticmethod
    def util_calc_euler_angles(origin:   Point | list | np.ndarray,
                               aimpoint: Point | list | np.ndarray,
                               zrot:     float):
        return math_utils.euler_angles(origin, aimpoint, zrot)
    
    @staticmethod
    def util_transform_to_local(posref:    np.ndarray,
                                cosref:    np.ndarray,
                                origin:    np.ndarray,
                                rreftoloc: np.ndarray):
        return math_utils.transform_to_local(posref, cosref, origin, rreftoloc)
    
    @staticmethod
    def util_transform_to_ref(posloc:    np.ndarray,
                              cosloc:    np.ndarray,
                              origin:    np.ndarray,
                              rloctoref: np.ndarray):
        return math_utils.util_transform_to_ref(posloc, cosloc, origin, rloctoref)
    
    @staticmethod
    def util_matrix_vector_mult(m: np.ndarray, v: np.ndarray):
        return math_utils.matrix_vector_mult(m, v)
    
    @staticmethod
    def util_calc_transforms(euler: Point | list): 
        return math_utils.euler_transforms(euler)
    
    @staticmethod
    def util_matrix_transpose(m: np.ndarray): 
        return m.T
    
    @staticmethod
    def util_rotation_arbitrary(theta: float,
                                axis:  Point | list,
                                axloc: Point | list,
                                pt:    Point | list):
        return math_utils.arbitrary_rotation(theta, axis, axloc, pt)

    @staticmethod
    def util_calc_unitvect(vect: Point | list): 
        return math_utils.unitize(vect)

    @staticmethod
    def util_calc_zrot_azel(vect: Point | list): 
        return math_utils.zrot_from_azel(vect)

# end class legacy ---------------------------------------------

# # ----------------------------------------------------------------------
class PySolTrace_v1:
    """
    see https://github.com/uwsbel/SolTrace/blob/5e08741c1da9e6e68cef5232d9da401586674618/app/deploy/api/pysoltrace.py for docstrings

    maintainance checklist: (x: done, T: have TODO, empty: not done)
    [x] - PySolTrace.Optics.Face.__init__
    [x] - PySolTrace.Optics.Face.copy
    [x] - PySolTrace.Optics.__init__
    [x] - PySolTrace.Optics.copy
    [x] - PySolTrace.Optics.Create
    [x] - PySolTrace.Sun.__init__
    [x] - PySolTrace.Sun.copy
    [x] - PySolTrace.Sun.Create
    [T] - PySolTrace.Sun.calc_sun_vector
    [x] - PySolTrace.Stage.Element.__init__
    [x] - PySolTrace.Stage.Element.copy
    [T] - PySolTrace.Stage.Element.Create
    [x] - PySolTrace.Stage.Element.surface_spherical
    [x] - PySolTrace.Stage.Element.surface_parabolic
    [x] - PySolTrace.Stage.Element.surface_flat
    [x] - PySolTrace.Stage.Element.surface_hypellip
    [x] - PySolTrace.Stage.Element.surface_conical
    [x] - PySolTrace.Stage.Element.surface_cylindrical
    [x] - PySolTrace.Stage.Element.surface_toroid
    [x] - PySolTrace.Stage.Element.surface_zernicke
    [x] - PySolTrace.Stage.Element.surface_polynomialrev
    [x] - PySolTrace.Stage.Element.surface_cubicspline
    [x] - PySolTrace.Stage.Element.surface_finiteelement
    [x] - PySolTrace.Stage.Element.surface_vshot
    [x] - PySolTrace.Stage.Element.aperture_circle
    [x] - PySolTrace.Stage.Element.aperture_hexagon
    [x] - PySolTrace.Stage.Element.aperture_triangle
    [x] - PySolTrace.Stage.Element.aperture_rectangle
    [x] - PySolTrace.Stage.Element.aperture_annulus
    [x] - PySolTrace.Stage.Element.aperture_singleax_curve
    [x] - PySolTrace.Stage.Element.aperture_irr_triangle
    [x] - PySolTrace.Stage.Element.aperture_quadrilateral
    [x] - PySolTrace.Stage.__init__
    [x] - PySolTrace.Stage.copy
    [ ] - PySolTrace.Stage.Create
    [x] - PySolTrace.Stage.add_element
    [x] - PySolTrace.__init__
    [x] - PySolTrace.copy
    [ ] - PySolTrace.Create
    [x] - PySolTrace.add_optic
    [ ] - PySolTrace.delete_optic
    [x] - PySolTrace.add_sun
    [x] - PySolTrace.add_stage
    [ ] - PySolTrace.delete_stage
    [ ] - PySolTrace.__load_dll
    [ ] - PySolTrace.run
    [x] - PySolTrace.__get_num_intersections
    [x] - PySolTrace.__get_sun_stats     -> PySolTrace.get_sun_stats
    [x] - PySolTrace.__get_ray_dataframe -> PySolTrace.get_ray_dataframe
    [ ] - PySolTrace.plot_trace
    [ ] - PySolTrace.plot_flux
    [ ] - PySolTrace.bin_rays
    [x] - PySolTrace.util_calc_euler_angles
    [x] - PySolTrace.util_transform_to_local
    [x] - PySolTrace.util_transform_to_ref
    [x] - PySolTrace.util_matrix_vector_mult
    [x] - PySolTrace.util_calc_transforms
    [x] - PySolTrace.util_matrix_transpose
    [x] - PySolTrace.util_rotation_arbitrary
    [x] - PySolTrace.util_calc_unitvect
    [x] - PySolTrace.util_calc_zrot_azel
    [x] - PySolTrace.write_soltrace_input_file
    """

    class Optics:
        class Face:
            def __init__(self):
                # ## Distribution type for surface interactions. One of:
                # # {'g':Gaussian, 'p':Pillbox, 'd':Diffuse }
                # self.dist_type = 'g'     #One of 'g'->Gauss 'p'->Pillbox 'd'->Diffuse
                # ## Real component of the refraction index
                # self.refraction_real = 1.1         #real component of the refraction index
                # ## [0..1] Surface reflectivity
                # self.reflectivity = 0.96         #reflectivity
                # ## [0..1] Surface transmissivity
                # self.transmissivity = 0.         #transmissivity
                # ## [mrad] Surface RMS slope error, half-angle
                # self.slope_error = 0.95          #RMS slope error [mrad]
                # ## [mrad] Surface specularity error, half-angle
                # self.spec_error = 0.2            #RMS specularity error [mrad]
                # ## Flag specifying use of user reflectivity table to modify reflectivity as a function of incidence angle
                # self.userefltable = False             #Flag [bool] use reflectivity table
                # ## [mrad,0..1] 2D list containing pairs of [angle,reflectivity] values.
                # self.refltable = []  #[[angle1,refl1],[...]]
                # ## Flag specifying use of user transmissivity table to modify transmissivity as a function of incidence angle
                # self.usetranstable = False             #Flag [bool] use transmissivity table
                # ## [mrad,0..1] 2D list containing pairs of [angle,transmissivity] values.
                # self.transtable = []  #[[angle1,trans1],[...]]
                pass

            def copy(self, fnew):
                # c = self.__dict__.copy()
                # for attr in self.__dict__.keys():
                #     fnew.__setattr__(attr, copy.deepcopy(c[attr]))
                # return
                pass

        # -------- methods of the Optics class -----------------------------------------
        def __init__(self, id : int):
            # ## Unique name for the optical property set
            # self.name = "new optic"
            # ## Identifying integer associated with the property set
            # self.id = id

            # ## properties associated with the front of the optical surface
            # self.front = PySolTrace.Optics.Face()
            # ## properties associated with the back of the optical surface
            # self.back = PySolTrace.Optics.Face()
            pass

        def copy(self, onew):
            # """
            # Deep copy of the current Optics instance

            # Inputs
            # ---------
            # onew : Optics.Face
            #     Reference to new Optics object to which data will be copied
            # """
            # onew.name = copy.copy(self.name)
            # onew.id = copy.copy(self.id)

            # self.front.copy(onew.front)
            # self.back.copy(onew.back)

            # return
            pass

        def Create(self, pdll, p_data) -> int:
            # """
            # Create Optics instance in the SolTrace context.

            # Returns
            # ----------
            # int
            #     1 if successful, 0 otherwise
            # """
            # pdll.st_add_optic.restype = c_int
            # pdll.st_add_optic(c_void_p(p_data), c_char_p(self.name.encode()))

            # pdll.st_optic.restype = c_int

            # dummy_grating = (c_number*3)()

            # resok = True

            # # for each face -- front or back
            # for i,opt in enumerate([self.front, self.back]):

            #     user_refl_angles = (c_number*len(opt.refltable))()
            #     user_refls = (c_number*len(opt.refltable))()
            #     if len(opt.refltable) > 1:
            #         user_refl_angles[:] = list(list(zip(*opt.refltable))[0])
            #         user_refls[:] = list(list(zip(*opt.refltable))[1])

            #     user_trans_angles = (c_number*len(opt.transtable))()
            #     user_trans = (c_number*len(opt.transtable))()
            #     if len(opt.transtable) > 1:
            #         user_trans_angles[:] = list(list(zip(*opt.transtable))[0])
            #         user_trans[:] = list(list(zip(*opt.transtable))[1])

            #     # Create surface optic
            #     resok = resok and pdll.st_optic( \
            #         c_void_p(p_data),
            #         c_uint32(self.id),
            #         c_int(i+1),   #front
            #         c_wchar(opt.dist_type[0]),
            #         c_int(1), #optical surface number
            #         c_int(3), #aperture grating
            #         c_int(4), #Diffraction order
            #         c_number(opt.refraction_real),
            #         c_number(0.), #imaginary component of refraction
            #         c_number(opt.reflectivity),
            #         c_number(opt.transmissivity),
            #         pointer(dummy_grating),
            #         c_number(opt.slope_error),
            #         c_number(opt.spec_error),
            #         c_int(1 if opt.userefltable else 0),
            #         c_int(len(opt.refltable)),
            #         pointer(user_refl_angles),
            #         pointer(user_refls),
            #         c_int(1 if opt.usetranstable else 0),
            #         c_int(len(opt.transtable)),
            #         pointer(user_trans_angles),
            #         pointer(user_trans),
            #         )

            # return 1 if resok else 0
            pass
    # ========end Optics class =================================================================

    # ==========================================================================================
    class Sun:
        def __init__(self):
            # ## Flag indicating whether the sun is modeled as a point source at a finite distance.
            # self.point_source = False
            # ## Sun shape model. One of: {'p':Pillbox, 'g':Gaussian, 'd':data table, 'f':gray diffuse}
            # self.shape = 'p'
            # ## [mrad] Half-width or std. dev. of the error distribution
            # self.sigma = 4.65
            # ## Location of the sun/sun vector in global coordinates
            # self.position = Point()
            # self.position.z = 100.

            # ## [mrad, 0..1] 2D list containing pairs of
            # # angle deviation from sun vector and irradiation intensity.
            # # A typical table will have angles spanning 0->~5mrad, and inten-
            # # sities starting at 1 and decreasing to zero. The table must
            # # contain at least 2 entries.
            # self.user_intensity_table = []
            pass

        def copy(self, snew):
            # snew.position = self.position.copy()
            # c = self.__dict__.copy()
            # for attr in self.__dict__.keys():
            #     if attr in ['_pdll','_p_data','position']:
            #         continue
            #     snew.__setattr__(attr, copy.deepcopy(c[attr]))
            # return
            pass

        def Create(self, pdll, p_data):
            # """
            # Create Sun instance in the SolTrace context.

            # Returns
            # ----------
            # int
            #     1 if successful, 0 otherwise
            # """

            # pdll.st_sun.restype = c_int
            # pdll.st_sun_xyz.restype = c_int

            # pdll.st_sun(c_void_p(p_data), c_int(int(self.point_source)), c_wchar(self.shape[0]), c_number(self.sigma))
            # pdll.st_sun_xyz(c_void_p(p_data), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))

            # # If a user intensity table is provided, and the shape is specified accordingly as 'd', load the data table into context
            # if len(self.user_intensity_table) > 2 and self.shape.lower()[0] == 'd':
            #     user_angles = (c_number*len(self.user_intensity_table))()
            #     user_ints = (c_number*len(self.user_intensity_table))()
            #     user_angles[:] = list(list(zip(*self.user_intensity_table))[0])
            #     user_ints[:] = list(list(zip(*self.user_intensity_table))[1])

            #     pdll.st_sun_userdata.restype = c_int
            #     return pdll.st_sun_userdata(c_void_p(p_data), c_uint32(len(self.user_intensity_table)), pointer(user_angles), pointer(user_ints))

            # return 1
            pass

        def calc_sun_vector(self, hour, day, lat):
            # """
            # Computes the sun vector associated with a given latitude, hour, and day. The coordinate system
            # follows the left hand rule:
            #     +x = east
            #     +y = north
            #     +z = zenith

            # Note that this sun position algorithm is not especially accurate compared to NREL/Solpos or the like.
            # If a very accurate method is required, use another calculation.

            # Parameters
            # ===========
            # hour : float
            #     Hour of the day, can be fractional. 12 corresponds to solar noon.
            # day : float
            #     Day of the year. 1 is january 1st
            # lat : float
            #     [rad] Lattitude (+ north, - south)

            # Returns
            # ========
            #     numpy.array([x,y,z]) unit vector in the direction of the sun
            # """

            # Declination = numpy.arcsin(0.39795 * numpy.cos(0.01720248870643171 * (day - 173)))   #[rad]
            # HourAngle = (hour/12 - 1)*numpy.pi  #[rad]
            # cos_Declination = numpy.cos(Declination)
            # sin_Declination = numpy.sin(Declination)
            # cos_HourAngle = numpy.cos(HourAngle)
            # sin_lat = numpy.sin(lat)
            # cos_lat = numpy.cos(lat)
            # Elevation = numpy.arcsin(sin_Declination * sin_lat + cos_Declination * cos_HourAngle * cos_lat)
            # Azimuth = numpy.arccos((sin_Declination * cos_lat - cos_Declination * sin_lat * cos_HourAngle) / numpy.cos(Elevation) + 0.0000000001)
            # if (numpy.sin(HourAngle) > 0.0):
            #     Azimuth = 2*numpy.pi - Azimuth
            # x = numpy.sin(Azimuth) * numpy.cos(Elevation)
            # y = numpy.cos(Azimuth) * numpy.cos(Elevation)
            # z = numpy.sin(Elevation)

            # return Point(x,y,z)
            pass

    # ===========end of the Sun class===========================================================

    # ==========================================================================================
    class Stage:

        class Element:
            # STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file);
            def __init__(self, parent_stage, element_id : int):
                # ## Identifying integer associated with the containing stage
                # self.stage_id = parent_stage.id
                # ## Identifying integer associated with element
                # self.id = element_id
                # ## Flag indicating whether the element is included in the model
                # self.enabled = True
                # ## Element location in stage coordinates
                # self.position = Point()
                # ## Element coordinate system aim point in stage coordinates
                # self.aim = Point()
                # self.aim.z = 1
                # ## [deg] Rotation of coordinate system around z-axis
                # self.zrot = 0.
                # ## Charater indicating aperture type.
                # self.aperture = 'r'
                # ## Up to 8 coefficients defining aperture -- values depend on selection for 'aperture'
                # self.aperture_params = [0. for i in range(8)]
                # ## Character indicating surface type.
                # self.surface = 'f'
                # ## Up to 8 coefficients defining surface -- values depend on selection for 'surface'
                # self.surface_params = [0. for i in range(8)]
                # ## Name for surface file, if using compatible type.
                # self.surface_file = None
                # ## Flag indicating optical interaction type. {1:refraction, 2:reflection}
                # self.interaction = 2        #1=refract, 2=reflect
                # ## Reference to *Optics* instance associated with this element
                # self.optic = None
                pass

            def copy(self, enew):
                # """
                # Deep copy of the current Element instance

                # Inputs
                # ---------
                # enew : Stage.Element
                #     Reference to new Element object to which data will be copied
                # """
                # c = self.__dict__.copy()
                # for attr in self.__dict__.keys():
                #     if attr in ['_pdll','_p_data','optic']:
                #         continue
                #     else:
                #         enew.__setattr__(attr, copy.deepcopy(c[attr]))
                # return
                pass

            def Create(self, pdll, p_data) -> int:
                # """
                # Create Element instance in the SolTrace context.

                # Returns
                # ----------
                # int
                #     1 if successful, 0 otherwise
                # """

                # pdll.st_add_element.restype = c_int
                # pdll.st_add_element(c_void_p(p_data), c_uint32(self.stage_id))

                # pdll.st_element_enabled.restype = c_int
                # pdll.st_element_xyz.restype = c_int
                # pdll.st_element_aim.restype = c_int
                # pdll.st_element_zrot.restype = c_int
                # pdll.st_element_aperture.restype = c_int
                # pdll.st_element_aperture_params.restype = c_int
                # pdll.st_element_surface.restype = c_int
                # pdll.st_element_surface_params.restype = c_int
                # pdll.st_element_surface_file.restype = c_int
                # pdll.st_element_interaction.restype = c_int
                # pdll.st_element_optic.restype = c_int

                # aperture_params = (c_number*len(self.aperture_params))()
                # surface_params = (c_number*len(self.surface_params))()
                # aperture_params[:] = self.aperture_params
                # surface_params[:] = self.surface_params

                # pdll.st_element_enabled(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(int(self.enabled)));
                # pdll.st_element_xyz(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.position.x),  c_number(self.position.y), c_number(self.position.z));
                # pdll.st_element_aim(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.aim.x),  c_number(self.aim.y), c_number(self.aim.z));
                # pdll.st_element_zrot(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.zrot) );
                # pdll.st_element_aperture(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.aperture[0]));
                # pdll.st_element_aperture_params(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(aperture_params));
                # pdll.st_element_surface(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.surface[0]));
                # pdll.st_element_surface_params(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(surface_params));
                # if self.surface_file:
                #     pdll.st_element_surface_file(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.surface_file.encode()));
                # pdll.st_element_interaction(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(self.interaction)); #/* 1=refract, 2=reflect */
                # pdll.st_element_optic(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.optic.name.encode()));

                # return 1
                pass

            def surface_spherical(self, radius):
                # """
                # Set up the surface as spherical type.

                # Surface centroid is at x=0, y=0, z=radius.

                # Parameters
                # ==========
                # radius : float
                #     Radius of the spherical surface
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = 1. / radius
                # self.surface = 's'
                # return True
                pass

            def surface_parabolic(self, focal_len_x, focal_len_y):
                # """
                # Set up the surface as parabolic.

                # Surface function is:
                #     Z(x,y) = 1/2 * (c_x * x^2 + c_y * y^2)
                #     where
                #     c_x = 1 / (2 * focal_len_x)
                #     c_y = 1 / (2 * focal_len_y)

                # The surface value is z=0 at x=y=0.

                # Parameters
                # ==========
                # focal_len_x : float
                #     Focal length of the surface in the x-direction. If infinite, use float('inf')
                # focal_len_y : float
                #     Focal length of the surface in the y-direction. If infinite, use float('inf')
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = 1. / (2.*focal_len_x)
                # if focal_len_y != 0.0:
                #     self.surface_params[1] = 1. / (2.*focal_len_y)
                # else:
                #     self.surface_params[1] = 0.0
                # self.surface = 'p'
                # return True
                pass

            def surface_flat(self):
                # """
                # Set up the surface as flat
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'f'
                # return True
                pass

            def surface_hypellip(self, vertex_curv, kappa):
                # """
                # Set up the surface described by equation:
                #     Z(x,y) = ( vertex_curv*(x^2 + y^2) ) /
                #                 (1 + sqrt(1-kappa*vertex_curv^2*(x^2 + y^2)))
                # Parameters
                # ----------
                # vertex_curv
                #     Curvature parameter
                # kappa
                #     Form parameter. Value of parameter determines geometry as follows:
                #     kappa < 0 --> tall hyperboloid
                #     kappa 0..1 --> ellipsoid
                #     kappa > 1 --> stout ellipsoid
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = vertex_curv
                # self.surface_params[1] = kappa
                # self.surface = 'o'
                # return True
                pass

            def surface_conical(self, theta):
                # """
                # Set up the surface described by cone with half-angle theta.

                # The axis of the cone coincides with the z-axis. The function of the surface is:
                #     Z(x,y) = sqrt(x^2 + y^2)/tan(theta)

                # Parameters
                # ----------
                # theta : float
                #     (degrees) half-angle of cone
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = theta
                # self.surface = 'c'
                # return True
                pass

            def surface_cylindrical(self, radius):
                # """
                # Set up the surface as cylindrical.

                # The surface centroid is located at x=0, y=0, z=radius. The cylinder's axis
                # is parallel to the Y-axis.

                # Parameters
                # ----------
                # radius
                #     Radius of the cylinder
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = 1./radius
                # self.surface = 't'
                # return True
                pass

            def surface_toroid(self, rad_annulus, rad_ring):
                # """
                # Set up the surface as a toroid "donut".

                # Parameters
                # ----------
                # rad_annulus
                #     Radius of the 'tube', the distance between the min and max radii of the torus
                # rad_ring
                #     The radius of the centerpoint of the annular tube
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface_params[0] = rad_annulus
                # self.surface_params[1] = rad_ring
                # self.surface = 'd'
                # return True
                pass

            def surface_zernicke(self, file_path):
                # """
                # Set up the surface from a file as a Zernicke surface, where the surface is described by the equation:
                # Z(x,y) = sum_i=0^N
                #             sum_j=0^i  Bi,j * x^j * y^(i-j)

                # Accepts *mon file extension specifying the Zernicke coefficients.
                # File format should be a single data column:
                #     N
                #     B0,0
                #     B1,0
                #     B1,1
                #     B2,1
                #     B2,2
                #     B2,3
                #     ...
                #     BN,N

                # Parameters
                # ----------
                # file_path
                #     Path to the file containing the data.
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'm'
                # self.surface_file = file_path
                # return True
                pass

            def surface_polynomialrev(self, file_path):
                # """
                # Set up the surface from a file as a rotationally symmetric polynomial, where the surface is described by
                # the equation:
                # Z(r) = sum_i=0^N  C_i * r^i,  where r=sqrt(x^2 + y^2)

                # Accepts *ply file extension specifying equation coefficients.
                # File format should be a single data column:
                #     N
                #     C0
                #     C1
                #     C2
                #     ...
                #     C,N

                # Parameters
                # ----------
                # file_path
                #     Path to the file containing the data.
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'r'
                # self.surface_file = file_path
                # return True
                pass

            def surface_cubicspline(self, file_path):
                # """
                # Set up the surface from a file as a rotationally symmetric cubic spline. Accepts *csi file extension.
                # File format should be two tab-separated columns:
                #     N
                #     r1      Z1
                #     r2      Z2
                #     r3      Z3
                #     ...
                #     rN      ZN
                #     dZ/dr1  dZ/drN

                # Parameters
                # ----------
                # file_path
                #     Path to the file containing the data.
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'i'
                # self.surface_file = file_path
                # return True
                pass

            def surface_finiteelement(self, file_path):
                # """
                # Set up the surface from a file using finite element data specifying the vertices of the elements in
                # x,y,z coordinates.

                # Accepts the *.fed file extension. File format should be 3 tab-separated
                # columns:
                #     N
                #     x1      y1      z1
                #     x2      y2      z2
                #     x3      y3      z3
                #     ...
                #     xN      yN      zN

                # Parameters
                # ----------
                # file_path
                #     Path to the file containing the data.
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'e'
                # self.surface_file = file_path
                # return True
                pass

            def surface_vshot(self, file_path):
                # """
                # Set up the surface from a file using VSHOT data specifying matrix coefficients generated by a VSHOT test.

                # Accepts the *.sht file extension. File format should be:
                #     First line - file name (skipped)
                #     Radius      Focal length        Target-dist
                #     0           order               num points
                #     rmsslope    rmsscale
                #     b00
                #     b10
                #     b11
                #     b20
                #     b21
                #     b22
                #     ...
                #     bDD    || where 'D' is order
                #     a1      b1      c1      d1      e1
                #     a2      b2      c2      d2      e2
                #     a3      b3      c3      d3      e3
                #     ...
                #     aN      bN      cN      dN      eN  || where 'N' is num points

                # Parameters
                # ----------
                # file_path
                #     Path to the file containing the data.
                # """
                # self.surface_params = [0. for i in range(8)]
                # self.surface = 'v'
                # self.surface_file = file_path
                # return True
                pass

            # ---------------------
            def aperture_circle(self, diameter):
                # """
                # Set up the aperture as circular with 'diameter'.

                # Aim: The X and Y directions lie in the plane of the circle. Z is normal to the plane.

                # Parameters
                # ----------
                # diameter
                #     Diameter of the circle
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = diameter
                # self.aperture = 'c'
                # return True
                pass

            def aperture_hexagon(self, diameter):
                # """
                # Set up the aperture as a hexagon centered at x=0,y=0. The hexagon is circumscribed
                # by a circle of 'diameter'.

                # Aim: The X and Y directions lie in the plane of the hexagon. X crosses through a vertex
                # between two segments, while Y bisects an edge segment. Z is normal to the plane.

                #  y^
                #  __
                # /  \ x->
                # \__/

                # Parameters
                # ----------
                # diameter
                #     Diameter of the circumscribing circle.
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = diameter
                # self.aperture = 'h'
                # return True
                pass

            def aperture_triangle(self, diameter):
                # """
                # Set up the aperture as a equilateral triangle with centroid at x=0,y=0. The triangle
                # is circumscribed by a circle of 'diameter'.

                # Aim: The X and Y directions lie in the plane of the triangle. Y crosses through a vertex between
                # two segments, while X crosses at an intermediate position on one leg of the triangle. Z is normal
                # to the plane. The coordinates are centered at the middle of a circle of diameter 'D' that circumscribes
                # the isoceles triangle

                #   y^
                #   /\
                #  /  \  x ->
                # /____\

                # Parameters
                # ----------
                # diameter
                #     Diameter of the circumscribing circle.
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = diameter
                # self.aperture = 't'
                # return True
                pass

            def aperture_rectangle(self, length_x, length_y):
                # """
                # Set up the aperture as a rectangle.

                # Aim: The X and Y directions lie in the plane of the rectangle. Y crosses bisects a horzontal leg
                # of width 'W', while X bisects a vertical leg of height 'H'. Z is normal to the plane. The
                # coordinates are centered x=W/2, y=H/2.

                # Parameters
                # ----------
                # length_x
                #     Width in x-coordinate direction
                # length_y
                #     Height in y-coordinate direction
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = length_x
                # self.aperture_params[1] = length_y
                # self.aperture = 'r'
                # return True
                pass

            def aperture_annulus(self, r_inner, r_outer, theta):
                # """
                # Set up the aperture as annular, where aperture is the annulus between to specified radii
                # and within an angular slice 'theta' which is centered around the x-axis.

                # Aim: The X and Y directions lie in the plane of the annulus. Z is normal to the plane.

                # Parameters
                # ----------
                # r_inner
                #     Inner radius of annular region
                # r_outer
                #     Outer radius of annular region
                # theta : deg
                #     Slice of the circle contained, centered around x-axis
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = r_inner
                # self.aperture_params[1] = r_outer
                # self.aperture_params[2] = theta
                # self.aperture = 'a'
                # return True
                pass

            def aperture_singleax_curve(self, x1, x2, L):
                # """
                # Set up the aperture as revolved around a single axis. Revolved window is between two
                # coordinates x1->x1, both non-negative and with x2 > x1. The aperture has
                # length 'L' in the y-direction.

                # This aperture is often used with a cylindrical surface. In this case,
                # both x1 and x2 should be zero, and the cylinder height specified with 'L'.

                # Aim: X and Z follow radial lines and cross through the curvature section. Y lies along the
                # centerline/axis of the cylindrical section at X=0, Z=0. The radial positions are with
                # respect to the X and Z coordinates.

                # ^ y
                # |    ___  ....L
                # |   |   |
                # |---|---|---> X
                # |   |___| ....
                # |   x1  x2

                # Parameters
                # ----------
                # x1
                #     inner coordinate of revolved section
                # x2
                #     outer coordinate of revolved section
                # L
                #     length of revolved section along axis of revolution
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = x1
                # self.aperture_params[1] = x2
                # self.aperture_params[2] = L
                # self.aperture = 'l'
                # return True
                pass
            
            def aperture_irr_triangle(self, x1, y1, x2, y2, x3, y3):
                # """
                # Set up the aperture as a triangle given by three (x,y) coordinate pairs.

                # Aim: X and Y are in the plane containing the coordinates. Z is normal to the plane.

                # Parameters
                # ----------
                # x1
                #     x-coordinate, point 1
                # y1
                #     y-coordinate, point 1
                # x2
                #     x-coordinate, point 2
                # y2
                #     y-coordinate, point 2
                # x3
                #     x-coordinate, point 3
                # y3
                #     y-coordinate, point 3
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = x1
                # self.aperture_params[1] = y1
                # self.aperture_params[2] = x2
                # self.aperture_params[3] = y2
                # self.aperture_params[4] = x3
                # self.aperture_params[5] = y3
                # self.aperture = 'i'
                # return True
                pass

            def aperture_quadrilateral(self, x1, y1, x2, y2, x3, y3, x4, y4):
                # """
                # Set up the aperture as a quadrilateral given by four (x,y) coordinate pairs.

                # Aim: X and Y are in the plane containing the coordinates. Z is normal to the plane.

                # Parameters
                # ----------
                # x1
                #     x-coordinate, point 1
                # y1
                #     y-coordinate, point 1
                # x2
                #     x-coordinate, point 2
                # y2
                #     y-coordinate, point 2
                # x3
                #     x-coordinate, point 3
                # y3
                #     y-coordinate, point 3
                # x4
                #     x-coordinate, point 4
                # y4
                #     y-coordinate, point 4
                # """
                # self.aperture_params = [0. for i in range(8)]
                # self.aperture_params[0] = x1
                # self.aperture_params[1] = y1
                # self.aperture_params[2] = x2
                # self.aperture_params[3] = y2
                # self.aperture_params[4] = x3
                # self.aperture_params[5] = y3
                # self.aperture_params[6] = x4
                # self.aperture_params[7] = y4
                # self.aperture = 'q'
                # return True
                pass
        # -------------------------- end Element class ---------------------------------


        # -----------methods of the 'Stage' class --------------------------------------
        def __init__(self, id : int):
            # """
            # """
            # ## Identifying integer associated with the stage
            # self.id = id
            # ## Stage location in global coordinates
            # self.position = Point()
            # ## Coordinate system aim point in global coordinates
            # self.aim = Point()
            # self.aim.z = 1
            # ## [deg] Rotation of coordinate system around z-axis
            # self.zrot = 0.
            # ## Flag indicating virtual stage
            # self.is_virtual = False
            # ## Flag indicating that rays can have multiple interactions within a single stage.
            # self.is_multihit = True
            # ## Flag indicating the stage is in trace-through mode
            # self.is_tracethrough = False
            # ## Descriptive name for this stage
            # self.name = "stage_{:d}".format(id)

            # ## list of all elements in the stage
            # self.elements = []
            # return
            pass

        def copy(self, snew):
            # """
            # Deep copy of the current Stage instance

            # Inputs
            # ---------
            # snew : Stage
            #     Reference to new Stage object to which data will be copied
            # """
            # c = self.__dict__.copy()
            # for attr in self.__dict__.keys():
            #     if attr in ['elements', '_pdll', '_p_data']:
            #         continue
            #     elif attr in ['position','aim']:  #points
            #         snew.__setattr__(attr, c[attr].copy())
            #     else:
            #         snew.__setattr__(attr, copy.deepcopy(c[attr]))

            # snew.elements = [PySolTrace.Stage.Element(snew, el.id) for el in self.elements]
            # for i in range(len(self.elements)):
            #     self.elements[i].copy(snew.elements[i])

            # return
            pass

        def Create(self, pdll, p_data) -> int:
            # """
            # Create Stage instance in the SolTrace context.
            # Note: This does not create any associated Elements, which must have their Create method called separately.

            # Returns
            # ----------
            # int
            #     1 if successful, 0 otherwise
            # """

            # pdll.st_add_stage.restype = c_int
            # pdll.st_add_stage(c_void_p(p_data) )

            # pdll.st_stage_xyz.restype = c_int
            # pdll.st_stage_aim.restype = c_int
            # pdll.st_stage_zrot.restype = c_int
            # pdll.st_stage_flags.restype = c_int

            # pdll.st_stage_xyz(c_void_p(p_data), c_uint32(self.id), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))
            # pdll.st_stage_aim(c_void_p(p_data), c_uint32(self.id), c_number(self.aim.x), c_number(self.aim.y), c_number(self.aim.z))
            # pdll.st_stage_zrot(c_void_p(p_data), c_uint32(self.id), c_number(self.zrot))
            # pdll.st_stage_flags(c_void_p(p_data), c_uint32(self.id), c_int(int(self.is_virtual)), c_int(int(self.is_multihit)), c_int(int(self.is_tracethrough)))

            # for element in self.elements:
            #     element.Create(pdll,p_data)

            # return 1
            pass

        def add_element(self) -> int:
            # """
            # Add one element to the stage. This method appends an Element object to the
            # stage's Stage.elements list.
            # To update element properties and settings, call the Element.Create method
            # on each element.

            # Returns
            # ----------
            # PySolTrace.Stage.Element
            #     Reference to the newly created element
            # """

            # new_e = PySolTrace.Stage.Element(self, len(self.elements) )
            # self.elements.append( new_e )
            # return new_e
            pass

    # ---------- methods of the PySolTrace class --------------------------------------------
    def __init__(self):
        # # Initialize lists for optics and stage instances
        # ## List of Optics instances
        # self.optics = []
        # ## List of Stage instances
        # self.stages = []
        # ## sun
        # self.sun = None

        # # Simulation settings
        # ## Minimum number of simulation ray hits
        # self.num_ray_hits = int(1e5)
        # ## Maximum number of ray hits in a simulation
        # self.max_rays_traced = self.num_ray_hits*100
        # ## Flag indicating whether sunshape should be included
        # self.is_sunshape = True
        # ## Flag indicating whether surface errors should be included
        # self.is_surface_errors = True
        # # Placeholder for output ray data
        # self.raydata = None
        # # Placeholder for sunstats data
        # self.sunstats = None
        # # Placeholder for power per ray
        # self.powerperray = None
        # # Direct normal irradince for calculations
        # self.dni = 1000.  #w/m^2
        pass

    def copy(self):
        # """
        # Deep copy of the current PySolTrace instance

        # Returns:
        # ========
        # Copy of the current PySolTrace instance

        # """
        # psnew = PySolTrace()
        # for attr in ['num_ray_hits','max_rays_traced','is_sunshape','is_surface_errors']:
        #     psnew.__setattr__(attr, copy.deepcopy(self.__getattribute__(attr)))
        # psnew.stages = [PySolTrace.Stage(st.id) for st in self.stages]
        # for i in range(len(self.stages)):
        #     self.stages[i].copy(psnew.stages[i])
        # psnew.optics = [PySolTrace.Optics(op.id) for op in self.optics]
        # for i in range(len(self.optics)):
        #     self.optics[i].copy(psnew.optics[i])
        # psnew.sun = PySolTrace.Sun()
        # self.sun.copy(psnew.sun)

        # #re-link optics and elements
        # opt_map = {}
        # for opt in psnew.optics:
        #     opt_map[opt.id] = opt

        # for i,stage in enumerate(self.stages):
        #     for j,element in enumerate(stage.elements):
        #         psnew.stages[i].elements[j].optic = opt_map[element.optic.id]

        # return psnew
        pass

    def Create(self, pdll, p_data):
        # """
        # Create soltrace context from data structures
        # """

        # self.sun.Create(pdll, p_data)
        # for opt in self.optics:
        #     opt.Create(pdll, p_data)
        # for stage in self.stages:
        #     stage.Create(pdll, p_data)
        pass

    def add_optic(self, optic_name : str):
        # """
        # Instantiates a new PySolTrace.Optics object, adding it to the optics list.
        # This method does not set optics properties, which instead is done using the Optics.Create method.

        # Parameters
        # ----------
        # optic_name : string
        #     Unique name for this Optics instance.

        # Returns
        # ----------
        # Optics
        #     Reference to the Optics object that was just created.
        # """

        # new_opt_id = len(self.optics)

        # self.optics.append( PySolTrace.Optics(new_opt_id ) )
        # self.optics[-1].name = optic_name

        # return self.optics[-1]  #return the last object in the list, which was the one just created
        pass

    def delete_optic(self, optic_id : int) -> int:
        # """
        # Delete Optics instance. The optics object is removed from the PySolTrace.optics list and
        # from the SolTrace context.

        # Parameters
        # ----------
        # optic_id : int
        #     ID associated with the optics to be deleted

        # Returns
        # ----------
        # int
        #     1 if successful, 0 otherwise
        # """

        # # find the appropriate optic
        # for opt in self.optics:
        #     if optic_id == opt.id:
        #         # clear it from the optics array
        #         self.optics.remove(opt)
        #         # Remove from the soltrace context
        #         # self._pdll.st_delete_optic.restype = c_int
        #         # return self._pdll.st_delete_optic(c_void_p(self._p_data), c_uint32(optic_id) )

        # # If reaching this point, the optic id was not found
        # return 0
        pass

    def add_sun(self):
        # """
        # Instantiates a PySolTrace.Sun object and associates it with the PySolTrace.sun member.
        # This does not create or modify the Sun data in the SolTrace context.

        # Returns
        # ----------
        # PySolTrace.Sun
        #     Reference to newly created Sun instance.
        # """
        # ## Object containing Sun class data
        # self.sun = PySolTrace.Sun()
        # return self.sun
        pass

    def add_stage(self):
        # """
        # Adds a new Stage instance to the PySolTrace.stages list. The Stage ID is automatically generated based on the number
        # of current stages.

        # Returns
        # ----------
        # PySolTrace.Stage
        #     Reference to the newly created Stage object.
        # """
        # new_st_id = len(self.stages)

        # self.stages.append( PySolTrace.Stage( new_st_id ) )

        # return self.stages[-1]
        pass

    def delete_stage(self, stage_id : int) -> int:
        # """
        # Delete Stage instance. The stage object is removed from the PySolTrace.stages list and
        # from the SolTrace context.

        # Parameters
        # ----------
        # stage_id : int
        #     ID associated with the stage to be deleted

        # Returns
        # ----------
        # int
        #     1 if successful, 0 otherwise
        # """
        # # find the appropriate optic
        # for st in self.stages:
        #     if stage_id == st.id:
        #         # clear it from the optics array
        #         self.stages.remove(st)
        #         # Remove from the soltrace context
        #         # self._pdll.st_delete_stage.restype = c_int
        #         # return self._pdll.st_delete_stage(c_void_p(self._p_data), c_uint32(stage_id) )

        # # If reaching this point, the stage id was not found
        # return 0
        pass

    def __load_dll(self):
        # cwd = os.getcwd()
        # if sys.platform == 'win32' or sys.platform == 'cygwin':
        #     ## loaded SolTrace library of exported functions
        #     pdll = CDLL(cwd + "/coretrace_api.dll")
        #     # print("Loaded win32")
        #     #pdll = CDLL(cwd + "/coretraced.dll") # for debugging
        # elif sys.platform == 'darwin':
        #     pdll = CDLL(cwd + "/coretrace_api.dylib")  # Never tested
        # elif sys.platform.startswith('linux'):
        #     pdll = CDLL(cwd +"/coretrace_api.so")
        # else:
        #     print( 'Platform not supported ', sys.platform)
        # return pdll
        pass

    def run(self, seed : int = -1, as_power_tower = False, nthread=1, thread_id=0, no_callback=False):
        # """
        # Run SolTrace simulation.

        # If calling this function in multithread mode, note that the run() function
        # **must** be called inside an import guard, e.g.:
        # > if __name__ == "__main__":
        # >    mypst_obj.run(...)
        # Otherwise, you'll receive an error.

        # Parameters
        # ----------
        # seed : int
        #     Seed for random number generator. [-1] for random seed. Seeding happens
        #     differently for single vs multi-thread modes.
        #         * If nthreads == 1 and seed < 0: a random int is chosen as the seed value.
        #         * If nthreads > 1 and seed < 0: a random int is chosen for the first
        #           thread seed value. Other threads i=1..(nthreads-1) are assigned
        #           (first value) + i*123.
        # as_power_tower : bool
        #     Flag indicating simulation should be processed as power
        #     tower / central receiver type, with corresponding efficiency adjustments.
        # nthread : int
        #     Number of threads to execute. Will be limited by the method to the number
        #     available on the machine.
        #         * If nthreads > 1, the function will call recursively while setting
        #           nthreads=0 for each thread spawned.
        #         * If nthreads == 1, the function will run in single-thread mode. Seed
        #           values are checked.
        #         * If nthreads == 0, the function will run in single-thread mode. Seed
        #           values are not checked and should be handled prior to calling in
        #           this mode.
        # thread_id : int
        #     Argument used by the multi-threading call. Do not manually specify this value.

        # Returns
        # ----------
        # int
        #     Simulation return value
        # """

        # pdll = self.__load_dll()

        # if seed<0:
        #     runseed = random.randint(1,int(1e9))
        # else:
        #     runseed = seed

        # if nthread in [0,1]:

        #     # Create an instance of soltrace in memory
        #     pdll.st_create_context.restype = c_void_p
        #     p_data = pdll.st_create_context()

        #     self.Create(pdll, p_data)

        #     pdll.st_sim_errors.restype = c_int
        #     pdll.st_sim_errors(c_void_p(p_data), c_int(1 if self.is_sunshape else 0), c_int(1 if self.is_surface_errors else 0))

        #     pdll.st_sim_params.restype = c_int
        #     pdll.st_sim_params(c_void_p(p_data), c_int(int(self.num_ray_hits)), c_int(int(self.max_rays_traced)), c_int(as_power_tower))

        #     if thread_id == 0:
        #         tstart = time.time()

        #     pdll.st_sim_run.restype = c_int
        #     if no_callback:
        #         res = pdll.st_sim_run( c_void_p(p_data), c_uint16(runseed), no_api_callback, thread_id)
        #     else:
        #         res = pdll.st_sim_run( c_void_p(p_data), c_uint16(runseed), api_callback, thread_id)
        #         if thread_id == 0:
        #             print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

        #     # Collect simulation output, including raw ray data and sunbox stats
        #     self.raydata = self.__get_ray_dataframe(pdll,p_data)
        #     self.sunstats = self.__get_sun_stats(pdll, p_data)
        #     # Compute and save power per ray
        #     self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni

        #     pdll.st_free_context.restype = c_bool
        #     pdll.st_free_context(c_void_p(p_data))

        #     return res
        # else:
        #     seeds = [seed + i*123 for i in range(nthread)]

        #     P = [[self.copy(), as_power_tower, seeds[i], i+1, no_callback] for i in range(nthread)]

        #     # modify the number of rays to match the required totals
        #     nrpt = int(float(self.num_ray_hits)/float(nthread))
        #     mrpt = int(float(self.max_rays_traced)/float(nthread))

        #     for p in P:
        #         p[0].num_ray_hits = nrpt
        #         p[0].max_rays_traced = mrpt
        #         if p == P[0]:
        #             p[0].num_ray_hits += int(float(self.num_ray_hits) % float(nthread))
        #             p[0].max_rays_traced += int(float(self.max_rays_traced) % float(nthread))

        #     pool = multiprocessing.Pool(nthread)
        #     if not no_callback:
        #         print("Launching {:d} threads...".format(nthread))
        #     tstart = time.time()
        #     res = pool.starmap_async(_thread_func, P)
        #     pool.close()
        #     pool.join()
        #     if not no_callback:
        #         print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

        #     # Modify the ray number for threads 2+ to avoid duplication
        #     try:
        #         dfs = [r[0] for r in res.get()]
        #         rstart = int(dfs[0].iloc[-1].number)
        #     except:
        #         print("Unknown error caused the simulation to fail. Try re-running.")
        #         return

        #     if len(dfs)>1:
        #         for d in dfs[1:]:
        #             d.number = d.number+rstart
        #             rstart = d.number.iloc[-1]

        #     self.raydata = pd.concat(dfs)
        #     self.raydata.reset_index(inplace=True)
        #     self.sunstats = res.get()[0][1]  #take the first thread result
        #     # add up all the sunrays from all threads
        #     srct = 0
        #     for r in res.get():
        #         srct += r[1]['nsunrays']
        #     self.sunstats['nsunrays'] = srct

        #     # Compute and save power per ray [W]
        #     self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni

        #     return 1
        pass

    def __get_num_intersections(self, pdll, p_data) -> int:
        # """
        # [Post simulation] Get the number of ray intersections detected in the simulation.

        # Returns
        # ----------
        # int
        #     Number of intersections
        # """

        # if p_data == 0:
        #     return 0

        # pdll.st_num_intersections.restype = c_int
        # return pdll.st_num_intersections(c_void_p(p_data))
        pass

    def __get_sun_stats(self, pdll, p_data):
        # """
        # Get information on the sun box.

        # Returns
        # ----------
        # dict
        #     Keys in the return dictionary are:
        #     'xmin' --> Minimum x extent of the bounding box for hit testing
        #     'xmax' --> Maximum x extent of the bounding box for hit testing
        #     'ymin' --> Minimum y extent of the bounding box for hit testing
        #     'ymax' --> Maximum y extent of the bounding box for hit testing
        #     'nsunrays' --> Number of sun rays simulated
        # """
        # if p_data == 0:
        #     raise "SolTrace context not assigned"

        # xmin = (c_number)()
        # xmax = (c_number)()
        # ymin = (c_number)()
        # ymax = (c_number)()
        # nsunrays = (c_int)()

        # pdll.st_sun_stats.restype = c_int
        # pdll.st_sun_stats(c_void_p(p_data), pointer(xmin), pointer(xmax), pointer(ymin), pointer(ymax), pointer(nsunrays))

        # return {
        #     'xmin':float(xmin.value),
        #     'xmax':float(xmax.value),
        #     'ymin':float(ymin.value),
        #     'ymax':float(ymax.value),
        #     'nsunrays':int(nsunrays.value),
        # }
        pass

    def __get_ray_dataframe(self, pdll, p_data):
        # """
        # Get a pandas dataframe with all of the ray data from the simulation.

        # Returns
        # ----------
        # Pandas.DataFrame
        #     with columns:
        #     loc_x   | Ray hit location, x-coordinate
        #     loc_y   | Ray hit location, y-coordinate
        #     loc_z   | Ray hit location, z-coordinate
        #     cos_x   | Ray directional vector, x-component
        #     cos_y   | Ray directional vector, y-component
        #     cos_z   | Ray directional vector, z-component
        #     element | Element associated with ray hit
        #     stage   | Stage associated with ray hit
        #     number  | Ray number
        # """
        # if p_data == 0:
        #     raise "SolTrace context not assigned"

        # data = {}

        # n_int = self.__get_num_intersections(pdll, p_data)
        # #print("Returning  {:d} intersections...".format(n_int))
        # data['loc_x'] = (c_number*n_int)()
        # data['loc_y'] = (c_number*n_int)()
        # data['loc_z'] = (c_number*n_int)()

        # pdll.st_locations.restype = c_int
        # pdll.st_locations(c_void_p(p_data), pointer(data['loc_x']), pointer(data['loc_y']), pointer(data['loc_z']))

        # data['cos_x'] = (c_number*n_int)()
        # data['cos_y'] = (c_number*n_int)()
        # data['cos_z'] = (c_number*n_int)()

        # pdll.st_cosines.restype = c_int
        # pdll.st_cosines(c_void_p(p_data), pointer(data['cos_x']), pointer(data['cos_y']), pointer(data['cos_z']))

        # data['element'] = (c_int*n_int)()

        # pdll.st_elementmap.restype = c_int
        # pdll.st_elementmap(c_void_p(p_data), pointer(data['element']))

        # data['stage'] = (c_int*n_int)()

        # pdll.st_stagemap.restype = c_int
        # pdll.st_stagemap(c_void_p(p_data), pointer(data['stage']))

        # data['number'] = (c_int*n_int)()

        # pdll.st_raynumbers.restype = c_int
        # pdll.st_raynumbers(c_void_p(p_data), pointer(data['number']))

        # for key in data.keys():
        #     data[key] = list(data[key])

        # df = pd.DataFrame(data)

        # return df
        pass

    def plot_trace(self, nrays:int = 100000, ntrace:int=100, show_sun_vector:bool=True):
        # """
        # Creates and (optionally) displays a 3D scatter and trace plot. This
        # function requires that the Python package `plotly` be installed.

        # Parameters
        # ------------
        # nrays : int
        #     Number of individual rays to include in the scatter plot. Very
        #     large values may render slowly.
        # ntrace : int
        #     Number of rays for which traces will be displayed. Large values
        #     may render slowly
        # show_sun_vector : bool
        #     Flag indicating whether the sun vector should be rendered on the plot
        # """

        # print("Generating 3D trace plots")
        # # Plotting with plotly
        # try:
        #     import plotly.graph_objects as go
        # except:
        #     raise RuntimeError("Missing library: plotly. \n Trace plotting requires the Plotly library to be installed. [$ pip install plotly]")

        # df = self.raydata

        # # Choose how many points to plot.
        # nn = min(nrays, len(df))
        # inds = numpy.random.choice(range(len(df)), size=nn, replace=False)

        # # Data for a three-dimensional line. Randomly choose points if fewer than the full amount are desired.
        # loc_x = df.loc_x.values[inds]
        # loc_y = df.loc_y.values[inds]
        # loc_z = df.loc_z.values[inds]
        # stage = df.stage.values[inds]
        # raynum = df.number.values[inds]

        # # Generate the 3D scatter plot
        # layout = go.Layout(scene=dict(aspectmode='data'))

        # if len(list(set(stage))) > 1:
        #     md = dict( size=0.75, color=stage, colorscale='jet', opacity=0.7, )
        # else:
        #     md = dict( size=0.75, color='black', opacity=0.7, )

        # fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=md ), layout=layout )

        # # Generate line traces for a subset of randomly selected rays
        # for i in numpy.random.choice(raynum, size=50, replace=False):
        #     dfr = df[df.number == i]    #find all rays numbered 'i'
        #     ray_x = dfr.loc_x
        #     ray_y = dfr.loc_y
        #     ray_z = dfr.loc_z
        #     fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))
        # # Add a trace for the sun vector
        # if show_sun_vector:
        #     tmp = df[df.stage==1].iloc[0]  #sun is coming from the cos vector of the elements in the first stage. Just take the first.
        #     sun_vec = numpy.array([-tmp.cos_x,-tmp.cos_y,-tmp.cos_z])  #negative of the vector
        #     # scale the vector based on the overall size of the sun bounding box
        #     sunrange = numpy.array([df.loc_x.max()-df.loc_x.min(), df.loc_y.max()-df.loc_y.min(), df.loc_z.max()-df.loc_z.min()])
        #     # sun_scale = ((self.sunstats['xmax']-self.sunstats['xmin'])**2 + (self.sunstats['ymax']-self.sunstats['ymin'])**2)**.5 *0.75
        #     # sun_scale = min([sun_scale, df.loc_x.max()])
        #     sun_vec *= (sunrange*sun_vec).max()
        #     fig.add_trace(go.Scatter3d(x=[0,sun_vec[0]], y=[0,sun_vec[1]], z=[0,sun_vec[2]], mode='lines', line=dict(color='orange', width=3)))
        #     fig.add_trace(go.Scatter3d(x=[0,sun_vec[0]], y=[0,sun_vec[1]], z=[0,0], mode='lines', line=dict(color='gray', width=2)))

        # fig.update_layout(showlegend=False)
        # fig.show()

        # return
        pass

    def plot_flux(self, element, nx:int = 25, ny:int = 25, figpath:str=None, display=True, figsize=(9,6),
                  absorbed_only:bool = True, levels=25, dpi:int=300, xlabel:str=None, ylabel:str = None):
        # """
        # Creates and (optionally) displays a flux plot for a given stage element.

        # Parameters
        # ----------
        # element : PySolTrace:Stage:Element
        #     Reference to the element for which the plot will be generated
        # nx : int (default 25)
        #     Number of flux bins along the aperture x-coordinate
        # ny : int (default 25)
        #     Number of flux bins along the aperture y-coordinate
        # figpath : str (default None)
        #     Path to file location where figure will be saved. If None, figure is not saved.
        # display : bool (default True)
        #     Flag indicating whether the figure should be displayed at runtime
        # figsize : tuple (default (9,6))
        #     Figure size in inches
        # absorbed_only : bool (default True)
        #     Only include rays that are absorbed by the element, omitting reflected rays
        # levels : int (default 25)
        #     Number of contour levels to include in the flux map
        # dpi : int (default 300)
        #     Resolution of the saved image
        # xlabel : str (default None)
        #     String specifying label to use on x-axis of plot
        # ylabel : str (default None)
        #     String specifying label to use on y-axis of plot

        # Returns
        # ------------
        # None
        # """

        # flux_st = self.bin_rays(element, nx, ny, absorbed_only)

        # el_id = element.id+1
        # st_id = element.stage_id+1

        # # plotting specifics for each surface type
        # if element.surface == 'f':
        #     # Flat
        #     W,H = element.aperture_params[0:2]
        #     x_rec = numpy.arange(0, W, W/nx)
        #     y_rec = numpy.arange(0, H, H/ny)
        #     xlabtemp = "X-axis position"

        # elif element.surface == 't':
        #     # Cylindrical
        #     D = 2./element.surface_params[0]
        #     H = element.aperture_params[2]
        #     x_rec = numpy.arange(0, numpy.pi*D, numpy.pi*D/nx)
        #     y_rec = numpy.arange(-H/2,H/2, H/ny)
        #     # plot label for later
        #     xlabtemp = "Circumferential position"

        # # check labels
        # if ylabel != None:
        #     ylabtemp = ylabel
        # else:
        #     ylabtemp = "Y-axis position"
        # if xlabel != None:
        #     xlabtemp = xlabel

        # # Generate new plot
        # plt.figure(figsize=figsize)
        # plt.title(f"Flux intensity: element {el_id}, stage {st_id}")
        # Xr,Yr = numpy.meshgrid(y_rec, x_rec)
        # plt.contourf(Yr, Xr, flux_st, levels=levels)
        # plt.colorbar()
        # plt.title(f"Stage {st_id}/Elem. {el_id} | Max {flux_st.max():.0f} | Mean {flux_st.mean():.1f}")
        # plt.xlabel(xlabtemp)
        # plt.ylabel(ylabtemp)
        # plt.tight_layout()
        # if figpath:
        #     plt.savefig(figpath, dpi=dpi)
        # if display:
        #     plt.show()
        # return
        pass

    def bin_rays(self, element, nx:int = 25, ny:int = 25, absorbed_only:bool = True):
        # """
        # Bins rays for plotting flux maps.

        # Parameters
        # ----------
        # element : PySolTrace:Stage:Element
        #     Reference to the element for which the plot will be generated
        # nx : int (default 25)
        #     Number of flux bins along the aperture x-coordinate
        # ny : int (default 25)
        #     Number of flux bins along the aperture y-coordinate
        # absorbed_only : bool (default True)
        #     Only include rays that are absorbed by the element, omitting reflected rays

        # Returns
        # ------------
        # flux_map : numpy.Array(nx, ny)
        #     Flux map of element
        # """
        # # Get a pandas dataframe with all of the ray data
        # df = self.raydata

        # if self.raydata.empty:
        #     raise(RuntimeError("Flux plot not created: no ray data available"))

        # # Check if surface type is supported
        # if element.surface not in ['f', 't']:
        #     raise(RuntimeError(f"Surface type {element.surface} is not supported for flux plot generation. Must be one of 'f', 't'."))

        # el_id = element.id+1
        # st_id = element.stage_id+1

        # dfr = df[df.stage==st_id]
        # if absorbed_only:
        #     dfr = dfr[dfr.element==-el_id]  #absorbed rays
        # else:
        #     dfr = dfr[(dfr.element==-el_id) & (dfr.element==el_id)]  #absorbed and reflected rays

        # dfr = dfr.copy()

        # # Compute the euler angles for the target element
        # eu_angles = self.util_calc_euler_angles(numpy.array(element.position.as_list()), numpy.array(element.aim.as_list()), element.zrot)
        # # Compute the transform matrix
        # transforms = self.util_calc_transforms(eu_angles)

        # # initialize
        # loc = dfr[['loc_x','loc_y','loc_z']].to_numpy()
        # e_pos = numpy.array(element.position.as_list())

        # pos_t = self.util_transform_to_local(loc, numpy.array([0,0,1]), e_pos, transforms['rreftoloc'])['posloc']
        # dfr['loc_xt'] = pos_t.T[0]
        # dfr['loc_yt'] = pos_t.T[1]
        # dfr['loc_zt'] = pos_t.T[2]

        # flux_st = numpy.zeros((ny,nx))
        # # handle mapping differently for each surface type
        # if element.surface == 'f':
        #     # Flat
        #     W,H = element.aperture_params[0:2]
        #     raybins_x = numpy.floor((dfr.loc_xt + W/2)/W*nx).astype(int)
        #     raybins_y = numpy.floor((dfr.loc_yt + H/2)/H*ny).astype(int)

        #     dx = W / nx
        #     dy = H / ny

        # elif element.surface == 't':
        #     # Cylindrical
        #     D = 2./element.surface_params[0]
        #     H = element.aperture_params[2]

        #     # bin the rays circumferentially and vertically
        #     raybins_x = numpy.floor((numpy.arctan2(dfr.loc_xt, dfr.loc_zt-D/2)+math.pi)*nx/(2*math.pi)).astype(int)
        #     raybins_y = numpy.floor((dfr.loc_yt + H/2)/H*ny).astype(int)

        #     # Create the coordinate meshes
        #     dx = D*numpy.pi / nx
        #     dy = H / ny

        # # Compute power per ray (ppr) based on node area
        # anode = dx*dy
        # ppr = self.powerperray / anode

        # for r in range(len(raybins_x)):
        #     flux_st[raybins_x.values[r], raybins_y.values[r]] += ppr

        # return flux_st
        pass


    # /* utility transform/math functions */
    # def util_calc_euler_angles(self, origin : numpy.array, aimpoint : numpy.array, zrot) -> numpy.array:
        # """
        # Calculate the Euler angles associated with a given origin, aimpoint, and z-axis rotation.

        # Parameters
        # ----------
        # origin : [float,*3]
        #     Origin of the coordinate system
        # aimpoint : [float,*3]
        #     Aimpoint of the vector originating at the origin
        # zrot : float
        #     Rotation around the z-axis coordinate (degr)

        # Returns
        # ----------
        # list
        #     Calculated Euler angles (rad)
        # """

        # # This duplicates the built-in function but uses numpy operators directly
        # dv = aimpoint - origin
        # d = math.sqrt(sum(dv**2))
        # if d == 0:
        #     return
        # dv /= d
        # euler = numpy.array([
        #     math.atan2(dv[0],dv[2]),
        #     math.asin(dv[1]),
        #     zrot*0.017453292519943295 # acos(-1)/180.0
        # ])
        # return euler
        pass

    # def util_transform_to_local(self, posref : numpy.array, cosref : numpy.array, origin : numpy.array, rreftoloc : numpy.array):
        # """
        # Perform coordinate transformation from reference system to local system.

        # Parameters
        # ----------
        # PosRef : numpy.array([float,]*3)
        #     X,Y,Z coordinates of ray point in reference system
        # CosRef : numpy.array([float,]*3)
        #     Direction cosines of ray in reference system
        # Origin : numpy.array([float,]*3)
        #     X,Y,Z coordinates of origin of local system as measured in reference system
        # RRefToLoc : numpy.array([float,]*3)
        #     Rotation matrices required for coordinate transform from reference to local

        # Returns
        # ----------
        # (dict)  Keys in return dictionary include:
        #     posloc : ([float,]*3) X,Y,Z coordinates of ray point in local system
        #     cosloc : ([float,]*3) Direction cosines of ray in local system
        # """
        # assert type(rreftoloc) == type(numpy.array([]))
        # assert rreftoloc.shape == (3,3)

        # # This duplicates the built-in function but uses numpy operators directly
        # posdum = posref - origin
        # rreftoloc_m = rreftoloc.reshape((3,3))
        # posloc = numpy.dot(rreftoloc_m, posdum.T).T
        # cosloc = numpy.dot(rreftoloc_m, cosref.T).T

        # return {'cosloc':cosloc, 'posloc':posloc}
        pass

    def util_transform_to_ref(self, posloc, cosloc, origin, rloctoref):
        # """
        # Perform coordinate transformation from local system to reference system.

        # Parameters
        # ----------
        # PosLoc : [float,]*3
        #     X,Y,Z coordinates of ray point in local system
        # CosLoc : [float,]*3
        #     Direction cosines of ray in local system
        # Origin : [float,]*3
        #     X,Y,Z coordinates of origin of local system as measured in reference system
        # RLocToRef
        #     Rotation matrices required for coordinate transform from local to reference
        #     -- inverse of reference to local transformation

        # Returns
        # ----------
        # dict
        #     Keys in return dictionary include:
        #     posref : ([float,]*3) X,Y,Z coordinates of ray point in reference system
        #     cosref : ([float,]*3) Direction cosines of ray in reference system
        # """
        # assert type(rloctoref) == type(numpy.array([]))
        # assert rloctoref.shape == (3,3)

        # posdum = numpy.dot(rloctoref, posloc)
        # cosref = numpy.dot(rloctoref, cosloc, cosref)
        # posref = posdum + origin

        # return {'cosref':cosref, 'posref':posref}
        pass

    def util_matrix_vector_mult(self, m, v):
        # """
        # Perform multiplication of a 3x3 matrix and a length-3 vector, returning the result vector.

        # Parameters
        # ----------
        # m : array
        #     m[3][3] - a 3x3 matrix
        # v : array
        #     v[3] - a list, length 3

        # Returns
        # ----------
        # list
        #     m x v [3]
        # """
        # assert type(m) == type(numpy.array([]))
        # assert m.shape == (3,3)
        # assert type(v) == type(numpy.array([]))
        # assert m.shape == (3,)

        # return numpy.dot(m,v)
        pass

    def util_calc_transforms(self, euler):
        # """
        # Calculate matrix transforms

        # Parameters
        # ----------
        # euler : [float,]*3
        #     Euler angles

        # Returns
        # ----------
        # (dict) A dictionary containing the keys:
        #     rreftoloc : Transformation matrix from Reference to Local system
        #     rloctoref : Transformation matrix from Local to Reference system
        # """

        # pdll = self.__load_dll()

        # a_euler = (c_number*3)()
        # rreftoloc = (c_number*9)()
        # rloctoref = (c_number*9)()

        # a_euler[:] = euler

        # pdll.st_calc_transform_matrices.restype = c_void_p
        # pdll.st_calc_transform_matrices(pointer(a_euler), pointer(rreftoloc), pointer(rloctoref))

        # # reshape
        # a_rreftoloc = numpy.zeros((3,3))
        # a_rloctoref = numpy.zeros((3,3))
        # for i in range(3):
        #     for j in range(3):
        #         a_rreftoloc[i,j] = rreftoloc[i*3+j]
        #         a_rloctoref[i,j] = rloctoref[i*3+j]

        # return {'rreftoloc':a_rreftoloc, 'rloctoref':a_rloctoref}
        pass

    def util_matrix_transpose(self, m):
        # """
        # Calculate matrix transpose

        # Parameters
        # ----------
        # m : [[float,]*3]*3
        #     Square matrix 3x3 to be transposed.

        # Returns
        # ----------
        # [[float,]*3]*3
        #     Square matrix 3x3, transpose of 'm'
        # """

        # assert type(m) == type(numpy.array([]))
        # return m.T
        pass

    def util_rotation_arbitrary(self, theta, axis, axloc, pt):
        # """
        # Rotation of a point 'pt' about an arbitrary axis with direction 'axis' centered at point 'axloc'.
        # The point is rotated through 'theta' radians.

        # Parameters
        # ----------
        # theta : float
        #     Angle of rotation (rad)
        # axis : Point()
        #     Vector (x=i,y=j,z=k) indicating direction of axis for rotation
        # axloc : Point()
        #     Location of the axis origin
        # pt : Point()
        #     Location of the point that is to be rotated

        # Returns
        # -----------
        # Point
        #     Point after rotation
        # """

        # assert type(axis) == type(Point())
        # assert type(axloc) == type(Point())
        # assert type(pt) == type(Point())

        # a = axloc.x     #Point through which the axis passes
        # b = axloc.y
        # c = axloc.z
        # x = pt.x        #Point that we're rotating
        # y = pt.y
        # z = pt.z
        # u = axis.x        #Direction of the axis that we're rotating about
        # v = axis.y
        # w = axis.z


        # sinth = math.sin(theta)
        # costh = math.cos(theta)

        # fin = Point()

        # fin.x = (a*(v*v+w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1.-costh) + x*costh + (-c*v + b*w - w*y + v*z)*sinth
        # fin.y = (b*(u*u+w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1.-costh) + y*costh + (c*u - a*w + w*x - u*z)*sinth
        # fin.z = (c*(u*u+v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1.-costh) + z*costh + (-b*u + a*v - v*x + u*y)*sinth

        # return fin
        pass

    def util_calc_unitvect(self, vect):
        # """
        # Scales a vector to have total magnitude of 1

        # Parameters
        # ----------
        # vect : list | Point
        #     list or Point containing the vector

        # Returns
        # ----------
        # list | Point
        #     Unitized vector of type list or Point, depending on input type
        # """
        # if type(vect) == list:
        #     v = Point()
        #     v.x = vect[0]
        #     v.y = vect[1]
        #     v.z = vect[2]
        # else:
        #     v = vect

        # vect_mag = math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
        # nvect = Point()
        # nvect.x = v.x / vect_mag
        # nvect.y = v.y / vect_mag
        # nvect.z = v.z / vect_mag

        # if type(vect) == list:
        #     return [nvect.x, nvect.y, nvect.z]
        # else:
        #     return nvect
        pass

    def util_calc_zrot_azel(self, vect) -> float:
        # """
        # Compute the z-rotation of a vector, assuming the vector's deviation from (0,0,1)
        # has been realized using azimuth-elevation transforms.

        # Parameters
        # ----------
        # vect : (list OR Point)
        #     i,j,k components of a vector

        # Returns
        # ----------
        # float
        #     Computed z-rotation (degrees)
        # """
        # if isinstance(vect, List):
        #     vect_i, vect_j, vect_k = vect
        # elif isinstance(vect, Point):
        #     vect_i = vect.x
        #     vect_j = vect.y
        #     vect_k = vect.z
        # else:
        #     raise TypeError("Function expects 'vect' of type List or Point")

        # az = math.atan2(vect_i,vect_j)
        # az = (az + 2.*math.pi) if az < 0. else az

        # el = math.asin(vect_k)

        # #Calculate Euler angles
        # alpha = math.atan2(vect_i, vect_k);        #Rotation about the Y axis
        # bsign = 1 if vect_j > 0. else -1
        # beta = -bsign*math.acos( ( math.pow(vect_i,2) + math.pow(vect_k,2) )/
        #                 max(math.sqrt(math.pow(vect_i,2) + math.pow(vect_k,2)), 1.e-8) )    #Rotation about the modified X axis

        # #Calculate the modified axis vector
        # modax = Point(math.cos(alpha), 0., - math.sin(alpha))

        # #Rotation references - axis point. Set as origin
        # axpos = Point(0., 0., 0.)
        # #sp_point to rotate. lower edge of heliostat
        # pbase = Point(0., -1., 0.)

        # #Rotated point
        # protv = self.util_rotation_arbitrary(beta, modax, axpos, pbase).unitize()

        # #Azimuth/elevation reference vector (vector normal to where the base of the heliostat should be)
        # azelref = Point()
        # azelref.x = math.sin(az)*math.sin(el)
        # azelref.y = math.cos(az)*math.sin(el)
        # azelref.z = -math.cos(el)

        # # the sign of the rotation angle is determined by whether the 'k' component of the cross product
        # # vector is positive or negative.
        # cp = Point()
        # cp.x = protv.y*azelref.z - protv.z*azelref.y
        # cp.y = protv.z*azelref.x - protv.x*azelref.z
        # cp.z = protv.x*azelref.y - protv.y*azelref.x

        # gamma = math.asin( cp.radius() )
        # gsign = (1 if cp.z > 0. else -1.) * (1 if vect_j > 0. else -1.)

        # return gamma * gsign * 180./math.pi
        pass

    def write_soltrace_input_file(self, path : str):
        # """
        # Write a SolTrace input file (.stinput) based on the currently created API objects. This file is written
        # using the objects and data in the PySolTrace instance, not necessarily on what has been created in the
        # coretrace 'context' data space. The 'context' may not match the PySolTrace instance if not all 'Create()'
        # methods have been called.

        # Parameters
        # ==========
        # path : str
        #     Path and file name to be used to write the resulting .stinput file

        # Returns
        # =======
        # None
        # """
        # with open(path, 'w') as fout:

        #     # Header
        #     dt = datetime.now()
        #     fout.write(
        #         "# SOLTRACE VERSION 2012.7.6 INPUT FILE -- GENERATED BY soltrace-api v{:s} on {:02d}/{:02d}/{:04d} at {:02d}:{:02d}:{:02d}\n".format(
        #             "v000", dt.day, dt.month, dt.year, dt.hour, dt.minute, dt.second)
        #     )

        #     #------------------- Sun shape
        #     fout.write( "SUN\tPTSRC\t{:d}\tSHAPE\t{:1s}\tSIGMA\t{:f}\tHALFWIDTH\t{:f}\n".format( 0, self.sun.shape, self.sun.sigma, self.sun.sigma) )
        #     fout.write( "XYZ\t{:f}\t{:f}\t{:f}\tUSELDH\t{:d}\tLDH\t{:f}\t{:f}\t{:f}\n".format( self.sun.position.x, self.sun.position.y, self.sun.position.z, 0, 0., 0., 0.) )
        #     if( self.sun.shape == 'd' ):
        #         np = len(self.sun.user_intensity_table)
        #         fout.write( "USER SHAPE DATA\t{:d}\n".format( np) )
        #         for i in range(np):
        #             fout.write( "{:f}\t{:f}\n".format( *self.sun.user_intensity_table[i] ) )
        #     else:
        #         fout.write( "USER SHAPE DATA\t{:d}\n".format( 0) )

        #     #------------------- Optics list
        #     fout.write( "OPTICS LIST COUNT\t{:d}\n".format( len(self.optics) ) )
        #     for optics in self.optics:
        #         fout.write( "OPTICAL PAIR\t{:s}\n".format( optics.name ) )
        #         for opt in [optics.front, optics.back]:
        #             fout.write( "OPTICAL\t{:1s}\t{:d}\t{:d}\t{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\n".format(
        #                     opt.dist_type,3, 1, 4, opt.reflectivity, opt.transmissivity, opt.slope_error, opt.spec_error , opt.refraction_real, 0., 0., 0., 0., 0. )
        #                 )

        #     #------------------- loop through Stage list
        #     fout.write( "STAGE LIST COUNT\t{:d}\n".format( len(self.stages)) )
        #     for stage in self.stages:
        #         fout.write(
        #             "STAGE\tXYZ\t{:f}\t{:f}\t{:f}\tAIM\t{:f}\t{:f}\t{:f}\tZROT\t{:f}\tVIRTUAL\t{:d}\tMULTIHIT\t{:d}\tELEMENTS\t{:d}\tTRACETHROUGH\t{:d}\n".format(
        #                 stage.position.x, stage.position.y,stage.position.z, stage.aim.x, stage.aim.y, stage.aim.z, stage.zrot, int(stage.is_virtual), int(stage.is_multihit), len(stage.elements), int(stage.is_tracethrough))
        #                 )

        #         fout.write( "{:s}\n".format( stage.name ) )

        #         #------------------- loop through element list
        #         for el in stage.elements:
        #             # format string
        #             fmt = '\t'.join(['{:d}'] + ['{:f}']*7 + (['{:1s}'] + ['{:f}']*8)*2 + ['{:s}', '{:s}', '{:d}', '{:s}\n'])
        #             fout.write(fmt.format(
        #                 int(el.enabled),
        #                 el.position.x, el.position.y, el.position.z,  #origin
        #                 el.aim.x, el.aim.y, el.aim.z,  #aim
        #                 el.zrot,
        #                 el.aperture,
        #                 # Ap_A, Ap_B, Ap_C, Ap_D, Ap_E, Ap_F, Ap_G, Ap_H,
        #                 *el.aperture_params,
        #                 el.surface,
        #                 # Su_A, Su_B, Su_C, Su_D, Su_E, Su_F, Su_G, Su_H,
        #                 *el.surface_params,
        #                 "" if not el.surface_file else el.surface_file, #Surface geometry file
        #                 el.optic.name, el.interaction,
        #                 "" ) )
        # return
        pass