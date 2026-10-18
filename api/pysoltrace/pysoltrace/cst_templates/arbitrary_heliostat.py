from __future__ import annotations
# import subprocess, os
import numpy as np
# import pvlib
# import pandas as pd
# import orjson
from typing import Literal
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from functools import reduce
# from datetime import datetime, timedelta
# from timer import timer
import pysoltrace.soltrace_json as stJSON
import pysoltrace.soltrace_structs as _STS
from pysoltrace.math_utils import arbitrary_rotation, zrot_from_azel
from pysoltrace import api


class Facet():
    def __init__(self, name: str, opt_id: int | str,
                 pos: np.ndarray, canting: np.ndarray,
                 ap_type: str, ap: np.ndarray, 
                 surf_type: str, surf: np.ndarray, 
                 canting_err: float = 0.0):
        self.name = name
        self.opt_id = opt_id
        self.position = pos # local to heliostat centroid
        self.canting_err = canting_err
        self.aperture_type = ap_type
        self.aperture = ap
        self.surface_type = surf_type
        self.surface = surf

        self.canting_err_x = np.random.normal(0, canting_err)
        self.canting_err_y = np.random.normal(0, canting_err)

        # canting direction relative to heliostat in face up orientation
        # ie relative to (0, 0, 1)
        self.canting = np.array([
            canting[0] + self.canting_err_x,
            canting[1] + self.canting_err_y,
            canting[2],
        ])

    def __repr__(self):
        string = self.name + '\n' \
                + 'Position:  ' + str(self.position) + '\n' \
                + 'Canting:   ' + str(self.canting) + '\n' \
                + 'Aperture: ' + str(self.aperture_type) + ' ' + str(self.aperture) + '\n' \
                + 'Surface:   ' + str(self.surface_type) + ' ' + str(self.surface) + '\n'
        return string

    @property
    def xyzuvw(self): return np.append(self.position, self.canting)
    
    @property
    def x(self): return self.position[0]

    @property
    def y(self): return self.position[1]

    @property
    def z(self): return self.position[2]

    @property
    def u(self): return self.canting[0]

    @property
    def v(self): return self.canting[1]

    @property
    def w(self): return self.canting[2]

    def aim(self, az, el, el_axis, pos_h):
        Z = np.array([0.0,0.0,1.0])
        pos_g = pos_h + self.position # global position, face up orientation
        aim_g = pos_g + self.canting # global aim point, face up orientation

        # rotate global postion around z axis then around the axis defined by el_axis
        # TODO: write np.array version of util_rotation_arbitrary
        position = arbitrary_rotation(-az, Z, pos_h, pos_g)
        position = arbitrary_rotation(np.pi/2 - el, el_axis, pos_h, position)

        # rotate global aim point around z axis then around the axis defined by el_axis
        aim = arbitrary_rotation(-az, Z, pos_h, aim_g)
        aim = arbitrary_rotation(np.pi/2 - el, el_axis, pos_h, aim)

        # normal vector of facet aim point relative to facet position
        norm = (aim - position).unitize()
        zrot = zrot_from_azel(norm)

        return position.as_list(), aim.as_list(), float(zrot)
    
    def to_dict(self, heliostat_name: str, optical_registry: stJSON.OpticalPropertyRegistry,
                azimuth: float, elevation: float, 
                elevation_axis: np.ndarray, pos_h: np.ndarray,
                group: int = -1) -> dict:
        """
        Convert the facet state to a dictionary representation.
        """

        pos, aim, zrot = self.aim(azimuth, elevation, elevation_axis, pos_h)

        return stJSON.element_json(
            f'{heliostat_name} {self.name}',
            pos,
            aim,
            zrot,
            stJSON.APERTURE_JSON[self.aperture_type](self.aperture, self.aperture),
            stJSON.SURFACE_JSON[self.surface_type](*self.surface),
            self.opt_id,
            optical_registry,
            group
        )

    def add(self, stapi: api, heliostat_name: str,
            azimuth: float, elevation: float, 
            elevation_axis: np.ndarray, pos_h: np.ndarray,
            group: int = -1):
        pos, aim, zrot = self.aim(azimuth, elevation, elevation_axis, pos_h)

        el_args = _STS.element([*map(float, pos)],
                               [*map(float, aim)],
                               zrot, 
                               True, False,
                               group)
        geom_args = _STS.geometry(self.aperture.tolist(),
                                  self.surface.tolist(),
                                  self.aperture_type,
                                  self.surface_type)
        el_id = stapi.data.element.add(el_args,
                                       self.opt_id,
                                       geom_args)
        return (el_id, el_args)

class Heliostat():

    def __init__(
            self, 
            name: str,
            pos: np.ndarray[Literal[3], float],
            target: np.ndarray[Literal[3], float],
            sun: np.ndarray[Literal[3], float],
            tracking_err: float = 0.0,
            group: int = -1
        ):
        _pos = pos
        _tar = target
        _sun = sun
        if not isinstance(_pos, np.ndarray): _pos = np.array(_pos)
        if not isinstance(_tar, np.ndarray): _tar = np.array(_tar)
        if not isinstance(_sun, np.ndarray): _sun = np.array(_sun)
        
        self.name = name
        self.position = _pos        # in global coordiante system
        self.facets_local = []
        self._num_facets = 0
        self._current_facet = 0
        self.target = _tar - _pos # in local coordiante system
        self.sun = _sun             # treating sun as completely parallel
        self.tracking_err = tracking_err
        self.group = group

        # possible aim point error strategy
        # self.aim = Point().set_cartesian(*self.tracking).rotate_spherical(
        #     np.random.uniform(0, 2* np.pi),     # random azimuth to deviate from
        #     np.random.normal(0, tracking_err)   # random amount of deviation in altitude
        # ).cartesian                             # average angular deviation = tracking_err

        # from NSTTF_Solar 1 flux target_multifacet.lk
        self.tracking_err_x = np.random.normal(0, tracking_err)
        self.tracking_err_y = np.random.normal(0, tracking_err)
        # aim direction relative to heliostat centroid
        self._set_tracking()

    def __repr__(self):
        string = self.name + '\n' \
            + 'Position: ' + str(self.position) + '\n' \
            + 'Facets:   ' + str(len(self.facets_local)) + '\n'
        return string

    def __iter__(self) -> Heliostat:
        # reset iteration each time
        self._current_facet = 0 
        return self

    def __next__(self) -> Facet:
        if self._current_facet < self._num_facets:
            self._current_facet += 1
            return self.facets_local[self._current_facet - 1]
        raise StopIteration
    
    def __getitem__(self, i) -> Facet:
        return self.facets_local[i]

    def __len__(self):
        return self._num_facets

    @property
    def xyzuvw_local(self): return np.array([f.xyzuvw for f in self.facets_local])
    
    @property
    def xyz_local(self): return np.array([f.position for f in self.facets_local])
    
    @property
    def x(self): return np.array([f.x for f in self.facets_local])

    @property
    def y(self): return np.array([f.y for f in self.facets_local])

    @property
    def z(self): return np.array([f.z for f in self.facets_local])
    
    @property
    def uvw_local(self): return np.array([f.canting for f in self.facets_local])

    @property
    def u(self): return np.array([f.u for f in self.facets_local])

    @property
    def v(self): return np.array([f.v for f in self.facets_local])

    @property
    def w(self): return np.array([f.w for f in self.facets_local])

    @property
    def area(self): return sum([f.area for f in self])
    
    def _set_tracking(self) -> None:
        """
        Compute the heliostat tracking and aim direction directions.

        The tracking direction is computed relative to the heliostat. This
        represents the ideal direction the heliostat should point.

        The aim direction represents the actual pointing direction of the
        heliostat, incorporating tracking error offsets in the x and y directions.

        Updates
        -------
        self.tracking : numpy.ndarray
            Unit vector representing the ideal tracking direction.
        self.aim : numpy.ndarray
            Aim direction vector including tracking error deviations.

        Returns
        -------
        None
        """

        self.tracking = self.target / np.linalg.norm(self.target) + self.sun / np.linalg.norm(self.sun)
        self.tracking /= np.linalg.norm(self.tracking)
        self.aim = np.array([
            self.tracking[0] + self.tracking_err_x,
            self.tracking[1] + self.tracking_err_y,
            self.tracking[2],
        ])

    def set_optical_property(self, opt_id: int):
        for f in self: f.opt_id = opt_id

    def set_target_global(self, new_target: np.ndarray[Literal[3], float]):
        _tar = new_target
        if not isinstance(_tar, np.ndarray): _tar = np.array(_tar)

        self.target = _tar - self.position # in local coordiante system

    def reaim(self, new_sun: np.ndarray) -> None: 
        """
        Update the heliostat aiming direction based on a new sun direction.

        Parameters
        ----------
        new_sun : numpy.ndarray
            Updated sun position.

        Updates
        -------
        self.sun : numpy.ndarray
            Updated sun position.

        Returns
        -------
        None
        """
        _sun = new_sun
        if not isinstance(_sun, np.ndarray): _sun = np.array(_sun)

        self.sun = _sun
        self._set_tracking()
    
    def add_facets_local(self, canting_info: np.ndarray, opt_id: int | str,
                         ap_type: str, ap: np.ndarray, 
                         surf_type: str, surf: np.ndarray, 
                         canting_err: float = 0.0) -> None:
        """
        Expects canting_info to be np.array with rows corresponding to each facet
        and columns in the order id, x, y, z, aim x, aim y, aim z.
        """
        _ap = ap
        _surf = surf
        if not isinstance(_ap, np.ndarray): _ap = np.array(_ap)
        if not isinstance(_surf, np.ndarray): _surf = np.array(_surf)

        # set up facets in local coords
        for i in range(0, canting_info.shape[0]):
            self.facets_local.append(
                Facet(
                    str(int(canting_info[i, 0])),
                    opt_id,
                    canting_info[i, 1:4],
                    canting_info[i, 4:],
                    ap_type, _ap,
                    surf_type, _surf,
                    canting_err,
                )
            )
        self._num_facets = len(self.facets_local)
    
    def save_to_json(self, optical_registry: stJSON.OpticalPropertyRegistry, g_or_l: str = 'g'):
        """
        Convert the heliostat state to a dictionary representation in JSON-compatible format.

        Parameters
        ----------
        g_or_l : str, optional
            Coordinate system specifier: 'g' for global or 'l' for local.

        Returns
        -------
        dict
            Dictionary containing the heliostat's configuration and state in JSON-compatible format.
        """

        if not (g_or_l == 'g' or g_or_l == 'l'): raise ValueError('chose local or global')
        if not len(self): raise ValueError('add facets to heliostat')
 
        az = np.arctan2(self.aim[0], self.aim[1]) - np.pi # [rad]
        el = np.arcsin(self.aim[2]) # [rad] elevation angle not zenith
        # print(f'\n{self.name}\taim: {self.aim}\taz: {az}\tel: {el}')

        elevation_axis = arbitrary_rotation(-az, np.array([0.0,0.0,1.0]), np.array([0.0,0.0,0.0]), np.array([1.0, 0.0, 0.0]))

        return [f.to_dict(self.name, optical_registry, az, el, elevation_axis, self.position, self.group) for f in self]

    def add(self, stapi: api):
        az = np.arctan2(self.aim[0], self.aim[1]) - np.pi # [rad]
        el = np.arcsin(self.aim[2]) # [rad] elevation angle not zenith

        elevation_axis = arbitrary_rotation(-az,
                                            np.array([0.0, 0.0, 1.0]),
                                            np.array([0.0, 0.0, 0.0]),
                                            np.array([1.0, 0.0, 0.0]))

        return [f.add(stapi, self.name,
                      az, el, elevation_axis,
                      self.position, self.group)
                for f in self]
