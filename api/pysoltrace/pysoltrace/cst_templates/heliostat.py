"""
@author: whamilton

Creates a heliostat template made of multiple mirror panels.
Heliostat template updates positioning based on sun position.
"""

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..', '..'))

from pysoltrace import PySolTrace, Point
from parabolic_trough import compute_euler_angles, sun_vector, trace, sun_position, meinel_clearsky

import matplotlib.pyplot as plt
import numpy as np

from typing import Union


class heliostat:
    
    def __init__(self, 
                 position: Point, 
                 aperture_size: tuple, 
                 number_panels: tuple, 
                 gaps: tuple, 
                 focal_length: Union[float,tuple],
                 aim_point: Point,
                 onaxis_canting_distance: float = np.nan, 
                 offaxis_canting_sun_position: tuple = (np.nan, np.nan)):
        """
        :param position: [m] Position of the heliostat
        :param aperture_size: [m] Size of the heliostat aperture (width, height)
        :param number_panels: [-] Number of panels (width, height)
        :param gaps: [m] Gaps between the panels (width, height)
        :param focal_length: [m] Focal length of mirror panels
        :param aim_point: [m] Aimpoint of heliostat on the receiver in global coordinates (x, y, z)
        :param onaxis_canting_distance (optional): [m] Distance from the heliostat to cant facets
        :param offaxis_canting_sun_position (optional): [deg] Solar position at canting (azimuth, zenith)
        """
        self.position = position

        if len(aperture_size) != 2:
            raise ValueError("Aperture size should be a tuple with 2 elements")
        if any(ap_dim < 0.0 for ap_dim in aperture_size):
            raise ValueError("Aperture dimensions must be positive values")
        self.aperture_size = aperture_size

        if len(number_panels) != 2:
            raise ValueError("Number of panels should be a tuple with 2 elements")
        if any(panels < 0 for panels in number_panels):
            raise ValueError("Number of panels must be positive values")
        self.number_panels = number_panels

        if len(gaps) != 2:
            raise ValueError("Gaps should be a tuple with 2 elements")
        if any(gap < 0.0 for gap in gaps):
            raise ValueError("Gaps must be positive values")
        self.gaps = gaps

        if isinstance(focal_length, float):            
            if focal_length < 0.0:
                raise ValueError("Focal length must be greater than or equal to zero")
            self.focal_length = (focal_length, focal_length)
        if isinstance(focal_length, tuple):
            if any(flen < 0.0 for flen in focal_length):
                raise ValueError("Focal lengths must be greater than or equal to zero")
            self.focal_length = focal_length

        self.aim_point = aim_point

        self.onaxis_canting_distance = np.nan
        self.offaxis_canting_sun_position = (np.nan, np.nan)
        self.canting_method = 'no_canting'  # We could use an enum type
        if not np.isnan(onaxis_canting_distance):
            if not any(np.isnan(np.array(offaxis_canting_sun_position, dtype=float))):
                raise ValueError("Only specify either on-axis canting distance or off-axis canting sun position not both")
            if onaxis_canting_distance <= 0.0:
                raise ValueError("On-axis canting distance must be a positive value")
            self.onaxis_canting_distance = onaxis_canting_distance
            self.canting_method = 'on-axis'
        elif not any(np.isnan(np.array(offaxis_canting_sun_position, dtype=float))):
            if offaxis_canting_sun_position[0] < 0 or offaxis_canting_sun_position[0] > 360:
                raise ValueError("Off-axis canting sun azimuth angle must be between 0 and 360 degrees")
            if offaxis_canting_sun_position[1] < 0 or offaxis_canting_sun_position[1] > 90:
                raise ValueError("Off-axis canting sun zenith angle must be between 0 and 90 degrees")
            self.offaxis_canting_sun_position = offaxis_canting_sun_position
            self.canting_method = 'off-axis'
        
        # Calculated values
        self.heliostat_area = 0.0
        self.tracking_azimuth = 0.0
        self.tracking_elevation = 90.0              # Pointing straight up
        self.normal = Point(0.0, 0.0, 1.0)
        self.elevation_axis = Point(1.0, 0.0, 0.0)

        self.efficiencies = dict()
        self.power2receiver = dict()

        self.receiver = PySolTrace.Stage.Element
        self.facets = list()

        # For testing only, relative positions and aimpoints
        num_el = number_panels[0] * number_panels[1]
        self.rel_positions = np.zeros((num_el, num_el))
        self.rel_aimpoints = np.zeros((num_el, num_el))

    def get_elements(self):
        return self.facets
    
    def set_receiver(self, receiver: PySolTrace.Stage.Element):
        self.receiver = receiver

    def create_geometry(self, PT:PySolTrace, facet_optics: PySolTrace.Optics):
        """
        Creates the elements of the heliostat in the PySolTrace object (PT). 
        NOTE: This assumes the stage has an origin of (0.0, 0.0, 0.0), Aimpoint of (0.0, 0.0, 1.0), and a z-rotation of 0.0 degrees
        
        :param PT: PySolTrace object
        :param facet_optics: facet surface optics (mirror)
        """

        # Determine the aperture bounds of each panel
        panel_width = (self.aperture_size[0] - self.gaps[0] * (self.number_panels[0] - 1)) / self.number_panels[0]
        panel_height = (self.aperture_size[1] - self.gaps[1] * (self.number_panels[1] - 1)) / self.number_panels[1]

        if self.canting_method == 'off-axis':
            # Off-axis canting calculations that are not dependent on specific panel
            sun_vec = sun_vector(self.offaxis_canting_sun_position[0], 90.0 - self.offaxis_canting_sun_position[1]) # TODO: zenith and elevation is mixed up.
            helio2aim = (self.aim_point - self.position).unitize()
            tracking = (helio2aim + sun_vec).unitize()                  # Heliostat tracking vector (i.e., normal)

            # Calculate the tracking azimuth/zenith based on the tracking vector
            track_az = np.atan2(tracking.x, tracking.y)
            track_el = np.asin(tracking.z)

            # Calculate the panel's actual x-y-z location w/r/t the global coordinates
            delta_azimuth = track_az - np.pi
            delta_elevation = track_el - np.pi/2
            elevation_axis = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), Point(), Point(1.0, 0.0, 0.0))

        # Build heliostat about the origin then translate to position
        stage = PT.stages[0]        # NOTE: assuming everything is in the first stage (removing stages)
        # Create mirror panels
        self.facets.clear()
        self.heliostat_area = 0.0
        panel_y = - self.aperture_size[1] / 2 + panel_height / 2
        for _ in range(self.number_panels[1]):          # height
            panel_x = - self.aperture_size[0] / 2 + panel_width / 2
            for _ in range(self.number_panels[0]):      # width
                element = stage.add_element()

                if self.canting_method == 'no_canting':
                    element.position = Point(panel_x, panel_y, 0.0)
                    element.aim = Point(panel_x, panel_y, 100)  # point straight up
                elif self.canting_method == 'on-axis':
                    c = 1 / (2. * self.onaxis_canting_distance)
                    z = 0.5 * c *(panel_x**2 + panel_y**2)
                    element.position = Point(panel_x, panel_y, z)     # NOTE: Z-position shift the panels slightly to follow canting curvature
                    # element.position = Point(panel_x, panel_y, 0.0)   # This option center points all fall on the same plane.
                    element.aim = Point(0.0, 0.0, 2 * self.onaxis_canting_distance)
                    # element.aim = Point(0.0, 0.0, self.onaxis_canting_distance)
                elif self.canting_method == 'off-axis':
                    element.position = Point(panel_x, panel_y, 0.0)   # This option panel center points all fall on the same plane.

                    # Calculate the panel's position within the global coordinates
                    p_pos = element.position
                    p_pos = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), Point(), p_pos)
                    p_pos = PT.util_rotation_arbitrary(-delta_elevation, elevation_axis, Point(), p_pos)
                    p_pos_global = p_pos + self.position
                    
                    # Determine the vector from the panel centroid to the aim point
                    pref = (self.aim_point - p_pos_global).unitize()
                    pnorm = (pref + sun_vec).unitize()      # Panel normal

                    # Translate back to stow position
                    pnorm = PT.util_rotation_arbitrary(delta_elevation, elevation_axis, Point(), pnorm)
                    pnorm = PT.util_rotation_arbitrary(delta_azimuth, Point(0.0,0.0,1.0), Point(), pnorm)

                    # Scale aim to target and translate to panel position
                    scale = 2.0 * (self.aim_point - p_pos_global).radius()
                    element.aim = pnorm * scale + element.position
                else:
                    raise ValueError("Unexpected canting method value")

                element.zrot = 0.0
                element.aperture_rectangle(panel_width, panel_height)
                self.heliostat_area += panel_width * panel_height

                if self.focal_length[0] == 0.0:
                    element.surface_flat()
                else:
                    element.surface_parabolic(self.focal_length[0], self.focal_length[1])
                
                element.optic = facet_optics
                element.enabled = True

                self.facets.append(element)
                # Increment the x position
                panel_x += panel_width + self.gaps[0]
                
            # Increment the y position
            panel_y += panel_height + self.gaps[1]

        # NOTE: for TESTING only, calculate relative positions and aimpoints
        for i, ref_el in enumerate(self.get_elements()):
            for j, el in enumerate(self.get_elements()):
                self.rel_positions[i,j] = (el.position - ref_el.position).radius()
                self.rel_aimpoints[i,j] = (el.aim - ref_el.aim).radius()

        # Aiming heliostat at aim point but in stow (straight up)
        helio2aim = self.aim_point - self.position
        y_axis = Point(-helio2aim.x, -helio2aim.y, 0.0).unitize()
        assert abs(np.dot(self.normal.as_list(), y_axis.as_list())) <= 1.e-6
        self.tracking_azimuth = np.atan2(helio2aim.x, helio2aim.y) * 180.0 / np.pi

        elevation_axis = np.cross(y_axis.as_list(), self.normal.as_list())
        self.elevation_axis = Point(elevation_axis[0], elevation_axis[1], elevation_axis[2])
        assert abs(np.dot(self.elevation_axis.as_list(), y_axis.as_list())) <= 1.e-6

        # Rotation matrix from stage to local coordinates
        rotation_matrix = np.array([[self.elevation_axis.x, y_axis.x, self.normal.x], 
                                    [self.elevation_axis.y, y_axis.y, self.normal.y],
                                    [self.elevation_axis.z, y_axis.z, self.normal.z]])
        
        assert abs(np.linalg.det(rotation_matrix) - 1) <= 1.e-8
        euler_angles = compute_euler_angles(np.linalg.matrix_transpose(rotation_matrix))   # NOTE: Stage coordinates are assumed to be Identity matrix
        # NOTE: Transpose of a rotation matrix is equal to its inversve.
        #print("Euler angles (initial) {}".format(euler_angles))

        # Updating element positions, aimpoints, and z-rotations based on aimpoint and positon
        azimuth_rads = (180.0 + self.tracking_azimuth) * np.pi / 180.0
        for element in self.get_elements(): 
            # negative azimuth due to right-hand rule
            element.position = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(), element.position)
            element.aim = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(), element.aim)

            element.zrot = euler_angles[2]

            # Translate to specific location
            element.position += self.position
            element.aim += self.position            

    def update_geometry(self, PT:PySolTrace, azimuth: float, elevation: float):
        """
        Update heliostat geometry based on solar position (azimuth and elevation)

        :param PT: PySolTrace instance
        :param azimuth: [deg] Solar azimuth angle
        :param elevation: [deg] Solar elevation angle
        """
        # Calculate and update solar position
        sun_vec = sun_vector(azimuth, elevation)
        PT.sun.position = sun_vec * 1000.0

        # Unit vect from heliostat to aim point
        helio2aim = (self.aim_point - self.position).unitize()
        self.normal = (helio2aim + sun_vec).unitize()

        # Calculate updated tracking azimuth and elevation angles
        prev_tracking_azimuth = self.tracking_azimuth
        self.tracking_azimuth = np.atan2(self.normal.x, self.normal.y) * 180.0 / np.pi
        delta_azimuth = (self.tracking_azimuth - prev_tracking_azimuth) * np.pi / 180.0

        prev_tracking_elevation = self.tracking_elevation
        self.tracking_elevation = np.asin(self.normal.z) * 180.0 / np.pi
        delta_elevation = (self.tracking_elevation - prev_tracking_elevation) * np.pi / 180.0

        # Update elevation axis
        self.elevation_axis = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), Point(), self.elevation_axis)

        # Updating element positions, aimpoints, and z-rotations based on tilt and azimuth
        for element in self.get_elements():
            element.position = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, element.position)
            element.position = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, element.position)

            element.aim = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, element.aim)
            element.aim = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, element.aim)

            # Updating each element's z-rotation
            norm = (element.aim - element.position).unitize()
            y_axis = np.cross(norm.as_list(), self.elevation_axis.as_list())
            y_axis = Point(y_axis[0], y_axis[1], y_axis[2]).unitize()

            element_el_axis = np.cross(y_axis.as_list(), norm.as_list())
            element_el_axis = Point(element_el_axis[0], element_el_axis[1], element_el_axis[2])

            rotation_matrix = np.array([[element_el_axis.x, y_axis.x, norm.x], 
                                        [element_el_axis.y, y_axis.y, norm.y],
                                        [element_el_axis.z, y_axis.z, norm.z]])
            
            assert abs(np.linalg.det(rotation_matrix) - 1) <= 1.e-8
            euler_angles = compute_euler_angles(np.linalg.matrix_transpose(rotation_matrix)) 
            element.zrot = euler_angles[2]  

        # Testing relative positions and aimpoints
        for i, ref_el in enumerate(self.get_elements()):
            for j, el in enumerate(self.get_elements()):
                rel_pos = (el.position - ref_el.position).radius()
                assert abs(self.rel_positions[i,j] - rel_pos) <= 1e-6

                rel_aim = (el.aim - ref_el.aim).radius()
                assert abs(self.rel_aimpoints[i,j] - rel_aim) <= 1e-6

    def calculate_receiver_power(self, PT: PySolTrace):
        """
        Calculates the power absorbed by the receiver tubes

        :param PT: PySolTrace instance

        :return total absorbed power: [kWt]
        """
        rayData = PT.raydata
        if rayData.empty:
            raise LookupError("RayData is empty when trying to calculate receiver power")
        
        # Filter ray data by mirror element
        power = 0.0
        for el in self.get_elements():
            hit_rays = rayData[rayData['element'] == el.id+1]['number']
            hitRayData = rayData[rayData['number'].isin(hit_rays)]
            power += len(hitRayData[hitRayData['element'] == -(self.receiver.id+1)]) * PT.powerperray / 1.e3
        
        self.power2receiver = power
        return power    

    def calculate_efficiencies(self, PT: PySolTrace):
        """
        Calculates heliostat efficiencies including cosine (includes shadowing), reflect (includes attenuation), 
        intercept, blocking, absorption (on receiver), heliostat total efficiency (reflect * intercept * blocking), 
        and total optical efficiency (cosine * reflect * intercept * blocking * absorption)

        :param PT: PySolTrace instance

        :return efficiencies: dict [-] with keys (cosine, reflect, intercept, blocking, absorption, helio_total, total)
        """
        rayData = PT.raydata
        if rayData.empty:
            raise LookupError("RayData is empty when trying to calculate receiver power")
        
        # Filter ray data to only include rays directly from sun
        sun_vec = PT.sun.position.unitize()
        sunRayData = rayData[(abs(sun_vec.x + rayData['cos_x']) <= 1e-6)           # This filter is typically sufficient 
                          & (abs(sun_vec.y + rayData['cos_y']) <= 1e-6)
                          & (abs(sun_vec.z + rayData['cos_z']) <= 1e-6)]

        Abox = (PT.sunstats['xmax']-PT.sunstats['xmin'])*(PT.sunstats['ymax'] - PT.sunstats['ymin']) 
        nsunrays = PT.sunstats['nsunrays']
        N_hin = 0
        N_habs = 0
        N_hblock = 0
        N_hout = 0
        N_rin = 0
        N_rabs = 0
        for el in self.get_elements():
            N_hin += len(sunRayData[sunRayData['element'].abs() == el.id+1]['number'].unique())
            N_habs += len(sunRayData[sunRayData['element'] == -(el.id+1)]['number'].unique())

            # Rays that direct hit the heliostat element
            raysHitEl = sunRayData[sunRayData['element'].abs() == el.id+1]['number']
            elRays = rayData[rayData['number'].isin(raysHitEl)]
            N_rin += len(elRays[elRays['element'].abs() == self.receiver.id+1]['number'].unique())
            N_rabs += len(elRays[elRays['element'] == -(self.receiver.id+1)]['number'].unique())

            # Rays that are blocked - remove rays that hit or miss the receiver
            raysHitOrMissRec = elRays[(elRays['element'].abs() == self.receiver.id+1)
                                         | (elRays['element'] == 0)]['number']
            raysDontHitOrMissRec = elRays[~elRays['number'].isin(raysHitOrMissRec)]
            blockRaysNums = raysDontHitOrMissRec[raysDontHitOrMissRec['element'].abs() != el.id+1]['number']
            blockRays = elRays[elRays['number'].isin(blockRaysNums)]
            N_hblock += len(blockRays['number'].unique())

        N_hout = N_hin - N_habs - N_hblock

        self.efficiencies = dict()
        self.efficiencies['cosine'] = N_hin / nsunrays * Abox / self.heliostat_area
        self.efficiencies['reflect'] = (N_hin - N_habs) / N_hin
        self.efficiencies['intercept'] = N_rin / N_hout
        self.efficiencies['blocking'] = 1. - N_hblock / (N_hin - N_habs)
        self.efficiencies['absorption'] = (N_rabs / N_rin)
        self.efficiencies['helio_total'] = (N_rin / N_hin)
        self.efficiencies['total'] = (self.efficiencies['cosine'] 
                                      * self.efficiencies['reflect'] 
                                      * self.efficiencies['intercept'] 
                                      * self.efficiencies['blocking'] 
                                      * self.efficiencies['absorption'])
        
        assert abs(self.efficiencies['helio_total'] - (self.efficiencies['reflect'] 
                                                       * self.efficiencies['intercept'] 
                                                       * self.efficiencies['blocking'])) <= 1e-3
        # TODO: 
        # - Ray counting should be validated or confirmed
        # - Getting attenuation and shadowing seperated would be nice...
        #   - shadowing would require the ray to be "extended" to see if intersects with another heliostat. 
        #           This is a typical workflow for ray tracers capaturing shadowing effects
        #   - Attenuation would require a probabilistic approach based on the distance of a reflected ray. 
        # - max and mean flux on receiver
        # - image size on receiver
        return self.efficiencies


if __name__ == "__main__":
    run_basic_setup = False
    catch_plane_on = False

    solar_tracking_test = False
    multi_rows_heliostats = True
    plot_efficiencies = True

    run_canting_tests = True       # FIXME: off-axis canting center panel...


    PT = PySolTrace()
    PT.add_stage()

    # Create sun
    sun = PT.add_sun()
    sun.position = Point(0.0, 0.0, 100.0)
    sun.shape = 'd'
    sun.sigma = 4.65

    # Limb-darkened sun-shape
    pts = np.linspace(0, 4.65, 26)
    intensity = np.maximum(1.0 - 0.5138*((pts/4.65)**4), 0.0)  
    intensity[-1] = 0.0
    sun.user_intensity_table = [[pts[i], intensity[i]] for i in range(len(pts))]

    # Create optics
    mirror_optic = PT.add_optic('mirror')
    mirror_optic.front.reflectivity = 0.95
    mirror_optic.front.slope_error = 1.5
    mirror_optic.front.spec_error = 0.5
    mirror_optic.back.reflectivity = 0.0
    mirror_optic.back.slope_error = 100.
    mirror_optic.back.spec_error = 0.0

    target_optic = PT.add_optic('target')
    target_optic.front.reflectivity = 0.05
    target_optic.front.slope_error = 1e-5
    target_optic.front.spec_error = 1e-5
    target_optic.back.reflectivity = 0.05
    target_optic.back.slope_error = 1e-5
    target_optic.back.spec_error = 1e-5

    if run_basic_setup:
        # Set-up the heliostat
        position = Point(100.0, 66.0, 0.0)
        # position = Point(-40.0, 40.0, 0.0)
        #position = Point(0.0, 40.0, 0.0)
        aperture_size = (12., 12.)
        number_panels = (3, 4)
        gaps = (0.1, 0.1)
        focal_length = 156.06 #114.89 #115
        aim_point = Point(0.0, 0.0, 100.)

        onaxis_canting_distance = 156.06
        offaxis_canting_sun_position = (180.0, 45.0)

        # helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, offaxis_canting_sun_position = offaxis_canting_sun_position)
        helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, onaxis_canting_distance)
        # helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point)
        helio.create_geometry(PT, mirror_optic)

        # Create target
        stage = PT.stages[0]
        #stage = PT.add_stage()
        receiver = stage.add_element()
        receiver.position = aim_point
        receiver.aim = receiver.position + Point(0.0, 1.0, 0.0)
        receiver.zrot = 0.0
        receiver.aperture_rectangle(12.0, 12.0)
        receiver.surface_flat()
        receiver.optic = target_optic
        receiver.enabled = True
        helio.set_receiver(receiver)

        PT.write_soltrace_input_file('heliostat_create.stinput')

        helio.update_geometry(PT, 90.0, 25.0)
        PT.write_soltrace_input_file('heliostat_morning.stinput')
        trace(PT, nthreads=12)
        power = helio.calculate_receiver_power(PT)
        print("Power on receiver: {:.2f}".format(power))

        # flux_map = PT.bin_rays(receiver, 50, 50, True)
        # PT.plot_flux(receiver, 50, 50)

        helio.update_geometry(PT, 180.0, 45.0)
        PT.write_soltrace_input_file('heliostat_noon.stinput')
        trace(PT, nthreads=12)
        power = helio.calculate_receiver_power(PT)
        print("Power on receiver: {:.2f}".format(power))

        helio.update_geometry(PT, 225.0, 35.0)
        PT.write_soltrace_input_file('heliostat_afternoon.stinput')
        trace(PT, nthreads=12)
        power = helio.calculate_receiver_power(PT)
        print("Power on receiver: {:.2f}".format(power))


    if catch_plane_on:
        stage = PT.add_stage()
        element = stage.add_element()
        element.position = aim_point
        element.aim = element.position + Point(0.0, 1.0, 0.0)
        element.zrot = 0.0
        element.aperture_rectangle(200.0, 200.0)
        element.surface_flat()
        element.optic = target_optic
        element.enabled = True
        PT.write_soltrace_input_file('heliostat_two_stage.stinput')

    # Sun tracking - with multiple heliostat
    if solar_tracking_test:
        latitude = 34.8639          # Daggett, CA
        altitude = 610.0 / 1000.0   # [km]
        day_of_year = 80

        if multi_rows_heliostats:
            # This should be driven by a layout pattern, e.g., radial stagger
            delta_azi = 12.0
            row_radius = 75.0
            angle_bounds = (-60, 60)
            positions = list()
            for theta in [angle_bounds[0] + i*delta_azi for i in range(int((angle_bounds[1]-angle_bounds[0])/delta_azi)+1)]:
                positions.append(Point(row_radius * np.sin(theta * np.pi / 180), 
                                    row_radius * np.cos(theta * np.pi / 180),
                                    0.0))

            row_radius = 95.0
            angle_bounds = (-65, 65)
            for theta in [angle_bounds[0] + i*delta_azi for i in range(int((angle_bounds[1]-angle_bounds[0])/delta_azi)+1)]:
                positions.append(Point(row_radius * np.sin(theta * np.pi / 180), 
                                    row_radius * np.cos(theta * np.pi / 180),
                                    0.0))
                
            row_radius = 115.0
            angle_bounds = (-60, 60)
            for theta in [angle_bounds[0] + i*delta_azi for i in range(int((angle_bounds[1]-angle_bounds[0])/delta_azi)+1)]:
                positions.append(Point(row_radius * np.sin(theta * np.pi / 180), 
                                    row_radius * np.cos(theta * np.pi / 180),
                                    0.0))
            plot_eff_heliostats = [22, 27, 32]
        else:
            # 3 heliostats
            positions = [Point(-60, 60.0, 0.0), Point(0.0, 84.85, 0.0), Point(60.0, 60.0, 0.0)]
            plot_eff_heliostats = [0, 1, 2]

        # Heliostat set-up
        aperture_size = (12., 12.)
        number_panels = (3, 3)  #(3,4)
        gaps = (0.1, 0.1)
        focal_length = 138.0
        aim_point = Point(0.0, 0.0, 100.)

        onaxis_canting_distance = 138.0
        offaxis_canting_sun_position = (180.0, 45.0)

        PT.stages.clear()
        PT.add_stage()
        heliostats = list()
        for i, position in enumerate(positions):
            # helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point)
            helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, onaxis_canting_distance)
            # helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, offaxis_canting_sun_position = offaxis_canting_sun_position)
            helio.create_geometry(PT, mirror_optic)
            heliostats.append(helio)

        # PT.write_soltrace_input_file('heliostat_TS_create.stinput')

        # Adding receiver
        stage = PT.stages[0]    # Within the first stage
        #stage = PT.add_stage()
        receiver = stage.add_element()
        receiver.position = aim_point
        receiver.aim = receiver.position + Point(0.0, 1.0, 0.0)
        receiver.zrot = 0.0
        receiver.aperture_rectangle(12.0, 12.0)
        receiver.surface_flat()
        receiver.optic = target_optic
        receiver.enabled = True

        for helio in heliostats:
            helio.set_receiver(receiver)

        # Time-simulation
        rec_power = list()
        dni = list()

        helio_power2rec = list()
        cosine_eff = list()
        intercept_eff = list()
        block_eff = list()
        total_eff = list()
        for helio in heliostats:
            helio_power2rec.append(list())
            cosine_eff.append(list())
            intercept_eff.append(list())
            block_eff.append(list())
            total_eff.append(list())
        hours = [5 + h*0.5 for h in range(29)]
        #hours = [5 + h for h in range(15)]
        for h in hours:
            print("Simulating hour {:}...".format(h))
            azimuth, elevation = sun_position(latitude, day_of_year, h)    
            zenith = (90.0 - elevation) * np.pi / 180.0
            dni_hour = meinel_clearsky(day_of_year, zenith, altitude)   # Clear sky DNI values
            dni.append(dni_hour)

            for helio in heliostats:
                helio.update_geometry(PT, azimuth, elevation)

            if h in [7, 8, 9, 10, 11, 12]:
                PT.write_soltrace_input_file('heliostat_{}_hour.stinput'.format(h))

            if dni_hour > 0.0:
                trace(PT, nrays=1e6, dni=dni_hour, nthreads=12)
                
                # Calculate power from only reflected rays
                ray_counts = PT.raydata['number'].value_counts()
                multi_hits = ray_counts[ray_counts > 1].index
                multiHitsRayData = PT.raydata[PT.raydata['number'].isin(multi_hits)]
                tot_power = len(multiHitsRayData[multiHitsRayData['element'] == -(receiver.id+1)]) * PT.powerperray / 1.e3
                
                power = 0.0
                for i, helio in enumerate(heliostats):
                    helio_power = helio.calculate_receiver_power(PT)
                    helio_power2rec[i].append(helio_power)
                    power += helio_power

                    # Calculate efficiencies
                    efficiencies = helio.calculate_efficiencies(PT)
                    cosine_eff[i].append(efficiencies['cosine'])
                    intercept_eff[i].append(efficiencies['intercept'])
                    block_eff[i].append(efficiencies['blocking'])
                    total_eff[i].append(efficiencies['total'])

                assert abs(tot_power - power) <= 1e-6
                rec_power.append(power)
                assert power > 100.0    # Non-trivial power
            else:
                rec_power.append(0.0)
                for i, helio in enumerate(heliostats):
                    helio_power2rec[i].append(0.0)
                    cosine_eff[i].append(0.0)
                    intercept_eff[i].append(0.0)
                    block_eff[i].append(0.0)
                    total_eff[i].append(0.0)

        # Plot power and dni over time
        fig = plt.figure()
        plt.plot(hours, rec_power, 'k', label = "Total Power")
        # for i, helio in enumerate(heliostats):
        #     plt.plot(hours, helio_power2rec[i], label= "Heliostat = {:}".format(i))
        plt.ylabel('Receiver Power [kWt]')
        plt.xlabel('Hour of Day')
        # plt.legend()

        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(hours, dni, 'r')
        ax2.set_ylabel(r"Direct Normal Irradiance [W/m$^2$]", color = 'r')
        ax2.tick_params('y', colors = 'r')
        plt.tight_layout()
        plt.show()

        if plot_efficiencies:
            fig = plt.figure()
            for i in plot_eff_heliostats:
                plt.plot(hours, cosine_eff[i], label="Heliostat = {:}".format(i))
            plt.ylabel('Cosine Efficiency [%]')
            plt.xlabel('Hour of Day')
            plt.legend()
            plt.tight_layout()
            plt.show()

            fig = plt.figure()
            for i in plot_eff_heliostats:
                plt.plot(hours, intercept_eff[i], label="Heliostat = {:}".format(i))
            plt.ylabel('Intercept Efficiency [%]')
            plt.xlabel('Hour of Day')
            plt.legend()
            plt.tight_layout()
            plt.show()

            fig = plt.figure()
            for i in plot_eff_heliostats:
                plt.plot(hours, block_eff[i], label="Heliostat = {:}".format(i))
            plt.ylabel('Blocking Efficiency [%]')
            plt.xlabel('Hour of Day')
            plt.legend()
            plt.tight_layout()
            plt.show()

            fig = plt.figure()
            for i in plot_eff_heliostats:
                plt.plot(hours, total_eff[i], label="Heliostat = {:}".format(i))
            plt.ylabel('Total Efficiency [%]')
            plt.xlabel('Hour of Day')
            plt.legend()
            plt.tight_layout()
            plt.show()


    # Canting tests
    if run_canting_tests:
        ###############################
        # No canting 
        ###############################
        PT.stages.clear()
        PT.add_stage()
        # Set-up the heliostat
        position = Point(10.0, -10.0, 5.0)  # arbitrary 
        aperture_size = (12., 8.)
        number_panels = (5, 5)
        gaps = (0.1, 0.1)
        focal_length = 0.0  # zero focus
        aim_point = position + Point(0.0, -1.0, 0.0)  # Straight south

        helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point)
        helio.create_geometry(PT, mirror_optic)
        helio.update_geometry(PT, 180.0, 0.0)   # Straight south

        # Create target - additional stage
        stage = PT.add_stage()
        target = stage.add_element()
        target.position = position + Point(0.0, -15.0, 0.0) 
        target.aim = target.position + Point(0.0, 1.0, 0.0)
        target.zrot = 0.0
        target.aperture_rectangle(aperture_size[0] + 5.0, aperture_size[1] + 5.0)
        target.surface_flat()
        target.optic = target_optic
        target.enabled = True
        helio.set_receiver(target)

        # PT.write_soltrace_input_file('no_canting_test.stinput')

        PT.num_ray_hits = 1.e6
        PT.max_rays_traced = 1.e6*200
        # Turn of errors
        PT.is_sunshape = False
        PT.is_surface_errors = False
        PT.dni = 1000.0
        res = PT.run(123, nthread=12, no_callback=True)

        nx = 250
        ny = 250
        flux_profile = PT.bin_rays(target, nx, ny)
        # FIXME: flux_profile is rotated by 90.0
        # plt.imshow(flux_profile, aspect='equal', cmap='viridis')
        # plt.show()

        target_width, target_height = target.aperture_params[0:2]
        dx = target_width / nx
        dy = target_height / ny

        max_height_idx = -1
        min_height_idx = ny
        max_width_idx = -1
        min_width_idx = nx

        for i, row in enumerate(flux_profile):
            for j, element in enumerate(row):
                if element > 0.0:
                    max_height_idx = max(max_height_idx, j)
                    min_height_idx = min(min_height_idx, j)
                    max_width_idx = max(max_width_idx, i)
                    min_width_idx = min(min_width_idx, i)

        image_height = (max_height_idx - min_height_idx + 1) * dy
        image_width = (max_width_idx - min_width_idx + 1) * dx
        
        assert abs(aperture_size[1] - image_height) < 2*dy
        assert abs(aperture_size[0] - image_width) < 2*dx

        ###############################
        # On-axis canting
        ###############################
        PT.stages.clear()
        PT.add_stage()
        # Set-up the heliostat
        onaxis_canting_distance = 100.0
        helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, onaxis_canting_distance=onaxis_canting_distance)
        helio.create_geometry(PT, mirror_optic)
        helio.update_geometry(PT, 180.0, 0.0)   # Straight south

        # Create target - additional stage
        stage = PT.add_stage()
        target = stage.add_element()
        target.position = position + Point(0.0, -onaxis_canting_distance, 0.0) 
        target.aim = target.position + Point(0.0, 1.0, 0.0)
        target.zrot = 0.0
        target.aperture_rectangle(aperture_size[0] + 5.0, aperture_size[1] + 5.0)
        target.surface_flat()
        target.optic = target_optic
        target.enabled = True
        helio.set_receiver(target)

        # PT.write_soltrace_input_file('onaxis_canting_test.stinput')

        PT.num_ray_hits = 1.e6
        PT.max_rays_traced = 1.e6*200
        # Turn of errors
        PT.is_sunshape = False
        PT.is_surface_errors = False
        PT.dni = 1000.0
        res = PT.run(123, nthread=12, no_callback=True)
        # PT.plot_trace()

        nx = 250
        ny = 250
        flux_profile = PT.bin_rays(target, nx, ny)
        # FIXME: flux_profile is rotated by 90.0
        # plt.imshow(flux_profile, aspect='equal', cmap='viridis')
        # plt.show()

        target_width, target_height = target.aperture_params[0:2]
        dx = target_width / nx
        dy = target_height / ny

        max_height_idx = -1
        min_height_idx = ny
        max_width_idx = -1
        min_width_idx = nx

        for i, row in enumerate(flux_profile):
            for j, element in enumerate(row):
                if element > 0.0:
                    max_height_idx = max(max_height_idx, j)
                    min_height_idx = min(min_height_idx, j)
                    max_width_idx = max(max_width_idx, i)
                    min_width_idx = min(min_width_idx, i)

        image_height = (max_height_idx - min_height_idx + 1) * dy
        image_width = (max_width_idx - min_width_idx + 1) * dx

        # Panel height and width
        expected_height = (aperture_size[1] - gaps[1] * (number_panels[1] - 1)) / number_panels[1]
        expected_width = (aperture_size[0] - gaps[0] * (number_panels[0] - 1)) / number_panels[0]
        
        assert abs(expected_height - image_height) < 2*dy
        assert abs(expected_width - image_width) < 2*dx

        ###############################
        # Off-axis canting - Straight on sun position - should be same as on-axis canting
        ###############################
        PT.stages.clear()
        PT.add_stage()
        # Set-up the heliostat
        offaxis_canting_sun_position = (180.0, 90.0)
        aim_point = position + Point(0.0, -100.0, 0.0)  # Straight south
        helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, offaxis_canting_sun_position=offaxis_canting_sun_position)
        helio.create_geometry(PT, mirror_optic)
        helio.update_geometry(PT, 180.0, 0.0)   # Straight south

        # Create target - additional stage
        stage = PT.add_stage()
        target = stage.add_element()
        target.position = aim_point
        target.aim = target.position + Point(0.0, 1.0, 0.0)
        target.zrot = 0.0
        target.aperture_rectangle(aperture_size[0] + 5.0, aperture_size[1] + 5.0)
        target.surface_flat()
        target.optic = target_optic
        target.enabled = True
        helio.set_receiver(target)

        PT.write_soltrace_input_file('offaxis_canting_test.stinput')
 
        PT.num_ray_hits = 1.e6
        PT.max_rays_traced = 1.e6*200
        # Turn of errors
        PT.is_sunshape = False
        PT.is_surface_errors = False

        PT.dni = 1000.0
        # TODO: This corrects the center panel rotation issue
        # Without rounding, center panel is rotated by 90 degrees in API
        # GUI from write_soltrace_input_file() does not have this issue because the stinput files are written to 6 decimal places
        # TODO: Create a unit test that reproduces this issue
        for mir in helio.get_elements():
            digits = 13
            position = mir.position
            if mir.id == 12:
                print("Element Id: {:}".format(mir.id))
                print("Position Error: {}, {}, {}".format(position.x - round(position.x, digits), position.y - round(position.y, digits), position.z - round(position.z, digits)))
            mir.position = Point(round(position.x, digits), round(position.y, digits), round(position.z, digits))

            aim = mir.aim
            if mir.id == 12:
                print("Aim Error: {}, {}, {}".format(aim.x - round(aim.x, digits), aim.y - round(aim.y, digits), aim.z - round(aim.z, digits)))
            mir.aim = Point(round(aim.x, digits), round(aim.y, digits), round(aim.z, digits))

            zrot = mir.zrot
            mir.zrot = round(zrot, digits)

        res = PT.run(123, nthread=12, no_callback=True)
        # PT.plot_trace()

        nx = 250
        ny = 250
        flux_profile = PT.bin_rays(target, nx, ny)
        # FIXME: flux_profile is rotated by 90.0
        # plt.imshow(flux_profile, aspect='equal', cmap='viridis')
        # plt.show()

        target_width, target_height = target.aperture_params[0:2]
        dx = target_width / nx
        dy = target_height / ny

        max_height_idx = -1
        min_height_idx = ny
        max_width_idx = -1
        min_width_idx = nx

        for i, row in enumerate(flux_profile):
            for j, element in enumerate(row):
                if element > 0.0:
                    max_height_idx = max(max_height_idx, j)
                    min_height_idx = min(min_height_idx, j)
                    max_width_idx = max(max_width_idx, i)
                    min_width_idx = min(min_width_idx, i)

        image_height = (max_height_idx - min_height_idx + 1) * dy
        image_width = (max_width_idx - min_width_idx + 1) * dx

        # Panel height and width
        expected_height = (aperture_size[1] - gaps[1] * (number_panels[1] - 1)) / number_panels[1]
        expected_width = (aperture_size[0] - gaps[0] * (number_panels[0] - 1)) / number_panels[0]
        
        assert abs(expected_height - image_height) < 2*dy
        assert abs(expected_width - image_width) < 2*dx

        ###############################
        # Off-axis canting - 45 degrees reflection angle
        ###############################
        PT.stages.clear()
        PT.add_stage()
        # Set-up the heliostat
        position = Point(50.0, 50.0, 5.0)  # at a 45 degree from target
        offaxis_canting_sun_position = (135., 90.0)
        aim_point = Point(0.0, 0.0, 5.0)
        helio = heliostat(position, aperture_size, number_panels, gaps, focal_length, aim_point, offaxis_canting_sun_position=offaxis_canting_sun_position)
        helio.create_geometry(PT, mirror_optic)
        helio.update_geometry(PT, offaxis_canting_sun_position[0], 90.0 - offaxis_canting_sun_position[1])

        # Create target - additional stage
        stage = PT.add_stage()
        target = stage.add_element()
        target.position = aim_point
        target.aim = target.position + Point(0.0, 1.0, 0.0)
        target.zrot = 0.0
        target.aperture_rectangle(aperture_size[0] + 5.0, aperture_size[1] + 5.0)
        target.surface_flat()
        target.optic = target_optic
        target.enabled = True
        helio.set_receiver(target)

        # PT.write_soltrace_input_file('offaxis_canting_45deg_test.stinput')

        PT.num_ray_hits = 1.e6
        PT.max_rays_traced = 1.e6*200
        # Turn of errors
        PT.is_sunshape = False
        PT.is_surface_errors = False
        PT.dni = 1000.0
        # TODO: Fix this in Pysoltrace
        for mir in helio.get_elements():
            digits = 13
            position = mir.position
            mir.position = Point(round(position.x, digits), round(position.y, digits), round(position.z, digits))
            aim = mir.aim
            mir.aim = Point(round(aim.x, digits), round(aim.y, digits), round(aim.z, digits))
            zrot = mir.zrot
            mir.zrot = round(zrot, digits)

        res = PT.run(123, nthread=12, no_callback=True)
        # PT.plot_trace()

        nx = 250
        ny = 250
        flux_profile = PT.bin_rays(target, nx, ny)
        # FIXME: flux_profile is rotated by 90.0
        # plt.imshow(flux_profile, aspect='equal', cmap='viridis')
        # plt.show()

        target_width, target_height = target.aperture_params[0:2]
        dx = target_width / nx
        dy = target_height / ny

        max_height_idx = -1
        min_height_idx = ny
        max_width_idx = -1
        min_width_idx = nx

        for i, row in enumerate(flux_profile):
            for j, element in enumerate(row):
                if element > 0.0:
                    max_height_idx = max(max_height_idx, j)
                    min_height_idx = min(min_height_idx, j)
                    max_width_idx = max(max_width_idx, i)
                    min_width_idx = min(min_width_idx, i)

        image_height = (max_height_idx - min_height_idx + 1) * dy
        image_width = (max_width_idx - min_width_idx + 1) * dx

        # Panel height and width
        expected_height = (aperture_size[1] - gaps[1] * (number_panels[1] - 1)) / number_panels[1]
        expected_width = (aperture_size[0] - gaps[0] * (number_panels[0] - 1)) / number_panels[0]
        
        assert abs(expected_height - image_height) < 2*dy
        assert abs(expected_width - image_width) < 2*dx



# TODO:
#   - Add tracking angle limits?
#   - Add canting and tracking error (elevation and azimuth)
