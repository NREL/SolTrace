"""
@author: whamilton

Creates a linear Fresnel template made of multiple mirror panels.
Linear Fresnel template updates positioning based on sun position with limits on tracking angle.
"""

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..', '..'))

from pysoltrace import PySolTrace, Point
from parabolic_trough import compute_euler_angles, sun_vector, trace, sun_position, meinel_clearsky

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 15})


class linear_fresnel:
    
    def __init__(self, 
                 position: Point, 
                 aperture_size: tuple, 
                 number_panels: tuple, 
                 gaps: tuple,
                 azimuth: float, 
                 tilt: float,
                 receiver_height: float, 
                 receiver_dimensions: tuple,
                 focused_panels: bool = True,
                 tracking_limits: tuple = (10.0, 170.0)):
        """
        :param position: [m] Position of the linear fresnel
        :param aperture_size: [m] Size of the aperture (width, length)
        :param number_panels: [-] Number of panels (width, length)
        :param gaps: [m] Gaps between the panels (width, length, center*) *optional
        :param azimuth: [deg] Azimuth angle of the linear fresnel (0 = North, 90 = East)
        :param tilt: [deg] Tilt angle of the linear fresnel (0 = Horizontal, 90 = Vertical)
        :param receiver_height: [m] Focal length of the linear fresnel
        :param receiver_dimenisons: [m] Absorber diameter, envelope diameter, envelope thickness, length
        :param focused_panels: [-] True if the panels are focused, False if they are flat
        :param tracking_limits (optional): [deg] minimum and maximum tracking angles of the linear fresnel panels
        TODO: simplified receiver (flat plate)
        """
        self.position = position

        if len(aperture_size) != 2:
            raise ValueError("Aperture size should be a tuple with 2 elements")
        if any(ap_dim < 0.0 for ap_dim in aperture_size):
            raise ValueError("Aperture dimensions must be positive values")
        self.aperture_size = aperture_size

        if len(number_panels) != 2:
            raise ValueError("Number of panels should be a tuple with 2 elements")
        if number_panels[0] % 2 != 0 and len(gaps) == 3:
            raise ValueError("Number of panels in the width should be an even number if center gap is provided")
        if any(panels < 0.0 for panels in number_panels):
            raise ValueError("Number of panels must be positive values")
        self.number_panels = number_panels

        if len(gaps) != 2 and len(gaps) != 3:
            raise ValueError("Gaps should be a tuple with 2 or 3 elements")
        if any(gap < 0.0 for gap in gaps):
            raise ValueError("Gaps must be positive values")
        self.gaps = gaps

        if receiver_height < 0.0:
            raise ValueError("Receiver height must be a positive value")
        self.receiver_height = receiver_height

        if azimuth > 180.0 or azimuth < -180.0:
            raise ValueError("Azimuth angle should be between -180 and 180 degrees")
        self.azimuth = azimuth

        if tilt > 90.0 or tilt < 0.0:
            raise ValueError("Tilt angle should be between 0 and 90 degrees")
        self.tilt = tilt

        if len(receiver_dimensions) != 3:
            raise ValueError("Receiver dimensions parameter should be a tuple with 3 elements")
        self.absorber_diameter = receiver_dimensions[0]
        self.envelope_diameter = receiver_dimensions[1]
        self.envelope_thickness = receiver_dimensions[2]
        if self.envelope_diameter + 2 * self.envelope_thickness < self.absorber_diameter:
            raise ValueError("Receiver envelope must be larger than absorber diameter")
        
        self.focused_panels = focused_panels

        if len(tracking_limits) != 2:
            raise ValueError("Tracking limits parameter should be a tuple with 2 elements")
        if tracking_limits[1] < tracking_limits[0]:
            raise ValueError("Maximum tracking limit angle must be greater than minimum tracking limit angle")
        self.tracking_limits = tracking_limits
        
        # Calculated values
        self.mirrors = list()
        self.absorbers = list()
        self.envelope = list()

        self.tracking_angles = list()
        self.aperture_normal = list()   # List of panel normals (per number of panel in width direction)
        self.panel2receiver = list()    # List of panel to receiver unit vectors (per number of panel in width direction) - calculate once during create
        self.track_rotation_axis = Point()      # Constant for all mirrors

    def get_elements(self):
        return self.mirrors + self.absorbers + self.envelope

    def create_geometry(self, PT:PySolTrace, optics: list):
        """
        Creates the elements of the linear Fresnel in the PySolTrace object. 
        NOTE: This assumes the stage has an origin of (0.0, 0.0, 0.0), Aimpoint of (0.0, 0.0, 1.0), and a z-rotation of 0.0 degrees
        
        :param PT: PySolTrace object
        :param optics: list of surface optics (mirror, absorber, outer envelope, inner envelope)
        """
        # Determine panel width and length
        panel_width = (self.aperture_size[0] - self.gaps[0] * (self.number_panels[0] - 1))
        if len(self.gaps) == 3:
            panel_width += self.gaps[0] - self.gaps[2] # add back a gap and subtract the center gap
        panel_width /= self.number_panels[0]
        panel_length = (self.aperture_size[1] - self.gaps[1] * (self.number_panels[1] - 1)) / self.number_panels[1]

        # Build linear Fresnel about the origin then translate to position
        stage = PT.stages[0]        # NOTE: assuming everything is in the first stage (removing stages)
        # Create mirror panels
        self.mirrors.clear()
        self.tracking_angles.clear()
        self.aperture_normal.clear()
        self.panel2receiver.clear()
        panel_x = - self.aperture_size[0] / 2 + panel_width / 2
        for i in range(self.number_panels[0]):          # width
            panel_y = - self.aperture_size[1] / 2 + panel_length / 2
            self.tracking_angles.append(90.0)           # Pointing straight up
            self.aperture_normal.append(Point(0.0, 0.0, 1.0))   # Pointing straight up
            self.panel2receiver.append(Point(-panel_x, 0.0, self.receiver_height).unitize())
            for j in range(self.number_panels[1]):      # length - mirrors that rotate the same are together
                element = stage.add_element()
                element.position = Point(panel_x, panel_y, 0.0)
                element.aim = Point(panel_x, panel_y, 1.0)          # point straight up
                element.zrot = 0.0
                element.aperture_rectangle(panel_width, panel_length)
                if self.focused_panels:
                    panel_focal_length = np.sqrt(panel_x**2 + self.receiver_height**2)
                    element.surface_parabolic(panel_focal_length, 0.0)
                else:
                    element.surface_flat()
                element.optic = optics[0]
                element.enabled = True

                self.mirrors.append(element)

                panel_y += panel_length + self.gaps[1]

            panel_x += panel_width 
            if self.number_panels[0] % 2 == 0 and j == self.number_panels[0] // 2 - 1:
                panel_x += self.gaps[2]
            else:
                panel_x += self.gaps[0]

        # Adding receiver
        # TODO: Assuming a single receiver per length of aperture - this should be broken up into multiple receiver tubes.
        # Absorber tube
        self.absorbers.clear()
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.receiver_height - self.absorber_diameter/2)
        element.aim = element.position + Point(0.0, 0.0, 1.0)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical(self.absorber_diameter/2)
        element.optic = optics[1]
        element.enabled = True
        self.absorbers.append(element)

        # Outer Envelope
        self.envelope.clear()
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.receiver_height - self.envelope_diameter/2)
        element.aim = element.position + Point(0.0, 0.0, 1.0)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical(self.envelope_diameter/2)
        element.optic = optics[2]
        element.interaction = 1
        element.enabled = True
        self.envelope.append(element)

        # Inner Envelope
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.receiver_height - ((self.envelope_diameter / 2) - self.envelope_thickness))
        element.aim = element.position + Point(0.0, 0.0, 1.0)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical((self.envelope_diameter / 2) - self.envelope_thickness)
        element.optic = optics[3]
        element.interaction = 1
        element.enabled = True
        self.envelope.append(element)

        # Azimuth and tilt convert to radians
        tilt_rads = self.tilt * (np.pi / 180.0)
        azimuth_rads = self.azimuth * (np.pi / 180.0)

        # Calculating aperture normal
        norm = Point(0.0, 0.0, 1.0)
        norm = PT.util_rotation_arbitrary( tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), norm)           # Rotation about the x-axis
        norm = PT.util_rotation_arbitrary( -azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), norm)       # Rotation about the z-axis
        for i,_ in enumerate(self.aperture_normal):
            self.aperture_normal[i] = norm

        # Updating panel to receiver vector
        for i, p2r in enumerate(self.panel2receiver):
            p2r = PT.util_rotation_arbitrary( tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), p2r)           # Rotation about the x-axis
            p2r = PT.util_rotation_arbitrary( -azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), p2r)       # Rotation about the z-axis
            self.panel2receiver[i] = p2r

        # Calculating the rotation axis - all panels have the same tracking axis
        track_axis = Point(0.0, 1.0, 0.0)
        track_axis = PT.util_rotation_arbitrary( tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), track_axis)           # Rotation about the x-axis
        track_axis = PT.util_rotation_arbitrary( -azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), track_axis)       # Rotation about the z-axis
        self.track_rotation_axis = track_axis
        
        for norm in self.aperture_normal:
            assert abs(np.dot(norm.as_list(), self.track_rotation_axis.as_list())) <= 1.e-6    # This should always be true        

        x_axis = np.cross(self.track_rotation_axis.as_list(), self.aperture_normal[0].as_list())
        x_axis = Point(x_axis[0], x_axis[1], x_axis[2]).unitize()

        # Rotation matrix from stage to local coordinates
        rotation_matrix = np.array([[x_axis.x, self.track_rotation_axis.x, self.aperture_normal[0].x], 
                                    [x_axis.y, self.track_rotation_axis.y, self.aperture_normal[0].y],
                                    [x_axis.z, self.track_rotation_axis.z, self.aperture_normal[0].z]])
        
        assert abs(np.linalg.det(rotation_matrix) - 1) <= 1.e-8
        euler_angles = compute_euler_angles(np.linalg.matrix_transpose(rotation_matrix))   # NOTE: Stage coordinates are assumed to be Identity matrix
        # NOTE: Transpose of a rotation matrix is equal to its inversve.
        #print("Euler angles (initial) {}".format(euler_angles))

        # Updating element positions, aimpoints, and z-rotations based on tilt and azimuth
        for element in self.get_elements(): 
            # Tilt first then azimuth, negative azimuth due to right-hand rule
            element.position = PT.util_rotation_arbitrary(tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), element.position)
            element.position = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), element.position)

            element.aim = PT.util_rotation_arbitrary(tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), element.aim)
            element.aim = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), element.aim)

            element.zrot = euler_angles[2]

            # Translate to specific location
            element.position += self.position
            element.aim += self.position


    def update_geometry(self, PT:PySolTrace, azimuth: float, elevation: float):
        """
        Update linear Fresnel geometry based on solar position (azimuth and elevation)

        :param PT: PySolTrace instance
        :param azimuth: [deg] Solar azimuth angle
        :param elevation: [deg] Solar elevation angle
        """
        # Calculate and update solar position
        sun_vec = sun_vector(azimuth, elevation)
        PT.sun.position = sun_vec * 1000.0

        # Determine change in tracking angle
        # Project the sun location onto the plane with a norm of the track rotation vector
        distance = np.dot((sun_vec * -1).as_list(), self.track_rotation_axis.as_list())
        sun_proj = sun_vec + (self.track_rotation_axis * distance)
        sun_proj = sun_proj.unitize()

        mirror_idx = 0
        for i in range(self.number_panels[0]):          # width
            # Calculate new normal
            prev_norm = self.aperture_normal[i]
            self.aperture_normal[i] = (self.panel2receiver[i] + sun_proj).unitize()
            delta_tracking_angle = np.arccos(np.dot(prev_norm.as_list(), self.aperture_normal[i].as_list()) 
                                             / (np.linalg.norm(prev_norm.as_list()) * np.linalg.norm(self.aperture_normal[i].as_list())))
            
            # Determine the direction of rotation
            direction = np.cross(prev_norm.as_list(), self.aperture_normal[i].as_list())
            direction = np.dot(direction, self.track_rotation_axis.as_list())
            if direction < 0:
                delta_tracking_angle *= -1.0

            # Update tracking angle and check limits, adjust if exceeds limits
            self.tracking_angles[i] += delta_tracking_angle * 180.0 / np.pi
            if self.tracking_angles[i] < self.tracking_limits[0]:
                delta_tracking_angle += (self.tracking_limits[0] - self.tracking_angles[i]) * np.pi / 180.0
                self.tracking_angles[i] = self.tracking_limits[0]
            elif self.tracking_angles[i] > self.tracking_limits[1]:
                delta_tracking_angle -= (self.tracking_angles[i] - self.tracking_limits[1]) * np.pi / 180.0
                self.tracking_angles[i] = self.tracking_limits[1]

            assert abs(np.dot(self.aperture_normal[i].as_list(), self.track_rotation_axis.as_list())) <= 1.e-5    # This should always be true

            # calculate x-axis by taking the cross product
            x_axis = np.cross(self.track_rotation_axis.as_list(), self.aperture_normal[i].as_list())
            x_axis = Point(x_axis[0], x_axis[1], x_axis[2]).unitize()

            x_loc = np.array([[x_axis.x, self.track_rotation_axis.x, self.aperture_normal[i].x], 
                              [x_axis.y, self.track_rotation_axis.y, self.aperture_normal[i].y],
                              [x_axis.z, self.track_rotation_axis.z, self.aperture_normal[i].z]])
            
            assert abs(np.linalg.det(x_loc) - 1) <= 1.e-8
            euler_angles = compute_euler_angles(np.linalg.matrix_transpose(x_loc))   # Stage coordinates is Indentity
            
            # Update aimpoint and z-rotation of mirrors
            for j in range(self.number_panels[1]):      # length - mirrors that rotate the same are together
                element = self.mirrors[mirror_idx]
                element.aim = PT.util_rotation_arbitrary(delta_tracking_angle, self.track_rotation_axis, element.position, element.aim)
                element.zrot = euler_angles[2]
                mirror_idx += 1
            

    def calculate_receiver_power(self, PT: PySolTrace):
        """
        Calculates the power absorbed by the receiver tubes

        :param PT: PySolTrace instance

        :return total absorbed power: [kWt]
        """
        # TODO: We could expand post processing to include receiver losses, optical efficiency
        # Shading and Cosine, Reflection, Intercept, Absorption, Overall Optical Efficiency
        # These calculations would require the full field of as a "linear Fresnel System" information.
        rayData = PT.raydata

        absorbed_power = 0.0
        for rec in self.absorbers:
            absorbed_power += len(rayData[rayData['element'] == -(rec.id+1)]) * PT.powerperray / 1.e3

        return absorbed_power
    

if __name__ == "__main__":

    solar_tracking_test = True

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

    absorber_optic = PT.add_optic('absorber')
    absorber_optic.front.reflectivity = 0.04
    absorber_optic.front.slope_error = 1e-5
    absorber_optic.front.spec_error = 1e-5
    absorber_optic.back.reflectivity = 0.04
    absorber_optic.back.slope_error = 1e-5
    absorber_optic.back.spec_error = 1e-5

    outer_env_optic = PT.add_optic('outer_env')
    outer_env_optic.front.reflectivity = 0.0
    outer_env_optic.front.transmissivity = 0.98
    outer_env_optic.front.refraction_real = 1.46
    outer_env_optic.front.slope_error = 1e-4
    outer_env_optic.front.spec_error = 1e-4
    outer_env_optic.back.reflectivity = 0.0
    outer_env_optic.back.transmissivity = 0.98
    outer_env_optic.back.refraction_real = 1.0
    outer_env_optic.back.slope_error = 1e-4
    outer_env_optic.back.spec_error = 1e-4

    inner_env_optic = PT.add_optic('inner_env')
    inner_env_optic.front.reflectivity = 0.0
    inner_env_optic.front.transmissivity = 0.965
    inner_env_optic.front.refraction_real = 1.0
    inner_env_optic.front.slope_error = 1e-4
    inner_env_optic.front.spec_error = 1e-4
    inner_env_optic.back.reflectivity = 0.0
    inner_env_optic.back.transmissivity = 0.965
    inner_env_optic.back.refraction_real = 1.46
    inner_env_optic.back.slope_error = 1e-4
    inner_env_optic.back.spec_error = 1e-4

    optics = [mirror_optic, absorber_optic, outer_env_optic, inner_env_optic]

    # Set-up the parabolic trough
    position = Point(10.0, -10.0, 0.0)
    aperture_size = (6.0, 12.0)
    number_panels = (10, 4)
    gaps = (0.15, 0.02, 0.15)
    receiver_height = 2.0
    azimuth = 25.0 #45.0
    tilt = 15.0#15.0
    receiver_dimensions = (0.07, 0.115, 0.003)
    # TODO: figure out handling of multiple receivers

    linearFresnel = linear_fresnel(position, aperture_size, number_panels, gaps, azimuth, tilt, receiver_height, receiver_dimensions)
    linearFresnel.create_geometry(PT, optics)
    # PT.write_soltrace_input_file('linear_fresnel.stinput')

    linearFresnel.update_geometry(PT, 135.0, 60.)
    # PT.write_soltrace_input_file('linear_fresnel_update.stinput')
    trace(PT, nthreads=12)
    rec_power = linearFresnel.calculate_receiver_power(PT)
    print("Power on the receiver: {:.2f} [kWt]".format(rec_power))

    linearFresnel.update_geometry(PT, 180.0, 90.0)
    # PT.write_soltrace_input_file('linear_fresnel_overhead.stinput')
    trace(PT, nthreads=12)
    rec_power = linearFresnel.calculate_receiver_power(PT)
    print("Power on the receiver: {:.2f} [kWt]".format(rec_power))


    # Sun tracking - with multiple Linear Fresnel
    if solar_tracking_test:
        latitude = 34.8639          # Daggett, CA
        altitude = 610.0 / 1000.0   # [km]
        day_of_year = 80

        # Linear Fresnel set-up
        positions = [Point(10.0, 0.0, 0.0), Point(0.0, 0.0, 0.0), Point(-10.0, 0.0, 0.0)]    # north-south
        aperture_size = (6.0, 12.0)
        number_panels = (10, 1)
        gaps = (0.15, 0.02, 0.15)
        receiver_height = 2.0
        azimuth = 0.0
        tilt = 0.0
        receiver_dimensions = (0.07, 0.115, 0.003)

        PT.stages.clear()
        PT.add_stage()
        linearFresnels = list()
        for position in positions:
            lf = linear_fresnel(position, aperture_size, number_panels, gaps, azimuth, tilt, receiver_height, receiver_dimensions)
            lf.create_geometry(PT, optics)
            linearFresnels.append(lf)

        # Time-simulation
        rec_power = list()
        dni = list()

        # hours = [10 + h*0.2 for h in range(11)]
        # hours = [5 + h*0.5 for h in range(29)]
        hours = [5 + h for h in range(15)]
        for h in hours:
            print("Simulating hour {:}...".format(h))
            azimuth, elevation = sun_position(latitude, day_of_year, h)    
            zenith = (90.0 - elevation) * np.pi / 180.0
            dni_hour = meinel_clearsky(day_of_year, zenith, altitude)   # Clear sky DNI values
            dni.append(dni_hour)

            for lf in linearFresnels:
                lf.update_geometry(PT, azimuth, elevation)

            if h in [6, 8, 10, 12, 17]:
                PT.write_soltrace_input_file('linear_fresnel_{}_hour.stinput'.format(h))

            if dni_hour > 0.0:
                trace(PT, dni=dni_hour, nthreads=12)
                # Sum power on receivers
                power = 0.0
                for lf in linearFresnels:
                    power += lf.calculate_receiver_power(PT)
                rec_power.append(power)
            else:
                rec_power.append(0.0)

        # Plot power and dni over time
        fig = plt.figure()
        plt.plot(hours, rec_power, 'k')
        plt.ylabel('Receiver Power [kWt]')
        plt.xlabel('Hour of Day')

        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(hours, dni, 'r')
        ax2.set_ylabel(r"Direct Normal Irradiance [W/m$^2$]", color = 'r')
        ax2.tick_params('y', colors = 'r')
        plt.tight_layout()
        plt.show()
