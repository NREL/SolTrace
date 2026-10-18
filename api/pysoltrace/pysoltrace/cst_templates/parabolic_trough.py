"""
@author: whamilton

Creates a parabolic trough template made of multiple mirror panels.
Trough template updates positioning based on sun position with limits on tracking angle.
"""

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..', '..'))

from pysoltrace import PySolTrace, Point

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 15})

def meinel_clearsky(day: int, zenith: float, altitude: float):
    """
    The Meinel model calculates solar intensity as a function of extraterrestrial radiation, site altitude, and solar zenith angle.
    NOTE: Taken from SolarPILOT Ambient.cpp

    :param day: [-] day of the year
    :param zenith: [radians] solar zenith angle
    :param altitude: [km] altitude of location
    
    :return dni [W/m^2]
    """
    czen = np.cos(zenith)
    s0 = 1.353*(1.+.0335*np.cos(2.*np.pi*(day+10.)/365.))  #[kW/m^2]      # NOTE: This is a different equation than the one in SolarPILOT's help documentation 
    if czen <= 0.0:
        dni = 0.0
    else:
        dni = (1.-.14*altitude)*np.exp(-.357/pow(czen,.678))+.14*altitude
    return s0 * dni * 1000.0 # [W/m^2]

def sun_position(latitude: float, day: int, hour: float)->tuple:
    """
    Computes the sun vector xyz given arguments (Copied from st_sun_position from stapi.cpp)
    TODO: This is just a place holder function. We should probably use NREL's Solar Position Algorithm (SPA)

	:param latitude: [deg] latitude 
	:param day: [-] day of the year 
	:param hour: [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

    :return (azimuth, elevation): [degrees]
    """
    D2R = np.pi / 180
    R2D = 180 / np.pi
    declination = R2D * np.asin(0.39795 * np.cos(0.98563 * D2R * (day - 173)))
    hourAngle = 15 * (hour - 12)
    elevation = R2D * np.asin(np.sin(declination * D2R) * np.sin(latitude * D2R) + np.cos(declination * D2R) * np.cos(hourAngle * D2R) * np.cos(latitude * D2R))
    azimuth = R2D * np.acos((np.sin(D2R * declination) * np.cos(D2R * latitude) - np.cos(D2R * declination) * np.sin(D2R * latitude) * np.cos(D2R * hourAngle)) / np.cos(D2R * elevation) + 0.0000000001)
    if (np.sin(hourAngle * D2R) > 0.0):
        azimuth = 360 - azimuth
    return (azimuth, elevation)

def sun_vector(azimuth: float, elevation: float)->Point:
    """
    Computes the sun vector xyz given arguments (Copied from st_sun_position from stapi.cpp)

    Assumes xyz coordinate system:
		x: +east
		y: +north
		z: +zenith

	:param azimuth: [deg] solar azimuth 
	:param elevation: [deg] solar elevation 

    :return solar position (x, y, z):
    """
    D2R = np.pi / 180
    x = np.sin(azimuth * D2R) * np.cos(elevation * D2R)
    y = np.cos(azimuth * D2R) * np.cos(elevation * D2R)
    z = np.sin(elevation * D2R)
    return Point(x, y, z)

def sun_vector_from_latitude(latitude: float, day: int, hour: float)->Point:
    """
    Computes the sun vector xyz given arguments (Copied from st_sun_position from stapi.cpp)

    Assumes xyz coordinate system:
		x: +east
		y: +north
		z: +zenith

	:param latitude: [deg] latitude 
	:param day: [-] day of the year 
	:param hour: [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

    :return solar position (x, y, z):
    """
    azimuth, elevation = sun_position(latitude, day, hour)
    return sun_vector(azimuth, elevation)

def compute_euler_angles(rotation: np.array):
    """
    Computes euler angles based on rotation matrix. This is consistant with the assumed rotation matrax within SolTrace.
    Source: https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf

    :param rotation: [3x3] matrix

    :return euler angles: [deg] array [beta, alpha, gamma]
    """
    if rotation[2,1] != 1.0 and rotation[2,1] != -1.0:
        beta = np.asin(rotation[2,1])
        alpha = np.atan2(rotation[2,0] / np.cos(beta), rotation[2,2] / np.cos(beta))
        gamma = - np.atan2(rotation[0,1] / np.cos(beta), rotation[1,1] / np.cos(beta))
    else:
        gamma = 0.0
        if rotation[2,1] == 1.0:
            beta = np.pi / 2
            alpha = np.atan2(rotation[0,2], rotation[1,2])
        else:
            beta = - np.pi / 2
            alpha = np.atan2(-rotation[0,2], -rotation[1,2])

    return np.array([beta, alpha, gamma])*180/np.pi

def simpsons_rule(f, a, b, n):
    if n % 2 == 1:  # Simpson's rule requires even number of intervals
        n += 1
    x = np.linspace(a, b, n + 1)
    y = f(x)
    h = (b - a) / n
    return (h / 3) * (y[0] + 4 * np.sum(y[1:n:2]) + 2 * np.sum(y[2:n-1:2]) + y[n])

class parabolic_trough:
    
    def __init__(self, 
                 position: Point, 
                 aperture_size: tuple, 
                 number_panels: tuple, 
                 gaps: tuple, 
                 focal_length: float, 
                 azimuth: float, 
                 tilt: float,
                 receiver_dimensions: tuple,
                 tracking_limits: tuple = (10.0, 170.0)):
        """
        :param position: [m] Position of the parabolic trough
        :param aperture_size: [m] Size of the aperture (width, length)
        :param number_panels: [-] Number of panels (width, length)
        :param gaps: [m] Gaps between the panels (width, length, center*) *optional
        :param focal_length: [m] Focal length of the parabolic trough
        :param azimuth: [deg] Azimuth angle of the parabolic trough (0 = North, 90 = East)
        :param tilt: [deg] Tilt angle of the parabolic trough (0 = Horizontal, 90 = Vertical)
        :param receiver_dimenisons: [m] Absorber diameter, envelope diameter, envelope thickness, length
        :param tracking_limits (optional): [deg] minimum and maximum tracking angles of the trough
        """
        self.position = position

        if len(aperture_size) != 2:
            raise ValueError("Aperture size should be a tuple with 2 elements")
        if any(ap_dim < 0.0 for ap_dim in aperture_size):
            raise ValueError("Aperture dimensions must be positive values")
        self.aperture_size = aperture_size

        if len(number_panels) != 2:
            raise ValueError("Number of panels should be a tuple with 2 elements")
        if number_panels[0] % 2 != 0 and number_panels[0] != 1:
            raise ValueError("Number of panels in the width should be an even number if greater than 1")
        if any(panels < 0.0 for panels in number_panels):
            raise ValueError("Number of panels must be positive values")
        self.number_panels = number_panels

        if len(gaps) != 2 and len(gaps) != 3:
            raise ValueError("Gaps should be a tuple with 2 or 3 elements")
        if any(gap < 0.0 for gap in gaps):
            raise ValueError("Gaps must be positive values")
        self.gaps = gaps

        if focal_length < 0.0:
            raise ValueError("Focal length must be a positive value")
        self.focal_length = focal_length

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
        
        if len(tracking_limits) != 2:
            raise ValueError("Tracking limits parameter should be a tuple with 2 elements")
        if tracking_limits[1] < tracking_limits[0]:
            raise ValueError("Maximum tracking limit angle must be greater than minimum tracking limit angle")
        self.tracking_limits = tracking_limits
        
        # Calculated values
        self.flip_element = dict()

        self.mirrors = list()
        self.absorbers = list()
        self.envelope = list()

        self.tracking_angle = 90.0      # pointing straight up
        self.x_axis = Point()
        self.track_rotation_axis = Point()
        self.aperture_normal = Point()

    def get_elements(self):
        return self.mirrors + self.absorbers + self.envelope

    def parabolic_arc_length_equation(self, x: float):
        """
        Equation for integrating the arc length of a 2D parabola.

        :param x: x-coordinate
        :return y: y-coordinate
        """
        c_x = 1 / (2 * self.focal_length)
        return (1 + (c_x * x)**2)**(1/2)

    def determine_end_x_coordinate(self, start:float, distance:float, dx:float = 1e-6) -> float:
        """
        Determines the ending x-coordinate based on the start x-coordinate and the distance along the parabolic arc (using trapezoidal rule).
        NOTE: we could use simpson's rule instead

        :param start: starting x-coordinate
        :param distance: distance along the parabolic arc to travel
        :param dx: step size for numerical integration
        :return: x-coordinate at the end of distance traveled along the parabolic arc
        """
        distance_traveled = 0
        # Starting step
        x = start
        y = self.parabolic_arc_length_equation(x)
        distance_traveled += 0.5 * y * dx
        while True:
            x += dx
            y = self.parabolic_arc_length_equation(x)
            distance_traveled += 0.5 * y * dx
            if distance_traveled < distance:
                distance_traveled += 0.5 * y * dx
            else:
                break

        error = distance_traveled - distance
        if abs(error) > 1e-4:
            raise ValueError("Error in determining the x-coordinate of the parabolic arc length is too high, error = {}".format(error))
        return x

    def create_geometry(self, PT:PySolTrace, optics: list):
        """
        Creates the elements of the parabolic trough in the PySolTrace object. 
        NOTE: This assumes the stage has an origin of (0.0, 0.0, 0.0), Aimpoint of (0.0, 0.0, 1.0), and a z-rotation of 0.0 degrees
        
        :param PT: PySolTrace object
        :param optics: list of surface optics (mirror, absorber, outer envelope, inner envelope)
        """
        # Determine the aperture bounds of each panel
        panel_aperture_bounds = []
        half_width = self.aperture_size[0] / 2
        if self.number_panels[0] == 1:
            panel_aperture_bounds = [(-half_width, half_width)]
        elif self.number_panels[0] % 2 == 0:    # Even number of panels
            # Analytical solution of arc length of parabola
            # Problem: Integral of sqrt(1 + f'(x)^2) dx from -half_width to half_width
            # where f(x) = 1/2 * c_x * x^2
            c_x = 1 / (2 * self.focal_length)
            arc_length = half_width * ((half_width * c_x)**2 + 1)**(1/2) + np.arcsinh(half_width * c_x) / c_x
            # Confirm result with numerical integration
            check_x = self.determine_end_x_coordinate(-half_width, arc_length)
            assert abs(check_x - half_width) < 1e-4, "Error in determining the x-coordinate of the parabolic arc length, x = {}".format(check_x)

            panel_arc_length = arc_length - self.gaps[0] * (self.number_panels[0] - 1)
            if len(self.gaps) == 3: # Center gap provided
                panel_arc_length += self.gaps[0] - self.gaps[2] # add back width gap and subtract center gap                
            panel_arc_length /= self.number_panels[0]

            # start from the center and work outwards
            positive_x_panels = []
            for i in range(self.number_panels[0] // 2):
                if i == 0:
                    half_center_gap = self.gaps[0] / 2
                    if len(self.gaps) == 3: # Center gap provided
                        half_center_gap = self.gaps[2] / 2
                    x_start = self.determine_end_x_coordinate(0, half_center_gap)
                else:
                    x_start = self.determine_end_x_coordinate(x_end, self.gaps[0])
                
                x_end = self.determine_end_x_coordinate(x_start, panel_arc_length)
                positive_x_panels.append((round(x_start, 5), round(x_end, 5)))

            # create complete list of panels - TODO: This might not be needed because of the SolTrace implementation
            for panel_bounds in reversed(positive_x_panels):
                panel_aperture_bounds.append((-panel_bounds[1], -panel_bounds[0]))
            panel_aperture_bounds.extend(positive_x_panels)
        else:
            raise ValueError("Number of panels in the width direction should be either an even number or equal to 1")

        panel_length = (self.aperture_size[1] - self.gaps[1] * (self.number_panels[1] - 1)) / self.number_panels[1]

        # NOTE: calculating the arc lengths of LS-3 example panels... they are not equal.        
        #arc1 = simpsons_rule(self.parabolic_arc_length_equation, 0.04, 1.62, 1000)
        #arc2 = simpsons_rule(self.parabolic_arc_length_equation, 1.6385, 2.887, 1000)

        # Build trough about the origin then translate to position
        stage = PT.stages[0]        # NOTE: assuming everything is in the first stage (removing stages)
        # Create mirror panels
        self.mirrors.clear()
        panel_y = - self.aperture_size[1] / 2 + panel_length / 2
        for i in range(self.number_panels[1]):          # length
            for j in range(self.number_panels[0]):      # width
                element = stage.add_element()
                element.position = Point(0.0, panel_y, 0.0)
                element.aim = Point(0.0, panel_y, 100)  # point straight up

                if self.number_panels[0] == 1:
                    element.zrot = 0.0
                    element.aperture_rectangle(self.aperture_size[0], panel_length)
                    self.flip_element[element.id] = False
                else:
                    if panel_aperture_bounds[j][0] < 0:
                        # NOTE: This is because aperture single axis curvature requires positive x1, x2 parameters
                        element.zrot = 180.0
                        element.aperture_singleax_curve(-panel_aperture_bounds[j][1], -panel_aperture_bounds[j][0], panel_length)
                        self.flip_element[element.id] = True
                    else:
                        element.zrot = 0.0
                        element.aperture_singleax_curve(panel_aperture_bounds[j][0], panel_aperture_bounds[j][1], panel_length)
                        self.flip_element[element.id] = False

                element.surface_parabolic(self.focal_length, 0.0)
                element.optic = optics[0]
                element.enabled = True

                self.mirrors.append(element)
                
            panel_y += panel_length + self.gaps[1]

        # Adding receiver
        # TODO: Assuming a single receiver per length of aperture - this should be broken up into multiple receiver tubes.
        # Absorber tube
        self.absorbers.clear()
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.focal_length - self.absorber_diameter/2)
        element.aim = Point(0.0, 0.0, 100)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical(self.absorber_diameter/2)
        element.optic = optics[1]
        element.enabled = True
        self.absorbers.append(element)
        self.flip_element[element.id] = False

        # Outer Envelope
        self.envelope.clear()
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.focal_length - self.envelope_diameter/2)
        element.aim = Point(0.0, 0.0, 100)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical(self.envelope_diameter/2)
        element.optic = optics[2]
        element.interaction = 1
        element.enabled = True
        self.envelope.append(element)
        self.flip_element[element.id] = False

        # Inner Envelope
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.focal_length - ((self.envelope_diameter / 2) - self.envelope_thickness))
        element.aim = Point(0.0, 0.0, 100)  # point straight up
        element.zrot = 0.0
        element.aperture_singleax_curve(0.0, 0.0, self.aperture_size[1])
        element.surface_cylindrical((self.envelope_diameter / 2) - self.envelope_thickness)
        element.optic = optics[3]
        element.interaction = 1
        element.enabled = True
        self.envelope.append(element)
        self.flip_element[element.id] = False

        # Azimuth and tilt convert to radians
        tilt_rads = self.tilt * (np.pi / 180.0)
        azimuth_rads = self.azimuth * (np.pi / 180.0)

        # Calculating aperture normal
        norm = Point(0.0, 0.0, 1.0)
        norm = PT.util_rotation_arbitrary( tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), norm)           # Rotation about the x-axis
        norm = PT.util_rotation_arbitrary( -azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), norm)       # Rotation about the z-axis
        self.aperture_normal = norm

        # Calculating the rotation axis
        track_axis = Point(0.0, 1.0, 0.0)
        track_axis = PT.util_rotation_arbitrary( tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), track_axis)           # Rotation about the x-axis
        track_axis = PT.util_rotation_arbitrary( -azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), track_axis)       # Rotation about the z-axis
        self.track_rotation_axis = track_axis
        
        assert abs(np.dot(self.aperture_normal.as_list(), self.track_rotation_axis.as_list())) <= 1.e-6    # This should always be true        

        x_axis = np.cross(self.track_rotation_axis.as_list(), self.aperture_normal.as_list())
        self.x_axis = Point(x_axis[0], x_axis[1], x_axis[2])

        # Rotation matrix from stage to local coordinates
        rotation_matrix = np.array([[self.x_axis.x, self.track_rotation_axis.x, self.aperture_normal.x], 
                                    [self.x_axis.y, self.track_rotation_axis.y, self.aperture_normal.y],
                                    [self.x_axis.z, self.track_rotation_axis.z, self.aperture_normal.z]])
        
        assert abs(np.linalg.det(rotation_matrix) - 1) <= 1.e-8
        euler_angles = compute_euler_angles(np.linalg.matrix_transpose(rotation_matrix))   # NOTE: Stage coordinates are assumed to be Identity matrix
        # NOTE: Transpose of a rotation matrix is equal to its inversve.
        #print("Euler angles (initial) {}".format(euler_angles))

        # Updating element positions, aimpoints, and z-rotations based on tilt and azimuth
        for element in self.get_elements(): 
            # Tilt first then azimuth, negative azimuth due to right-hand rule
            # TODO: We could use the rotation matrix to modify position and aim
            element.position = PT.util_rotation_arbitrary(tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), element.position)
            element.position = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), element.position)

            element.aim = PT.util_rotation_arbitrary(tilt_rads, Point(1,0,0), Point(0.0,0.0,0.0), element.aim)
            element.aim = PT.util_rotation_arbitrary(-azimuth_rads, Point(0,0,1), Point(0.0,0.0,0.0), element.aim)

            element.zrot = 180.0 if self.flip_element[element.id] else 0.0
            element.zrot += euler_angles[2]

            # Translate to specific location
            element.position += self.position
            element.aim += self.position
            
            # Compute the euler angles for the target element -> can be used to compare rotation_matrix and rreftoloc
            # eu_angles = PT.util_calc_euler_angles(np.array(el.position.as_list()), np.array(el.aim.as_list()), el.zrot)
            # # Compute the transform matrix
            # transforms = PT.util_calc_transforms(eu_angles)
            # rreftoloc = transforms['rreftoloc']

    def update_geometry(self, PT:PySolTrace, azimuth: float, elevation: float):
        """
        Update trough geometry based on solar position (azimuth and elevation)

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
        delta_tracking_angle = np.arccos(np.dot(self.aperture_normal.as_list(), sun_proj.as_list()) 
                                   / (np.linalg.norm(self.aperture_normal.as_list()) * np.linalg.norm(sun_proj.as_list())))
        
        # Determine the direction of rotation
        direction = np.cross(self.aperture_normal.as_list(), sun_proj.as_list())
        direction = np.dot(direction, self.track_rotation_axis.as_list())
        if direction < 0:
            delta_tracking_angle *= -1.0
        
        # Update tracking angle and check limits, adjust if exceeds limits
        self.tracking_angle += delta_tracking_angle * 180.0 / np.pi
        if self.tracking_angle < self.tracking_limits[0]:
            delta_tracking_angle += (self.tracking_limits[0] - self.tracking_angle) * np.pi / 180.0
            self.tracking_angle = self.tracking_limits[0]
        elif self.tracking_angle > self.tracking_limits[1]:
            delta_tracking_angle -= (self.tracking_angle - self.tracking_limits[1]) * np.pi / 180.0
            self.tracking_angle = self.tracking_limits[1]
        # print("Tracking angle: {}".format(self.tracking_angle))

        # x_loc_prev = np.array([[self.x_axis.x, self.track_rotation_axis.x, self.aperture_normal.x], 
        #                        [self.x_axis.y, self.track_rotation_axis.y, self.aperture_normal.y],
        #                        [self.x_axis.z, self.track_rotation_axis.z, self.aperture_normal.z]])

        # Update aperture normal
        self.aperture_normal = PT.util_rotation_arbitrary(delta_tracking_angle, self.track_rotation_axis, Point(), self.aperture_normal)

        assert abs(np.dot(self.aperture_normal.as_list(), self.track_rotation_axis.as_list())) <= 1.e-5    # This should always be true        

        # Update x-axis by taking the cross product
        x_axis = np.cross(self.track_rotation_axis.as_list(), self.aperture_normal.as_list())
        self.x_axis = Point(x_axis[0], x_axis[1], x_axis[2])

        x_loc = np.array([[self.x_axis.x, self.track_rotation_axis.x, self.aperture_normal.x], 
                          [self.x_axis.y, self.track_rotation_axis.y, self.aperture_normal.y],
                          [self.x_axis.z, self.track_rotation_axis.z, self.aperture_normal.z]])
        
        assert abs(np.linalg.det(x_loc) - 1) <= 1.e-8
        euler_angles = compute_euler_angles(np.linalg.matrix_transpose(x_loc))   # Stage coordinates is Indentity
        # print("Euler angles (from stage) {}".format(euler_angles))

        # NOTE: we could calculate the rotation matrix, then update position and aim using the matrix directly.
        # rotation_matrix = x_loc @ np.linalg.inv(x_loc_prev)      # rotation from previous state

        # Updating element positions, aimpoints, and z-rotations based on tilt and azimuth
        for element in self.get_elements():
            # test_position = rotation_matrix @ np.array(element.position.as_list())
            # test_aim = rotation_matrix @ np.array(element.aim.as_list())

            element.position = PT.util_rotation_arbitrary(delta_tracking_angle, self.track_rotation_axis, self.position, element.position)
            element.aim = PT.util_rotation_arbitrary(delta_tracking_angle, self.track_rotation_axis, self.position, element.aim)

            element.zrot = 180.0 if self.flip_element[element.id] else 0.0
            element.zrot += euler_angles[2]
            
            # assert abs(test_position[0] - element.position.x) <= 1e-6
            # assert abs(test_position[1] - element.position.y) <= 1e-6
            # assert abs(test_position[2] - element.position.z) <= 1e-6

            # assert abs(test_aim[0] - element.aim.x) <= 1e-6
            # assert abs(test_aim[1] - element.aim.y) <= 1e-6
            # assert abs(test_aim[2] - element.aim.z) <= 1e-6

    def calculate_receiver_power(self, PT: PySolTrace):
        """
        Calculates the power absorbed by the receiver tubes

        :param PT: PySolTrace instance

        :return total absorbed power: [kWt]
        """
        # TODO: We could expand post processing to include receiver losses, optical efficiency
        # Shading and Cosine, Reflection, Intercept, Absorption, Overall Optical Efficiency
        # These calculations would require the full field of as a "Trough System" information.
        rayData = PT.raydata

        absorbed_power = 0.0
        for rec in self.absorbers:
            absorbed_power += len(rayData[rayData['element'] == -(rec.id+1)]) * PT.powerperray / 1.e3

        return absorbed_power
    
def trace(PT: PySolTrace, dni: float = 1000.0, nrays: int = 1e6, plot_trace: bool = False, nthreads: int = 1):
    PT.num_ray_hits = nrays
    PT.max_rays_traced = nrays*200
    PT.is_sunshape = True
    PT.is_surface_errors = True
    PT.dni = dni
    
    res = PT.run(123, nthread=nthreads, no_callback=True)
    if plot_trace:
        PT.plot_trace()

if __name__ == "__main__":

    solar_tracking_test = True
    run_sun_positions_test = False
    run_azi_tilt_test = False

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
    position = Point(20.0, -20.0, 30.0)
    aperture_size = (5.774, 11.96)
    number_panels = (4, 7) # (1, 7)
    gaps = (0.02, 0.01, 0.08)
    focal_length = 1.71
    azimuth = 0.0
    tilt = 0.0
    receiver_dimensions = (0.07, 0.115, 0.003)
    # TODO: figure out handling of multiple receivers

    trough = parabolic_trough(position, aperture_size, number_panels, gaps, focal_length, azimuth, tilt, receiver_dimensions)
    trough.create_geometry(PT, optics)
    trough.update_geometry(PT, 135.0, 60.)
    PT.write_soltrace_input_file('parabolic_trough.stinput')

    trace(PT)
    rec_power = trough.calculate_receiver_power(PT)
    print("Power on the receiver: {:.2f} [kWt]".format(rec_power))

    # Sun tracking - with multiple troughs
    if solar_tracking_test:
        latitude = 34.8639          # Daggett, CA
        altitude = 610.0 / 1000.0   # [km]
        day_of_year = 80

        # Trough set-up
        #positions = [Point(10.0, -13.0, 0.0), Point(10.0, 0.0, 0.0), Point(10.0, 13.0, 0.0)]    # north-south
        positions = [Point(10.0, 0.0, 0.0), Point(0.0, 0.0, 0.0), Point(-10.0, 0.0, 0.0)]    # north-south
        # positions = [Point(-13.0, 10.0, 0.0), Point(0.0, 10.0, 0.0), Point(13.0, 10.0, 0.0)]      # east-west
        aperture_size = (6.0, 12.0)
        number_panels = (4, 7)
        gaps = (0.02, 0.01, 0.08)
        focal_length = 1.71
        azimuth = 0.0
        tilt = 0.0
        receiver_dimensions = (0.07, 0.115, 0.003)

        PT.stages.clear()
        PT.add_stage()
        troughs = list()
        for position in positions:
            trough = parabolic_trough(position, aperture_size, number_panels, gaps, focal_length, azimuth, tilt, receiver_dimensions)
            trough.create_geometry(PT, optics)
            troughs.append(trough)

        # Time-simulation
        rec_power = list()
        dni = list()
        hours = [5 + h*0.5 for h in range(29)]
        #hours = [5 + h for h in range(15)]
        for h in hours:
            print("Simulating hour {:}...".format(h))
            azimuth, elevation = sun_position(latitude, day_of_year, h)    
            zenith = (90.0 - elevation) * np.pi / 180.0
            dni_hour = meinel_clearsky(day_of_year, zenith, altitude)   # Clear sky DNI values
            dni.append(dni_hour)

            for trough in troughs:
                trough.update_geometry(PT, azimuth, elevation)

            if h in [6, 8, 10, 12, 17]:
                PT.write_soltrace_input_file('parabolic_trough_{}_hour.stinput'.format(h))

            if dni_hour > 0.0:
                trace(PT, dni=dni_hour, nthreads=12)
                # Sum power on receivers
                power = 0.0
                for trough in troughs:
                    power += trough.calculate_receiver_power(PT)
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

    # Moves sun and updates tracking
    if run_sun_positions_test:
        position = Point(0.0, 0.0, 0.0)
        aperture_size = (5.774, 11.96)
        number_panels = (2, 2) # (4, 7)   # (1, 7)
        gaps = (0.02, 0.01, 0.08)
        focal_length = 1.71
        azimuth = 0.0
        tilt = 0.0
        receiver_dimensions = (0.07, 0.115, 0.003)

        PT.stages.clear()
        PT.add_stage()
        trough = parabolic_trough(position, aperture_size, number_panels, gaps, focal_length, azimuth, tilt, receiver_dimensions)
        trough.create_geometry(PT, optics)

        sun_azimuths = [70.0, 120.0, 135.0, 180.0, 225.0, 270.0]
        sun_elevation = [10.0, 30.0, 60.0, 90.0]
        for sun_azi in sun_azimuths:
            for sun_el in sun_elevation:
                print("Sun azimuth {}, elevation {}".format(sun_azi, sun_el))
                trough.update_geometry(PT, sun_azi, sun_el)
                print("Tracking angle: {:.2f} [deg]".format(trough.tracking_angle))

                # Moving the sun directly overhead
                #PT.sun.position = trough.aperture_normal * 1000.0
                # TODO: This is not a good approach because it would not capture if the aperture normal was calculated incorrectly
                 
                #PT.write_soltrace_input_file('parabolic_trough_sunAzi_{}_sunEl_{}.stinput'.format(int(sun_azi), int(sun_el)))
                trace(PT)
                rec_power = trough.calculate_receiver_power(PT)
                print("Power on the receiver: {:.2f} [kWt]".format(rec_power))
                if sun_azi == 180.0 and sun_el == 10.0:
                    continue    # Skiping because of very low power
                else: 
                    assert rec_power >= 10.0    # Non-trivial power

    # Tests for azimuth and tilt moving sun directly overhead of trough
    if run_azi_tilt_test:
        position = Point(100.0, -100.0, 0.0)
        aperture_size = (5.774, 11.96)
        number_panels = (4, 7)
        gaps = (0.02, 0.01, 0.08)
        focal_length = 1.71
        receiver_dimensions = (0.07, 0.115, 0.003)

        azimuths = [0.0, 30.0, 45.0, 90.0]
        tilts = [0.0, 15.0, 30.0]
        for azimuth in azimuths:
            for tilt in tilts:
                sun_vec = sun_vector(azimuth + 180.0, 90. - tilt)
                PT.sun.position = sun_vec * 1000.0      # + position # Is not needed because sun position is relative to all geometry

                PT.stages.clear()
                PT.add_stage()
                trough = parabolic_trough(position, aperture_size, number_panels, gaps, focal_length, azimuth, tilt, receiver_dimensions)
                trough.create_geometry(PT, optics)

                trace(PT, nrays=1e5)
                overhead_rec_power = trough.calculate_receiver_power(PT)
                # PT.write_soltrace_input_file('parabolic_trough_azi_{}_tilt_{}.stinput'.format(azimuth,tilt))
                print("Azimuth {:}, Tilt {:}".format(azimuth, tilt))
                print("Sun Position: {}".format(PT.sun.position))
                print("Power on the receiver (sun directly overhead) {:.2f} [kWt]".format(overhead_rec_power))
                assert abs(overhead_rec_power - 54.5) <= 0.2 

    # Testing solar position calculation:
    if False:
        latitude = 39.733056 # [deg]
        day = 182 

        hours = [6 + h*0.5 for h in range(25)]
        sun_x = list()
        sun_y = list()
        sun_z = list()
        azimuths = list()
        elevation = list()
        for h in hours:
            sun_vect = sun_vector_from_latitude(latitude, day, h)
            sun_x.append(sun_vect.x)
            sun_y.append(sun_vect.y)
            sun_z.append(sun_vect.z)

            sun_pos = sun_position(latitude, day, h)
            azimuths.append(sun_pos[0])
            elevation.append(sun_pos[1])

        fig = plt.figure()
        plt.plot(hours, sun_x, label= "x-component")
        plt.plot(hours, sun_y, label= "y-component")
        plt.plot(hours, sun_z, label= "z-component")
        plt.legend()
        plt.show()

        fig = plt.figure()
        plt.plot(hours, azimuths, label= "azimuth")
        plt.plot(hours, elevation, label= "elevation")
        plt.legend()
        plt.show()

    # TODO:
    #   - Add tracking error to troughs
    #   - Add more fidelity to receiver objects
    #   - Create a trough system class that is made of multiple troughs objects
    #       - Calculate field opical parameters
    #   - Add defocus options?

    # NOTE: 
    #   - Inner Envelope does not transmit rays when both sun shape and optical errors are turned off
