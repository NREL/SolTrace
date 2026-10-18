"""
@author: whamilton

Creates a parabolic dish template made of multiple mirror panels.
Dish template updates positioning based on sun position following a tilt-and-rotate two-axis tracking method.
"""

import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..', '..'))

from pysoltrace import PySolTrace, Point
from parabolic_trough import compute_euler_angles, sun_vector, trace, sun_position, meinel_clearsky

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 15})


class parabolic_dish:
    
    def __init__(self, 
                 position: Point, 
                 aperture_size: float, 
                 number_panels: tuple, 
                 gaps: tuple, 
                 focal_length: float, 
                 receiver_dimensions: tuple,
                 ):
        """
        :param position: [m] Position of the parabolic dish
        :param aperture_size: [m] Diameter of the aperture
        :param number_panels: [-] Number of panels (radial, angular)
        :param gaps: [m] Gaps between the panels (radial, angular, center radius*) *optional
        :param focal_length: [m] Focal length of the parabolic dish
        :param receiver_dimenisons: [m] Absorber diameter, distance from dish vertex
        """
        self.position = position

        if aperture_size < 0.0:
            raise ValueError("Aperture diameter must be a positive value")
        self.aperture_size = aperture_size

        if len(number_panels) != 2:
            raise ValueError("Number of panels should be a tuple with 2 elements")
        if not any(isinstance(panels, int) for panels in number_panels):
            raise ValueError("Number of panels must be integers")
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

        if len(receiver_dimensions) != 2:
            raise ValueError("Receiver dimensions parameter should be a tuple with 2 elements")
        if any(rec_dim < 0.0 for rec_dim in receiver_dimensions):
            raise ValueError("Receiver dimensions must be positive values")
        self.absorber_diameter = receiver_dimensions[0]
        self.absorber_distance = receiver_dimensions[1]
        
        # Calculated values
        self.mirrors = list()
        self.absorbers = list()

        self.tracking_elevation = 90.0   # pointing straight up
        self.tracking_azimuth = 180.0     # pointing straight up

        self.elevation_axis = Point(1.0, 0.0, 0.0)

    def get_elements(self):
        return self.mirrors + self.absorbers

    def parabolic_arc_length_equation(self, x: float):
        """
        Equation for integrating the arc length of a 2D parabola.

        :param x: x-coordinate
        :return y: y-coordinate
        """
        c_x = 1 / (2 * self.focal_length)
        return (1 + (c_x * x)**2)**(1/2)

    def determine_end_x_coordinate(self, start:float, distance:float, dx:float = 1e-5) -> float:
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
        Creates the elements of the parabolic dish in the PySolTrace object. 
        NOTE: This assumes the stage has an origin of (0.0, 0.0, 0.0), Aimpoint of (0.0, 0.0, 1.0), and a z-rotation of 0.0 degrees
        
        :param PT: PySolTrace object
        :param optics: list of surface optics (mirror, absorber)
        """
        # Determine the aperture bounds of each panel
        panel_aperture_bounds = []
        radius = self.aperture_size / 2
        if self.number_panels[0] == 1 and len(self.gaps) == 2:
            panel_aperture_bounds = [(0, radius)]
        else:
            # Analytical solution of arc length of parabola
            # Problem: Integral of sqrt(1 + f'(x)^2) dx from -radius to radius
            # where f(x) = 1/2 * c_x * x^2
            c_x = 1 / (2 * self.focal_length)
            arc_length = radius * ((radius * c_x)**2 + 1)**(1/2) + np.arcsinh(radius * c_x) / c_x
            # Confirm result with numerical integration
            check_x = self.determine_end_x_coordinate(-radius, arc_length)
            assert abs(check_x - radius) < 1e-4, "Error in determining the x-coordinate of the parabolic arc length, x = {}".format(check_x)

            # Calculate the arc length the panels
            panel_arc_length = arc_length / 2 - self.gaps[0] * (self.number_panels[0] - 1)
            if len(self.gaps) == 3: # Center gap provided
                # Calculate the arc length of the center gap
                dr = 1e-5   # step size for numerical integration
                center_arc_length = self.gaps[2] # initial guess
                while True:
                    r = self.determine_end_x_coordinate(0.0, center_arc_length)
                    if r < self.gaps[2]:
                        center_arc_length += dr
                    else:
                        break
                assert abs(r - self.gaps[2]) < 1e-4

                panel_arc_length -= center_arc_length # subtract center gap                
            panel_arc_length /= self.number_panels[0]

            # start from the center and work outwards
            for i in range(self.number_panels[0]):
                if i == 0:
                    x_start = self.gaps[0] / 2
                    if len(self.gaps) == 3: # Center gap provided
                        x_start = self.gaps[2]
                else:
                    x_start = self.determine_end_x_coordinate(x_end, self.gaps[0])
                
                x_end = self.determine_end_x_coordinate(x_start, panel_arc_length)
                panel_aperture_bounds.append((round(x_start, 5), round(x_end, 5)))

        # FIXME: There is a limitation on defining panel angular gaps. We cannot mantain a constant gap between panels because
        # panels are defined by an angle from the center of the dish.
        gap_angle = 0.0
        panel_angle = 360.0
        if self.number_panels[1] != 1:
            gap_angle = self.gaps[1] * 360. / (2 * np.pi * radius/2.)   # [deg] Using the half radius of the dish
            panel_angle = (360. - gap_angle * self.number_panels[1]) / self.number_panels[1]
        self.panel_zrot_offset = panel_angle + gap_angle

        # Build dish about the origin then translate to position
        stage = PT.stages[0]        # NOTE: assuming everything is in the first stage (removing stages)
        # Create mirror panels
        self.mirrors.clear()

        z_rot = 0.0
        for i in range(self.number_panels[1]):          # angular
            for j in range(self.number_panels[0]):      # radial
                element = stage.add_element()
                element.position = Point(0.0, 0.0, 0.0)
                element.aim = Point(0.0, 0.0, 100.0)  # point straight up

                if self.number_panels[0] == 1 and self.number_panels[1] == 1 and len(self.gaps) == 2:
                    # Single panel with no center gap
                    element.zrot = 0.0
                    element.aperture_circle(self.aperture_size)
                else:
                    element.zrot = z_rot
                    element.aperture_annulus(panel_aperture_bounds[j][0], panel_aperture_bounds[j][1], panel_angle)

                element.surface_parabolic(self.focal_length, self.focal_length)
                element.optic = optics[0]
                element.enabled = True

                self.mirrors.append(element)
                
            z_rot += self.panel_zrot_offset

        # Absorber surface
        self.absorbers.clear()
        element = stage.add_element()
        element.position = Point(0.0, 0.0, self.absorber_distance)  # distance from the dish vertex
        element.aim = element.position + Point(0.0, 0.0, -1)        # point straight down
        element.zrot = 0.0
        element.aperture_circle(self.absorber_diameter)
        element.surface_flat()
        element.optic = optics[1]
        element.enabled = True
        self.absorbers.append(element)

        # Updating element positions and aimpoints based on global position
        # Translate to specific location
        for element in self.get_elements(): 
            element.position += self.position
            element.aim += self.position

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

        # Calculate updated tracking azimuth and elevation angles
        prev_tracking_azimuth = self.tracking_azimuth
        self.tracking_azimuth = azimuth
        delta_azimuth = (self.tracking_azimuth - prev_tracking_azimuth) * np.pi / 180.0

        prev_tracking_elevation = self.tracking_elevation
        self.tracking_elevation = elevation
        delta_elevation = (self.tracking_elevation - prev_tracking_elevation) * np.pi / 180.0

        # Update elevation axis
        self.elevation_axis = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), Point(), self.elevation_axis)

        # Updating element positions and aimpoints
        # all mirror have the same position and aimpoint - do rotations once up front
        mirror_position = self.mirrors[0].position
        mirror_position = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, mirror_position)
        mirror_position = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, mirror_position)

        mirror_aim = self.mirrors[0].aim
        mirror_aim = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, mirror_aim)
        mirror_aim = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, mirror_aim)

        for element in self.mirrors:
            element.position = mirror_position
            element.aim = mirror_aim
        
        # z-rotation - this ensures the mirrors maintain a tilt-and-rotation orientation
        norm = (mirror_aim - mirror_position).unitize()
        y_axis = np.cross(norm.as_list(), self.elevation_axis.as_list())
        y_axis = Point(y_axis[0], y_axis[1], y_axis[2]).unitize()

        element_el_axis = np.cross(y_axis.as_list(), norm.as_list())
        element_el_axis = Point(element_el_axis[0], element_el_axis[1], element_el_axis[2])

        rotation_matrix = np.array([[element_el_axis.x, y_axis.x, norm.x], 
                                    [element_el_axis.y, y_axis.y, norm.y],
                                    [element_el_axis.z, y_axis.z, norm.z]])
            
        assert abs(np.linalg.det(rotation_matrix) - 1) <= 1.e-8
        euler_angles = compute_euler_angles(np.linalg.matrix_transpose(rotation_matrix)) 
        
        z_rot = euler_angles[2]
        mirror_idx = 0
        for i in range(self.number_panels[1]):          # angular
            for j in range(self.number_panels[0]):      # radial
                element = self.mirrors[mirror_idx]
                element.zrot = z_rot
                mirror_idx += 1
            z_rot += self.panel_zrot_offset

        # Updating absorber position and aimpoint
        for element in self.absorbers:
            element.position = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, element.position)
            element.position = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, element.position)

            element.aim = PT.util_rotation_arbitrary(-delta_azimuth, Point(0.0,0.0,1.0), self.position, element.aim)
            element.aim = PT.util_rotation_arbitrary(-delta_elevation, self.elevation_axis, self.position, element.aim)

            # # Updating each element's z-rotation
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

    def calculate_receiver_power(self, PT: PySolTrace):
        """
        Calculates the power absorbed by the receiver surface

        :param PT: PySolTrace instance

        :return total absorbed power: [kWt]
        """
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

    optics = [mirror_optic, absorber_optic]

    # Set-up the parabolic dish parameters
    position = Point(20.0, -20.0, 30.0)
    aperture_size = 10.0
    number_panels = (4, 7) #(4, 7) # (1, 7)
    gaps = (0.02, 0.01, 0.5)
    focal_length = 7.5
    receiver_dimensions = (0.25, 7.25)

    dish = parabolic_dish(position, aperture_size, number_panels, gaps, focal_length, receiver_dimensions)
    dish.create_geometry(PT, optics)
    PT.write_soltrace_input_file('parabolic_dish.stinput')

    dish.update_geometry(PT, 135.0, 60.)
    PT.write_soltrace_input_file('parabolic_dish_update.stinput')

    dish.update_geometry(PT, 215.0, 25.)
    PT.write_soltrace_input_file('parabolic_dish_update_2.stinput')


    # Sun tracking - with multiple heliostat
    if solar_tracking_test:
        latitude = 34.8639          # Daggett, CA
        altitude = 610.0 / 1000.0   # [km]
        day_of_year = 80

        # Create multi dish system -> clear stages first
        PT.stages.clear()
        PT.add_stage()
        dishs = list()
        # Set-up the parabolic dish parameters
        position = Point(-15.0, 5.0, 10.0)
        aperture_size = 10.0
        number_panels = (4, 7) #(4, 7) # (1, 7)
        gaps = (0.02, 0.01, 0.5)
        focal_length = 7.5
        receiver_dimensions = (0.25, 7.5)
        dish = parabolic_dish(position, aperture_size, number_panels, gaps, focal_length, receiver_dimensions)
        dish.create_geometry(PT, optics)
        dishs.append(dish)

        # Create another dish
        position = Point(0.0, -10.0, 10.0)
        aperture_size = 5.0
        number_panels = (1, 1) #(4, 7) # (1, 7)
        gaps = (0.0, 0.0)
        focal_length = 5.0
        receiver_dimensions = (0.1, 5.0)
        dish = parabolic_dish(position, aperture_size, number_panels, gaps, focal_length, receiver_dimensions)
        dish.create_geometry(PT, optics)
        dishs.append(dish)

        # Create another dish
        position = Point(15.0, 20.0, 10.0)
        aperture_size = 15.0
        number_panels = (2, 6) #(4, 7) # (1, 7)
        gaps = (0.05, 0.05, 0.5)
        focal_length = 10.0
        receiver_dimensions = (0.35, 10.0)
        dish = parabolic_dish(position, aperture_size, number_panels, gaps, focal_length, receiver_dimensions)
        dish.create_geometry(PT, optics)
        dishs.append(dish)

        # Time-simulation
        rec_power = list()
        dni = list()
        dish_power2rec = list()
        for dish in dishs:
            dish_power2rec.append(list())

        # hours = [5 + h*0.5 for h in range(29)]
        hours = [5 + h for h in range(15)]
        for h in hours:
            print("Simulating hour {:}...".format(h))
            azimuth, elevation = sun_position(latitude, day_of_year, h)    
            zenith = (90.0 - elevation) * np.pi / 180.0
            dni_hour = meinel_clearsky(day_of_year, zenith, altitude)   # Clear sky DNI values
            dni.append(dni_hour)

            for dish in dishs:
                dish.update_geometry(PT, azimuth, elevation)

            if h in [7, 10, 12]:
                PT.write_soltrace_input_file('parabolic_dish_{}_hour.stinput'.format(h))

            if dni_hour > 0.0:
                trace(PT, nrays=1e6, dni=dni_hour, nthreads=12)
                
                power = 0.0
                for i, dish in enumerate(dishs):
                    dish_power = dish.calculate_receiver_power(PT)
                    dish_power2rec[i].append(dish_power)
                    power += dish_power

                rec_power.append(power)
                assert power > 50.0    # Non-trivial power
            else:
                rec_power.append(0.0)
                for i, dish in enumerate(dishs):
                    dish_power2rec[i].append(0.0)

        # Plot power and dni over time
        fig = plt.figure()
        plt.plot(hours, rec_power, 'k', label = "Total Power")
        for i, dish in enumerate(dishs):
            plt.plot(hours, dish_power2rec[i], label= "dish = {:}".format(i))
        plt.ylabel('Receiver Power [kWt]')
        plt.xlabel('Hour of Day')
        plt.legend()

        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(hours, dni, 'r')
        ax2.set_ylabel(r"Direct Normal Irradiance [W/m$^2$]", color = 'r')
        ax2.tick_params('y', colors = 'r')
        plt.tight_layout()
        plt.show()
