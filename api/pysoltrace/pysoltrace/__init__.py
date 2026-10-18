__version__ = "0.0.1"

# find chedder, either in pyshared during development or in pysoltrace.chedder during production
import os, sys
from pathlib import Path


# try:
from pysoltrace.chedder import dot_h, found_in # pyright: ignore[reportMissingImports]
# except (ImportError, ModuleNotFoundError):
#     pyshared_path = Path(__file__).parent.parent / 'pyshared'
#     sys.path.insert(1, str(pyshared_path))
#     from pyshared.chedder import dot_h, found_in

# no outside of pysoltrace dependencies
from pysoltrace.point import Point
from pysoltrace import soltrace_json
# depend on other pysoltrace modules
from pysoltrace import math_utils # Point
from pysoltrace import soltrace_constants # dot_h
from pysoltrace import soltrace_structs # dot_h
from pysoltrace.api import api, STAPIv2Exception # soltrace_constants, soltrace_structs

from pysoltrace.legacy import legacy as PySolTrace # api, dot_h, soltrace_json, math_utils, Point

from pysoltrace import cst_templates as cst # PySolTrace

__all__ = [
    'dot_h',
    'found_in',
    'math_utils',
    'Point',
    'PySolTrace',
    'soltrace_constants',
    'soltrace_json',
    'soltrace_structs',
    'api',
    'STAPIv2Exception',
    'cst'
]

"""
Create context handle for API                      -> st_create_context
Reset data at handle                               -> st_reset_context
Free handle                                        -> st_free_context
Import SolTrace data from string                   -> st_read_input_json
Import SolTrace data from json file                -> st_read_input_json_file
Set all simulation parameters                      -> st_set_simulation_parameters
Set ray counts and power tower parameters (legacy) -> st_sim_params
Set ray counts parameters                          -> st_sim_rays
Set power tower parameter                          -> st_sim_power_tower
Set sun and optical error parameters               -> st_sim_errors
Set latitude and longitude parameters              -> st_sim_location
Set convergence tolerence parameter                -> st_sim_tolerance
Get current parameters                             -> st_get_simulation_parameters
Get number of optical property sets                -> st_num_optics
Add optical property set                           -> st_add_optical_properties_set
Get optical property set                           -> st_get_optical_properties_set
Delete optical property set                        -> st_delete_optic
Delete all optical property sets                   -> st_clear_optics
Get number of elements                             -> st_num_elements
Add element                                        -> st_add_element
Get element                                        -> st_get_element
Delete element                                     -> st_delete_element
Delete all elements                                -> st_clear_elements
Set if element is enabled                          -> st_element_enabled
Set if element is virtual                          -> st_element_virtual
Set element origin/position                        -> st_element_xyz
Set element aim point                              -> st_element_aim
Set element z rotation                             -> st_element_zrot
Set element aperture                               -> st_element_aperture
Set element surface                                -> st_element_surface
Set element optical property set                   -> st_element_optic
Set element group                                  -> st_element_group
Add sun                                            -> st_add_sun
Get sun                                            -> st_get_sun
Set sun shape                                      -> st_sun_shape
Set sun position                                   -> st_sun_xyz
Get sun position (legacy)                          -> st_sun_position
Set user-defined sun                               -> st_sun_userdata
Get sun azimuth/zenith angles                      -> st_get_sun_az_zen
Get sun azimuth/elevation angles                   -> st_get_sun_az_el
Get sun position                                   -> st_get_sun_vector
Write SolTrace JSON file                           -> st_export_json_file
Get which runners are installed (contextless)      -> st_get_installed_runners
Get if a runner is installed (contextless)         -> st_is_runner_installed
Set up SolTrace runner                             -> st_sim_setup
Run SolTrace simulation                            -> st_sim_run_v2
Report SolTrace simulation                         -> st_sim_report
Write ray record result to CSV                     -> st_write_results_csv
Write group results to JSON                        -> st_write_group_results_json
Get number of intersections in simulation          -> st_num_intersections
Get ray intersection locations                     -> st_locations
Get ray intersection directions                    -> st_cosines
Get elements intersected                           -> st_elementmap
Get stages intersected (legacy)                    -> st_stagemap
Get ray ids                                        -> st_raynumbers
Get sun plane size and number of generated rays    -> st_sun_stats
Get struct of results                              -> st_get_results_data
Use an array to call API functions                 -> st_batch
"""