from selectors import PollSelector
import sys, os
import pandas as pd
from ctypes import *
c_number = c_double   #must be either c_double or c_float depending on copilot.h definition

# Callback to print command line progress messages
@CFUNCTYPE(c_int, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32, c_void_p)
def api_callback(ntracedtotal, ntraced, ntotrace, curstage, nstages, data):
    print("\tProgress: Stage ({:d}/{:d}) - Complete {:.2f}%".format(curstage, nstages, 100.*float(ntraced)/float(ntotrace)) )
    return 1

# ==========================================================================================
class Point:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z 
        return 

# ----------------------------------------------------------------------
class PySolTrace:
    """
    A class to access PySolTrace (SolTrace's Python API)

    Attributes
    ----------
    pdll : class ctypes.CDLL
        loaded SolTrace library of exported functions

    Methods
    -------

    """

    class Optics:
        class Face:
            def __init__(self):
                self.dist_type = 'g'     #One of 'g'->Gauss 'p'->Pillbox 'd'->Diffuse
                self.refraction_real = 1.1         #real component of the refraction index 
                self.reflectivity = 0.96         #reflectivity
                self.transmissivity = 0.         #transmissivity
                self.slope_error = 0.95          #RMS slope error [mrad] 
                self.spec_error = 0.2            #RMS specularity error [mrad]
                self.userefltable = False             #Flag [bool] use reflectivity table 
                self.refltable = []  #[[angle1,refl1],[...]] 

        def __init__(self, p_dll, p_data : int, id : int):
            self.p_dll = p_dll
            self.p_data = p_data #pointer to data context
            
            self.name = "new optic"
            self.id = id 

            self.front = PySolTrace.Optics.Face()
            self.back = PySolTrace.Optics.Face()


        def Create(self) -> int:
            self.p_dll.st_optic.restype = c_int

            dummy_grating = (c_number*3)()
            
            for i,opt in enumerate([self.front, self.back]):

                user_angles = (c_number*len(opt.refltable))()
                user_refls = (c_number*len(opt.refltable))()
                if len(opt.refltable) > 1:
                    user_angles[:] = list(list(zip(*opt.refltable))[0])
                    user_refls[:] = list(list(zip(*opt.refltable))[1])

                # Create surface optic
                self.p_dll.st_optic( \
                    c_void_p(self.p_data),
                    c_uint32(self.id),
                    c_int(i+1),   #front
                    c_wchar(opt.dist_type[0]),
                    c_int(1), #optical surface number
                    c_int(3), #aperture grating
                    c_int(4), #Diffraction order
                    c_number(opt.refraction_real),
                    c_number(0.), #imaginary component of refraction
                    c_number(opt.reflectivity),
                    c_number(opt.transmissivity),
                    pointer(dummy_grating),
                    c_number(opt.slope_error),
                    c_number(opt.spec_error),
                    c_int(1 if opt.userefltable else 0),
                    c_int(len(opt.refltable)),
                    pointer(user_angles),
                    pointer(user_refls),
                    )

            return 1
    # ==========================================================================================
            

    # ==========================================================================================
    # STCORE_API int st_sun(st_context_t pcxt, int point_source, char shape, double sigma_halfwidth);
    # STCORE_API int st_sun_xyz(st_context_t pcxt,  double x, double y, double z );
    # STCORE_API int st_sun_userdata(st_context_t pcxt,  st_uint_t npoints, double angle[], double intensity[]);

    class Sun:
        def __init__(self, pdll, p_data):
            self.pdll = pdll
            self.p_data = p_data 

            self.point_source = False 
            self.shape = 'p'
            self.sigma = 4.65
            self.position = Point()
            self.position.z = 100.
            self.user_intensity_table = []   

        def Create(self):

            self.pdll.st_sun.restype = c_int 
            self.pdll.st_sun_xyz.restype = c_int 

            self.pdll.st_sun(c_void_p(self.p_data), c_int(int(self.point_source)), c_wchar(self.shape[0]), c_number(self.sigma))
            self.pdll.st_sun_xyz(c_void_p(self.p_data), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))

            if len(self.user_intensity_table) > 2:
                user_angles = (c_number*len(self.user_intensity_table))()
                user_ints = (c_number*len(self.user_intensity_table))()
                user_angles[:] = list(list(zip(*self.user_intensity_table))[0])
                user_ints[:] = list(list(zip(*self.user_intensity_table))[1])
                self.pdll.st_sun_userdata.restype = c_int 
                return self.pdll.st_sun_userdata(c_void_p(self.p_data), c_uint32(len(self.user_intensity_table)), pointer(user_angles), pointer(user_ints))

            return 1

    # ==========================================================================================


    # ==========================================================================================
    class Stage:

        class Element:
            # STCORE_API int st_element_enabled(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int enabled);
            # STCORE_API int st_element_xyz(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double x, double y, double z);
            # STCORE_API int st_element_aim(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double ax, double ay, double az);
            # STCORE_API int st_element_zrot(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double zrot);
            # STCORE_API int st_element_aperture(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char ap);
            # STCORE_API int st_element_aperture_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8]);
            # STCORE_API int st_element_surface(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char surf);
            # STCORE_API int st_element_surface_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8]);
            # STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file);
            # STCORE_API int st_element_interaction(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int type); /* 1=refract, 2=reflect */
            # STCORE_API int st_element_optic(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *name);
            def __init__(self, parent_stage, element_id : int):
                self.pdll = parent_stage.pdll
                self.p_data = parent_stage.p_data
                self.stage_id = parent_stage.id
                self.id = element_id
                self.enabled = True 
                self.position = Point()
                self.aim = Point()
                self.aim.z = 1
                self.zrot = 0.
                self.aperture = 'r'
                self.aperture_params = [1, 1, 0, 0, 0, 0, 0, 0]
                self.surface = 'f'
                self.surface_params = [0. for i in range(8)]
                self.surface_file = None 
                self.interaction = 2        #1=refract, 2=reflect
                self.optic = None

            def Create(self) -> int:
                self.pdll.st_element_enabled.restype = c_int
                self.pdll.st_element_xyz.restype = c_int
                self.pdll.st_element_aim.restype = c_int
                self.pdll.st_element_zrot.restype = c_int
                self.pdll.st_element_aperture.restype = c_int
                self.pdll.st_element_aperture_params.restype = c_int
                self.pdll.st_element_surface.restype = c_int
                self.pdll.st_element_surface_params.restype = c_int
                self.pdll.st_element_surface_file.restype = c_int
                self.pdll.st_element_interaction.restype = c_int
                self.pdll.st_element_optic.restype = c_int

                aperture_params = (c_number*len(self.aperture_params))()
                surface_params = (c_number*len(self.surface_params))()
                aperture_params[:] = self.aperture_params
                surface_params[:] = self.surface_params

                self.pdll.st_element_enabled(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(int(self.enabled)));
                self.pdll.st_element_xyz(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.position.x),  c_number(self.position.y), c_number(self.position.z));
                self.pdll.st_element_aim(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.aim.x),  c_number(self.aim.y), c_number(self.aim.z));
                self.pdll.st_element_zrot(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.zrot) );
                self.pdll.st_element_aperture(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.aperture[0]));
                self.pdll.st_element_aperture_params(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(aperture_params));
                self.pdll.st_element_surface(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.surface[0]));
                self.pdll.st_element_surface_params(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(surface_params));
                if self.surface_file:
                    self.pdll.st_element_surface_file(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.surface_file.encode()));
                self.pdll.st_element_interaction(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(self.interaction)); #/* 1=refract, 2=reflect */
                self.pdll.st_element_optic(c_void_p(self.p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.optic.name.encode()));

                return 1
        # -------------------------- end Element class ---------------------------------


        # -----------methods of the 'Stage' class --------------------------------------
        def __init__(self, pdll, p_data : int, id : int):
            """
            """
            self.pdll = pdll
            self.p_data = p_data 
            self.id = id
            self.position = Point()
            self.aim = Point()
            self.aim.z = 1
            self.zrot = 0.
            self.is_virtual = False
            self.is_multihit = True 
            self.is_tracethrough = False 

            self.elements = []
            return 

        def Create(self) -> int:
            # STCORE_API int st_stage_xyz(st_context_t pcxt, st_uint_t idx, double x, double y, double z);
            # STCORE_API int st_stage_aim(st_context_t pcxt, st_uint_t idx, double ax, double ay, double az);
            # STCORE_API int st_stage_zrot(st_context_t pcxt, st_uint_t idx, double zrot);
            # STCORE_API int st_stage_flags(st_context_t pcxt, st_uint_t idx, int virt, int multihit, int tracethrough);

            self.pdll.st_stage_xyz.restype = c_int 
            self.pdll.st_stage_aim.restype = c_int 
            self.pdll.st_stage_zrot.restype = c_int 
            self.pdll.st_stage_flags.restype = c_int 

            self.pdll.st_stage_xyz(c_void_p(self.p_data), c_uint32(self.id), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))
            self.pdll.st_stage_aim(c_void_p(self.p_data), c_uint32(self.id), c_number(self.aim.x), c_number(self.aim.y), c_number(self.aim.z))
            self.pdll.st_stage_zrot(c_void_p(self.p_data), c_uint32(self.id), c_number(self.zrot))
            self.pdll.st_stage_flags(c_void_p(self.p_data), c_uint32(self.id), c_int(int(self.is_virtual)), c_int(int(self.is_multihit)), c_int(int(self.is_tracethrough)))

            return 1;

        # STCORE_API int st_num_elements(st_context_t pcxt, st_uint_t stage);
        def get_num_elements(self)->int:
            self.pdll.st_num_elements.restype = c_int
            return self.pdll.st_num_elements(c_void_p(self.p_data), c_uint32(self.id))

        # STCORE_API int st_add_element(st_context_t pcxt, st_uint_t stage);
        # STCORE_API int st_add_elements(st_context_t pcxt, st_uint_t stage, st_uint_t num);
        def add_elements(self, n_to_add = 1) -> int:

            for i in range(n_to_add):
                self.elements.append(  PySolTrace.Stage.Element(self, len(self.elements) ) )

            if n_to_add > 1:
                self.pdll.st_add_elements.restype = c_int 
                return self.pdll.st_add_elements(c_void_p(self.p_data), c_uint32(self.id), c_int(n_to_add))
            else:
                self.pdll.st_add_element.restype = c_int 
                return self.pdll.st_add_element(c_void_p(self.p_data), c_uint32(self.id))

        # STCORE_API int st_delete_element(st_context_t pcxt, st_uint_t stage, st_uint_t idx);
        # STCORE_API int st_clear_elements(st_context_t pcxt, st_uint_t stage);

        

    # ---------- methods of the PySolTrace class --------------------------------------------
    def __init__(self):
        cwd = os.getcwd()
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            self.pdll = CDLL(cwd + "/coretrace.dll")
            print("Loaded win32")
            #self.pdll = CDLL(cwd + "/coretraced.dll") # for debugging
        elif sys.platform == 'darwin':
            self.pdll = CDLL(cwd + "/coretrace.dylib")  # Never tested
        elif sys.platform.startswith('linux'):
            self.pdll = CDLL(cwd +"/coretrace.so")  # Never tested
        else:
            print( 'Platform not supported ', sys.platform)

        self.p_data = None

        self.optics = []
        self.stages = []

        self.num_ray_hits = int(1e5)
        self.max_rays_traced = self.num_ray_hits*100
        self.is_sunshape = True 
        self.is_surface_errors = True

    #     return

    # def data_create(self) -> int:
        """Creates an instance of SolarPILOT in memory

        Returns
        -------
        int
            memory address of SolarPILOT instance 
        """

        self.pdll.st_create_context.restype = c_void_p
        self.p_data = self.pdll.st_create_context()
        # return self.p_data

    # def data_free(self, p_data : int) -> bool:
    def clear_context(self,*args):
        """Frees SolarPILOT instance from memory

        Parameters
        ----------
        p_data : int
            memory address of SolarPILOT instance 

        Returns
        -------
        bool
            True if successful, False otherwise
        """

        self.pdll.st_free_context.restype = c_bool
        self.pdll.st_free_context(c_void_p(self.p_data))

        return True

    def get_num_optics(self) -> int:
        self.pdll.st_num_optics.restype = c_int
        return self.pdll.st_num_optics(c_void_p(self.p_data))


    def add_optic(self, optic_name : str) -> int:
        """
        """

        new_opt_id = len(self.optics)

        self.optics.append( PySolTrace.Optics(self.pdll, self.p_data, new_opt_id ) )
        self.optics[-1].name = optic_name

        self.pdll.st_add_optic.restype = c_int
        self.pdll.st_add_optic(c_void_p(self.p_data), c_char_p(optic_name.encode()))

        return self.optics[-1]

    def delete_optic(self, optic_id : int) -> int:
        # find the appropriate optic
        for opt in self.optics: 
            if optic_id == opt.id:
                # clear it from the optics array 
                self.optics.remove(opt)
                # Remove from the soltrace context
                self.pdll.st_delete_optic.restype = c_int 
                return self.pdll.st_delete_optic(c_void_p(self.p_data), c_uint32(optic_id) )
        
        # If reaching this point, the optic id was not found
        return -1

    def add_sun(self):
        self.sun = PySolTrace.Sun(self.pdll, self.p_data)
        return self.sun


    # STCORE_API int st_num_stages(st_context_t pcxt);
    def get_num_stages(self) -> int:
        self.pdll.st_num_stages.restype = c_int
        return self.pdll.st_num_stages(c_void_p(self.p_data))

    # STCORE_API int st_add_stage(st_context_t pcxt);
    # STCORE_API int st_add_stages(st_context_t pcxt, st_uint_t num);
    def add_stage(self) -> int:
        """
        """
        new_st_id = len(self.stages)

        self.stages.append( PySolTrace.Stage(self.pdll, self.p_data, new_st_id ) )

        self.pdll.st_add_stage.restype = c_int
        self.pdll.st_add_stage(c_void_p(self.p_data) )

        return self.stages[-1]

    # STCORE_API int st_delete_stage(st_context_t pcxt, st_uint_t idx);
    # STCORE_API int st_clear_stages(st_context_t pcxt);
    def delete_stage(self, stage_id : int) -> int:
        # find the appropriate optic
        for st in self.stages: 
            if stage_id == st.id:
                # clear it from the optics array 
                self.stages.remove(st)
                # Remove from the soltrace context
                self.pdll.st_delete_stage.restype = c_int 
                return self.pdll.st_delete_stage(c_void_p(self.p_data), c_uint32(stage_id) )
        
        # If reaching this point, the stage id was not found
        return -1

    # STCORE_API int st_sim_params(st_context_t pcxt, int raycount, int maxcount);
    # def sim_params(self, raycount : int, maxcount : int) -> int:
    # STCORE_API int st_sim_errors(st_context_t pcxt, int include_sun_shape, int include_optics);
    # STCORE_API int st_sim_run( st_context_t pcxt, unsigned int seed, bool AsPowerTower, 
    # int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data), void *data);
    def run(self, seed : int = -1, as_power_tower = False):

        self.pdll.st_sim_errors.restype = c_int 
        self.pdll.st_sim_errors(c_void_p(self.p_data), c_int(1 if self.is_sunshape else 0), c_int(1 if self.is_surface_errors else 0))

        self.pdll.st_sim_params.restype = c_int 
        self.pdll.st_sim_params(c_void_p(self.p_data), c_int(int(self.num_ray_hits)), c_int(int(self.max_rays_traced)))

        self.pdll.st_sim_run.restype = c_int 
        self.pdll.st_sim_run( c_void_p(self.p_data), c_uint16(seed), c_bool(as_power_tower), api_callback)


    # STCORE_API int st_num_intersections(st_context_t pcxt);
    def get_num_intersections(self) -> int:
        self.pdll.st_num_intersections.restype = c_int 
        return self.pdll.st_num_intersections(c_void_p(self.p_data))

    # STCORE_API int st_locations(st_context_t pcxt, double *loc_x, double *loc_y, double *loc_z);
    def get_intersect_locations(self):
        n_int = self.get_num_intersections()

        loc_x = (c_number*n_int)()
        loc_y = (c_number*n_int)()
        loc_z = (c_number*n_int)()

        self.pdll.st_locations.restype = c_int 
        self.pdll.st_locations(c_void_p(self.p_data), pointer(loc_x), pointer(loc_y), pointer(loc_z))

        return list(map(list,zip(loc_x, loc_y, loc_z)))

    # STCORE_API int st_cosines(st_context_t pcxt, double *cos_x, double *cos_y, double *cos_z);
    def get_intersect_cosines(self):
        n_int = self.get_num_intersections()

        cos_x = (c_number*n_int)()
        cos_y = (c_number*n_int)()
        cos_z = (c_number*n_int)()

        self.pdll.st_cosines.restype = c_int 
        self.pdll.st_cosines(c_void_p(self.p_data), pointer(cos_x), pointer(cos_y), pointer(cos_z))

        return list(map(list,zip(cos_x, cos_y, cos_z)))

    # STCORE_API int st_elementmap(st_context_t pcxt, int *element_map);
    def get_intersect_elementmap(self):
        n_int = self.get_num_intersections()
        element_map = (c_int*n_int)()

        self.pdll.st_elementmap.restype = c_int 
        self.pdll.st_elementmap(c_void_p(self.p_data), pointer(element_map))

        return list(element_map)

    # STCORE_API int st_stagemap(st_context_t pcxt, int *stage_map);
    def get_intersect_stagemap(self):
        n_int = self.get_num_intersections()
        stage_map = (c_int*n_int)()

        self.pdll.st_stagemap.restype = c_int 
        self.pdll.st_stagemap(c_void_p(self.p_data), pointer(stage_map))

        return list(stage_map)

    # STCORE_API int st_raynumbers(st_context_t pcxt, int *ray_numbers);
    def get_intersect_raynumbers(self):
        n_int = self.get_num_intersections()
        raynumbers = (c_int*n_int)()

        self.pdll.st_raynumbers.restype = c_int 
        self.pdll.st_raynumbers(c_void_p(self.p_data), pointer(raynumbers))

        return list(raynumbers)

    # STCORE_API int st_sun_stats(st_context_t pcxt, double *xmin, double *xmax, double *ymin, double *ymax, int *nsunrays );
    def get_sun_stats(self):
        xmin = c_number
        xmax = c_number
        ymin = c_number
        ymax = c_number
        nsunrays = c_int

        self.pdll.st_sun_stats.restype = c_int 
        self.pdll.st_sun_stats(c_void_p(self.p_data), pointer(xmin), pointer(xmax), pointer(ymin), pointer(ymax), pointer(nsunrays))

        return {
            'xmin':float(xmin),
            'xmax':float(xmax),
            'ymin':float(ymin),
            'ymax':float(ymax),
            'nsunrays':int(nsunrays),
        }

    def get_ray_dataframe(self, intersect_only=False):
        
        data = {}

        n_int = self.get_num_intersections()
        data['loc_x'] = (c_number*n_int)()
        data['loc_y'] = (c_number*n_int)()
        data['loc_z'] = (c_number*n_int)()

        self.pdll.st_locations.restype = c_int 
        self.pdll.st_locations(c_void_p(self.p_data), pointer(data['loc_x']), pointer(data['loc_y']), pointer(data['loc_z']))

        data['cos_x'] = (c_number*n_int)()
        data['cos_y'] = (c_number*n_int)()
        data['cos_z'] = (c_number*n_int)()

        self.pdll.st_cosines.restype = c_int 
        self.pdll.st_cosines(c_void_p(self.p_data), pointer(data['cos_x']), pointer(data['cos_y']), pointer(data['cos_z']))

        data['element_map'] = (c_int*n_int)()

        self.pdll.st_elementmap.restype = c_int 
        self.pdll.st_elementmap(c_void_p(self.p_data), pointer(data['element_map']))

        data['stage_map'] = (c_int*n_int)()

        self.pdll.st_stagemap.restype = c_int 
        self.pdll.st_stagemap(c_void_p(self.p_data), pointer(data['stage_map']))

        data['raynumbers'] = (c_int*n_int)()

        self.pdll.st_raynumbers.restype = c_int 
        self.pdll.st_raynumbers(c_void_p(self.p_data), pointer(data['raynumbers']))

        for key in data.keys():
            data[key] = list(data[key])

        return pd.DataFrame(data)


    # /* utility transform/math functions */
    def util_calc_euler_angles(self, origin, aimpoint, zrot):

        a_origin = (c_number*3)()
        a_aimpoint = (c_number*3)()
        a_euler = (c_number*3)()
        a_origin[:] = origin
        a_aimpoint[:] = aimpoint

        self.pdll.st_calc_euler_angles.restype = c_void_p
        self.pdll.st_calc_euler_angles(pointer(a_origin), pointer(a_aimpoint), c_number(zrot), pointer(a_euler))

        return list(a_euler)

    # STCORE_API void st_transform_to_local( double posref[3], double cosref[3], double origin[3], double rreftoloc[3][3], double posloc[3], double cosloc[3]);
    def util_transform_to_local(self, posref, cosref, origin, rreftoloc):
        a_posref = (c_number*3)()
        a_cosref = (c_number*3)()
        a_origin = (c_number*3)() 
        a_rreftoloc = (c_number*9)() 
        #output
        posloc = (c_number*3)()
        cosloc = (c_number*3)()
        
        a_posref[:] = posref
        a_cosref[:] = cosref
        a_origin[:] = origin
        a_rreftoloc[:] = sum([], rreftoloc)

        self.pdll.st_transform_to_local.restype = c_void_p
        self.pdll.st_transform_to_local(pointer(a_posref), pointer(a_cosref), pointer(a_origin), pointer(a_rreftoloc), pointer(posloc), pointer(cosloc))

        return {'cosloc':list(cosloc), 'posloc':list(posloc)}

    # STCORE_API void st_transform_to_reference( double posloc[3], double cosloc[3], double origin[3], double rloctoref[3][3], double posref[3], double cosref[3]);
    def util_transform_to_reference(self, posloc, cosloc, origin, rloctoref):
        a_posloc = (c_number*3)()
        a_cosloc = (c_number*3)()
        a_origin = (c_number*3)() 
        a_rloctoref = (c_number*9)() 
        #output
        posref = (c_number*3)()
        cosref = (c_number*3)()
        
        a_posloc[:] = posloc
        a_cosloc[:] = cosloc
        a_origin[:] = origin
        a_rloctoref[:] = sum([], rloctoref)

        self.pdll.st_transform_to_reference.restype = c_void_p
        self.pdll.st_transform_to_reference(pointer(a_posloc), pointer(a_cosloc), pointer(a_origin), pointer(a_rloctoref), pointer(posref), pointer(cosref))

        return {'cosref':list(cosref), 'posref':list(posref)}

    # STCORE_API void st_matrix_vector_mult( double m[3][3], double v[3], double mxv[3] );
    def util_matrix_vector_mult(self, m, v):
        """
        Arguments:
            m[3][3] - a 3x3 list
            v[3] - a list, length 3

        Returns:
            m x v [3]
        """

        a_m = (c_number*9)()
        a_v = (c_number*3)()
        a_mv = (c_number*3)()
        a_m[:] = sum([], m)
        a_v[:] = v 

        self.pdll.st_matrix_vector_mult.restype = c_void_p
        self.pdll.st_matrix_vector_mult(pointer(a_m), pointer(a_v), pointer(a_mv))

        return list(a_mv)

    # STCORE_API void st_calc_transform_matrices( double euler[3], double rreftoloc[3][3], double rloctoref[3][3] );
    def util_calc_transform_matrices(self, euler):
        """
        input:  Euler = Euler angles
        output: RRefToLoc = Transformation matrix from Reference to Local system
        RLocToRef = ""             ""     ""   Local to Reference system (transpose of above)}
        """

        a_euler = (c_number*3)()
        rreftoloc = (c_number*9)()
        rloctoref = (c_number*9)()

        a_euler[:] = euler

        self.pdll.st_calc_transform_matrices.restype = c_void_p
        self.pdll.st_calc_transform_matrices(pointer(a_euler), pointer(rreftoloc), pointer(rloctoref))

        # reshape
        a_rreftoloc = []
        a_rloctoref = []
        for i in range(0,10,3):
            a_rreftoloc.append(rreftoloc[i:i+3])
            a_rloctoref.append(rloctoref[i:i+3])

        return {'rreftoloc':a_rreftoloc, 'rloctoref':a_rloctoref}

    # STCORE_API void st_matrix_transpose( double input[3][3], double output[3][3] );
    def util_matrix_transpose(self, m):

        a_m = (c_number*9)()
        a_m[:] = sum([], m)
        output = (c_number*9)()

        self.pdll.st_matrix_transpose.restype = c_void_p
        self.pdll.st_matrix_transpose(pointer(a_m), pointer(output))

        a_output = []
        for i in range(0,10,3):
            a_output.append(output[i:i+3])

        return a_output

        