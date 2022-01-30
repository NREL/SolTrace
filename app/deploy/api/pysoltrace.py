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

    Attributes:
    ----------
    pdll            | (ctypes.CDLL) loaded SolTrace library of exported functions
    p_data          | Memory location for soltrace context
    optics          | ([PySolTrace.Optics,]) List of Optics instances
    stages          | ([PySolTrace.Stage,]) List of Stage instances
    num_ray_hits    | (int) Minimum number of simulation ray hits 
    max_rays_traced | (int) Maximum number of ray hits in a simulation
    is_sunshape     | (bool) Flag indicating whether sunshape should be included
    is_surface_errors| (bool) Flag indicating whether surface errors should be included

    Methods:
    -------
    clear_context   | Frees SolTrace instance from memory at p_data
    get_num_optics  | Get number of optical elements in the SolTrace context
    add_optics      | Instantiates a new PySolTrace.Optics object
    delete_optic    | Delete Optics instance
    add_sun         | Adds Sun instance
    get_num_stages  | Get number of Stages
    add_stage       | Adds Stage instance
    validate        | Detect common errors in simulation setup 
    run             | Runs SolTrace simulation
    --------------------------------------------------------------------------
    get_num_intersections   | Get the number of simulation ray intersections
    get_ray_dataframe       | Get all simulation ray data in pandas dataframe
    get_intersect_locations | Get ray intersection locations
    get_intersect_cosines   | Get ray intersection cosines / direction vectors
    get_intersect_elementmap| Get ray intersection associated element
    get_intersect_stagemap  | Get ray intersection associated stage
    get_intersect_raynumbers| Get ray intersection number
    get_sun_stats           | Get information on sun box
    --------------------------------------------------------------------------
    util_calc_euler_angles  | Calculate Euler angles for a position, aimpoint, and rotation
    util_transform_to_local | Transform a coordinate system from reference to a local system
    util_transform_to_ref   | Transform a coordinate system from local to reference system
    util_matrix_vector_mult | Calculate product of a square matrix and a vector
    util_calc_transforms    | Calculate matrix transforms
    util_matrix_transpose   |
    """

    class Optics:
        """ *Optics* is a subclass of PySolTrace, and represents an optical property set. 
        A PySolTrace instance may have multiple Optics member instances, which are stored in 
        the PySolTrace.optics list.
        
        Optics contains a subclass *Face*, which collects properties associated with the front 
        or back face of an optical surface.

        Attributes:
            p_dll       | Reference for API DLL, managed by PySolTrace
            p_data      | Memory location for soltrace context, managed by PySolTrace
            name        | (str) Unique name for the optical property set
            id          | (int) Identifying integer associated with the property set
            front       | (Face) properties associated with the front of the optical surface
            back        | (Face) properties associated with the back of the optical surface

        Methods:
            Create      | Calls methods to instantiate and construct optical surface in the
                        | SolTrace context.
        """
        class Face:
            """
            Subclass of Optics, contains properties associated with one of the optical faces.

            Attributes:
                dist_type       | (char) Distribution type for surface interactions. One of:
                                | {'g':Gaussian, 'p':Pillbox, 'd':Diffuse }
                refraction_real | (float) Real component of the refraction index
                reflectivity    | (float) [0..1] Surface reflectivity
                transmissivity  | (float) [0..1] Surface transmissivity
                slope_error     | (float) [mrad] Surface RMS slope error, half-angle
                spec_error      | (float) [mrad] Surface specularity error, half-angle
                userefltable    | (bool) Flag specifying use of user reflectivity table to  
                                |        modify reflectivity as a function of incidence angle
                refltable       | ([[float,float],]) [mrad,0..1] 2D list containing pairs of
                                |        [angle,reflectivity] values.              
            """
            def __init__(self):
                self.dist_type = 'g'     #One of 'g'->Gauss 'p'->Pillbox 'd'->Diffuse
                self.refraction_real = 1.1         #real component of the refraction index 
                self.reflectivity = 0.96         #reflectivity
                self.transmissivity = 0.         #transmissivity
                self.slope_error = 0.95          #RMS slope error [mrad] 
                self.spec_error = 0.2            #RMS specularity error [mrad]
                self.userefltable = False             #Flag [bool] use reflectivity table 
                self.refltable = []  #[[angle1,refl1],[...]] 

        # -------- methods of the Optics class -----------------------------------------
        def __init__(self, p_dll, p_data : int, id : int):
            self.p_dll = p_dll
            self.p_data = p_data #pointer to data context
            
            self.name = "new optic"
            self.id = id 

            self.front = PySolTrace.Optics.Face()
            self.back = PySolTrace.Optics.Face()

        def Create(self) -> int:
            """
            Create Optics instance in the SolTrace context.

            Returns:
                int     | 1 if successful, 0 otherwise
            """
            self.p_dll.st_optic.restype = c_int

            dummy_grating = (c_number*3)()

            resok = True
            
            # for each face -- front or back
            for i,opt in enumerate([self.front, self.back]):

                user_angles = (c_number*len(opt.refltable))()
                user_refls = (c_number*len(opt.refltable))()
                if len(opt.refltable) > 1:
                    user_angles[:] = list(list(zip(*opt.refltable))[0])
                    user_refls[:] = list(list(zip(*opt.refltable))[1])

                # Create surface optic
                resok = resok and self.p_dll.st_optic( \
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

            return 1 if resok else 0
    # ========end Optics class =================================================================
            

    # ==========================================================================================
    class Sun:
        """ *Sun* is a subclass of PySolTrace, and represents a sun property set. 
        A PySolTrace instance may have a single Sun member instance, which is stored as the 
        PySolTrace.sun member.

        Attributes:
            p_dll           | Reference for API DLL, managed by PySolTrace
            p_data          | Memory location for soltrace context, managed by PySolTrace
            point_source    | (bool) Flag indicating whether the sun is modeled as a point source
                            |   at a finite distance.
            shape           | (char) Sun shape model. One of:
                            |   {'p':Pillbox, 'g':Gaussian, 'd':data table, 'f':gray diffuse}
            sigma           | (float) [mrad] Half-width or std. dev. of the error distribution
            position        | (Point) Location of the sun/sun vector in global coordinates
            user_intensity..|
              .._table      | ([[float,float],]) [mrad, 0..1] 2D list containing pairs of
                            |   angle deviation from sun vector and irradiation intensity.
                            |   A typical table will have angles spanning 0->~5mrad, and inten-
                            |   sities starting at 1 and decreasing to zero. The table must
                            |   contain at least 2 entries.    
        Methods:
            Create          | Calls methods to instantiate and construct optical surface in the
                            |   SolTrace context.

        """
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
            """
            Create Sun instance in the SolTrace context.

            Returns:
                int     | 1 if successful, 0 otherwise
            """

            self.pdll.st_sun.restype = c_int 
            self.pdll.st_sun_xyz.restype = c_int 

            self.pdll.st_sun(c_void_p(self.p_data), c_int(int(self.point_source)), c_wchar(self.shape[0]), c_number(self.sigma))
            self.pdll.st_sun_xyz(c_void_p(self.p_data), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))

            # If a user intensity table is provided, and the shape is specified accordingly as 'd', load the data table into context
            if len(self.user_intensity_table) > 2 and self.shape.lower()[0] == 'd':
                user_angles = (c_number*len(self.user_intensity_table))()
                user_ints = (c_number*len(self.user_intensity_table))()
                user_angles[:] = list(list(zip(*self.user_intensity_table))[0])
                user_ints[:] = list(list(zip(*self.user_intensity_table))[1])

                self.pdll.st_sun_userdata.restype = c_int 
                return self.pdll.st_sun_userdata(c_void_p(self.p_data), c_uint32(len(self.user_intensity_table)), pointer(user_angles), pointer(user_ints))

            return 1
    # ===========end of the Sun class===========================================================

    # ==========================================================================================
    class Stage:
        """ *Stage* is a subclass of PySolTrace, and represents a grouping of elements. 
        A PySolTrace instance may have multiple Stage member instances, which are stored in 
        the PySolTrace.stages list.

        Stage contains a subclass *Element*, which collects properties and geometry associated
        with individual geometric elements.

        Attributes:
            p_dll           | Reference for API DLL, managed by PySolTrace
            p_data          | Memory location for soltrace context, managed by PySolTrace
            id              | (int) Identifying integer associated with the stage
            position        | (Point) Stage location in global coordinates 
            aim             | (Point) Coordinate system aim point in global coordinates
            zrot            | (float) [deg] Rotation of coordinate system around z-axis
            is_virtual      | (bool) Flag indicating virtual stage
            is_multihit     | (bool) Flag indicating that rays can have multiple interactions
                            |   within a single stage.
            is_tracethrough | (bool) Flag indicating the stage is in trace-through mode
            elements        | ([Stage.Element,]) list of all elements in the stage
            
        Methods:
            Create          | Calls methods to instantiate and construct a stage in the context.
            get_num_elements| Returns number of elements associated with this stage (from context).
            add_elements    | Creates new element in Stage.element[] list, and adds to context
        """

        class Element:
            """ *Element* is a subclass of PySolTrace.Stage, and represents a set of properties and
            geometric settings related to a single optical element in SolTrace. Elements are associ-
            ated with a single stage, and are stored in the respective stage's Stage.elements[] list.

            Attributes:
                p_dll           | Reference for API DLL, managed by PySolTrace
                p_data          | Memory location for soltrace context, managed by PySolTrace
                stage_id        | (int) Identifying integer associated with the containing stage
                id              | (int) Identifying integer associated with element
                enabled         | (bool) Flag indicating whether the element is included in the model
                position        | (Point) Element location in stage coordinates 
                aim             | (Point) Element coordinate system aim point in stage coordinates
                zrot            | (float) [deg] Rotation of coordinate system around z-axis
                interaction     | (int) Flag indicating optical interaction type. {1:refraction, 2:reflection}
                optic           | (Optics) Reference to *Optics* instance associated with this element
                aperture        | (char) Charater indicating aperture type. One of:
                                | {'c':circle, 'h':hexagon, 't':triangle, 'r':rectangle, 'a':annulus,
                                |  'l':single-axis curvature, 'i':irregular triangle, 'q':quadrilateral}
                aperture_params | ([float,]) Up to 8 coefficients defining aperture -- values depend on 
                                |   selection for 'aperture'
                surface         | (char) Character indicating surface type. One of:
                                | {'s':spherical, 'p':parabolic, 'f':flat plane, 'o':hyperboloid/ellipsoid,
                                |  'c':conical, 't':cylindrical, 'd':toroid, 'm':Zernicke monomial, 
                                |  'r':Polynomial revolution, 'i':cubic spline interpolation, 
                                |  'e':finite element data, 'v':VSHOT data}
                surface_params  | ([float,]) Up to 8 coefficients defining surface -- values depend on 
                                |   selection for 'surface'
                surface_file    | (string) Name for surface file, if using compatible type. File extension:
                                |       *.mon --> 'm' / Zernicke monomial 
                                |       *.sht --> 'v' / VSHOT data
                                |       *.ply --> 'r' / Polynomial revolution
                                |       *.csi --> 'i' / Cubic spline interpolation
                                |       *.fed --> 'e' / Finite element data
            Methods:
                Create          | Calls methods to instantiate and construct element in the SolTrace context
            """
            
            # STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file);
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
                """
                Create Element instance in the SolTrace context.

                Returns:
                    int     | 1 if successful, 0 otherwise
                """
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
            """
            Create Stage instance in the SolTrace context. 
            Note: This does not create any associated Elements, which must have their Create method called separately. 

            Returns:
                int     | 1 if successful, 0 otherwise
            """

            self.pdll.st_stage_xyz.restype = c_int 
            self.pdll.st_stage_aim.restype = c_int 
            self.pdll.st_stage_zrot.restype = c_int 
            self.pdll.st_stage_flags.restype = c_int 

            self.pdll.st_stage_xyz(c_void_p(self.p_data), c_uint32(self.id), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))
            self.pdll.st_stage_aim(c_void_p(self.p_data), c_uint32(self.id), c_number(self.aim.x), c_number(self.aim.y), c_number(self.aim.z))
            self.pdll.st_stage_zrot(c_void_p(self.p_data), c_uint32(self.id), c_number(self.zrot))
            self.pdll.st_stage_flags(c_void_p(self.p_data), c_uint32(self.id), c_int(int(self.is_virtual)), c_int(int(self.is_multihit)), c_int(int(self.is_tracethrough)))

            return 1;

        def get_num_elements(self)->int:
            """
            Retrieve the number of elements associated with this stage that are currently created 
            in the SolTrace context. 

            Returns:
                (int) Number of elements
            """
            self.pdll.st_num_elements.restype = c_int
            return self.pdll.st_num_elements(c_void_p(self.p_data), c_uint32(self.id))

        def add_elements(self, n_to_add = 1) -> int:
            """
            Add one or more elements to the stage. This method appends an Element object to the
            stage's Stage.elements list and adds the element(s) to the SolTrace context. 
            To update element properties and settings, call the Element.Create method on each element.

            Arguments:
                n_to_add <= 1>      | (int) Number of elements to add

            Returns: 
                (int) 1 if successful, 0 otherwise
            """

            for i in range(n_to_add):
                self.elements.append(  PySolTrace.Stage.Element(self, len(self.elements) ) )

            # Call the appropriate method depending on number of elements to add
            if n_to_add > 1:
                self.pdll.st_add_elements.restype = c_int 
                return self.pdll.st_add_elements(c_void_p(self.p_data), c_uint32(self.id), c_int(n_to_add))
            else:
                self.pdll.st_add_element.restype = c_int 
                return self.pdll.st_add_element(c_void_p(self.p_data), c_uint32(self.id))

    # ---------- methods of the PySolTrace class --------------------------------------------
    def __init__(self):
        cwd = os.getcwd()
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            self.pdll = CDLL(cwd + "/coretrace.dll")
            # print("Loaded win32")
            #self.pdll = CDLL(cwd + "/coretraced.dll") # for debugging
        elif sys.platform == 'darwin':
            self.pdll = CDLL(cwd + "/coretrace.dylib")  # Never tested
        elif sys.platform.startswith('linux'):
            self.pdll = CDLL(cwd +"/coretrace.so")  # Never tested
        else:
            print( 'Platform not supported ', sys.platform)

        # memory location of soltrace data context
        self.p_data = None

        # Initialize lists for optics and stage instances
        self.optics = []
        self.stages = []

        # Simulation settings
        self.num_ray_hits = int(1e5)
        self.max_rays_traced = self.num_ray_hits*100
        self.is_sunshape = True 
        self.is_surface_errors = True

        # Create an instance of soltrace in memory
        self.pdll.st_create_context.restype = c_void_p
        self.p_data = self.pdll.st_create_context()


    # def data_free(self, p_data : int) -> bool:
    def clear_context(self,*args):
        """Frees SolTrace instance from memory

        Parameters
        ----------
        p_data : int
            memory address of SolTrace instance 

        Returns
        -------
        bool
            True if successful, False otherwise
        """

        self.pdll.st_free_context.restype = c_bool
        self.pdll.st_free_context(c_void_p(self.p_data))

        return True

    def get_num_optics(self) -> int:
        """
        Get number of optical elements in the SolTrace context
        
        Returns:
            (int) number of elements
        """
        self.pdll.st_num_optics.restype = c_int
        return self.pdll.st_num_optics(c_void_p(self.p_data))


    def add_optic(self, optic_name : str):
        """
        Instantiates a new PySolTrace.Optics object, adding it to the optics list, and
        creating the optics in the SolTrace context. This method does not set optics
        properties, which instead is done using the Optics.Create method.

        Arguments:
        optic_name      | (string) Unique name for this Optics instance.

        Returns:
            (Optics) Reference to the Optics object that was just created.
        """

        new_opt_id = len(self.optics)

        self.optics.append( PySolTrace.Optics(self.pdll, self.p_data, new_opt_id ) )
        self.optics[-1].name = optic_name

        self.pdll.st_add_optic.restype = c_int
        self.pdll.st_add_optic(c_void_p(self.p_data), c_char_p(optic_name.encode()))

        return self.optics[-1]  #return the last object in the list, which was the one just created

    def delete_optic(self, optic_id : int) -> int:
        """
        Delete Optics instance. The optics object is removed from the PySolTrace.optics list and 
        from the SolTrace context.

        Arguments:
            optic_id        | (int) ID associated with the optics to be deleted

        Returns:
            (int) 1 if successful, 0 otherwise
        """

        # find the appropriate optic
        for opt in self.optics: 
            if optic_id == opt.id:
                # clear it from the optics array 
                self.optics.remove(opt)
                # Remove from the soltrace context
                self.pdll.st_delete_optic.restype = c_int 
                return self.pdll.st_delete_optic(c_void_p(self.p_data), c_uint32(optic_id) )
        
        # If reaching this point, the optic id was not found
        return 0

    def add_sun(self):
        """
        Instantiates a PySolTrace.Sun object and associates it with the PySolTrace.sun member.
        This does not create or modify the Sun data in the SolTrace context.

        Returns:
            (PySolTrace.Sun) Reference to newly created Sun instance.
        """
        self.sun = PySolTrace.Sun(self.pdll, self.p_data)
        return self.sun


    def get_num_stages(self) -> int:
        """
        Retrieve the number of stages that have been created in the SolTrace context. 

        Returns:
            (int) Number of stages.
        """
        self.pdll.st_num_stages.restype = c_int
        return self.pdll.st_num_stages(c_void_p(self.p_data))

    def add_stage(self):
        """
        Adds a new Stage instance to the PySolTrace.stages list and creates the stage in
        the SolTrace context. The Stage ID is automatically generated based on the number 
        of current stages.

        Returns:
            (PySolTrace.Stage) Reference to the newly created Stage object. 
        """
        new_st_id = len(self.stages)

        self.stages.append( PySolTrace.Stage(self.pdll, self.p_data, new_st_id ) )

        self.pdll.st_add_stage.restype = c_int
        self.pdll.st_add_stage(c_void_p(self.p_data) )

        return self.stages[-1]

    def delete_stage(self, stage_id : int) -> int:
        """
        Delete Stage instance. The stage object is removed from the PySolTrace.stages list and 
        from the SolTrace context.

        Arguments:
            stage_id       | (int) ID associated with the stage to be deleted

        Returns:
            (int) 1 if successful, 0 otherwise 
        """
        # find the appropriate optic
        for st in self.stages: 
            if stage_id == st.id:
                # clear it from the optics array 
                self.stages.remove(st)
                # Remove from the soltrace context
                self.pdll.st_delete_stage.restype = c_int 
                return self.pdll.st_delete_stage(c_void_p(self.p_data), c_uint32(stage_id) )
        
        # If reaching this point, the stage id was not found
        return 0

    def run(self, seed : int = -1, as_power_tower = False):
        """
        Run SolTrace simulation.

        Arguments:
            seed            | (int) Seed for random number generator. [-1] for random seed.
            as_power_tower  | (bool) Flag indicating simulation should be processed as power 
                            | tower / central receiver type, with corresponding efficiency adjustments.
        Returns:
            (int) Simulation return value
        """

        self.pdll.st_sim_errors.restype = c_int 
        self.pdll.st_sim_errors(c_void_p(self.p_data), c_int(1 if self.is_sunshape else 0), c_int(1 if self.is_surface_errors else 0))

        self.pdll.st_sim_params.restype = c_int 
        self.pdll.st_sim_params(c_void_p(self.p_data), c_int(int(self.num_ray_hits)), c_int(int(self.max_rays_traced)))

        self.pdll.st_sim_run.restype = c_int 
        return self.pdll.st_sim_run( c_void_p(self.p_data), c_uint16(seed), c_bool(as_power_tower), api_callback)


    def get_num_intersections(self) -> int:
        """
        [Post simulation] Get the number of ray intersections detected in the simulation.

        Returns:
            (int) Number of intersections
        """
        self.pdll.st_num_intersections.restype = c_int 
        return self.pdll.st_num_intersections(c_void_p(self.p_data))

    def get_intersect_locations(self):
        """
        [Post simulation] Get the ray intersection locations in global coordinates.

        Returns:
            ([[x,y,z],]) 2D array of ray positions
        """
        n_int = self.get_num_intersections()

        # Initialize result arrays
        loc_x = (c_number*n_int)()
        loc_y = (c_number*n_int)()
        loc_z = (c_number*n_int)()

        # Retrieve raw data
        self.pdll.st_locations.restype = c_int 
        self.pdll.st_locations(c_void_p(self.p_data), pointer(loc_x), pointer(loc_y), pointer(loc_z))
        
        # Combine data and transpose
        return list(map(list,zip(loc_x, loc_y, loc_z)))

    def get_intersect_cosines(self):
        """
        [Post simulation] Get the ray intersection cosines (direction vectors) in global coordinates.

        Returns:
            ([[cx,cy,cz],]) 2D array of ray vectors associated with the intersection at the same index.
        """
        n_int = self.get_num_intersections()

        # Initialize data vectors
        cos_x = (c_number*n_int)()
        cos_y = (c_number*n_int)()
        cos_z = (c_number*n_int)()

        self.pdll.st_cosines.restype = c_int 
        self.pdll.st_cosines(c_void_p(self.p_data), pointer(cos_x), pointer(cos_y), pointer(cos_z))

        # Combine data and transpose
        return list(map(list,zip(cos_x, cos_y, cos_z)))

    def get_intersect_elementmap(self):
        """
        [Post simulation] Get the elements associated with the ray intersection. The element numbers are
        interpreted as follows:
            0     | Ray did not hit an element in this stage 
            >= +1 | Ray hit this element ID and was reflected / refracted 
            <= -1 | Ray hit this element ID and was absorbed

        Returns:
            ([int,]) 2D array of integers associated with the intersection at the same index.
        """
        n_int = self.get_num_intersections()
        element_map = (c_int*n_int)()

        self.pdll.st_elementmap.restype = c_int 
        self.pdll.st_elementmap(c_void_p(self.p_data), pointer(element_map))

        return list(element_map)

    def get_intersect_stagemap(self):
        """
        [Post simulation] Get the stage associated with the ray intersection at the same positional index.

        Returns:
            ([int,]) List of stages associated with the intersection at the same index.
        """
        n_int = self.get_num_intersections()
        stage_map = (c_int*n_int)()

        self.pdll.st_stagemap.restype = c_int 
        self.pdll.st_stagemap(c_void_p(self.p_data), pointer(stage_map))

        return list(stage_map)

    def get_intersect_raynumbers(self):
        """
        [Post simulation] Get the ray intersection numbers. 

        Returns:
            ([int,]) List of ray intersection numbers.
        """
        n_int = self.get_num_intersections()
        raynumbers = (c_int*n_int)()

        self.pdll.st_raynumbers.restype = c_int 
        self.pdll.st_raynumbers(c_void_p(self.p_data), pointer(raynumbers))

        return list(raynumbers)

    def get_sun_stats(self):
        """
        Get information on the sun box. 

        Returns:
            (dict)      | Keys in the return dictionary are:
                        |   'xmin' --> Minimum x extent of the bounding box for hit testing 
                        |   'xmax' --> Maximum x extent of the bounding box for hit testing 
                        |   'ymin' --> Minimum y extent of the bounding box for hit testing 
                        |   'ymax' --> Maximum y extent of the bounding box for hit testing 
                        |   'nsunrays' --> Number of sun rays simulated
        """
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

    def get_ray_dataframe(self):
        """
        Get a pandas dataframe with all of the ray data from the simulation. 

        Returns:
            (Pandas.DataFrame) with columns:
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

        data['element'] = (c_int*n_int)()

        self.pdll.st_elementmap.restype = c_int 
        self.pdll.st_elementmap(c_void_p(self.p_data), pointer(data['element']))

        data['stage'] = (c_int*n_int)()

        self.pdll.st_stagemap.restype = c_int 
        self.pdll.st_stagemap(c_void_p(self.p_data), pointer(data['stage']))

        data['number'] = (c_int*n_int)()

        self.pdll.st_raynumbers.restype = c_int 
        self.pdll.st_raynumbers(c_void_p(self.p_data), pointer(data['number']))

        for key in data.keys():
            data[key] = list(data[key])

        df = pd.DataFrame(data)

        return df

    def validate(self) -> bool:
        """
        Validate that the current SolTrace context has been configured correctly, and
        that commonly missed inputs are specified.

        Returns:
            bool | False if error is identified, True otherwise
        """



        return True


    # /* utility transform/math functions */
    def util_calc_euler_angles(self, origin, aimpoint, zrot):
        """
        
        """

        a_origin = (c_number*3)()
        a_aimpoint = (c_number*3)()
        a_euler = (c_number*3)()
        a_origin[:] = origin
        a_aimpoint[:] = aimpoint

        self.pdll.st_calc_euler_angles.restype = c_void_p
        self.pdll.st_calc_euler_angles(pointer(a_origin), pointer(a_aimpoint), c_number(zrot), pointer(a_euler))

        return list(a_euler)

    def util_transform_to_local(self, posref, cosref, origin, rreftoloc):
        """
        Perform coordinate transformation from reference system to local system.
        
        Arguments:
            PosRef      | ([float,]*3) X,Y,Z coordinates of ray point in reference system
            CosRef      | ([float,]*3) Direction cosines of ray in reference system
            Origin      | ([float,]*3) X,Y,Z coordinates of origin of local system as measured
                        |   in reference system
            RRefToLoc   | Rotation matrices required for coordinate transform from reference to local
        
        Returns:
            (dict)  Keys in return dictionary include:
                posloc : ([float,]*3) X,Y,Z coordinates of ray point in local system
                cosloc : ([float,]*3) Direction cosines of ray in local system
        """
        # Allocate space for input 
        a_posref = (c_number*3)()
        a_cosref = (c_number*3)()
        a_origin = (c_number*3)() 
        a_rreftoloc = (c_number*9)() 
        #output
        posloc = (c_number*3)()
        cosloc = (c_number*3)()
        # Assign input
        a_posref[:] = posref
        a_cosref[:] = cosref
        a_origin[:] = origin
        a_rreftoloc[:] = sum([], rreftoloc)

        self.pdll.st_transform_to_local.restype = c_void_p
        self.pdll.st_transform_to_local(pointer(a_posref), pointer(a_cosref), pointer(a_origin), pointer(a_rreftoloc), pointer(posloc), pointer(cosloc))

        return {'cosloc':list(cosloc), 'posloc':list(posloc)}

    def util_transform_to_ref(self, posloc, cosloc, origin, rloctoref):
        """
        Perform coordinate transformation from local system to reference system.
        
        Arguments:
            PosLoc      | ([float,]*3) X,Y,Z coordinates of ray point in local system
            CosLoc      | ([float,]*3) Direction cosines of ray in local system
            Origin      | ([float,]*3) X,Y,Z coordinates of origin of local system as measured
                        |   in reference system
            RLocToRef   | Rotation matrices required for coordinate transform from local to reference
                        |   -- inverse of reference to local transformation
        
        Returns:
            (dict)  Keys in return dictionary include:
                posref : ([float,]*3) X,Y,Z coordinates of ray point in reference system
                cosref : ([float,]*3) Direction cosines of ray in reference system
        """

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

    def util_matrix_vector_mult(self, m, v):
        """
        Perform multiplication of a 3x3 matrix and a length-3 vector, returning the result vector.

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

    def util_calc_transforms(self, euler):
        """
        Calculate matrix transforms

        Arguments:
        euler       | ([float,]*3) Euler angles

        Returns:
            (dict) A dictionary containing the keys:
                rreftoloc : Transformation matrix from Reference to Local system
                rloctoref : Transformation matrix from Local to Reference system
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

    def util_matrix_transpose(self, m):
        """
        Calculate matrix transpose 

        Arguments:
            m       | ([[float,]*3]*3) Square matrix 3x3 to be transposed.

        Returns:
            ([[float,]*3]*3)    | Square matrix 3x3, transpose of 'm'
        """

        a_m = (c_number*9)()
        a_m[:] = sum([], m)
        output = (c_number*9)()

        self.pdll.st_matrix_transpose.restype = c_void_p
        self.pdll.st_matrix_transpose(pointer(a_m), pointer(output))

        a_output = []
        for i in range(0,10,3):
            a_output.append(output[i:i+3])

        return a_output

        