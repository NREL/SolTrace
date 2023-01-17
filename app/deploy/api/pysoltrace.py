from datetime import datetime
import sys, os, copy
from typing import List
import pandas as pd
from ctypes import *
c_number = c_double   #must be either c_double or c_float depending on coretrace definition
import multiprocessing
import time
import math
import random

# Callback to print command line progress messages
@CFUNCTYPE(c_int, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32, c_void_p)
def api_callback(ntracedtotal, ntraced, ntotrace, curstage, nstages, thread_id):
    if thread_id:
        if thread_id > 1:
            return 1

    w = 50
    prog = 100.*float(ntraced)/float(ntotrace)
    pbprog = int(prog * w / 100.)
    pbar = pbprog *"▮" + (w-pbprog)*"▯" 

    print("{:s}  | Stage ({:d}/{:d}) - Complete {:.2f}%".format(pbar, curstage, nstages, prog), end='\r')
    return 1


def _thread_func(pobj, seed, id):
    pobj.run(seed,True, 1, id)
    return copy.deepcopy(pobj.raydata), copy.copy(pobj.sunstats)


# ==========================================================================================
class Point:
    """
    A simple class to manage points in Cartesian coordinates.
    """
    def __init__(self, x=0, y=0, z=0):
        """
        Parameters
        ==========
        x : float
            x-coordinate
        y : float
            y-coordinate
        z : float
            z-coordinate
        """
        ## (float) x-coordinate
        self.x = x  
        ## (float) y-coordinate
        self.y = y  
        ## (float) z-coordinate
        self.z = z  
        return 
    def copy(self):
        pnew = Point()
        c = self.__dict__.copy()
        for attr in self.__dict__.keys():
            pnew.__setattr__(attr, c[attr])
        return pnew
    def __str__(self):
        return "[{:f}, {:f}, {:f}]".format(self.x, self.y, self.z)
    def __add__(self, obj):
        """
        Add to the current point coordinate values. 

        Parameters
        ==========
        obj : variant
            If obj = (Point), adds component-wise to the current point 
            If obj = (float), adds obj to each component
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, Point):
            return Point(self.x + obj.x, self.y + obj.y, self.z + obj.z)
        elif isinstance(obj, float):
            return Point(self.x + obj, self.y + obj, self.z + obj)
        else:
            raise ValueError("Invalid addition operator object")
        return self
    def __sub__(self, obj):
        """
        Subtract from the current point coordinate values. 

        Parameters
        ==========
        obj : variant
            If obj = (Point), subtracts component-wise from the current point 
            If obj = (float), subtracts obj from each component
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, Point):
            return Point(self.x - obj.x, self.y - obj.y, self.z - obj.z)
        elif isinstance(obj, float):
            return Point(self.x - obj, self.y - obj, self.z - obj)
        else:
            raise ValueError("Invalid subtraction operator object")
        return self
    def __mul__(self, obj):
        """
        Multiplies the current point coordinate values. 

        Parameters
        ==========
        obj : variant
            If obj = (Point), multiplies the current point component-wise by obj 
            If obj = (float), multiplies each component of the current point by obj
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, (float,int)):
            return Point(self.x*obj, self.y*obj, self.z*obj)
        else:
            raise ValueError("Invalid multiplication operator object")
    def __floordiv__(self, obj):
        """
        Divides the current point coordinate values, taking the floor of the result. 

        Parameters
        ==========
        obj : variant
            If obj = (Point), divides current point component-wise 
            If obj = (float), divides components of current point by obj
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, (float,int)):
            return Point(self.x//obj, self.y//obj, self.z//obj)
        else:
            raise ValueError("Invalid division operator object")
    def __truediv__(self, obj):
        """
        Divides the current point coordinate values/ 

        Parameters
        ==========
        obj : variant
            If obj = (Point), divides current point component-wise 
            If obj = (float), divides components of current point by obj
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, (float,int)):
            return Point(self.x/obj, self.y/obj, self.z/obj)
        else:
            raise ValueError("Invalid division operator object")
    def radius(self):
        """
        Computes distance between current point and the origin (0,0,0)

        Returns
        =======
        float
            Calculated radius
        """
        mag = (self.x*self.x + self.y*self.y + self.z*self.z)**0.5
        return mag
    def unitize(self, inplace : bool = False):
        """
        Converts current point into a unit vector with magnitude 1

        Parameters
        ==========
        inplace : bool = False
            Specifies whether the current point is converted to a unit vector in place, or whether the 
            current point remains unchanged and a unitized copy of the vector is returned.
        """
        mag = self.radius()
        if mag > 0:
            if inplace:
                self.x /= mag
                self.y /= mag 
                self.z /= mag
            else:
                return Point(self.x/mag, self.y/mag, self.z/mag)
# ----------------------------------------------------------------------
class PySolTrace:
    """
    A class to access PySolTrace (SolTrace's Python API)

    Attributes
    ----------
    pdll : ctypes.CDLL
        loaded SolTrace library of exported functions
    p_data 
        Memory location for soltrace context
    optics : [PySolTrace.Optics,]
        List of Optics instances
    stages : [PySolTrace.Stage,]
        List of Stage instances
    num_ray_hits : int
        Minimum number of simulation ray hits 
    max_rays_traced : int
        Maximum number of ray hits in a simulation
    is_sunshape : bool
        Flag indicating whether sunshape should be included
    is_surface_errors : bool
        Flag indicating whether surface errors should be included

    Methods
    -------
    clear_context 
        Frees SolTrace instance from memory at p_data
    add_optics 
        Instantiates a new PySolTrace.Optics object
    delete_optic 
        Delete Optics instance
    add_sun 
        Adds Sun instance
    add_stage 
        Adds Stage instance
    validate 
        Detect common errors in simulation setup 
    run 
        Runs SolTrace simulation
    get_num_intersections 
        Get the number of simulation ray intersections
    get_ray_dataframe 
        Get all simulation ray data in pandas dataframe
    get_intersect_locations 
        Get ray intersection locations
    get_intersect_cosines 
        Get ray intersection cosines / direction vectors
    get_intersect_elementmap 
        Get ray intersection associated element
    get_intersect_stagemap 
        Get ray intersection associated stage
    get_intersect_raynumbers 
        Get ray intersection number
    get_sun_stats 
        Get information on sun box
    util_calc_euler_angles 
        Calculate Euler angles for a position, aimpoint, and rotation
    util_transform_to_local 
        Transform a coordinate system from reference to a local system
    util_transform_to_ref 
        Transform a coordinate system from local to reference system
    util_matrix_vector_mult 
        Calculate product of a square matrix and a vector
    util_calc_transforms 
        Calculate matrix transforms
    util_matrix_transpose 
        Compute the transpose of a matrix
    """

    class Optics:
        """ *Optics* is a subclass of PySolTrace, and represents an optical property set. 
        A PySolTrace instance may have multiple Optics member instances, which are stored in 
        the PySolTrace.optics list.
        
        Optics contains a subclass *Face*, which collects properties associated with the front 
        or back face of an optical surface.

        Attributes
        ----------
        p_dll 
            Reference for API DLL, managed by PySolTrace
        p_data 
            Memory location for soltrace context, managed by PySolTrace
        name : str
            Unique name for the optical property set
        id : int
            Identifying integer associated with the property set
        front : Face
            properties associated with the front of the optical surface
        back : Face
            properties associated with the back of the optical surface

        Methods
        ----------
        Create 
            Calls methods to instantiate and construct optical surface in the SolTrace context.
        """
        class Face:
            """
            Subclass of Optics, contains properties associated with one of the optical faces.

            Attributes
            ----------
            dist_type : char
                Distribution type for surface interactions. One of:
                {'g':Gaussian, 'p':Pillbox, 'd':Diffuse }
            refraction_real : float
                Real component of the refraction index
            reflectivity : float
                [0..1] Surface reflectivity
            transmissivity : float
                [0..1] Surface transmissivity
            slope_error : float
                [mrad] Surface RMS slope error, half-angle
            spec_error : float
                [mrad] Surface specularity error, half-angle
            userefltable : bool
                Flag specifying use of user reflectivity table to modify reflectivity as a function of incidence angle
            refltable : [[float,float],]
                [mrad,0..1] 2D list containing pairs of [angle,reflectivity] values.              
            usetranstable : bool
                Flag specifying use of user transmissivity table to modify transmissivity as a function of incidence angle
            transtable : [[float,float],]
                [mrad,0..1] 2D list containing pairs of [angle,transmissivity] values.    
            """
            def __init__(self):
                ## Distribution type for surface interactions. One of:
                # {'g':Gaussian, 'p':Pillbox, 'd':Diffuse }
                self.dist_type = 'g'     #One of 'g'->Gauss 'p'->Pillbox 'd'->Diffuse
                ## Real component of the refraction index
                self.refraction_real = 1.1         #real component of the refraction index 
                ## [0..1] Surface reflectivity
                self.reflectivity = 0.96         #reflectivity
                ## [0..1] Surface transmissivity
                self.transmissivity = 0.         #transmissivity
                ## [mrad] Surface RMS slope error, half-angle
                self.slope_error = 0.95          #RMS slope error [mrad] 
                ## [mrad] Surface specularity error, half-angle
                self.spec_error = 0.2            #RMS specularity error [mrad]
                ## Flag specifying use of user reflectivity table to modify reflectivity as a function of incidence angle
                self.userefltable = False             #Flag [bool] use reflectivity table 
                ## [mrad,0..1] 2D list containing pairs of [angle,reflectivity] values.      
                self.refltable = []  #[[angle1,refl1],[...]] 
                ## Flag specifying use of user transmissivity table to modify transmissivity as a function of incidence angle
                self.usetranstable = False             #Flag [bool] use transmissivity table 
                ## [mrad,0..1] 2D list containing pairs of [angle,transmissivity] values.      
                self.transtable = []  #[[angle1,trans1],[...]] 

            def copy(self, fnew):
                c = self.__dict__.copy()
                for attr in self.__dict__.keys():
                    fnew.__setattr__(attr, copy.deepcopy(c[attr]))
                return

        # -------- methods of the Optics class -----------------------------------------
        def __init__(self, id : int):
            
            ## Unique name for the optical property set
            self.name = "new optic"
            ## Identifying integer associated with the property set
            self.id = id 

            ## properties associated with the front of the optical surface
            self.front = PySolTrace.Optics.Face()
            ## properties associated with the back of the optical surface
            self.back = PySolTrace.Optics.Face()

        def copy(self, onew):
            onew.name = copy.copy(self.name)
            onew.id = copy.copy(self.id)

            self.front.copy(onew.front)
            self.back.copy(onew.back)

            return 

        def Create(self, pdll, p_data) -> int:
            """
            Create Optics instance in the SolTrace context.

            Returns
            ----------
            int 
                1 if successful, 0 otherwise
            """
            pdll.st_add_optic.restype = c_int
            pdll.st_add_optic(c_void_p(p_data), c_char_p(self.name.encode()))

            pdll.st_optic.restype = c_int

            dummy_grating = (c_number*3)()

            resok = True
            
            # for each face -- front or back
            for i,opt in enumerate([self.front, self.back]):

                user_refl_angles = (c_number*len(opt.refltable))()
                user_refls = (c_number*len(opt.refltable))()
                if len(opt.refltable) > 1:
                    user_refl_angles[:] = list(list(zip(*opt.refltable))[0])
                    user_refls[:] = list(list(zip(*opt.refltable))[1])

                user_trans_angles = (c_number*len(opt.transtable))()
                user_trans = (c_number*len(opt.transtable))()
                if len(opt.transtable) > 1:
                    user_trans_angles[:] = list(list(zip(*opt.transtable))[0])
                    user_trans[:] = list(list(zip(*opt.transtable))[1])

                # Create surface optic
                resok = resok and pdll.st_optic( \
                    c_void_p(p_data),
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
                    pointer(user_refl_angles),
                    pointer(user_refls),
                    c_int(1 if opt.usetranstable else 0),
                    c_int(len(opt.transtable)),
                    pointer(user_trans_angles),
                    pointer(user_trans),
                    )

            return 1 if resok else 0
    # ========end Optics class =================================================================
            

    # ==========================================================================================
    class Sun:
        """ *Sun* is a subclass of PySolTrace, and represents a sun property set. 
        A PySolTrace instance may have a single Sun member instance, which is stored as the 
        PySolTrace.sun member.

        Attributes
        ----------
        p_dll 
            Reference for API DLL, managed by PySolTrace
        p_data 
            Memory location for soltrace context, managed by PySolTrace
        point_source : bool
            Flag indicating whether the sun is modeled as a point source at a finite distance.
        shape : char
            Sun shape model. One of: {'p':Pillbox, 'g':Gaussian, 'd':data table, 'f':gray diffuse}
        sigma : float
            [mrad] Half-width or std. dev. of the error distribution
        position : Point
            Location of the sun/sun vector in global coordinates
        user_intensity_table : [[float,float],]
            [mrad, 0..1] 2D list containing pairs of
            angle deviation from sun vector and irradiation intensity.
            A typical table will have angles spanning 0->~5mrad, and inten-
            sities starting at 1 and decreasing to zero. The table must
            contain at least 2 entries.    
        Methods
        ----------
        Create 
            Calls methods to instantiate and construct optical surface in the SolTrace context.

        """
        def __init__(self):

            ## Flag indicating whether the sun is modeled as a point source at a finite distance.
            self.point_source = False 
            ## Sun shape model. One of: {'p':Pillbox, 'g':Gaussian, 'd':data table, 'f':gray diffuse}
            self.shape = 'p'
            ## [mrad] Half-width or std. dev. of the error distribution
            self.sigma = 4.65
            ## Location of the sun/sun vector in global coordinates
            self.position = Point()
            self.position.z = 100.

            ## [mrad, 0..1] 2D list containing pairs of
            # angle deviation from sun vector and irradiation intensity.
            # A typical table will have angles spanning 0->~5mrad, and inten-
            # sities starting at 1 and decreasing to zero. The table must
            # contain at least 2 entries.    
            self.user_intensity_table = []   

        def copy(self, snew):
            snew.position = self.position.copy() 
            c = self.__dict__.copy()
            for attr in self.__dict__.keys():
                if attr in ['_pdll','_p_data','position']:
                    continue
                snew.__setattr__(attr, copy.deepcopy(c[attr]))
            return 

        def Create(self, pdll, p_data):
            """
            Create Sun instance in the SolTrace context.

            Returns
            ----------
            int 
                1 if successful, 0 otherwise
            """

            pdll.st_sun.restype = c_int 
            pdll.st_sun_xyz.restype = c_int 

            pdll.st_sun(c_void_p(p_data), c_int(int(self.point_source)), c_wchar(self.shape[0]), c_number(self.sigma))
            pdll.st_sun_xyz(c_void_p(p_data), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))

            # If a user intensity table is provided, and the shape is specified accordingly as 'd', load the data table into context
            if len(self.user_intensity_table) > 2 and self.shape.lower()[0] == 'd':
                user_angles = (c_number*len(self.user_intensity_table))()
                user_ints = (c_number*len(self.user_intensity_table))()
                user_angles[:] = list(list(zip(*self.user_intensity_table))[0])
                user_ints[:] = list(list(zip(*self.user_intensity_table))[1])

                pdll.st_sun_userdata.restype = c_int 
                return pdll.st_sun_userdata(c_void_p(p_data), c_uint32(len(self.user_intensity_table)), pointer(user_angles), pointer(user_ints))

            return 1
    # ===========end of the Sun class===========================================================

    # ==========================================================================================
    class Stage:
        """ *Stage* is a subclass of PySolTrace, and represents a grouping of elements. 
        A PySolTrace instance may have multiple Stage member instances, which are stored in 
        the PySolTrace.stages list.

        Stage contains a subclass *Element*, which collects properties and geometry associated
        with individual geometric elements.

        Attributes
        ----------
        p_dll 
            Reference for API DLL, managed by PySolTrace
        p_data 
            Memory location for soltrace context, managed by PySolTrace
        id : int
            Identifying integer associated with the stage
        position : Point
            Stage location in global coordinates 
        aim : Point
            Coordinate system aim point in global coordinates
        zrot : float
            [deg] Rotation of coordinate system around z-axis
        is_virtual : bool
            Flag indicating virtual stage
        is_multihit : bool
            Flag indicating that rays can have multiple interactions within a single stage.
        is_tracethrough : bool
            Flag indicating the stage is in trace-through mode
        name : str 
            Descriptive name for this stage
        elements : [Stage.Element,]
            list of all elements in the stage
            
        Methods
        ----------
        Create 
            Calls methods to instantiate and construct a stage in the context.
        add_elements 
            Creates new element in Stage.element[] list
        """

        class Element:
            """ *Element* is a subclass of PySolTrace.Stage, and represents a set of properties and
            geometric settings related to a single optical element in SolTrace. Elements are associ-
            ated with a single stage, and are stored in the respective stage's Stage.elements[] list.

            Attributes
            ----------
            p_dll 
                Reference for API DLL, managed by PySolTrace
            p_data 
                Memory location for soltrace context, managed by PySolTrace
            stage_id : int
                Identifying integer associated with the containing stage
            id : int
                Identifying integer associated with element
            enabled : bool
                Flag indicating whether the element is included in the model
            position : Point
                Element location in stage coordinates 
            aim : Point
                Element coordinate system aim point in stage coordinates
            zrot : float
                [deg] Rotation of coordinate system around z-axis
            interaction : int
                Flag indicating optical interaction type. {1:refraction, 2:reflection}
            optic : Optics
                Reference to *Optics* instance associated with this element
            aperture : char
                Charater indicating aperture type. One of:
                {'c':circle, 'h':hexagon, 't':triangle, 'r':rectangle, 'a':annulus,
                'l':single-axis curvature, 'i':irregular triangle, 'q':quadrilateral}
            aperture_params : [float,]
                Up to 8 coefficients defining aperture -- values depend on selection for 'aperture'
            surface : char
                Character indicating surface type. One of:
                {'s':spherical, 'p':parabolic, 'f':flat plane, 'o':hyperboloid/ellipsoid,
                'c':conical, 't':cylindrical, 'd':toroid, 'm':Zernicke monomial, 
                'r':Polynomial revolution, 'i':cubic spline interpolation, 
                'e':finite element data, 'v':VSHOT data}
            surface_params : [float,]
                Up to 8 coefficients defining surface -- values depend on selection for 'surface'
            surface_file : string
                Name for surface file, if using compatible type. File extension:
                *.mon --> 'm' / Zernicke monomial 
                *.sht --> 'v' / VSHOT data
                *.ply --> 'r' / Polynomial revolution
                *.csi --> 'i' / Cubic spline interpolation
                *.fed --> 'e' / Finite element data
            
            Methods
            ----------
            Create 
                Calls methods to instantiate and construct element in the SolTrace context
            surface_xxxx 
                Family of methods that compute surface coefficients

            """
            
            # STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file);
            def __init__(self, parent_stage, element_id : int):
                ## Identifying integer associated with the containing stage
                self.stage_id = parent_stage.id
                ## Identifying integer associated with element
                self.id = element_id
                ## Flag indicating whether the element is included in the model
                self.enabled = True 
                ## Element location in stage coordinates 
                self.position = Point()
                ## Element coordinate system aim point in stage coordinates
                self.aim = Point()
                self.aim.z = 1
                ## [deg] Rotation of coordinate system around z-axis
                self.zrot = 0.
                ## Charater indicating aperture type. 
                self.aperture = 'r'
                ## Up to 8 coefficients defining aperture -- values depend on selection for 'aperture'
                self.aperture_params = [0. for i in range(8)]
                ## Character indicating surface type. 
                self.surface = 'f'
                ## Up to 8 coefficients defining surface -- values depend on selection for 'surface'
                self.surface_params = [0. for i in range(8)]
                ## Name for surface file, if using compatible type. 
                self.surface_file = None 
                ## Flag indicating optical interaction type. {1:refraction, 2:reflection}
                self.interaction = 2        #1=refract, 2=reflect
                ## Reference to *Optics* instance associated with this element
                self.optic = None

            def copy(self, enew):
                c = self.__dict__.copy()
                for attr in self.__dict__.keys():
                    if attr in ['_pdll','_p_data','optic']:
                        continue
                    else:
                        enew.__setattr__(attr, copy.deepcopy(c[attr]))
                return

            def Create(self, pdll, p_data) -> int:
                """
                Create Element instance in the SolTrace context.

                Returns
                ----------
                int 
                    1 if successful, 0 otherwise
                """

                pdll.st_add_element.restype = c_int 
                pdll.st_add_element(c_void_p(p_data), c_uint32(self.stage_id))

                pdll.st_element_enabled.restype = c_int
                pdll.st_element_xyz.restype = c_int
                pdll.st_element_aim.restype = c_int
                pdll.st_element_zrot.restype = c_int
                pdll.st_element_aperture.restype = c_int
                pdll.st_element_aperture_params.restype = c_int
                pdll.st_element_surface.restype = c_int
                pdll.st_element_surface_params.restype = c_int
                pdll.st_element_surface_file.restype = c_int
                pdll.st_element_interaction.restype = c_int
                pdll.st_element_optic.restype = c_int

                aperture_params = (c_number*len(self.aperture_params))()
                surface_params = (c_number*len(self.surface_params))()
                aperture_params[:] = self.aperture_params
                surface_params[:] = self.surface_params

                pdll.st_element_enabled(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(int(self.enabled)));
                pdll.st_element_xyz(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.position.x),  c_number(self.position.y), c_number(self.position.z));
                pdll.st_element_aim(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.aim.x),  c_number(self.aim.y), c_number(self.aim.z));
                pdll.st_element_zrot(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_number(self.zrot) );
                pdll.st_element_aperture(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.aperture[0]));
                pdll.st_element_aperture_params(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(aperture_params));
                pdll.st_element_surface(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_wchar(self.surface[0]));
                pdll.st_element_surface_params(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), pointer(surface_params));
                if self.surface_file:
                    pdll.st_element_surface_file(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.surface_file.encode()));
                pdll.st_element_interaction(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_int(self.interaction)); #/* 1=refract, 2=reflect */
                pdll.st_element_optic(c_void_p(p_data), c_uint32( self.stage_id ), c_uint32( self.id ), c_char_p(self.optic.name.encode()));

                return 1

            def surface_spherical(self, radius):
                """
                Set up the surface as spherical type

                Parameters
                ==========
                radius : float
                    Radius of the spherical surface
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = 1. / radius
                self.surface = 's'
                return True
            def surface_parabolic(self, focal_len_x, focal_len_y):
                """
                Set up the surface as parabolic

                Parameters
                ==========
                focal_len_x : float
                    Focal length of the surface in the x-direction. If infinite, use float('inf')
                focal_len_y : float
                    Focal length of the surface in the y-direction. If infinite, use float('inf')
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = 1. / (2.*focal_len_x)
                self.surface_params[1] = 1. / (2.*focal_len_y)
                self.surface = 'p'
                return True
            def surface_flat(self):
                """
                Set up the surface as flat
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'f'
                return True
            def surface_hypellip(self, vertex_curv, kappa):
                """
                Set up the surface described by equation:
                    Z(x,y) = ( vertex_curv*(x^2 + y^2) ) / 
                                (1 + sqrt(1-kappa*vertex_curv^2*(x^2 + y^2)))
                Parameters
                ----------
                vertex_curv 
                    Curvature parameter
                kappa 
                    Form parameter. Value of parameter determines geometry as follows:
                    kappa < 0 --> tall hyperboloid
                    kappa 0..1 --> ellipsoid 
                    kappa > 1 --> stout ellipsoid
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = vertex_curv
                self.surface_params[1] = kappa
                self.surface = 'o'
                return True
            def surface_conical(self, theta):
                """
                Set up the surface described by cone with half-angle theta. 

                Parameters
                ----------
                theta : deg
                    half-angle of cone
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = theta
                self.surface = 'c'
                return True
            def surface_cylindrical(self, radius):
                """
                Set up the surface as cylindrical.

                Parameters
                ----------
                radius 
                    Radius of the cylinder
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = 1./radius
                self.surface = 't'
                return True
            def surface_toroid(self, rad_annulus, rad_ring):
                """
                Set up the surface as a toroid "donut".

                Parameters
                ---------- 
                rad_annulus 
                    Radius of the 'tube', the distance between the min and max radii of the torus
                rad_ring 
                    The radius of the centerpoint of the annular tube
                """
                self.surface_params = [0. for i in range(8)]
                self.surface_params[0] = rad_annulus
                self.surface_params[1] = rad_ring
                self.surface = 'd'
                return True
            def surface_zernicke(self, file_path):
                """
                Set up the surface from a file as a Zernicke surface, where the surface is described by the equation:
                Z(x,y) = sum_i=0^N
                            sum_j=0^i  Bi,j * x^j * y^(i-j) 

                Accepts *mon file extension specifying the Zernicke coefficients. 
                File format should be a single data column:
                    N
                    B0,0
                    B1,0
                    B1,1
                    B2,1
                    B2,2
                    B2,3
                    ...
                    BN,N

                Parameters
                ----------
                file_path
                    Path to the file containing the data.
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'm'
                self.surface_file = file_path
                return True
            def surface_polynomialrev(self, file_path):
                """
                Set up the surface from a file as a rotationally symmetric polynomial, where the surface is described by
                the equation:
                Z(r) = sum_i=0^N  C_i * r^i,  where r=sqrt(x^2 + y^2)

                Accepts *ply file extension specifying equation coefficients.
                File format should be a single data column:
                    N
                    C0
                    C1
                    C2
                    ...
                    C,N
                
                Parameters
                ----------
                file_path
                    Path to the file containing the data.
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'r'
                self.surface_file = file_path
                return True
            def surface_cubicspline(self, file_path):
                """
                Set up the surface from a file as a rotationally symmetric cubic spline. Accepts *csi file extension. 
                File format should be two tab-separated columns:
                    N
                    r1      Z1
                    r2      Z2
                    r3      Z3
                    ...
                    rN      ZN
                    dZ/dr1  dZ/drN

                Parameters
                ----------
                file_path
                    Path to the file containing the data.
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'i'
                self.surface_file = file_path
                return True
            def surface_finiteelement(self, file_path):
                """
                Set up the surface from a file using finite element data specifying the vertices of the elements in 
                x,y,z coordinates. 

                Accepts the *.fed file extension. File format should be 3 tab-separated
                columns:
                    N
                    x1      y1      z1
                    x2      y2      z2
                    x3      y3      z3
                    ...
                    xN      yN      zN

                Parameters
                ----------
                file_path
                    Path to the file containing the data.
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'e'
                self.surface_file = file_path
                return True
            def surface_vshot(self, file_path):
                """
                Set up the surface from a file using VSHOT data specifying matrix coefficients generated by a VSHOT test.

                Accepts the *.sht file extension. File format should be:
                    First line - file name (skipped)
                    Radius      Focal length        Target-dist
                    0           order               num points
                    rmsslope    rmsscale
                    b00
                    b10
                    b11
                    b20
                    b21
                    b22
                    ...
                    bDD    || where 'D' is order
                    a1      b1      c1      d1      e1
                    a2      b2      c2      d2      e2
                    a3      b3      c3      d3      e3
                    ...
                    aN      bN      cN      dN      eN  || where 'N' is num points

                Parameters
                ----------
                file_path
                    Path to the file containing the data.
                """
                self.surface_params = [0. for i in range(8)]
                self.surface = 'v'
                self.surface_file = file_path
                return True
            # ---------------------
            def aperture_circle(self, diameter):
                """
                Set up the aperture as circular with 'diameter'. 

                Parameters
                ----------
                diameter 
                    Diameter of the circle
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = diameter
                self.aperture = 'c'
                return True
            def aperture_hexagon(self, diameter):
                """
                Set up the aperture as a hexagon centered at x=0,y=0. The hexagon is circumscribed
                by a circle of 'diameter'.

                Parameters
                ----------
                diameter 
                    Diameter of the circumscribing circle.
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = diameter
                self.aperture = 'h'
                return True
            def aperture_triangle(self, diameter):
                """
                Set up the aperture as a equilateral triangle with centroid at x=0,y=0. The triangle
                is circumscribed by a circle of 'diameter'.

                Parameters
                ----------
                diameter 
                    Diameter of the circumscribing circle.
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = diameter
                self.aperture = 't'
                return True
            def aperture_rectangle(self, length_x, length_y):
                """
                Set up the aperture as a rectangle with the centroid at x=0,y=0. 

                Parameters
                ----------
                length_x 
                    Width in x-coordinate direction
                length_y 
                    Height in y-coordinate direction
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = length_x
                self.aperture_params[1] = length_y
                self.aperture = 'r'
                return True
            def aperture_annulus(self, r_inner, r_outer, theta):
                """
                Set up the aperture as annular, where aperture is the annulus between to specified radii
                and within an angular slice 'theta' which is centered around the x-axis.

                Parameters
                ----------
                r_inner 
                    Inner radius of annular region
                r_outer 
                    Outer radius of annular region
                theta : deg
                    Slice of the circle contained, centered around x-axis
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = r_inner
                self.aperture_params[1] = r_outer
                self.aperture_params[2] = theta
                self.aperture = 'a'
                return True
            def aperture_singleax_curve(self, x1, x2, L):
                """
                Set up the aperture as revolved around a single axis. Revolved window is between two
                coordinates x1->x1, both non-negative and with x2 > x1. The aperture has 
                length 'L' in the y-direction. 

                This aperture is often used with a cylindrical surface. In this case, 
                both x1 and x2 should be zero, and the cylinder height specified with 'L'.

                Parameters
                ----------
                x1 
                    inner coordinate of revolved section
                x2 
                    outer coordinate of revolved section
                L 
                    length of revolved section along axis of revolution
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = x1
                self.aperture_params[1] = x2
                self.aperture_params[2] = L
                self.aperture = 'l'
                return True
            def aperture_irr_triangle(self, x1, y1, x2, y2, x3, y3):
                """
                Set up the aperture as a triangle given by three (x,y) coordinate pairs.

                Parameters
                ----------
                x1 
                    x-coordinate, point 1
                y1 
                    y-coordinate, point 1
                x2 
                    x-coordinate, point 2
                y2 
                    y-coordinate, point 2
                x3 
                    x-coordinate, point 3
                y3 
                    y-coordinate, point 3
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = x1
                self.aperture_params[1] = y1
                self.aperture_params[2] = x2
                self.aperture_params[3] = y2
                self.aperture_params[4] = x3
                self.aperture_params[5] = y3
                self.aperture = 'i'
                return True
            def aperture_quadrilateral(self, x1, y1, x2, y2, x3, y3, x4, y4):
                """
                Set up the aperture as a quadrilateral given by four (x,y) coordinate pairs.

                Parameters
                ----------
                x1 
                    x-coordinate, point 1
                y1 
                    y-coordinate, point 1
                x2 
                    x-coordinate, point 2
                y2 
                    y-coordinate, point 2
                x3 
                    x-coordinate, point 3
                y3 
                    y-coordinate, point 3
                x4 
                    x-coordinate, point 4
                y4 
                    y-coordinate, point 4
                """
                self.aperture_params = [0. for i in range(8)]
                self.aperture_params[0] = x1
                self.aperture_params[1] = y1
                self.aperture_params[2] = x2
                self.aperture_params[3] = y2
                self.aperture_params[4] = x3
                self.aperture_params[5] = y3
                self.aperture_params[6] = x4
                self.aperture_params[7] = y4
                self.aperture = 'q'
                return True
        # -------------------------- end Element class ---------------------------------


        # -----------methods of the 'Stage' class --------------------------------------
        def __init__(self, id : int):
            """
            """
            ## Identifying integer associated with the stage
            self.id = id
            ## Stage location in global coordinates 
            self.position = Point()
            ## Coordinate system aim point in global coordinates
            self.aim = Point()
            self.aim.z = 1
            ## [deg] Rotation of coordinate system around z-axis
            self.zrot = 0.
            ## Flag indicating virtual stage
            self.is_virtual = False
            ## Flag indicating that rays can have multiple interactions within a single stage.
            self.is_multihit = True 
            ## Flag indicating the stage is in trace-through mode
            self.is_tracethrough = False 
            ## Descriptive name for this stage
            self.name = "stage_{:d}".format(id)

            ## list of all elements in the stage
            self.elements = []
            return 

        def copy(self, snew):
            c = self.__dict__.copy()
            for attr in self.__dict__.keys():
                if attr in ['elements', '_pdll', '_p_data']:
                    continue
                elif attr in ['position','aim']:  #points
                    snew.__setattr__(attr, c[attr].copy())
                else:
                    snew.__setattr__(attr, copy.deepcopy(c[attr]))

            snew.elements = [PySolTrace.Stage.Element(snew, el.id) for el in self.elements]
            for i in range(len(self.elements)):
                self.elements[i].copy(snew.elements[i])

            return 


        def Create(self, pdll, p_data) -> int:
            """
            Create Stage instance in the SolTrace context. 
            Note: This does not create any associated Elements, which must have their Create method called separately. 

            Returns
            ----------
            int 
                1 if successful, 0 otherwise
            """

            pdll.st_add_stage.restype = c_int
            pdll.st_add_stage(c_void_p(p_data) )

            pdll.st_stage_xyz.restype = c_int 
            pdll.st_stage_aim.restype = c_int 
            pdll.st_stage_zrot.restype = c_int 
            pdll.st_stage_flags.restype = c_int 

            pdll.st_stage_xyz(c_void_p(p_data), c_uint32(self.id), c_number(self.position.x), c_number(self.position.y), c_number(self.position.z))
            pdll.st_stage_aim(c_void_p(p_data), c_uint32(self.id), c_number(self.aim.x), c_number(self.aim.y), c_number(self.aim.z))
            pdll.st_stage_zrot(c_void_p(p_data), c_uint32(self.id), c_number(self.zrot))
            pdll.st_stage_flags(c_void_p(p_data), c_uint32(self.id), c_int(int(self.is_virtual)), c_int(int(self.is_multihit)), c_int(int(self.is_tracethrough)))

            for element in self.elements:
                element.Create(pdll,p_data)

            return 1

        def add_element(self) -> int:
            """
            Add one element to the stage. This method appends an Element object to the
            stage's Stage.elements list. 
            To update element properties and settings, call the Element.Create method 
            on each element.

            Returns
            ----------
            PySolTrace.Stage.Element 
                Reference to the newly created element
            """
            
            new_e = PySolTrace.Stage.Element(self, len(self.elements) )
            self.elements.append( new_e )
            return new_e

    # ---------- methods of the PySolTrace class --------------------------------------------
    def __init__(self):

        # Initialize lists for optics and stage instances
        ## List of Optics instances
        self.optics = []
        ## List of Stage instances
        self.stages = []
        ## sun
        self.sun = None

        # Simulation settings
        ## Minimum number of simulation ray hits 
        self.num_ray_hits = int(1e5)
        ## Maximum number of ray hits in a simulation
        self.max_rays_traced = self.num_ray_hits*100
        ## Flag indicating whether sunshape should be included
        self.is_sunshape = True 
        ## Flag indicating whether surface errors should be included
        self.is_surface_errors = True
        # Placeholder for output ray data
        self.raydata = None
        # Placeholder for sunstats data
        self.sunstats = None
        # Placeholder for power per ray
        self.powerperray = None
        # Direct normal irradince for calculations
        self.dni = 1000.  #w/m^2

    def copy(self):
        psnew = PySolTrace()
        for attr in ['num_ray_hits','max_rays_traced','is_sunshape','is_surface_errors']:
            psnew.__setattr__(attr, copy.deepcopy(self.__getattribute__(attr)))
        psnew.stages = [PySolTrace.Stage(st.id) for st in self.stages]
        for i in range(len(self.stages)):
            self.stages[i].copy(psnew.stages[i])
        psnew.optics = [PySolTrace.Optics(op.id) for op in self.optics]
        for i in range(len(self.optics)):
            self.optics[i].copy(psnew.optics[i])
        psnew.sun = PySolTrace.Sun()
        self.sun.copy(psnew.sun)

        #re-link optics and elements
        opt_map = {}
        for opt in psnew.optics:
            opt_map[opt.id] = opt

        for i,stage in enumerate(self.stages):
            for j,element in enumerate(stage.elements):
                psnew.stages[i].elements[j].optic = opt_map[element.optic.id]

        return psnew

    def Create(self, pdll, p_data):
        """
        Create soltrace context from data structures
        """

        self.sun.Create(pdll, p_data)
        for opt in self.optics:
            opt.Create(pdll, p_data)
        for stage in self.stages:
            stage.Create(pdll, p_data) 

    def add_optic(self, optic_name : str):
        """
        Instantiates a new PySolTrace.Optics object, adding it to the optics list. 
        This method does not set optics properties, which instead is done using the Optics.Create method.

        Parameters
        ----------
        optic_name : string
            Unique name for this Optics instance.

        Returns
        ----------
        Optics 
            Reference to the Optics object that was just created.
        """

        new_opt_id = len(self.optics)

        self.optics.append( PySolTrace.Optics(new_opt_id ) )
        self.optics[-1].name = optic_name

        return self.optics[-1]  #return the last object in the list, which was the one just created

    def delete_optic(self, optic_id : int) -> int:
        """
        Delete Optics instance. The optics object is removed from the PySolTrace.optics list and 
        from the SolTrace context.

        Parameters
        ----------
        optic_id : int
            ID associated with the optics to be deleted

        Returns
        ----------
        int 
            1 if successful, 0 otherwise
        """

        # find the appropriate optic
        for opt in self.optics: 
            if optic_id == opt.id:
                # clear it from the optics array 
                self.optics.remove(opt)
                # Remove from the soltrace context
                # self._pdll.st_delete_optic.restype = c_int 
                # return self._pdll.st_delete_optic(c_void_p(self._p_data), c_uint32(optic_id) )
        
        # If reaching this point, the optic id was not found
        return 0

    def add_sun(self):
        """
        Instantiates a PySolTrace.Sun object and associates it with the PySolTrace.sun member.
        This does not create or modify the Sun data in the SolTrace context.

        Returns
        ----------
        PySolTrace.Sun 
            Reference to newly created Sun instance.
        """
        ## Object containing Sun class data
        self.sun = PySolTrace.Sun()
        return self.sun

    def add_stage(self):
        """
        Adds a new Stage instance to the PySolTrace.stages list. The Stage ID is automatically generated based on the number 
        of current stages.

        Returns
        ----------
        PySolTrace.Stage 
            Reference to the newly created Stage object. 
        """
        new_st_id = len(self.stages)

        self.stages.append( PySolTrace.Stage( new_st_id ) )

        return self.stages[-1]

    def delete_stage(self, stage_id : int) -> int:
        """
        Delete Stage instance. The stage object is removed from the PySolTrace.stages list and 
        from the SolTrace context.

        Parameters
        ----------
        stage_id : int
            ID associated with the stage to be deleted

        Returns
        ----------
        int 
            1 if successful, 0 otherwise 
        """
        # find the appropriate optic
        for st in self.stages: 
            if stage_id == st.id:
                # clear it from the optics array 
                self.stages.remove(st)
                # Remove from the soltrace context
                # self._pdll.st_delete_stage.restype = c_int 
                # return self._pdll.st_delete_stage(c_void_p(self._p_data), c_uint32(stage_id) )
        
        # If reaching this point, the stage id was not found
        return 0

    def __load_dll(self):
        cwd = os.getcwd()
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            ## loaded SolTrace library of exported functions
            pdll = CDLL(cwd + "/coretrace_api.dll")
            # print("Loaded win32")
            #pdll = CDLL(cwd + "/coretraced.dll") # for debugging
        elif sys.platform == 'darwin':
            pdll = CDLL(cwd + "/coretrace_api.dylib")  # Never tested
        elif sys.platform.startswith('linux'):
            pdll = CDLL(cwd +"/coretrace_api.so")  
        else:
            print( 'Platform not supported ', sys.platform)
        return pdll

    def run(self, seed : int = -1, as_power_tower = False, nthread=1, thread_id=0):
        """
        Run SolTrace simulation. 

        If calling this function in multithread mode, note that the run() function
        **must** be called inside an import guard, e.g.:
        > if __name__ == "__main__":
        >    mypst_obj.run(...)
        Otherwise, you'll receive an error.

        Parameters
        ----------
        seed : int
            Seed for random number generator. [-1] for random seed.
        as_power_tower : bool
            Flag indicating simulation should be processed as power 
            tower / central receiver type, with corresponding efficiency adjustments.
        
        Returns
        ----------
        int 
            Simulation return value
        """

        pdll = self.__load_dll()
        
        if nthread == 1:
            
            # Create an instance of soltrace in memory
            pdll.st_create_context.restype = c_void_p
            p_data = pdll.st_create_context()

            self.Create(pdll, p_data)

            pdll.st_sim_errors.restype = c_int 
            pdll.st_sim_errors(c_void_p(p_data), c_int(1 if self.is_sunshape else 0), c_int(1 if self.is_surface_errors else 0))

            pdll.st_sim_params.restype = c_int 
            pdll.st_sim_params(c_void_p(p_data), c_int(int(self.num_ray_hits)), c_int(int(self.max_rays_traced)))

            if thread_id == 0:
                tstart = time.time()

            pdll.st_sim_run.restype = c_int 
            res = pdll.st_sim_run( c_void_p(p_data), c_uint16(seed), 
                    c_bool(as_power_tower), api_callback, thread_id)
            
            if thread_id == 0:
                print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

            # Collect simulation output, including raw ray data and sunbox stats
            self.raydata = self.__get_ray_dataframe(pdll,p_data)
            self.sunstats = self.__get_sun_stats(pdll, p_data)
            # Compute and save power per ray
            self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni

            pdll.st_free_context.restype = c_bool
            pdll.st_free_context(c_void_p(p_data))

            return res
        else:
            seeds = [seed + i if seed>0 else random.randint(1,int(1e9)) for i in range(nthread)]

            P = [[self.copy(), seeds[i], i+1] for i in range(nthread)]

            # modify the number of rays to match the required totals
            nrpt = int(float(self.num_ray_hits)/float(nthread))
            mrpt = int(float(self.max_rays_traced)/float(nthread))

            for p in P:
                p[0].num_ray_hits = nrpt
                p[0].max_rays_traced = mrpt
                if p == P[0]:
                    p[0].num_ray_hits += int(float(self.num_ray_hits) % float(nthread))
                    p[0].max_rays_traced += int(float(self.max_rays_traced) % float(nthread))

            pool = multiprocessing.Pool(nthread)
            print("Launching {:d} threads...".format(nthread))
            tstart = time.time()
            res = pool.starmap_async(_thread_func, P)
            pool.close()
            pool.join()
            print("\nSimulation complete. Total simulation time {:.2f} seconds.".format(time.time()-tstart))

            # Modify the ray number for threads 2+ to avoid duplication
            dfs = [r[0] for r in res.get()]
            rstart = int(dfs[0].iloc[-1].number)
            if len(dfs)>1:
                for d in dfs[1:]:
                    d.number = d.number+rstart
                    rstart = d.number.iloc[-1]

            self.raydata = pd.concat(dfs)
            self.sunstats = res.get()[0][1]  #take the first thread result
            # add up all the sunrays from all threads
            srct = 0
            for r in res.get():
                srct += r[1]['nsunrays']
            self.sunstats['nsunrays'] = srct

            # Compute and save power per ray [W]
            self.powerperray = (self.sunstats['xmax']-self.sunstats['xmin'])*(self.sunstats['ymax'] - self.sunstats['ymin']) / self.sunstats['nsunrays'] * self.dni
            
            return 1

    def get_num_intersections(self, pdll, p_data) -> int:
        """
        [Post simulation] Get the number of ray intersections detected in the simulation.

        Returns
        ----------
        int 
            Number of intersections
        """

        if p_data == 0:
            return 0

        pdll.st_num_intersections.restype = c_int
        return pdll.st_num_intersections(c_void_p(p_data))

    def __get_sun_stats(self, pdll, p_data):
        """
        Get information on the sun box. 

        Returns
        ----------
        dict      
            Keys in the return dictionary are:
            'xmin' --> Minimum x extent of the bounding box for hit testing 
            'xmax' --> Maximum x extent of the bounding box for hit testing 
            'ymin' --> Minimum y extent of the bounding box for hit testing 
            'ymax' --> Maximum y extent of the bounding box for hit testing 
            'nsunrays' --> Number of sun rays simulated
        """
        if p_data == 0:
            raise "SolTrace context not assigned"
        
        xmin = (c_number)()
        xmax = (c_number)()
        ymin = (c_number)()
        ymax = (c_number)()
        nsunrays = (c_int)()

        pdll.st_sun_stats.restype = c_int 
        pdll.st_sun_stats(c_void_p(p_data), pointer(xmin), pointer(xmax), pointer(ymin), pointer(ymax), pointer(nsunrays))

        return {
            'xmin':float(xmin.value),
            'xmax':float(xmax.value),
            'ymin':float(ymin.value),
            'ymax':float(ymax.value),
            'nsunrays':int(nsunrays.value),
        }

    def __get_ray_dataframe(self, pdll, p_data):
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
        if p_data == 0:
            raise "SolTrace context not assigned"
        
        data = {}

        n_int = self.get_num_intersections(pdll, p_data)
        data['loc_x'] = (c_number*n_int)()
        data['loc_y'] = (c_number*n_int)()
        data['loc_z'] = (c_number*n_int)()

        pdll.st_locations.restype = c_int 
        pdll.st_locations(c_void_p(p_data), pointer(data['loc_x']), pointer(data['loc_y']), pointer(data['loc_z']))

        data['cos_x'] = (c_number*n_int)()
        data['cos_y'] = (c_number*n_int)()
        data['cos_z'] = (c_number*n_int)()

        pdll.st_cosines.restype = c_int 
        pdll.st_cosines(c_void_p(p_data), pointer(data['cos_x']), pointer(data['cos_y']), pointer(data['cos_z']))

        data['element'] = (c_int*n_int)()

        pdll.st_elementmap.restype = c_int 
        pdll.st_elementmap(c_void_p(p_data), pointer(data['element']))

        data['stage'] = (c_int*n_int)()

        pdll.st_stagemap.restype = c_int 
        pdll.st_stagemap(c_void_p(p_data), pointer(data['stage']))

        data['number'] = (c_int*n_int)()

        pdll.st_raynumbers.restype = c_int 
        pdll.st_raynumbers(c_void_p(p_data), pointer(data['number']))

        for key in data.keys():
            data[key] = list(data[key])

        df = pd.DataFrame(data)

        return df

    def validate(self) -> bool:
        """
        Validate that the current SolTrace context has been configured correctly, and
        that commonly missed inputs are specified.

        Returns
        ----------
        bool 
            False if error is identified, True otherwise
        """



        return True


    # /* utility transform/math functions */
    def util_calc_euler_angles(self, origin, aimpoint, zrot):
        """
        Calculate the Euler angles associated with a given origin, aimpoint, and z-axis rotation. 

        Parameters
        ----------
        origin : [float,*3]
            Origin of the coordinate system
        aimpoint : [float,*3]
            Aimpoint of the vector originating at the origin
        zrot : float
            Rotation around the z-axis coordinate (degr)

        Returns
        ----------
        list
            Calculated Euler angles (rad)
        """

        a_origin = (c_number*3)()
        a_aimpoint = (c_number*3)()
        a_euler = (c_number*3)()
        a_origin[:] = origin
        a_aimpoint[:] = aimpoint

        pdll = self.__load_dll()

        pdll.st_calc_euler_angles.restype = c_void_p
        pdll.st_calc_euler_angles(pointer(a_origin), pointer(a_aimpoint), c_number(zrot), pointer(a_euler))

        return list(a_euler)

    def util_transform_to_local(self, posref, cosref, origin, rreftoloc):
        """
        Perform coordinate transformation from reference system to local system.
        
        Parameters
        ----------
        PosRef : [float,]*3
            X,Y,Z coordinates of ray point in reference system
        CosRef : [float,]*3
            Direction cosines of ray in reference system
        Origin : [float,]*3
            X,Y,Z coordinates of origin of local system as measured in reference system
        RRefToLoc 
            Rotation matrices required for coordinate transform from reference to local
        
        Returns
        ----------
        (dict)  Keys in return dictionary include:
            posloc : ([float,]*3) X,Y,Z coordinates of ray point in local system
            cosloc : ([float,]*3) Direction cosines of ray in local system
        """
        pdll = self.__load_dll()
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

        pdll.st_transform_to_local.restype = c_void_p
        pdll.st_transform_to_local(pointer(a_posref), pointer(a_cosref), pointer(a_origin), pointer(a_rreftoloc), pointer(posloc), pointer(cosloc))

        return {'cosloc':list(cosloc), 'posloc':list(posloc)}

    def util_transform_to_ref(self, posloc, cosloc, origin, rloctoref):
        """
        Perform coordinate transformation from local system to reference system.
        
        Parameters
        ----------
        PosLoc : [float,]*3
            X,Y,Z coordinates of ray point in local system
        CosLoc : [float,]*3
            Direction cosines of ray in local system
        Origin : [float,]*3
            X,Y,Z coordinates of origin of local system as measured in reference system
        RLocToRef 
            Rotation matrices required for coordinate transform from local to reference
            -- inverse of reference to local transformation
        
        Returns
        ----------
        dict  
            Keys in return dictionary include:
            posref : ([float,]*3) X,Y,Z coordinates of ray point in reference system
            cosref : ([float,]*3) Direction cosines of ray in reference system
        """
        pdll = self.__load_dll()

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

        pdll.st_transform_to_reference.restype = c_void_p
        pdll.st_transform_to_reference(pointer(a_posloc), pointer(a_cosloc), pointer(a_origin), pointer(a_rloctoref), pointer(posref), pointer(cosref))

        return {'cosref':list(cosref), 'posref':list(posref)}

    def util_matrix_vector_mult(self, m, v):
        """
        Perform multiplication of a 3x3 matrix and a length-3 vector, returning the result vector.

        Parameters
        ----------
        m : list
            m[3][3] - a 3x3 matrix
        v : list
            v[3] - a list, length 3

        Returns
        ----------
        list
            m x v [3]
        """
        pdll = self.__load_dll()

        a_m = (c_number*9)()
        a_v = (c_number*3)()
        a_mv = (c_number*3)()
        a_m[:] = sum([], m)
        a_v[:] = v 

        pdll.st_matrix_vector_mult.restype = c_void_p
        pdll.st_matrix_vector_mult(pointer(a_m), pointer(a_v), pointer(a_mv))

        return list(a_mv)

    def util_calc_transforms(self, euler):
        """
        Calculate matrix transforms

        Parameters
        ----------
        euler : [float,]*3
            Euler angles

        Returns
        ----------
        (dict) A dictionary containing the keys:
            rreftoloc : Transformation matrix from Reference to Local system
            rloctoref : Transformation matrix from Local to Reference system
        """

        pdll = self.__load_dll()

        a_euler = (c_number*3)()
        rreftoloc = (c_number*9)()
        rloctoref = (c_number*9)()

        a_euler[:] = euler

        pdll.st_calc_transform_matrices.restype = c_void_p
        pdll.st_calc_transform_matrices(pointer(a_euler), pointer(rreftoloc), pointer(rloctoref))

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

        Parameters
        ----------
        m : [[float,]*3]*3
            Square matrix 3x3 to be transposed.

        Returns
        ----------
        [[float,]*3]*3
            Square matrix 3x3, transpose of 'm'
        """

        pdll = self.__load_dll()

        a_m = (c_number*9)()
        a_m[:] = sum([], m)
        output = (c_number*9)()

        pdll.st_matrix_transpose.restype = c_void_p
        pdll.st_matrix_transpose(pointer(a_m), pointer(output))

        a_output = []
        for i in range(0,10,3):
            a_output.append(output[i:i+3])

        return a_output

    def util_rotation_arbitrary(self, theta, axis, axloc, pt):
        """
        Rotation of a point 'pt' about an arbitrary axis with direction 'axis' centered at point 'axloc'. 
        The point is rotated through 'theta' radians.

        Parameters
        ----------
        theta : float
            Angle of rotation (rad)
        axis : Point()
            Vector (x=i,y=j,z=k) indicating direction of axis for rotation
        axloc : Point()
            Location of the axis origin
        pt : Point()
            Location of the point that is to be rotated

        Returns
        -----------
        Point
            Point after rotation

        """

        a = axloc.x     #Point through which the axis passes
        b = axloc.y
        c = axloc.z
        x = pt.x        #Point that we're rotating
        y = pt.y
        z = pt.z
        u = axis.x        #Direction of the axis that we're rotating about
        v = axis.y
        w = axis.z

        
        sinth = math.sin(theta)
        costh = math.cos(theta)

        fin = Point()

        fin.x = (a*(v*v+w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1.-costh) + x*costh + (-c*v + b*w - w*y + v*z)*sinth
        fin.y = (b*(u*u+w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1.-costh) + y*costh + (c*u - a*w + w*x - u*z)*sinth
        fin.z = (c*(u*u+v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1.-costh) + z*costh + (-b*u + a*v - v*x + u*y)*sinth

        return fin

    def util_calc_unitvect(self, vect):
        if type(vect) == list:
            v = Point()
            v.x = vect[0]
            v.y = vect[1]
            v.z = vect[2]
        else:
            v = vect

        vect_mag = math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
        nvect = Point()
        nvect.x = v.x / vect_mag
        nvect.y = v.y / vect_mag
        nvect.z = v.z / vect_mag

        if type(vect) == list:
            return [nvect.x, nvect.y, nvect.z]
        else:
            return nvect

    def util_calc_zrot_azel(self, vect) -> float:
        """
        Compute the z-rotation of a vector, assuming the vector's deviation from (0,0,1) 
        has been realized using azimuth-elevation transforms.

        Parameters
        ----------
        vect : (list OR Point)
            i,j,k components of a vector

        Returns
        ----------
        float
            Computed z-rotation (degrees)
        """
        if isinstance(vect, List):
            vect_i, vect_j, vect_k = vect
        elif isinstance(vect, Point):
            vect_i = vect.x 
            vect_j = vect.y 
            vect_k = vect.z 
        else:
            raise TypeError("Function expects 'vect' of type List or Point")

        az = math.atan2(vect_i,vect_j)
        az = (az + 2.*math.pi) if az < 0. else az 

        el = math.asin(vect_k)

        #Calculate Euler angles
        alpha = math.atan2(vect_i, vect_k);        #Rotation about the Y axis
        bsign = 1 if vect_j > 0. else -1
        beta = -bsign*math.acos( ( math.pow(vect_i,2) + math.pow(vect_k,2) )/ 
                        max(math.sqrt(math.pow(vect_i,2) + math.pow(vect_k,2)), 1.e-8) )    #Rotation about the modified X axis

        #Calculate the modified axis vector
        modax = Point(math.cos(alpha), 0., - math.sin(alpha))  

        #Rotation references - axis point. Set as origin
        axpos = Point(0., 0., 0.)
        #sp_point to rotate. lower edge of heliostat
        pbase = Point(0., -1., 0.)
        
        #Rotated point
        protv = self.util_rotation_arbitrary(beta, modax, axpos, pbase).unitize()

        #Azimuth/elevation reference vector (vector normal to where the base of the heliostat should be)
        azelref = Point()   
        azelref.x = math.sin(az)*math.sin(el)
        azelref.y = math.cos(az)*math.sin(el)
        azelref.z = -math.cos(el)

        # the sign of the rotation angle is determined by whether the 'k' component of the cross product
        # vector is positive or negative. 
        cp = Point()
        cp.x = protv.y*azelref.z - protv.z*azelref.y
        cp.y = protv.z*azelref.x - protv.x*azelref.z
        cp.z = protv.x*azelref.y - protv.y*azelref.x

        gamma = math.asin( cp.radius() )
        gsign = (1 if cp.z > 0. else -1.) * (1 if vect_j > 0. else -1.)

        return gamma * gsign * 180./math.pi

    def write_soltrace_input_file(self, path : str):
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
            fout.write(
                "# SOLTRACE VERSION 2012.7.6 INPUT FILE -- GENERATED BY soltrace-api v{:s} on {:02d}/{:02d}/{:04d} at {:02d}:{:02d}:{:02d}\n".format( 
                    "v000", dt.day, dt.month, dt.year, dt.hour, dt.minute, dt.second)
            )

            #------------------- Sun shape
            fout.write( "SUN\tPTSRC\t{:d}\tSHAPE\t{:1s}\tSIGMA\t{:f}\tHALFWIDTH\t{:f}\n".format( 0, self.sun.shape, self.sun.sigma, self.sun.sigma) )
            fout.write( "XYZ\t{:f}\t{:f}\t{:f}\tUSELDH\t{:d}\tLDH\t{:f}\t{:f}\t{:f}\n".format( self.sun.position.x, self.sun.position.y, self.sun.position.z, 0, 0., 0., 0.) )
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
                            opt.dist_type,3, 1, 4, opt.reflectivity, opt.transmissivity, opt.slope_error, opt.spec_error , 0., 0., 0., 0., 0., 0. ) 
                        )

            #------------------- loop through Stage list
            fout.write( "STAGE LIST COUNT\t{:d}\n".format( len(self.stages)) )
            for stage in self.stages:
                fout.write( 
                    "STAGE\tXYZ\t{:f}\t{:f}\t{:f}\tAIM\t{:f}\t{:f}\t{:f}\tZROT\t{:f}\tVIRTUAL\t{:d}\tMULTIHIT\t{:d}\tELEMENTS\t{:d}\tTRACETHROUGH\t{:d}\n".format( 
                        stage.position.x, stage.position.y,stage.position.z, stage.aim.x, stage.aim.y, stage.aim.z, stage.zrot, int(stage.is_virtual), int(stage.is_multihit), len(stage.elements), int(stage.is_tracethrough)) 
                        ) 

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
# -----------------------------------------------------------------------------------------------------------------------------

def loaddll():
    cwd = os.getcwd()
    pdll = CDLL(cwd + "/coretrace_api.dll")


if __name__ == "__main__":
    # p1 = PySolTrace()
    # p1.add_sun()

    # p2 = p1.copy()

    # print(p1,p2)

    import time
    tstart = time.time()

    for i in range(1000):
        loaddll()

    print(time.time() - tstart)
