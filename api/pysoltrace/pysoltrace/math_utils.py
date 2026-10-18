import math
from collections.abc import Callable
import numpy as np

from pysoltrace import Point

# Freedman–Diaconis rule
def freedman_diaconis_np(arr, fudge: float = 1.0):
    min_val, iqr25, iqr75, max_val = np.quantile(arr, [0, .25, .75, 1])
    bin_width = 2 * (iqr75 - iqr25) * len(arr) ** (-1/3)
    return int(np.ceil(fudge * (max_val - min_val) / bin_width))

def euler_angles(origin:   Point | list | np.ndarray,
                 aimpoint: Point | list | np.ndarray,
                 zrot:     float) -> Point | list | np.ndarray:
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
        _arg_origin_type = type(origin)
        if _arg_origin_type == np.ndarray:
            _origin = origin
        elif _arg_origin_type == list:
            _origin = np.array(origin)
        elif _arg_origin_type == Point:
            _origin = np.array(origin.as_list())
        else:
            raise ValueError(f'Can\'t convert {_arg_origin_type} to numpy array.')

        _arg_aimpoint_type = type(aimpoint)
        if _arg_aimpoint_type == np.ndarray:
            _aimpoint = aimpoint
        elif _arg_aimpoint_type == list:
            _aimpoint = np.array(aimpoint)
        elif _arg_aimpoint_type == Point:
            _aimpoint = np.array(aimpoint.as_list())
        else:
            raise ValueError(f'Can\'t convert {_arg_aimpoint_type} to numpy array.')

        # cast to float64 to avoid issues with numpy division
        _origin = _origin.astype(np.float64)
        _aimpoint = _aimpoint.astype(np.float64)

        # This duplicates the built-in function but uses numpy operators directly
        dv = _aimpoint - _origin
        d = math.sqrt(sum(dv**2))
        if d == 0:
            return
        dv /= d
            # math.atan2(-dv[0], math.sqrt(dv[1]**2 + dv[2]**2)),
            # math.atan2(-dv[1], dv[2]),
        euler = np.array([
            math.atan2(dv[0], dv[2]),
            math.asin(dv[1]),
            zrot * 0.017453292519943295 # acos(-1)/180.0
        ])

        if _arg_origin_type == list:
            return euler.tolist()
        elif _arg_origin_type == Point:
            return Point.from_list(euler)
        return euler

def transform_to_local(posref:    np.ndarray,
                       cosref:    np.ndarray,
                       origin:    np.ndarray,
                       rreftoloc: np.ndarray):
    """
    Perform coordinate transformation from reference system to local system.

    Parameters
    ----------
    PosRef : numpy.array([float,]*3)
        X,Y,Z coordinates of ray point in reference system
    CosRef : numpy.array([float,]*3)
        Direction cosines of ray in reference system
    Origin : numpy.array([float,]*3)
        X,Y,Z coordinates of origin of local system as measured in reference system
    RRefToLoc : numpy.array([float,]*3)
        Rotation matrices required for coordinate transform from reference to local

    Returns
    ----------
    (dict)  Keys in return dictionary include:
        posloc : ([float,]*3) X,Y,Z coordinates of ray point in local system
        cosloc : ([float,]*3) Direction cosines of ray in local system
    """
    assert type(rreftoloc) == type(np.array([]))
    assert rreftoloc.shape == (3,3)

    # This duplicates the built-in function but uses numpy operators directly
    posdum = posref - origin
    rreftoloc_m = rreftoloc.reshape((3,3))
    posloc = np.dot(rreftoloc_m, posdum.T).T
    cosloc = np.dot(rreftoloc_m, cosref.T).T

    return {'cosloc': cosloc, 'posloc': posloc}

def util_transform_to_ref(posloc:    np.ndarray,
                          cosloc:    np.ndarray,
                          origin:    np.ndarray,
                          rloctoref: np.ndarray):
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
    assert type(rloctoref) == type(np.array([]))
    assert rloctoref.shape == (3,3)

    posdum = np.dot(rloctoref, posloc)
    cosref = np.dot(rloctoref, cosloc, cosref)
    posref = posdum + origin

    return {'cosref': cosref, 'posref': posref}
        

def matrix_vector_mult(m: np.ndarray, v: np.ndarray):
    """
    Perform multiplication of a 3x3 matrix and a length-3 vector, returning the result vector.

    Parameters
    ----------
    m : array
        m[3][3] - a 3x3 matrix
    v : array
        v[3] - a list, length 3

    Returns
    ----------
    list
        m x v [3]
    """
    assert type(m) == type(np.array([]))
    assert m.shape == (3,3)
    assert type(v) == type(np.array([]))
    assert v.shape == (3,)

    return np.dot(m, v)

def euler_transforms(euler: Point | list):
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

    assert isinstance(euler, (Point, list)) or hasattr(euler, '__getitem__'), 'Euler angle values must be a Point or list.'
    if isinstance(euler, list):       assert len(euler) == 3,                 f'Must have 3 Euler angles, not {len(euler)}.'
    if isinstance(euler, np.ndarray): assert euler.shape == (3,),             f'Must have 3 Euler angles, not {euler.shape}.'

    Alpha = euler[0]
    Beta  = euler[1]
    Gamma = euler[2]
    CosAlpha = np.cos(Alpha)
    CosBeta  = np.cos(Beta)
    CosGamma = np.cos(Gamma)
    SinAlpha = np.sin(Alpha)
    SinBeta  = np.sin(Beta)
    SinGamma = np.sin(Gamma)

    # Fill in elements of the transformation matrix as per 
    # Spencer and Murty paper page 673 equation (2)
    rreftoloc = np.array([[CosAlpha * CosGamma + SinAlpha * SinBeta * SinGamma,
                           CosAlpha * SinGamma - SinAlpha * SinBeta * CosGamma,
                           SinAlpha * CosBeta],
                          [-CosBeta * SinGamma, CosBeta * CosGamma, SinBeta],
                          [-SinAlpha * CosGamma + CosAlpha * SinBeta * SinGamma,
                           -SinAlpha * SinGamma - CosAlpha * SinBeta * CosGamma,
                            CosAlpha * CosBeta]])

    return {'rreftoloc':rreftoloc, 'rloctoref':rreftoloc.T}

def arbitrary_rotation(theta: float,
                       axis:  Point | list,
                       axloc: Point | list,
                       pt:    Point | list):
    """    
    Rotation of a point 'pt' about an arbitrary axis with direction 'axis' centered at point 'axloc'.
    The point is rotated through 'theta' radians.

    Parameters
    ----------
    theta : float
        Angle of rotation (rad)
    axis : list or np.ndarray
        Vector (x=i,y=j,z=k) indicating direction of axis for rotation
    axloc : list or np.ndarray
        Location of the axis origin
    pt : list or np.ndarray
        Location of the point that is to be rotated

    Returns
    -----------
    np.ndarray
        Point after rotation
    """

    a = axloc[0]     #Point through which the axis passes
    b = axloc[1]
    c = axloc[2]
    x = pt[0]        #Point that we're rotating
    y = pt[1]
    z = pt[2]
    u = axis[0]        #Direction of the axis that we're rotating about
    v = axis[1]
    w = axis[2]

    sinth = math.sin(theta)
    costh = math.cos(theta)

    return Point((a*(v*v+w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1.-costh) + x*costh + (-c*v + b*w - w*y + v*z)*sinth,
                 (b*(u*u+w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1.-costh) + y*costh + ( c*u - a*w + w*x - u*z)*sinth,
                 (c*(u*u+v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1.-costh) + z*costh + (-b*u + a*v - v*x + u*y)*sinth)

def unitize(vect: Point | list):
    """
    Scales a vector to have total magnitude of 1

    Parameters
    ----------
    vect : list | Point
        list or Point containing the vector

    Returns
    ----------
    list | Point
        Unitized vector of type list or Point, depending on input type
    """
    _arg_type = type(vect)
    if _arg_type == list or _arg_type == np.ndarray:
        v = Point.from_list(vect)
    elif _arg_type == Point:
        v = vect
    else:
        raise ValueError(f'Can\'t convert {_arg_type} to Point.')

    v.unitize(inplace = True)

    if _arg_type == list:
        return v.as_list()
    elif _arg_type == np.ndarray:
        return np.array(v.as_list())
    else:
        return v

def zrot_from_azel(vect: Point | list) -> float:
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
    if isinstance(vect, (Point, list)) or hasattr(vect, '__getitem__'):
        vect_i, vect_j, vect_k = vect[0], vect[1], vect[2]
    else:
        raise TypeError("Function expects 'vect' of type List or Point")

    az = math.atan2(vect_i, vect_j)
    az = (az + 2.*math.pi) if az < 0. else az

    el = math.asin(vect_k)

    #Calculate Euler angles
    alpha = math.atan2(vect_i, vect_k);                                        #Rotation about the Y axis
    bsign = 1 if vect_j > 0. else -1
    beta = -bsign * math.acos( 
            (math.pow(vect_i,2) + math.pow(vect_k,2)) /
            max(math.sqrt(math.pow(vect_i,2) + math.pow(vect_k,2)), 1.e-8))    #Rotation about the modified X axis

    #Calculate the modified axis vector
    modax = Point(math.cos(alpha), 0., - math.sin(alpha))

    #Rotation references - axis point. Set as origin
    axpos = Point(0., 0., 0.)
    #sp_point to rotate. lower edge of heliostat
    pbase = Point(0., -1., 0.)

    #Rotated point
    protv = arbitrary_rotation(beta, modax, axpos, pbase).unitize()

    #Azimuth/elevation reference vector (vector normal to where the base of the heliostat should be)
    azelref = Point(
        math.sin(az)*math.sin(el),
        math.cos(az)*math.sin(el),
        -math.cos(el)
    )
    # the sign of the rotation angle is determined by whether the 'k' component of the cross product
    # vector is positive or negative.
    cp = protv @ azelref

    gamma = math.asin(cp.radius())
    gsign = (1 if cp.z > 0. else -1.) * (1 if vect_j > 0. else -1.)

    return gamma * gsign * 180./math.pi

def get_unstager(pos: Point | list | np.ndarray, aim: Point | list | np.ndarray, zrot: float) -> Callable | None:
    assert isinstance(pos, (Point, list, np.ndarray)), f'pos must be Point, list, or array, not {type(pos)}'
    assert isinstance(aim, (Point, list, np.ndarray)), f'aim must be Point, list, or array, not {type(aim)}'
    assert len(pos) == 3, f'pos must have 3 elements, not {len(pos)}'
    assert len(aim) == 3, f'aim must have 3 elements, not {len(aim)}'

    _pos = pos
    if type(_pos) == list:
        _pos = np.array(_pos)
    elif type(_pos) == Point:
        _pos = _pos.as_array()

    _aim = aim
    if type(_aim) == list:
        _aim = np.array(_aim)
    elif type(_aim) == Point:
        _aim = _aim.as_array()

    euler = euler_angles(_pos, _aim, zrot)
    significant_rotation = np.linalg.norm(euler) > 1e-7
    significant_translation = np.linalg.norm(_pos) > 1e-7

    if significant_rotation:
        transform = euler_transforms(euler)['rloctoref'] # ['rreftoloc']
        if significant_translation:
            unstager = lambda v: transform @ v + _pos
        else:
            unstager = lambda v: transform @ v
    elif significant_translation:
        unstager = lambda v: v + _pos
    else:
        unstager = None
    return unstager
