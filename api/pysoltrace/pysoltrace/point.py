from __future__ import annotations

import functools
from typing import Literal
from numpy import array, ndarray

class Point:
    """
    A simple class to manage points in Cartesian coordinates.
    """

    __slots__ = ('x', 'y', 'z', '__i',)

    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
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

        self.__i = 0
        return
    
    def copy(self) -> Point:
        pnew = Point()
        pnew.x = self.x
        pnew.y = self.y
        pnew.z = self.z
        return pnew

    def __repr__(self) -> str:
        return f"[{self.x:.2f}, {self.y:.2f}, {self.z:.2f}]"
        # return "[{:f}, {:f}, {:f}]".format(self.x, self.y, self.z)

    def __format__(self, format_spec):
        if not format_spec: format_spec = '.2f'
        chunk = '{:' + format_spec + '}'
        return f"[{chunk}, {chunk}, {chunk}]".format(self.x, self.y, self.z)

    def __bool__(self):
        return self.x != 0 or self.y != 0 or self.z != 0

    @property
    def no_zeros(self):
        return self.x != 0 and self.y != 0 and self.z != 0

    def __float__(self):
        return self.radius()

    def __int__(self):
        return int(self.radius())

    def __iter__(self) -> Point:
        # reset iteration each time
        self.__i = 0 
        return self

    def __next__(self) -> float | int:
        if self.__i < 3:
            self.__i += 1
            return self[self.__i - 1]
        raise StopIteration
    
    def __getitem__(self, i) -> float | int:
        if i == 0:   return self.x
        elif i == 1: return self.y
        elif i == 2: return self.z
        else:        raise IndexError

    def __setitem__(self, i, value: float | int) -> None:
        if i == 0:   self.x = value
        elif i == 1: self.y = value
        elif i == 2: self.z = value
        else:        raise IndexError

    def __len__(self):
        return 3

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)

    def __abs__(self):
        return self.radius()

    def __add__(self, obj: Point | float | int | list) -> Point:
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
            return Point(self.x + obj.x, 
                         self.y + obj.y,
                         self.z + obj.z)
        elif isinstance(obj, (float, int)):
            return Point(self.x + obj, 
                         self.y + obj,
                         self.z + obj)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.x + obj[0], 
                         self.y + obj[1],
                         self.z + obj[2])
        return NotImplemented

    def __radd__(self, obj: Point | float | int | list) -> Point:
        return self + obj

    def __iadd__(self, obj: Point | float | int | list) -> Point:
        if isinstance(obj, Point):
            self.x += obj.x
            self.y += obj.y
            self.z += obj.z
        elif isinstance(obj, (float, int)):
            self.x += obj
            self.y += obj
            self.z += obj
        elif isinstance(obj, (list, ndarray)):
            self.x += obj[0]
            self.y += obj[1]
            self.z += obj[2]
        else: return NotImplemented
        return self

    def __sub__(self, obj: Point | float | int | list) -> Point:
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
            return Point(self.x - obj.x, 
                         self.y - obj.y, 
                         self.z - obj.z)
        elif isinstance(obj, (float, int)):
            return Point(self.x - obj, 
                         self.y - obj, 
                         self.z - obj)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.x - obj[0], 
                         self.y - obj[1],
                         self.z - obj[2])
        return NotImplemented

    def __rsub__(self, obj: Point | float | int | list) -> Point:
        if isinstance(obj, Point):
            return Point(obj.x - self.x, 
                         obj.y - self.y, 
                         obj.z - self.z)
        elif isinstance(obj, (float, int)):
            return Point(obj - self.x, 
                         obj - self.y, 
                         obj - self.z)
        elif isinstance(obj, (list, ndarray)):
            return Point(obj[0] - self.x, 
                         obj[1] - self.y,
                         obj[2] - self.z)
        return NotImplemented

    def __isub__(self, obj: Point | float | int | list) -> Point:
        if isinstance(obj, Point):
            self.x -= obj.x
            self.y -= obj.y
            self.z -= obj.z
        elif isinstance(obj, (float, int)):
            self.x -= obj
            self.y -= obj
            self.z -= obj
        elif isinstance(obj, (list, ndarray)):
            self.x -= obj[0]
            self.y -= obj[1]
            self.z -= obj[2]
        else: return NotImplemented
        return self

    def __mul__(self, obj: Point | list | float | int) -> Point:
        """
        Multiplies the current point coordinate values.

        Parameters
        ==========
        obj : variant
            If obj = (Point | list), multiplies each component of the current point by obj
            If obj = (float), multiplies each component of the current point by obj
        Returns
        =======
        Point
            Reference to this point
        """
        if isinstance(obj, Point):
            return Point(self.x * obj.x,
                         self.y * obj.y,
                         self.z * obj.z)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.x * obj[0],
                         self.y * obj[1],
                         self.z * obj[2])
        elif isinstance(obj, (float, int)):
            return Point(self.x * obj, 
                         self.y * obj, 
                         self.z * obj)  
        return NotImplemented

    def __rmul__(self, obj: Point | list | float | int) -> Point:
        return self * obj

    def __imul__(self, obj: Point | list | float | int) -> Point:
        if isinstance(obj, Point):
            self.x *= obj.x
            self.y *= obj.y
            self.z *= obj.z
        elif isinstance(obj, (list, ndarray)):
            self.x *= obj[0]
            self.y *= obj[1]
            self.z *= obj[2]
        elif isinstance(obj, (float, int)):
            self.x *= obj
            self.y *= obj
            self.z *= obj
        else: return NotImplemented
        return self

    def __can_div(self, obj: Point | list | float | int) -> bool:
        if isinstance(obj, Point): 
            return obj.no_zeros
        elif isinstance(obj, (list, ndarray)):
            return obj[0] != 0 and obj[1] != 0 and obj[2] != 0
        elif isinstance(obj, (float, int)):
            return obj != 0
        return False
        
    def __floordiv__(self, obj: Point | list | float | int) -> Point:
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
        if not self.__can_div(obj): raise ValueError(f'No divide by zero. {obj}')
        if isinstance(obj, Point):
            return Point(self.x // obj.x,
                         self.y // obj.y,
                         self.z // obj.z)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.x // obj[0],
                         self.y // obj[1],
                         self.z // obj[2])
        elif isinstance(obj, (float, int)):
            return Point(self.x // obj, 
                         self.y // obj, 
                         self.z // obj)
        return NotImplemented

    def __ifloordiv__(self, obj: Point | list | float | int) -> Point:
        if not self.__can_div(obj): raise ValueError(f'No divide by zero. {obj}')
        if isinstance(obj, Point):
            self.x //= obj.x
            self.y //= obj.y
            self.z //= obj.z
        elif isinstance(obj, (list, ndarray)):
            self.x //= obj[0]
            self.y //= obj[1]
            self.z //= obj[2]
        elif isinstance(obj, (float, int)):
            self.x //= obj
            self.y //= obj
            self.z //= obj
        else: return NotImplemented
        return self

    def __truediv__(self, obj: Point | list | float | int) -> Point:
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
        if not self.__can_div(obj): raise ValueError(f'No divide by zero. {obj}')
        if isinstance(obj, Point):
            return Point(self.x / obj.x,
                         self.y / obj.y,
                         self.z / obj.z)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.x / obj[0],
                         self.y / obj[1],
                         self.z / obj[2])
        elif isinstance(obj, (float, int)):
            return Point(self.x / obj, 
                         self.y / obj, 
                         self.z / obj)
        return NotImplemented

    def __itruediv__(self, obj: Point | list | float | int) -> Point:
        if not self.__can_div(obj): raise ValueError(f'No divide by zero. {obj}')
        if isinstance(obj, Point):
            self.x /= obj.x
            self.y /= obj.y
            self.z /= obj.z
        elif isinstance(obj, (list, ndarray)):
            self.x /= obj[0]
            self.y /= obj[1]
            self.z /= obj[2]
        elif isinstance(obj, (float, int)):
            self.x /= obj
            self.y /= obj 
            self.z /= obj
        else: return NotImplemented
        return self

    def __matmul__(self, obj: Point | list) -> Point:
        """
        implemented as cross porduct between self and obj

        Parameters
        ==========
        obj : Point | list
            calculates cross project
            
        Returns
        =======
        Point
        """
        if isinstance(obj, Point):
            return Point(self.y * obj.z - self.z * obj.y,
                         self.z * obj.x - self.x * obj.z,
                         self.x * obj.y - self.y * obj.x)
        elif isinstance(obj, (list, ndarray)):
            return Point(self.y * obj[2] - self.z * obj[1],
                         self.z * obj[0] - self.x * obj[2],
                         self.x * obj[1] - self.y * obj[0])
        return NotImplemented

    def __rmatmul__(self, obj: ndarray[Literal[3, 3], float]):
        """
        right side matmul implemented as matrix multiplication 
        for transformation matricies
        
        Parameters
        ==========
        obj : np.ndarray(3, 3)
            must be 3x3 matrix, calculates matrix multiplication

        Returns
        =======
        Point
        """
        if isinstance(obj, ndarray) and obj.shape == (3, 3):
            return Point(*(obj @ self.as_array()))
        return NotImplemented

    def __imatmul__(self, obj: Point | list) -> Point:
        """
        see Point.__matmul__
        """
        tmp = self.copy()
        if isinstance(obj, Point):
            self.x = tmp.y * obj.z - tmp.z * obj.y
            self.y = tmp.z * obj.x - tmp.x * obj.z
            self.z = tmp.x * obj.y - tmp.y * obj.x
        elif isinstance(obj, (list, ndarray)):
            self.x = tmp.y * obj[2] - tmp.z * obj[1]
            self.y = tmp.z * obj[0] - tmp.x * obj[2]
            self.z = tmp.x * obj[1] - tmp.y * obj[0]
        else: return NotImplemented
        return self

    def __eq__(self, obj: Point | list) -> bool:
        if isinstance(obj, Point):
            return self.x == obj.x and self.y == obj.y and self.z == obj.z
        elif isinstance(obj, (list, ndarray)):
            return self.x == obj[0] and self.y == obj[1] and self.z == obj[2]
        return NotImplemented

    def __ne__(self, obj: Point | list) -> bool:
        return not self == obj

    def __lt__(self, obj: Point | list) -> bool:
        if isinstance(obj, (list, ndarray)):
            obj = self.from_list(obj)
        if isinstance(obj, Point):
            return self.radius() < obj.radius()
        return NotImplemented

    def __gt__(self, obj: Point | list) -> bool:
        if isinstance(obj, (list, ndarray)):
            obj = self.from_list(obj)
        if isinstance(obj, Point):
            return self.radius() > obj.radius()
        return NotImplemented

    def __le__(self, obj: Point | list) -> bool:
        if isinstance(obj, (list, ndarray)):
            obj = self.from_list(obj)
        if isinstance(obj, Point):
            return self.radius() <= obj.radius()
        return NotImplemented

    def __ge__(self, obj: Point | list) -> bool:
        if isinstance(obj, (list, ndarray)):
            obj = self.from_list(obj)
        if isinstance(obj, Point):
            return self.radius() >= obj.radius()
        return NotImplemented

    def dot(self, obj: Point | list) -> float:
        if isinstance(obj, Point):
            return self.x * obj.x + self.y * obj.y + self.z * obj.z
        elif isinstance(obj, (list, ndarray)):
            return self.x * obj[0] + self.y * obj[1] + self.z * obj[2]
        raise NotImplementedError

    def reduce(self) -> float:
        return self.x + self.y + self.z

    def radius(self) -> float:
        """
        Computes distance between current point and the origin (0,0,0)

        Returns
        =======
        float
            Calculated radius
        """
        # return (self.x*self.x + self.y*self.y + self.z*self.z)**0.5
        return (self.dot(self))**0.5

    def unitize(self, inplace: bool = False) -> Point:
        """
        Converts current point into a unit vector with magnitude 1

        Parameters
        ==========
        inplace : bool = False
            Specifies whether the current point is converted to a unit vector in place, or whether the
            current point remains unchanged and a unitized copy of the vector is returned.
        """
        if (mag := self.radius()) > 0:
            if inplace: self /= mag
            else:       return self / mag
        return self

    def as_list(self) -> list:
        """
        Returns:
        --------
        coordinates as a python list [x,y,z]
        """
        return [self.x, self.y, self.z]

    def as_array(self) -> ndarray[Literal[3,], float]:
        """
        Returns:
        --------
        coordinates as a numpy array [x,y,z]
        """
        return array([self.x, self.y, self.z])

    @staticmethod
    def from_list(l: list) -> Point:
        return Point(l[0], l[1], l[2])

def unitizer(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        func(self, *args, **kwargs)
        self.unitize(inplace=True)
    return wrapper

class UnitVector(Point):
    @unitizer
    def __init__(self, x = 0, y = 0, z = 0):
        super().__init__(x, y, z)

    @unitizer
    def __setitem__(self, i, value: float | int) -> None:
        super().__setitem__(i, value)

    @unitizer
    def __add__(self, obj: Point | float | int | list) -> Point:
        super().__add__(obj)
    
    @unitizer
    def __iadd__(self, obj: Point | float | int | list) -> Point:
        super().__iadd__(obj)

    @unitizer
    def __sub__(self, obj: Point | float | int | list) -> Point:
        super().__sub__(obj)
    
    @unitizer
    def __rsub__(self, obj: Point | float | int | list) -> Point:
        super().__rsub__(obj)

    @unitizer
    def __isub__(self, obj: Point | float | int | list) -> Point:
        super().__isub__(obj)

    @unitizer
    def __mul__(self, obj: Point | list | float | int) -> Point:
        super().__mul__(obj)
    
    @unitizer
    def __imul__(self, obj: Point | list | float | int) -> Point:
        super().__imul__(obj)

    @unitizer
    def __floordiv__(self, obj: Point | list | float | int) -> Point:
        super().__floordiv__(obj)
    
    @unitizer
    def __ifloordiv__(self, obj: Point | list | float | int) -> Point:
        super().__ifloordiv__(obj)

    @unitizer
    def __truediv__(self, obj: Point | list | float | int) -> Point:
        super().__truediv__(obj)
    
    @unitizer
    def __itruediv__(self, obj: Point | list | float | int) -> Point:
        super().__itruediv__(obj)

    @unitizer
    def __matmul__(self, obj: Point | list) -> Point:
        super().__matmul__(obj)
    
    @unitizer
    def __imatmul__(self, obj: Point | list) -> Point:
        super().__imatmul__(obj)

    @staticmethod
    def from_list(l: list) -> Point:
        return UnitVector(l[0], l[1], l[2])