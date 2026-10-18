"""
This file demonstrates stuff that can be done with the Point class
"""

# Load the pysoltrace api from the parent directory ---
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

from pysoltrace import Point

a = Point(6, 0, 8)

# __repr__
print('Examples of printing:')
print(f'a -> {a}')
print(f'a -> {a:.4e}')

# __bool__
zeros = Point()
print('\nExamples of boolean casting:')
print(f'bool(a) -> {bool(a)}')
print(f'bool(zeros) -> {bool(zeros)}')

# __float__
print('\nExample of float casting:')
print(f'float(a) -> {float(a)}')

# __int__
print('\nExample of int casting:')
print(f'int(a) -> {int(a)}')

# __iter__
# __next__
print('\nExample of iterating:')
sum = 0
for v in a: sum += v
print(f'for v in a: sum += v -> {sum}')

# __getitem__
# __setitem__
print('\nExamples of list-ish interface:')
print(f'a[1] -> {a[1]}')
a[1] = -17
print(f'a[1] = -17 -> {a[1]}')

# __len__
print('\nExample of length:')
print(f'len(a) -> {len(a)}')

# __neg__
print('\nExample of negation:')
print(f'-a -> {-a}')

# __abs__
print('\nExample of absolute value:')
print(f'abs(a) -> {abs(a)}')

b = Point(3, 0, 3)
# __add__
print('\nExamples of adding:')
print(f'a + b         -> {a + b}')
print(f'a + 1         -> {a + 1}')
print(f'a + 2.1       -> {a + 2.1}')
print(f'a + [3, 0, 3] -> {a + [3, 0, 3]}')

# __radd__
print('\nExamples of right adding:')
print(f'1 + b         -> {1 + b}')
print(f'2.1 + b       -> {2.1 + b}')
print(f'[6, 0, 8] + b -> {[6, 0, 8] + b}')

# __iadd__
print('\nExamples of inplace adding:')
c = Point()
c += a + b
print(f'c += a + b     -> {c}')
c += 1
print(f'c += 1         -> {c}')
c += 2.1
print(f'c += 2.1       -> {c}')
c += [6, 0, 8]
print(f'c += [6, 0, 8] -> {c}')

# __sub__
print('\nExamples of subtracting:')
print(f'a - b         -> {a - b}')
print(f'a - 1         -> {a - 1}')
print(f'a - 2.1       -> {a - 2.1}')
print(f'a - [3, 0, 3] -> {a - [3, 0, 3]}')

# __rsub__
print('\nExamples of right subtracting:')
print(f'1 - b         -> {1 - b}')
print(f'2.1 - b       -> {2.1 - b}')
print(f'[6, 0, 8] - b -> {[6, 0, 8] - b}')

# __isub__
print('\nExamples of inplace subtracting:')
c -= a + b
print(f'c -= a + b     -> {c}')
c -= 1
print(f'c -= 1         -> {c}')
c -= 2.1
print(f'c -= 2.1       -> {c}')
c -= [6, 0, 8]
print(f'c -= [6, 0, 8] -> {c}')

# __mul__
print('\nExamples of multiplying:')
print(f'a * b         -> {a * b}')
print(f'a * 4         -> {a * 4}')
print(f'a * 2.1       -> {a * 2.1}')
print(f'a * [3, 0, 3] -> {a * [3, 0, 3]}')

# __rmul__
print('\nExamples of right multiplying:')
print(f'4 * b         -> {4 * b}')
print(f'2.1 * b       -> {2.1 * b}')
print(f'[6, 0, 8] * b -> {[6, 0, 8] * b}')

# __imul__
print('\nExamples of inplace multiplying:')
c += a
c *= b
print(f'c *= b         -> {c}')
c = a.copy()
c *= 4
print(f'c *= 4         -> {c}')
c *= 2.1
print(f'c *= 2.1       -> {c}')
c *= [6, 0, 8]
print(f'c *= [6, 0, 8] -> {c}')

# __floordiv__
print('\nExamples of floordiv:')
d = Point(3.2, 5.0, 2.4)
print(f'a // d         -> {a // d}')
print(f'a // 4         -> {a // 4}')
print(f'a // 2.1       -> {a // 2.1}')
print(f'a // [3, 10, 3] -> {a // [3, 10, 3]}')

# __ifloordiv__
print('\nExamples of inplace floordiv:')
c = a.copy()
c //= d
print(f'c //= d         -> {c}')
c = a.copy()
c //= 4
print(f'c //= 4         -> {c}')
c //= 2.1
print(f'c //= 2.1       -> {c}')
c //= [6, 2, 8]
print(f'c //= [6, 2, 8] -> {c}')

# __truediv__
print('\nExamples of truediv:')
print(f'a / d         -> {a / d}')
print(f'a / 4         -> {a / 4}')
print(f'a / 2.1       -> {a / 2.1}')
print(f'a / [3, 10, 3] -> {a / [3, 10, 3]}')

# __itruediv__
print('\nExamples of inplace truediv:')
c = a.copy()
c /= d
print(f'c /= d         -> {c}')
c = a.copy()
c /= 4
print(f'c /= 4         -> {c}')
c /= 2.1
print(f'c /= 2.1       -> {c}')
c /= [6, 2, 8]
print(f'c /= [6, 2, 8] -> {c}')

# __matmul__
# __imatmul__
print('\nExamples of matmul (cross product):')
print(a)
print(b)
print(f'a @ b          -> {a @ b}')
print(f'a.dot(a @ b)   -> {a.dot(a @ b)}')
print(f'a @ [3, 10, 3] -> {a @ [3, 10, 3]}')

print('\nExamples of inplace matmul (cross product):')
c = a.copy()
c @= b
print(f'c @= b         -> {c}')
c = a.copy()
c @= [6, 0, 8]
print(f'c @= [6, 0, 8] -> {c}')

# __eq__
e = Point(6, -17, 8)
print('\nExamples of equality:')
print(f'a == e           -> {a == e}')
print(f'a == e + 1       -> {a == e + 1}')
print(f'a == [6, -17, 8] -> {a == [6, -17, 8]}')
print(f'a == [7, -16, 9] -> {a == [7, -16, 9]}')

# __ne__
print('\nExamples of not equality:')
print(f'a != e           -> {a != e}')
print(f'a != e + 1       -> {a != e + 1}')
print(f'a != [6, -17, 8] -> {a != [6, -17, 8]}')
print(f'a != [7, -16, 9] -> {a != [7, -16, 9]}')

# __lt__
print('\nExamples of less than:')
print(f'a < e           -> {a < e}')
print(f'a < e + 1       -> {a < e + 1}')
print(f'a < e - 1       -> {a < e - 1}')
print(f'a < [6, -17, 8] -> {a < [6, -17, 8]}')
print(f'a < [7, -16, 9] -> {a < [7, -16, 9]}')
print(f'a < [5, -18, 7] -> {a < [5, -18, 7]}')

# __gt__
print('\nExamples of greater than:')
print(f'a > e           -> {a > e}')
print(f'a > e + 1       -> {a > e + 1}')
print(f'a > e - 1       -> {a > e - 1}')
print(f'a > [6, -17, 8] -> {a > [6, -17, 8]}')
print(f'a > [7, -16, 9] -> {a > [7, -16, 9]}')
print(f'a > [5, -18, 7] -> {a > [5, -18, 7]}')

# __le__
print('\nExamples of less than or equal:')
print(f'a <= e           -> {a <= e}')
print(f'a <= e + 1       -> {a <= e + 1}')
print(f'a <= e - 1       -> {a <= e - 1}')
print(f'a <= [6, -17, 8] -> {a <= [6, -17, 8]}')
print(f'a <= [7, -16, 9] -> {a <= [7, -16, 9]}')
print(f'a <= [5, -18, 7] -> {a <= [5, -18, 7]}')

# __ge__
print('\nExamples of greater than or equal:')
print(f'a >= e           -> {a >= e}')
print(f'a >= e + 1       -> {a >= e + 1}')
print(f'a >= e - 1       -> {a >= e - 1}')
print(f'a >= [6, -17, 8] -> {a >= [6, -17, 8]}')
print(f'a >= [7, -16, 9] -> {a >= [7, -16, 9]}')
print(f'a >= [5, -18, 7] -> {a >= [5, -18, 7]}')

# reduce
print('\nExamples of Point.reduce:')
print(f'a.reduce() -> {a.reduce()} from {a}')
print(f'b.reduce() -> {b.reduce()} from {b}')
print(f'c.reduce() -> {c.reduce()} from {c}')
print(f'd.reduce() -> {d.reduce()} from {d}')
print(f'e.reduce() -> {e.reduce()} from {e}')

# radius
print('\nExamples of Point.radius:')
print(f'a.radius()     -> {a.radius()} from {a}')
print(f'b.radius()     -> {b.radius()} from {b}')
print(f'c.radius()     -> {c.radius()} from {c}')
print(f'd.radius()     -> {d.radius()} from {d}')
print(f'e.radius()     -> {e.radius()} from {e}')
print(f'zeros.radius() -> {zeros.radius()} from {zeros}')

# unitize
print('\nExamples of Point.unitize:')
print(f'a.unitize()     -> {a.unitize()} from {a}')
print(f'b.unitize()     -> {b.unitize()} from {b}')
print(f'c.unitize()     -> {c.unitize()} from {c}')
print(f'd.unitize()     -> {d.unitize()} from {d}')
print(f'e.unitize()     -> {e.unitize()} from {e}')
print(f'zeros.unitize() -> {zeros.unitize()} from {zeros}')

# as_list
print('\nExamples of Point.as_list:')
print(f'a.as_list()     -> {a.as_list()} from {a}')
print(f'b.as_list()     -> {b.as_list()} from {b}')
print(f'c.as_list()     -> {c.as_list()} from {c}')
print(f'd.as_list()     -> {d.as_list()} from {d}')
print(f'e.as_list()     -> {e.as_list()} from {e}')
print(f'zeros.as_list() -> {zeros.as_list()} from {zeros}')

# from_list
print('\nExamples of Point.from_list:')
print(f'Point.from_list([6, 0, 8]) -> {Point.from_list([6, 0, 8])}')
print(f'Point.from_list([3, 0, 3]) -> {Point.from_list([3, 0, 3])}')