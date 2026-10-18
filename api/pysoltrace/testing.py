import ctypes, unittest, os, sys
from pathlib import Path

# pysoltrace_path = Path(__file__).parent / 'pysoltrace'
# sys.path.insert(1, str(pysoltrace_path))

from math import sin, cos, pi, sqrt
import orjson
import numpy as np

# TODO: make a script to blow up the build dir, run the wheel build, make the build dir, make .venv and pip install wheel, run this file

from pysoltrace import dot_h, found_in, \
                       api, \
                       PySolTrace, \
                       Point, \
                       STAPIv2Exception, \
                       soltrace_constants as _STC, \
                       soltrace_structs as _STS

locate = lambda rel: (Path(__file__).parent / rel).resolve()

class PointTests(unittest.TestCase):
    def setUp(self):
        self.a = Point(6, 0, 8)
        self.b = Point(3, 0, 3)
        self.d = Point(3.2, 5.0, 2.4)
        self.e = Point(6, 0, 8)
        self.zeros = Point()
        return super().setUp()
    
    def test_repr(self):
        self.assertEqual(f'{self.a}', '[6.00, 0.00, 8.00]')
        self.assertEqual(f'{self.a:.4e}', '[6.0000e+00, 0.0000e+00, 8.0000e+00]')

    def test_bool(self):
        self.assertEqual(bool(self.a), True)
        self.assertEqual(bool(self.zeros), False)

    def test_float(self):
        self.assertEqual(float(self.a), 10.0)

    def test_int(self):
        self.assertEqual(int(self.a), 10)

    def test_iter(self):
        sum = 0
        for v in self.a: sum += v
        self.assertEqual(sum, 14)

    def test_get_set_item(self):
        self.assertEqual(self.a[1], 0)
        self.a[1] = -17
        self.assertEqual(self.a[1], -17)

    def test_len(self):
        self.assertEqual(len(self.a), 3)

    def test_negation(self):
        self.assertEqual(-self.a, Point(-6, 0, -8))

    def test_abs(self):
        self.a[1] = -17
        self.assertEqual(abs(self.a), 19.72308292331602)

    def test_add(self):
        self.assertEqual(self.a + self.b,              Point(9, 0, 11))
        self.assertEqual(self.a + 1,                   Point(7, 1, 9))
        self.assertEqual(self.a + 2.1,                 Point(8.1, 2.1, 10.1))
        self.assertEqual(self.a + [3, 0, 3],           Point(9, 0, 11))
        self.assertEqual(self.a + np.array([3, 0, 3]), Point(9, 0, 11))
        # self.assertEqual(self.a + 'not implemented',   NotImplemented)

    def test_radd(self):
        self.assertEqual(1 + self.b,                   Point(4, 1, 4))
        self.assertEqual(2.1 + self.b,                 Point(5.1, 2.1, 5.1))
        self.assertEqual([6, 0, 8] + self.b,           Point(9, 0, 11))
        # right add with np.array returns np.array
        temp = np.array([6, 0, 8]) + self.b
        self.assertAlmostEqual(Point(9, 0, 11) - temp, 0)

    def test_iadd(self):
        c = Point()
        c += self.a + self.b
        self.assertEqual(c, Point(9, 0, 11))
        c += 1
        self.assertEqual(c, Point(10, 1, 12))
        c += 2.1
        self.assertEqual(c, Point(12.1, 3.1, 14.1))
        c += [2, 2, 2]
        self.assertEqual(c, Point(14.1, 5.1, 16.1))
        c += np.array([-2, -2, -2])
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(c - Point(12.1, 3.1, 14.1), 0)

    def test_sub(self):
        self.assertEqual(self.a - self.b,              Point(3, 0, 5))
        self.assertEqual(self.a - 1,                   Point(5, -1, 7))
        self.assertEqual(self.a - 2.1,                 Point(3.9, -2.1, 5.9))
        self.assertEqual(self.a - [3, 0, 3],           Point(3, 0, 5))
        self.assertEqual(self.a - np.array([3, 0, 3]), Point(3, 0, 5))

    def test_rsub(self):
        self.assertEqual(1 - self.b, Point(-2, 1, -2))
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(2.1 - self.b - Point(-.9, 2.1, -.9), 0)
        self.assertEqual([6, 0, 8] - self.b,           Point(3, 0, 5))
        # right subtraction with np.array returns np.array
        temp = np.array([6, 0, 8]) - self.b
        self.assertAlmostEqual(Point(3, 0, 5) - temp, 0)

    def test_isub(self):
        c = Point()
        c -= self.a + self.b
        self.assertEqual(c, Point(-9, 0, -11))
        c -= 1
        self.assertEqual(c, Point(-10, -1, -12))
        c -= 2.1
        self.assertEqual(c, Point(-12.1, -3.1, -14.1))
        c -= [2, 2, 2]
        self.assertEqual(c, Point(-14.1, -5.1, -16.1))
        c -= np.array([-2, -2, -2])
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(c - Point(-12.1, -3.1, -14.1), 0)

    def test_mul(self):
        self.assertEqual(self.a * self.b, Point(18, 0, 24))
        self.assertEqual(self.a * 4,      Point(24, 0, 32))
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(self.a * 2.1 - Point(12.6, 0, 16.8), 0)
        self.assertEqual(self.a * [3, 0, 3],           Point(18, 0, 24))
        self.assertEqual(self.a * np.array([3, 0, 3]), Point(18, 0, 24))

    def test_rmul(self):
        self.assertEqual(4 * self.b, Point(12, 0, 12))
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(2.1 * self.b - Point(6.3, 0, 6.3), 0)
        self.assertEqual([6, 0, 8] * self.b,           Point(18, 0, 24))
        # right subtraction with np.array returns np.array
        temp = np.array([6, 0, 8]) * self.b
        self.assertAlmostEqual(Point(18, 0, 24) - temp, 0)

    def test_imul(self):
        c = self.a.copy()
        c *= self.b
        self.assertEqual(c, Point(18, 0, 24))
        c *= 4
        self.assertEqual(c, Point(72, 0, 96))
        c *= 2.1
        # floating point error makes assertEqual fail
        self.assertAlmostEqual(c - Point(151.2, 0, 201.6), 0)
        c *= [2, 2, 2]
        self.assertAlmostEqual(c - Point(302.4, 0, 403.2), 0)
        c *= np.array([2, 2, 2])
        self.assertAlmostEqual(c - Point(604.8, 0, 806.4), 0)

    def test_floordiv(self):
        self.assertEqual(self.a // self.d,              Point(1, 0, 3))
        self.assertEqual(self.a // 4,                   Point(1, 0, 2))
        self.assertEqual(self.a // 2.1,                 Point(2, 0, 3))
        self.assertEqual(self.a // [3, 1, 3],           Point(2, 0, 2))
        self.assertEqual(self.a // np.array([3, 1, 3]), Point(2, 0, 2))

    def test_ifloordiv(self):
        c = self.a.copy()
        c //= self.d
        self.assertEqual(c, Point(1, 0, 3))
        c = self.a.copy()
        c //= 4
        self.assertEqual(c, Point(1, 0, 2))
        c = self.a.copy()
        c //= 2.1
        self.assertEqual(c, Point(2, 0, 3))
        c = self.a.copy()
        c //= [3, 1, 3]
        self.assertEqual(c, Point(2, 0, 2))
        c = self.a.copy()
        c //= np.array([3, 1, 3])
        self.assertEqual(c, Point(2, 0, 2))

    def test_truediv(self):
        self.assertEqual(self.a / self.d,              Point(6 / 3.2, 0, 8 / 2.4))
        self.assertEqual(self.a / 4,                   Point(3 / 2, 0, 2))
        self.assertEqual(self.a / 2.1,                 Point(6 / 2.1, 0, 8 / 2.1))
        self.assertEqual(self.a / [3, 1, 3],           Point(2, 0, 8 / 3))
        self.assertEqual(self.a / np.array([3, 1, 3]), Point(2, 0, 8 / 3))

    def test_itruediv(self):
        c = self.a.copy()
        c /= self.d
        self.assertEqual(c, Point(6 / 3.2, 0, 8 / 2.4))
        c = self.a.copy()
        c /= 4
        self.assertEqual(c, Point(3 / 2, 0, 2))
        c = self.a.copy()
        c /= 2.1
        self.assertEqual(c, Point(6 / 2.1, 0, 8 / 2.1))
        c = self.a.copy()
        c /= [3, 1, 3]
        self.assertEqual(c, Point(2, 0, 8 / 3))
        c = self.a.copy()
        c /= np.array([3, 1, 3])
        self.assertEqual(c, Point(2, 0, 8 / 3))

    def test_cross(self):
        self.a[1] = -17
        self.assertEqual(self.a @ self.b,               Point(-51, 6, 51))
        self.assertEqual(self.a.dot(self.a @ self.b),   0)
        self.assertEqual(self.a @ [3, 10, 3],           Point(-131, 6, 111))
        self.assertEqual(self.a @ np.array([3, 10, 3]), Point(-131, 6, 111))

    def test_icross(self):
        self.a[1] = -17
        c = self.a.copy()
        c @= self.b
        self.assertEqual(c, Point(-51, 6, 51))
        c = self.a.copy()
        c @= [3, 10, 3]
        self.assertEqual(c, Point(-131, 6, 111))
        c = self.a.copy()
        c @= np.array([3, 10, 3])
        self.assertEqual(c, Point(-131, 6, 111))

    def test_eq(self):
        self.assertEqual(self.a,    self.e)
        self.assertNotEqual(self.a, self.e + 1)
        self.assertEqual(self.a,    [6, 0, 8])
        self.assertNotEqual(self.a, [6, 1, 8])
        self.assertEqual(self.a,    np.array([6, 0, 8]))
        self.assertNotEqual(self.a, np.array([6, 1, 8]))

    def test_neq(self):
        self.assertEqual(self.a != self.e,              False)
        self.assertEqual(self.a != self.e + 1,          True)
        self.assertEqual(self.a != [6, 0, 8],           False)
        self.assertEqual(self.a != [6, 1, 8],           True)
        self.assertEqual(self.a != np.array([6, 0, 8]), False)
        self.assertEqual(self.a != np.array([6, 1, 8]), True)

    def test_lt(self):
        self.assertEqual(self.a < self.e,               False)
        self.assertEqual(self.a < self.e + 1,           True)
        self.assertEqual(self.a < self.e - 1,           False)
        self.assertEqual(self.a < [6, 0, 8],            False)
        self.assertEqual(self.a < [6, 1, 8],            True)
        self.assertEqual(self.a < [6, -1, 8],           True)
        self.assertEqual(self.a < [5, -1, 7],           False)
        self.assertEqual(self.a < np.array([6, 0, 8]),  False)
        self.assertEqual(self.a < np.array([6, 1, 8]),  True)
        self.assertEqual(self.a < np.array([6, -1, 8]), True)
        self.assertEqual(self.a < np.array([5, -1, 7]), False)

    def test_gt(self):
        self.assertEqual(self.a > self.e,               False)
        self.assertEqual(self.a > self.e + 1,           False)
        self.assertEqual(self.a > self.e - 1,           True)
        self.assertEqual(self.a > [6, 0, 8],            False)
        self.assertEqual(self.a > [6, 1, 8],            False)
        self.assertEqual(self.a > [6, -1, 8],           False)
        self.assertEqual(self.a > [5, -1, 7],           True)
        self.assertEqual(self.a > np.array([6, 0, 8]),  False)
        self.assertEqual(self.a > np.array([6, 1, 8]),  False)
        self.assertEqual(self.a > np.array([6, -1, 8]), False)
        self.assertEqual(self.a > np.array([5, -1, 7]), True)

    def test_lte(self):
        self.assertEqual(self.a <= self.e,               True)
        self.assertEqual(self.a <= self.e + 1,           True)
        self.assertEqual(self.a <= self.e - 1,           False)
        self.assertEqual(self.a <= [6, 0, 8],            True)
        self.assertEqual(self.a <= [6, 1, 8],            True)
        self.assertEqual(self.a <= [6, -1, 8],           True)
        self.assertEqual(self.a <= [5, -1, 7],           False)
        self.assertEqual(self.a <= np.array([6, 0, 8]),  True)
        self.assertEqual(self.a <= np.array([6, 1, 8]),  True)
        self.assertEqual(self.a <= np.array([6, -1, 8]), True)
        self.assertEqual(self.a <= np.array([5, -1, 7]), False)

    def test_gte(self):
        self.assertEqual(self.a >= self.e,               True)
        self.assertEqual(self.a >= self.e + 1,           False)
        self.assertEqual(self.a >= self.e - 1,           True)
        self.assertEqual(self.a >= [6, 0, 8],            True)
        self.assertEqual(self.a >= [6, 1, 8],            False)
        self.assertEqual(self.a >= [6, -1, 8],           False)
        self.assertEqual(self.a >= [5, -1, 7],           True)
        self.assertEqual(self.a >= np.array([6, 0, 8]),  True)
        self.assertEqual(self.a >= np.array([6, 1, 8]),  False)
        self.assertEqual(self.a >= np.array([6, -1, 8]), False)
        self.assertEqual(self.a >= np.array([5, -1, 7]), True)

    def test_reduce(self):
        self.assertEqual(self.a.reduce(),     14)
        self.assertEqual(self.b.reduce(),     6)
        self.assertEqual(self.d.reduce(),     10.6)
        self.assertEqual(self.e.reduce(),     14)
        self.assertEqual(self.zeros.reduce(), 0)

    def test_radius(self):
        self.assertEqual(self.a.radius(),     10)
        self.assertEqual(self.b.radius(),     4.242640687119285)
        self.assertEqual(self.d.radius(),     6.4031242374328485)
        self.assertEqual(self.e.radius(),     10)
        self.assertEqual(self.zeros.radius(), 0)

    def test_unitize(self):
        self.assertEqual(self.a.unitize(),     Point(.6, 0, .8))
        self.assertEqual(self.b.unitize(),     Point(3 / sqrt(18), 0, 3 / sqrt(18)))
        self.assertEqual(self.d.unitize(),     Point(self.d.x / self.d.radius(),
                                                     self.d.y / self.d.radius(),
                                                     self.d.z / self.d.radius()))
        self.assertEqual(self.e.unitize(),     Point(.6, 0, .8))
        self.assertEqual(self.zeros.unitize(), Point())

    def test_as_list(self):
        self.assertEqual(self.a.as_list(),     [6, 0, 8])
        self.assertEqual(self.b.as_list(),     [3, 0, 3])
        self.assertEqual(self.d.as_list(),     [3.2, 5, 2.4])
        self.assertEqual(self.e.as_list(),     [6, 0, 8])
        self.assertEqual(self.zeros.as_list(), [0, 0, 0])

    def test_from_list(self):
        self.assertEqual(Point.from_list([6, 0, 8]),     Point(6, 0, 8))
        self.assertEqual(Point.from_list([3, 0, 3]),     Point(3, 0, 3))
        self.assertEqual(Point.from_list([3.2, 5, 2.4]), Point(3.2, 5, 2.4))
        self.assertEqual(Point.from_list(np.array([6, 0, 8])),     Point(6, 0, 8))
        self.assertEqual(Point.from_list(np.array([3, 0, 3])),     Point(3, 0, 3))
        self.assertEqual(Point.from_list(np.array([3.2, 5, 2.4])), Point(3.2, 5, 2.4))

class STAPIv2TestCase(unittest.TestCase):
    def setUp(self):
        self.stapi = api(testing=True)
        return super().setUp()

    def assertStructEqual(self, value, check):
        for field in value._fields_:
            self.assertEqual(getattr(value, field[0]),
                             getattr(check, field[0]))

class ConstantsTests(STAPIv2TestCase):
    def test_error_code_msg(self):
        self.assertIn(dot_h.st_return_code.FAILURE, _STC.ST_RETURN_CODE_ERROR_MSG)
        self.assertEqual(len(_STC.ST_RETURN_CODE_ERROR_MSG[dot_h.st_return_code.FAILURE]), 0)

        for i in range(dot_h.st_return_code.CANCEL, 
                       dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE):
            self.assertIn(i, _STC.ST_RETURN_CODE_ERROR_MSG)
            self.assertGreater(len(_STC.ST_RETURN_CODE_ERROR_MSG[i]), 0)

        _, _, check_return_code = self.stapi.sneak()
        for i in range(dot_h.st_return_code.CANCEL, 
                       dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE):
            with self.assertRaises(STAPIv2Exception) as ex:
                check_return_code(i)
            self.assertEqual(ex.exception.code, i)

        for i in range(dot_h.st_return_code.RETURN_COUNT, 
                       dot_h.st_return_code.RETURN_COUNT + 2):
            with self.assertRaises(STAPIv2Exception) as ex:
                check_return_code(i)
            self.assertEqual(ex.exception.code, i)

    def test_warning_code_msg(self):
        for i in range(dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE, 
                       dot_h.st_return_code.RETURN_COUNT):
            self.assertIn(i, _STC.ST_RETURN_CODE_WARNING_MSG)
            self.assertGreater(len(_STC.ST_RETURN_CODE_WARNING_MSG[i]), 0)

        _, _, check_return_code = self.stapi.sneak()
        for i in range(dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE, 
                       dot_h.st_return_code.RETURN_COUNT):
            with self.assertWarns(UserWarning):
                check_return_code(i)

    def test_args_simulation_parameters(self):
        params = _STS.simulation_parameters(35.962278, -106.5122622, 1, 100, .1, True, True, False)
        c_params = dot_h.args_simulation_parameters(1, 100, .1, 35.962278, -106.5122622, True, True, False)
        self.assertStructEqual(params.ctype, c_params)

    def test_args_optical_properties_face(self):
        args = (.5, .5, 5, 5, b'g')
        face = _STS.optical_properties_face(*args)
        c_face = dot_h.args_optical_properties_face(*args)
        self.assertStructEqual(face.ctype, c_face)

    def test_args_optical_properties_set(self):
        args = (b'test', 1.1, 1.1, 2)
        opt_set = _STS.optical_properties_set(*args)
        c_opt_set = dot_h.args_optical_properties_set(*args)
        self.assertStructEqual(opt_set.ctype, c_opt_set)

    def test_args_element(self):
        args_list = ([2, 2, 2], [2, 2, 2], 2, False, True)
        args_no_list = (2, 2, 2, 2, 2, 2, 2, False, True)
        el = _STS.element(*args_list)
        c_el = dot_h.args_element(*args_no_list)
        self.assertStructEqual(el.ctype, c_el)

    def test_args_geometry(self):
        args_no_array = ([2], [2, 2], b'c', b'p')
        args_array = (_STS.make_c_double_8([2]), _STS.make_c_double_8([2, 2]), b'c', b'p')
        geom_fields = _STS.geometry(*args_no_array).ctype._fields_
        c_geom_fields = dot_h.args_geometry(*args_array)._fields_
        # check individually be point to different array locations
        # aperture params
        self.assertEqual(geom_fields[0][1][0], c_geom_fields[0][1][0])
        # surface params
        self.assertEqual(geom_fields[1][1][0], c_geom_fields[1][1][0])
        self.assertEqual(geom_fields[1][1][1], c_geom_fields[1][1][1])
        # aperture
        self.assertEqual(geom_fields[2][1], c_geom_fields[2][1])
        # surf
        self.assertEqual(geom_fields[3][1], c_geom_fields[3][1])

    def test_args_sun(self):
        args = (0, 2, 2, 2, .5, b'b')
        sun = _STS.sun(*args)
        c_sun = dot_h.args_sun(*args)
        self.assertStructEqual(sun.ctype, c_sun)

    def test_args_sun_location(self):
        args = (40.0, -105.0, -7.0)
        loc = _STS.sun_location(*args)
        c_loc = dot_h.args_sun_location(*args)
        self.assertStructEqual(loc.ctype, c_loc)

    def test_args_sun_datetime(self):
        args = (2025, 6, 20)
        dt = _STS.sun_datetime(*args)
        c_dt = dot_h.args_sun_datetime(*args)
        self.assertStructEqual(dt.ctype, c_dt)

class ParametersTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        self.params = _STS.simulation_parameters(35.962278, -106.5122622, 1, 100, .1, True, True, False)

    def test_set_parameters(self):
        self.stapi.parameters.set(self.params)

    def test_get_parameters(self):
        self.stapi.parameters.set(self.params)

        rt_params = self.stapi.parameters.get()
        self.assertDictEqual(rt_params, self.params.ctype.value)

    def test_set_ray_parameters(self):
        self.stapi.parameters.rays(10, 1000)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['number_of_rays'], 10)
        self.assertEqual(rt_params['max_number_of_rays'], 1000)

    def test_set_power_tower_parameter(self):
        self.stapi.parameters.power_tower(True)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['as_power_tower'], True)

    def test_set_error_parameters(self):
        self.stapi.parameters.errors(True, True)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['include_sun_shape_errors'], True)
        self.assertEqual(rt_params['include_optical_errors'], True)

    def test_set_location_parameters(self):
        self.stapi.parameters.location(35.962278, -106.5122622)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['latitude'], 35.962278)
        self.assertEqual(rt_params['longitude'], -106.5122622)

    def test_set_tolerance_parameter(self):
        self.stapi.parameters.tolerance(.1)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['tolerance'], .1)

    def test_legacy_params(self):
        self.stapi.legacy.sim_params(10, 1000, True)

        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['number_of_rays'], 10)
        self.assertEqual(rt_params['max_number_of_rays'], 1000)
        self.assertEqual(rt_params['as_power_tower'], True)

class DataJSONTests(STAPIv2TestCase):
    def test_load_json_str(self):
        f = open(locate('sample.json'), mode='rb')
        _json = f.read()
        f.close()
        self.stapi.data.json.load(_json)
        self.assertEqual(self.stapi.data.element.num(), 126)

        with self.assertRaises(STAPIv2Exception) as ex:
            f = open(locate('errors.json'), mode='rb')
            _bad_json = f.read()
            f.close()
            self.stapi.data.json.load(_bad_json)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.EXCEPTION)

    def test_load_json_dict(self):
        f = open(locate('sample.json'), mode='rb')
        _json = orjson.loads(f.read())
        f.close()
        self.stapi.data.json.load(_json)
        self.assertEqual(self.stapi.data.element.num(), 126)

        with self.assertRaises(STAPIv2Exception) as ex:
            f = open(locate('errors.json'), mode='rb')
            _json = orjson.loads(f.read())
            f.close()
            self.stapi.data.json.load(_json)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.EXCEPTION)

    def test_load_json_file(self):
        good_path = locate('sample.json')
        self.stapi.data.json.load(good_path)
        self.assertEqual(self.stapi.data.element.num(), 126)

        bad_path = locate('errors.json')
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.json.load(bad_path)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.EXCEPTION)

class OpticalPropertiesTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        self.opt_set = _STS.optical_properties_set(b'test', 1.1, 1.1, 2)
        self.front   = _STS.optical_properties_face(.5, .5, 5, 5, b'g')
        self.back    = _STS.optical_properties_face(.25, .25, 2, 2, b'g')

    def test_add_optic(self):
        optical_id = self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.assertEqual(optical_id, 0)

        with self.assertRaises(STAPIv2Exception) as ex:
            bad_front = _STS.optical_properties_face(.5, .5, 5, 5, b'z')
            self.stapi.data.optic.add(self.opt_set, bad_front, self.back)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

        with self.assertRaises(STAPIv2Exception) as ex:
            bad_back = _STS.optical_properties_face(.25, .25, 2, 2, b'z')
            self.stapi.data.optic.add(self.opt_set, self.front, bad_back)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

    def test_get_optic(self):
        optical_id = self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.assertEqual(optical_id, 0)

        _opt_set, _front, _back = self.stapi.data.optic.get(optical_id)

        self.assertDictEqual(_opt_set, self.opt_set.ctype.value)
        self.assertDictEqual(_front, self.front.ctype.value)
        self.assertDictEqual(_back, self.back.ctype.value)

        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.optic.get(optical_id + 1)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.DATA_VALUE_NOT_FOUND)

    def test_delete_optic(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)

        self.assertEqual(self.stapi.data.optic.num(), 2)

        with self.assertWarns(UserWarning):
            self.stapi.data.optic.delete(2)

        self.stapi.data.optic.delete(0)
        self.assertEqual(self.stapi.data.optic.num(), 1)

        self.stapi.data.optic.delete(1)
        self.assertEqual(self.stapi.data.optic.num(), 0)

    def test_clear_optics(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.assertEqual(self.stapi.data.optic.num(), 2)

        self.stapi.data.optic.clear()
        self.assertEqual(self.stapi.data.optic.num(), 0)

class ElementTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        # set up dummy optical set
        opt_set = _STS.optical_properties_set(b'dummy', 1.1, 1.1, 0)
        front   = _STS.optical_properties_face(.5, .5, 5, 5, b'g')
        back    = _STS.optical_properties_face(.25, .25, 2, 2, b'g')
        
        self.stapi.data.optic.add(opt_set, front, back)
        self.el_args = _STS.element([2, 2, 2], [2, 2, 2], 2, False, True)
        self.opt_id = 0
        self.geom_args = _STS.geometry([2], [2, 2], b'c', b'p')

    def test_add_element(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # optical property not set
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.add(self.el_args, self.opt_id + 1, self.geom_args)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.DATA_VALUE_NOT_FOUND)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # bad aperture char
        with self.assertRaises(STAPIv2Exception) as ex:
            bad_ap = _STS.geometry([2], [2, 2], b'z', b'p')
            # self.el_args.ap = b'z'
            self.stapi.data.element.add(self.el_args, self.opt_id, bad_ap)
            # self.el_args.ap = b'c'
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # bad surface char
        with self.assertRaises(STAPIv2Exception) as ex:
            bad_surf = _STS.geometry([2], [2, 2], b'c', b'z')
            # self.el_args.surf = b'z'
            self.stapi.data.element.add(self.el_args, self.opt_id, bad_surf)
            # self.el_args.surf = b'p'
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # bad aperture params
        with self.assertRaises(STAPIv2Exception) as ex:
            bad_ap_params = _STS.geometry([-2], [2, 2], b'c', b'p')
            self.stapi.data.element.add(self.el_args, self.opt_id, bad_ap_params)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # bad surface params
        with self.assertRaises(STAPIv2Exception) as ex:
            bad_surf_params = _STS.geometry([2], [2, float('nan')], b'c', b'p')
            self.stapi.data.element.add(self.el_args, self.opt_id, bad_surf_params)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)
        self.assertEqual(self.stapi.data.element.num(), 1)

    def test_get_element(self):
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(id, 1)

        args, optic_id, geom = self.stapi.data.element.get(id)

        self.assertDictEqual(args, self.el_args.ctype.value)
        self.assertEqual(optic_id, self.opt_id)
        self.assertEqual(geom['a_params'][0], self.geom_args.a_params[0])
        self.assertEqual(geom['s_params'][0], 1 / (2 * self.geom_args.s_params[0]))
        self.assertEqual(geom['s_params'][1], 1 / (2 * self.geom_args.s_params[1]))

        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.get(id + 1)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.DATA_VALUE_NOT_FOUND)

    def test_delete_element(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 2)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.delete(0)
        
        with self.assertWarns(UserWarning):
            self.stapi.data.element.delete(3)

        self.stapi.data.element.delete(1)
        self.assertEqual(self.stapi.data.element.num(), 1)
        self.stapi.data.element.delete(2)
        self.assertEqual(self.stapi.data.element.num(), 0)

    def test_clear_elements(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 2)

        self.stapi.data.element.clear()
        self.assertEqual(self.stapi.data.element.num(), 0)

    def test_element_enabled(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.enabled(0, True)

        self.stapi.data.element.enabled(1, True)
        
    def test_element_virtual(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.virtual(0, False)

        self.stapi.data.element.virtual(1, False)

    def test_element_xyz(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.xyz(0, [3, 3, 3])

        with self.assertRaises(AssertionError):
            self.stapi.data.element.xyz(0, [3, 3])

        with self.assertRaises(AssertionError):
            self.stapi.data.element.xyz(0, [3, 3, 3, 3])

        self.stapi.data.element.xyz(1, [3, 3, 3])

    def test_element_aim(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.aim(0, [3, 3, 3])

        with self.assertRaises(AssertionError):
            self.stapi.data.element.aim(0, [3, 3])

        with self.assertRaises(AssertionError):
            self.stapi.data.element.aim(0, [3, 3, 3, 3])

        self.stapi.data.element.aim(1, [3, 3, 3])

    def test_element_zrot(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        with self.assertWarns(UserWarning):
            self.stapi.data.element.zrot(0, 3)

        self.stapi.data.element.zrot(1, 3)

    def test_element_aperture(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)
        new_ap = 'r'
        new_params = [3, 3]

        with self.assertWarns(UserWarning):
            self.stapi.data.element.aperture(0, new_ap, new_params)

        # bad aperture char
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.aperture(1, 'z', new_params)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

        # bad aperture params
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.aperture(1, new_ap, [-1, 2])
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

        self.stapi.data.element.aperture(1, new_ap, new_params)

    def test_element_surface(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)
        new_surf = 's'
        new_params = [3]

        with self.assertWarns(UserWarning):
            self.stapi.data.element.surface(0, new_surf, new_params)

        # bad surface char
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.surface(1, 'z', new_params)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

        # bad surface params
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.surface(1, new_surf, [0])
        self.assertEqual(ex.exception.code, dot_h.st_return_code.INVALID_ARGUMENTS)

        self.stapi.data.element.surface(1, new_surf, new_params)

    def test_element_optic(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # optical property not set
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.element.optic(1, self.opt_id + 1)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.DATA_VALUE_NOT_FOUND)

        self.stapi.data.element.optic(1, self.opt_id)

    def test_element_group(self):
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        self.assertEqual(self.stapi.data.element.num(), 1)

        # ignored groups number
        with self.assertWarns(UserWarning):
            self.stapi.data.element.group(1, -2)

        self.stapi.data.element.group(1, 1)

class SunTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        self.good_angles      = [0, 1, 2]
        self.good_intensities = [0, 1, 2]
        self.bad_intensities  = [0, -1]
        self.args_sun = _STS.sun(3, 608, 303, 1000, 0, b'd')

        self.calc = dot_h.SolarPositionCalculationMethod.SPA
        self.loc = _STS.sun_location(40.0, -105.0, -7.0)
        self.dt  = _STS.sun_datetime(2025, 6, 20)
        self.az = 178.61128380
        self.el = 73.439035265
        self.zen = 90 - self.el
        self.sun_vector = Point(0.006908, -0.2849516, 0.9585169)

    def test_add_sun(self):
        # test bad intensities
        with self.assertRaises(AssertionError) as ex:
            self.stapi.data.sun.add(self.args_sun, self.good_angles, self.bad_intensities)

        # test bad shape
        self.args_sun.npoints = 0
        self.args_sun.shape = b'z'
        with self.assertWarns(UserWarning):
            self.stapi.data.sun.add(self.args_sun, [], [])
        
        # test bad value
        self.args_sun.shape = b'b'
        self.args_sun.sigma_halfwidth_csr = 1
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.sun.add(self.args_sun, [], [])
        self.assertEqual(ex.exception.code, dot_h.st_return_code.EXCEPTION)

        self.args_sun.npoints = 3
        self.stapi.data.sun.add(self.args_sun, self.good_angles, self.good_intensities)

    def test_get_sun(self):
        sun_args = _STS.sun(0, 608, 303, 1000, .5, b'b')
        self.stapi.data.sun.add(sun_args)

        rt_sun_args, rt_angle, rt_intensity = self.stapi.data.sun.get()
        self.assertDictEqual(rt_sun_args, sun_args.ctype.value)
    
    def test_get_sun_userdata(self):
        self.stapi.data.sun.add(self.args_sun, self.good_angles, self.good_intensities)

        rt_sun_args, rt_angle, rt_intensity = self.stapi.data.sun.get(3)
        self.assertDictEqual(rt_sun_args, self.args_sun.ctype.value)
        self.assertEqual(rt_angle, self.good_angles)
        self.assertEqual(rt_intensity, self.good_intensities)
        # TODO: test userdata

    def test_sun_shape(self):
        # test bad shape
        with self.assertWarns(UserWarning):
            self.stapi.data.sun.shape('z', 4.65)
        
        # test bad value
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.data.sun.shape('g', -4.65)
        self.assertEqual(ex.exception.code, dot_h.st_return_code.EXCEPTION)

        self.stapi.data.sun.shape('g', 5)

    def test_sun_xyz(self):
        self.stapi.data.sun.xyz(1, 1, 1)

    def test_sun_userdata(self):
        # test bad intensities
        with self.assertRaises(AssertionError) as ex:
            self.stapi.data.sun.userdata(3, self.good_angles, self.bad_intensities)

        self.stapi.data.sun.userdata(3, self.good_angles, self.good_intensities)

    def test_sun_az_zen(self):
        az, zen = self.stapi.data.sun.az_zen(self.calc, self.loc, self.dt)
        self.assertAlmostEqual(az, self.az, 7)
        self.assertAlmostEqual(zen, self.zen, 7)

    def test_sun_az_el(self):
        az, el = self.stapi.data.sun.az_el(self.calc, self.loc, self.dt)
        self.assertAlmostEqual(az, self.az, 7)
        self.assertAlmostEqual(el, self.el, 7)

    def test_sun_vector(self):
        v = self.stapi.data.sun.vector(self.calc, self.loc, self.dt)
        self.assertAlmostEqual((v - self.sun_vector).radius(), 0, 6)

class RunnerTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        self.stapi.data.json.load(locate('sample.json'))

    def test_set_up_native(self):
        # if this is called after sim_setup call, leads to OS error 
        # bc I think dll handle gets cleaned up before call completes
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.runner.setup(dot_h.st_runner_type_t.NATIVE, 8, [608, 303])
        self.assertEqual(ex.exception.code, dot_h.st_return_code.RUNNER_NUMBER_THREADS_SEEDS_MISMATCH_FAILURE)

        self.stapi.runner.setup(dot_h.st_runner_type_t.NATIVE)

    def test_set_up_embree(self):
        if self.stapi.runner.is_installed(dot_h.st_runner_type_t.EMBREE):
            self.stapi.runner.setup(dot_h.st_runner_type_t.EMBREE)
        else:
            with self.assertWarns(UserWarning):
                self.stapi.runner.setup(dot_h.st_runner_type_t.EMBREE)
        
    def test_set_up_optix(self):
        if self.stapi.runner.is_installed(dot_h.st_runner_type_t.OPTIX):
            self.stapi.runner.setup(dot_h.st_runner_type_t.OPTIX)

            with self.assertWarns(UserWarning):
                self.stapi.runner.setup(dot_h.st_runner_type_t.OPTIX, 1)
        else:
            with self.assertWarns(UserWarning):
                self.stapi.runner.setup(dot_h.st_runner_type_t.OPTIX)

    def runner_run(self, runner_type):
        if not self.stapi.runner.is_installed(runner_type):
            runner_name = dot_h.st_runner_type_t(runner_type).name.title()
            self.skipTest(f'{runner_name} runner not installed, skip testing running {runner_name}Runner.')
        
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.runner.run()
        self.assertEqual(ex.exception.code, dot_h.st_return_code.RUNNER_NOT_FOUND)

        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(runner_type)
        self.stapi.runner.run()

    def test_run_native(self): self.runner_run(dot_h.st_runner_type_t.NATIVE)
    def test_run_embree(self): self.runner_run(dot_h.st_runner_type_t.EMBREE)
    def test_run_optix(self):  self.runner_run(dot_h.st_runner_type_t.OPTIX)

    def runner_report(self, runner_type):
        if not self.stapi.runner.is_installed(runner_type):
            runner_name = dot_h.st_runner_type_t(runner_type).name.title()
            self.skipTest(f'{runner_name} runner not installed, skip testing running {runner_name}Runner.')
        
        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.runner.report()
        self.assertEqual(ex.exception.code, dot_h.st_return_code.RUNNER_NOT_FOUND)

        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(runner_type)

        with self.assertRaises(STAPIv2Exception) as ex:
            self.stapi.runner.report()
        self.assertEqual(ex.exception.code, dot_h.st_return_code.RUNNER_NOT_READY_TO_REPORT)

        self.stapi.runner.run()
        self.stapi.runner.report()

    def test_report_native(self): self.runner_report(dot_h.st_runner_type_t.NATIVE)
    def test_report_embree(self): self.runner_report(dot_h.st_runner_type_t.EMBREE)
    def test_report_optix(self):  self.runner_report(dot_h.st_runner_type_t.OPTIX)

class Results:
    def test_csv_dump(self):
        self.stapi.result.csv.dump('./sample.csv')

    def test_locations(self):
        loc_x, loc_y, loc_z = self.stapi.result.locations()
        self.assertEqual(len(loc_x), self.n_intersections)
        self.assertEqual(len(loc_y), self.n_intersections)
        self.assertEqual(len(loc_z), self.n_intersections)

    def test_cosines(self):
        coz_x, coz_y, coz_z = self.stapi.result.cosines()
        self.assertEqual(len(coz_x), self.n_intersections)
        self.assertEqual(len(coz_y), self.n_intersections)
        self.assertEqual(len(coz_z), self.n_intersections)

    def test_elementmap(self):
        element_map = self.stapi.result.elementmap()
        self.assertEqual(len(element_map), self.n_intersections)

    def test_stagemap(self):
        stage_map = self.stapi.result.stagemap()
        self.assertEqual(len(stage_map), self.n_intersections)

    def test_raynumbers(self):
        ray_numbers = self.stapi.result.raynumbers()
        self.assertEqual(len(ray_numbers), self.n_intersections)

    def test_sun_stats(self):
        width, height, area, nsunrays = self.stapi.result.sun_stats()

        if self.runner_type != dot_h.st_runner_type_t.OPTIX: self.assertGreater(width, 0)
        else:                                                self.assertEqual(width, 0)
        if self.runner_type != dot_h.st_runner_type_t.OPTIX: self.assertGreater(height, 0)
        else:                                                self.assertEqual(height, 0)
        
        self.assertGreater(area, 0)
        self.assertGreater(nsunrays, 0)

    def test_get_results_data(self):
        res = self.stapi.result.get()
        for k in res.keys():
            self.assertEqual(len(res[k][:self.n_intersections]), self.n_intersections)

class ResultsNativeTests(STAPIv2TestCase, Results):
    def setUp(self):
        super().setUp()
        self.runner_type = dot_h.st_runner_type_t.NATIVE
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(self.runner_type)
        self.stapi.runner.run()
        self.stapi.runner.report()
        self.n_intersections = self.stapi.result.num()

@unittest.skipIf(not api.is_runner_installed(dot_h.st_runner_type_t.EMBREE),
                 "Embree runner not installed, skip testing EmbreeRunner results.")
class ResultsEmbreeTests(STAPIv2TestCase, Results):
    def setUp(self):
        super().setUp()
        self.runner_type = dot_h.st_runner_type_t.EMBREE
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(self.runner_type)
        self.stapi.runner.run()
        self.stapi.runner.report()
        self.n_intersections = self.stapi.result.num()

@unittest.skipIf(not api.is_runner_installed(dot_h.st_runner_type_t.OPTIX),
                 "Optix runner not installed, skip testing OptixRunner results.")
class ResultsOptixTests(STAPIv2TestCase, Results):
    def setUp(self):
        super().setUp()
        self.runner_type = dot_h.st_runner_type_t.OPTIX
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(self.runner_type)
        self.stapi.runner.run()
        self.stapi.runner.report()
        self.n_intersections = self.stapi.result.num()

class BatchTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()
        # set up dummy values
        self.sim_params = _STS.simulation_parameters(35.962278, -106.5122622, 1, 100, .1, True, True, False)

        self.args_sun_buie = _STS.sun(0, 2, 2, 2, .5, b'b')
        self.args_sun_user = _STS.sun(3, 608, 303, 1000, 5, b' ')
        self.good_angles      = [0, 1, 2]
        self.good_intensities = [0, 1, 2]
        
        self.opt_set = _STS.optical_properties_set(b'dummy', 1.1, 1.1, 2)
        self.front   = _STS.optical_properties_face(.5, .5, 5, 5, b'g')
        self.back    = _STS.optical_properties_face(.25, .25, 2, 2, b'g')
        
        self.el_args  = _STS.element([2, 2, 2], [2, 2, 2], 2, False, True)
        self.opt_id   = 0
        self.geom_args = _STS.geometry([2], [2, 2], b'c', b'p')

        self.calc = dot_h.SolarPositionCalculationMethod.SPA
        self.loc = _STS.sun_location(40.0, -105.0, -7.0)
        self.dt  = _STS.sun_datetime(2025, 6, 20)
        self.az = 178.61128380
        self.el = 73.439035265
        self.zen = 90 - self.el
        self.sun_vector = Point(0.006908, -0.2849516, 0.9585169)

    def test_simple(self):
        self.stapi.batch.data.json.load(locate('sample.json'))
        el_rec = self.stapi.batch.data.element.num()
        self.stapi.batch.runner.setup(dot_h.st_runner_type_t.NATIVE)
        self.stapi.batch.runner.run()
        self.stapi.batch.runner.report()
        inter_rec = self.stapi.batch.result.num()

        # test value is None before batching
        self.assertEqual(el_rec.value, None)
        self.assertEqual(inter_rec.value, None)

        self.stapi.batch()
        self.assertEqual(el_rec.value, 126)
        self.assertEqual(inter_rec.value, self.stapi.result.num())

    # functions for simulation data management thru json strings
    def test_call_st_read_input_json(self):
        f = open(locate('sample.json'), mode='rb')
        _json = f.read()
        f.close()
        self.stapi.batch.data.json.load(_json)
        self.stapi.batch()
        self.assertEqual(self.stapi.data.element.num(), 126)
    
    def test_call_st_read_input_json_file(self):
        self.stapi.batch.data.json.load(locate('sample.json'))
        self.stapi.batch()
        self.assertEqual(self.stapi.data.element.num(), 126)
    
    # functions for simulation data management directly
    def test_call_st_set_simulation_parameters(self):
        self.stapi.batch.parameters.set(self.sim_params)

        # batch call and get/test record value 
        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertDictEqual(rt_params, self.sim_params.ctype.value)

    def test_call_st_sim_rays(self):
        self.stapi.batch.parameters.rays(10, 1000)

        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['number_of_rays'], 10)
        self.assertEqual(rt_params['max_number_of_rays'], 1000)

    def test_call_st_sim_power_tower(self):
        self.stapi.batch.parameters.power_tower(True)

        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['as_power_tower'], True)

    def test_call_st_sim_errors(self):
        self.stapi.batch.parameters.errors(True, True)

        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['include_sun_shape_errors'], True)
        self.assertEqual(rt_params['include_optical_errors'], True)

    def test_call_st_sim_location(self): 
        self.stapi.batch.parameters.location(35.962278, -106.5122622)

        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['latitude'], 35.962278)
        self.assertEqual(rt_params['longitude'], -106.5122622)

    def test_call_st_sim_tolerance(self):
        self.stapi.batch.parameters.tolerance(.1)
        
        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['tolerance'], .1)

    def test_call_st_sim_params(self):
        self.stapi.batch.legacy.sim_params(10, 1000, True)
        
        self.stapi.batch()
        rt_params = self.stapi.parameters.get()
        self.assertEqual(rt_params['number_of_rays'], 10)
        self.assertEqual(rt_params['max_number_of_rays'], 1000)
        self.assertEqual(rt_params['as_power_tower'], True)
        
    # functions to add/remove/set optical properties
    def test_call_st_num_optics(self):
        rec = self.stapi.batch.data.optic.num()

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        self.assertEqual(rec.value, 0)

        self.stapi.batch.data.optic.add(self.opt_set, self.front, self.back)
        rec = self.stapi.batch.data.optic.num()
        self.stapi.batch()
        self.assertEqual(rec.value, 1)

    def test_call_st_add_optical_properies_set(self):
        rec = self.stapi.batch.data.optic.add(self.opt_set, self.front, self.back)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        self.assertEqual(rec.value, 0)

    def test_call_st_delete_optic(self):
        id = self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.batch.data.optic.delete(id)

        self.stapi.batch()
        self.assertEqual(self.stapi.data.optic.num(), 0)
        
    def test_call_st_clear_optics(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.batch.data.optic.clear()

        self.stapi.batch()
        self.assertEqual(self.stapi.data.optic.num(), 0)

    # functions to add/remove elements
    def test_call_st_num_elements(self):
        rec = self.stapi.batch.data.element.num()
        
        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        self.assertEqual(rec.value, 0)

        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)
        rec = self.stapi.batch.data.element.num()
        self.stapi.batch()
        self.assertEqual(rec.value, 1)

    def test_call_st_add_element(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        rec = self.stapi.batch.data.element.add(self.el_args, self.opt_id, self.geom_args)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        self.assertEqual(rec.value, 1)

    def test_call_st_delete_element(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.delete(id)
        
        self.stapi.batch()
        self.assertEqual(self.stapi.data.element.num(), 0)

    def test_call_st_clear_elements(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.clear()
        
        self.stapi.batch()
        self.assertEqual(self.stapi.data.element.num(), 0)
    
    # functions to modify elements
    def test_call_st_element_enabled(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.enabled(id, True)

        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['enabled_flag'], True)
        
    def test_call_st_element_virtual(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.virtual(id, False)

        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['virtual_flag'], False)

    def test_call_st_element_zyx(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.xyz(id, [4, 4, 4])
        
        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['x'], 4)
        self.assertEqual(el['y'], 4)
        self.assertEqual(el['z'], 4)

    def test_call_st_element_aim(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.aim(id, [4, 4, 4])
        
        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['ax'], 4)
        self.assertEqual(el['ay'], 4)
        self.assertEqual(el['az'], 4)

    def test_call_st_element_zrot(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.zrot(id, 4)
        
        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['zrot'], 4)

    def test_call_st_element_aperture(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        new_ap = b'r'
        new_params = [3, 3]

        self.stapi.batch.data.element.aperture(id, new_ap, new_params)

        self.stapi.batch()
        _, _, geom = self.stapi.data.element.get(id)
        self.assertEqual(geom['ap'], new_ap)
        self.assertEqual(geom['a_params'][0], new_params[0])
        self.assertEqual(geom['a_params'][1], new_params[1])

    def test_call_st_element_surface(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        new_surf = b's'
        new_params = [3]

        self.stapi.batch.data.element.surface(id, new_surf, new_params)

        self.stapi.batch()
        _, _, geom = self.stapi.data.element.get(id)
        self.assertEqual(geom['surf'], new_surf)
        self.assertEqual(geom['s_params'][0], new_params[0])

    def test_call_st_element_optic(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        other_opt = _STS.optical_properties_set(b'other', 1.1, 1.1, 0)
        other_id = self.stapi.data.optic.add(other_opt, self.front, self.back)

        self.stapi.batch.data.element.optic(id, other_id)

        self.stapi.batch()
        _, test_id, _ = self.stapi.data.element.get(id)
        self.assertEqual(test_id, other_id)

    def test_call_st_element_group(self):
        self.stapi.data.optic.add(self.opt_set, self.front, self.back)
        id = self.stapi.data.element.add(self.el_args, self.opt_id, self.geom_args)

        self.stapi.batch.data.element.group(id, 1)

        self.stapi.batch()
        el, *_ = self.stapi.data.element.get(id)
        self.assertEqual(el['group'], 1)

    # sun functions
    def test_call_st_add_sun(self):
        self.stapi.batch.data.sun.add(self.args_sun_user, self.good_angles, self.good_intensities)

        self.stapi.batch()
        sun, *_ = self.stapi.data.sun.get(3)
        self.assertEqual(sun['npoints'], 3)
        self.assertEqual(sun['x'], 608)
        self.assertEqual(sun['y'], 303)
        self.assertEqual(sun['z'], 1000)
        self.assertEqual(sun['shape'], b'd')

        self.stapi.batch.data.sun.add(self.args_sun_buie)

        self.stapi.batch()
        sun, *_ = self.stapi.data.sun.get()
        self.assertDictEqual(sun, self.args_sun_buie.ctype.value)
    
    def test_call_st_sun_shape(self):
        self.stapi.batch.data.sun.shape('g', 5)

        self.stapi.batch()

        sun, *_ = self.stapi.data.sun.get()
        self.assertEqual(sun['shape'], b'g')
        self.assertEqual(sun['sigma_halfwidth_csr'], 5)

    def test_call_st_sun_xyz(self):
        self.stapi.data.sun.add(self.args_sun_buie)
        
        self.stapi.batch.data.sun.xyz(608, 303, 1000)

        self.stapi.batch()
        sun, *_ = self.stapi.data.sun.get()
        self.assertEqual(sun['x'], 608)
        self.assertEqual(sun['y'], 303)
        self.assertEqual(sun['z'], 1000)

    def test_call_st_sun_userdata(self):
        self.stapi.data.sun.add(self.args_sun_buie)

        self.stapi.batch.data.sun.userdata(3, self.good_angles, self.good_intensities)

        self.stapi.batch()
        sun, *_ = self.stapi.data.sun.get(3)
        self.assertEqual(sun['npoints'], 3)
        self.assertEqual(sun['shape'], b'd')

    def test_call_st_get_sun_az_zen(self):
        rec = self.stapi.batch.data.sun.az_zen(self.calc, self.loc, self.dt)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        self.stapi.batch()
        az, zen = rec.value
        self.assertAlmostEqual(az, self.az, 7)
        self.assertAlmostEqual(zen, self.zen, 7)

    def test_call_st_get_sun_az_el(self):
        rec = self.stapi.batch.data.sun.az_el(self.calc, self.loc, self.dt)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        self.stapi.batch()
        az, el = rec.value
        self.assertAlmostEqual(az, self.az, 7)
        self.assertAlmostEqual(el, self.el, 7)

    def test_call_st_get_sun_vector(self):
        rec = self.stapi.batch.data.sun.vector(self.calc, self.loc, self.dt)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        self.stapi.batch()
        v = rec.value
        self.assertAlmostEqual((v - self.sun_vector).radius(), 0, 6)
    
    # functions for SolTrace runner management
    def test_call_st_sim_setup(self):
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.batch.runner.setup(dot_h.st_runner_type_t.NATIVE)
        self.stapi.batch()

    def test_call_st_sim_run_v2(self):
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(dot_h.st_runner_type_t.NATIVE)
        self.stapi.batch.runner.run()
        self.stapi.batch()

    def test_call_st_sim_report(self):
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(dot_h.st_runner_type_t.NATIVE)
        self.stapi.runner.run()
        self.stapi.batch.runner.report()
        self.stapi.batch()
        
    # functions for SolTrace results management
    def set_up_run_report(self):
        self.stapi.data.json.load(locate('sample.json'))
        self.stapi.parameters.rays(1000, 10000)
        self.stapi.runner.setup(dot_h.st_runner_type_t.NATIVE)
        self.stapi.runner.run()
        self.stapi.runner.report()
        return self.stapi.result.num()
    
    def test_call_st_write_results_csv(self):
        self.set_up_run_report()
        self.stapi.batch.result.csv.dump('./batch_sample.csv')
        self.stapi.batch()

    # functions to get results directly
    def test_call_st_locations(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.locations(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        loc_x, loc_y, loc_z = rec.value
        self.assertEqual(len(loc_x), n)
        self.assertEqual(len(loc_y), n)
        self.assertEqual(len(loc_z), n)

    def test_call_st_cosines(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.cosines(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        cos_x, cos_y, cos_z = rec.value
        self.assertEqual(len(cos_x), n)
        self.assertEqual(len(cos_y), n)
        self.assertEqual(len(cos_z), n)

    def test_call_st_elementmap(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.elementmap(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        el = rec.value
        self.assertEqual(len(el), n)
        
    def test_call_st_stagemap(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.stagemap(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        stage = rec.value
        self.assertEqual(len(stage), n)

    def test_call_st_raynumbers(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.raynumbers(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        rays = rec.value
        self.assertEqual(len(rays), n)

    def test_call_st_sun_stats(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.sun_stats()

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        width, height, area, nsunrays = rec.value
        self.assertGreater(width, 0.)
        self.assertGreater(height, 0.)
        self.assertGreater(area, 0.)
        self.assertGreater(nsunrays, 0)

    def test_call_st_get_results_data(self):
        n = self.set_up_run_report()
        rec = self.stapi.batch.result.get(n)

        # test value is None before batching
        self.assertEqual(rec.value, None)

        # batch call and get/test record value 
        self.stapi.batch()
        results = rec.value
        for k in results.keys():
            self.assertEqual(len(results[k][:n]), n)

class LegacyTests(STAPIv2TestCase):
    def setUp(self):
        super().setUp()

        self.legacy = PySolTrace(testing = True)
        # set up dummy values
        # sun
        self.sun = self.legacy.add_sun()
        self.sun.position = Point(608, 303, 1000)
        self.sun.shape = 'b'
        self.sun.sigma = .5
        self.sun_args = _STS.sun(0, 608, 303, 1000, .5, 'b'.encode())

        # optical property set
        self.opt = self.legacy.add_optic('dummy')

        # front face object (f is 6th letter)
        self.opt.front.dist_type       = 'f'
        self.opt.front.refraction_real = 1.6
        self.opt.front.reflectivity    = .6
        self.opt.front.transmissivity  = .6
        self.opt.front.slope_error     = 6
        self.opt.front.spec_error      = 6

        # back face object (b is 2nd letter)
        self.opt.back.dist_type       = 'p'
        self.opt.back.refraction_real = 1.2
        self.opt.back.reflectivity    = .2
        self.opt.back.transmissivity  = .2
        self.opt.back.slope_error     = 2
        self.opt.back.spec_error      = 2

        self.opt_args   = _STS.optical_properties_set("dummy".encode(), 1.6, 1.2, 2)
        self.front_args = _STS.optical_properties_face(.6, .6, 6, 6, 'f'.encode())
        self.back_args  = _STS.optical_properties_face(.2, .2, 2, 2, 'p'.encode())

        # dummy element (e is 5th letter)
        self.stage = self.legacy.add_stage()
        self.el = self.stage.add_element()

        self.el.optic    = self.opt
        self.el.position = Point(5, 5, 5)
        self.el.aim      = Point(5, 5, 5)
        self.el.zrot     = 5
        self.el.aperture_circle(5)
        self.el.surface_parabolic(5, 5)

        self.el_args = _STS.element([5, 5, 5], [5, 5, 5], 5, True, False)
        self.opt_id = 0
        self.geom_args = _STS.geometry([5], [5, 5], b'c', b'p')

    def test_legacy_sun(self):
        self.sun.Create(self.stapi, _do = True)

        rt_sun_args, *_ = self.stapi.data.sun.get()
        self.assertDictEqual(rt_sun_args, self.sun_args.ctype.value)

    def test_legacy_optic(self):
        self.opt.Create(self.stapi, _do = True)

        _opt_set, _front, _back = self.stapi.data.optic.get(0)
        self.assertDictEqual(_opt_set, self.opt_args.ctype.value)
        self.assertDictEqual(_front, self.front_args.ctype.value)
        self.assertDictEqual(_back, self.back_args.ctype.value)

    def test_legacy_element(self):
        self.opt.Create(self.stapi, _do = True)
        self.el.Create(self.stapi, _do = True)

        args, optic_id, geom = self.stapi.data.element.get(1)
        
        self.assertDictEqual(args, self.el_args.ctype.value)
        self.assertEqual(optic_id, self.opt_id)
        self.assertEqual(geom['a_params'][0], self.geom_args.a_params[0])
        self.assertEqual(geom['s_params'][0], 1 / (2 * self.geom_args.s_params[0]))
        self.assertEqual(geom['s_params'][1], 1 / (2 * self.geom_args.s_params[1]))

class TowerDemoTest(STAPIv2TestCase):
    def setUp(self):
        super().setUp()

    def test(self):
        PT = PySolTrace(testing = True)
        # this object has data from previous legacy tests?
        PT.clear_optics()
        PT.clear_stages()
        # Create two optics types - one for reflector, and one for absorber.
        opt_ref = PT.add_optic("Reflector")
        opt_ref.front.reflectivity = 1.
        opt_abs = PT.add_optic("Absorber")
        opt_abs.front.reflectivity = 0.

        # Sun
        sun = PT.add_sun()
        # Give sun an arbitrary position
        sun.position.x = 1.
        sun.position.y = -1.
        sun.position.z = 99.

        # Reflector stage
        st = PT.add_stage()

        #absorber element height
        abs_pos = Point(0., 0., 10.)

        # Create a heliostat at some random x,y position, reflecting to the receiver
        for i in range(-1,2):
            hpos = [sin(i*pi/2)*5, cos(i*pi/2)*5]
            # hpos = [random.uniform(-10,10), random.uniform(-10,10)]
            el = st.add_element()
            el.optic = opt_ref
            el.position.x = hpos[0]
            el.position.y = hpos[1]
            # calculate the vectors - receiver, sun, and aim
            rvec = (abs_pos - el.position).unitize()
            svec = sun.position.unitize()
            avec = (rvec + svec)/2.
            # assign the aim vector. scale by a large number
            el.aim = el.position + avec*100.
            # compute surface z rotation to align with plane of the ground
            el.zrot = PT.util_calc_zrot_azel(avec)
            # Set surface and aperture characteristics
            el.surface_flat()
            el.aperture_rectangle(1.0,1.95)

        sta = PT.add_stage()
        # flat absorber element
        ela = sta.add_element()
        ela.position = abs_pos
        ela.aim.z = 0.
        ela.aim.x = 0.
        ela.aim.y = 5.
        ela.optic = opt_abs
        ela.surface_flat()
        ela.aperture_rectangle(2,2)  

        # set simulation parameters
        PT.num_ray_hits = 1e3
        PT.max_rays_traced = PT.num_ray_hits*100
        PT.is_sunshape = True 
        PT.is_surface_errors = True

        PT.run(-1, True, 1)

if __name__ == '__main__':
    # print(f'\n\n\n\n{found_in(dot_h)}\n\n\n\n')
    unittest.main()