import unittest
import numpy as np
import calc_grain_size


class TestCalcGrainSize(unittest.TestCase):
    def test_G_numgrain(self):
        self.assertEqual(calc_grain_size.G_numgrain(5), 19.299297468431693)

    def test_G_meanbarA(self):
        self.assertEqual(calc_grain_size.G_meanbarA(5), 14.65544127865697)

    def test_G_meanintl(self):
        self.assertEqual(calc_grain_size.G_meanintl(5), 11.999999999999998)

    def test_get_polygon_area_shoelace_formula(self):
        self.x = np.array([4.87336245, 4.87336245, 196.87336245, 196.87336245])
        self.y = np.array([4.87336245, 196.87336245, 196.87336245, 4.87336245])

        self.assertEqual(calc_grain_size.get_polygon_area_shoelace_formula(self.x, self.y), 36864.00000000001)


if __name__ == '__main__':
    unittest.main()
