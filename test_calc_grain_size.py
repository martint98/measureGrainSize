import unittest
import calc_grain_size


class TestCalcGrainSize(unittest.TestCase):
    def test_G_numgrain(self):
        self.assertEqual(calc_grain_size.G_numgrain(5), 19.299297468431693)

    def test_G_meanbarA(self):
        self.assertEqual(calc_grain_size.G_meanbarA(5), 14.65544127865697)

    def test_G_meanintl(self):
        self.assertEqual(calc_grain_size.G_meanintl(5), 11.999999999999998)


if __name__ == '__main__':
    unittest.main()
