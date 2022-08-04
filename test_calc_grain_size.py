import numpy as np
import pickle
import calc_grain_size
from calc_grain_size import EBSD
import pytest


class TestClass:
    def test_G_numgrain(self):
        assert calc_grain_size.G_numgrain(5) == 19.299297468431693

    def test_G_meanbarA(self):
        assert calc_grain_size.G_meanbarA(5) == 14.65544127865697

    def test_G_meanintl(self):
        assert calc_grain_size.G_meanintl(5) == 11.999999999999998

    def test_get_polygon_area_shoelace_formula(self):
        self.x = np.array([4.87336245, 4.87336245, 196.87336245, 196.87336245])
        self.y = np.array([4.87336245, 196.87336245, 196.87336245, 4.87336245])

        assert calc_grain_size.get_polygon_area_shoelace_formula(self.x, self.y) == 36864.00000000001

    # def test_grainsize_linint_random(self):
    #     with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
    #         self.ebsd = pickle.load(pickle_in)
    #
    #     expected_results =
    #     actual_results = calc_grain_size.grainsize_linint_random(self.ebsd)
    #
    #     np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E112_SaltikovPlanimetric(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = 8.092847, 0.002115885, 44
        actual_results = calc_grain_size.GrainSize_E112_SaltikovPlanimetric(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E112_JeffriesPlanimetric(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = 8.305933, 0.00245266, 36
        actual_results = calc_grain_size.GrainSize_E112_JeffriesPlanimetric(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E2627_AsWritten(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = 7.847896085958176, 560.0738607832295, 48, 0.0036676302599370552, 720.0612244897765, \
            np.array([[309.68135619], [498.84632253], [623.93928415], [139.58543887], [546.90032608], [378.32993269],
                      [1129.65046433], [295.95164089], [549.95137392], [610.97233081], [189.16496634], [589.6149959],
                      [476.72622566], [970.23321447], [305.86754639], [234.16792205], [320.36002365], [496.55803665],
                      [704.79205202], [291.37506913], [262.3901146], [770.38958067], [732.25148262], [1025.91483763],
                      [865.7348258], [675.04433554], [554.52794569], [694.11338457], [596.47985355], [504.18565626],
                      [895.48254229], [563.68108922], [1112.86970119], [858.10720619], [234.16792205], [318.07173776],
                      [620.12547434], [154.07791613], [1301.27190557], [384.43202837], [356.97259778], [643.00833317],
                      [406.55212525], [454.60612879], [170.09591732], [793.2724395],  [716.23348144], [556.81623157]])
        actual_results = calc_grain_size.GrainSize_E2627_AsWritten(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results[0:5], expected_results[0:5])
        np.testing.assert_array_almost_equal(len(actual_results[5]), len(expected_results[5]))
        np.testing.assert_array_almost_equal(np.shape(actual_results[5]), np.shape(expected_results[5]))

    def test_GrainSize_E2627_CustomMinGS(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = 7.876088800924972, 549.2353116758231, 49, 0.0037440392236857438, 720.0612244897765, \
            np.array([[28.98495452], [309.68135619], [498.84632253], [623.93928415], [139.58543887], [546.90032608],
                      [378.32993269], [1129.65046433], [295.95164089], [549.95137392], [610.97233081], [189.16496634],
                      [589.6149959], [476.72622566], [970.23321447], [305.86754639], [234.16792205], [320.36002365],
                      [496.55803665], [704.79205202], [291.37506913], [262.3901146], [770.38958067], [732.25148262],
                      [1025.91483763], [865.7348258], [675.04433554], [554.52794569], [694.11338457], [596.47985355],
                      [504.18565626], [895.48254229], [563.68108922], [1112.86970119], [858.10720619], [234.16792205],
                      [318.07173776], [620.12547434], [154.07791613], [1301.27190557], [384.43202837], [356.97259778],
                      [643.00833317], [406.55212525], [454.60612879], [170.09591732], [793.2724395], [716.23348144],
                      [556.81623157]])
        actual_results = calc_grain_size.GrainSize_E2627_CustomMinGS(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results[0:5], expected_results[0:5])
        np.testing.assert_array_almost_equal(len(actual_results[5]), len(expected_results[5]))
        np.testing.assert_array_almost_equal(np.shape(actual_results[5]), np.shape(expected_results[5]))

    # def test_GrainSize_E112_HeynRandomLineMLI(self):
    #     with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
    #         self.ebsd = pickle.load(pickle_in)
    #
    #     expected_results =
    #     actual_results = calc_grain_size.GrainSize_E112_HeynRandomLineMLI(self.ebsd)
    #
    #     np.testing.assert_array_almost_equal(actual_results, expected_results)

    # def test_GrainSize_E112_HeynRandomLinePL(self):
    #     with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
    #         self.ebsd = pickle.load(pickle_in)
    #
    #     expected_results =
    #     actual_results = calc_grain_size.GrainSize_E112_HeynRandomLinePL(self.ebsd)
    #
    #     np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E112_Hilliard(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = 6.148001162887553, 38.0, [15.87331025], [603.18578949]
        actual_results = calc_grain_size.GrainSize_E112_Hilliard(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E112_Abrams(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        expected_results = np.array([24.85820742]), 70, np.array([17.2338797]), np.array([1206.37157898])
        actual_results = calc_grain_size.GrainSize_E112_Abrams(self.ebsd)

        np.testing.assert_array_almost_equal(actual_results, expected_results)

    def test_GrainSize_E930_ALA(self):
        with open(f"C:\\git\\measureGrainSize\\myHighResEBSD.pickle", "rb") as pickle_in:
            self.ebsd = pickle.load(pickle_in)

        G_S, _, _ = calc_grain_size.GrainSize_E112_SaltikovPlanimetric(self.ebsd)

        expected_results = (np.array([6.63166264]), np.array([0]))
        actual_results = calc_grain_size.GrainSize_E930_ALA(self.ebsd, G_S)

        np.testing.assert_array_almost_equal(actual_results, expected_results)


if __name__ == '__main__':
    pytest.main()
