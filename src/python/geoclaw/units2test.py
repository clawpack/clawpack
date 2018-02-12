from units2 import convertUnits
import unittest

class TestTempConversion(unittest.TestCase):
    def test_F_to_C(self):   #assertAlmostEqual checks accuracy within 7 decimal places
        self.assertAlmostEqual(convertUnits("t", 'F', 'C', 190.00004), 87.7778)
        self.assertAlmostEqual(convertUnits("t", 'F', 'C', -11.99992), -24.4444)

    def test_C_to_F(self):
        self.assertAlmostEqual(convertUnits("t", 'C', 'F', 120), 248)
        self. assertAlmostEqual(convertUnits("t", 'C', 'F', -12442.222), -22363.9996)

    def test_K_to_C(self):
        self.assertAlmostEqual(convertUnits("t", 'K', 'C', 150), -123.15)
        self.assertAlmostEqual(convertUnits("t", 'K', 'C', -129), -402.15)

    def test_K_to_F(self):
        self.assertAlmostEqual(convertUnits("t", 'K', 'F', 120), -243.67)
        self.assertAlmostEqual(convertUnits("t", 'K', 'F', -122), -679.27)

    def test_C_to_K(self):
        self.assertAlmostEqual(convertUnits("t", 'C', 'K', 400), 673.15)
        self.assertAlmostEqual(convertUnits("t", 'C', 'K', -223), 50.15)

    def test_F_to_K(self):
        self.assertAlmostEqual(convertUnits("t", 'F', 'K', 212), 373.15)
        self.assertAlmostEqual(convertUnits("t", 'F', 'K', -122.0008), 187.594)

class TestPressureConversion(unittest.TestCase):
    def test_Pa_to_kPa(self):
        self.assertAlmostEqual(convertUnits("p", 'Pa', 'kPa', .01), 1e-5)
        self.assertAlmostEqual(convertUnits("p", 'Pa', 'kPa', 9124), 9.124)
    
    def test_GPa_to_mbar(self):
        self.assertAlmostEqual(convertUnits("p", 'GPa', 'mbar', 912), 9120000000)
        self.assertAlmostEqual(convertUnits("p", 'GPa', 'mbar', .012), 120000) 

    def test_dyne_per_cm2_to_kPa(self):
        self.assertAlmostEqual(convertUnits("p", 'dyne/cm^2', 'kPa', 1234), .1234)
        self.assertAlmostEqual(convertUnits("p", 'dyne/cm^2', 'kPa', .01234), 1.234e-6)

class TestLengthConversion(unittest.TestCase):
    def test_m_to_miles(self):
        self.assertAlmostEqual(convertUnits("l", 'm', 'miles', 1.99999), 0.001242742)
    def test_km_to_m(self):
        self.assertEqual(convertUnits("l", 'km', 'm', 12345.0000512), 12345000.0512)

class TestSpeedConversion(unittest.TestCase):
    def test_m_per_s_to_knots(self):
        self.assertAlmostEqual(convertUnits("s", 'm/s', 'knots', 0.999997688), 1.9438399999515)
    def test_knots_to_m_per_s(self):
        self.assertAlmostEqual(convertUnits("s", 'knots', 'm/s', 3.9999963533499998114), 2.0577759000011681678)

class TestMomentConversion(unittest.TestCase):
    def test_N_m_to_dyne_cm(self):
        self.assertAlmostEqual(convertUnits("m", 'N-m', 'dyne-cm', 1234), 12340000000)
    def test_dyne_cm_to_N_m(self):
        self.assertAlmostEqual(convertUnits("m", 'dyne-cm', 'N-m', 12340000), 1.2340)

if __name__ == '__main__':
    unittest.main()
