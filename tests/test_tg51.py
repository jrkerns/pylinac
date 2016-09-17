import unittest

from pylinac import tg51


class TestFunctions(unittest.TestCase):

    def test_p_tp(self):
        temps = (22, 25, 19)
        presss = (760, 770, 740)
        expected_ptp = (1.0, 0.997, 1.0165)
        for temp, press, exp in zip(temps, presss, expected_ptp):
            self.assertAlmostEqual(tg51.p_tp(temp, press), exp, delta=0.001)

    def test_p_ion(self):
        low_vals = (20, 20.05)
        high_vals = (20, 20.1)
        expected_pion = (1.0, 1.0025)
        for low, high, exp in zip(low_vals, high_vals, expected_pion):
            self.assertAlmostEqual(tg51.p_ion(300, 150, high, low), exp, delta=0.001)

    def test_p_pol(self):
        m_ref = (20, -20.2, 19.8)
        m_opp = (-20, 19.8, -20.1)
        expected_ppol = (1.0, 0.99, 1.0075)
        for ref, opp, exp in zip(m_ref, m_opp, expected_ppol):
            self.assertAlmostEqual(tg51.p_pol(ref, opp), exp, delta=0.001)

    def test_pddx(self):
        pdds = (66.4, 73.3, 76.7, 79.3)
        energies = (6, 10, 15, 18)
        pddxs = (66.4, 73.3, 77.18, 80.47)
        for pdd, energy, pddx in zip(pdds, energies, pddxs):
            self.assertAlmostEqual(tg51.pddx(pdd, energy), pddx, delta=0.01)

    def test_kq(self):
        models = ('30010', 'A12')
        pddxs = (66.4, 76.7)
        kqs = (0.9927, 0.976)
        for model, pddx, kq in zip(models, pddxs, kqs):
            self.assertAlmostEqual(tg51.kq(model, pddx), kq, delta=0.001)


class TestTG51(unittest.TestCase):

    def test_tg51(self):
        tg51_6x = tg51.TG51Photon(temp=22, press=751.2, model='30013', n_dw=5.45,
                                  m_raw=(24.82,), m_opp=-24.83, m_low=24.79, measured_pdd=66.8,
                                  clinical_pdd=66.8)
        self.assertAlmostEqual(tg51_6x.dose_mu_10, 0.679, delta=0.001)
        self.assertAlmostEqual(tg51_6x.dose_mu_dmax, 1.0167, delta=0.001)
