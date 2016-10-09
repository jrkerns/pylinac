import unittest

from pylinac import tg51


class TestFunctions(unittest.TestCase):

    def test_p_tp(self):
        temps = (22, 25, 19)
        presss = (760, 770, 740)
        expected_ptp = (1.0, 0.997, 1.0165)
        for temp, press, exp in zip(temps, presss, expected_ptp):
            self.assertAlmostEqual(tg51.p_tp(temp, press), exp, delta=0.001)

    def test_p_pol(self):
        m_ref = (20, -20.2, 19.8)
        m_opp = (-20, 19.8, -20.1)
        expected_ppol = (1.0, 0.99, 1.0075)
        for ref, opp, exp in zip(m_ref, m_opp, expected_ppol):
            self.assertAlmostEqual(tg51.p_pol(ref, opp), exp, delta=0.001)

    def test_p_ion(self):
        low_vals = (20, 20.05)
        high_vals = (20, 20.1)
        expected_pion = (1.0, 1.0025)
        for low, high, exp in zip(low_vals, high_vals, expected_pion):
            self.assertAlmostEqual(tg51.p_ion(300, 150, high, low), exp, delta=0.001)

    def test_dref(self):
        i50s = (3, 5, 7)
        drefs = (1.72, 2.96, 4.19)
        for i50, dref in zip(i50s, drefs):
            self.assertAlmostEqual(tg51.d_ref(i50), dref, delta=0.01)

    def test_r50(self):
        i50s = (3.5, 5.5, 12)
        r50s = (3.54, 5.60, 18.71)
        for i50, r50 in zip(i50s, r50s):
            self.assertAlmostEqual(tg51.r_50(i50), r50, delta=0.01)

    def test_kp_r50(self):
        r50s = (3.5, 5, 8)
        kpr50s = (1.017, 1, 1)
        for r50, kpr50 in zip(r50s, kpr50s):
            self.assertAlmostEqual(tg51.kp_r50(r50), kpr50, delta=0.01)

        with self.assertRaises(ValueError):
            tg51.kp_r50(9.1)
        with self.assertRaises(ValueError):
            tg51.kp_r50(1.8)

    def test_pq_gr(self):
        m_refs = (20, 25, 30)
        m_ref_pluss = (19.9, 23, 30)
        pq_grs = (0.995, 0.92, 1.0)
        for m_ref, m_ref_plus, pq_gr in zip(m_refs, m_ref_pluss, pq_grs):
            self.assertAlmostEqual(tg51.pq_gr(m_ref_plus, m_ref), pq_gr, delta=0.01)

    def test_m_corr(self):
        exp = 20.225
        res = tg51.m_corrected(1.01, 0.995, 1, 1.005, (20, 20.05))
        self.assertAlmostEqual(exp, res, delta=0.002)

    def test_pddx(self):
        pdds = (66.4, 70.5, 72.8, 73.3, 76.7, 77.1, 77.1, 79.3)
        energies = (6, 10, 10, 10, 15, 15, 15, 18)
        pddxs = (66.4, 70.5, 72.8, 72.87, 77.18, 77.57, 78.27, 80.47)
        foils = (None, '30cm', '50cm', None, None, '50cm', '30cm', None)
        for pdd, energy, pddx, foil in zip(pdds, energies, pddxs, foils):
            self.assertAlmostEqual(tg51.pddx(pdd, energy, foil), pddx, delta=0.01)

    def test_kq(self):
        # Test via PDDs
        models = ('30010', 'A12')
        pddxs = (66.4, 76.7)
        kqs = (0.9927, 0.976)
        for model, pddx, kq in zip(models, pddxs, kqs):
            self.assertAlmostEqual(tg51.kq(model, pddx), kq, delta=0.001)

        # test via TPRs
        tprs = (0.65, 0.76)
        kqs = (0.994, 0.975)
        for model, tpr, kq in zip(models, tprs, kqs):
            self.assertAlmostEqual(tg51.kq(model, tpr=tpr), kq, delta=0.001)

        # neither TPR or PDD passed
        with self.assertRaises(ValueError):
            tg51.kq()

        # both defined
        with self.assertRaises(ValueError):
            tg51.kq(pddx=0.66, tpr=0.72)

        # PDD too low
        with self.assertRaises(ValueError):
            tg51.kq(pddx=61)

        # TPR too high
        with self.assertRaises(ValueError):
            tg51.kq(tpr=81)


class TestTG51(unittest.TestCase):

    def test_tg51_photon(self):
        # MDA TB2 2015
        tg51_15x = tg51.TG51Photon(temp=20.5, press=760, model='30013', n_dw=5.444,
                                   p_elec=1.002, m_raw=29.28, m_opp=-29.33, m_low=29.10,
                                   measured_pdd=76.9, clinical_pdd=77.4)
        self.assertAlmostEqual(tg51_15x.dose_mu_10, 0.780, delta=0.001)
        self.assertAlmostEqual(tg51_15x.dose_mu_dmax, 1.007, delta=0.001)

        # MDA TB1 2015
        tg51_10x = tg51.TG51Photon(temp=21, press=763, model='30013', n_dw=5.393,
                                   p_elec=1.003, m_raw=27.727, m_opp=-27.784, m_low=27.635,
                                   measured_pdd=73.42, clinical_pdd=73.5)
        self.assertAlmostEqual(tg51_10x.dose_mu_10, 0.734, delta=0.001)
        self.assertAlmostEqual(tg51_10x.dose_mu_dmax, 0.998, delta=0.001)

        # ACB5 2011 pre-adjust 9/4
        tg51_6x = tg51.TG51Photon(temp=22, press=751.2, model='30013', n_dw=5.450,
                                  m_raw=(24.82,), m_opp=-24.83, m_low=24.79, measured_pdd=66.8,
                                  clinical_pdd=66.8, tissue_correction=0.99)
        self.assertAlmostEqual(tg51_6x.dose_mu_10, 0.672, delta=0.001)
        self.assertAlmostEqual(tg51_6x.dose_mu_dmax, 1.0066, delta=0.001)

        # ACB5 2012
        tg51_6x = tg51.TG51Photon(temp=21.7, press=757.2, model='30013', n_dw=5.446,
                                  m_raw=25.27, m_opp=-25.19, m_low=25.17, measured_pdd=66.8,
                                  clinical_pdd=66.8, tissue_correction=0.99)
        self.assertAlmostEqual(tg51_6x.dose_mu_10, 0.679, delta=0.001)
        self.assertAlmostEqual(tg51_6x.dose_mu_dmax, 1.0159, delta=0.001)

        # ACB5 2012
        tg51_18x = tg51.TG51Photon(temp=21.7, press=757.2, model='30013', n_dw=5.446,
                                  m_raw=30.67, m_opp=-30.65, m_low=30.50, energy=18, measured_pdd=79.5,
                                  clinical_pdd=79.7, tissue_correction=0.99, lead_foil=None)
        self.assertAlmostEqual(tg51_18x.dose_mu_10, 0.8059, delta=0.001)
        self.assertAlmostEqual(tg51_18x.dose_mu_dmax, 1.011, delta=0.001)

        # IMMC TB 15x
        tg51_15x = tg51.TG51Photon(temp=22.4, press=748.1, model='30013', n_dw=5.394,
                                  m_raw=14.307, m_opp=-14.323, m_low=14.22, energy=15, measured_pdd=76.79,
                                  clinical_pdd=76.7, lead_foil='30cm', volt_high=-300,
                                   volt_low=-150, mu=100)
        self.assertAlmostEqual(tg51_15x.dose_mu_10, 0.769, delta=0.001)
        self.assertAlmostEqual(tg51_15x.dose_mu_dmax, 1.002, delta=0.001)

    def test_tg51_electron(self):
        # ACB5 2011 pre-adjust 9/4
        tg51_9e = tg51.TG51Electron(temp=21.6, press=751.9, model='30013', n_dw=5.450,
                                     m_raw=39.79, m_opp=-39.71, m_low=39.33, i_50=3.87,
                                     clinical_pdd=100, k_ecal=0.897, m_plus=39.79, tissue_correction=0.99)
        self.assertAlmostEqual(tg51_9e.dose_mu_dref, 0.997, delta=0.001)
        self.assertAlmostEqual(tg51_9e.dose_mu_dmax, 0.998, delta=0.001)

        # ACB5 2012
        tg51_16e = tg51.TG51Electron(temp=21.5, press=758, model='30013', n_dw=5.446,
                                  m_raw=(40.71,), m_opp=-40.71, m_low=40.22, i_50=6.42,
                                  clinical_pdd=99.5, k_ecal=0.897, m_plus=40.71, tissue_correction=0.99)
        self.assertAlmostEqual(tg51_16e.dose_mu_dref, 1.000, delta=0.001)
        self.assertAlmostEqual(tg51_16e.dose_mu_dmax, 1.005, delta=0.001)