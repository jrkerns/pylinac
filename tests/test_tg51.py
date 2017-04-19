from unittest import TestCase

from pylinac import tg51


class TestFunctions(TestCase):

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


class TestTG51Base:
    temperature = 22
    pressure = 760
    model = '30013'
    nd_w = 5.555
    p_elec = 1.000
    volt_high = -300
    volt_low = -150
    m_raw = (20, 20, 20)
    m_opp = (20, 20, 20)
    m_low = (20, 20, 20)
    clinical_pdd = 66
    dose_mu_dmax = 1.000
    tissue_correction = 1.000
    mu = 200

    def test_dose_dmax(self):
        self.assertAlmostEqual(self.dose_mu_dmax, self.tg51.dose_mu_dmax, delta=0.0005)


class TestTG51Photon(TestTG51Base):
    energy = 6
    measured_pdd = 66
    lead_foil = None
    dose_mu_10 = 1.000

    def setUp(self):
        self.tg51 = tg51.TG51Photon(temp=self.temperature, press=self.pressure,
                                    model=self.model, n_dw=self.nd_w, p_elec=self.p_elec,
                                    measured_pdd=self.measured_pdd, lead_foil=self.lead_foil,
                                    clinical_pdd=self.clinical_pdd, energy=self.energy,
                                    volt_high=self.volt_high, volt_low=self.volt_low,
                                    m_raw=self.m_raw, m_opp=self.m_opp, m_low=self.m_low,
                                    mu=self.mu, tissue_correction=self.tissue_correction)

    def test_dose_10(self):
        self.assertAlmostEqual(self.dose_mu_10, self.tg51.dose_mu_10, delta=0.0005)


class TestTG51Electron(TestTG51Base):
    k_ecal = 1.000
    i_50 = 7.5
    dose_mu_dref = 1.000

    def setUp(self):
        self.tg51 = tg51.TG51Electron(temp=self.temperature, press=self.pressure,
                                      model=self.model, n_dw=self.nd_w, p_elec=self.p_elec,
                                      clinical_pdd=self.clinical_pdd,
                                      volt_high=self.volt_high, volt_low=self.volt_low,
                                      m_raw=self.m_raw, m_opp=self.m_opp, m_low=self.m_low,
                                      mu=self.mu, tissue_correction=self.tissue_correction,
                                      i_50=self.i_50)

    def test_dose_dref(self):
        self.assertAlmostEqual(self.dose_mu_dref, self.tg51.dose_mu_dref, delta=0.0005)


class MDATB2_2015(TestTG51Photon, TestCase):
    temperature = 20.5
    pressure = 760
    energy = 15
    nd_w = 5.444
    p_elec = 1.002
    m_raw = 29.28
    m_opp = -29.33
    m_low = 29.10
    measured_pdd = 76.9
    clinical_pdd = 77.4
    dose_mu_10 = 0.779
    dose_mu_dmax = 1.007


class MDATB1_2015(TestTG51Photon, TestCase):
    temperature = 21
    pressure = 763
    nd_w = 5.393
    energy = 10
    p_elec = 1.003
    m_raw = 27.727
    m_opp = 27.784
    m_low = 27.635
    measured_pdd = 73.42
    clinical_pdd = 73.5
    dose_mu_10 = 0.734
    dose_mu_dmax = 0.999


class ACB5_2011(TestTG51Photon, TestCase):
    temperature = 22
    pressure = 751.2
    nd_w = 5.450
    m_raw = 24.82
    m_opp = -24.83
    m_low = 24.79
    measured_pdd = 66.8
    clinical_pdd = 66.8
    tissue_correction = 0.99
    dose_mu_10 = 0.672
    dose_mu_dmax = 1.0066


class ACB5_2012_6X(TestTG51Photon, TestCase):
    temperature = 21.7
    pressure = 757.2
    nd_w = 5.446
    m_raw = 25.27
    m_opp = -25.19
    m_low = 25.17
    measured_pdd = 66.8
    clinical_pdd = 66.8
    tissue_correction = 0.99
    dose_mu_10 = 0.679
    dose_mu_dmax = 1.0159


class ACB5_2012_18X(TestTG51Photon, TestCase):
    temperature = 21.7
    pressure = 757.2
    nd_w = 5.446
    m_raw = 30.67
    m_opp = -30.65
    m_low = 30.50
    energy = 18
    measured_pdd = 79.5
    clinical_pdd = 79.7
    tissue_correction = 0.99
    lead_foil = None
    dose_mu_10 = 0.8059
    dose_mu_dmax = 1.011


class IMMCTB_15X(TestTG51Photon, TestCase):
    temperature = 22.4
    pressure = 748.1
    nd_w = 5.394
    m_raw = 14.307
    m_opp = -14.323
    m_low = 14.22
    energy = 15
    measured_pdd = 76.79
    clinical_pdd = 76.7
    lead_foil = '30cm'
    mu = 100
    dose_mu_10 = 0.769
    dose_mu_dmax = 1.002


class ACB5_2011_9E(TestTG51Electron, TestCase):
    temperature = 21.6
    pressure = 751.9
    nd_w = 5.45
    m_raw = 39.79
    m_opp = -39.71
    m_low = 39.33
    i_50 = 3.87
    clinical_pdd = 100
    tissue_correction = 0.99
    dose_mu_dref = 0.997
    dose_mu_dmax = 0.997


class ACB5_2012_16E(TestTG51Electron, TestCase):
    temperature = 21.5
    pressure = 758
    nd_w = 5.446
    m_raw = 40.71
    m_opp = -40.71
    m_low = 40.22
    i_50 = 6.42
    clinical_pdd = 99.5
    k_ecal = 0.897
    m_plus = 40.71
    tissue_correction = 0.99
    dose_mu_dref = 1.000
    dose_mu_dmax = 1.005


class IMMC_TB_9E(TestTG51Electron, TestCase):
    mu = 100
    temperature = 22
    pressure = 748.2
    p_elec = 0.999
    nd_w = 5.394
    m_raw = 19.877
    m_opp = 19.933
    m_low = 19.643
    i_50 = 3.55
    clinical_pdd = 100
    tissue_correction = 1.0
    dose_mu_dref = 1.006
    dose_mu_dmax = 1.006