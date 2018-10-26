from unittest import TestCase

from pylinac.calibration import trs398
from tests_basic.utils import save_file


class TestFunctions(TestCase):

    def test_k_s(self):
        low_vals = (20, 20.05)
        high_vals = (20, 20.1)
        expected_pion = (1.0, 1.0025)
        for low, high, exp in zip(low_vals, high_vals, expected_pion):
            self.assertAlmostEqual(trs398.k_s(voltage_reference=300, voltage_reduced=150, m_reference=high, m_reduced=low), exp,
                                   delta=0.001)

    def test_m_corrected(self):
        exp = 20.225
        res = trs398.m_corrected(k_s=1.01, k_tp=0.995, k_elec=1, k_pol=1.005, m_reference=(20, 20.05))
        self.assertAlmostEqual(exp, res, delta=0.002)

    def test_kq_photon(self):
        models = ('30010', 'A12')
        tprs = (0.65, 0.75)
        kqs = (0.994, 0.983)
        for model, tpr, kq in zip(models, tprs, kqs):
            self.assertAlmostEqual(trs398.kq_photon(chamber=model, tpr=tpr), kq, delta=0.001)

    def test_kq_electron(self):
        models = ('30013', '2571')
        r_50s = (4.5, 8.2)
        kqs = (0.909, 0.905)
        for model, r_50, kq in zip(models, r_50s, kqs):
            self.assertAlmostEqual(trs398.kq_electron(chamber=model, r_50=r_50), kq, delta=0.001)


class TestTRS398Base:
    temp = float
    press = float
    chamber = '30013'
    nd_w = float
    k_elec = 1.000
    voltage_reference = -300
    voltage_reduced = -150
    m_reference = tuple
    m_opposite = tuple
    m_reduced = tuple
    dose_mu_zmax = 1.000
    tissue_correction = 1.000
    mu = 200
    open_pdf = False
    print_data = False

    def test_dose_zmax(self):
        self.assertAlmostEqual(self.dose_mu_zmax, self.trs398.dose_mu_zmax, delta=0.0005)

    def test_dose_zref(self):
        self.assertAlmostEqual(self.dose_mu_zref, self.trs398.dose_mu_zref, delta=0.0005)

    def test_pdf(self):
        save_file(self.trs398.publish_pdf)
        if self.open_pdf:
            self.trs398.publish_pdf('testtrs.pdf', open_file=True)

    def print_results(self):
        print('kQ determined', self.trs398.kq)
        print('Pion', self.trs398.k_s)
        print('Ppol', self.trs398.k_pol)
        print('Ptp', self.trs398.k_tp)


class TestTRS398Photon(TestTRS398Base):
    energy = 6
    setup_condition = 'SSD'
    # clinical_pdd_zref = None
    clinical_tmr_zref = None
    fff = False

    def setUp(self):
        self.trs398 = trs398.TRS398Photon(
            setup=self.setup_condition,
            tpr2010=self.tpr2010,
            temp=self.temp, press=self.press,
            chamber=self.chamber, n_dw=self.nd_w, k_elec=self.k_elec,
            clinical_pdd_zref=self.clinical_pdd_zref,
            voltage_reference=self.voltage_reference, voltage_reduced=self.voltage_reduced,
            m_reference=self.m_reference, m_opposite=self.m_opposite, m_reduced=self.m_reduced,
            clinical_tmr_zref=self.clinical_tmr_zref,
            mu=self.mu, tissue_correction=self.tissue_correction, fff=self.fff, energy=self.energy
        )
        if self.print_data:
            self.print_results()


class TestTRS398Electron(TestTRS398Base):
    energy = None
    k_ecal = None
    i_50 = 7.5
    cone = '15x15'
    dose_mu_10 = 1.000

    def setUp(self):
        self.trs398 = trs398.TRS398Electron(temp=self.temp, press=self.press, energy=self.energy,
                                    chamber=self.chamber, n_dw=self.nd_w, k_elec=self.k_elec,
                                    clinical_pdd_zref=self.clinical_pdd_zref, i_50=self.i_50,
                                    voltage_reference=self.voltage_reference, voltage_reduced=self.voltage_reduced,
                                    m_reference=self.m_reference, m_opposite=self.m_opposite, m_reduced=self.m_reduced,
                                    mu=self.mu, tissue_correction=self.tissue_correction, cone=self.cone)


class MDA_TB2_2015_15x(TestTRS398Photon, TestCase):
    energy = 15
    temp = 20.5
    press = trs398.mmHg2kPa(760)
    nd_w = 5.444
    k_elec = 1.002
    m_reference = 29.28
    m_opposite = -29.33
    m_reduced = 29.10
    dose_mu_zref = 0.779
    dose_mu_zmax = 1.007
    clinical_pdd_zref = 77.4
    tpr2010 = 0.762


class MDA_TB1_2015_10x(TestTRS398Photon, TestCase):
    energy = 10
    temp = 21
    press = trs398.mmHg2kPa(763)
    nd_w = 5.393
    k_elec = 1.003
    m_reference = 27.727
    m_opposite = 27.784
    m_reduced = 27.635
    clinical_pdd_zref = 73.5
    dose_mu_zref = 0.734
    dose_mu_zmax = 0.998
    tpr2010 = (73.42/73.7) * trs398.tpr2010_from_pdd2010(pdd2010=46.3/73.7)
    # open_pdf = True
    # print_data = True


class ACB5_2011_6x(TestTRS398Photon, TestCase):
    temp = 22
    press = trs398.mmHg2kPa(751.2)
    nd_w = 5.450
    tpr2010 = trs398.tpr2010_from_pdd2010(pdd2010=38.4/66.8)
    m_reference = 24.82
    m_opposite = -24.83
    m_reduced = 24.79
    clinical_pdd_zref = 66.8
    tissue_correction = 0.99
    dose_mu_zref = 0.673
    dose_mu_zmax = 1.007

    def test_zmax_adjusted(self):
        self.trs398.m_reference_adjusted = 24.65
        self.assertAlmostEqual(self.trs398.dose_mu_zmax_adjusted, 1.000, delta=0.0005)

    def test_zref_adjusted(self):
        self.trs398.m_reference_adjusted = 24.65
        self.assertAlmostEqual(self.trs398.dose_mu_zref_adjusted, 0.668, delta=0.0005)


class ACB5_2012_6X(TestTRS398Photon, TestCase):
    temp = 21.7
    press = trs398.mmHg2kPa(757.2)
    nd_w = 5.446
    m_reference = 25.27
    m_opposite = -25.19
    m_reduced = 25.17
    clinical_pdd_zref = 66.8
    tpr2010 = trs398.tpr2010_from_pdd2010(pdd2010=38.4/66.8)
    tissue_correction = 0.99
    dose_mu_zref = 0.679
    dose_mu_zmax = 1.0159


class ACB5_2012_18X(TestTRS398Photon, TestCase):
    energy = 18
    temp = 21.7
    press = trs398.mmHg2kPa(757.2)
    tpr2010 = trs398.tpr2010_from_pdd2010(pdd2010=52.5/79.4)
    nd_w = 5.446
    m_reference = 30.67
    m_opposite = -30.65
    m_reduced = 30.50
    clinical_pdd_zref = 79.7
    tissue_correction = 0.99
    dose_mu_zref = 0.807
    dose_mu_zmax = 1.0125


class IMMCTB_6FFF(TestTRS398Photon, TestCase):
    energy = 6
    fff = True
    temp = 22.5
    press = trs398.mmHg2kPa(749)
    tpr2010 = (64.16 / 63.6) * trs398.tpr2010_from_pdd2010(pdd2010=34.5 / 63.6)
    nd_w = 5.394
    m_reference = 11.610
    m_opposite = -11.613
    m_reduced = 11.533
    clinical_pdd_zref = 63.5
    mu = 100
    dose_mu_zref = 0.638
    dose_mu_zmax = 1.005
    print_data = True


class IMMCTB_10FFF(TestTRS398Photon, TestCase):
    energy = 10
    fff = True
    temp = 22.4
    press = trs398.mmHg2kPa(748.1)
    nd_w = 5.394
    m_reference = 13.00067
    m_opposite = -13.013
    m_reduced = 12.867
    tpr2010 = trs398.tpr2010_from_pdd2010(pdd2010=(43/71.2))
    clinical_pdd_zref = 71.1
    mu = 100
    dose_mu_zref = 0.712
    dose_mu_zmax = 1.0005
    # open_pdf = True


class IMMCTB_15X(TestTRS398Photon, TestCase):
    energy = 15
    temp = 22.4
    press = trs398.mmHg2kPa(748.1)
    nd_w = 5.394
    m_reference = 14.307
    m_opposite = -14.323
    m_reduced = 14.220
    clinical_pdd_zref = 76.7
    tpr2010 = trs398.tpr2010_from_pdd2010(pdd2010=(49.9/76.9)) * (76.79/76.9)
    mu = 100
    dose_mu_zref = 0.770
    dose_mu_zmax = 1.004
    # print_data = True
    # open_pdf = True


class IMMC_TB_20E(TestTRS398Electron, TestCase):
    energy = 20
    cone = '15x15'
    mu = 100
    temp = 22.1
    press = trs398.mmHg2kPa(748.2)
    k_elec = 0.999
    nd_w = 5.394
    m_reference = 19.670 * 0.99354  # * Pgradient because TRS398 is done at dref+0.5cm
    m_opposite = 19.707 * 0.99354
    m_reduced = 19.437 * 0.99354
    i_50 = 8.22
    clinical_pdd_zref = 96.8
    dose_mu_zref = 0.972
    dose_mu_zmax = 1.004
