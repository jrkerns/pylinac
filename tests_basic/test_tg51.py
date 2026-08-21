import copy
from dataclasses import FrozenInstanceError
from unittest import TestCase

from argue import BoundsError
from parameterized import parameterized

from pylinac.calibration import tg51
from tests_basic.utils import save_file


class TestFunctions(TestCase):
    @parameterized.expand([(22, 101.33, 1.0), (25, 102.66, 0.997), (19, 98.66, 1.0165)])
    def test_p_tp(self, temperature, pressure, ptp):
        self.assertAlmostEqual(
            tg51.p_tp(temp=temperature, press=pressure), ptp, delta=0.001
        )

    def test_override_p_tp(self):
        original_temp = copy.copy(tg51.MAX_TEMP)
        tg51.MAX_TEMP = 20
        temp = 22
        press = 101.33
        with self.assertRaises(BoundsError):
            tg51.p_tp(temp=temp, press=press)
        tg51.MAX_TEMP = original_temp  # set back so other tests don't fail

    @parameterized.expand([(20, -20, 1.0), (-20.2, 19.8, 0.99), (19.8, -20.1, 1.0075)])
    def test_p_pol(self, m_ref, m_opposite, ppol):
        self.assertAlmostEqual(
            tg51.p_pol(m_reference=m_ref, m_opposite=m_opposite), ppol, delta=0.001
        )

    @parameterized.expand([(20, 20, 1.0), (20.05, 20.1, 1.0025)])
    def test_p_ion(self, m_low, m_high, pion):
        self.assertAlmostEqual(
            tg51.p_ion(
                voltage_reference=300,
                voltage_reduced=150,
                m_reference=m_high,
                m_reduced=m_low,
            ),
            pion,
            delta=0.001,
        )

    @parameterized.expand([(3, 1.72), (5, 2.96), (7, 4.19)])
    def test_dref(self, i50, dref):
        self.assertAlmostEqual(tg51.d_ref(i_50=i50), dref, delta=0.01)

    @parameterized.expand([(3.5, 3.54), (5.5, 5.60), (12, 18.71)])
    def test_r50(self, i50, r50):
        self.assertAlmostEqual(tg51.r_50(i_50=i50), r50, delta=0.01)

    @parameterized.expand(
        [(0.38 / 0.663, 0.6662), (0.385 / 0.667, 0.6713), (0.527 / 0.793, 0.7819)]
    )
    def test_pdd_to_tpr(self, pdd, tpr):
        self.assertAlmostEqual(tg51.tpr2010_from_pdd2010(pdd2010=pdd), tpr, delta=0.01)

    def test_m_corr(self):
        exp = 20.225
        res = tg51.m_corrected(
            p_ion=1.01, p_tp=0.995, p_elec=1, p_pol=1.005, m_reference=(20, 20.05)
        )
        self.assertAlmostEqual(exp, res, delta=0.002)

    @parameterized.expand(
        [
            (66.4, 6, 66.4, None),
            (70.5, 10, 70.5, "30cm"),
            (72.8, 10, 72.8, "50cm"),
            (73.3, 10, 73.3, None),
            (76.7, 15, 77.18, None),
            (77.1, 15, 77.57, "50cm"),
            (77.1, 15, 78.27, "30cm"),
            (79.3, 18, 80.47, None),
        ]
    )
    def test_pddx(self, pdd, energy, pddx, foil):
        self.assertAlmostEqual(
            tg51.pddx(pdd=pdd, energy=energy, lead_foil=foil), pddx, delta=0.01
        )

    def test_pddx_accepts_lead_foil_enum(self):
        enum_result = tg51.pddx(pdd=77.1, energy=15, lead_foil=tg51.LeadFoil.CM30)
        string_result = tg51.pddx(pdd=77.1, energy=15, lead_foil="30cm")
        self.assertEqual(enum_result, string_result)

    @parameterized.expand([("30010", 66.4, 0.9927), ("A12", 76.7, 0.976)])
    def test_kq_photon_pdd(self, chamber, pddx, kq):
        self.assertAlmostEqual(
            tg51.kq_photon_pddx(chamber=chamber, pddx=pddx), kq, delta=0.001
        )

    @parameterized.expand([("30010", 0.666, 0.9927)])
    def test_kq_photon_tpr(self, chamber, tpr, kq):
        self.assertAlmostEqual(
            tg51.kq_photon_tpr(chamber=chamber, tpr=tpr), kq, delta=0.001
        )

    @parameterized.expand([("30010", 3, 0.926), ("A12", 5, 0.915)])
    def test_kq_electron(self, chamber, r_50, kq):
        # Test via PDDs
        self.assertAlmostEqual(
            tg51.kq_electron(chamber=chamber, r_50=r_50), kq, delta=0.001
        )


class TestIonChamber(TestCase):
    def test_registry_values_match_chamber_table_keys(self):
        chamber_names = {chamber.name for chamber in tg51.IonChambers()}
        table_keys = set(tg51.KQ_PHOTONS) | set(tg51.KQ_ELECTRONS)
        self.assertEqual(chamber_names, table_keys)

    def test_iteration_returns_chambers(self):
        chambers = list(tg51.IonChambers())
        self.assertEqual(len(chambers), 34)
        self.assertTrue(
            all(isinstance(chamber, tg51.IonChamber) for chamber in chambers)
        )

    def test_from_name_returns_matching_chamber(self):
        chamber = tg51.IonChambers.from_name("30013")
        self.assertIs(chamber, tg51.IonChambers.PTW_30013)

    def test_from_name_builds_custom_photon_chamber_from_legacy_table(self):
        name = "Custom photon chamber"
        tg51.KQ_PHOTONS[name] = {
            "a": 1,
            "b": 0.001,
            "c": -0.00001,
            "a'": 2,
            "b'": -3,
            "c'": 4,
            "d'": -2,
        }
        self.addCleanup(tg51.KQ_PHOTONS.pop, name, None)

        chamber = tg51.IonChambers.from_name(name)

        self.assertEqual(chamber.name, name)
        self.assertIsNone(chamber.electron_coefficients)
        self.assertEqual(
            tg51.kq_photon_pddx(chamber=name, pddx=70),
            1 + 0.001 * 70 - 0.00001 * 70**2,
        )

    def test_from_name_builds_custom_electron_chamber_from_legacy_table(self):
        name = "Custom electron chamber"
        tg51.KQ_ELECTRONS[name] = {
            "kQ,ecal": 0.9,
            "a": 0.95,
            "b": 0.12,
            "c": 0.5,
        }
        self.addCleanup(tg51.KQ_ELECTRONS.pop, name, None)

        chamber = tg51.IonChambers.from_name(name)

        self.assertEqual(chamber.name, name)
        self.assertIsNone(chamber.photon_coefficients)
        self.assertEqual(
            tg51.kq_electron(chamber=name, r_50=3),
            (0.95 + 0.12 * 3**-0.5) * 0.9,
        )

    def test_from_name_combines_custom_photon_and_electron_coefficients(self):
        name = "Custom dual-energy chamber"
        tg51.KQ_PHOTONS[name] = {
            "a": 1,
            "b": 0.001,
            "c": -0.00001,
            "a'": 2,
            "b'": -3,
            "c'": 4,
            "d'": -2,
        }
        tg51.KQ_ELECTRONS[name] = {
            "kQ,ecal": 0.9,
            "a": 0.95,
            "b": 0.12,
            "c": 0.5,
        }
        self.addCleanup(tg51.KQ_PHOTONS.pop, name, None)
        self.addCleanup(tg51.KQ_ELECTRONS.pop, name, None)

        chamber = tg51.IonChambers.from_name(name)

        self.assertIsNotNone(chamber.photon_coefficients)
        self.assertIsNotNone(chamber.electron_coefficients)

    def test_from_name_identifies_missing_custom_coefficient_key(self):
        name = "Incomplete custom chamber"
        tg51.KQ_PHOTONS[name] = {"a": 1}
        self.addCleanup(tg51.KQ_PHOTONS.pop, name, None)

        with self.assertRaisesRegex(
            ValueError,
            "Photon coefficients for ion chamber 'Incomplete custom chamber' "
            "are missing required key 'b'",
        ):
            tg51.IonChambers.from_name(name)

    def test_from_name_error_is_clear(self):
        with self.assertRaisesRegex(
            ValueError,
            "No matching ion chamber exists for the name 'not-a-chamber'",
        ):
            tg51.IonChambers.from_name("not-a-chamber")

    def test_lookup_tables_retain_string_keys(self):
        self.assertTrue(all(isinstance(key, str) for key in tg51.KQ_PHOTONS))
        self.assertTrue(all(isinstance(key, str) for key in tg51.KQ_ELECTRONS))
        self.assertEqual(tg51.KQ_PHOTONS["30013"]["a"], 0.9652)
        self.assertEqual(tg51.KQ_ELECTRONS["30013"]["kQ,ecal"], 0.901)

    def test_chamber_and_coefficients_are_frozen(self):
        chamber = tg51.IonChambers.PTW_30013
        with self.assertRaises(FrozenInstanceError):
            chamber.name = "changed"
        with self.assertRaises(FrozenInstanceError):
            chamber.photon_coefficients.a = 2

    def test_lookup_tables_project_dataclass_coefficients(self):
        chamber = tg51.IonChambers.PTW_30013
        self.assertEqual(tg51.KQ_PHOTONS["30013"]["a"], chamber.photon_coefficients.a)
        self.assertEqual(
            tg51.KQ_ELECTRONS["30013"]["kQ,ecal"],
            chamber.electron_coefficients.kq_ecal,
        )

    def test_custom_chamber_can_be_used(self):
        coefficients = tg51.PhotonChamberCoefficients(
            a=1,
            b=0.001,
            c=-0.00001,
            a_prime=2,
            b_prime=-3,
            c_prime=4,
            d_prime=-2,
        )
        chamber = tg51.IonChamber(
            name="Custom",
            photon_coefficients=coefficients,
            electron_coefficients=None,
        )
        result = tg51.kq_photon_pddx(chamber=chamber, pddx=70)
        self.assertEqual(
            result, coefficients.a + coefficients.b * 70 + coefficients.c * 70**2
        )

    def test_kq_photon_pddx_accepts_chamber(self):
        chamber_result = tg51.kq_photon_pddx(
            chamber=tg51.IonChambers.PTW_30010, pddx=66.4
        )
        string_result = tg51.kq_photon_pddx(chamber="30010", pddx=66.4)
        self.assertEqual(chamber_result, string_result)

    def test_kq_photon_tpr_accepts_chamber(self):
        chamber_result = tg51.kq_photon_tpr(
            chamber=tg51.IonChambers.PTW_30010, tpr=0.666
        )
        string_result = tg51.kq_photon_tpr(chamber="30010", tpr=0.666)
        self.assertEqual(chamber_result, string_result)

    def test_kq_electron_accepts_chamber(self):
        chamber_result = tg51.kq_electron(chamber=tg51.IonChambers.PTW_30010, r_50=3)
        string_result = tg51.kq_electron(chamber="30010", r_50=3)
        self.assertEqual(chamber_result, string_result)

    def test_pddx_rejects_electron_only_chamber(self):
        with self.assertRaisesRegex(ValueError, "does not have photon coefficients"):
            tg51.kq_photon_pddx(chamber=tg51.IonChambers.PTW_31013, pddx=66.4)

    def test_kq_rejects_photon_only_chamber(self):
        with self.assertRaisesRegex(ValueError, "does not have electron coefficients"):
            tg51.kq_electron(chamber=tg51.IonChambers.EXRADIN_A2, r_50=3)

    def test_invalid_string_raises_clear_error(self):
        with self.assertRaisesRegex(
            ValueError,
            "No matching ion chamber exists for the name 'not-a-chamber'",
        ):
            tg51.kq_photon_pddx(chamber="not-a-chamber", pddx=66.4)

    @staticmethod
    def _common_class_arguments() -> dict:
        return {
            "unit": "Test Linac",
            "temp": 22,
            "press": 101.33,
            "chamber": "30013",
            "n_dw": 5.4,
            "p_elec": 1.0,
            "energy": 6,
            "voltage_reference": -300,
            "voltage_reduced": -150,
            "m_reference": 20,
            "m_opposite": -20,
            "m_reduced": 19.9,
            "mu": 100,
        }

    def test_photon_class_resolves_string_to_chamber(self):
        calibration = tg51.TG51Photon(
            **self._common_class_arguments(),
            measured_pdd10=66,
            clinical_pdd10=66,
        )
        self.assertIs(calibration.chamber, tg51.IonChambers.PTW_30013)

    def test_modern_electron_class_resolves_string_to_chamber(self):
        calibration = tg51.TG51ElectronModern(
            **self._common_class_arguments(),
            clinical_pdd=100,
            cone="10x10",
            i_50=3,
            tissue_correction=1,
        )
        self.assertIs(calibration.chamber, tg51.IonChambers.PTW_30013)

    def test_legacy_electron_class_resolves_known_string_to_chamber(self):
        calibration = tg51.TG51ElectronLegacy(
            **self._common_class_arguments(),
            clinical_pdd=100,
            cone="10x10",
            i_50=3,
            tissue_correction=1,
            k_ecal=0.9,
            m_gradient=20,
        )
        self.assertIs(calibration.chamber, tg51.IonChambers.PTW_30013)

    def test_legacy_electron_class_still_accepts_arbitrary_string(self):
        arguments = self._common_class_arguments()
        arguments["chamber"] = "custom chamber"
        calibration = tg51.TG51ElectronLegacy(
            **arguments,
            clinical_pdd=100,
            cone="10x10",
            i_50=3,
            tissue_correction=1,
            k_ecal=0.9,
            m_gradient=20,
        )
        self.assertEqual(calibration.chamber, "custom chamber")


class TG51TestBase:
    temperature = 22
    pressure = 760
    chamber = "30013"
    unit = "TBtest"
    nd_w = 5.555
    p_elec = 1.000
    volt_high = -300
    volt_low = -150
    clinical_pdd = 66
    dose_mu_dmax = 1.000
    dose_mu_dmax_adjusted = None
    tissue_correction = 1.000
    mu = 200
    print_data = False
    open_pdf = False
    m_reference: float
    m_reduced: float
    m_adjusted: float = None
    tg51: tg51.TG51Base

    def test_dose_dmax(self):
        self.assertAlmostEqual(self.dose_mu_dmax, self.tg51.dose_mu_dmax, delta=0.0005)

    def test_dose_dmax_adjusted(self):
        if self.m_adjusted is not None:
            self.assertAlmostEqual(
                self.dose_mu_dmax_adjusted,
                self.tg51.dose_mu_dmax_adjusted,
                delta=0.0005,
            )

    def test_pdf(self):
        save_file(self.tg51.publish_pdf)
        if self.open_pdf:
            self.tg51.publish_pdf("testtg51.pdf", open_file=True)


class TG51PhotonTestBase(TG51TestBase):
    energy = 6
    measured_pdd = 66
    lead_foil = None
    dose_mu_10 = 1.000
    fff = False
    tg51 = tg51.TG51Photon

    def setUp(self):
        self.tg51 = tg51.TG51Photon(
            temp=self.temperature,
            press=self.pressure,
            unit=self.unit,
            chamber=self.chamber,
            n_dw=self.nd_w,
            p_elec=self.p_elec,
            measured_pdd10=self.measured_pdd10,
            lead_foil=self.lead_foil,
            clinical_pdd10=self.clinical_pdd10,
            energy=self.energy,
            voltage_reference=self.volt_high,
            voltage_reduced=self.volt_low,
            m_reference=self.m_reference,
            m_opposite=self.m_opposite,
            m_reduced=self.m_reduced,
            m_reference_adjusted=self.m_adjusted,
            mu=self.mu,
            tissue_correction=self.tissue_correction,
            fff=self.fff,
        )
        if self.print_data:
            self.print_results()

    def test_dose_10(self):
        self.assertAlmostEqual(self.dose_mu_10, self.tg51.dose_mu_10, delta=0.0005)

    def test_lead_foil_is_enum(self):
        self.assertIsInstance(self.tg51.lead_foil, tg51.LeadFoil)

    def print_results(self):
        print("kQ determined", self.tg51.kq)
        print("Pion", self.tg51.p_ion)
        print("Ppol", self.tg51.p_pol)
        print("Ptp", self.tg51.p_tp)


class TG51ElectronLegacyTestBase(TG51TestBase):
    k_ecal = None
    i_50 = 7.5
    dose_mu_dref = 1.000
    tissue_correction = 1.0
    m_gradient = 0
    m_opposite: float
    energy = 0
    cone = ""
    tg51 = tg51.TG51ElectronLegacy

    def setUp(self):
        self.tg51 = tg51.TG51ElectronLegacy(
            temp=self.temperature,
            press=self.pressure,
            chamber=self.chamber,
            n_dw=self.nd_w,
            p_elec=self.p_elec,
            clinical_pdd=self.clinical_pdd,
            voltage_reference=self.volt_high,
            voltage_reduced=self.volt_low,
            m_reference=self.m_reference,
            m_opposite=self.m_opposite,
            m_reduced=self.m_reduced,
            m_reference_adjusted=self.m_adjusted,
            mu=self.mu,
            tissue_correction=self.tissue_correction,
            i_50=self.i_50,
            k_ecal=self.k_ecal,
            m_gradient=self.m_gradient,
            energy=self.energy,
            cone=self.cone,
        )

    def test_dose_dref(self):
        self.assertAlmostEqual(self.dose_mu_dref, self.tg51.dose_mu_dref, delta=0.0005)


class TG51ElectronModernTestBase(TG51TestBase):
    i_50 = 7.5
    dose_mu_dref = 1.000
    tissue_correction = 1.0
    energy = 0
    cone = ""
    tg51 = tg51.TG51ElectronModern

    def setUp(self):
        self.tg51 = tg51.TG51ElectronModern(
            temp=self.temperature,
            press=self.pressure,
            chamber=self.chamber,
            n_dw=self.nd_w,
            p_elec=self.p_elec,
            clinical_pdd=self.clinical_pdd,
            voltage_reference=self.volt_high,
            voltage_reduced=self.volt_low,
            m_reference=self.m_reference,
            m_opposite=self.m_reduced,
            m_reduced=self.m_reduced,
            m_reference_adjusted=self.m_adjusted,
            mu=self.mu,
            tissue_correction=self.tissue_correction,
            i_50=self.i_50,
            energy=self.energy,
            cone=self.cone,
        )

    def test_dose_dref(self):
        self.assertAlmostEqual(self.dose_mu_dref, self.tg51.dose_mu_dref, delta=0.0005)


class MDA_TB2_2015_15x(TG51PhotonTestBase, TestCase):
    energy = 15
    temperature = 20.5
    pressure = tg51.mmHg2kPa(760)
    nd_w = 5.444
    p_elec = 1.002
    m_adjusted = 29.28
    m_reference = 29.28
    m_opposite = -29.33
    m_reduced = 29.10
    measured_pdd10 = 76.9
    clinical_pdd10 = 77.4
    dose_mu_10 = 0.779
    dose_mu_dmax = 1.007
    dose_mu_dmax_adjusted = 1.007


class MDA_TB1_2015_10x(TG51PhotonTestBase, TestCase):
    energy = 10
    temperature = 21
    pressure = tg51.mmHg2kPa(763)
    nd_w = 5.393
    p_elec = 1.003
    m_reference = 27.727
    m_opposite = 27.784
    m_reduced = 27.635
    measured_pdd10 = 73.42
    clinical_pdd10 = 73.5
    dose_mu_10 = 0.733
    dose_mu_dmax = 0.998
    # print_data = True


class ACB5_2011_6x(TG51PhotonTestBase, TestCase):
    temperature = 22
    pressure = tg51.mmHg2kPa(751.2)
    nd_w = 5.450
    m_reference = 24.82
    m_opposite = -24.83
    m_reduced = 24.79
    measured_pdd10 = 66.8
    clinical_pdd10 = 66.8
    tissue_correction = 0.99
    dose_mu_10 = 0.672
    dose_mu_dmax = 1.0064


class ACB5_2012_6X(TG51PhotonTestBase, TestCase):
    temperature = 21.7
    pressure = tg51.mmHg2kPa(757.2)
    nd_w = 5.446
    m_reference = 25.27
    m_opposite = -25.19
    m_reduced = 25.17
    measured_pdd10 = 66.8
    clinical_pdd10 = 66.8
    tissue_correction = 0.99
    dose_mu_10 = 0.679
    dose_mu_dmax = 1.0159


class ACB5_2012_18X(TG51PhotonTestBase, TestCase):
    temperature = 21.7
    pressure = tg51.mmHg2kPa(757.2)
    nd_w = 5.446
    m_reference = 30.67
    m_opposite = -30.65
    m_reduced = 30.50
    energy = 18
    measured_pdd10 = 79.5
    clinical_pdd10 = 79.7
    tissue_correction = 0.99
    lead_foil = None
    dose_mu_10 = 0.8059
    dose_mu_dmax = 1.011


class IMMCTB_6FFF(TG51PhotonTestBase, TestCase):
    energy = 6
    fff = True
    temperature = 22.5
    pressure = tg51.mmHg2kPa(749)
    nd_w = 5.394
    m_reference = 11.610
    m_opposite = -11.613
    m_reduced = 11.533
    measured_pdd10 = 64.16
    clinical_pdd10 = 63.5
    mu = 100
    dose_mu_10 = 0.637
    dose_mu_dmax = 1.0033
    # print_data = True


class IMMCTB_10FFF(TG51PhotonTestBase, TestCase):
    energy = 10
    fff = True
    temperature = 22.4
    pressure = tg51.mmHg2kPa(748.1)
    nd_w = 5.394
    m_reference = 13.00067
    m_opposite = -13.013
    m_reduced = 12.867
    measured_pdd10 = 71.386
    clinical_pdd10 = 71.1
    lead_foil = "30cm"
    mu = 100
    dose_mu_10 = 0.710
    dose_mu_dmax = 0.9985
    # open_pdf = True


class IMMCTB_15X(TG51PhotonTestBase, TestCase):
    energy = 15
    temperature = 22.4
    pressure = tg51.mmHg2kPa(748.1)
    nd_w = 5.394
    m_reference = 14.307
    m_opposite = -14.323
    m_reduced = 14.220
    measured_pdd10 = 76.79
    clinical_pdd10 = 76.7
    mu = 100
    dose_mu_10 = 0.770
    dose_mu_dmax = 1.0036
    # print_data = True
    # open_pdf = True


class IMMC_TB_6E(TG51ElectronLegacyTestBase, TestCase):
    energy = 6
    cone = "15x15"
    mu = 100
    temperature = 22
    pressure = tg51.mmHg2kPa(748.2)
    k_ecal = 0.897
    p_elec = 0.999
    nd_w = 5.394
    m_reference = 19.730
    m_opposite = 19.797
    m_reduced = 19.497
    m_gradient = 19.710
    i_50 = 2.35
    clinical_pdd = 100
    tissue_correction = 1.0
    dose_mu_dref = 1.0085
    dose_mu_dmax = 1.0085
    # open_pdf = True


class IMMC_TB_9E(TG51ElectronLegacyTestBase, TestCase):
    energy = 9
    cone = "15x15"
    mu = 100
    temperature = 22
    pressure = tg51.mmHg2kPa(748.2)
    k_ecal = 0.897
    p_elec = 0.999
    nd_w = 5.394
    m_reference = 19.877
    m_opposite = 19.933
    m_reduced = 19.643
    m_gradient = 19.877
    i_50 = 3.55
    clinical_pdd = 100
    tissue_correction = 1.0
    dose_mu_dref = 1.006
    dose_mu_dmax = 1.006


class IMMC_TB_12E(TG51ElectronLegacyTestBase, TestCase):
    energy = 12
    cone = "15x15"
    mu = 100
    temperature = 22.1
    pressure = tg51.mmHg2kPa(748.2)
    k_ecal = 0.897
    p_elec = 0.999
    nd_w = 5.394
    m_reference = 20.080
    m_opposite = 20.143
    m_reduced = 19.850
    m_gradient = 20.047
    i_50 = 4.96
    clinical_pdd = 99.9
    tissue_correction = 1.0
    dose_mu_dref = 1.006
    dose_mu_dmax = 1.0068


class IMMC_TB_20E(TG51ElectronLegacyTestBase, TestCase):
    energy = 20
    cone = "15x15"
    mu = 100
    temperature = 22.1
    pressure = tg51.mmHg2kPa(748.2)
    k_ecal = 0.897
    p_elec = 0.999
    nd_w = 5.394
    m_reference = 19.670
    m_opposite = 19.707
    m_reduced = 19.437
    i_50 = 8.22
    m_gradient = 19.543
    clinical_pdd = 96.8
    dose_mu_dref = 0.970
    dose_mu_dmax = 1.002
    # open_pdf = True


class IMMC_TB_20E_Modern(TG51ElectronModernTestBase, TestCase):
    energy = 20
    cone = "15x15"
    mu = 100
    temperature = 22.1
    pressure = tg51.mmHg2kPa(748.2)
    p_elec = 0.999
    nd_w = 5.394
    m_reference = 19.670
    m_opposite = 19.707
    m_reduced = 19.437
    i_50 = 8.22
    clinical_pdd = 96.8
    dose_mu_dref = 0.974
    dose_mu_dmax = 1.006
