"""
The TG-51 module contains a number of helper functions and classes that can calculate parameters for performing the
TG-51 absolute linac dose calibration although there are some modifications from the original TG-51. The modifications 
include updated kQ and kecal values from Muir and Rodgers' set of papers.
Functions include all relevant calculations for TG-51 including PDDx, kQ,
Dref, and chamber reading corrections. Where Muir & Rodgers' values/equations are used they are specified in the documentation.

Classes include photon and electron calibrations using cylindrical chambers. Pass all the relevant raw measurements
and the class will compute all corrections and corrected readings and dose at 10cm and dmax/dref.
"""
from datetime import datetime
from typing import Optional

import argue
import numpy as np

from ..core.typing import NumberLike, NumberOrArray
from ..core.utilities import Structure, open_path
from ..core.pdf import PylinacCanvas


KQ_PHOTONS = {
    # Exradin
    'A12': {"a": 1.0146, "b": 0.777e-3, "c": -1.666e-5, "a'": 2.6402, "b'": -7.2304, "c'": 10.7573, "d'": -5.4294},
    'A19': {"a": 0.9934, "b": 1.384e-3, "c": -2.125e-5, "a'": 3.0907, "b'": -9.1930, "c'": 13.5957, "d'": -6.7969},
    'A2': {"a": 0.9819, "b": 1.609e-3, "c": -2.184e-5, "a'": 2.8458, "b'": -8.1619, "c'": 12.1411, "d'": -6.1041},
    'T2': {"a": 1.0173, "b": 0.854e-3, "c": -1.941e-5, "a'": 3.3433, "b'": -10.2649, "c'": 15.1247, "d'": -7.5415},
    'A12S': {"a": 0.9692, "b": 1.974e-3, "c": -2.448e-5, "a'": 2.9597, "b'": -8.6777, "c'": 12.9155, "d'": -6.4903},
    'A18': {"a": 0.9944, "b": 1.286e-3, "c": -1.980e-5, "a'": 2.5167, "b'": -6.7567, "c'": 10.1519, "d'": -5.1709},
    'A1': {"a": 1.0029, "b": 1.023e-3, "c": -1.803e-5, "a'": 2.0848, "b'": -4.9174, "c'": 7.5446, "d'": -3.9441},
    'T1': {"a": 1.0552, "b": -0.196e-3, "c": -1.275e-5, "a'": 2.8060, "b'": -7.9273, "c'": 11.7541, "d'": -5.9263},
    'A1SL': {"a": 0.9896, "b": 1.410e-3, "c": -2.049e-5, "a'": 2.8029, "b'": -7.9648, "c'": 11.8445, "d'": -5.9568},
    'A14': {"a": 0.9285, "b": 2.706e-3, "c": -2.599e-5, "a'": 5.4677, "b'": -19.1795, "c'": 27.4542, "d'": -13.1336},
    'T14': {"a": 0.9622, "b": 2.009e-3, "c": -2.401e-5, "a'": 4.9690, "b'": -17.1074, "c'": 24.6292, "d'": -11.8877},
    'A14SL': {"a": 0.9017, "b": 3.454e-3, "c": -3.083e-5, "a'": 5.1205, "b'": -17.7884, "c'": 25.6123, "d'": -12.3232},
    'A16': {"a": 0.8367, "b": 4.987e-3, "c": -3.877e-5, "a'": 6.0571, "b'": -21.7829, "c'": 31.2289, "d'": -14.9168},

    # PTW
    '30010': {"a": 1.0093, "b": 0.926e-3, "c": -1.771e-5, "a'": 2.5318, "b'": -6.7948, "c'": 10.1779, "d'": -5.1746},
    '30011': {"a": 0.9676, "b": 2.061e-3, "c": -2.528e-5, "a'": 2.9044, "b'": -8.4576, "c'": 12.6339, "d'": -6.3742},
    '30012': {"a": 0.9537, "b": 2.440e-3, "c": -2.750e-5, "a'": 3.2836, "b'": -10.0610, "c'": 14.8867, "d'": -7.4212},
    '30013': {"a": 0.9652, "b": 2.141e-3, "c": -2.623e-5, "a'": 3.2012, "b'": -9.7211, "c'": 14.4211, "d'": -7.2184},
    '31010': {"a": 0.9590, "b": 2.265e-3, "c": -2.684e-5, "a'": 3.1578, "b'": -9.5422, "c'": 14.1676, "d'": -7.0964},
    '31016': {"a": 1.0085, "b": 1.028e-3, "c": -1.968e-5, "a'": 2.9524, "b'": -8.6054, "c'": 12.7757, "d'": -6.4265},
    '31014': {"a": 1.0071, "b": 1.048e-3, "c": -1.967e-5, "a'": 3.0178, "b'": -8.8735, "c'": 13.1372, "d'": -6.5867},

    # IBA
    'CC25': {"a": 0.9551, "b": 2.353e-3, "c": -2.687e-5, "a'": 2.4567, "b'": -6.5932, "c'": 10.0471, "d'": -5.1775},
    'CC13': {"a": 0.9515, "b": 2.455e-3, "c": -2.768e-5, "a'": 3.1982, "b'": -9.7182, "c'": 14.4210, "d'": -7.2121},
    'CC08': {"a": 0.9430, "b": 2.637e-3, "c": -2.884e-5, "a'": 3.7328, "b'": -11.9800, "c'": 17.5884, "d'": -8.6843},
    'CC04': {"a": 0.9714, "b": 1.938e-3, "c": -2.432e-5, "a'": 3.0054, "b'": -8.8633, "c'": 13.1704, "d'": -6.6075},
    'CC01': {"a": 0.9116, "b": 3.358e-3, "c": -3.177e-5, "a'": 4.3376, "b'": -14.4935, "c'": 21.0293, "d'": -10.2208},
    'FC65-G': {"a": 0.9708, "b": 1.972e-3, "c": -2.480e-5, "a'": 3.3221, "b'": -10.2012, "c'": 15.0497, "d'": -7.4872},
    'FC65-P': {"a": 0.9828, "b": 1.664e-3, "c": -2.296e-5, "a'": 3.0872, "b'": -9.1919, "c'": 13.6137, "d'": -6.8118},
    'FC23-C': {"a": 0.9820, "b": 1.579e-3, "c": -2.166e-5, "a'": 3.0511, "b'": -9.0243, "c'": 13.3378, "d'": -6.6559},

    # Other
    'NE2581': {"a": 1.0318, "b": 0.488e-3, "c": -1.731e-5, "a'": 2.9190, "b'": -8.4561, "c'": 12.5690, "d'": -6.3468},
    'NE2571': {"a": 0.9882, "b": 1.486e-3, "c": -2.140e-5, "a'": 2.2328, "b'": -5.5779, "c'": 8.5325, "d'": -4.4352},
    'NE2561': {"a": 1.0200, "b": 0.596e-3, "c": -1.551e-5, "a'": 2.4235, "b'": -6.3179, "c'": 9.4737, "d'": -4.8307},
    'PR06C/G': {"a": 0.9519, "b": 2.432e-3, "c": -2.704e-5, "a'": 2.9110, "b'": -8.4916, "c'": 12.6817, "d'": -6.3874},
}

KQ_ELECTRONS = {
    # Exradin
    'A12': {'kQ,ecal': 0.907, 'a': 0.965, 'b': 0.119, 'c': 0.607},
    'A19': {'kQ,ecal': 0.904, 'a': 0.957, 'b': 0.119, 'c': 0.505},
    'A12S': {'kQ,ecal': 0.907, 'a': 0.937, 'b': 0.136, 'c': 0.378},
    'A18': {'kQ,ecal': 0.914, 'a': 0.352, 'b': 0.711, 'c': 0.046},
    'A1SL': {'kQ,ecal': 0.914, 'a': 0.205, 'b': 0.854, 'c': 0.036},

    # PTW
    '30010': {'kQ,ecal': 0.904, 'a': 0.980, 'b': 0.119, 'c': 0.891},
    '30011': {'kQ,ecal': 0.901, 'a': 0.976, 'b': 0.120, 'c': 0.793},
    '30012': {'kQ,ecal': 0.908, 'a': 0.972, 'b': 0.121, 'c': 0.728},
    '30013': {'kQ,ecal': 0.901, 'a': 0.978, 'b': 0.112, 'c': 0.816},
    '31013': {'kQ,ecal': 0.902, 'a': 0.945, 'b': 0.133, 'c': 0.441},

    # IBA
    'FC65-G': {'kQ,ecal': 0.904, 'a': 0.971, 'b': 0.113, 'c': 0.680},
    'FC65-P': {'kQ,ecal': 0.902, 'a': 0.973, 'b': 0.110, 'c': 0.692},
    'FC23-C': {'kQ,ecal': 0.904, 'a': 0.971, 'b': 0.097, 'c': 0.591},
    'CC25': {'kQ,ecal': 0.904, 'a': 0.964, 'b': 0.105, 'c': 0.539},
    'CC13': {'kQ,ecal': 0.904, 'a': 0.926, 'b': 0.129, 'c': 0.279},

    # Other
    'PR06C/G': {'kQ,ecal': 0.906, 'a': 0.972, 'b': 0.122, 'c': 0.729},
    'NE2571': {'kQ,ecal': 0.903, 'a': 0.977, 'b': 0.117, 'c': 0.817},
    'NE2611': {'kQ,ecal': 0.896, 'a': 0.979, 'b': 0.120, 'c': 0.875},
}

LEAD_OPTIONS = {
    'None': None,
    '30cm': '30cm',
    '50cm': '50cm'
}


def mmHg2kPa(mmHg: float) -> float:
    """Utility function to convert from mmHg to kPa."""
    return mmHg*101.33/760


def mbar2kPa(mbar: float) -> float:
    """Utility function to convert from millibars to kPa."""
    return mbar/10


def fahrenheit2celsius(f: float) -> float:
    """Utility function to convert from Fahrenheit to Celsius."""
    return (f - 32) * 5/9


@argue.bounds(pdd2010=(0.5, 1))
def tpr2010_from_pdd2010(*, pdd2010: float) -> float:
    """Calculate TPR20,10 from PDD20,10. From TRS-398 pg 62 and Followill et al 1998 eqn 1."""
    return 1.2661*pdd2010 - 0.0595


@argue.bounds(temp=(17, 27), message="Temperature {:2.2f} out of range. Did you use Fahrenheit? Consider using the utility function fahrenheit2celsius()")
@argue.bounds(press=(91,  111), message="Pressure {:2.2f} out of range. Did you use kPa? Consider using the utility functions mmHg2kPa() or mbar2kPa()")
def p_tp(*, temp: NumberLike, press: NumberLike) -> float:
    """Calculate the temperature & pressure correction.

    Parameters
    ----------
    temp : float (17-27)
        The temperature in degrees Celsius.
    press : float (91-111)
        The value of pressure in kPa. Can be converted from mmHg and mbar; see :func:`~pylinac.calibration.tg51.mmHg2kPa` and :func:`~pylinac.calibration.tg51.mbar2kPa`.
    """
    return ((273.2+temp)/295.2)*(101.33/press)


def p_pol(*, m_reference: NumberOrArray, m_opposite: NumberOrArray) -> float:
    """Calculate the polarity correction.

    Parameters
    ----------
    m_reference : number, array
        The readings of the ion chamber at the reference polarity and voltage.
    m_opposite : number, array
        The readings of the ion chamber at the polarity opposite the reference. The sign does not make a difference.

    Raises
    ------
    BoundsError if calculated Ppol is >1% from 1.0.
    """
    mref_avg = np.mean(m_reference)
    mopp_avg = np.mean(m_opposite)
    polarity = (abs(mref_avg) + abs(mopp_avg))/abs(2*mref_avg)
    argue.verify_bounds(polarity, bounds=(0.99, 1.01), message="Polarity correction {:2.2f} out of range (+/-2%). Verify inputs")
    return float(polarity)


def p_ion(*, voltage_reference: int, voltage_reduced: int, m_reference: NumberOrArray, m_reduced: NumberOrArray) -> float:
    """Calculate the ion chamber collection correction.

    Parameters
    ----------
    voltage_reference : int
        The "high" voltage; same as the TG51 measurement voltage.
    voltage_reduced : int
        The "low" voltage; usually half of the high voltage.
    m_reference : float, iterable
        The readings of the ion chamber at the "high" voltage.
    m_reduced : float, iterable
        The readings of the ion chamber at the "low" voltage.

    Raises
    ------
    BoundsError if calculated Pion is outside the range 1.00-1.05.
    """
    ion = (1 - voltage_reference / voltage_reduced) / (np.mean(m_reference) / np.mean(m_reduced) - voltage_reference / voltage_reduced)
    argue.verify_bounds(ion, bounds=(1, 1.05), message="Pion out of range (1.00-1.05). Check inputs or chamber")
    return float(ion)


@argue.bounds(i_50=argue.POSITIVE, message="i50 should be positive")
def d_ref(*, i_50: float) -> float:
    """Calculate the dref of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    r50 = r_50(i_50=i_50)
    return 0.6*r50-0.1


@argue.bounds(i_50=argue.POSITIVE, message="i50 should be positive")
def r_50(*, i_50: float) -> float:
    """Calculate the R50 depth of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    if i_50 < 10:
        r50 = 1.029 * i_50 - 0.06
    else:
        r50 = 1.59 * i_50 - 0.37
    return r50


@argue.bounds(r_50=(2, 9))
def kp_r50(*, r_50: float) -> float:
    """Calculate k'R50 for Farmer-like chambers.

    Parameters
    ----------
    r_50 : float (2-9)
        The R50 value in cm.
    """
    return 0.9905+0.071*np.exp(-r_50/3.67)


def pq_gr(*, m_dref_plus: NumberOrArray, m_dref: NumberOrArray) -> float:
    """Calculate PQ_gradient for a cylindrical chamber.

    Parameters
    ----------
    m_dref_plus : float, iterable
        The readings of the ion chamber at dref + 0.5rcav.
    m_dref : float, iterable
        The readings of the ion chamber at dref.
    """
    return float(np.mean(m_dref_plus) / np.mean(m_dref))


@argue.bounds(p_ion=(1, 1.05), p_tp=(0.92, 1.08), p_elec=(0.98, 1.02), p_pol=(0.98, 1.02))
def m_corrected(*, p_ion: float, p_tp: float, p_elec: float, p_pol: float, m_reference: NumberOrArray) -> float:
    """Calculate M_corrected, the ion chamber reading with all corrections applied.

    Parameters
    ----------
    p_ion : float (1.00-1.05)
        The ion collection correction.
    p_tp : float (0.92-1.08)
        The temperature & pressure correction.
    p_elec : float (0.98-1.02)
        The electrometer correction.
    p_pol : float (0.98-1.02)
        The polarity correction.
    m_reference : float, iterable
        The raw ion chamber reading(s).

    Returns
    -------
    float
    """
    return float(p_ion*p_tp*p_elec*p_pol*np.mean(m_reference))


@argue.bounds(pdd=(62.7, 89.0))
@argue.options(lead_foil=LEAD_OPTIONS.values())
def pddx(*, pdd: float, energy: int, lead_foil: Optional[str]=None) -> float:
    """Calculate PDDx based on the PDD.

    Parameters
    ----------
    pdd : {>62.7, <89.0}
        The measured PDD. If lead foil was used, this assumes the pdd as measured with the lead in place.
    energy : int
        The nominal energy in MV.
    lead_foil : {None, '30cm', '50cm'}
        Applicable only for energies >10MV.
        Whether a lead foil was used to acquire the pdd.
        Use ``None`` if no lead foil was used and the interim equation should be used. This is the default
        Use ``50cm`` if the lead foil was set to 50cm from the phantom surface.
        Use ``30cm`` if the lead foil was set to 30cm from the phantom surface.
    """
    if energy < 10:
        return pdd
    elif energy >= 10:
        if lead_foil is None:
            return 1.267*pdd-20
        elif lead_foil == LEAD_OPTIONS['50cm']:
            if pdd < 73:
                return pdd
            else:
                return (0.8905+0.0015*pdd)*pdd
        elif lead_foil == LEAD_OPTIONS['30cm']:
            if pdd < 71:
                return pdd
            else:
                return (0.8116+0.00264*pdd)*pdd


@argue.bounds(pddx=(63.0, 86.0))
@argue.options(chamber=KQ_PHOTONS.keys())
def kq_photon_pddx(*, chamber: str, pddx: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements of PDD(10)x. This will calculate kQ for photons
    for *CYLINDRICAL* chambers only.

    Parameters
    ----------
    chamber : str
        The chamber of the chamber. Valid values are those listed in
        Table III of Muir and Rodgers and Table I of the TG-51 Addendum.
    pddx : {>63.0, <86.0}
        The **PHOTON-ONLY** PDD measurement at 10cm depth for a 10x10cm2 field.

        .. note:: Use the :func:`~pylinac.calibration.tg51.pddx` function to convert PDD to PDDx as needed.

        .. note:: Muir and Rogers state limits of 0.627 - 0.861. The TG-51 addendum states them as 0.63 and 0.86.
                  The TG-51 addendum limits are used here.
    """
    ch = KQ_PHOTONS[chamber]
    return ch["a"] + ch["b"] * pddx + ch["c"] * (pddx ** 2)


@argue.bounds(tpr=(0.623, 0.805))
@argue.options(chamber=KQ_PHOTONS.keys())
def kq_photon_tpr(*, chamber: str, tpr: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements of TPR20,10. This will calculate kQ for photons
    for *CYLINDRICAL* chambers only.

    Parameters
    ----------
    chamber : str
        The chamber of the chamber. Valid values are those listed in
        Table III of Muir and Rodgers and Table I of the TG-51 Addendum.
    tpr : {>0.630, <0.860}
        The TPR(20,10) value.

        .. note::
         Use the :func:`~pylinac.calibration.tg51.tpr2010_from_pdd2010` function to convert from PDD without needing to take TPR measurements.
    """
    ch = KQ_PHOTONS[chamber]
    return ch["a'"] + ch["b'"] * tpr + ch["c'"] * (tpr ** 2) + ch["d'"] * (tpr ** 3)


@argue.options(chamber=KQ_ELECTRONS.keys())
def kq_electron(*, chamber: str, r_50: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements. This will calculate kQ for electrons
    for *CYLINDRICAL* chambers only according to Muir & Rodgers.

    Parameters
    ----------
    chamber : str
        The chamber of the chamber. Valid values are those listed in
        Tables VI and VII of Muir and Rodgers 2014.
    r_50 : float
        The R50 value in cm of an electron beam.
    """
    ch = KQ_ELECTRONS[chamber]
    return (ch['a'] + ch['b'] * r_50 ** -ch['c']) * ch['kQ,ecal']


class TG51Base(Structure):

    @property
    def p_tp(self) -> float:
        """Temperature/Pressure correction."""
        return p_tp(temp=self.temp, press=self.press)

    @property
    def p_ion(self) -> float:
        """Ionization collection correction."""
        return p_ion(voltage_reference=self.voltage_reference, voltage_reduced=self.voltage_reduced,
                     m_reference=self.m_reference,
                     m_reduced=self.m_reduced)

    @property
    def p_pol(self) -> float:
        """Polarity correction."""
        return p_pol(m_reference=self.m_reference, m_opposite=self.m_opposite)

    @property
    def m_corrected(self) -> float:
        """Corrected chamber reading."""
        return m_corrected(p_ion=self.p_ion, p_tp=self.p_tp, p_elec=self.p_elec, p_pol=self.p_pol,
                           m_reference=self.m_reference)

    @property
    def m_corrected_adjustment(self) -> float:
        """Corrected chamber reading after adjusting the output."""
        if self.m_reference_adjusted is not None:
            return m_corrected(p_ion=self.p_ion, p_tp=self.p_tp, p_elec=self.p_elec, p_pol=self.p_pol,
                               m_reference=self.m_reference_adjusted)

    @property
    def output_was_adjusted(self) -> float:
        """Boolean specifiying if output was adjusted."""
        return self.m_reference_adjusted is not None


class TG51Photon(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in a photon beam.

    Parameters
    ----------
    institution : str
        Institution name.
    physicist : str
        Physicist performing calibration.
    unit : str
        Unit name; e.g. TrueBeam1.
    measurement_date : str
        Date of measurement. E.g. 10/22/2018.
    temp : float
        The temperature in Celsius. Use :func:`~pylinac.calibration.tg51.fahrenheit2celsius` to convert if necessary.
    press : float
        The value of pressure in kPa. Can be converted from mmHg and mbar; see :func:`~pylinac.calibration.tg51.mmHg2kPa` and :func:`~pylinac.calibration.tg51.mbar2kPa`.
    energy : float
        Nominal energy of the beam in MV.
    chamber : str
        Chamber model. Must be one of the listed chambers in TG-51 Addendum.
    n_dw : float
        NDW value in Gy/nC.
    p_elec : float
        Electrometer correction factor; given by the calibration laboratory.
    measured_pdd10 : float
        The measured value of PDD(10); will be converted to PDDx(10) and used for calculating kq.
    lead_foil : {None, '50cm', '30cm'}
        Whether a lead foil was used to acquire PDD(10)x and where its position was. Used to calculate kq.
    clinical_pdd10 : float
        The PDD used to correct the dose at 10cm back to dmax. Usually the TPS PDD(10) value.
    voltage_reference : int
        Reference voltage; i.e. voltage when taking the calibration measurement.
    voltage_reduced : int
        Reduced voltage; usually half of the reference voltage.
    m_reference : float, tuple
        Ion chamber reading(s) at the reference voltage.
    m_opposite : float, tuple
        Ion chamber reading(s) at the opposite voltage of reference.
    m_reduced : float, tuple
        Ion chamber reading(s) at the reduced voltage.
    mu : int
        The MU delivered to measure the reference reading. E.g. 200.
    fff : bool
        Whether the beam is FFF or flat.
    tissue_correction : float
        Correction value to calibration to, e.g., muscle. A value of 1.0 means no correction (i.e. water).
    """

    @argue.options(chamber=KQ_PHOTONS.keys(), lead_foil=LEAD_OPTIONS.values())
    def __init__(self, *,
                 institution: str='',
                 physicist: str='',
                 unit: str,
                 measurement_date: str='',
                 temp: NumberLike,
                 press: NumberLike,
                 chamber: str,
                 n_dw: float,
                 p_elec: float,
                 electrometer: str='',
                 measured_pdd10: Optional[float]=None,
                 lead_foil: Optional[str]=None,
                 clinical_pdd10: float,
                 energy: int,
                 fff: bool=False,
                 voltage_reference: int,
                 voltage_reduced: int,
                 m_reference: NumberOrArray,
                 m_opposite: NumberOrArray,
                 m_reduced: NumberOrArray,
                 mu: int,
                 tissue_correction: float=1.0,
                 m_reference_adjusted: Optional[NumberOrArray]=None):
        super().__init__(temp=temp, press=press, chamber=chamber, n_dw=n_dw, p_elec=p_elec, measured_pdd10=measured_pdd10,
                         energy=energy, voltage_reference=voltage_reference, voltage_reduced=voltage_reduced,
                         m_reference=m_reference, m_opposite=m_opposite, m_reduced=m_reduced, clinical_pdd10=clinical_pdd10, mu=mu,
                         tissue_correction=tissue_correction, lead_foil=lead_foil, electrometer=electrometer,
                         m_reference_adjusted=m_reference_adjusted, institution=institution, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, fff=fff)
            # add check for tpr vs pdd

    @property
    def pddx(self) -> float:
        """The photon-only PDD(10) value."""
        return pddx(pdd=self.measured_pdd10, energy=self.energy, lead_foil=self.lead_foil)

    @property
    def kq(self) -> float:
        """The chamber-specific beam quality correction factor."""
        return kq_photon_pddx(chamber=self.chamber, pddx=self.pddx)

    @property
    def dose_mu_10(self) -> float:
        """cGy/MU at a depth of 10cm."""
        return self.tissue_correction * self.m_corrected * self.kq * self.n_dw / self.mu

    @property
    def dose_mu_dmax(self) -> float:
        """cGy/MU at a depth of dmax."""
        return self.dose_mu_10 / (self.clinical_pdd10 / 100)

    @property
    def dose_mu_10_adjusted(self) -> float:
        """The dose/mu at 10cm depth after adjustment."""
        return self.tissue_correction*self.m_corrected_adjustment*self.kq*self.n_dw/self.mu

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """The dose/mu at dmax depth after adjustment."""
        return self.dose_mu_10_adjusted / (self.clinical_pdd10 / 100)

    def publish_pdf(self, filename: str, notes: Optional[list]=None, open_file: bool=False,
                    metadata: Optional[dict]=None):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : str, file-like object
            The file to write the results to.
        notes : str, list
            Any notes to be added to the report. If a string, adds everything as one line.
            If a list, must be a list of strings; each string item will be a new line.
        open_file : bool
            Whether to open the file after creation. Will use the default PDF program.
        metadata : dict
            Any data that should be appended to every page of the report. This differs from notes in that
            metadata is at the top of every page while notes is at the bottom of the report.
        """
        was_adjusted = 'Yes' if self.output_was_adjusted else 'No'
        title = [
            'TG-51 Photon Report',
            f"{self.unit} - {self.energy} MV{' FFF' if self.fff else ''}"
        ]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            'Site Data:',
            f'Institution: {self.institution}',
            f'Performed by: {self.physicist}',
            f'Measurement Date: {self.measurement_date}',
            f'Date of Report: {datetime.now().strftime("%A, %B %d, %Y")}',
            f'Unit: {self.unit}',
            f"Energy: {self.energy} MV {'FFF' if self.fff else ''}",
            '',
            'Instrumentation:',
            f'Chamber: {self.chamber}',
            f'Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}',
            f'Electrometer: {self.electrometer}',
            f'Pelec: {self.p_elec:2.3f}',
            f'MU: {self.mu}',
            '',
            'Beam Quality:',
            f"Lead foil: {'No' if self.lead_foil is None else self.lead_foil}",
            f"Measured PDD(10){'' if self.lead_foil is None else 'Pb'} {self.measured_pdd10:2.2f}",
            f'Calculated PDD(10)x: {self.pddx:2.2f}',
            f'Determined kQ: {self.kq:2.3f}',
            '',
            'Chamber Corrections/Measurements:',
            f'Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}',
            f'Pressure (kPa): {self.press:2.1f}',
            f'Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}',
            f'Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}',
            f'Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}',
            f'Ptp: {self.p_tp:2.3f}',
            f'Pion: {self.p_ion:2.3f}',
            f'Ppol: {self.p_pol:2.3f}',
            '',
            'Dose Determination:',
            f'Fully corrected M (nC): {self.m_corrected:2.3f}',
            f'Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}',
            f'Dose/MU @ 10cm depth (cGy): {self.dose_mu_10:2.3f}',
            f'Clinical PDD (%): {self.clinical_pdd10:2.2f}',
            f'Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}',
            '',
            f'Output Adjustment?: {was_adjusted}',
        ]
        if was_adjusted == 'Yes':
            text.append(f'Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjusted}')
            text.append(f'Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}')
            text.append(f'Adjusted Dose/MU @ 10cm depth (cGy): {self.dose_mu_10_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}')
        canvas.add_text(text=text, location=(2, 25.5), font_size=12)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            open_path(filename)


class TG51ElectronLegacy(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in an electron beam.

    Parameters
    ----------
    institution : str
        Institution name.
    physicist : str
        Physicist performing calibration.
    unit : str
        Unit name; e.g. TrueBeam1.
    measurement_date : str
        Date of measurement. E.g. 10/22/2018.
    temp : float (17-27)
        The temperature in degrees Celsius.
    press : float (91-111)
        The value of pressure in kPa. Can be converted from mmHg and mbar; see :func:`~pylinac.calibration.tg51.mmHg2kPa` and :func:`~pylinac.calibration.tg51.mbar2kPa`.
    chamber : str
        Chamber model; only for bookkeeping.
    n_dw : float
        NDW value in Gy/nC. Given by the calibration laboratory.
    k_ecal : float
        Kecal value which is chamber specific. This value is the major difference between the legacy class and modern class where no kecal is needed.
    p_elec : float
        Electrometer correction factor; given by the calibration laboratory.
    clinical_pdd : float
        The PDD used to correct the dose back to dref.
    voltage_reference : float
        Reference voltage; i.e. voltage when taking the calibration measurement.
    voltage_reduced : float
        Reduced voltage; usually half of the reference voltage.
    m_reference : float, tuple
        Ion chamber reading(s) at the reference voltage.
    m_opposite : float, tuple
        Ion chamber reading(s) at the opposite voltage of reference.
    m_reduced : float, tuple
        Ion chamber reading(s) at the reduced voltage.
    mu : int
        The MU delivered to measure the reference reading. E.g. 200.
    i_50 : float
        Depth of 50% ionization.
    tissue_correction : float
        Correction value to calibration to, e.g., muscle. A value of 1.0 means no correction (i.e. water).
    """

    def __init__(self, *,
                 institution: str='',
                 physicist: str='',
                 unit: str='',
                 measurement_date: str='',
                 energy: int,
                 temp: NumberLike,
                 press: NumberLike,
                 chamber: str,
                 k_ecal: float,
                 n_dw: float,
                 electrometer: str='',
                 p_elec: float,
                 clinical_pdd: float,
                 voltage_reference: int,
                 voltage_reduced: int,
                 m_reference: NumberOrArray,
                 m_opposite: NumberOrArray,
                 m_reduced: NumberOrArray,
                 m_gradient: NumberOrArray,
                 cone: str,
                 mu: int,
                 i_50: float,
                 tissue_correction: float=1.0,
                 m_reference_adjusted=None):
        super().__init__(temp=temp, press=press, chamber=chamber, n_dw=n_dw, p_elec=p_elec,
                         voltage_reference=voltage_reference, voltage_reduced=voltage_reduced, m_reference=m_reference,
                         m_opposite=m_opposite, m_reduced=m_reduced, clinical_pdd=clinical_pdd, mu=mu,
                         i_50=i_50, tissue_correction=tissue_correction,
                         institution=institution, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, electrometer=electrometer,
                         m_reference_adjusted=m_reference_adjusted, cone=cone, energy=energy, k_ecal=k_ecal,
                         m_gradient=m_gradient)

    @property
    def r_50(self) -> float:
        """Depth of the 50% dose value."""
        return r_50(i_50=self.i_50)

    @property
    def dref(self) -> float:
        """Depth of the reference point."""
        return d_ref(i_50=self.i_50)

    @property
    def pq_gr(self):
        """Gradient factor"""
        return pq_gr(m_dref_plus=self.m_gradient, m_dref=self.m_reference)

    @property
    def kq(self) -> float:
        """The kQ value using classic TG-51"""
        return self.k_ecal * kp_r50(r_50=self.r_50)

    @property
    def dose_mu_dref(self) -> float:
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.m_corrected * self.kq * self.n_dw * self.pq_gr / self.mu

    @property
    def dose_mu_dmax(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref / (self.clinical_pdd / 100)

    @property
    def dose_mu_dref_adjusted(self) -> float:
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.m_corrected_adjustment * self.kq * self.n_dw * self.pq_gr / self.mu

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref_adjusted / (self.clinical_pdd / 100)

    def publish_pdf(self, filename: str, notes: Optional[list]=None, open_file: bool=False,
                    metadata: Optional[dict]=None):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : str, file-like object
            The file to write the results to.
        notes : str, list
            Any notes to be added to the report. If a string, adds everything as one line.
            If a list, must be a list of strings; each string item will be a new line.
        open_file : bool
            Whether to open the file after creation. Will use the default PDF program.
        metadata : dict
            Any data that should be appended to every page of the report. This differs from notes in that
            metadata is at the top of every page while notes is at the bottom of the report.
        """
        was_adjusted = 'Yes' if self.output_was_adjusted else 'No'
        title = [
            'TG-51 Electron Report (Legacy)',
            f'{self.unit} - {self.energy} MeV'
        ]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            'Site Data:',
            f'Institution: {self.institution}',
            f'Performed by: {self.physicist}',
            f'Measurement Date: {self.measurement_date}',
            f'Date of Report: {datetime.now().strftime("%A, %B %d, %Y")}',
            f'Unit: {self.unit}',
            f'Energy: {self.energy} MeV',
            f'Cone: {self.cone}',
            f'MU: {self.mu}',
            '',
            'Instrumentation:',
            f'Chamber chamber: {self.chamber}',
            f'Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}',
            f'Electrometer: {self.electrometer}',
            f'Pelec: {self.p_elec:2.2f}',
            '',
            'Beam Quality:',
            f'I50 (cm): {self.i_50:2.2f}',
            f'R50 (cm): {self.r_50:2.2f}',
            f'Dref (cm): {self.dref:2.2f}',
            f'Kecal: {self.k_ecal:2.3f}',
            f"kQ: {self.kq:2.3f}",
            '',
            'Chamber Corrections/Measurements:',
            f'Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}',
            f'Pressure (kPa): {self.press:2.1f}',
            f'Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}',
            f'Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}',
            f'Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}',
            f'Ptp: {self.p_tp:2.3f}',
            f'Pion: {self.p_ion:2.3f}',
            f'Ppol: {self.p_pol:2.3f}',
            f'Mraw @ Dref + 0.5rcav (nC): {self.m_gradient}',
            '',
            'Dose Determination:',
            f'Fully corrected M (nC): {self.m_corrected:2.3f}',
            f'Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}',
            f'Dose/MU @ Dref depth (cGy): {self.dose_mu_dref:2.3f}',
            f'Clinical PDD (%): {self.clinical_pdd:2.2f}',
            f'Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}',
            '',
            f'Output Adjustment?: {was_adjusted}',
        ]
        if was_adjusted == 'Yes':
            text.append(f'Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjustment}')
            text.append(f'Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}')
            text.append(f'Adjusted Dose/MU @ dref depth (cGy): {self.dose_mu_dref_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}')
        canvas.add_text(text=text, location=(2, 25.5), font_size=11)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            open_path(filename)


class TG51ElectronModern(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in an electron beam.

    .. warning::
        This class uses the values of Muir & Rodgers. These values are likely to be included in the new TG-51
        addendum, but are not official. The results can be up to 1% different. Physicists should use their own
        judgement when deciding which class to use. To use a manual kecal value, Pgradient and the classic TG-51 equations use
        the :class:`~pylinac.calibration.tg51.TG51ElectronLegacy` class.

    Parameters
    ----------
    institution : str
        Institution name.
    physicist : str
        Physicist performing calibration.
    unit : str
        Unit name; e.g. TrueBeam1.
    measurement_date : str
        Date of measurement. E.g. 10/22/2018.
    press : float
        The value of pressure in kPa. Can be converted from mmHg and mbar; see :func:`~pylinac.calibration.tg51.mmHg2kPa` and :func:`~pylinac.calibration.tg51.mbar2kPa`.
    temp : float
        The temperature in Celsius.
    voltage_reference : int
        The reference voltage; i.e. the voltage for the calibration reading (e.g. 300V).
    voltage_reduced : int
        The reduced voltage, usually a fraction of the reference voltage (e.g. 150V).
    m_reference : array, float
        The reading(s) of the chamber at reference voltage.
    m_reduced : array, float
        The reading(s) of the chamber at the reduced voltage.
    m_opposite : array, float
        The reading(s) of the chamber at the opposite voltage from reference. Sign of the reading does not matter.
    k_elec : float
        The electrometer correction value given by the calibration laboratory. jyh,lykllp;ljljuhyk nmdrzj
    chamber : str
        Ion chamber model.
    n_dw : float
        NDW value in Gy/nC
    p_elec : float
        Electrometer correction given by the calibration laboratory.
    clinical_pdd : float
        The PDD used to correct the dose back to dref.
    mu : int
        MU delivered.
    i_50 : float
        Depth of 50% ionization
    tissue_correction : float
        Correction value to calibration to, e.g., muscle. A value of 1.0 means no correction (i.e. water).
    """

    def __init__(self, *,
                 institution: str='',
                 physicist: str='',
                 unit: str='',
                 measurement_date: str='',
                 energy: int,
                 temp: NumberLike,
                 press: NumberLike,
                 chamber: str,
                 n_dw: float,
                 electrometer: str='',
                 p_elec: float,
                 clinical_pdd: float,
                 voltage_reference: int,
                 voltage_reduced: int,
                 m_reference: NumberOrArray,
                 m_opposite: NumberOrArray,
                 m_reduced: NumberOrArray,
                 cone: str,
                 mu: int,
                 i_50: float,
                 tissue_correction: float,
                 m_reference_adjusted=None):
        super().__init__(temp=temp, press=press, chamber=chamber, n_dw=n_dw, p_elec=p_elec,
                         voltage_reference=voltage_reference, voltage_reduced=voltage_reduced, m_reference=m_reference,
                         m_opposite=m_opposite, m_reduced=m_reduced, clinical_pdd=clinical_pdd, mu=mu,
                         i_50=i_50, tissue_correction=tissue_correction,
                         institution=institution, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, electrometer=electrometer,
                         m_reference_adjusted=m_reference_adjusted, cone=cone, energy=energy,
                         )

    @property
    def r_50(self) -> float:
        """Depth of the 50% dose value."""
        return r_50(i_50=self.i_50)

    @property
    def dref(self) -> float:
        """Depth of the reference point."""
        return d_ref(i_50=self.i_50)

    @property
    def kq(self) -> float:
        """The kQ value using the updated Muir & Rodgers values from their 2014 paper, equation 11, or classically
        if kecal is passed."""
        return kq_electron(chamber=self.chamber, r_50=self.r_50)

    @property
    def dose_mu_dref(self) -> float:
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.m_corrected * self.kq * self.n_dw / self.mu

    @property
    def dose_mu_dmax(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref / (self.clinical_pdd / 100)

    @property
    def dose_mu_dref_adjusted(self) -> float:
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.m_corrected_adjusted * self.kq * self.n_dw / self.mu

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref_adjusted / (self.clinical_pdd / 100)

    def publish_pdf(self, filename: str, notes: Optional[list]=None, open_file: bool=False,
                    metadata: Optional[dict]=None):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : str, file-like object
            The file to write the results to.
        notes : str, list
            Any notes to be added to the report. If a string, adds everything as one line.
            If a list, must be a list of strings; each string item will be a new line.
        open_file : bool
            Whether to open the file after creation. Will use the default PDF program.
        metadata : dict
            Any data that should be appended to every page of the report. This differs from notes in that
            metadata is at the top of every page while notes is at the bottom of the report.
        """
        was_adjusted = 'Yes' if self.output_was_adjusted else 'No'
        title = [
            'TG-51 Electron Report (Modern)',
            f'{self.unit} - {self.energy} MeV'
        ]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            'Site Data:',
            f'Institution: {self.institution}',
            f'Performed by: {self.physicist}',
            f'Measurement Date: {self.measurement_date}',
            f'Date of Report: {datetime.now().strftime("%A, %B %d, %Y")}',
            f'Unit: {self.unit}',
            f'Energy: {self.energy} MeV',
            f'Cone: {self.cone}',
            f'MU: {self.mu}',
            '',
            'Instrumentation:',
            f'Chamber: {self.chamber}',
            f'Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}',
            f'Electrometer: {self.electrometer}',
            f'Pelec: {self.p_elec:2.2f}',
            '',
            'Beam Quality:',
            f'I50 (cm): {self.i_50:2.2f}',
            f'R50 (cm): {self.r_50:2.2f}',
            f'Dref (cm): {self.dref:2.2f}',
            f"Calculated kQ: {self.kq:2.3f}",
            '',
            'Chamber Corrections/Measurements:',
            f'Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}',
            f'Pressure (kPa): {self.press:2.1f}',
            f'Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}',
            f'Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}',
            f'Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}',
            f'Ptp: {self.p_tp:2.3f}',
            f'Pion: {self.p_ion:2.3f}',
            f'Ppol: {self.p_pol:2.3f}',
            '',
            'Dose Determination:',
            f'Fully corrected M (nC): {self.m_corrected:2.3f}',
            f'Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}',
            f'Dose/MU @ Dref depth (cGy): {self.dose_mu_dref:2.3f}',
            f'Clinical PDD (%): {self.clinical_pdd:2.2f}',
            f'Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}',
            '',
            f'Output Adjustment?: {was_adjusted}',
        ]
        if was_adjusted == 'Yes':
            text.append(f'Adjusted corrected M @ reference voltage (nC): {self.m_corrected_adjustment}')
            text.append(f'Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}')
            text.append(f'Adjusted Dose/MU @ dref depth (cGy): {self.dose_mu_dref_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}')
        canvas.add_text(text=text, location=(2, 25.5), font_size=11)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            open_path(filename)
