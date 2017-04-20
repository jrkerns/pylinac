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

import numpy as np
from reportlab.lib.units import cm

from .core.utilities import Structure
from .core import pdf


CHAMBERS_PHOTONS = {
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

CHAMBERS_ELECTRONS = {
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
    'FC65G': {'kQ,ecal': 0.904, 'a': 0.971, 'b': 0.113, 'c': 0.680},
    'FC65P': {'kQ,ecal': 0.902, 'a': 0.973, 'b': 0.110, 'c': 0.692},
    'FC23C': {'kQ,ecal': 0.904, 'a': 0.971, 'b': 0.097, 'c': 0.591},
    'CC25': {'kQ,ecal': 0.904, 'a': 0.964, 'b': 0.105, 'c': 0.539},
    'CC13': {'kQ,ecal': 0.904, 'a': 0.926, 'b': 0.129, 'c': 0.279},

    # Other
    'PR06C/G': {'kQ,ecal': 0.906, 'a': 0.972, 'b': 0.122, 'c': 0.729},
    '2571': {'kQ,ecal': 0.903, 'a': 0.977, 'b': 0.117, 'c': 0.817},
    '2611': {'kQ,ecal': 0.896, 'a': 0.979, 'b': 0.120, 'c': 0.875},
}


def p_tp(temp=22, press=760):
    """Calculate the temperature & pressure correction.

    Parameters
    ----------
    temp : float
        The temperature in degrees Celsius.
    press : float
        The pressure in mmHg.
    """
    return (760/press)*((273.2+temp)/295.2)


def p_pol(m_reference=(1, 2), m_opposite=(-3, -4)):
    """Calculate the polarity correction.

    Parameters
    ----------
    m_reference : iterable
        The readings of the ion chamber at the reference polarity and voltage.
    m_opposite : iterable
        The readings of the ion chamber at the polarity opposite the reference.
        This value should be of the opposite sign of the M reference value.
        If it's not, its sign will automatically be flipped.
    """
    mref_avg = np.mean(m_reference)
    mopp_avg = np.mean(m_opposite)
    # if same sign given, flip one.
    # Technically, they're opposite charges, but most physicists pass positive values for both
    if np.sign(mref_avg) == np.sign(mopp_avg):
        mopp_avg = -mopp_avg
    return (mref_avg - mopp_avg)/(2*mref_avg)


def p_ion(volt_high=300, volt_low=150, m_high=(1, 2), m_low=(3, 4)):
    """Calculate the ion chamber collection correction.

    Parameters
    ----------
    volt_high : int
        The "high" voltage; same as the TG51 measurement voltage.
    volt_low : int
        The "low" voltage; usually half of the high voltage.
    m_high : float, iterable
        The readings of the ion chamber at the "high" voltage.
    m_low : float, iterable
        The readings of the ion chamber at the "low" voltage.
    """
    return (1 - volt_high/volt_low)/(np.mean(m_high)/np.mean(m_low) - volt_high/volt_low)


def d_ref(i_50):
    """Calculate the dref of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    r50 = r_50(i_50)
    return 0.6*r50-0.1


def r_50(i_50):
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


def kp_r50(r_50):
    """Calculate k'R50 for Farmer-like chambers.

    Parameters
    ----------
    r_50 : float
        The R50 value in cm.
    """
    if r_50 >= 9 or r_50 <= 2:
        raise ValueError("Cannot calculate k prime R50 with an R50 value of <=2 or >=9cm")
    return 0.9905+0.071*np.exp(-r_50/3.67)


def pq_gr(m_dref_plus=(1, 2), m_dref=(3, 4)):
    """Calculate PQ_gradient for a cylindrical chamber.

    Parameters
    ----------
    m_dref_plus : float, iterable
        The readings of the ion chamber at dref + 0.5rcav.
    m_dref : float, iterable
        The readings of the ion chamber at dref.
    """
    return np.mean(m_dref_plus) / np.mean(m_dref)


def m_corrected(p_ion=1.0, p_tp=1.0, p_elec=1.0, p_pol=1.0, m_raw=(1.1, 2.2)):
    """Calculate M_corrected, the ion chamber reading with all corrections applied.

    Parameters
    ----------
    p_ion : float
        The ion collection correction.
    p_tp : float
        The temperature & pressure correction.
    p_elec : float
        The electrometer correction.
    p_pol : float
        The polarity correction.
    m_raw : float, iterable
        The raw ion chamber readings.

    Returns
    -------
    float
    """
    return p_ion*p_tp*p_elec*p_pol*np.mean(m_raw)


def pddx(pdd=66.4, energy=6, lead_foil=None):
    """Calculate PDDx based on the PDD.

    Parameters
    ----------
    pdd : {>0.627, <0.890}
        The measured PDD. If lead foil was used, this assumes the pdd as measured with the lead in place.
    energy : int
        The nominal energy in MV.
    lead_foil : {None, '30cm', '50cm'}
        Applicable only for energies >10MV.
        Whether a lead foil was used to acquire the pdd.
        Use ``None`` if no lead foil was used and the interim equation should be used.
        Use ``50cm`` if the lead foil was set to 50cm from the phantom surface.
        Use ``30cm`` if the lead foil was set to 30cm from the phantom surface.
    """
    if energy < 10:
        return pdd
    elif energy >= 10:
        if lead_foil is None:
            return 1.267*pdd-20
        elif lead_foil == '50cm':
            if pdd < 73:
                return pdd
            else:
                return (0.8905+0.0015*pdd)*pdd
        elif lead_foil == '30cm':
            if pdd < 71:
                return pdd
            else:
                return (0.8116+0.00264*pdd)*pdd


def kq(model='30010', pddx=None, tpr=None, r_50=None):
    """Calculate kQ based on the model and clinical measurements. This will calculate kQ for both photons and electrons 
    for *CYLINDRICAL* chambers only.

    Parameters
    ----------
    model : str
        The model of the chamber. Valid values are those listed in
        Table III of Muir and Rodgers and Table I of the TG-51 Addendum.
    pddx : {>=0.627, <=0.861}
        The *PHOTON-ONLY* PDD measurement at 10cm depth for a 10x10cm2 field.
    tpr : {>=0.623, <=0.805}
        The TPR ratio of the 20cm measurement divided by the 10cm measurement.
    r_50 : float
        The R50 value in cm of an electron beam.

    .. warning::
        Only 1 of ``pddx``, ``tpr`` or ``r_50`` can be defined.
    """
    PDD_LOW = 62.7
    PDD_HIGH = 86.1
    TPR_LOW = 0.623
    TPR_HIGH = 0.805

    # error checking
    if not any((pddx, tpr, r_50)):
        raise ValueError("At least one of the parameters pddx, tpr, or r_50 must be defined.")
    if pddx is not None and tpr is not None:
        raise ValueError("Only the PDD or TPR parameter can be defined, not both.")
    if any((pddx, tpr)) and r_50 is not None:
        raise ValueError("Cannot define both a photon component (PDDx, TPR) and an electron component (R50)")

    if pddx is not None:
        if pddx > PDD_HIGH or pddx < PDD_LOW:
            raise ValueError("Measured PDD is out of range; must be between {:2.2} and {:2.2}.".format(PDD_LOW, PDD_HIGH))
        else:
            ch = CHAMBERS_PHOTONS[model]
            return ch["a"] + ch["b"]*pddx + ch["c"]*(pddx**2)

    if tpr is not None:
        if tpr > TPR_HIGH or tpr < TPR_LOW:
            raise ValueError("Measured TPR is out of range; must be between {:2.2} and {:2.2}.".format(TPR_LOW, TPR_HIGH))
        else:
            ch = CHAMBERS_PHOTONS[model]
            return ch["a'"] + ch["b'"]*tpr + ch["c'"]*(tpr**2) + ch["d'"]*(tpr**3)

    if r_50 is not None:
        ch = CHAMBERS_ELECTRONS[model]
        return (ch['a'] + ch['b'] * r_50**-ch['c']) * ch['kQ,ecal']


class TG51Base(Structure):

    @property
    def p_tp(self):
        """Temperature/Pressure correction."""
        return p_tp(self.temp, self.press)

    @property
    def p_ion(self):
        """Ionization collection correction."""
        return p_ion(self.volt_high, self.volt_low, self.m_raw, self.m_low)

    @property
    def p_pol(self):
        """Polarity correction."""
        return p_pol(self.m_raw, self.m_opp)

    @property
    def m_corrected(self):
        """Corrected chamber reading."""
        return m_corrected(self.p_ion, self.p_tp, self.p_elec, self.p_pol, self.m_raw)

    @property
    def adjusted_m_corrected(self):
        """Corrected chamber reading after adjusting the output."""
        return m_corrected(self.p_ion, self.p_tp, self.p_elec, self.p_pol, self.adjusted_m_raw)

    @property
    def output_was_adjusted(self):
        """Boolean specifiying if output was adjusted."""
        return self.adjusted_m_raw is not None


class TG51Photon(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in a photon beam.

    Attributes
    ----------
    temp : float
    press : float
    energy : float
        Nominal energy of the beam in MV.
    model : str
        Chamber model
    n_dw : float
        NDW value in Gy/nC
    p_elec : float
    measured_pdd : float
        The measured value of PDD(10); used for calculating kq.
    lead_foil : {None, '50cm', '30cm'}
        Whether a lead foil was used to acquire PDD(10)x and where its position was. Used to calculate kq.
    clinical_pdd : float
        The PDD used to correct the dose at 10cm back to dmax. Usually the TPS PDD(10) value.
    volt_high : float
    volt_low : float
    m_raw : float, tuple
    m_opp : float, tuple
    m_low : float, tuple
    mu : float
    tissue_correction : float
        Correction value to calibration to, e.g., muscle. A value of 1.0 means no correction (i.e. water).
    """

    def __init__(self,
                 institution='',
                 physicist='',
                 unit='',
                 measurement_date='',
                 temp=22,
                 press=760,
                 model='30010',
                 n_dw=5.9,
                 p_elec=1.0,
                 electrometer='',
                 measured_pdd=66.4,
                 lead_foil=None,
                 clinical_pdd=66.4,
                 energy=6,
                 fff=False,
                 volt_high=300,
                 volt_low=150,
                 m_raw=(1, 2),
                 m_opp=(1, 2),
                 m_low=(1, 2),
                 mu=200,
                 tissue_correction=1.0,
                 adjusted_m_raw=None):
        super().__init__(temp=temp, press=press, model=model, n_dw=n_dw, p_elec=p_elec, measured_pdd=measured_pdd,
                         energy=energy, volt_high=volt_high, volt_low=volt_low, m_raw=m_raw,
                         m_opp=m_opp, m_low=m_low, clinical_pdd=clinical_pdd, mu=mu,
                         tissue_correction=tissue_correction, lead_foil=lead_foil, electrometer=electrometer,
                         adjusted_m_raw=adjusted_m_raw, institution=institution, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, fff=fff)

    @property
    def pddx(self):
        """The photon-only PDD(10) value."""
        return pddx(self.measured_pdd, self.energy, self.lead_foil)

    @property
    def kq(self):
        """The chamber-specific beam quality correction factor."""
        return kq(self.model, self.pddx)

    @property
    def dose_mu_10(self):
        """cGy/MU at a depth of 10cm."""
        return self.tissue_correction * self.m_corrected * self.kq * self.n_dw / self.mu

    @property
    def dose_mu_dmax(self):
        """cGy/MU at a depth of dmax."""
        return self.dose_mu_10 / (self.clinical_pdd / 100)

    @property
    def adjusted_dose_mu_10(self):
        """The dose/mu at 10cm depth after adjustment."""
        return self.tissue_correction*self.adjusted_m_corrected*self.kq*self.n_dw/self.mu

    @property
    def adjusted_dose_mu_dmax(self):
        """The dose/mu at dmax depth after adjustment."""
        return self.adjusted_dose_mu_10/(self.clinical_pdd/100)

    def publish_pdf(self, filename, notes=None, open_file=False):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : str, file-like object
            The file to write the results to.
        notes : str, list
            Any notes to be added to the report. If a string, adds everything as one line.
            If a list, must be a list of strings; each string item will be a new line.
        """
        was_adjusted = 'Yes' if self.output_was_adjusted else 'No'
        title = 'TG-51 Photon Report - {} MV'.format(self.energy)
        if self.fff:
            title += ' FFF'

        canvas = pdf.create_pylinac_page_template(filename, analysis_title=title)
        text = ['Site Data:',
                'Institution: {}'.format(self.institution),
                'Performed by: {}'.format(self.physicist),
                'Measurement Date: {}'.format(self.measurement_date),
                'Date of Report: {}'.format(datetime.now().strftime("%A, %B %d, %Y")),
                'Unit: {}'.format(self.unit),
                'Energy: {} MV {}'.format(self.energy, 'FFF' if self.fff else ''),
                '',
                'Instrumentation:',
                'Chamber model: {}'.format(self.model),
                'Chamber Calibration Factor Ndw (cGy/nC): {:2.3f}'.format(self.n_dw),
                'Electrometer: {}'.format(self.electrometer),
                'Pelec: {}'.format(self.p_elec),
                'MU: {}'.format(self.mu),
                '',
                'Beam Quality:',
                'Lead foil used?: {}'.format('No' if self.lead_foil is None else self.lead_foil),
                'Measured %dd(10) (this is %dd(10)Pb if lead was used): {:2.2f}'.format(self.measured_pdd),
                'Calculated %dd(10)x: {:2.2f}'.format(self.pddx),
                'Determined kQ: {:2.3f}'.format(self.kq),
                '',
                'Chamber Corrections/Measurements:',
                'Temperature (\N{DEGREE SIGN}C): {:2.1f}'.format(self.temp),
                'Pressure (mmHg): {:2.1f}'.format(self.press),
                'Ptp: {:2.3f}'.format(self.p_tp),
                'Reference voltage (V): {}'.format(self.volt_high),
                'Mraw @ reference voltage (nC): {}'.format(self.m_raw),
                '"Lower" voltage (V): {}'.format(self.volt_low),
                'Mraw @ "lower" voltage (nC): {}'.format(self.m_low),
                'Opposite voltage (V): {}'.format(-self.volt_high),
                'Mraw @ opposite voltage (nC): {}'.format(self.m_opp),
                'Pion: {:2.3f}'.format(self.p_ion),
                'Ppol: {:2.3f}'.format(self.p_pol),
                '',
                'Dose Determination:',
                'Fully corrected M (nC): {:2.3f}'.format(self.m_corrected),
                'Tissue correction (e.g. muscle): {:2.3f}'.format(self.tissue_correction),
                'Dose/MU @ 10cm depth (cGy): {:2.3f}'.format(self.dose_mu_10),
                'Clinical PDD (%): {:2.2f}'.format(self.clinical_pdd),
                'Dose/MU @ dmax (cGy): {:2.3f}'.format(self.dose_mu_dmax),
                '',
                'Output Adjustment?: {}'.format(was_adjusted),
                ]
        if was_adjusted == 'Yes':
            text.append('Adjusted Mraw @ reference voltage (nC): {}'.format(self.adjusted_m_raw))
            text.append('Adjusted fully corrected M (nC): {:2.3f}'.format(self.adjusted_m_corrected))
            text.append('Adjusted Dose/MU @ 10cm depth (cGy): {:2.3f}'.format(self.adjusted_dose_mu_10))
            text.append('Adjusted Dose/MU @ dmax (cGy): {:2.3f}'.format(self.adjusted_dose_mu_dmax))
        pdf.draw_text(canvas, x=2 * cm, y=25.5 * cm, fontsize=12,
                      text=text)
        if notes is not None:
            pdf.draw_text(canvas, x=12 * cm, y=6.5 * cm, fontsize=14, text="Notes:")
            pdf.draw_text(canvas, x=12 * cm, y=6 * cm, text=notes)
        pdf.finish(canvas, open_file=open_file, filename=filename)


class TG51Electron(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in an electron beam.

    Attributes
    ----------
    temp : float
    press : float
    model : str
        Chamber model
    n_dw : float
        NDW value in Gy/nC
    p_elec : float
    clinical_pdd : float
        The PDD used to correct the dose back to dref.
    volt_high : float
    volt_low : float
    m_raw : float, tuple
    m_opp : float, tuple
    m_low : float, tuple
    mu : float
    i_50 : float
        Depth of 50% ionization
    tissue_correction : float
        Correction value to calibration to, e.g., muscle. A value of 1.0 means no correction (i.e. water).
    """

    def __init__(self,
                 institution='',
                 physicist='',
                 unit='',
                 measurement_date='',
                 energy=9,
                 temp=22,
                 press=760,
                 model='30010',
                 n_dw=5.9,
                 electrometer='',
                 p_elec=1.0,
                 clinical_pdd=99.0,
                 volt_high=300,
                 volt_low=150,
                 m_raw=(1, 2),
                 m_opp=(1, 2),
                 m_low=(1, 2),
                 cone='15x15',
                 mu=200,
                 i_50=4,
                 tissue_correction=1.0,
                 adjusted_m_raw=None):
        super().__init__(temp=temp, press=press, model=model, n_dw=n_dw, p_elec=p_elec,
                         volt_high=volt_high, volt_low=volt_low, m_raw=m_raw,
                         m_opp=m_opp, m_low=m_low, clinical_pdd=clinical_pdd, mu=mu,
                         i_50=i_50, tissue_correction=tissue_correction,
                         institution=institution, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, electrometer=electrometer,
                         adjusted_m_raw=adjusted_m_raw, cone=cone, energy=energy)

    @property
    def r_50(self):
        """Depth of the 50% dose value."""
        return r_50(self.i_50)

    @property
    def dref(self):
        """Depth of the reference point."""
        return d_ref(self.i_50)

    @property
    def kq(self):
        """The kQ value using the updated Muir & Rodgers values from their 2014 paper, equation 11."""
        return kq(self.model, r_50=self.r_50)

    @property
    def dose_mu_dref(self):
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.m_corrected * self.kq * self.n_dw / self.mu

    @property
    def dose_mu_dmax(self):
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref / (self.clinical_pdd / 100)

    @property
    def adjusted_dose_mu_dref(self):
        """cGy/MU at the depth of Dref."""
        return self.tissue_correction * self.adjusted_m_corrected * self.kq * self.n_dw / self.mu

    @property
    def adjusted_dose_mu_dmax(self):
        """cGy/MU at the depth of dmax."""
        return self.adjusted_dose_mu_dref / (self.clinical_pdd / 100)

    def publish_pdf(self, filename, notes=None):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : str, file-like object
            The file to write the results to.
        notes : str, list
            Any notes to be added to the report. If a string, adds everything as one line.
            If a list, must be a list of strings; each string item will be a new line.
        """
        was_adjusted = 'Yes' if self.output_was_adjusted else 'No'
        title = 'TG-51 Photon Report - {} MeV'.format(self.energy)

        canvas = pdf.create_pylinac_page_template(filename, analysis_title=title)
        text = ['Site Data:',
                'Institution: {}'.format(self.institution),
                'Performed by: {}'.format(self.physicist),
                'Measurement Date: {}'.format(self.measurement_date),
                'Date of Report: {}'.format(datetime.now().strftime("%A, %B %d, %Y")),
                'Unit: {}'.format(self.unit),
                'Energy: {} MeV'.format(self.energy),
                'Cone: {}'.format(self.cone),
                'MU: {}'.format(self.mu),
                '',
                'Instrumentation:',
                'Chamber model: {}'.format(self.model),
                'Chamber Calibration Factor Ndw (cGy/nC): {:2.3f}'.format(self.n_dw),
                'Electrometer: {}'.format(self.electrometer),
                'Pelec: {}'.format(self.p_elec),
                '',
                'Beam Quality:',
                'I50 (cm): {:2.2f}'.format(self.i_50),
                'R50 (cm): {:2.2f}'.format(self.r_50),
                'Dref (cm): {:2.2f}'.format(self.dref),
                "kQ: {:2.3f}".format(self.kq),
                '',
                'Chamber Corrections/Measurements:',
                'Temperature (\N{DEGREE SIGN}C): {:2.1f}'.format(self.temp),
                'Pressure (mmHg): {:2.1f}'.format(self.press),
                'Ptp: {:2.3f}'.format(self.p_tp),
                'Reference voltage (V): {}'.format(self.volt_high),
                'Mraw @ reference voltage (nC): {}'.format(self.m_raw),
                '"Lower" voltage (V): {}'.format(self.volt_low),
                'Mraw @ "lower" voltage (nC): {}'.format(self.m_low),
                'Opposite voltage (V): {}'.format(-self.volt_high),
                'Mraw @ opposite voltage (nC): {}'.format(self.m_opp),
                'Pion: {:2.3f}'.format(self.p_ion),
                'Ppol: {:2.3f}'.format(self.p_pol),
                'Mraw @ Dref + 0.5rcav (nC): {}'.format(self.m_plus),
                '',
                'Dose Determination:',
                'Fully corrected M (nC): {:2.3f}'.format(self.m_corrected),
                'Tissue correction (e.g. muscle): {:2.3f}'.format(self.tissue_correction),
                'Dose/MU @ 10cm depth (cGy): {:2.3f}'.format(self.dose_mu_dref),
                'Clinical PDD (%): {:2.2f}'.format(self.clinical_pdd),
                'Dose/MU @ dmax (cGy): {:2.3f}'.format(self.dose_mu_dmax),
                '',
                'Output Adjustment?: {}'.format(was_adjusted),
                ]
        if was_adjusted == 'Yes':
            text.append('Adjusted Mraw @ reference voltage (nC): {}'.format(self.adjusted_m_raw))
            text.append('Adjusted fully corrected M (nC): {:2.3f}'.format(self.adjusted_m_corrected))
            text.append('Adjusted Dose/MU @ dref depth (cGy): {:2.3f}'.format(self.adjusted_dose_mu_dref))
            text.append('Adjusted Dose/MU @ dmax (cGy): {:2.3f}'.format(self.adjusted_dose_mu_dmax))
        pdf.draw_text(canvas, x=2 * cm, y=25.5 * cm, fontsize=11,
                      text=text)
        if notes is not None:
            pdf.draw_text(canvas, x=12 * cm, y=6.5 * cm, fontsize=14, text="Notes:")
            pdf.draw_text(canvas, x=12 * cm, y=6 * cm, text=notes)
        canvas.showPage()
        canvas.save()
