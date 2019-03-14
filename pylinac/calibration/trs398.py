from datetime import datetime
from typing import Union, List, Optional

import argue
import numpy as np

from pylinac.core.pdf import PylinacCanvas
from . import tg51 as _tg51
from .tg51 import mmHg2kPa, mbar2kPa, fahrenheit2celsius, tpr2010_from_pdd2010  # make available to module
from ..core.utilities import is_close, Structure, open_path
from ..core.typing import NumberOrArray


V1_V2_FITS = {
    2.0: {'a0': 2.337, 'a1': -3.636, 'a2': 2.299},
    2.5: {'a0': 1.474, 'a1': -1.587, 'a2': 1.114},
    3.0: {'a0': 1.198, 'a1': -0.875, 'a2': 0.677},
    3.5: {'a0': 1.080, 'a1': -0.542, 'a2': 0.463},
    4.0: {'a0': 1.022, 'a1': -0.363, 'a2': 0.341},
    5.0: {'a0': 0.975, 'a1': -0.188, 'a2': 0.214}
}

KQ_PHOTON_TPRS = (0.50, 0.53, 0.56, 0.59, 0.62, 0.65, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84)

# chamber kQ values from table 6.III
KQ_PHOTON_CHAMBERS = {

    # Capintec
    'PR-06C/G': (1.001, 1.001, 1.000, 0.998, 0.998, 0.995, 0.992, 0.990, 0.988, 0.984, 0.980, 0.972, 0.965,
                 0.956, 0.944),

    # Exradin
    'A12': (1.001, 1.001, 1.000, 1.000, 0.999, 0.997, 0.994, 0.992, 0.990, 0.986, 0.981, 0.974, 0.966, 0.957, 0.944),

    # Nuclear Associates
    '30-751': (1.002, 1.002, 1.000, 0.999, 0.997, 0.994, 0.991, 0.989, 0.985, 0.981, 0.977, 0.969, 0.961, 0.953, 0.940),
    '30-752': (1.004, 1.003, 1.001, 1.000, 0.998, 0.996, 0.993, 0.991, 0.989, 0.985, 0.981, 0.974, 0.967, 0.959, 0.947),

    # NE
    '2505': (1.001, 1.001, 1.000, 0.999, 0.997, 0.994, 0.991, 0.988, 0.984, 0.980, 0.975, 0.967, 0.959, 0.950, 0.937),
    '2505/A': (1.005, 1.003, 1.001, 0.997, 0.995, 0.990, 0.985, 0.982, 0.978, 0.974, 0.969, 0.962, 0.955, 0.947, 0.936),
    '2505/3, 3A': (1.005, 1.004, 1.002, 1.000, 0.998, 0.995, 0.993, 0.991, 0.989, 0.986, 0.982, 0.975, 0.969, 0.961,
                   0.949),
    '2505/3, 3B': (1.006, 1.004, 1.001, 0.999, 0.996, 0.991, 0.987, 0.984, 0.980, 0.976, 0.971, 0.964, 0.957,
                   0.950, 0.938),
    '2571': (1.005, 1.004, 1.002, 1.000, 0.998, 0.995, 0.993, 0.991, 0.989, 0.986, 0.982, 0.975, 0.969, 0.961, 0.949),
    '2581': (1.005, 1.003, 1.001, 0.998, 0.995, 0.991, 0.986, 0.983, 0.980, 0.975, 0.970, 0.963, 0.956, 0.949, 0.937),

    # PTW
    '30001': (1.004, 1.003, 1.001, 0.999, 0.997, 0.994, 0.990, 0.988, 0.985, 0.981, 0.976, 0.969, 0.962, 0.955,
              0.943),
    '30010': (1.004, 1.003, 1.001, 0.999, 0.997, 0.994, 0.990, 0.988, 0.985, 0.981, 0.976, 0.969, 0.962, 0.955,
              0.943),
    '30002': (1.006, 1.004, 1.001, 0.999, 0.997, 0.994, 0.992, 0.990, 0.987, 0.984, 0.980, 0.973, 0.967, 0.959,
              0.948),
    '30011': (1.006, 1.004, 1.001, 0.999, 0.997, 0.994, 0.992, 0.990, 0.987, 0.984, 0.980, 0.973, 0.967, 0.959,
              0.948),
    '30004': (1.006, 1.005, 1.002, 1.000, 0.999, 0.996, 0.994, 0.992, 0.989, 0.986, 0.982, 0.976, 0.969, 0.962,
              0.950),
    '30012': (1.006, 1.005, 1.002, 1.000, 0.999, 0.996, 0.994, 0.992, 0.989, 0.986, 0.982, 0.976, 0.969, 0.962,
              0.950),
    '30006': (1.002, 1.002, 1.000, 0.999, 0.997, 0.994, 0.990, 0.988, 0.984, 0.980, 0.975, 0.968,
              0.960, 0.952, 0.940),
    '30013': (1.002, 1.002, 1.000, 0.999, 0.997, 0.994, 0.990, 0.988, 0.984, 0.980, 0.975, 0.968,
              0.960, 0.952, 0.940),
}

KQ_ELECTRON_R50S = (4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0, 13.0, 16.0, 20.0)

KQ_ELECTRON_CHAMBERS = {

    # Capintec
    'PR06C': (0.916, 0.914, 0.912, 0.911, 0.909, 0.906, 0.904, 0.899, 0.891, 0.884, 0.874),

    # Exradin
    'A12': (0.921, 0.919, 0.918, 0.916, 0.914, 0.911, 0.909, 0.903, 0.896, 0.888, 0.878),

    # NE
    '2571': (0.918, 0.916, 0.915, 0.913, 0.911, 0.909, 0.906, 0.901, 0.893, 0.886, 0.876),
    '2581': (0.899, 0.898, 0.896, 0.894, 0.893, 0.890, 0.888, 0.882, 0.875, 0.868, 0.859),

    # PTW
    '30001': (0.911, 0.909, 0.907, 0.905, 0.904, 0.901, 0.898, 0.893, 0.885, 0.877, 0.868),
    '30010': (0.911, 0.909, 0.907, 0.905, 0.904, 0.901, 0.898, 0.893, 0.885, 0.877, 0.868),
    '30002': (0.916, 0.914, 0.912, 0.910, 0.909, 0.906, 0.903, 0.897, 0.890, 0.882, 0.873),
    '30011': (0.916, 0.914, 0.912, 0.910, 0.909, 0.906, 0.903, 0.897, 0.890, 0.882, 0.873),
    '30004': (0.920, 0.918, 0.916, 0.915, 0.913, 0.910, 0.907, 0.902, 0.894, 0.887, 0.877),
    '30012': (0.920, 0.918, 0.916, 0.915, 0.913, 0.910, 0.907, 0.902, 0.894, 0.887, 0.877),
    '30006': (0.911, 0.909, 0.907, 0.906, 0.904, 0.901, 0.898, 0.893, 0.885, 0.878, 0.868),
    '30013': (0.911, 0.909, 0.907, 0.906, 0.904, 0.901, 0.898, 0.893, 0.885, 0.878, 0.868),

    # Wellhofer
    'FC65-P': (0.914, 0.912, 0.911, 0.909, 0.907, 0.904, 0.902, 0.896, 0.889, 0.881, 0.872),
    'FC65-G': (0.920, 0.918, 0.916, 0.914, 0.913, 0.910, 0.907, 0.902, 0.894, 0.887, 0.877),
}

# Rename common functions from TG-51
k_tp = _tg51.p_tp
k_pol = _tg51.p_pol
z_ref = _tg51.d_ref
r_50 = _tg51.r_50


def k_s(*, voltage_reference: int, voltage_reduced: int,
        m_reference: NumberOrArray, m_reduced: NumberOrArray) -> float:
    """Calculate the ion recombination effect using readings at two voltages. The voltages should have
    a ratio of 2, 2.5, 3, 3.5, 4, or 5.

    Parameters
    ----------
    voltage_reference : int
        The voltage at which calibration will be performed (e.g. 300V)
    voltage_reduced : int
        The voltage which is lower than reference (e.g. 150V)
    m_reference : array, float
        The reading(s) at the reference voltage.
    m_reduced : array, float
        The reading(s) at the reduced voltage.

    Returns
    -------
    k_s : float
        The ion recombination factor.

    Raises
    ------
    ValueError
        If the voltage ratio is not valid.
    ValueError
        If the calculated ks value is outside the range (1.0, 1.05).
    """
    v_ratio = voltage_reference / voltage_reduced
    _verify_voltage_ratio_is_valid(v_ratio)
    a = V1_V2_FITS[v_ratio]
    m_ratio = np.mean(m_reference) / np.mean(m_reduced)
    argue.verify_bounds(m_ratio, bounds=(1.0, 1.05), message="Ks is out of bounds. Verify inputs or check chamber")
    return float(a['a0'] + a['a1']*m_ratio + a['a2']*(m_ratio**2))


def _verify_voltage_ratio_is_valid(voltage_ratio):
    """Helper function to verify the voltage ratio of high/low is one of a valid set as defined by TRS-398 Table 4.VII"""
    if not is_close(voltage_ratio, target=(2, 2.5, 3, 3.5, 4, 5), delta=0.001):
        raise ValueError("voltage_reference and voltage_reduced are not a valid ratio. Valid ratios are: 2, 2.5, 3, 3.5, 4, 5")


@argue.options(chamber=KQ_PHOTON_CHAMBERS.keys())
@argue.bounds(tpr=(KQ_PHOTON_TPRS[0], KQ_PHOTON_TPRS[-1]))
def kq_photon(*, chamber: str, tpr: float) -> float:
    """Calculate the kQ factor for a photon beam given the chamber model and TPR20/10 using Table 6.III.
    Linear interpolation is used between given TPR ratios.

    Parameters
    ----------
    chamber : str
        Allowable chambers are those listed in Table 6.III that are also Farmer-type (e.g. Exradin A14 Farmer).
    tpr : float
        The ratio of measured TPR(20cm) / TPR(10cm). Note that this can also be calculated from PDD. See :func:`~pylinac.calibration.tg51.tpr2010_from_pdd2010`.

    Returns
    -------
    kQ : float
        The calculated kQ given table Table 6.III

    Raises
    ------
    KeyError
        If the passed chamber is not within the acceptable list.
    ValueError
        If the TPR is not within the range defined by Table 6.III
    """
    return np.interp([tpr], KQ_PHOTON_TPRS, KQ_PHOTON_CHAMBERS[chamber])[0]


@argue.options(chamber=KQ_ELECTRON_CHAMBERS.keys())
@argue.bounds(r_50=(KQ_ELECTRON_R50S[0], KQ_ELECTRON_R50S[-1]))
def kq_electron(*, chamber: str, r_50: float) -> float:
    """Calculate the kQ factor for an electron beam given the chamber model and R50 using Table 7.III.
    Linear interpolation is used between given R50 values.

    Parameters
    ----------
    chamber : str
        The Farmer-type chambers listed in Table 7.III (e.g. PTW 30004/30012).
    r_50 : float
        The depth of R50 in cm in water.

    Returns
    -------
    kQ : float
        The calculated kQ from Table 7.III

    Raises
    ------
    KeyError
        If the passed chamber is not within the acceptable list.
    ValueError
        If the R50 is not within the range defined by Table 7.III
    """
    return np.interp([r_50], KQ_ELECTRON_R50S, KQ_ELECTRON_CHAMBERS[chamber])[0]


@argue.bounds(k_tp=(0.9, 1.1), k_elec=(0.95, 1.05), k_pol=(0.95, 1.05), k_s=(1.0, 1.05))
def m_corrected(*, m_reference, k_tp, k_elec, k_pol, k_s) -> float:
    """The fully corrected chamber reading.

    Parameters
    ----------
    m_reference : array, float
        The chamber reading(s) at the calibration position.
    k_tp : float
        Temperature/Pressure correction. See :func:`~pylinac.calibration.tg51.p_tp`.
    k_elec : float
        Electrometer correction; given by the calibration laboratory.
    k_pol : float
        Polarity correction. See :func:`~pylinac.calibration.tg51.p_pol`.
    k_s : float
        Ion recombination correction. See :func:`~pylinac.calibration.trs398.k_s`.

    Returns
    -------
    m : float
        The fully corrected chamber reading.
    """
    return float(np.mean(m_reference) * k_tp * k_elec * k_pol * k_s)


class TRS398Base(Structure):

    @property
    def k_tp(self):
        """Temperature/Pressure correction"""
        return k_tp(temp=self.temp, press=self.press)

    @property
    def k_pol(self):
        """Polarity correction"""
        return k_pol(m_reference=self.m_reference, m_opposite=self.m_opposite)

    @property
    def k_s(self):
        """Ionization collection correction."""
        return k_s(voltage_reference=self.voltage_reference, voltage_reduced=self.voltage_reduced,
                   m_reference=self.m_reference, m_reduced=self.m_reduced)

    @property
    def m_corrected(self):
        """Corrected chamber reading."""
        return m_corrected(m_reference=self.m_reference, k_tp=self.k_tp, k_elec=self.k_elec,
                           k_pol=self.k_pol, k_s=self.k_s)

    @property
    def dose_mu_zref(self):
        """cGy/MU at a zref depth."""
        return self.tissue_correction * self.m_corrected * self.n_dw * self.kq / self.mu

    @property
    def m_corrected_adjusted(self):
        """Corrected chamber reading after adjusting the output."""
        # TODO: add check for m adjusted != None
        return m_corrected(m_reference=self.m_reference_adjusted, k_tp=self.k_tp, k_elec=self.k_elec,
                           k_pol=self.k_pol, k_s=self.k_s)

    @property
    def dose_mu_zref_adjusted(self):
        """The dose/mu at 10cm depth after adjustment."""
        return self.tissue_correction * self.m_corrected_adjusted * self.n_dw * self.kq / self.mu

    @property
    def output_was_adjusted(self):
        """Boolean specifiying if output was adjusted."""
        return self.m_reference_adjusted is not None


class TRS398Photon(TRS398Base):
    """Calculation of dose to water at zmax and zref from a high energy photon beam. Setup can be SSD or SAD.

    Parameters
    ----------
    setup : {'SSD', 'SAD'}
        The physical setup of the calibration.
    institution : str
        Institution name.
    physicist : str
        Physicist performing calibration.
    unit : str
        Unit name; e.g. TrueBeam1.
    measurement_date : str
        Date of measurement. E.g. 10/22/2018.
    chamber : str
        Farmer-type chamber model from Table 6.III.
    n_dw : float
        NDw of the chamber given by the calibration laboratory.
    mu : float, int
        The number of MU given per reading
    energy : int
        Nominal energy of the linac in MV; e.g. 6. Bookkeeping only.
    fff : bool
        Whether the beam is FFF or flat. Bookkeeping only.
    tpr2010 : float
        The value of TPR(20)/TPR(10). Can be derived from PDD; see :func:`~pylinac.calibration.tg51.tpr2010_from_pdd2010`.
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
        The electrometer correction value given by the calibration laboratory.
    clinical_pdd_zref : optional, float
        The PDD at the depth of calibration. Use the actual percentage (e.g. 66.7 not 0.667). If not supplied the clinical_tmr_zref value must be supplied.
    clinical_tmr_zref : optional, float
        The TMR at the depth of calibration. If not supplied the clinical_pdd_zref value must be supplied.
    tissue_correction : float
        The correction of calibration to a medium other than water. Default value is 1 which is water. E.g. use 0.99 if calibrating to muscle.
    """

    @argue.options(setup=('SSD', 'SAD'))
    @argue.options(chamber=KQ_PHOTON_CHAMBERS.keys())
    @argue.bounds(tpr2010=(KQ_PHOTON_TPRS[0], KQ_PHOTON_TPRS[-1]))
    def __init__(self, *,
                 institution: str = '',
                 physicist: str = '',
                 unit: str = '',
                 measurement_date: str = '',
                 electrometer: str = '',
                 setup: str,
                 chamber: str,
                 n_dw: float,
                 mu: int,
                 tpr2010: float,
                 energy: int,
                 fff: bool,
                 press: float,
                 temp: float,
                 voltage_reference: int,
                 voltage_reduced: int,
                 m_reference: Union[tuple, float],
                 m_reduced: Union[tuple, float],
                 m_opposite: Union[tuple, float],
                 k_elec: float,
                 clinical_pdd_zref: Optional[float]=None,
                 clinical_tmr_zref: Optional[float]=None,
                 tissue_correction: float=1.0,
                 ):
        super().__init__(chamber=chamber, tpr2010=tpr2010, press=press, temp=temp,
                         voltage_reference=voltage_reference, voltage_reduced=voltage_reduced, m_reference=m_reference,
                         m_reduced=m_reduced, m_opposite=m_opposite, k_elec=k_elec, clinical_pdd_zref=clinical_pdd_zref,
                         clinical_tmr_zref=clinical_tmr_zref, n_dw=n_dw, setup=setup, mu=mu,
                         tissue_correction=tissue_correction, fff=fff, energy=energy, institution=institution,
                         physicist=physicist, unit=unit, measurement_date=measurement_date, electrometer=electrometer)
        self.m_reference_adjusted = None

    @property
    def kq(self):
        """kQ of the chamber and TPR."""
        return kq_photon(chamber=self.chamber, tpr=self.tpr2010)

    @property
    def dose_mu_zmax(self):
        """cGy/MU at a depth of zmax."""
        if self.setup == 'SSD':
            return (100 * self.dose_mu_zref) / self.clinical_pdd_zref
        else:
            return self.dose_mu_zref / self.clinical_tmr_zref

    @property
    def dose_mu_zmax_adjusted(self):
        """The dose/mu at dmax depth after adjustment."""
        if self.setup == 'SSD':
            return (100 * self.dose_mu_zref_adjusted) / self.clinical_pdd_zref
        else:
            return self.dose_mu_zref_adjusted / self.clinical_tmr_zref

    def publish_pdf(self, filename: str, notes: Optional[list]=None, open_file: bool=False,
                    metadata: Optional[dict] = None):
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
        title = f'TRS-398 Photon Report - {self.energy} MV'
        if self.fff:
            title += ' FFF'

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
            f'Kelec: {self.k_elec:2.3f}',
            f'MU: {self.mu}',
            '',
            'Beam Quality:',
            f'TPR(20)/TPR(10): {self.tpr2010:2.3f}',
            f'Determined kQ: {self.kq:2.3f}',
            '',
            'Chamber Corrections/Measurements:',
            f'Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}',
            f'Pressure (kPa): {self.press:2.1f}',
            f'Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}',
            f'Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}',
            f'Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}',
            f'Ktp: {self.k_tp:2.3f}',
            f'Ks: {self.k_s:2.3f}',
            f'Kpol: {self.k_pol:2.3f}',
            '',
            'Dose Determination:',
            f'Fully corrected M (nC): {self.m_corrected:2.3f}',
            f'Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}',
            f'Dose/MU @ zref depth (cGy): {self.dose_mu_zref:2.3f}',
        ]
        if self.setup == 'SSD':
            text.append(
                f'Clinical PDD (%): {self.clinical_pdd_zref:2.2f}',
            )
        else:
            text.append(
                f'Clinical TMR (%): {self.clinical_tmr_zref:2.2f}',
            )
        text.append(f'Dose/MU @ zmax (cGy): {self.dose_mu_zmax:2.3f}')
        text.append("")
        text.append(f'Output Adjustment?: {was_adjusted}')
        if was_adjusted == 'Yes':
            text.append(f'Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjusted}')
            text.append(f'Adjusted fully corrected M (nC): {self.m_corrected_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ zref depth (cGy): {self.dose_mu_zref_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ zmax (cGy): {self.dose_mu_zmax_adjusted:2.3f}')
        canvas.add_text(text=text, location=(2, 25.5), font_size=12)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            open_path(filename)


class TRS398Electron(TRS398Base):
    """Calculation of dose to water at zmax and zref from a high energy electron beam.

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
    chamber : str
        Farmer-type chamber model from Table 6.III.
    n_dw : float
        NDw of the chamber given by the calibration laboratory.
    mu : float, int
        The number of MU given per reading.
    i_50 : float
        The depth of ionization 50% in cm.
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
        The electrometer correction value given by the calibration laboratory.
    pdd_zref : optional, float
        The PDD at the depth of calibration. Use the actual percentage (e.g. 66.7 not 0.667). If not supplied the tmr_zref value should be supplied.
    tissue_correction : float
        The correction of calibration to a medium other than water. Default value is 1 which is water. E.g. use 0.99 if calibrating to muscle.
    """

    def __init__(self, *,
                 institution: str = '',
                 physicist: str = '',
                 unit: str = '',
                 measurement_date: str = '',
                 electrometer: str = '',
                 energy: str,
                 cone: str,
                 chamber: str,
                 n_dw: float,
                 mu: int,
                 i_50: float,
                 press: float,
                 temp: float,
                 voltage_reference: int,
                 voltage_reduced: int,
                 m_reference: tuple,
                 m_reduced: tuple,
                 m_opposite: tuple,
                 k_elec: float,
                 clinical_pdd_zref: float,
                 tissue_correction: float=1.0,
                 ):
        super().__init__(chamber=chamber, i_50=i_50, press=press, temp=temp, energy=energy, institution=institution,
                         voltage_reference=voltage_reference, voltage_reduced=voltage_reduced, m_reference=m_reference,
                         m_reduced=m_reduced, m_opposite=m_opposite, k_elec=k_elec, clinical_pdd_zref=clinical_pdd_zref,
                         n_dw=n_dw, mu=mu, tissue_correction=tissue_correction, physicist=physicist, unit=unit,
                         measurement_date=measurement_date, electrometer=electrometer, cone=cone)
        self.m_reference_adjusted = None

    @property
    def r_50(self) -> float:
        """The depth of R50 in cm, derived from I50."""
        return r_50(i_50=self.i_50)

    @property
    def zref(self) -> float:
        """Depth of the reference point."""
        return z_ref(i_50=self.i_50)

    @property
    def kq(self):
        """kQ given the chamber and R50."""
        return kq_electron(chamber=self.chamber, r_50=self.r_50)

    @property
    def dose_mu_zmax(self):
        """cGy/MU at a depth of zmax."""
        return (100 * self.dose_mu_zref) / self.clinical_pdd_zref

    @property
    def dose_mu_zmax_adjusted(self):
        """The dose/mu at dmax depth after adjustment."""
        return (100 * self.dose_mu_zref_adjusted) / self.clinical_pdd_zref

    def publish_pdf(self, filename: str, notes: Optional[list]=None, open_file: bool=False,
                    metadata: Optional[dict] = None):
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
        title = 'TRS-398 Electron Report - {} MV'.format(self.energy)

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
            f'Kelec: {self.k_elec:2.3f}',
            '',
            'Beam Quality:',
            f'I50 (cm): {self.i_50:2.2f}',
            f'R50 (cm): {self.r_50:2.2f}',
            f'Zref (cm): {self.zref:2.2f}',
            f"Calculated kQ: {self.kq:2.3f}",
            '',
            'Chamber Corrections/Measurements:',
            f'Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}',
            f'Pressure (kPa): {self.press:2.1f}',
            f'Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}',
            f'Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}',
            f'Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}',
            f'Ktp: {self.k_tp:2.3f}',
            f'Ks: {self.k_s:2.3f}',
            f'Kpol: {self.k_pol:2.3f}',
            '',
            'Dose Determination:',
            f'Fully corrected M (nC): {self.m_corrected:2.3f}',
            f'Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}',
            f'Dose/MU @ zref (cGy): {self.dose_mu_zref:2.3f}',
            f'Clinical PDD (%): {self.clinical_pdd_zref:2.2f}',
            f'Dose/MU @ zmax (cGy): {self.dose_mu_zmax:2.3f}',
            "",
            f'Output Adjustment?: {was_adjusted}',
        ]
        if was_adjusted == 'Yes':
            text.append(f'Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjusted}')
            text.append(f'Adjusted fully corrected M (nC): {self.m_corrected_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ zref depth (cGy): {self.dose_mu_zref_adjusted:2.3f}')
            text.append(f'Adjusted Dose/MU @ zmax (cGy): {self.dose_mu_zmax_adjusted:2.3f}')
        canvas.add_text(text=text, location=(2, 25.5), font_size=12)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            open_path(filename)
