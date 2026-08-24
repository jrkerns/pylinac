"""
The TG-51 module contains a number of helper functions and classes that can calculate parameters for performing the
TG-51 absolute linac dose calibration although there are some modifications from the original TG-51. The modifications
include updated kQ and kecal values from Muir and Rogers' set of papers.
Functions include all relevant calculations for TG-51 including PDDx, kQ,
Dref, and chamber reading corrections. Where Muir & Rogers' values/equations are used they are specified in the documentation.

Classes include photon and electron calibrations using cylindrical chambers. Pass all the relevant raw measurements
and the class will compute all corrections and corrected readings and dose at 10cm and dmax/dref.
"""

from __future__ import annotations

import webbrowser
from abc import abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import cast

import argue
import numpy as np

from ..core.pdf import PylinacCanvas
from ..core.typing import NumberOrArray
from ..core.utilities import Structure, convert_to_enum

MIN_TEMP = 15
MAX_TEMP = 35
MIN_PRESSURE = 90
MAX_PRESSURE = 115
MIN_PION = 1
MAX_PION = 1.05
MIN_PTP = 0.9
MAX_PTP = 1.1
MIN_PELEC = 0.98
MAX_PELEC = 1.02
MIN_PPOL = 0.98
MAX_PPOL = 1.02


class LeadFoil(Enum):
    """Lead foil positions supported by the TG-51 photon calculations."""

    NONE = None
    CM30 = "30cm"
    CM50 = "50cm"

    @classmethod
    def options(cls) -> tuple[LeadFoil | str | None, ...]:
        """Return enum and legacy values for argument validation."""
        return *tuple(cls), *(option.value for option in cls)

    def __str__(self) -> str:
        """Return the historical lead foil value."""
        return "None" if self.value is None else self.value


@dataclass(frozen=True)
class PhotonChamberCoefficients:
    """Photon beam-quality correction coefficients for an ion chamber.

    The unprimed coefficients are used with PDD(10)x. Their ``b`` and ``c``
    values are stored in the units published in Table I of the TG-51 Addendum;
    Equation 1 applies factors of :math:`10^{-3}` and :math:`10^{-5}` to them,
    respectively. The primed coefficients are optional and, when present, are
    used with TPR(20,10) as an alternative to calculate kQ. See Muir & Rogers
    2010 Table III and Equation 10: https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.3495537
    """

    a: float
    b: float
    c: float
    a_prime: float | None = None
    b_prime: float | None = None
    c_prime: float | None = None
    d_prime: float | None = None

    def __post_init__(self) -> None:
        """Require TPR coefficients to be provided as a complete set."""
        tpr_coefficients = (
            self.a_prime,
            self.b_prime,
            self.c_prime,
            self.d_prime,
        )
        if any(value is None for value in tpr_coefficients) and not all(
            value is None for value in tpr_coefficients
        ):
            raise ValueError(
                "TPR coefficients a_prime, b_prime, c_prime, and d_prime "
                "must all be provided or all be omitted"
            )

    def as_dict(self) -> dict[str, float | None]:
        """Return the coefficients using the historical dictionary keys."""
        return {
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "a'": self.a_prime,
            "b'": self.b_prime,
            "c'": self.c_prime,
            "d'": self.d_prime,
        }


@dataclass(frozen=True)
class ElectronChamberCoefficients:
    """Electron beam-quality correction coefficients for an ion chamber."""

    kq_ecal: float
    a: float
    b: float
    c: float

    def as_dict(self) -> dict[str, float]:
        """Return the coefficients using the historical dictionary keys."""
        return {
            "kQ,ecal": self.kq_ecal,
            "a": self.a,
            "b": self.b,
            "c": self.c,
        }


@dataclass(frozen=True)
class IonChamber:
    """An ion chamber and its available TG-51 correction coefficients.

    Predefined chambers are from the TG-51 Addendum (McEwen, 2014).
    """

    name: str
    photon_coefficients: PhotonChamberCoefficients | None
    electron_coefficients: ElectronChamberCoefficients | None
    # The vast majority of kQ factors are for the range 63-86, per McEwen Table I.
    pddx_range: tuple[float, float] = (63.0, 86.0)

    def __post_init__(self) -> None:
        """Validate the supported PDD(10)x range."""
        if len(self.pddx_range) != 2:
            raise ValueError("pddx_range must contain exactly two bounds")
        lower, upper = self.pddx_range
        if not isinstance(lower, (int, float)) or not isinstance(upper, (int, float)):
            raise TypeError("pddx_range bounds must be numeric")
        if lower >= upper:
            raise ValueError("pddx_range lower bound must be less than upper bound")

    def __str__(self) -> str:
        """Return the historical chamber model name."""
        return self.name


class IonChambers:
    """Predefined ion chambers supported by the TG-51 calculations."""

    def __iter__(self) -> Iterator[IonChamber]:
        """Iterate over the predefined ion chambers."""
        for chamber in vars(type(self)).values():
            if isinstance(chamber, IonChamber):
                yield chamber

    @classmethod
    def from_name(cls, name: str) -> IonChamber:
        """Retrieve a chamber by its historical model name.

        Predefined chambers are returned directly. Chambers added dynamically to
        :data:`KQ_PHOTONS` or :data:`KQ_ELECTRONS` are converted to an
        :class:`IonChamber` with frozen coefficient objects. This is for historical compatibility
        for users who modified ``KQ_PHOTONS`` or ``KQ_ELECTRONS`` on the fly in previous version.

        Parameters
        ----------
        name
            The chamber model name, such as ``"30013"`` or ``"A12"``.

        Raises
        ------
        ValueError
            If no ion chamber has the given name or a custom coefficient table is
            missing a required key.
        """
        # predefined list
        for chamber in cls():
            if chamber.name == name:
                return chamber

        # for historical compatibility, construct any chamber that was added to
        # KQ_PHOTONS/ELECTRONS on the fly and check them.
        photon_coefficients = cls._photon_coefficients_from_table(name)
        electron_coefficients = cls._electron_coefficients_from_table(name)
        if photon_coefficients is not None or electron_coefficients is not None:
            return IonChamber(
                name=name,
                photon_coefficients=photon_coefficients,
                electron_coefficients=electron_coefficients,
            )

        valid_names = ", ".join(
            sorted(
                {chamber.name for chamber in cls()}
                | KQ_PHOTONS.keys()
                | KQ_ELECTRONS.keys()
            )
        )
        raise ValueError(
            f"No matching ion chamber exists for the name {name!r}. "
            f"Valid ion chamber names are: {valid_names}"
        )

    @staticmethod
    def _photon_coefficients_from_table(
        name: str,
    ) -> PhotonChamberCoefficients | None:
        """Historical compatibility shim."""
        coefficients = KQ_PHOTONS.get(name)
        if coefficients is None:
            return None
        try:
            a = coefficients["a"]
            b = coefficients["b"]
            c = coefficients["c"]
        except KeyError as error:
            raise ValueError(
                f"Photon coefficients for ion chamber {name!r} are missing "
                f"required key {error.args[0]!r}"
            ) from error

        tpr_keys = ("a'", "b'", "c'", "d'")
        present_tpr_keys = tuple(key for key in tpr_keys if key in coefficients)
        if present_tpr_keys and len(present_tpr_keys) != len(tpr_keys):
            missing_key = next(key for key in tpr_keys if key not in coefficients)
            raise ValueError(
                f"Photon coefficients for ion chamber {name!r} are missing "
                f"required key {missing_key!r}"
            )

        return PhotonChamberCoefficients(
            a=a,
            b=b,
            c=c,
            a_prime=coefficients.get("a'"),
            b_prime=coefficients.get("b'"),
            c_prime=coefficients.get("c'"),
            d_prime=coefficients.get("d'"),
        )

    @staticmethod
    def _electron_coefficients_from_table(
        name: str,
    ) -> ElectronChamberCoefficients | None:
        """Historical compatibility shim."""
        coefficients = KQ_ELECTRONS.get(name)
        if coefficients is None:
            return None
        try:
            return ElectronChamberCoefficients(
                kq_ecal=coefficients["kQ,ecal"],
                a=coefficients["a"],
                b=coefficients["b"],
                c=coefficients["c"],
            )
        except KeyError as error:
            raise ValueError(
                f"Electron coefficients for ion chamber {name!r} are missing "
                f"required key {error.args[0]!r}"
            ) from error

    # Exradin
    EXRADIN_A12 = IonChamber(
        name="A12",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0146,
            b=0.777,
            c=-1.666,
            a_prime=2.6402,
            b_prime=-7.2304,
            c_prime=10.7573,
            d_prime=-5.4294,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.907, a=0.965, b=0.119, c=0.607
        ),
    )
    EXRADIN_A19 = IonChamber(
        name="A19",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9934,
            b=1.384,
            c=-2.125,
            a_prime=3.0907,
            b_prime=-9.193,
            c_prime=13.5957,
            d_prime=-6.7969,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.957, b=0.119, c=0.505
        ),
    )
    EXRADIN_A2 = IonChamber(
        name="A2",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9819,
            b=1.609,
            c=-2.184,
            a_prime=2.8458,
            b_prime=-8.1619,
            c_prime=12.1411,
            d_prime=-6.1041,
        ),
        electron_coefficients=None,
    )
    EXRADIN_T2 = IonChamber(
        name="T2",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0173,
            b=0.854,
            c=-1.941,
            a_prime=3.3433,
            b_prime=-10.2649,
            c_prime=15.1247,
            d_prime=-7.5415,
        ),
        electron_coefficients=None,
    )
    EXRADIN_A12S = IonChamber(
        name="A12S",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9692,
            b=1.974,
            c=-2.448,
            a_prime=2.9597,
            b_prime=-8.6777,
            c_prime=12.9155,
            d_prime=-6.4903,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.907, a=0.937, b=0.136, c=0.378
        ),
    )
    EXRADIN_A18 = IonChamber(
        name="A18",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9944,
            b=1.286,
            c=-1.980,
            a_prime=2.5167,
            b_prime=-6.7567,
            c_prime=10.1519,
            d_prime=-5.1709,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.914, a=0.352, b=0.711, c=0.046
        ),
    )
    EXRADIN_A1 = IonChamber(
        name="A1",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0029,
            b=1.023,
            c=-1.803,
            a_prime=2.0848,
            b_prime=-4.9174,
            c_prime=7.5446,
            d_prime=-3.9441,
        ),
        electron_coefficients=None,
    )
    EXRADIN_T1 = IonChamber(
        name="T1",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0552,
            b=-0.196,
            c=-1.275,
            a_prime=2.806,
            b_prime=-7.9273,
            c_prime=11.7541,
            d_prime=-5.9263,
        ),
        electron_coefficients=None,
    )
    EXRADIN_A1SL = IonChamber(
        name="A1SL",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9896,
            b=1.410,
            c=-2.049,
            a_prime=2.8029,
            b_prime=-7.9648,
            c_prime=11.8445,
            d_prime=-5.9568,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.914, a=0.205, b=0.854, c=0.036
        ),
    )
    EXRADIN_A14 = IonChamber(
        name="A14",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9285,
            b=2.706,
            c=-2.599,
            a_prime=5.4677,
            b_prime=-19.1795,
            c_prime=27.4542,
            d_prime=-13.1336,
        ),
        electron_coefficients=None,
    )
    EXRADIN_T14 = IonChamber(
        name="T14",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9622,
            b=2.009,
            c=-2.401,
            a_prime=4.969,
            b_prime=-17.1074,
            c_prime=24.6292,
            d_prime=-11.8877,
        ),
        electron_coefficients=None,
    )
    EXRADIN_A14SL = IonChamber(
        name="A14SL",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9017,
            b=3.454,
            c=-3.083,
            a_prime=5.1205,
            b_prime=-17.7884,
            c_prime=25.6123,
            d_prime=-12.3232,
        ),
        electron_coefficients=None,
    )
    EXRADIN_A16 = IonChamber(
        name="A16",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.8367,
            b=4.987,
            c=-3.877,
            a_prime=6.0571,
            b_prime=-21.7829,
            c_prime=31.2289,
            d_prime=-14.9168,
        ),
        electron_coefficients=None,
    )

    # PTW
    PTW_30010 = IonChamber(
        name="30010",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0093,
            b=0.926,
            c=-1.771,
            a_prime=2.5318,
            b_prime=-6.7948,
            c_prime=10.1779,
            d_prime=-5.1746,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.98, b=0.119, c=0.891
        ),
    )
    PTW_30011 = IonChamber(
        name="30011",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9676,
            b=2.061,
            c=-2.528,
            a_prime=2.9044,
            b_prime=-8.4576,
            c_prime=12.6339,
            d_prime=-6.3742,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.901, a=0.976, b=0.12, c=0.793
        ),
    )
    PTW_30012 = IonChamber(
        name="30012",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9537,
            b=2.440,
            c=-2.750,
            a_prime=3.2836,
            b_prime=-10.061,
            c_prime=14.8867,
            d_prime=-7.4212,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.908, a=0.972, b=0.121, c=0.728
        ),
    )
    PTW_30013 = IonChamber(
        name="30013",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9652,
            b=2.141,
            c=-2.623,
            a_prime=3.2012,
            b_prime=-9.7211,
            c_prime=14.4211,
            d_prime=-7.2184,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.901, a=0.978, b=0.112, c=0.816
        ),
    )
    PTW_31010 = IonChamber(
        name="31010",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.959,
            b=2.265,
            c=-2.684,
            a_prime=3.1578,
            b_prime=-9.5422,
            c_prime=14.1676,
            d_prime=-7.0964,
        ),
        electron_coefficients=None,
    )
    # Per TG-51A1 3.E.i 31003 is identical to 31013.
    PTW_31003 = IonChamber(
        name="31003",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9725,
            b=1.957,
            c=-2.498,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.902, a=0.945, b=0.133, c=0.441
        ),
    )
    PTW_31013 = IonChamber(
        name="31013",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9725,
            b=1.957,
            c=-2.498,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.902, a=0.945, b=0.133, c=0.441
        ),
    )
    PTW_31014 = IonChamber(
        name="31014",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0071,
            b=1.048,
            c=-1.967,
            a_prime=3.0178,
            b_prime=-8.8735,
            c_prime=13.1372,
            d_prime=-6.5867,
        ),
        electron_coefficients=None,
    )
    PTW_31016 = IonChamber(
        name="31016",
        photon_coefficients=PhotonChamberCoefficients(
            a=1.0085,
            b=1.028,
            c=-1.968,
            a_prime=2.9524,
            b_prime=-8.6054,
            c_prime=12.7757,
            d_prime=-6.4265,
        ),
        electron_coefficients=None,
    )

    # IBA
    IBA_CC25 = IonChamber(
        name="CC25",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9551,
            b=2.353,
            c=-2.687,
            a_prime=2.4567,
            b_prime=-6.5932,
            c_prime=10.0471,
            d_prime=-5.1775,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.964, b=0.105, c=0.539
        ),
    )
    IBA_CC13 = IonChamber(
        name="CC13",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9515,
            b=2.455,
            c=-2.768,
            a_prime=3.1982,
            b_prime=-9.7182,
            c_prime=14.421,
            d_prime=-7.2121,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.926, b=0.129, c=0.279
        ),
    )
    # Per TG-51A1 3.E.ii the IC10 and CC13 are functionally identical
    IBA_IC10 = IonChamber(
        name="IC10",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9515,
            b=2.455,
            c=-2.768,
            a_prime=3.1982,
            b_prime=-9.7182,
            c_prime=14.421,
            d_prime=-7.2121,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.926, b=0.129, c=0.279
        ),
    )
    IBA_CC08 = IonChamber(
        name="CC08",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.943,
            b=2.637,
            c=-2.884,
            a_prime=3.7328,
            b_prime=-11.98,
            c_prime=17.5884,
            d_prime=-8.6843,
        ),
        electron_coefficients=None,
    )
    IBA_CC04 = IonChamber(
        name="CC04",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9714,
            b=1.938,
            c=-2.432,
            a_prime=3.0054,
            b_prime=-8.8633,
            c_prime=13.1704,
            d_prime=-6.6075,
        ),
        electron_coefficients=None,
    )
    IBA_CC01 = IonChamber(
        name="CC01",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9116,
            b=3.358,
            c=-3.177,
            a_prime=4.3376,
            b_prime=-14.4935,
            c_prime=21.0293,
            d_prime=-10.2208,
        ),
        electron_coefficients=None,
    )
    IBA_FC65_G = IonChamber(
        name="FC65-G",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9708,
            b=1.972,
            c=-2.480,
            a_prime=3.3221,
            b_prime=-10.2012,
            c_prime=15.0497,
            d_prime=-7.4872,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.971, b=0.113, c=0.68
        ),
    )
    IBA_FC65_P = IonChamber(
        name="FC65-P",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9828,
            b=1.664,
            c=-2.296,
            a_prime=3.0872,
            b_prime=-9.1919,
            c_prime=13.6137,
            d_prime=-6.8118,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.902, a=0.973, b=0.11, c=0.692
        ),
    )
    IBA_FC23_C = IonChamber(
        name="FC23-C",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.982,
            b=1.579,
            c=-2.166,
            a_prime=3.0511,
            b_prime=-9.0243,
            c_prime=13.3378,
            d_prime=-6.6559,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.904, a=0.971, b=0.097, c=0.591
        ),
    )

    # Nuclear Enterprises
    # Per TG-51A1 3.E.v the NE2581 is no longer recommended
    NE_2571 = IonChamber(
        name="NE2571",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9882,
            b=1.486,
            c=-2.140,
            a_prime=2.2328,
            b_prime=-5.5779,
            c_prime=8.5325,
            d_prime=-4.4352,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.903, a=0.977, b=0.117, c=0.817
        ),
    )
    NE_2561 = IonChamber(
        name="NE2561",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9722,
            b=1.977,
            c=-2.463,
            a_prime=2.4235,
            b_prime=-6.3179,
            c_prime=9.4737,
            d_prime=-4.8307,
        ),
        electron_coefficients=None,
    )
    NE_2611 = IonChamber(
        name="NE2611",
        photon_coefficients=None,
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.896, a=0.979, b=0.12, c=0.875
        ),
    )

    # Capintec
    CAPINTEC_PR06C_G = IonChamber(
        name="PR06C/G",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9519,
            b=2.432,
            c=-2.704,
            a_prime=2.911,
            b_prime=-8.4916,
            c_prime=12.6817,
            d_prime=-6.3874,
        ),
        electron_coefficients=ElectronChamberCoefficients(
            kq_ecal=0.906, a=0.972, b=0.122, c=0.729
        ),
    )

    # Sun Nuclear
    # https://www.sunnuclear.com/uploads/publications/Report-IRS-2066-Monte-Carlo-calculation-of-the-kQ-quality-conversion-factor-for-the-SNC600c-ionization-chamber-for-photon-beam-reference-dosimetry.pdf
    SNC_600C = IonChamber(
        name="SNC600c",
        photon_coefficients=PhotonChamberCoefficients(
            a=0.9992,
            b=1.254,
            c=-2.070,
        ),
        electron_coefficients=None,
    )


# Reconstructed dicts for historical compatibility in the event
# users modify/add to this on the fly in their scripts.
KQ_PHOTONS = {
    chamber.name: chamber.photon_coefficients.as_dict()
    for chamber in IonChambers()
    if (coefficients := chamber.photon_coefficients) is not None
}

KQ_ELECTRONS = {
    chamber.name: chamber.electron_coefficients.as_dict()
    for chamber in IonChambers()
    if (coefficients := chamber.electron_coefficients) is not None
}


def mmHg2kPa(mmHg: float) -> float:
    """Utility function to convert from mmHg to kPa."""
    return mmHg * 101.33 / 760


def mbar2kPa(mbar: float) -> float:
    """Utility function to convert from millibars to kPa."""
    return mbar / 10


def fahrenheit2celsius(f: float) -> float:
    """Utility function to convert from Fahrenheit to Celsius."""
    return (f - 32) * 5 / 9


@argue.bounds(pdd2010=(0.5, 1))
def tpr2010_from_pdd2010(*, pdd2010: float) -> float:
    """Calculate TPR20,10 from PDD20,10. From TRS-398 footnote 25, section 6.3.1, p.68 (https://www-pub.iaea.org/MTCD/Publications/PDF/TRS398_scr.pdf),
    and Followill et al 1998 eqn 1."""
    return 1.2661 * pdd2010 - 0.0595


def p_tp(*, temp: float, press: float) -> float:
    """Calculate the temperature & pressure correction.

    Parameters
    ----------
    temp : float (17-27)
        The temperature in degrees Celsius.
    press : float (91-111)
        The value of pressure in kPa. Can be converted from mmHg and mbar;
        see :func:`~pylinac.calibration.tg51.mmHg2kPa` and :func:`~pylinac.calibration.tg51.mbar2kPa`.
    """
    argue.verify_bounds(
        temp,
        bounds=(MIN_TEMP, MAX_TEMP),
        message="Temperature {:2.2f} out of range. Did you use Fahrenheit? Consider using the utility function fahrenheit2celsius()",
    )
    argue.verify_bounds(
        press,
        bounds=(MIN_PRESSURE, MAX_PRESSURE),
        message="Pressure {:2.2f} out of range. Did you use kPa? Consider using the utility functions mmHg2kPa() or mbar2kPa()",
    )
    return ((273.2 + temp) / 295.2) * (101.33 / press)


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
    polarity = (abs(mref_avg) + abs(mopp_avg)) / abs(2 * mref_avg)
    argue.verify_bounds(
        polarity,
        bounds=(MIN_PPOL, MAX_PPOL),
        message="Polarity correction {:2.2f} out of range (+/-2%). Verify inputs",
    )
    return float(polarity)


def p_ion(
    *,
    voltage_reference: int,
    voltage_reduced: int,
    m_reference: NumberOrArray,
    m_reduced: NumberOrArray,
) -> float:
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
    ion = (1 - voltage_reference / voltage_reduced) / (
        np.mean(m_reference) / np.mean(m_reduced) - voltage_reference / voltage_reduced
    )
    argue.verify_bounds(
        ion,
        bounds=(MIN_PION, MAX_PION),
        message="Pion out of range (1.00-1.05). Check inputs or chamber",
    )
    return float(ion)


def d_ref(*, i_50: float) -> float:
    """Calculate the dref of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    argue.verify_bounds(i_50, bounds=argue.POSITIVE, message="i50 should be positive")
    r50 = r_50(i_50=i_50)
    return 0.6 * r50 - 0.1


def r_50(*, i_50: float) -> float:
    """Calculate the R50 depth of an electron beam based on the I50 depth.

    Parameters
    ----------
    i_50 : float
        The value of I50 in cm.
    """
    argue.verify_bounds(i_50, bounds=argue.POSITIVE, message="i50 should be positive")
    if i_50 < 10:
        r50 = 1.029 * i_50 - 0.06
    else:
        r50 = 1.59 * i_50 - 0.37
    return r50


def kp_r50(*, r_50: float) -> float:
    """Calculate k'R50 for Farmer-like chambers.

    Parameters
    ----------
    r_50 : float (2-9)
        The R50 value in cm.
    """
    argue.verify_bounds(r_50, bounds=(2, 9))
    return 0.9905 + 0.071 * np.exp(-r_50 / 3.67)


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


def m_corrected(
    *,
    p_ion: float,
    p_tp: float,
    p_elec: float,
    p_pol: float,
    m_reference: NumberOrArray,
) -> float:
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
    argue.verify_bounds(p_ion, bounds=(MIN_PION, MAX_PION))
    argue.verify_bounds(p_tp, bounds=(MIN_PTP, MAX_PTP))
    argue.verify_bounds(p_elec, bounds=(MIN_PELEC, MAX_PELEC))
    argue.verify_bounds(p_pol, bounds=(MIN_PPOL, MAX_PPOL))
    return float(p_ion * p_tp * p_elec * p_pol * np.mean(m_reference))


@argue.bounds(pdd=(62.7, 89.0))
@argue.options(lead_foil=LeadFoil.options())
def pddx(*, pdd: float, energy: int, lead_foil: str | LeadFoil | None = None) -> float:
    """Calculate PDDx based on the PDD.

    Parameters
    ----------
    pdd : {>62.7, <89.0}
        The measured PDD. If lead foil was used, this assumes the pdd as measured with the lead in place.
    energy : int
        The nominal energy in MV.
    lead_foil : LeadFoil, {None, '30cm', '50cm'}
        Applicable only for energies >10MV.
        Whether a lead foil was used to acquire the pdd.
        Use ``None`` if no lead foil was used and the interim equation should be used. This is the default
        Use ``50cm`` if the lead foil was set to 50cm from the phantom surface.
        Use ``30cm`` if the lead foil was set to 30cm from the phantom surface.
    """
    lead_foil = cast(LeadFoil, convert_to_enum(lead_foil, LeadFoil))
    if energy < 10:
        return pdd
    elif energy >= 10:
        if lead_foil is LeadFoil.NONE:
            if pdd <= 75:
                return pdd
            elif 75 < pdd <= 89:
                return 1.267 * pdd - 20
            else:
                raise ValueError(f"PDD value of {pdd} was outside the bound of 89%")
        elif lead_foil is LeadFoil.CM50:
            if pdd < 73:
                return pdd
            else:
                return (0.8905 + 0.0015 * pdd) * pdd
        elif lead_foil is LeadFoil.CM30:
            if pdd < 71:
                return pdd
            else:
                return (0.8116 + 0.00264 * pdd) * pdd


def kq_photon_pddx(*, chamber: str | IonChamber, pddx: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements of PDD(10)x. This will calculate kQ for photons
    for *CYLINDRICAL* chambers only.

    Parameters
    ----------
    chamber : str or IonChamber
        The chamber of the chamber. Valid values are those listed in
        Table III of Muir and Rogers and Table I of the TG-51 Addendum.
    pddx : float
        The **PHOTON-ONLY** PDD measurement at 10cm depth for a 10x10cm2 field.
        The allowed range is defined by the selected chamber per TG-51 Addendum. The default range is
        63-86; chambers with a different validated range expose it through
        :attr:`~pylinac.calibration.tg51.IonChamber.pddx_range`.

        .. note:: Use the :func:`~pylinac.calibration.tg51.pddx` function to convert PDD to PDDx as needed.

    Raises
    ------
    argue.BoundsError
        If ``pddx`` is outside the selected chamber's supported range.
    """
    if isinstance(chamber, str):
        chamber = IonChambers.from_name(chamber)
    argue.verify_bounds(pddx, bounds=chamber.pddx_range)
    coef = chamber.photon_coefficients
    if coef is None:
        raise ValueError(f"Chamber {chamber} does not have photon coefficients")
    return coef.a + coef.b * 10**-3 * pddx + coef.c * 10**-5 * pddx**2


@argue.bounds(tpr=(0.623, 0.805))
def kq_photon_tpr(*, chamber: str | IonChamber, tpr: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements of TPR20,10. This will calculate kQ for photons
    for *CYLINDRICAL* chambers only.

    Notes
    -----
    This calculates kQ based on Muir & Rogers, 2010 Equation 10 and table III.
    https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.3495537

    Parameters
    ----------
    chamber : str or IonChamber
        The chamber of the chamber.
    tpr : {>=0.623, <=0.805}
        The TPR(20,10) value.

        .. note::
         Use the :func:`~pylinac.calibration.tg51.tpr2010_from_pdd2010` function to convert from PDD without needing to take TPR measurements.

    Raises
    ------
    ValueError
        If the selected chamber does not have a published TPR(20,10) fit.
    """
    if isinstance(chamber, str):
        chamber = IonChambers.from_name(chamber)
    coef = chamber.photon_coefficients
    if (
        coef is None
        or coef.a_prime is None
        or coef.b_prime is None
        or coef.c_prime is None
        or coef.d_prime is None
    ):
        raise ValueError(f"Chamber {chamber} does not have TPR photon coefficients")
    return (
        coef.a_prime
        + coef.b_prime * tpr
        + coef.c_prime * (tpr**2)
        + coef.d_prime * (tpr**3)
    )


def kq_electron(*, chamber: str | IonChamber, r_50: float) -> float:
    """Calculate kQ based on the chamber and clinical measurements. This will calculate kQ for electrons
    for *CYLINDRICAL* chambers only according to Muir & Rogers.

    Parameters
    ----------
    chamber : str or IonChamber
        The chamber of the chamber. Valid values are those listed in
        Tables VI and VII of Muir and Rogers 2014.
    r_50 : float
        The R50 value in cm of an electron beam.
    """
    if isinstance(chamber, str):
        chamber = IonChambers.from_name(chamber)
    coef = chamber.electron_coefficients
    if coef is None:
        raise ValueError(f"Chamber {chamber} does not have electron coefficients")
    return (coef.a + coef.b * r_50**-coef.c) * coef.kq_ecal


class TG51Base(Structure):
    institution: str
    physicist: str
    unit: str
    measurement_date: str
    temp: float
    press: float
    chamber: str | IonChamber
    n_dw: float
    p_elec: float
    electrometer: str
    energy: int
    voltage_reference: int
    voltage_reduced: int
    m_reference: NumberOrArray
    m_opposite: NumberOrArray
    m_reduced: NumberOrArray
    mu: int
    tissue_correction: float
    m_reference_adjusted: NumberOrArray | None = None

    @property
    def p_tp(self) -> float:
        """Temperature/Pressure correction."""
        return p_tp(temp=self.temp, press=self.press)

    @property
    def p_ion(self) -> float:
        """Ionization collection correction."""
        return p_ion(
            voltage_reference=self.voltage_reference,
            voltage_reduced=self.voltage_reduced,
            m_reference=self.m_reference,
            m_reduced=self.m_reduced,
        )

    @property
    def p_pol(self) -> float:
        """Polarity correction."""
        return p_pol(m_reference=self.m_reference, m_opposite=self.m_opposite)

    @property
    def m_corrected(self) -> float:
        """Corrected chamber reading."""
        return m_corrected(
            p_ion=self.p_ion,
            p_tp=self.p_tp,
            p_elec=self.p_elec,
            p_pol=self.p_pol,
            m_reference=self.m_reference,
        )

    @property
    def m_corrected_adjustment(self) -> float:
        """Corrected chamber reading after adjusting the output."""
        if self.m_reference_adjusted is not None:
            return m_corrected(
                p_ion=self.p_ion,
                p_tp=self.p_tp,
                p_elec=self.p_elec,
                p_pol=self.p_pol,
                m_reference=self.m_reference_adjusted,
            )

    @property
    def output_was_adjusted(self) -> float:
        """Boolean specifiying if output was adjusted."""
        return self.m_reference_adjusted is not None

    @abstractmethod
    def publish_pdf(self, *args, **kwargs):
        pass


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
    chamber : str or IonChamber
        Chamber model. Must be one of the listed chambers in TG-51 Addendum.
    n_dw : float
        NDW value in Gy/nC.
    p_elec : float
        Electrometer correction factor; given by the calibration laboratory.
    measured_pdd10 : float
        The measured value of PDD(10); will be converted to PDDx(10) and used for calculating kq.
    lead_foil : LeadFoil, {None, '50cm', '30cm'}
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

    fff: bool
    measured_pdd10: float | None
    clinical_pdd10: float
    lead_foil: LeadFoil

    @argue.options(lead_foil=LeadFoil.options())
    def __init__(
        self,
        *,
        institution: str = "",
        physicist: str = "",
        unit: str,
        measurement_date: str = "",
        temp: float,
        press: float,
        chamber: str | IonChamber,
        n_dw: float,
        p_elec: float,
        electrometer: str = "",
        measured_pdd10: float | None = None,
        lead_foil: str | LeadFoil | None = None,
        clinical_pdd10: float,
        energy: int,
        fff: bool = False,
        voltage_reference: int,
        voltage_reduced: int,
        m_reference: NumberOrArray,
        m_opposite: NumberOrArray,
        m_reduced: NumberOrArray,
        mu: int,
        tissue_correction: float = 1.0,
        m_reference_adjusted: NumberOrArray | None = None,
    ):
        super().__init__(
            temp=temp,
            press=press,
            chamber=(
                IonChambers.from_name(chamber) if isinstance(chamber, str) else chamber
            ),
            n_dw=n_dw,
            p_elec=p_elec,
            measured_pdd10=measured_pdd10,
            energy=energy,
            voltage_reference=voltage_reference,
            voltage_reduced=voltage_reduced,
            m_reference=m_reference,
            m_opposite=m_opposite,
            m_reduced=m_reduced,
            clinical_pdd10=clinical_pdd10,
            mu=mu,
            tissue_correction=tissue_correction,
            lead_foil=cast(LeadFoil, convert_to_enum(lead_foil, LeadFoil)),
            electrometer=electrometer,
            m_reference_adjusted=m_reference_adjusted,
            institution=institution,
            physicist=physicist,
            unit=unit,
            measurement_date=measurement_date,
            fff=fff,
        )
        # add check for tpr vs pdd

    @property
    def pddx(self) -> float:
        """The photon-only PDD(10) value."""
        return pddx(
            pdd=self.measured_pdd10, energy=self.energy, lead_foil=self.lead_foil
        )

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
        return (
            self.tissue_correction
            * self.m_corrected_adjustment
            * self.kq
            * self.n_dw
            / self.mu
        )

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """The dose/mu at dmax depth after adjustment."""
        return self.dose_mu_10_adjusted / (self.clinical_pdd10 / 100)

    def publish_pdf(
        self,
        filename: str,
        notes: list | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
    ):
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
        was_adjusted = "Yes" if self.output_was_adjusted else "No"
        title = [
            "TG-51 Photon Report",
            f"{self.unit} - {self.energy} MV{' FFF' if self.fff else ''}",
        ]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            "Site Data:",
            f"Institution: {self.institution}",
            f"Performed by: {self.physicist}",
            f"Measurement Date: {self.measurement_date}",
            f"Date of Report: {datetime.now().strftime('%A, %B %d, %Y')}",
            f"Unit: {self.unit}",
            f"Energy: {self.energy} MV {'FFF' if self.fff else ''}",
            "",
            "Instrumentation:",
            f"Chamber: {self.chamber}",
            f"Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}",
            f"Electrometer: {self.electrometer}",
            f"Pelec: {self.p_elec:2.3f}",
            f"MU: {self.mu}",
            "",
            "Beam Quality:",
            f"Lead foil: {'No' if self.lead_foil is LeadFoil.NONE else self.lead_foil}",
            f"Measured PDD(10){'' if self.lead_foil is LeadFoil.NONE else 'Pb'} {self.measured_pdd10:2.2f}",
            f"Calculated PDD(10)x: {self.pddx:2.2f}",
            f"Determined kQ: {self.kq:2.3f}",
            "",
            "Chamber Corrections/Measurements:",
            f"Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}",
            f"Pressure (kPa): {self.press:2.1f}",
            f"Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}",
            f"Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}",
            f"Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}",
            f"Ptp: {self.p_tp:2.3f}",
            f"Pion: {self.p_ion:2.3f}",
            f"Ppol: {self.p_pol:2.3f}",
            "",
            "Dose Determination:",
            f"Fully corrected M (nC): {self.m_corrected:2.3f}",
            f"Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}",
            f"Dose/MU @ 10cm depth (cGy): {self.dose_mu_10:2.3f}",
            f"Clinical PDD (%): {self.clinical_pdd10:2.2f}",
            f"Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}",
            "",
            f"Output Adjustment?: {was_adjusted}",
        ]
        if was_adjusted == "Yes":
            text.append(
                f"Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjusted}"
            )
            text.append(
                f"Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ 10cm depth (cGy): {self.dose_mu_10_adjusted:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}"
            )
        canvas.add_text(text=text, location=(2, 25.5), font_size=12)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            webbrowser.open(filename)


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
    chamber : str or IonChamber
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

    m_gradient: NumberOrArray
    cone: str
    clinical_pdd: float
    i_50: float
    k_ecal: float

    def __init__(
        self,
        *,
        institution: str = "",
        physicist: str = "",
        unit: str = "",
        measurement_date: str = "",
        energy: int,
        temp: float,
        press: float,
        chamber: str | IonChamber,
        k_ecal: float,
        n_dw: float,
        electrometer: str = "",
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
        tissue_correction: float = 1.0,
        m_reference_adjusted=None,
    ):
        if isinstance(chamber, str):
            try:
                chamber = IonChambers.from_name(chamber)
            except ValueError:
                pass  # legacy calibrations permit arbitrary bookkeeping labels
        super().__init__(
            temp=temp,
            press=press,
            chamber=chamber,
            n_dw=n_dw,
            p_elec=p_elec,
            voltage_reference=voltage_reference,
            voltage_reduced=voltage_reduced,
            m_reference=m_reference,
            m_opposite=m_opposite,
            m_reduced=m_reduced,
            clinical_pdd=clinical_pdd,
            mu=mu,
            i_50=i_50,
            tissue_correction=tissue_correction,
            institution=institution,
            physicist=physicist,
            unit=unit,
            measurement_date=measurement_date,
            electrometer=electrometer,
            m_reference_adjusted=m_reference_adjusted,
            cone=cone,
            energy=energy,
            k_ecal=k_ecal,
            m_gradient=m_gradient,
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
        return (
            self.tissue_correction
            * self.m_corrected
            * self.kq
            * self.n_dw
            * self.pq_gr
            / self.mu
        )

    @property
    def dose_mu_dmax(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref / (self.clinical_pdd / 100)

    @property
    def dose_mu_dref_adjusted(self) -> float:
        """cGy/MU at the depth of Dref."""
        return (
            self.tissue_correction
            * self.m_corrected_adjustment
            * self.kq
            * self.n_dw
            * self.pq_gr
            / self.mu
        )

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref_adjusted / (self.clinical_pdd / 100)

    def publish_pdf(
        self,
        filename: str,
        notes: list | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
    ):
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
        was_adjusted = "Yes" if self.output_was_adjusted else "No"
        title = ["TG-51 Electron Report (Legacy)", f"{self.unit} - {self.energy} MeV"]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            "Site Data:",
            f"Institution: {self.institution}",
            f"Performed by: {self.physicist}",
            f"Measurement Date: {self.measurement_date}",
            f"Date of Report: {datetime.now().strftime('%A, %B %d, %Y')}",
            f"Unit: {self.unit}",
            f"Energy: {self.energy} MeV",
            f"Cone: {self.cone}",
            f"MU: {self.mu}",
            "",
            "Instrumentation:",
            f"Chamber chamber: {self.chamber}",
            f"Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}",
            f"Electrometer: {self.electrometer}",
            f"Pelec: {self.p_elec:2.2f}",
            "",
            "Beam Quality:",
            f"I50 (cm): {self.i_50:2.2f}",
            f"R50 (cm): {self.r_50:2.2f}",
            f"Dref (cm): {self.dref:2.2f}",
            f"Kecal: {self.k_ecal:2.3f}",
            f"kQ: {self.kq:2.3f}",
            "",
            "Chamber Corrections/Measurements:",
            f"Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}",
            f"Pressure (kPa): {self.press:2.1f}",
            f"Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}",
            f"Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}",
            f"Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}",
            f"Ptp: {self.p_tp:2.3f}",
            f"Pion: {self.p_ion:2.3f}",
            f"Ppol: {self.p_pol:2.3f}",
            f"Mraw @ Dref + 0.5rcav (nC): {self.m_gradient}",
            "",
            "Dose Determination:",
            f"Fully corrected M (nC): {self.m_corrected:2.3f}",
            f"Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}",
            f"Dose/MU @ Dref depth (cGy): {self.dose_mu_dref:2.3f}",
            f"Clinical PDD (%): {self.clinical_pdd:2.2f}",
            f"Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}",
            "",
            f"Output Adjustment?: {was_adjusted}",
        ]
        if was_adjusted == "Yes":
            text.append(
                f"Adjusted Mraw @ reference voltage (nC): {self.m_reference_adjusted}"
            )
            text.append(
                f"Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ dref depth (cGy): {self.dose_mu_dref_adjusted:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}"
            )
        canvas.add_text(text=text, location=(2, 25.5), font_size=11)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            webbrowser.open(filename)


class TG51ElectronModern(TG51Base):
    """Class for calculating absolute dose to water using a cylindrical chamber in an electron beam.

    .. warning::
        This class uses the values of Muir & Rogers. These values are likely to be included in the new TG-51
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
    chamber : str or IonChamber
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

    clinical_pdd: float
    i_50: float
    cone: str

    def __init__(
        self,
        *,
        institution: str = "",
        physicist: str = "",
        unit: str = "",
        measurement_date: str = "",
        energy: int,
        temp: float,
        press: float,
        chamber: str | IonChamber,
        n_dw: float,
        electrometer: str = "",
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
        m_reference_adjusted=None,
    ):
        super().__init__(
            temp=temp,
            press=press,
            chamber=(
                IonChambers.from_name(chamber) if isinstance(chamber, str) else chamber
            ),
            n_dw=n_dw,
            p_elec=p_elec,
            voltage_reference=voltage_reference,
            voltage_reduced=voltage_reduced,
            m_reference=m_reference,
            m_opposite=m_opposite,
            m_reduced=m_reduced,
            clinical_pdd=clinical_pdd,
            mu=mu,
            i_50=i_50,
            tissue_correction=tissue_correction,
            institution=institution,
            physicist=physicist,
            unit=unit,
            measurement_date=measurement_date,
            electrometer=electrometer,
            m_reference_adjusted=m_reference_adjusted,
            cone=cone,
            energy=energy,
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
        """The kQ value using the updated Muir & Rogers values from their 2014 paper, equation 11, or classically
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
        return (
            self.tissue_correction
            * self.m_corrected_adjustment
            * self.kq
            * self.n_dw
            / self.mu
        )

    @property
    def dose_mu_dmax_adjusted(self) -> float:
        """cGy/MU at the depth of dmax."""
        return self.dose_mu_dref_adjusted / (self.clinical_pdd / 100)

    def publish_pdf(
        self,
        filename: str,
        notes: list | None = None,
        open_file: bool = False,
        metadata: dict | None = None,
    ):
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
        was_adjusted = "Yes" if self.output_was_adjusted else "No"
        title = ["TG-51 Electron Report (Modern)", f"{self.unit} - {self.energy} MeV"]

        canvas = PylinacCanvas(filename, page_title=title, metadata=metadata)
        text = [
            "Site Data:",
            f"Institution: {self.institution}",
            f"Performed by: {self.physicist}",
            f"Measurement Date: {self.measurement_date}",
            f"Date of Report: {datetime.now().strftime('%A, %B %d, %Y')}",
            f"Unit: {self.unit}",
            f"Energy: {self.energy} MeV",
            f"Cone: {self.cone}",
            f"MU: {self.mu}",
            "",
            "Instrumentation:",
            f"Chamber: {self.chamber}",
            f"Chamber Calibration Factor Ndw (cGy/nC): {self.n_dw:2.3f}",
            f"Electrometer: {self.electrometer}",
            f"Pelec: {self.p_elec:2.2f}",
            "",
            "Beam Quality:",
            f"I50 (cm): {self.i_50:2.2f}",
            f"R50 (cm): {self.r_50:2.2f}",
            f"Dref (cm): {self.dref:2.2f}",
            f"Calculated kQ: {self.kq:2.3f}",
            "",
            "Chamber Corrections/Measurements:",
            f"Temperature (\N{DEGREE SIGN}C): {self.temp:2.1f}",
            f"Pressure (kPa): {self.press:2.1f}",
            f"Mraw @ ({self.voltage_reference}V, Reference) (nC): {self.m_reference}",
            f"Mraw @ ({self.voltage_reduced}V, Reduced) (nC): {self.m_reduced}",
            f"Mraw @ ({-self.voltage_reference}V, Opposite) (nC): {self.m_opposite}",
            f"Ptp: {self.p_tp:2.3f}",
            f"Pion: {self.p_ion:2.3f}",
            f"Ppol: {self.p_pol:2.3f}",
            "",
            "Dose Determination:",
            f"Fully corrected M (nC): {self.m_corrected:2.3f}",
            f"Tissue correction (e.g. muscle): {self.tissue_correction:2.3f}",
            f"Dose/MU @ Dref depth (cGy): {self.dose_mu_dref:2.3f}",
            f"Clinical PDD (%): {self.clinical_pdd:2.2f}",
            f"Dose/MU @ dmax (cGy): {self.dose_mu_dmax:2.3f}",
            "",
            f"Output Adjustment?: {was_adjusted}",
        ]
        if was_adjusted == "Yes":
            text.append(
                f"Adjusted corrected M @ reference voltage (nC): {self.m_corrected_adjustment}"
            )
            text.append(
                f"Adjusted fully corrected M (nC): {self.m_corrected_adjustment:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ dref depth (cGy): {self.dose_mu_dref_adjusted:2.3f}"
            )
            text.append(
                f"Adjusted Dose/MU @ dmax (cGy): {self.dose_mu_dmax_adjusted:2.3f}"
            )
        canvas.add_text(text=text, location=(2, 25.5), font_size=11)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(12, 6.5), font_size=14)
            canvas.add_text(text=notes, location=(12, 6))

        canvas.finish()

        if open_file:
            webbrowser.open(filename)
