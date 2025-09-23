from __future__ import annotations

import datetime
import math
from abc import ABC, abstractmethod
from collections.abc import Iterable
from copy import deepcopy
from enum import Enum
from pathlib import Path
from typing import Literal

import numpy as np
import pydicom
from matplotlib.figure import Figure
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from ..core import scale
from ..core.image_generator.layers import ArrayLayer
from ..core.image_generator.simulators import Simulator
from ..core.scale import wrap360
from .fluence import generate_fluences, plot_fluences
from .mlc import MLCShaper


class GantryDirection(Enum):
    CLOCKWISE = "CW"
    COUNTER_CLOCKWISE = "CC"
    NONE = "NONE"


class GantrySpeedTransition(Enum):
    LEADING = "leading"
    TRAILING = "trailing"


class FluenceMode(Enum):
    STANDARD = "STANDARD"
    FFF = "FFF"
    SRS = "SRS"


class Stack(Enum):
    DISTAL = "distal"
    PROXIMAL = "proximal"
    BOTH = "both"


MLC_MILLENNIUM_BOUNDARIES = (
    list(np.arange(-200, -100 + 1, 10))
    + list(np.arange(-95, 95 + 1, 5))
    + list(np.arange(100, 200 + 1, 10))
)
MLC_120HDMIL_BOUNDARIES = (
    list(np.arange(-110, -40 + 1, 5))
    + list(np.arange(-37.5, 37.5 + 1, 2.5))
    + list(np.arange(40, 110 + 1, 10))
)
MLC_DISTAL_BOUNDARIES = list(np.arange(-140, 140 + 1, 10))
MLC_PROXIMAL_BOUNDARIES = list(np.arange(-145, 145 + 1, 10))


class _Beam(ABC):
    """Represents a DICOM beam dataset. Has methods for creating the dataset and adding control points.
    Generally not created on its own but rather under the hood as part of a PlanGenerator object.

    It contains enough independent logic steps that it's worth separating out from the PlanGenerator class.
    """

    ROUNDING_DECIMALS = 6

    meterset: float

    def __init__(
        self,
        beam_limiting_device_sequence: Sequence,
        beam_name: str,
        energy: float,
        fluence_mode: FluenceMode,
        dose_rate: int,
        metersets: list[float],
        gantry_angles: float | list[float],
        coll_angle: float,
        beam_limiting_device_positions: dict[str, list],
        couch_vrt: float,
        couch_lat: float,
        couch_lng: float,
        couch_rot: float,
    ):
        """
        Parameters
        ----------
        beam_limiting_device_sequence : Sequence
            The beam_limiting_device_sequence as defined in the template plan.
        beam_name : str
            The name of the beam. Must be less than 16 characters.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        metersets : list[float]
            The meter sets for each control point.
        gantry_angles : Union[float, list[float]]
            The gantry angle(s) of the beam. If a single number, it's assumed to be a static beam. If multiple numbers, it's assumed to be a dynamic beam.
        coll_angle : float
            The collimator angle.
        beam_limiting_device_positions : dict[str, list]
            The positions of the beam_limiting_device_positions for each control point,
            where key is the type of beam limiting device (e.g. "MLCX") and the value contains the positions.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation.
        """

        number_of_control_points = len(metersets)

        # The Meterset at a given Control Point is equal to Beam Meterset (300A,0086)
        # specified in the Referenced Beam Sequence (300C,0004) of the RT Fraction Scheme Module,
        # multiplied by the Cumulative Meterset Weight (300A,0134) for the Control Point,
        # divided by the Final Cumulative Meterset Weight (300A,010E)
        # https://dicom.innolitics.com/ciods/rt-plan/rt-beams/300a00b0/300a0111/300a0134
        metersets_weights = np.array(metersets) / metersets[-1]
        self.meterset = np.round(metersets[-1], self.ROUNDING_DECIMALS)

        if len(beam_name) > 16:
            raise ValueError("Beam name must be less than or equal to 16 characters")

        if not isinstance(gantry_angles, Iterable):
            # if it's just a single number (like for a static beam) set it to an array of that value
            gantry_angles = [gantry_angles] * number_of_control_points

        # Round all possible dynamic elements  to avoid floating point comparisons.
        # E.g. to evaluate is an axis is static, all elements should be equal to the first
        # Note: using np.isclose does not solve the problem since the tolerance should be the same
        # as Eclipse/Machine, and we don't know which tolerance they use.
        # Here we assume that their tolerance is tighter than ROUNDING_DECIMALS
        metersets_weights = np.round(metersets_weights, self.ROUNDING_DECIMALS)
        gantry_angles = np.round(gantry_angles, self.ROUNDING_DECIMALS)
        bld_positions = {
            k: np.round(v, self.ROUNDING_DECIMALS)
            for k, v in beam_limiting_device_positions.items()
        }

        # Infer gantry rotation from the gantry angles
        # It assumes the gantry cannot rotate over 180, so there is only one possible direction to go from A to B.
        ga_wrap180 = scale.wrap180(np.array(gantry_angles))
        # This dictionary is used for mapping the sign of the difference with the GantryDirection enum.
        gantry_direction_map = {
            0: GantryDirection.NONE,
            1: GantryDirection.CLOCKWISE,
            -1: GantryDirection.COUNTER_CLOCKWISE,
        }
        gantry_direction = [
            gantry_direction_map[s] for s in np.sign(np.diff(ga_wrap180))
        ]
        # The last GantryRotationDirection should always be 'NONE'
        gantry_direction += [GantryDirection.NONE]

        # Infer if a beam is static or dynamic from the control points
        gantry_is_static = len(set(gantry_direction)) == 1
        dict_bld_is_static = {
            k: np.all(pos == pos[0]) for k, pos in bld_positions.items()
        }
        blds_are_static = np.all(list(dict_bld_is_static.values()))
        beam_is_static = gantry_is_static and blds_are_static
        beam_type = "STATIC" if beam_is_static else "DYNAMIC"

        # Create dataset with basic beam info
        self.ds = self._create_basic_beam_info(
            beam_name,
            beam_type,
            fluence_mode,
            beam_limiting_device_sequence=beam_limiting_device_sequence,
            number_of_control_points=number_of_control_points,
        )

        # Add initial control point
        cp0 = Dataset()
        cp0.ControlPointIndex = 0
        cp0.NominalBeamEnergy = energy
        cp0.DoseRateSet = dose_rate
        beam_limiting_device_position_sequence = Sequence()
        for key, values in bld_positions.items():
            beam_limiting_device_position = Dataset()
            beam_limiting_device_position.RTBeamLimitingDeviceType = key
            beam_limiting_device_position.LeafJawPositions = list(values[0])
            beam_limiting_device_position_sequence.append(beam_limiting_device_position)
        cp0.BeamLimitingDevicePositionSequence = beam_limiting_device_position_sequence
        cp0.GantryAngle = gantry_angles[0]
        cp0.GantryRotationDirection = gantry_direction[0].value
        cp0.BeamLimitingDeviceAngle = coll_angle
        cp0.BeamLimitingDeviceRotationDirection = "NONE"
        cp0.PatientSupportAngle = couch_rot
        cp0.PatientSupportRotationDirection = "NONE"
        cp0.TableTopEccentricAngle = 0.0
        cp0.TableTopEccentricRotationDirection = "NONE"
        cp0.TableTopVerticalPosition = couch_vrt
        cp0.TableTopLongitudinalPosition = couch_lng
        cp0.TableTopLateralPosition = couch_lat
        cp0.IsocenterPosition = None
        cp0.CumulativeMetersetWeight = 0.0
        self.ds.ControlPointSequence.append(cp0)

        # Add rest of the control points
        for cp_idx in range(1, number_of_control_points):
            cp = Dataset()
            cp.ControlPointIndex = cp_idx
            cp.CumulativeMetersetWeight = metersets_weights[cp_idx]

            if not gantry_is_static:
                cp.GantryAngle = gantry_angles[cp_idx]
                cp.GantryRotationDirection = gantry_direction[cp_idx].value

            bld_position_sequence = Sequence()
            for bld, positions in bld_positions.items():
                if not dict_bld_is_static[bld]:
                    bld_position = Dataset()
                    bld_position.RTBeamLimitingDeviceType = bld
                    bld_position.LeafJawPositions = list(positions[cp_idx])
                    bld_position_sequence.append(bld_position)
            if len(bld_position_sequence) > 0:
                cp.BeamLimitingDevicePositionSequence = bld_position_sequence

            self.ds.ControlPointSequence.append(cp)

    def as_dicom(self) -> Dataset:
        """Return the beam as a DICOM dataset that represents a BeamSequence item."""
        return self.ds

    @staticmethod
    def _create_basic_beam_info(
        beam_name: str,
        beam_type: str,
        fluence_mode: FluenceMode,
        beam_limiting_device_sequence: Sequence,
        number_of_control_points: int,
    ) -> Dataset:
        beam = Dataset()
        beam.Manufacturer = "Radformation"
        beam.ManufacturerModelName = "RadMachine"
        beam.PrimaryDosimeterUnit = "MU"
        beam.SourceAxisDistance = 1000.0

        # Primary Fluence Mode Sequence
        primary_fluence_mode1 = Dataset()
        if fluence_mode == FluenceMode.STANDARD:
            primary_fluence_mode1.FluenceMode = "STANDARD"
        elif fluence_mode == FluenceMode.FFF:
            primary_fluence_mode1.FluenceMode = "NON_STANDARD"
            primary_fluence_mode1.FluenceModeID = "FFF"
        elif fluence_mode == FluenceMode.SRS:
            primary_fluence_mode1.FluenceMode = "NON_STANDARD"
            primary_fluence_mode1.FluenceModeID = "SRS"
        beam.PrimaryFluenceModeSequence = Sequence((primary_fluence_mode1,))

        # Beam Limiting Device Sequence
        beam.BeamLimitingDeviceSequence = beam_limiting_device_sequence

        # beam numbers start at 0 and increment from there.
        beam.BeamName = beam_name
        beam.BeamType = beam_type
        beam.RadiationType = "PHOTON"
        beam.TreatmentDeliveryType = "TREATMENT"
        beam.NumberOfWedges = 0
        beam.NumberOfCompensators = 0
        beam.NumberOfBoli = 0
        beam.NumberOfBlocks = 0
        beam.FinalCumulativeMetersetWeight = 1.0
        beam.NumberOfControlPoints = number_of_control_points

        # Control Point Sequence
        beam.ControlPointSequence = Sequence()
        return beam


class TrueBeamBeam(_Beam):
    """Represents a DICOM beam dataset for a TrueBeam. Has methods for creating the dataset and adding control points.
    Generally not created on its own but rather under the hood as part of a PlanGenerator object.

    It contains enough independent logic steps that it's worth separating out from the PlanGenerator class.
    """

    def __init__(
        self,
        is_mlc_hd: bool,
        beam_name: str,
        energy: float,
        fluence_mode: FluenceMode,
        dose_rate: int,
        metersets: list[float],
        gantry_angles: float | list[float],
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        mlc_positions: list[list[float]],
        coll_angle: float,
        couch_vrt: float,
        couch_lat: float,
        couch_lng: float,
        couch_rot: float,
    ):
        """
        Parameters
        ----------
        is_mlc_hd : bool
            Whether the MLC type is HD or Millennium
        beam_name : str
            The name of the beam. Must be less than 16 characters.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        metersets : list[float]
            The meter sets for each control point. The length must match the number of control points in mlc_positions.
        gantry_angles : Union[float, list[float]]
            The gantry angle(s) of the beam. If a single number, it's assumed to be a static beam. If multiple numbers, it's assumed to be a dynamic beam.
        x1 : float
            The left jaw position.
        x2 : float
            The right jaw position.
        y1 : float
            The bottom jaw position.
        y2 : float
            The top jaw position.
        mlc_positions : list[list[float]]
            The MLC positions for each control point. This is the x-position of each leaf for each control point.
        coll_angle : float
            The collimator angle.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation.
        """
        jaw_x = Dataset()
        jaw_x.RTBeamLimitingDeviceType = "X"
        jaw_x.NumberOfLeafJawPairs = 1
        jaw_y = Dataset()
        jaw_y.RTBeamLimitingDeviceType = "Y"
        jaw_y.NumberOfLeafJawPairs = 1
        jaw_asymx = Dataset()
        jaw_asymx.RTBeamLimitingDeviceType = "ASYMX"
        jaw_asymx.NumberOfLeafJawPairs = 1
        jaw_asymy = Dataset()
        jaw_asymy.RTBeamLimitingDeviceType = "ASYMX"
        jaw_asymy.NumberOfLeafJawPairs = 1
        mlc = Dataset()
        mlc.RTBeamLimitingDeviceType = "MLCX"
        mlc.NumberOfLeafJawPairs = 60
        if is_mlc_hd:
            mlc.LeafPositionBoundaries = MLC_120HDMIL_BOUNDARIES
        else:
            mlc.LeafPositionBoundaries = MLC_MILLENNIUM_BOUNDARIES
        bld_sequence = Sequence((jaw_x, jaw_y, jaw_asymx, jaw_asymy, mlc))

        beam_limiting_device_positions = {
            "ASYMX": [[x1, x2]],
            "ASYMY": [[y1, y2]],
            "MLCX": mlc_positions,
        }

        super().__init__(
            beam_limiting_device_sequence=bld_sequence,
            beam_name=beam_name,
            energy=energy,
            fluence_mode=fluence_mode,
            dose_rate=dose_rate,
            metersets=metersets,
            gantry_angles=gantry_angles,
            beam_limiting_device_positions=beam_limiting_device_positions,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
        )


class HalcyonBeam(_Beam):
    """A beam representing a Halcyon. Halcyons have dual MLC stacks, no X-jaws, no couch rotation, etc."""

    def __init__(
        self,
        beam_name: str,
        metersets: list[float],
        gantry_angles: float | list[float],
        distal_mlc_positions: list[list[float]],
        proximal_mlc_positions: list[list[float]],
        coll_angle: float,
        couch_vrt: float,
        couch_lat: float,
        couch_lng: float,
    ):
        """
        Parameters
        ----------
        beam_name : str
            The name of the beam. Must be less than 16 characters.
        metersets : list[float]
            The meter sets for each control point. The length must match the number of control points in mlc_positions.
        gantry_angles : Union[float, list[float]]
            The gantry angle(s) of the beam. If a single number, it's assumed to be a static beam. If multiple numbers, it's assumed to be a dynamic beam.
        distal_mlc_positions : list[list[float]]
            The distal MLC positions for each control point. This is the x-position of each leaf for each control point.
        proximal_mlc_positions : list[list[float]]
            The proximal MLC positions for each control point. This is the x-position of each leaf for each control point.
        coll_angle : float
            The collimator angle.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        """
        jaw_x = Dataset()
        jaw_x.RTBeamLimitingDeviceType = "X"
        jaw_x.NumberOfLeafJawPairs = 1
        jaw_y = Dataset()
        jaw_y.RTBeamLimitingDeviceType = "Y"
        jaw_y.NumberOfLeafJawPairs = 1
        mlc_x1 = Dataset()
        mlc_x1.RTBeamLimitingDeviceType = "MLCX1"
        mlc_x1.NumberOfLeafJawPairs = 28
        mlc_x1.LeafPositionBoundaries = MLC_DISTAL_BOUNDARIES
        mlc_x2 = Dataset()
        mlc_x2.RTBeamLimitingDeviceType = "MLCX2"
        mlc_x2.NumberOfLeafJawPairs = 29
        mlc_x2.LeafPositionBoundaries = MLC_PROXIMAL_BOUNDARIES
        bld_sequence = Sequence((jaw_x, jaw_y, mlc_x1, mlc_x2))

        beam_limiting_device_positions = {
            "X": [[-140, 140]],
            "Y": [[-140, 140]],
            "MLCX1": distal_mlc_positions,
            "MLCX2": proximal_mlc_positions,
        }

        super().__init__(
            beam_limiting_device_sequence=bld_sequence,
            beam_name=beam_name,
            energy=6,
            fluence_mode=FluenceMode.FFF,
            dose_rate=600,
            metersets=metersets,
            gantry_angles=gantry_angles,
            beam_limiting_device_positions=beam_limiting_device_positions,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=0,
        )


class PlanGenerator(ABC):
    """A tool for generating new QA RTPlan files based on an initial, somewhat empty RTPlan file.

    Attributes
    ----------
    machine_name : str
        The name of the machine
    """

    machine_name: str

    def __init__(
        self,
        ds: Dataset,
        plan_label: str,
        plan_name: str,
        patient_name: str | None,
        patient_id: str | None,
        max_mlc_position: float,
        max_mlc_speed: float,
        max_gantry_speed: float,
        max_overtravel_mm: float,
    ):
        """A tool for generating new QA RTPlan files based on an initial, somewhat empty RTPlan file.

        Parameters
        ----------
        ds : Dataset
              The RTPLAN dataset to base the new plan off of. The plan must already have MLC positions.
        plan_label : str
            The label of the new plan.
        plan_name : str
            The name of the new plan.
        patient_name : str, optional
            The name of the patient. If not provided, it will be taken from the RTPLAN file.
        patient_id : str, optional
            The ID of the patient. If not provided, it will be taken from the RTPLAN file.
        max_mlc_position : float
            The max mlc position in mm
        max_mlc_speed : float
            The maximum speed of the MLC leaves in mm/s
        max_gantry_speed : float
            The maximum speed of the gantry in degrees/s.
        max_overtravel_mm : float
            The maximum distance the MLC leaves can overtravel from each other as well as the jaw size (for tail exposure protection).
        """
        if ds.Modality != "RTPLAN":
            raise ValueError("File is not an RTPLAN file")
        self.max_overtravel_mm = max_overtravel_mm
        self.max_mlc_position = max_mlc_position
        self.max_mlc_speed = max_mlc_speed
        self.max_gantry_speed = max_gantry_speed
        patient_name = patient_name or getattr(ds, "PatientName", None)
        if not patient_name:
            raise ValueError(
                "RTPLAN file must have PatientName or pass it via `patient_name`"
            )
        patient_id = patient_id or getattr(ds, "PatientID", None)
        if not patient_id:
            raise ValueError(
                "RTPLAN file must have PatientID or pass it via `patient_id`"
            )
        if not hasattr(ds, "ToleranceTableSequence"):
            raise ValueError("RTPLAN file must have ToleranceTableSequence")
        if not hasattr(ds, "BeamSequence"):
            raise ValueError(
                "RTPLAN file must have at least one beam in the beam sequence"
            )
        has_mlc_data: bool = any(
            "MLC" in bld.RTBeamLimitingDeviceType
            for bs in ds.BeamSequence
            for bld in bs.BeamLimitingDeviceSequence
        )
        if not has_mlc_data:
            raise ValueError("RTPLAN file must have MLC data")

        ######  Create a copy of the template plan
        # A shallow copy won’t work because beam data is cleared.
        # The inherited classes require access to the original beam state to determine leaf boundaries
        self.ds = deepcopy(ds)

        ######  Clear/initialize the metadata for the new plan
        self.ds.PatientName = patient_name
        self.ds.PatientID = patient_id
        self.ds.RTPlanLabel = plan_label
        self.ds.RTPlanName = plan_name
        date = datetime.datetime.now().strftime("%Y%m%d")
        time = datetime.datetime.now().strftime("%H%M%S")

        self.ds.InstanceCreationDate = date
        self.ds.InstanceCreationTime = time
        self.ds.SOPInstanceUID = generate_uid()

        # Patient Setup Sequence
        patient_setup = Dataset()
        patient_setup.PatientPosition = "HFS"
        patient_setup.PatientSetupNumber = 0
        self.ds.PatientSetupSequence = Sequence((patient_setup,))

        # Dose Reference Sequence
        dose_ref1 = Dataset()
        dose_ref1.DoseReferenceNumber = 1
        dose_ref1.DoseReferenceUID = generate_uid()
        dose_ref1.DoseReferenceStructureType = "SITE"
        dose_ref1.DoseReferenceDescription = "PTV"
        dose_ref1.DoseReferenceType = "TARGET"
        dose_ref1.DeliveryMaximumDose = 20.0
        dose_ref1.TargetPrescriptionDose = 40.0
        dose_ref1.TargetMaximumDose = 20.0
        self.ds.DoseReferenceSequence = Sequence((dose_ref1,))

        # Fraction Group Sequence
        frxn_gp_sequence = Sequence()
        self.ds.FractionGroupSequence = frxn_gp_sequence
        # Fraction Group Sequence: Fraction Group 1
        frxn_gp1 = Dataset()
        frxn_gp1.FractionGroupNumber = 1
        frxn_gp1.NumberOfFractionsPlanned = 1
        frxn_gp1.NumberOfBeams = 0
        frxn_gp1.NumberOfBrachyApplicationSetups = 0
        frxn_gp1.ReferencedBeamSequence = Sequence()
        self.ds.FractionGroupSequence = Sequence((frxn_gp1,))

        # Clear beam sequence
        # This will be filled with the custom beams
        self.ds.BeamSequence = Sequence()

        # Machine name
        self.machine_name = ds.BeamSequence[0].TreatmentMachineName

        # Validate machine type
        self._validate_machine_type(ds.BeamSequence)

    @classmethod
    def from_rt_plan_file(cls, rt_plan_file: str | Path, **kwargs) -> PlanGenerator:
        """Load an existing RTPLAN file and create a new plan based on it.

        Parameters
        ----------
        rt_plan_file : str
            The path to the RTPLAN file.
        kwargs
            See the PlanGenerator constructor for details.
        """
        ds = pydicom.dcmread(rt_plan_file)
        return cls(ds, **kwargs)

    @abstractmethod
    def _validate_machine_type(self, beam_sequence: Sequence):
        pass

    def add_beam(self, beam: HalcyonBeam | TrueBeamBeam):
        """Add a beam to the plan using the Beam object. Although public,
        this is a low-level method that is used by the higher-level methods like add_open_field_beam.
        This handles the associated metadata like the referenced beam sequence and fraction group sequence.
        """
        beam_dataset = beam.as_dicom()

        # Update the beam
        beam_dataset.BeamNumber = len(self.ds.BeamSequence) + 1
        beam_dataset.TreatmentMachineName = self.machine_name
        patient_setup_nr = self.ds.PatientSetupSequence[0].PatientSetupNumber
        beam_dataset.ReferencedPatientSetupNumber = patient_setup_nr
        tolerance_table_nr = self.ds.ToleranceTableSequence[0].ToleranceTableNumber
        beam_dataset.ReferencedToleranceTableNumber = tolerance_table_nr
        self.ds.BeamSequence.append(beam_dataset)

        # increment number of beams
        fr = self.ds.FractionGroupSequence[0]
        fr.NumberOfBeams += 1

        # Update plan references
        referenced_beam = Dataset()
        referenced_beam.BeamDose = 1.0
        referenced_beam.BeamMeterset = beam.meterset
        referenced_beam.ReferencedBeamNumber = beam_dataset.BeamNumber
        dose_reference_uid = self.ds.DoseReferenceSequence[0].DoseReferenceUID
        referenced_beam.ReferencedDoseReferenceUID = dose_reference_uid
        self.ds.FractionGroupSequence[0].ReferencedBeamSequence.append(referenced_beam)

    def to_file(self, filename: str | Path) -> None:
        """Write the DICOM dataset to file"""
        self.ds.save_as(filename, write_like_original=False)

    def as_dicom(self) -> Dataset:
        """Return the new DICOM dataset."""
        return self.ds

    def plot_fluences(
        self,
        width_mm: float = 400,
        resolution_mm: float = 0.5,
        dtype: np.dtype = np.uint16,
    ) -> list[Figure]:
        """Plot the fluences of the beams generated

        See Also
        --------
        :func:`~pydicom_planar.PlanarImage.plot_fluences`
        """
        return plot_fluences(self.as_dicom(), width_mm, resolution_mm, dtype, show=True)

    def to_dicom_images(
        self, simulator: type[Simulator], invert: bool = True
    ) -> list[Dataset]:
        """Generate simulated DICOM images of the plan. This provides a way to
        generate an end-to-end simulation of the plan. The images will always be
        at 1000mm SID.

        Parameters
        ----------
        simulator : Simulator
            The simulator to use to generate the images. This provides the
            size of the image and the pixel size
        invert: bool
            Invert the fluence. Setting to True simulates EPID-style images where
            dose->lower pixel value.
        """
        image_ds = []
        fluences = generate_fluences(
            rt_plan=self.as_dicom(),
            width_mm=simulator.shape[1] * simulator.pixel_size,
            resolution_mm=simulator.pixel_size,
        )
        for beam, fluence in zip(self.ds.BeamSequence, fluences):
            beam_info = beam.ControlPointSequence[0]
            sim = simulator(sid=1000)
            sim.add_layer(ArrayLayer(fluence))
            ds = sim.as_dicom(
                gantry_angle=beam_info.GantryAngle,
                coll_angle=beam_info.BeamLimitingDeviceAngle,
                table_angle=beam_info.PatientSupportAngle,
                invert_array=invert,
            )
            image_ds.append(ds)
        return image_ds


class TrueBeamPlanGenerator(PlanGenerator):
    # Private fields
    _is_mlc_hd: bool  # Whether the MLC type id HD or Millennium.
    _leaf_boundaries: list[float]  # The leaf boundaries of the MLC.

    def __init__(
        self,
        ds: Dataset,
        plan_label: str,
        plan_name: str,
        patient_name: str | None = None,
        patient_id: str | None = None,
        max_mlc_position: float = 200,
        max_mlc_speed: float = 25,
        max_gantry_speed: float = 4.8,
        max_overtravel_mm: float = 140,
    ):
        super().__init__(
            ds,
            plan_label,
            plan_name,
            patient_name,
            patient_id,
            max_mlc_position,
            max_mlc_speed,
            max_gantry_speed,
            max_overtravel_mm,
        )

        # Fill Fields
        self._is_mlc_hd = any(
            bld.LeafPositionBoundaries[0] == -110
            for bs in ds.BeamSequence
            for bld in bs.BeamLimitingDeviceSequence
            if bld.RTBeamLimitingDeviceType == "MLCX"
        )
        self._leaf_boundaries = (
            MLC_120HDMIL_BOUNDARIES if self._is_mlc_hd else MLC_MILLENNIUM_BOUNDARIES
        )

    def _validate_machine_type(self, beam_sequence: Sequence):
        has_valid_mlc_data: bool = any(
            bld.RTBeamLimitingDeviceType == "MLCX"
            for bs in beam_sequence
            for bld in bs.BeamLimitingDeviceSequence
        )
        if not has_valid_mlc_data:
            raise ValueError(
                "The machine on the template plan does not seem to be a TrueBeam machine."
            )

    def _create_mlc(
        self, sacrifice_gap_mm: float = None, sacrifice_max_move_mm: float = None
    ) -> MLCShaper:
        """Utility to create MLC shaper instances."""
        return MLCShaper(
            leaf_y_positions=self._leaf_boundaries,
            max_mlc_position=self.max_mlc_position,
            sacrifice_gap_mm=sacrifice_gap_mm,
            sacrifice_max_move_mm=sacrifice_max_move_mm,
            max_overtravel_mm=self.max_overtravel_mm,
        )

    def add_picketfence_beam(
        self,
        strip_width_mm: float = 3,
        strip_positions_mm: tuple[float | int, ...] = (-45, -30, -15, 0, 15, 30, 45),
        y1: float = -100,
        y2: float = 100,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        energy: float = 6,
        gantry_angle: float = 0,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lng: float = 1000,
        couch_lat: float = 0,
        couch_rot: float = 0,
        mu: int = 200,
        jaw_padding_mm: float = 10,
        beam_name: str = "PF",
        max_sacrificial_move_mm: float = 50,
    ):
        """Add a picket fence beam to the plan.

        Parameters
        ----------
        strip_width_mm : float
            The width of the strips in mm.
        strip_positions_mm : tuple
            The positions of the strips in mm relative to the center of the image.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        energy : float
            The energy of the beam.
        gantry_angle : float
            The gantry angle of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lng : float
            The couch longitudinal position.
        couch_lat : float
            The couch lateral position.
        couch_rot : float
            The couch rotation.
        mu : int
            The monitor units of the beam.
        jaw_padding_mm : float
            The padding to add to the X jaws.
        beam_name : str
            The name of the beam.
        max_sacrificial_move_mm : float
            The maximum distance the sacrificial leaves can move in a given control point.
            Smaller values generate more control points and more back-and-forth movement.
            Too large of values may cause deliverability issues.
        """
        # check MLC overtravel; machine may prevent delivery if exposing leaf tail
        x1 = min(strip_positions_mm) - jaw_padding_mm
        x2 = max(strip_positions_mm) + jaw_padding_mm
        max_dist_to_jaw = max(
            max(abs(pos - x1), abs(pos + x2)) for pos in strip_positions_mm
        )
        if max_dist_to_jaw > self.max_overtravel_mm:
            raise ValueError(
                "Picket fence beam exceeds MLC overtravel limits. Lower padding, the number of pickets, or the picket spacing."
            )
        mlc = self._create_mlc(sacrifice_max_move_mm=max_sacrificial_move_mm)
        # create initial starting point; start under the jaws
        mlc.add_strip(
            position_mm=strip_positions_mm[0] - 2,
            strip_width_mm=strip_width_mm,
            meterset_at_target=0,
        )

        for strip in strip_positions_mm:
            # starting control point
            mlc.add_strip(
                position_mm=strip,
                strip_width_mm=strip_width_mm,
                meterset_at_target=1 / len(strip_positions_mm),
            )
        beam = TrueBeamBeam(
            beam_name=beam_name,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            fluence_mode=fluence_mode,
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)

    def add_mlc_transmission(
        self,
        bank: Literal["A", "B"],
        mu: int = 50,
        overreach: float = 10,
        beam_name: str = "MLC Tx",
        energy: int = 6,
        dose_rate: int = 600,
        x1: float = -50,
        x2: float = 50,
        y1: float = -100,
        y2: float = 100,
        gantry_angle: float = 0,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lat: float = 0,
        couch_lng: float = 1000,
        couch_rot: float = 0,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
    ):
        """Add a single-image MLC transmission beam to the plan.
        The beam is delivered with the MLCs closed and moved to one side underneath the jaws.

        Parameters
        ----------
        bank : str
            The MLC bank to move. Either "A" or "B".
        mu : int
            The monitor units to deliver.
        overreach : float
            The amount to tuck the MLCs under the jaws in mm.
        beam_name : str
            The name of the beam.
        energy : int
            The energy of the beam.
        dose_rate : int
            The dose rate of the beam.
        x1 : float
            The left jaw position. Usually negative. More negative is left.
        x2 : float
            The right jaw position. Usually positive. More positive is right.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        gantry_angle : float
            The gantry angle of the beam in degrees.
        coll_angle : float
            The collimator angle of the beam in degrees.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation in degrees.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        """
        mlc = self._create_mlc()
        if bank == "A":
            mlc_tips = x2 + overreach
        elif bank == "B":
            mlc_tips = x1 - overreach
        else:
            raise ValueError("Bank must be 'A' or 'B'")
        # test for overtravel
        if abs(x2 - x1) + overreach > self.max_overtravel_mm:
            raise OvertravelError(
                "The MLC overtravel is too large for the given jaw positions and overreach. Reduce the x-jaw opening size and/or overreach value."
            )
        mlc.add_strip(
            position_mm=mlc_tips,
            strip_width_mm=1,
            meterset_at_target=1,
        )
        beam = TrueBeamBeam(
            beam_name=f"{beam_name} {bank}",
            energy=energy,
            dose_rate=dose_rate,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            fluence_mode=fluence_mode,
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)

    def add_dose_rate_beams(
        self,
        dose_rates: tuple[int, ...] = (100, 300, 500, 600),
        default_dose_rate: int = 600,
        gantry_angle: float = 0,
        desired_mu: int = 50,
        energy: float = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lat: float = 0,
        couch_lng: float = 1000,
        couch_rot: float = 0,
        jaw_padding_mm: float = 5,
        roi_size_mm: float = 25,
        y1: float = -100,
        y2: float = 100,
        max_sacrificial_move_mm: float = 50,
    ):
        """Create a single-image dose rate test. Multiple ROIs are generated. A reference beam is also
        created where all ROIs are delivered at the default dose rate for comparison.
        The field names are generated automatically based on the min and max dose rates tested.

        Parameters
        ----------
        dose_rates : tuple
            The dose rates to test in MU/min. Each dose rate will have its own ROI.
        default_dose_rate : int
            The default dose rate. Typically, this is the clinical default. The reference beam
            will be delivered at this dose rate for all ROIs.
        gantry_angle : float
            The gantry angle of the beam.
        desired_mu : int
            The desired monitor units to deliver. It can be that based on the dose rates asked for,
            the MU required might be higher than this value.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation.
        jaw_padding_mm : float
            The padding to add to the X jaws. The X-jaws will close around the ROIs plus this padding.
        roi_size_mm : float
            The width of the ROIs in mm.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        max_sacrificial_move_mm : float
            The maximum distance the sacrificial leaves can move in a given control point.
            Smaller values generate more control points and more back-and-forth movement.
            Too large of values may cause deliverability issues.
        """
        if roi_size_mm * len(dose_rates) > self.max_overtravel_mm:
            raise ValueError(
                "The ROI size * number of dose rates must be less than the overall MLC allowable width"
            )
        # calculate MU
        mlc_transition_time = roi_size_mm / self.max_mlc_speed
        min_mu = mlc_transition_time * max(dose_rates) * len(dose_rates) / 60
        mu = max(desired_mu, math.ceil(min_mu))

        # create MLC sacrifices
        times_to_transition = [
            mu * 60 / (dose_rate * len(dose_rates)) for dose_rate in dose_rates
        ]
        sacrificial_movements = [tt * self.max_mlc_speed for tt in times_to_transition]

        mlc = self._create_mlc(sacrifice_max_move_mm=max_sacrificial_move_mm)
        ref_mlc = self._create_mlc()

        roi_centers = np.linspace(
            -roi_size_mm * len(dose_rates) / 2 + roi_size_mm / 2,
            roi_size_mm * len(dose_rates) / 2 - roi_size_mm / 2,
            len(dose_rates),
        )
        # we have a starting and ending strip
        ref_mlc.add_strip(
            position_mm=float(roi_centers[0] - roi_size_mm / 2),
            strip_width_mm=0,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=float(roi_centers[0] - roi_size_mm / 2),
            strip_width_mm=0,
            meterset_at_target=0,
            initial_sacrificial_gap_mm=5,
        )
        for sacrifice_distance, center in zip(sacrificial_movements, roi_centers):
            ref_mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,
                top_position=max(self._leaf_boundaries),
                bottom_position=min(self._leaf_boundaries),
                outer_strip_width=5,
                meterset_at_target=0,
                meterset_transition=0.5 / len(dose_rates),
                sacrificial_distance=0,
            )
            ref_mlc.add_strip(
                position_mm=center + roi_size_mm / 2,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=0.5 / len(dose_rates),
                sacrificial_distance_mm=0,
            )
            mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,  # not used
                top_position=max(self._leaf_boundaries),
                bottom_position=min(self._leaf_boundaries),
                outer_strip_width=5,  # not used
                meterset_at_target=0,
                meterset_transition=0.5 / len(dose_rates),
                sacrificial_distance=sacrifice_distance,
            )
            mlc.add_strip(
                position_mm=center + roi_size_mm / 2,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=0.5 / len(dose_rates),
                sacrificial_distance_mm=sacrifice_distance,
            )
        ref_beam = TrueBeamBeam(
            beam_name="DR Ref",
            energy=energy,
            dose_rate=default_dose_rate,
            x1=float(roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm),
            x2=float(roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm),
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=[mu * m for m in ref_mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(ref_beam)
        beam = TrueBeamBeam(
            beam_name=f"DR{min(dose_rates)}-{max(dose_rates)}",
            energy=energy,
            dose_rate=default_dose_rate,
            x1=float(roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm),
            x2=float(roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm),
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)

    def add_mlc_speed_beams(
        self,
        speeds: tuple[float | int, ...] = (5, 10, 15, 20),
        roi_size_mm: float = 20,
        mu: int = 50,
        default_dose_rate: int = 600,
        gantry_angle: float = 0,
        energy: float = 6,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lat: float = 0,
        couch_lng: float = 1000,
        couch_rot: float = 0,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        jaw_padding_mm: float = 5,
        y1: float = -100,
        y2: float = 100,
        beam_name: str = "MLC Speed",
        max_sacrificial_move_mm: float = 50,
    ):
        """Create a single-image MLC speed test. Multiple speeds are generated. A reference beam is also
        generated. The reference beam is delivered at the maximum MLC speed.

        Parameters
        ----------
        speeds : tuple[float]
            The speeds to test in mm/s. Each speed will have its own ROI.
        roi_size_mm : float
            The width of the ROIs in mm.
        mu : int
            The monitor units to deliver.
        default_dose_rate : int
            The dose rate used for the reference beam.
        gantry_angle : float
            The gantry angle of the beam.
        energy : int
            The energy of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        jaw_padding_mm : float
            The padding to add to the X jaws. The X-jaws will close around the ROIs plus this padding.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        beam_name : str
            The name of the beam. The reference beam will be called "MLC Sp Ref".
        max_sacrificial_move_mm : float
            The maximum distance the sacrificial leaves can move in a given control point.
            Smaller values generate more control points and more back-and-forth movement.
            Too large of values may cause deliverability issues.


        Notes
        -----

        The desired speed can be achieved through the following formula:

           speed = roi_size_mm * max dose rate / MU * 60

        We solve for MU with the desired speed. The 60 is for converting the dose rate as MU/min to MU/sec.
        Thus,

            MU = roi_size_mm * max dose rate / speed * 60

        MUs are calculated automatically based on the speed and the ROI size.

        """
        if max(speeds) > self.max_mlc_speed:
            raise ValueError(
                f"Maximum speed given {max(speeds)} is greater than the maximum MLC speed {self.max_mlc_speed}"
            )
        if min(speeds) <= 0:
            raise ValueError("Speeds must be greater than 0")
        if roi_size_mm * len(speeds) > self.max_overtravel_mm:
            raise ValueError(
                "The ROI size * number of speeds must be less than the overall MLC allowable width"
            )
        # create MLC positions
        times_to_transition = [roi_size_mm / speed for speed in speeds]
        sacrificial_movements = [tt * self.max_mlc_speed for tt in times_to_transition]

        mlc = self._create_mlc(sacrifice_max_move_mm=max_sacrificial_move_mm)
        ref_mlc = self._create_mlc()

        roi_centers = np.linspace(
            -roi_size_mm * len(speeds) / 2 + roi_size_mm / 2,
            roi_size_mm * len(speeds) / 2 - roi_size_mm / 2,
            len(speeds),
        )
        # we have a starting and ending strip
        ref_mlc.add_strip(
            position_mm=float(roi_centers[0] - roi_size_mm / 2),
            strip_width_mm=0,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=float(roi_centers[0] - roi_size_mm / 2),
            strip_width_mm=0,
            meterset_at_target=0,
            initial_sacrificial_gap_mm=5,
        )
        for sacrifice_distance, center in zip(sacrificial_movements, roi_centers):
            ref_mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,  # not relevant
                top_position=max(self._leaf_boundaries),
                bottom_position=min(self._leaf_boundaries),
                outer_strip_width=5,  # not relevant
                meterset_at_target=0,
                meterset_transition=0.5 / len(speeds),
                sacrificial_distance=0,
            )
            ref_mlc.add_strip(
                position_mm=center + roi_size_mm / 2,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=0.5 / len(speeds),
                sacrificial_distance_mm=0,
            )
            mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,  # not used
                top_position=max(self._leaf_boundaries),
                bottom_position=min(self._leaf_boundaries),
                outer_strip_width=5,  # not used
                meterset_at_target=0,
                meterset_transition=0.5 / len(speeds),
                sacrificial_distance=sacrifice_distance,
            )
            mlc.add_strip(
                position_mm=center + roi_size_mm / 2,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=0.5 / len(speeds),
                sacrificial_distance_mm=sacrifice_distance,
            )
        ref_beam = TrueBeamBeam(
            beam_name=f"{beam_name} Ref",
            energy=energy,
            dose_rate=default_dose_rate,
            x1=float(roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm),
            x2=float(roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm),
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=[mu * m for m in ref_mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(ref_beam)
        beam = TrueBeamBeam(
            beam_name=beam_name,
            energy=energy,
            dose_rate=default_dose_rate,
            x1=float(roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm),
            x2=float(roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm),
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)

    def add_winston_lutz_beams(
        self,
        x1: float = -10,
        x2: float = 10,
        y1: float = -10,
        y2: float = 10,
        defined_by_mlcs: bool = True,
        energy: float = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        axes_positions: Iterable[dict] = ({"gantry": 0, "collimator": 0, "couch": 0},),
        couch_vrt: float = 0,
        couch_lng: float = 1000,
        couch_lat: float = 0,
        mu: int = 10,
        padding_mm: float = 5,
    ):
        """Add Winston-Lutz beams to the plan. Will create a beam for each set of axes positions.
        Field names are generated automatically based on the axes positions.

        Parameters
        ----------
        x1 : float
            The left jaw position.
        x2 : float
            The right jaw position.
        y1 : float
            The bottom jaw position.
        y2 : float
            The top jaw position.
        defined_by_mlcs : bool
            Whether the field edges are defined by the MLCs or the jaws.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        axes_positions : Iterable[dict]
            The positions of the axes. Each dict should have keys 'gantry', 'collimator', 'couch', and optionally 'name'.
        couch_vrt : float
            The couch vertical position.
        couch_lng : float
            The couch longitudinal position.
        couch_lat : float
            The couch lateral position.
        mu : int
            The monitor units of the beam.
        padding_mm : float
            The padding to add. If defined by the MLCs, this is the padding of the jaws. If defined by the jaws,
            this is the padding of the MLCs.
        """
        for axes in axes_positions:
            if defined_by_mlcs:
                mlc_padding = 0
                jaw_padding = padding_mm
            else:
                mlc_padding = padding_mm
                jaw_padding = 0
            mlc = self._create_mlc()
            mlc.add_rectangle(
                left_position=x1 - mlc_padding,
                right_position=x2 + mlc_padding,
                top_position=y2 + mlc_padding,
                bottom_position=y1 - mlc_padding,
                outer_strip_width=5,
                meterset_at_target=1.0,
                x_outfield_position=x1 - mlc_padding - jaw_padding - 20,
            )
            beam_name = (
                axes.get("name")
                or f"G{axes['gantry']:g}C{axes['collimator']:g}P{axes['couch']:g}"
            )
            beam = TrueBeamBeam(
                beam_name=beam_name,
                energy=energy,
                dose_rate=dose_rate,
                x1=x1 - jaw_padding,
                x2=x2 + jaw_padding,
                y1=y1 - jaw_padding,
                y2=y2 + jaw_padding,
                gantry_angles=axes["gantry"],
                coll_angle=axes["collimator"],
                couch_vrt=couch_vrt,
                couch_lat=couch_lat,
                couch_lng=couch_lng,
                couch_rot=0,
                mlc_positions=mlc.as_control_points(),
                metersets=[mu * m for m in mlc.as_metersets()],
                fluence_mode=fluence_mode,
                is_mlc_hd=self._is_mlc_hd,
            )
            self.add_beam(beam)

    def add_gantry_speed_beams(
        self,
        speeds: tuple[float | int, ...] = (2, 3, 4, 4.8),
        max_dose_rate: int = 600,
        start_gantry_angle: float = 179,
        energy: float = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lat: float = 0,
        couch_lng: float = 1000,
        couch_rot: float = 0,
        beam_name: str = "GS",
        gantry_rot_dir: GantryDirection = GantryDirection.CLOCKWISE,
        jaw_padding_mm: float = 5,
        roi_size_mm: float = 30,
        y1: float = -100,
        y2: float = 100,
        mu: int = 120,
    ):
        """Create a single-image gantry speed test. Multiple speeds are generated. A reference beam is also
        generated. The reference beam is delivered without gantry movement.

        Parameters
        ----------
        speeds : tuple
            The gantry speeds to test. Each speed will have its own ROI.
        max_dose_rate : int
            The max dose rate clinically allowed for the energy.
        start_gantry_angle : float
            The starting gantry angle. The gantry will rotate around this point. It is up to the user
            to know what the machine's limitations are. (i.e. don't go through 180 for Varian machines).
            The ending gantry angle will be the starting angle + the sum of the gantry deltas generated
            by the speed ROIs. Slower speeds require more gantry angle to reach the same MU.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lat : float
            The couch lateral position.
        couch_lng : float
            The couch longitudinal position.
        couch_rot : float
            The couch rotation.
        beam_name : str
            The name of the beam.
        gantry_rot_dir : GantryDirection
            The direction of gantry rotation.
        jaw_padding_mm : float
            The padding to add to the X jaws. The X-jaws will close around the ROIs plus this padding.
        roi_size_mm : float
            The width of the ROIs in mm.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        mu : int
            The monitor units of the beam.

        Notes
        -----

        The gantry angle to cover can be determined via the following:

         gantry speed = gantry_range * max_dose_rate / (MU * 60)

         We can thus solve for the gantry range:

            gantry_range = gantry_speed * MU * 60 / max_dose_rate

        """
        if max(speeds) > self.max_gantry_speed:
            raise ValueError(
                f"Maximum speed given {max(speeds)} is greater than the maximum gantry speed {self.max_gantry_speed}"
            )
        if roi_size_mm * len(speeds) > self.max_overtravel_mm:
            raise ValueError(
                "The ROI size * number of speeds must be less than the overall MLC allowable width"
            )
        # determine sacrifices and gantry angles
        gantry_deltas = [speed * mu * 60 / max_dose_rate for speed in speeds]
        gantry_sign = -1 if gantry_rot_dir == GantryDirection.CLOCKWISE else 1
        g_angles_uncorrected = [start_gantry_angle] + (
            start_gantry_angle + gantry_sign * np.cumsum(gantry_deltas)
        ).tolist()
        gantry_angles = [round(wrap360(angle), 2) for angle in g_angles_uncorrected]

        if sum(gantry_deltas) >= 360:
            raise ValueError(
                "Gantry travel is >360 degrees. Lower the beam MU, use fewer speeds, or decrease the desired gantry speeds"
            )

        mlc = self._create_mlc()
        ref_mlc = self._create_mlc()

        roi_centers = np.linspace(
            -roi_size_mm * len(speeds) / 2 + roi_size_mm / 2,
            roi_size_mm * len(speeds) / 2 - roi_size_mm / 2,
            len(speeds),
        )
        # we have a starting and ending strip
        ref_mlc.add_strip(
            position_mm=float(roi_centers[0]),
            strip_width_mm=roi_size_mm,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=float(roi_centers[0]),
            strip_width_mm=roi_size_mm,
            meterset_at_target=0,
        )
        for center, gantry_angle in zip(roi_centers, gantry_angles):
            ref_mlc.add_strip(
                position_mm=center,
                strip_width_mm=roi_size_mm,
                meterset_at_target=0,
                meterset_transition=1 / len(speeds),
            )
            mlc.add_strip(
                position_mm=center,
                strip_width_mm=roi_size_mm,
                meterset_at_target=0,
                meterset_transition=1 / len(speeds),
            )

        beam = TrueBeamBeam(
            beam_name=beam_name,
            energy=energy,
            dose_rate=max_dose_rate,
            x1=min(roi_centers) - roi_size_mm - jaw_padding_mm,
            x2=max(roi_centers) + roi_size_mm + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angles,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)
        ref_beam = TrueBeamBeam(
            beam_name=f"{beam_name} Ref",
            energy=energy,
            dose_rate=max_dose_rate,
            x1=min(roi_centers) - roi_size_mm - jaw_padding_mm,
            x2=max(roi_centers) + roi_size_mm + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angles[-1],
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=[mu * m for m in ref_mlc.as_metersets()],
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(ref_beam)

    def add_open_field_beam(
        self,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        defined_by_mlcs: bool = True,
        energy: float = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        gantry_angle: float = 0,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lng: float = 1000,
        couch_lat: float = 0,
        couch_rot: float = 0,
        mu: int = 200,
        padding_mm: float = 5,
        beam_name: str = "Open",
        outside_strip_width_mm: float = 5,
    ):
        """Add an open field beam to the plan.

        Parameters
        ----------
        x1 : float
            The left jaw position.
        x2 : float
            The right jaw position.
        y1 : float
            The bottom jaw position.
        y2 : float
            The top jaw position.
        defined_by_mlcs : bool
            Whether the field edges are defined by the MLCs or the jaws.
        energy : float
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        gantry_angle : float
            The gantry angle of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lng : float
            The couch longitudinal position.
        couch_lat : float
            The couch lateral position.
        couch_rot : float
            The couch rotation.
        mu : int
            The monitor units of the beam.
        padding_mm : float
            The padding to add to the jaws or MLCs.
        beam_name : str
            The name of the beam.
        outside_strip_width_mm : float
            The width of the strip of MLCs outside the field. The MLCs will be placed to the
            left, under the X1 jaw by ~2cm.
        """
        if defined_by_mlcs:
            mlc_padding = 0
            jaw_padding = padding_mm
        else:
            mlc_padding = padding_mm
            jaw_padding = 0
        mlc = self._create_mlc()
        mlc.add_rectangle(
            left_position=x1 - mlc_padding,
            right_position=x2 + mlc_padding,
            top_position=y2 + mlc_padding,
            bottom_position=y1 - mlc_padding,
            outer_strip_width=outside_strip_width_mm,
            x_outfield_position=x1 - mlc_padding - jaw_padding - 20,
            meterset_at_target=1.0,
        )
        beam = TrueBeamBeam(
            beam_name=beam_name,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1 - jaw_padding,
            x2=x2 + jaw_padding,
            y1=y1 - jaw_padding,
            y2=y2 + jaw_padding,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=mlc.as_control_points(),
            metersets=[mu * m for m in mlc.as_metersets()],
            fluence_mode=fluence_mode,
            is_mlc_hd=self._is_mlc_hd,
        )
        self.add_beam(beam)


class HalcyonPlanGenerator(PlanGenerator):
    """A class to generate a plan with two beams stacked on top of each other such as the Halcyon. This
    also assumes no jaws."""

    _distal_leaf_boundaries: list[float] = MLC_DISTAL_BOUNDARIES
    _proximal_leaf_bondaries: list[float] = MLC_PROXIMAL_BOUNDARIES

    def __init__(
        self,
        ds: Dataset,
        plan_label: str,
        plan_name: str,
        patient_name: str | None = None,
        patient_id: str | None = None,
        max_mlc_position: float = 140,
        max_mlc_speed: float = 25,
        max_gantry_speed: float = 4.8,
        max_overtravel_mm: float = 140,
    ):
        super().__init__(
            ds,
            plan_label,
            plan_name,
            patient_name,
            patient_id,
            max_mlc_position,
            max_mlc_speed,
            max_gantry_speed,
            max_overtravel_mm,
        )

    def _validate_machine_type(self, beam_sequence: Sequence):
        has_valid_mlc_data: bool = any(
            bld.RTBeamLimitingDeviceType == "MLCX1"
            for bs in beam_sequence
            for bld in bs.BeamLimitingDeviceSequence
        )
        if not has_valid_mlc_data:
            raise ValueError(
                "The machine on the template plan does not seem to be a Halcyon machine."
            )

    def _create_mlc(self) -> tuple[MLCShaper, MLCShaper]:
        """Create 2 MLC shaper objects, one for each stack."""
        proximal_mlc = MLCShaper(
            leaf_y_positions=self._proximal_leaf_bondaries,
            max_mlc_position=self.max_mlc_position,
            max_overtravel_mm=self.max_overtravel_mm,
            sacrifice_gap_mm=None,
            sacrifice_max_move_mm=None,
        )
        distal_mlc = MLCShaper(
            leaf_y_positions=self._distal_leaf_boundaries,
            max_mlc_position=self.max_mlc_position,
            max_overtravel_mm=self.max_overtravel_mm,
            sacrifice_gap_mm=None,
            sacrifice_max_move_mm=None,
        )
        return proximal_mlc, distal_mlc

    def add_picketfence_beam(
        self,
        stack: Stack,
        strip_width_mm: float = 3,
        strip_positions_mm: tuple[float, ...] = (-45, -30, -15, 0, 15, 30, 45),
        gantry_angle: float = 0,
        coll_angle: float = 0,
        couch_vrt: float = 0,
        couch_lng: float = 1000,
        couch_lat: float = 0,
        mu: int = 200,
        beam_name: str = "PF",
    ):
        """Add a picket fence beam to the plan. The beam will be delivered with the MLCs stacked on top of each other.

        Parameters
        ----------
        stack: Stack
            Which MLC stack to use for the beam. The other stack will be parked.
        strip_width_mm : float
            The width of the strips in mm.
        strip_positions_mm : tuple
            The positions of the strips in mm.
        gantry_angle : float
            The gantry angle of the beam.
        coll_angle : float
            The collimator angle of the beam.
        couch_vrt : float
            The couch vertical position.
        couch_lng : float
            The couch longitudinal position.
        couch_lat : float
            The couch lateral position.
        mu : int
            The monitor units of the beam.
        beam_name : str
            The name of the beam.
        """
        prox_mlc, dist_mlc = self._create_mlc()

        # we prepend the positions with an initial starting position 2mm from the first strip
        # that way, each picket is the same cadence where the leaves move into position dynamically.
        # If you didn't do this, the first picket might be different as it has the advantage
        # of starting from a static position vs the rest of the pickets being dynamic.
        strip_positions = [strip_positions_mm[0] - 2, *strip_positions_mm]
        metersets = [0, *[1 / len(strip_positions_mm) for _ in strip_positions_mm]]

        for strip, meterset in zip(strip_positions, metersets):
            if stack in (Stack.DISTAL, Stack.BOTH):
                dist_mlc.add_strip(
                    position_mm=strip,
                    strip_width_mm=strip_width_mm,
                    meterset_at_target=meterset,
                )
                if stack == Stack.DISTAL:
                    prox_mlc.park(meterset=meterset)
            if stack in (Stack.PROXIMAL, Stack.BOTH):
                prox_mlc.add_strip(
                    position_mm=strip,
                    strip_width_mm=strip_width_mm,
                    meterset_at_target=meterset,
                )
                if stack == Stack.PROXIMAL:
                    dist_mlc.park(meterset=meterset)

        beam = HalcyonBeam(
            beam_name=beam_name,
            gantry_angles=gantry_angle,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            proximal_mlc_positions=prox_mlc.as_control_points(),
            distal_mlc_positions=dist_mlc.as_control_points(),
            # can use either MLC for metersets
            metersets=[mu * m for m in prox_mlc.as_metersets()],
        )
        self.add_beam(beam)

    def add_open_field_beam(
        self,
    ):
        """Add an open field beam to the plan. This can be defined by one or both MLC stacks.
        The x1, x2, y1, and y2 terms are colloquial. The MLCs define the field edges.
        The nearest MLC to the field edge in the y-direction will be the one used. It can round up or down.
        """
        raise NotImplementedError(
            "Open field beams are not yet implemented for Halcyon plans"
        )

    def add_dose_rate_beams(
        self,
    ):
        raise NotImplementedError(
            "Dose rate beams are not yet implemented for Halcyon plans"
        )

    def add_mlc_speed_beams(self):
        raise NotImplementedError(
            "MLC speed beams are not yet implemented for Halcyon plans"
        )

    def add_gantry_speed_beams(
        self,
    ):
        raise NotImplementedError(
            "Gantry speed beams are not yet implemented for Halcyon plans"
        )

    def add_winston_lutz_beams(
        self,
    ):
        raise NotImplementedError(
            "Winston-Lutz beams are not yet implemented for Halcyon plans"
        )


class OvertravelError(ValueError):
    pass
