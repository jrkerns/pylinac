from __future__ import annotations

import datetime
import math
from enum import Enum
from pathlib import Path
from typing import Iterable

import numpy as np
import pydicom
from matplotlib.figure import Figure
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from ..core.scale import wrap360
from .fluence import plot_fluences
from .mlc import MLCShaper


class GantryDirection(Enum):
    CLOCKWISE = "CC"
    COUNTER_CLOCKWISE = "CCW"
    NONE = "NONE"


class GantrySpeedTransition(Enum):
    LEADING = "leading"
    TRAILING = "trailing"


class BeamType(Enum):
    DYNAMIC = "DYNAMIC"
    STATIC = "STATIC"


class FluenceMode(Enum):
    STANDARD = "STANDARD"
    FFF = "FFF"
    SRS = "SRS"


class Beam:
    """Represents a DICOM beam dataset. Has methods for creating the dataset and adding control points.
    Generally not created on its own but rather under the hood as part of a PlanGenerator object.

    It contains enough independent logic steps that it's worth separating out from the PlanGenerator class.
    """

    ds: Dataset

    def __init__(
        self,
        plan_dataset: Dataset,
        beam_name: str,
        beam_type: BeamType,
        energy: float,
        dose_rate: int,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        machine_name: str,
        gantry_angles: float | list[float],
        gantry_direction: GantryDirection,
        coll_angle: float,
        couch_vrt: float,
        couch_lat: float,
        couch_lng: float,
        couch_rot: float,
        mlc_boundaries: list[float],
        mlc_positions: list[list[float]],
        metersets: list[float],
        fluence_mode: FluenceMode,
    ):
        """
        Parameters
        ----------
        plan_dataset
            The plan dataset. Used for dynamic links to other Sequences of the plan, such as Dose Reference and Tolerance Table.
        beam_name : str
            The name of the beam. Must be less than 16 characters.
        beam_type : BeamType
            The type of beam: dynamic or static.
        energy : float
            The energy of the beam.
        dose_rate : int
            The dose rate of the beam.
        x1 : float
            The left jaw position.
        x2 : float
            The right jaw position.
        y1 : float
            The bottom jaw position.
        y2 : float
            The top jaw position.
        machine_name : str
            The name of the machine.
        gantry_angles : Union[float, list[float]]
            The gantry angle(s) of the beam. If a single number, it's assumed to be a static beam. If multiple numbers, it's assumed to be a dynamic beam.
        gantry_direction : GantryDirection
            The direction of the gantry rotation. Only relevant if multiple gantry angles are specified.
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
        mlc_boundaries : list[float]
            The MLC boundaries. These are the same thing as the LeafPositionBoundaries in the DICOM file.
        mlc_positions : list[list[float]]
            The MLC positions for each control point. This is the x-position of each leaf for each control point.
        metersets : list[float]
            The meter sets for each control point. The length must match the number of control points in mlc_positions.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        """
        if len(beam_name) > 16:
            raise ValueError("Beam name must be less than or equal to 16 characters")

        if beam_type == BeamType.STATIC and len(mlc_positions) != 1:
            raise ValueError(
                f"Static beam can only have one MLC position change ({len(mlc_positions)})"
            )

        if len(metersets) != len(mlc_positions):
            raise ValueError(
                f"The number of meter sets ({len(metersets)}) "
                f"must match the number of MLC position changes ({len(mlc_positions)})"
            )

        if (
            isinstance(gantry_angles, Iterable)
            and gantry_direction == GantryDirection.NONE
        ):
            raise ValueError(
                "Cannot specify multiple gantry angles without a gantry direction"
            )

        if isinstance(gantry_angles, Iterable) and beam_type == BeamType.STATIC:
            raise ValueError("Cannot specify multiple gantry angles as a static beam")
        self.plan_ds = plan_dataset
        self.ds = self._create_basic_beam_info(
            beam_name,
            beam_type,
            fluence_mode,
            machine_name,
            num_leaves=len(mlc_positions[0]),
            mlc_boundaries=mlc_boundaries,
        )
        if not isinstance(
            gantry_angles, Iterable
        ):  # if it's just a single number (like for a static beam) set it to an array of that value
            gantry_angles = [gantry_angles] * len(metersets)
        self._append_initial_control_point(
            energy,
            dose_rate,
            x1,
            x2,
            y1,
            y2,
            gantry_angles[0],
            gantry_direction,
            coll_angle,
            couch_vrt,
            couch_lat,
            couch_lng,
            couch_rot,
            mlc_positions=mlc_positions[0],
        )
        for mlc_pos, meterset, gantry_angle in zip(
            mlc_positions[1:], metersets[1:], gantry_angles[1:]
        ):
            if beam_type == BeamType.DYNAMIC:
                self._append_secondary_dynamic_control_point(
                    mlc_positions=mlc_pos,
                    meterset=meterset,
                    gantry_angle=gantry_angle,
                    gantry_dir=gantry_direction,
                )
            else:
                self._append_secondary_static_control_point(meterset=meterset)

    def _create_basic_beam_info(
        self,
        beam_name: str,
        beam_type: BeamType,
        fluence_mode: FluenceMode,
        machine_name: str,
        num_leaves: int,
        mlc_boundaries: list[float],
    ) -> Dataset:
        beam = Dataset()
        beam.Manufacturer = "Radformation"
        beam.ManufacturerModelName = "RadMachine"
        beam.TreatmentMachineName = machine_name
        beam.PrimaryDosimeterUnit = "MU"
        beam.SourceAxisDistance = "1000.0"

        # Primary Fluence Mode Sequence
        primary_fluence_mode_sequence = Sequence()
        beam.PrimaryFluenceModeSequence = primary_fluence_mode_sequence

        primary_fluence_mode1 = Dataset()
        if fluence_mode == FluenceMode.STANDARD:
            primary_fluence_mode1.FluenceMode = "STANDARD"
        elif fluence_mode == FluenceMode.FFF:
            primary_fluence_mode1.FluenceMode = "NON_STANDARD"
            primary_fluence_mode1.FluenceModeID = "FFF"
        elif fluence_mode == FluenceMode.SRS:
            primary_fluence_mode1.FluenceMode = "NON_STANDARD"
            primary_fluence_mode1.FluenceModeID = "SRS"

        primary_fluence_mode_sequence.append(primary_fluence_mode1)

        # Beam Limiting Device Sequence
        beam_limiting_device_sequence = Sequence()
        beam.BeamLimitingDeviceSequence = beam_limiting_device_sequence

        # Beam Limiting Device Sequence: Beam Limiting Device 1
        beam_limiting_device1 = Dataset()
        beam_limiting_device1.RTBeamLimitingDeviceType = "ASYMX"
        beam_limiting_device1.NumberOfLeafJawPairs = "1"
        beam_limiting_device_sequence.append(beam_limiting_device1)

        # Beam Limiting Device Sequence: Beam Limiting Device 2
        beam_limiting_device2 = Dataset()
        # TODO: likely that Elekta will need tweaking for this
        beam_limiting_device2.RTBeamLimitingDeviceType = "ASYMY"
        beam_limiting_device2.NumberOfLeafJawPairs = "1"
        beam_limiting_device_sequence.append(beam_limiting_device2)

        # Beam Limiting Device Sequence: Beam Limiting Device 3
        beam_limiting_device3 = Dataset()

        beam_limiting_device3.RTBeamLimitingDeviceType = "MLCX"
        beam_limiting_device3.NumberOfLeafJawPairs = int(num_leaves / 2)
        beam_limiting_device3.LeafPositionBoundaries = mlc_boundaries

        beam_limiting_device_sequence.append(beam_limiting_device3)

        # beam numbers start at 0 and increment from there.
        beam.BeamNumber = len(self.plan_ds.BeamSequence)
        beam.BeamName = beam_name
        beam.BeamType = beam_type.value
        beam.RadiationType = "PHOTON"
        beam.TreatmentDeliveryType = "TREATMENT"
        beam.NumberOfWedges = "0"
        beam.NumberOfCompensators = "0"
        beam.NumberOfBoli = "0"
        beam.NumberOfBlocks = "0"
        beam.FinalCumulativeMetersetWeight = "1.0"
        beam.NumberOfControlPoints = "0"

        # Control Point Sequence
        cp_sequence = Sequence()
        beam.ControlPointSequence = cp_sequence
        return beam

    def _append_initial_control_point(
        self,
        energy: int,
        dose_rate: int,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        gantry_angle: float | list[float],
        gantry_rot: GantryDirection,
        coll_angle: float,
        couch_vrt: float,
        couch_lat: float,
        couch_lng: float,
        couch_rot: float,
        mlc_positions: list[float],
    ):
        # Control Point Sequence: Control Point 0
        cp0 = Dataset()
        cp0.ControlPointIndex = "0"
        cp0.NominalBeamEnergy = str(energy)
        cp0.DoseRateSet = str(dose_rate)

        # Beam Limiting Device Position Sequence
        beam_limiting_device_position_sequence = Sequence()
        cp0.BeamLimitingDevicePositionSequence = beam_limiting_device_position_sequence

        # Beam Limiting Device Position Sequence: Beam Limiting Device Position 1
        beam_limiting_device_position1 = Dataset()
        beam_limiting_device_position1.RTBeamLimitingDeviceType = "ASYMX"
        beam_limiting_device_position1.LeafJawPositions = [x1, x2]
        beam_limiting_device_position_sequence.append(beam_limiting_device_position1)

        # Beam Limiting Device Position Sequence: Beam Limiting Device Position 2
        beam_limiting_device_position2 = Dataset()
        beam_limiting_device_position2.RTBeamLimitingDeviceType = "ASYMY"
        beam_limiting_device_position2.LeafJawPositions = [y1, y2]
        beam_limiting_device_position_sequence.append(beam_limiting_device_position2)

        # Beam Limiting Device Position Sequence: Beam Limiting Device Position 3
        beam_limiting_device_position3 = Dataset()
        beam_limiting_device_position3.RTBeamLimitingDeviceType = "MLCX"
        beam_limiting_device_position3.LeafJawPositions = [
            f"{m:6f}" for m in mlc_positions
        ]  # convert to truncated string to fit VR limitations
        beam_limiting_device_position_sequence.append(beam_limiting_device_position3)

        cp0.GantryAngle = str(gantry_angle)
        cp0.GantryRotationDirection = gantry_rot.value
        cp0.BeamLimitingDeviceAngle = str(coll_angle)
        cp0.BeamLimitingDeviceRotationDirection = "NONE"
        cp0.PatientSupportAngle = str(couch_rot)
        cp0.PatientSupportRotationDirection = "NONE"
        cp0.TableTopEccentricAngle = "0.0"
        cp0.TableTopEccentricRotationDirection = "NONE"
        cp0.TableTopVerticalPosition = str(couch_vrt)
        cp0.TableTopLongitudinalPosition = str(couch_lng)
        cp0.TableTopLateralPosition = str(couch_lat)
        cp0.IsocenterPosition = None
        cp0.CumulativeMetersetWeight = "0.0"

        # Referenced Dose Reference Sequence
        refd_dose_ref_sequence = Sequence()
        cp0.ReferencedDoseReferenceSequence = refd_dose_ref_sequence

        # linked setup number and tolerance table
        # we *could* hardcode it but this makes it more flexible to changes upstream
        self.ds.ReferencedPatientSetupNumber = self.plan_ds.PatientSetupSequence[
            0
        ].PatientSetupNumber
        self.ds.ReferencedToleranceTableNumber = self.plan_ds.ToleranceTableSequence[
            0
        ].ToleranceTableNumber

        self.ds.NumberOfControlPoints = 1  # increment this
        self.ds.ControlPointSequence.append(cp0)

    def _append_secondary_static_control_point(self, meterset: float) -> None:
        # Control Point Sequence: Control Point 1
        cp1 = Dataset()
        cp1.ControlPointIndex = len(self.ds.ControlPointSequence)

        cp1.CumulativeMetersetWeight = (
            f"{meterset:.5f}"  # convert to truncated string to fit VR limitations
        )

        self.ds.NumberOfControlPoints = int(self.ds.NumberOfControlPoints) + 1
        self.ds.ControlPointSequence.append(cp1)

    def _append_secondary_dynamic_control_point(
        self,
        mlc_positions: list[float],
        meterset: float,
        gantry_angle: float,
        gantry_dir: GantryDirection,
    ) -> None:
        # Control Point Sequence: Control Point 1
        cp1 = Dataset()
        cp1.ControlPointIndex = len(self.ds.ControlPointSequence)

        if gantry_dir != GantryDirection.NONE:
            cp1.GantryAngle = gantry_angle
            cp1.GantryRotationDirection = gantry_dir.value
        cp1.CumulativeMetersetWeight = (
            f"{meterset:.5f}"  # convert to truncated string to fit VR limitations
        )

        # Beam Limiting Device Position Sequence
        beam_limiting_device_position_sequence = Sequence()
        cp1.BeamLimitingDevicePositionSequence = beam_limiting_device_position_sequence

        # Beam Limiting Device Position Sequence: Beam Limiting Device Position 1
        beam_limiting_device_position1 = Dataset()
        beam_limiting_device_position1.RTBeamLimitingDeviceType = "MLCX"
        beam_limiting_device_position1.LeafJawPositions = [
            f"{m:6f}" for m in mlc_positions
        ]  # convert to truncated string to fit VR limitations
        beam_limiting_device_position_sequence.append(beam_limiting_device_position1)

        # Referenced Dose Reference Sequence
        refd_dose_ref_sequence = Sequence()
        cp1.ReferencedDoseReferenceSequence = refd_dose_ref_sequence

        self.ds.NumberOfControlPoints = int(self.ds.NumberOfControlPoints) + 1
        self.ds.ControlPointSequence.append(cp1)

    def as_dicom(self) -> Dataset:
        """Return the beam as a DICOM dataset that represents a BeamSequence item."""
        return self.ds


class PlanGenerator:
    def __init__(
        self,
        ds: Dataset,
        plan_label: str,
        plan_name: str,
        x_width_mm: float = 400,
        max_mlc_speed: float = 25,
        max_gantry_speed: float = 4.8,
        sacrificial_gap_mm: float = 5,
        max_sacrificial_move_mm: float = 50,
        max_overtravel_mm: float = 140,
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
        x_width_mm : float
            The overall width of the MLC movement in the x-direction. Generally, this is the x field size.
        max_mlc_speed : float
            The maximum speed of the MLC leaves in mm/s
        max_gantry_speed : float
            The maximum speed of the gantry in degrees/s.
        sacrificial_gap_mm : float
            For certain dynamic beams, the top and bottom leaf pair are used to slow axes down. This is the gap
            between those leaves at any given time.
        max_sacrificial_move_mm : float
            The maximum distance the sacrificial leaves can move in a given control point.
            Smaller values generate more control points and more back-and-forth movement.
            Too large of values may cause deliverability issues.
        max_overtravel_mm : float
            The maximum distance the MLC leaves can overtravel from each other as well as the jaw size (for tail exposure protection).
        """
        if ds.Modality != "RTPLAN":
            raise ValueError("File is not an RTPLAN file")
        self.max_overtravel_mm = max_overtravel_mm
        self.x_width = x_width_mm
        self.max_mlc_speed = max_mlc_speed
        self.max_gantry_speed = max_gantry_speed
        self.sacrificial_gap = sacrificial_gap_mm
        self.max_sacrificial_move = max_sacrificial_move_mm
        if not hasattr(ds, "PatientName") or not hasattr(ds, "PatientID"):
            raise ValueError("RTPLAN file must have PatientName and PatientID")
        if not hasattr(ds, "ToleranceTableSequence"):
            raise ValueError("RTPLAN file must have ToleranceTableSequence")
        if not hasattr(ds, "BeamSequence") or not hasattr(
            ds.BeamSequence[0].BeamLimitingDeviceSequence[-1], "LeafPositionBoundaries"
        ):
            raise ValueError("RTPLAN file must have MLC data")

        self.ds = ds

        ######  Clear/initialize the metadata for the new plan

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
        patient_setup.PatientSetupNumber = "0"
        self.ds.PatientSetupSequence = Sequence((patient_setup,))

        # Dose Reference Sequence
        dose_ref_sequence = Sequence()
        self.ds.DoseReferenceSequence = dose_ref_sequence
        dose_ref1 = Dataset()
        dose_ref1.DoseReferenceNumber = "1"
        dose_ref1.DoseReferenceUID = generate_uid()
        dose_ref1.DoseReferenceStructureType = "SITE"
        dose_ref1.DoseReferenceDescription = "PTV"
        dose_ref1.DoseReferenceType = "TARGET"
        dose_ref1.DeliveryMaximumDose = "20.0"
        dose_ref1.TargetPrescriptionDose = "40.0"
        dose_ref1.TargetMaximumDose = "20.0"
        self.ds.DoseReferenceSequence.append(dose_ref1)

        # Fraction Group Sequence
        frxn_gp_sequence = Sequence()
        self.ds.FractionGroupSequence = frxn_gp_sequence
        # Fraction Group Sequence: Fraction Group 1
        frxn_gp1 = Dataset()
        frxn_gp1.FractionGroupNumber = "1"
        frxn_gp1.NumberOfFractionsPlanned = "1"
        frxn_gp1.NumberOfBeams = str(0)
        frxn_gp1.NumberOfBrachyApplicationSetups = "0"

        # Referenced Beam Sequence
        refd_beam_sequence = Sequence()
        frxn_gp1.ReferencedBeamSequence = refd_beam_sequence
        self.ds.FractionGroupSequence.append(frxn_gp1)

        # Beam Sequence
        self._old_beam = self.ds.BeamSequence[0]
        beam_sequence = Sequence()
        self.ds.BeamSequence = beam_sequence

    @classmethod
    def from_rt_plan_file(cls, rt_plan_file: str, **kwargs) -> PlanGenerator:
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

    @property
    def num_leaves(self) -> int:
        """The number of leaves in the MLC."""
        return self._old_beam.BeamLimitingDeviceSequence[-1].NumberOfLeafJawPairs * 2

    @property
    def leaf_config(self) -> list[float]:
        """The leaf boundaries of the MLC."""
        return self._old_beam.BeamLimitingDeviceSequence[-1].LeafPositionBoundaries

    @property
    def machine_name(self) -> str:
        """The name of the machine."""
        return self._old_beam.TreatmentMachineName

    def _create_mlc(self):
        """Utility to create MLC shaper instances."""
        return MLCShaper(
            leaf_y_positions=self.leaf_config,
            max_x_mm=self.x_width / 2,
            sacrifice_gap_mm=self.sacrificial_gap,
            sacrifice_max_move_mm=self.max_sacrificial_move,
            max_overtravel_mm=self.max_overtravel_mm,
        )

    def add_beam(self, beam_dataset: Dataset, mu: int):
        """Add a beam to the plan using the Beam object. Although public,
        this is a low-level method that is used by the higher-level methods like add_open_field_beam.
        This handles the associated metadata like the referenced beam sequence and fraction group sequence.
        """
        self.ds.BeamSequence.append(beam_dataset)
        self._update_references(mu=mu)

    def add_picketfence_beam(
        self,
        strip_width_mm: float = 3,
        strip_positions_mm: tuple[float] = (-45, -30, -15, 0, 15, 30, 45),
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
        mlc = self._create_mlc()
        # create initial starting point; start under the jaws
        mlc.add_strip(
            position_mm=strip_positions_mm[0] + 2,
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
        beam = Beam(
            plan_dataset=self.ds,
            beam_name=beam_name,
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1,
            x2=x2,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=mlc.as_control_points(),
            metersets=mlc.as_metersets(),
            fluence_mode=fluence_mode,
            mlc_boundaries=self.leaf_config,
            machine_name=self.machine_name,
        )
        self.add_beam(beam.as_dicom(), mu=mu)

    def add_dose_rate_beams(
        self,
        dose_rates: tuple[int] = (100, 300, 500, 600),
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

        mlc = self._create_mlc()
        ref_mlc = self._create_mlc()

        roi_centers = np.linspace(
            -roi_size_mm * len(dose_rates) / 2 + roi_size_mm / 2,
            roi_size_mm * len(dose_rates) / 2 - roi_size_mm / 2,
            len(dose_rates),
        )
        # we have a starting and ending strip
        ref_mlc.add_strip(
            position_mm=roi_centers[0] - roi_size_mm / 2,
            strip_width_mm=0,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=roi_centers[0] - roi_size_mm / 2,
            strip_width_mm=0,
            meterset_at_target=0,
            initial_sacrificial_gap_mm=5,
        )
        for sacrifice_distance, center in zip(sacrificial_movements, roi_centers):
            ref_mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,
                top_position=max(self.leaf_config),
                bottom_position=min(self.leaf_config),
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
                top_position=max(self.leaf_config),
                bottom_position=min(self.leaf_config),
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
        ref_beam = Beam(
            plan_dataset=self.ds,
            beam_name="DR Ref",
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=default_dose_rate,
            x1=roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm,
            x2=roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=ref_mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(ref_beam.as_dicom(), mu=mu)
        beam = Beam(
            plan_dataset=self.ds,
            beam_name=f"DR{min(dose_rates)}-{max(dose_rates)}",
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=default_dose_rate,
            x1=roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm,
            x2=roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(beam.as_dicom(), mu=mu)

    def add_mlc_speed_beams(
        self,
        speeds: tuple[float] = (5, 10, 15, 20),
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

        mlc = self._create_mlc()
        ref_mlc = self._create_mlc()

        roi_centers = np.linspace(
            -roi_size_mm * len(speeds) / 2 + roi_size_mm / 2,
            roi_size_mm * len(speeds) / 2 - roi_size_mm / 2,
            len(speeds),
        )
        # we have a starting and ending strip
        ref_mlc.add_strip(
            position_mm=roi_centers[0] - roi_size_mm / 2,
            strip_width_mm=0,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=roi_centers[0] - roi_size_mm / 2,
            strip_width_mm=0,
            meterset_at_target=0,
            initial_sacrificial_gap_mm=5,
        )
        for sacrifice_distance, center in zip(sacrificial_movements, roi_centers):
            ref_mlc.add_rectangle(
                left_position=center - roi_size_mm / 2,
                right_position=center + roi_size_mm / 2,
                x_outfield_position=-200,  # not relevant
                top_position=max(self.leaf_config),
                bottom_position=min(self.leaf_config),
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
                top_position=max(self.leaf_config),
                bottom_position=min(self.leaf_config),
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
        ref_beam = Beam(
            plan_dataset=self.ds,
            beam_name="MLC Sp Ref",
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=default_dose_rate,
            x1=roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm,
            x2=roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=ref_mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(ref_beam.as_dicom(), mu=mu)
        beam = Beam(
            plan_dataset=self.ds,
            beam_name=beam_name,
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=default_dose_rate,
            x1=roi_centers[0] - roi_size_mm / 2 - jaw_padding_mm,
            x2=roi_centers[-1] + roi_size_mm / 2 + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(beam.as_dicom(), mu=mu)

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
            The positions of the axes. Each dict should have keys 'gantry', 'collimator', and 'couch'.
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
            beam = Beam(
                plan_dataset=self.ds,
                beam_name=f"G{axes['gantry']:g}C{axes['collimator']:g}P{axes['couch']:g}",
                beam_type=BeamType.DYNAMIC,
                energy=energy,
                dose_rate=dose_rate,
                x1=x1 - jaw_padding,
                x2=x2 + jaw_padding,
                y1=y1 - jaw_padding,
                y2=y2 + jaw_padding,
                gantry_angles=axes["gantry"],
                gantry_direction=GantryDirection.NONE,
                coll_angle=axes["collimator"],
                couch_vrt=couch_vrt,
                couch_lat=couch_lat,
                couch_lng=couch_lng,
                couch_rot=0,
                mlc_positions=mlc.as_control_points(),
                metersets=mlc.as_metersets(),
                fluence_mode=fluence_mode,
                mlc_boundaries=self.leaf_config,
                machine_name=self.machine_name,
            )
            self.add_beam(beam.as_dicom(), mu=mu)

    def add_gantry_speed_beams(
        self,
        speeds: tuple = (2, 3, 4, 4.8),
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
        gantry_angles = [wrap360(angle) for angle in g_angles_uncorrected]

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
            position_mm=roi_centers[0],
            strip_width_mm=roi_size_mm,
            meterset_at_target=0,
        )
        mlc.add_strip(
            position_mm=roi_centers[0],
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

        beam = Beam(
            plan_dataset=self.ds,
            beam_name=beam_name,
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=max_dose_rate,
            x1=min(roi_centers) - roi_size_mm - jaw_padding_mm,
            x2=max(roi_centers) + roi_size_mm + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angles,
            gantry_direction=gantry_rot_dir,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=mlc.as_control_points(),
            metersets=mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(beam.as_dicom(), mu=mu)
        ref_beam = Beam(
            plan_dataset=self.ds,
            beam_name="G Sp Ref",
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=max_dose_rate,
            x1=min(roi_centers) - roi_size_mm - jaw_padding_mm,
            x2=max(roi_centers) + roi_size_mm + jaw_padding_mm,
            y1=y1,
            y2=y2,
            gantry_angles=gantry_angles[-1],
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            machine_name=self.machine_name,
            fluence_mode=fluence_mode,
            mlc_positions=ref_mlc.as_control_points(),
            metersets=ref_mlc.as_metersets(),
            mlc_boundaries=self.leaf_config,
        )
        self.add_beam(ref_beam.as_dicom(), mu=mu)

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
        beam = Beam(
            plan_dataset=self.ds,
            beam_name=beam_name,
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=dose_rate,
            x1=x1 - jaw_padding,
            x2=x2 + jaw_padding,
            y1=y1 - jaw_padding,
            y2=y2 + jaw_padding,
            gantry_angles=gantry_angle,
            gantry_direction=GantryDirection.NONE,
            coll_angle=coll_angle,
            couch_vrt=couch_vrt,
            couch_lat=couch_lat,
            couch_lng=couch_lng,
            couch_rot=couch_rot,
            mlc_positions=mlc.as_control_points(),
            metersets=mlc.as_metersets(),
            fluence_mode=fluence_mode,
            mlc_boundaries=self.leaf_config,
            machine_name=self.machine_name,
        )
        self.add_beam(beam.as_dicom(), mu=mu)

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

    def _update_references(self, mu: float) -> None:
        """Update the other sequences that reference the beam sequence."""
        referenced_beam = Dataset()
        referenced_beam.BeamDose = "1.0"
        referenced_beam.BeamMeterset = str(mu)
        referenced_beam.ReferencedDoseReferenceUID = self.ds.DoseReferenceSequence[
            0
        ].DoseReferenceUID
        self.ds.FractionGroupSequence[0].ReferencedBeamSequence.append(referenced_beam)
        referenced_beam.ReferencedBeamNumber = (
            len(self.ds.FractionGroupSequence[0].ReferencedBeamSequence) - 1
        )
        # increment number of beams
        self.ds.FractionGroupSequence[0].NumberOfBeams = (
            int(self.ds.FractionGroupSequence[0].NumberOfBeams) + 1
        )
