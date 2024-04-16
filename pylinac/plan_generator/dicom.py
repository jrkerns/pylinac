import datetime
from enum import Enum
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Union

import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

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
        beam_name: str,
        beam_type: BeamType,
        energy: int,
        dose_rate: int,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        machine_name: str,
        gantry_angles: Union[float, list[float]],
        gantry_direction: GantryDirection,
        coll_angle: int,
        couch_vrt: int,
        couch_lat: int,
        couch_lng: int,
        couch_rot: int,
        mlc_boundaries: list[float],
        mlc_positions: list[list[float]],
        meter_sets: list[float],
        fluence_mode: FluenceMode,
    ):
        """
        Parameters
        ----------

        beam_name : str
            The name of the beam. Must be less than 16 characters.
        beam_type : BeamType
            The type of beam: dynamic or static.
        energy : int
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
        coll_angle : int
            The collimator angle.
        couch_vrt : int
            The couch vertical position.
        couch_lat : int
            The couch lateral position.
        couch_lng : int
            The couch longitudinal position.
        couch_rot : int
            The couch rotation.
        mlc_boundaries : list[float]
            The MLC boundaries. These are the same thing as the LeafPositionBoundaries in the DICOM file.
        mlc_positions : list[list[float]]
            The MLC positions for each control point. This is the x-position of each leaf for each control point.
        meter_sets : list[float]
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

        if beam_type != BeamType.STATIC and len(meter_sets) != len(mlc_positions):
            raise ValueError(
                f"The number of meter sets ({len(meter_sets)}) "
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
            gantry_angles = [gantry_angles] * len(meter_sets)
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
        for mlc_pos, meter_set, gantry_angle in zip_longest(
            mlc_positions[1:], meter_sets[1:], gantry_angles[1:]
        ):
            if beam_type == BeamType.DYNAMIC:
                self._append_secondary_dynamic_control_point(
                    mlc_positions=mlc_pos,
                    meter_set=meter_set,
                    gantry_angle=gantry_angle,
                    gantry_dir=gantry_direction,
                )
            else:
                self._append_secondary_static_control_point(meter_set=meter_set)

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

        # TODO: use MLC config instead
        beam_limiting_device3.RTBeamLimitingDeviceType = "MLCX"
        beam_limiting_device3.NumberOfLeafJawPairs = num_leaves
        beam_limiting_device3.LeafPositionBoundaries = mlc_boundaries

        beam_limiting_device_sequence.append(beam_limiting_device3)

        # beam1.BeamNumber = len(self.ds.BeamSequence)  # set later in plan;
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
        gantry_angle: Union[float, list[float]],
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

        self.ds.NumberOfControlPoints = 1  # increment this
        self.ds.ControlPointSequence.append(cp0)

    def _append_secondary_static_control_point(self, meter_set: float) -> None:
        # Control Point Sequence: Control Point 1
        cp1 = Dataset()
        cp1.ControlPointIndex = len(self.ds.ControlPointSequence)

        cp1.CumulativeMetersetWeight = (
            f"{meter_set:.5f}"  # convert to truncated string to fit VR limitations
        )

        self.ds.NumberOfControlPoints = int(self.ds.NumberOfControlPoints) + 1
        self.ds.ControlPointSequence.append(cp1)

    def _append_secondary_dynamic_control_point(
        self,
        mlc_positions: list[float],
        meter_set: float,
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
            f"{meter_set:.5f}"  # convert to truncated string to fit VR limitations
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


class PlanGenerator:
    @classmethod
    def from_rt_plan_file(cls, rt_plan_file: str, plan_label: str, plan_name: str):
        """Load an existing RTPLAN file and create a new plan based on it.

        Parameters
        ----------
        rt_plan_file : str
            The path to the RTPLAN file.
        plan_label : str
            The label of the new plan.
        plan_name : str
            The name of the new plan.
        """
        ds = pydicom.dcmread(rt_plan_file)
        return cls(ds, plan_label, plan_name)

    def __init__(
        self, ds: Dataset, plan_label: str, plan_name: str, x_width: float = 400
    ):
        """A tool for generating new QA RTPlan files.

        Parameters
        ----------
        ds : Dataset
              The RTPLAN dataset to base the new plan off of. The plan must already have MLC positions.
        plan_label : str
            The label of the new plan.
        plan_name : str
            The name of the new plan.
        x_width : float
            The overall width of the MLC movement in the x-direction.
        """
        if ds.Modality != "RTPLAN":
            raise ValueError("File is not an RTPLAN file")
        self.x_width = x_width
        # check has patient name, id

        # check has tolerance table sequence

        # check has MLC data

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
        patient_setup_sequence = Sequence()
        self.ds.PatientSetupSequence = patient_setup_sequence

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

    def add_beam(self, beam: Beam, mu: int):
        """Add a beam to the plan using the Beam object. Although public,
        this is a low-level method that is used by the higher-level methods like add_open_field_beam.
        This handles the associated metadata like the referenced beam sequence and fraction group sequence.
        """
        self.ds.BeamSequence.append(beam.ds)
        self.ds.BeamSequence[-1].BeamNumber = len(self.ds.BeamSequence)
        self._add_referenced_beam(mu=mu)
        self.ds.FractionGroupSequence[0].NumberOfBeams = (
            int(self.ds.FractionGroupSequence[0].NumberOfBeams) + 1
        )

    def add_picketfence_beam(
        self,
        strip_width: float = 3,
        strip_positions: tuple = (-50, -30, -10, 10, 30, 50),
        y1: float = -200,
        y2: float = 200,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        energy: int = 6,
        gantry_angle: int = 0,
        coll_angle: int = 0,
        couch_vrt: int = 0,
        couch_lng: int = 100,
        couch_lat: int = 0,
        couch_rot: int = 0,
        beam_mu: int = 200,
        transition_dose: float = 0.01,
        jaw_padding: int = 20,
        beam_name: str = "PF",
    ):
        """Add a picket fence beam to the plan.

        Parameters
        ----------
        strip_width : float
            The width of the strips in mm.
        strip_positions : tuple
            The positions of the strips in mm.
        y1 : float
            The bottom jaw position. Usually negative. More negative is lower.
        y2 : float
            The top jaw position. Usually positive. More positive is higher.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        energy : int
            The energy of the beam.
        gantry_angle : int
            The gantry angle of the beam.
        coll_angle : int
            The collimator angle of the beam.
        couch_vrt : int
            The couch vertical position.
        couch_lng : int
            The couch longitudinal position.
        couch_lat : int
            The couch lateral position.
        couch_rot : int
            The couch rotation.
        beam_mu : int
            The monitor units of the beam.
        transition_dose : float
            The dose in between the strips. Usually you need *some* dose in between the strips.
        jaw_padding : int
            The padding to add to the X jaws.
        beam_name : str
            The name of the beam.
        """
        # create MU weighting; we create a [0, 1/N, ..., 1] where N is the number of strips sequence but also add in transition doses
        start_meter_sets = np.linspace(0, 1, len(strip_positions) + 1)
        end_meter_sets = np.linspace(0, 1, len(strip_positions) + 1)
        meter_sets = np.concatenate((start_meter_sets, end_meter_sets))[:-1]
        meter_sets.sort()
        meter_sets[1:-1:2] = [
            s + transition_dose for idx, s in enumerate(meter_sets[1:-1:2])
        ]
        mlc = MLCShaper(
            leaf_y_positions=self.leaf_config, max_x=self.x_width / 2, sacrifice_gap=5
        )
        # create initial starting point; start under the jaws
        mlc.add_strip(
            position=strip_positions[0] - jaw_padding - 10,
            strip_width=strip_width,
        )

        for strip in strip_positions:
            # starting control point
            mlc.add_strip(
                position=strip,
                strip_width=strip_width,
            )
            # ending control point
            mlc.add_strip(
                position=strip,
                strip_width=strip_width,
            )
        beam = Beam(
            beam_name=beam_name,
            beam_type=BeamType.DYNAMIC,
            energy=energy,
            dose_rate=dose_rate,
            x1=min(strip_positions) - jaw_padding,
            x2=max(strip_positions) + jaw_padding,
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
            meter_sets=meter_sets,
            beam_mu=beam_mu,
            fluence_mode=fluence_mode,
            mlc_boundaries=self.leaf_config,
            machine_name=self.machine_name,
        )
        self.add_beam(beam, mu=beam_mu)

    def add_dose_rate_beams(self):
        pass

    def add_mlc_speed_beams(self):
        pass

    def add_winston_lutz_beams(
        self,
        x1: float = -10,
        x2: float = 10,
        y1: float = -10,
        y2: float = 10,
        defined_by_mlcs: bool = True,
        energy: int = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        axes_positions: Iterable[dict] = ({"gantry": 0, "collimator": 0, "couch": 0},),
        couch_vrt: int = 0,
        couch_lng: int = 100,
        couch_lat: int = 0,
        beam_mu: int = 10,
        jaw_padding: int = 5,
    ):
        """Add Winston-Lutz beams to the plan. Will create a beam for each set of axes positions.

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
        energy : int
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        axes_positions : Iterable[dict]
            The positions of the axes. Each dict should have keys 'gantry', 'collimator', and 'couch'.
        couch_vrt : int
            The couch vertical position.
        couch_lng : int
            The couch longitudinal position.
        couch_lat : int
            The couch lateral position.
        beam_mu : int
            The monitor units of the beam.
        jaw_padding : int
            The padding to add. If defined by the MLCs, this is the padding of the jaws. If defined by the jaws,
            this is the padding of the MLCs.
        """
        for axes in axes_positions:
            if defined_by_mlcs:
                mlc_padding = 0
                jaw_padding = jaw_padding
            else:
                mlc_padding = jaw_padding
                jaw_padding = 0
            mlc = MLCShaper(
                leaf_y_positions=self.leaf_config, max_x=200, sacrifice_gap=5
            )
            mlc.add_rectangle(
                left_position=x1 - mlc_padding,
                right_position=x2 + mlc_padding,
                top_position=y2 + mlc_padding,
                bottom_position=y1 - mlc_padding,
                outer_strip_width=5,
                x_outfield_position=-200 + 5,
            )
            beam = Beam(
                beam_name=f"G{axes['gantry']:g}C{axes['collimator']:g}P{axes['couch']:g}",
                beam_type=BeamType.STATIC,
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
                meter_sets=[0, 1],
                fluence_mode=fluence_mode,
                mlc_boundaries=self.leaf_config,
                machine_name=self.machine_name,
            )
            self.add_beam(beam, mu=beam_mu)

    def add_gantry_speed_beams(self):
        pass

    def add_dlg_beam(self):
        pass

    def add_open_field_beam(
        self,
        x1: float,
        x2: float,
        y1: float,
        y2: float,
        defined_by_mlcs: bool = True,
        energy: int = 6,
        fluence_mode: FluenceMode = FluenceMode.STANDARD,
        dose_rate: int = 600,
        gantry_angle: int = 0,
        coll_angle: int = 0,
        couch_vrt: int = 0,
        couch_lng: int = 100,
        couch_lat: int = 0,
        couch_rot: int = 0,
        beam_mu: int = 200,
        padding: float = 5,
        beam_name: str = "Open",
        outside_strip_width: float = 5,
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
        energy : int
            The energy of the beam.
        fluence_mode : FluenceMode
            The fluence mode of the beam.
        dose_rate : int
            The dose rate of the beam.
        gantry_angle : int
            The gantry angle of the beam.
        coll_angle : int
            The collimator angle of the beam.
        couch_vrt : int
            The couch vertical position.
        couch_lng : int
            The couch longitudinal position.
        couch_lat : int
            The couch lateral position.
        couch_rot : int
            The couch rotation.
        beam_mu : int
            The monitor units of the beam.
        padding : float
            The padding to add to the jaws or MLCs.
        beam_name : str
            The name of the beam.
        outside_strip_width : float
            The width of the strip of MLCs outside the field.
        """
        if defined_by_mlcs:
            mlc_padding = 0
            jaw_padding = padding
        else:
            mlc_padding = padding
            jaw_padding = 0
        mlc = MLCShaper(leaf_y_positions=self.leaf_config, max_x=200, sacrifice_gap=5)
        mlc.add_rectangle(
            left_position=x1 - mlc_padding,
            right_position=x2 + mlc_padding,
            top_position=y2 + mlc_padding,
            bottom_position=y1 - mlc_padding,
            outer_strip_width=outside_strip_width,
            x_outfield_position=-self.x_width / 2 + outside_strip_width,
        )
        beam = Beam(
            beam_name=beam_name,
            beam_type=BeamType.STATIC,
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
            meter_sets=[0, 1],
            fluence_mode=fluence_mode,
            mlc_boundaries=self.leaf_config,
            machine_name=self.machine_name,
        )
        self.add_beam(beam, mu=beam_mu)

    def to_file(self, filename: str | Path) -> None:
        """Write the DICOM dataset to file"""
        self.ds.save_as(filename, write_like_original=False)

    def as_dicom(self) -> Dataset:
        """Return the new DICOM dataset."""
        return self.ds

    def _add_referenced_beam(self, mu: float) -> None:
        # Referenced Beam Sequence: Referenced Beam 1
        refd_beam1 = Dataset()
        refd_beam1.BeamDose = "1.0"
        refd_beam1.BeamMeterset = str(mu)
        self.ds.FractionGroupSequence[0].ReferencedBeamSequence.append(refd_beam1)
        refd_beam1.ReferencedBeamNumber = (
            len(self.ds.FractionGroupSequence[0].ReferencedBeamSequence) - 1
        )
