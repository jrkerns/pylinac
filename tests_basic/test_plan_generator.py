import tempfile
from unittest import TestCase

import numpy as np
import pydicom
from matplotlib.figure import Figure

from pylinac.plan_generator.dicom import (
    Beam,
    BeamType,
    FluenceMode,
    GantryDirection,
    PlanGenerator,
)
from pylinac.plan_generator.mlc import (
    MLCShaper,
    interpolate_control_points,
    next_sacrifice_shift,
    split_sacrifice_travel,
)
from tests_basic.utils import get_file_from_cloud_test_repo

RT_PLAN_FILE = get_file_from_cloud_test_repo(["plan_generator", "Murray-plan.dcm"])
RT_PLAN_DS = pydicom.dcmread(RT_PLAN_FILE)


class TestPlanGenerator(TestCase):
    def test_from_rt_plan_file(self):
        # shouldn't raise; happy path
        PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )

    def test_from_not_rt_plan_file(self):
        file = get_file_from_cloud_test_repo(["picket_fence", "AS500#2.dcm"])
        with self.assertRaises(ValueError):
            PlanGenerator.from_rt_plan_file(
                file, plan_label="label", plan_name="my name"
            )

    def test_to_file(self):
        pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )
        pg.add_mlc_speed_beams()
        with tempfile.NamedTemporaryFile(delete=False) as t:
            pg.to_file(t.name)
        # shouldn't raise; should be valid DICOM
        ds = pydicom.dcmread(t.name)
        self.assertEqual(ds.RTPlanLabel, "label")
        self.assertEqual(len(ds.BeamSequence), 2)

    def test_from_rt_plan_dataset(self):
        # happy path using a dicom dataset
        dataset = pydicom.dcmread(RT_PLAN_FILE)
        PlanGenerator(dataset, plan_label="label", plan_name="my name")

    def test_no_patient_id(self):
        ds = pydicom.dcmread(RT_PLAN_FILE)
        ds.pop("PatientID")
        with self.assertRaises(ValueError):
            PlanGenerator(ds, plan_label="label", plan_name="my name")

    def test_no_patient_name(self):
        ds = pydicom.dcmread(RT_PLAN_FILE)
        ds.pop("PatientName")
        with self.assertRaises(ValueError):
            PlanGenerator(ds, plan_label="label", plan_name="my name")

    def test_no_tolerance_table(self):
        ds = pydicom.dcmread(RT_PLAN_FILE)
        ds.pop("ToleranceTableSequence")
        with self.assertRaises(ValueError):
            PlanGenerator(ds, plan_label="label", plan_name="my name")

    def test_no_beam_sequence(self):
        ds = pydicom.dcmread(RT_PLAN_FILE)
        ds.pop("BeamSequence")
        with self.assertRaises(ValueError):
            PlanGenerator(ds, plan_label="label", plan_name="my name")

    def test_no_mlc_data(self):
        ds = pydicom.dcmread(RT_PLAN_FILE)
        # pop MLC part of the data; at this point it's just an open field
        ds.BeamSequence[0].BeamLimitingDeviceSequence.pop()
        with self.assertRaises(ValueError):
            PlanGenerator(ds, plan_label="label", plan_name="my name")

    def test_num_leaves(self):
        pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )
        self.assertEqual(pg.num_leaves, 120)

    def test_machine_name(self):
        pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )
        self.assertEqual(pg.machine_name, "TrueBeamSN5837")

    def test_leaf_config(self):
        pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )
        self.assertEqual(len(pg.leaf_config), 61)
        self.assertEqual(max(pg.leaf_config), 200)
        self.assertEqual(min(pg.leaf_config), -200)

    def test_instance_uid_changes(self):
        dcm = pydicom.dcmread(RT_PLAN_FILE)
        pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE, plan_label="label", plan_name="my name"
        )
        pg_dcm = pg.as_dicom()
        self.assertNotEqual(pg_dcm.SOPInstanceUID, dcm.SOPInstanceUID)


def create_beam(**kwargs) -> Beam:
    return Beam(
        plan_dataset=kwargs.get("plan_dataset", RT_PLAN_DS),
        beam_name=kwargs.get("beam_name", "name"),
        beam_type=kwargs.get("beam_type", BeamType.DYNAMIC),
        energy=kwargs.get("energy", 6),
        dose_rate=kwargs.get("dose_rate", 600),
        x1=kwargs.get("x1", -5),
        x2=kwargs.get("x2", 5),
        y1=kwargs.get("y1", -5),
        y2=kwargs.get("y2", 5),
        machine_name=kwargs.get("machine_name", "TrueBeam"),
        gantry_angles=kwargs.get("gantry_angles", 0),
        gantry_direction=kwargs.get("gantry_direction", GantryDirection.NONE),
        coll_angle=kwargs.get("coll_angle", 0),
        couch_vrt=kwargs.get("couch_vrt", 0),
        couch_lng=kwargs.get("couch_lng", 0),
        couch_lat=kwargs.get("couch_lat", 0),
        couch_rot=kwargs.get("couch_rot", 0),
        mlc_boundaries=kwargs.get("mlc_boundaries", [-200, -100, 0, 100, 200]),
        mlc_positions=kwargs.get("mlc_positions", [[0], [0]]),
        metersets=kwargs.get("metersets", [0, 1]),
        fluence_mode=kwargs.get("fluence_mode", FluenceMode.STANDARD),
    )


class TestBeam(TestCase):
    def test_beam_normal(self):
        # shouldn't raise; happy path
        beam = create_beam(
            gantry_angles=0,
        )
        beam_dcm = beam.as_dicom()
        self.assertEqual(beam_dcm.BeamName, "name")
        self.assertEqual(beam_dcm.BeamType, "DYNAMIC")
        self.assertEqual(beam_dcm.ControlPointSequence[0].GantryAngle, 0)

    def test_too_long_beam_name(self):
        with self.assertRaises(ValueError):
            create_beam(beam_name="superlongbeamname")

    def test_1_mlc_position_for_static(self):
        with self.assertRaises(ValueError):
            create_beam(
                mlc_positions=[[0], [0], [0], [0], [0]],
                beam_type=BeamType.STATIC,
                metersets=[0, 0, 0, 0, 0],
            )

        # valid
        beam = create_beam(
            mlc_positions=[[0]], metersets=[0], beam_type=BeamType.STATIC
        )
        self.assertEqual(beam.as_dicom().BeamType, "STATIC")

    def test_mlc_positions_match_metersets(self):
        with self.assertRaises(ValueError):
            create_beam(mlc_positions=[[0] * 5], metersets=[1, 2, 3])

    def test_gantry_must_be_dynamic_for_multiple_angles(self):
        with self.assertRaises(ValueError) as context:
            create_beam(
                gantry_angles=[0, 90],
                mlc_positions=[
                    [0],
                ],
                metersets=[0],
                beam_type=BeamType.STATIC,
            )
        self.assertIn("Cannot specify multiple gantry angles", str(context.exception))

        # valid
        beam = create_beam(
            gantry_angles=[0, 90],
            beam_type=BeamType.DYNAMIC,
            gantry_direction=GantryDirection.CLOCKWISE,
        )
        beam_dcm = beam.as_dicom()
        self.assertEqual(beam_dcm.BeamType, "DYNAMIC")
        self.assertEqual(beam_dcm.ControlPointSequence[0].GantryAngle, 0)
        self.assertEqual(beam_dcm.ControlPointSequence[0].GantryRotationDirection, "CC")
        self.assertEqual(beam_dcm.ControlPointSequence[1].GantryAngle, 90)

    def test_must_have_gantry_direction_if_dynamic(self):
        with self.assertRaises(ValueError):
            create_beam(
                gantry_angles=[0, 90],
                beam_type=BeamType.DYNAMIC,
                gantry_direction=GantryDirection.NONE,
            )

        # valid
        create_beam(
            gantry_angles=[0, 90],
            beam_type=BeamType.DYNAMIC,
            gantry_direction=GantryDirection.CLOCKWISE,
        )

    def test_jaw_positions(self):
        b = create_beam(x1=-5, x2=7, y1=-11, y2=13)
        dcm = b.as_dicom()
        self.assertEqual(
            len(dcm.ControlPointSequence[0].BeamLimitingDevicePositionSequence), 3
        )
        self.assertEqual(
            dcm.ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[0]
            .LeafJawPositions,
            [-5, 7],
        )
        self.assertEqual(
            dcm.ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions,
            [-11, 13],
        )

    def test_multiple_gantry_angles_static_not_allowed(self):
        with self.assertRaises(ValueError):
            create_beam(gantry_angles=[0, 90], beam_type=BeamType.STATIC)


class TestPlanGeneratorBeams(TestCase):
    """Test real workflow where beams are added"""

    def setUp(self) -> None:
        self.pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE,
            plan_label="label",
            plan_name="my name",
        )

    def test_add_beam_low_level(self):
        self.pg.add_beam(
            create_beam(plan_dataset=self.pg.as_dicom()).as_dicom(), mu=100
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 1)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "name")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 1)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 100
        )
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].ReferencedBeamNumber,
            0,
        )

    def test_add_2_beams(self):
        self.pg.add_beam(create_beam().as_dicom(), mu=100)
        self.pg.add_beam(create_beam(beam_name="beam2").as_dicom(), mu=200)
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 2)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 2)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 100
        )
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[1].BeamMeterset, 200
        )
        self.assertEqual(dcm.BeamSequence[1].BeamName, "beam2")
        self.assertEqual(dcm.BeamSequence[1].BeamNumber, 1)

    def test_plot_fluence(self):
        # just tests it works
        self.pg.add_open_field_beam(x1=-5, x2=5, y1=-5, y2=5, mu=100)
        figs = self.pg.plot_fluences()
        self.assertIsInstance(figs, list)
        self.assertIsInstance(figs[0], Figure)


class TestPlanPrefabs(TestCase):
    def setUp(self) -> None:
        self.pg = PlanGenerator.from_rt_plan_file(
            RT_PLAN_FILE,
            plan_label="label",
            plan_name="my name",
        )

    def test_create_open_field(self):
        self.pg.add_open_field_beam(
            x1=-100,
            x2=100,
            y1=-110,
            y2=110,
            mu=123,
            beam_name="Open Field",
            defined_by_mlcs=True,
            padding_mm=0,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 1)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "Open Field")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 1)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[0]
            .LeafJawPositions,
            [-100, 100],
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions,
            [-110, 110],
        )

    def test_open_field_jaws(self):
        self.pg.add_open_field_beam(
            x1=-100,
            x2=100,
            y1=-110,
            y2=110,
            mu=123,
            beam_name="Open Field",
            defined_by_mlcs=False,
            padding_mm=0,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 1)
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[0]
            .LeafJawPositions,
            [-100, 100],
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions,
            [-110, 110],
        )

    def test_create_picket_fence(self):
        self.pg.add_picketfence_beam(
            y1=-10,
            y2=10,
            mu=123,
            beam_name="Picket Fence",
            strip_positions_mm=(-50, -30, -10, 10, 30, 50),
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 1)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "Picket Fence")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 1)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[0]
            .LeafJawPositions,
            [-60, 60],
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions,
            [-10, 10],
        )

    def test_picket_fence_too_wide(self):
        with self.assertRaises(ValueError):
            self.pg.add_picketfence_beam(
                y1=-10,
                y2=10,
                mu=123,
                beam_name="Picket Fence",
                strip_positions_mm=(-100, 100),
            )

    def test_winston_lutz_beams(self):
        self.pg.add_winston_lutz_beams(
            axes_positions=(
                {"gantry": 0, "collimator": 0, "couch": 0},
                {"gantry": 90, "collimator": 0, "couch": 0},
                {"gantry": 180, "collimator": 0, "couch": 45},
            ),
            x1=-10,
            x2=10,
            y1=-10,
            y2=10,
            mu=123,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 3)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "G0C0P0")
        self.assertEqual(dcm.BeamSequence[2].BeamName, "G180C0P45")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 3)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )
        self.assertEqual(dcm.BeamSequence[0].ControlPointSequence[0].GantryAngle, 0)
        self.assertEqual(dcm.BeamSequence[1].ControlPointSequence[0].GantryAngle, 90)
        self.assertEqual(dcm.BeamSequence[2].ControlPointSequence[0].GantryAngle, 180)

    def test_winston_lutz_jaw_defined(self):
        self.pg.add_winston_lutz_beams(
            axes_positions=(
                {"gantry": 0, "collimator": 0, "couch": 0},
                {"gantry": 90, "collimator": 0, "couch": 0},
                {"gantry": 180, "collimator": 0, "couch": 45},
            ),
            x1=-10,
            x2=10,
            y1=-10,
            y2=10,
            mu=123,
            defined_by_mlcs=False,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 3)
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[0]
            .LeafJawPositions,
            [-10, 10],
        )
        self.assertEqual(
            dcm.BeamSequence[0]
            .ControlPointSequence[0]
            .BeamLimitingDevicePositionSequence[1]
            .LeafJawPositions,
            [-10, 10],
        )

    def test_dose_rate_beams(self):
        self.pg.add_dose_rate_beams(
            dose_rates=(100, 400, 600),
            y1=-10,
            y2=10,
            desired_mu=123,
            default_dose_rate=600,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 2)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "DR Ref")
        self.assertEqual(dcm.BeamSequence[1].BeamName, "DR100-600")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 2)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )

    def test_dose_rate_too_wide(self):
        with self.assertRaises(ValueError):
            self.pg.add_dose_rate_beams(
                dose_rates=(100, 150, 200, 250, 300, 350, 400, 600),
                roi_size_mm=30,
                y1=-10,
                y2=10,
                desired_mu=123,
                default_dose_rate=600,
            )

    def test_mlc_speed_beams(self):
        self.pg.add_mlc_speed_beams(
            speeds=(0.5, 1, 1.5, 2),
            y1=-100,
            y2=100,
            mu=123,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 2)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "MLC Sp Ref")
        self.assertEqual(dcm.BeamSequence[1].BeamName, "MLC Speed")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 2)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )

    def test_mlc_speed_too_fast(self):
        with self.assertRaises(ValueError):
            self.pg.add_mlc_speed_beams(
                speeds=(10, 20, 30, 40, 50),
                y1=-100,
                y2=100,
            )

    def test_mlc_speed_too_wide(self):
        with self.assertRaises(ValueError):
            self.pg.add_mlc_speed_beams(
                speeds=(0.5, 1, 1.5, 2),
                roi_size_mm=50,
                y1=-100,
                y2=100,
            )

    def test_0_mlc_speed(self):
        with self.assertRaises(ValueError):
            self.pg.add_mlc_speed_beams(
                speeds=(0, 1, 2),
                y1=-100,
                y2=100,
            )

    def test_gantry_speed_beams(self):
        # max speed is 2.5 by default
        self.pg.add_gantry_speed_beams(
            speeds=(1, 2, 3, 4),
            y1=-100,
            y2=100,
            mu=123,
        )
        dcm = self.pg.as_dicom()
        self.assertEqual(len(dcm.BeamSequence), 2)
        self.assertEqual(dcm.BeamSequence[0].BeamName, "GS")
        self.assertEqual(dcm.BeamSequence[1].BeamName, "G Sp Ref")
        self.assertEqual(dcm.BeamSequence[0].BeamNumber, 0)
        self.assertEqual(dcm.FractionGroupSequence[0].NumberOfBeams, 2)
        self.assertEqual(
            dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset, 123
        )

    def test_gantry_speed_too_fast(self):
        # max speed is 4.8 by default
        with self.assertRaises(ValueError):
            self.pg.add_gantry_speed_beams(
                speeds=(1, 2, 3, 4, 5),
                y1=-100,
                y2=100,
            )

    def test_gantry_speed_too_wide(self):
        with self.assertRaises(ValueError):
            self.pg.add_gantry_speed_beams(
                speeds=(1, 2, 3, 4),
                roi_size_mm=50,
                y1=-100,
                y2=100,
            )

    def test_gantry_range_over_360(self):
        with self.assertRaises(ValueError):
            self.pg.add_gantry_speed_beams(
                speeds=(4, 4, 4, 4),
                y1=-100,
                y2=100,
                mu=250,
            )


class TestMLCShaper(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # simplistic MLC setup
        cls.leaf_boundaries: list[float] = np.arange(
            start=-200, stop=201, step=5
        ).tolist()

    def test_init(self):
        MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )

    def test_num_leaves(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        self.assertEqual(shaper.num_leaves, 160)

    def test_meterset_over_1(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        with self.assertRaises(ValueError):
            shaper.add_strip(position_mm=-5, strip_width_mm=0, meterset_at_target=2)

    def test_sacrifice_without_transition_dose(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        with self.assertRaises(ValueError):
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=1,
                meterset_transition=0,
                sacrificial_distance_mm=50,
            )

    def test_initial_sacrificial_gap(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        shaper.add_strip(
            position_mm=-5,
            strip_width_mm=0,
            meterset_at_target=1,
            initial_sacrificial_gap_mm=10,
        )
        self.assertEqual(shaper.control_points[0][0], -10)

    def test_cant_add_sacrificial_gap_after_first_point(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        shaper.add_strip(
            position_mm=-5,
            strip_width_mm=0,
            meterset_at_target=0.2,
            initial_sacrificial_gap_mm=5,
        )
        with self.assertRaises(ValueError) as context:
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=0.2,
                initial_sacrificial_gap_mm=10,
            )
        self.assertIn("already control points", str(context.exception))

    def test_cant_have_initial_sacrifice_and_transition_dose(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        with self.assertRaises(ValueError):
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=1,
                initial_sacrificial_gap_mm=5,
            )

    def test_cant_have_meterset_transition_for_first_control_point(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        with self.assertRaises(ValueError) as context:
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=0,
                meterset_transition=1,
            )
        self.assertIn("Cannot have a transition", str(context.exception))

    def test_cant_have_initial_sacrificial_gap_and_sacrificial_distance(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        with self.assertRaises(ValueError) as context:
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=0.5,
                meterset_transition=0.1,
                sacrificial_distance_mm=5,
                initial_sacrificial_gap_mm=5,
            )
        self.assertIn("Cannot specify both", str(context.exception))

    def test_cannot_have_sacrifical_gap_on_secondary_control_point(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        shaper.add_strip(position_mm=-5, strip_width_mm=0, meterset_at_target=0.5)
        with self.assertRaises(ValueError) as context:
            shaper.add_strip(
                position_mm=-5,
                strip_width_mm=0,
                meterset_at_target=0.5,
                initial_sacrificial_gap_mm=10,
            )
        self.assertIn("already control points", str(context.exception))

    def test_split_sacrifices(self):
        res = split_sacrifice_travel(distance=33, max_travel=20)
        self.assertCountEqual(res, [20, 13])
        res = split_sacrifice_travel(distance=11, max_travel=20)
        self.assertCountEqual(res, [11])
        res = split_sacrifice_travel(distance=66, max_travel=20)
        self.assertCountEqual(res, [20, 20, 20, 6])

    def test_as_control_points(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        shaper.add_strip(position_mm=-5, strip_width_mm=0, meterset_at_target=1)
        cp = shaper.as_control_points()
        self.assertEqual(
            len(cp), 2
        )  # start and end positions given meterset increments
        self.assertEqual(cp[0][0], -5)

    def test_as_metersets(self):
        shaper = MLCShaper(
            leaf_y_positions=self.leaf_boundaries, max_x_mm=400, sacrifice_gap_mm=5
        )
        shaper.add_strip(position_mm=-5, strip_width_mm=0, meterset_at_target=1)
        metersets = shaper.as_metersets()
        self.assertEqual(metersets, [0, 1])


class TestNextSacrificeShift(TestCase):
    def test_easy(self):
        target = next_sacrifice_shift(
            current_position_mm=0,
            travel_mm=5,
            x_width_mm=400,
            other_mlc_position=0,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, -5)

    def test_toward_target_right(self):
        target = next_sacrifice_shift(
            current_position_mm=-5,
            travel_mm=50,
            x_width_mm=400,
            other_mlc_position=0,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, 50)

    def test_toward_target_left(self):
        target = next_sacrifice_shift(
            current_position_mm=45,
            travel_mm=50,
            x_width_mm=400,
            other_mlc_position=0,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, -50)

    def test_travel_too_large(self):
        with self.assertRaises(ValueError):
            next_sacrifice_shift(
                current_position_mm=0,
                travel_mm=200,
                x_width_mm=400,
                other_mlc_position=0,
                max_overtravel_mm=140,
            )

    def test_travel_can_be_over_max_overtravel_if_on_other_side(self):
        target = next_sacrifice_shift(
            current_position_mm=0,
            travel_mm=200,
            x_width_mm=400,
            other_mlc_position=100,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, 200)

    def test_at_edge_of_width(self):
        target = next_sacrifice_shift(
            current_position_mm=-180,
            travel_mm=30,
            x_width_mm=400,
            other_mlc_position=-190,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, 30)

        target = next_sacrifice_shift(
            current_position_mm=180,
            travel_mm=30,
            x_width_mm=400,
            other_mlc_position=190,
            max_overtravel_mm=140,
        )
        self.assertEqual(target, -30)

    def test_width_vs_overtravel(self):
        with self.assertRaises(ValueError):
            next_sacrifice_shift(
                current_position_mm=0,
                travel_mm=30,
                x_width_mm=100,
                other_mlc_position=-190,
                max_overtravel_mm=140,
            )


class TestInterpolateControlPoints(TestCase):
    """For these tests, we use a simplified version of a 3-pair MLC. The first and last pair are the sacrificial leaves."""

    def test_control_point_lengths_mismatch(self):
        with self.assertRaises(ValueError):
            interpolate_control_points(
                control_point_start=[0, 0, 0, 0, 0],
                control_point_end=[10, 10, 10, 10],
                interpolation_ratios=[0.5],
                sacrifice_chunks=[5],
                max_overtravel=140,
            )

    def test_no_interpolation(self):
        with self.assertRaises(ValueError):
            interpolate_control_points(
                control_point_start=[0, 0, 0, 0, 0],
                control_point_end=[10, 10, 10, 10, 10],
                interpolation_ratios=[],
                sacrifice_chunks=[5],
                max_overtravel=140,
            )

    def test_interpolate_simple(self):
        interp_cp = interpolate_control_points(
            control_point_start=[0, 0, 0, 0, 0, 0],
            control_point_end=[10, 10, 10, 10, 10, 10],
            interpolation_ratios=[0.5],
            sacrifice_chunks=[1],
            max_overtravel=140,
        )
        # the first, middle, and last values should be the sacrifice
        # the middle values should be the interpolated values
        self.assertEqual(interp_cp, [[-1, 5, -1, -1, 5, -1]])

    def test_interpolate_multiple(self):
        interp_cp = interpolate_control_points(
            control_point_start=[0, 0, 0, 0, 0, 0],
            control_point_end=[10, 10, 10, 10, 10, 10],
            interpolation_ratios=[0.25, 0.5, 0.75],
            sacrifice_chunks=[3, 5, 7],
            max_overtravel=140,
        )
        # the sacrifice goes 0 - 3 -> -3 + 5 -> 2 + 7 -> 9
        cp1 = [-3, 2.5, -3, -3, 2.5, -3]
        self.assertEqual(interp_cp[0], cp1)
        cp2 = [2, 5, 2, 2, 5, 2]
        self.assertEqual(interp_cp[1], cp2)
        cp3 = [9, 7.5, 9, 9, 7.5, 9]
        self.assertEqual(interp_cp[2], cp3)

    def test_overtravel(self):
        # 30 is over the max overtravel of 20
        with self.assertRaises(ValueError):
            interpolate_control_points(
                control_point_start=[0, 0, 0, 0, 0, 0],
                control_point_end=[10, 10, 10, 10, 10, 10],
                interpolation_ratios=[0.5],
                sacrifice_chunks=[30],
                max_overtravel=20,
            )

    def test_interpolation_over_1_or_0(self):
        with self.assertRaises(ValueError):
            interpolate_control_points(
                control_point_start=[0, 0, 0, 0, 0, 0],
                control_point_end=[10, 10, 10, 10, 10, 10],
                interpolation_ratios=[1.5],
                sacrifice_chunks=[5],
                max_overtravel=140,
            )
        with self.assertRaises(ValueError):
            interpolate_control_points(
                control_point_start=[0, 0, 0, 0, 0, 0],
                control_point_end=[10, 10, 10, 10, 10, 10],
                interpolation_ratios=[-0.5],
                sacrifice_chunks=[5],
                max_overtravel=140,
            )
