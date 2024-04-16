from unittest import TestCase

import pydicom

from pylinac.plan_generator.dicom import (
    Beam,
    BeamType,
    FluenceMode,
    GantryDirection,
    PlanGenerator,
)
from tests_basic.utils import get_file_from_cloud_test_repo


class TestPlanGenerator(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rt_plan_file = get_file_from_cloud_test_repo(
            ["plan_generator", "Murray-plan.dcm"]
        )

    def test_from_rt_plan_file(self):
        # shouldn't raise; happy path
        PlanGenerator.from_rt_plan_file(
            self.rt_plan_file, plan_label="label", plan_name="my name"
        )

    def test_from_not_rt_plan_file(self):
        file = get_file_from_cloud_test_repo(["picket_fence", "AS500#2.dcm"])
        with self.assertRaises(ValueError):
            PlanGenerator.from_rt_plan_file(
                file, plan_label="label", plan_name="my name"
            )

    def test_from_rt_plan_dataset(self):
        # happy path using a dicom dataset
        dataset = pydicom.dcmread(self.rt_plan_file)
        PlanGenerator(dataset, plan_label="label", plan_name="my name")

    def test_num_leaves(self):
        pg = PlanGenerator.from_rt_plan_file(
            self.rt_plan_file, plan_label="label", plan_name="my name"
        )
        self.assertEqual(pg.num_leaves, 120)

    def test_machine_name(self):
        pg = PlanGenerator.from_rt_plan_file(
            self.rt_plan_file, plan_label="label", plan_name="my name"
        )
        self.assertEqual(pg.machine_name, "TrueBeamSN5837")

    def test_leaf_config(self):
        pg = PlanGenerator.from_rt_plan_file(
            self.rt_plan_file, plan_label="label", plan_name="my name"
        )
        self.assertEqual(len(pg.leaf_config), 61)
        self.assertEqual(max(pg.leaf_config), 200)
        self.assertEqual(min(pg.leaf_config), -200)

    def test_instance_uid_changes(self):
        dcm = pydicom.dcmread(self.rt_plan_file)
        pg = PlanGenerator.from_rt_plan_file(
            self.rt_plan_file, plan_label="label", plan_name="my name"
        )
        pg_dcm = pg.as_dicom()
        self.assertNotEqual(pg_dcm.SOPInstanceUID, dcm.SOPInstanceUID)


def create_beam(**kwargs) -> Beam:
    return Beam(
        beam_name=kwargs.get("beam_name", "name"),
        beam_type=kwargs.get("beam_type", BeamType.STATIC),
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
        meter_sets=kwargs.get("meter_sets", [0, 1]),
        fluence_mode=kwargs.get("fluence_mode", FluenceMode.STANDARD),
    )


class TestBeam(TestCase):
    def test_beam_normal(self):
        # shouldn't raise; happy path
        create_beam()

    def test_too_long_beam_name(self):
        with self.assertRaises(ValueError):
            create_beam(beam_name="superlongbeamname")

    def test_1_mlc_position_for_static(self):
        with self.assertRaises(ValueError):
            create_beam(
                mlc_positions=[[0], [0], [0], [0], [0]], beam_type=BeamType.STATIC
            )

        # valid
        create_beam(mlc_positions=[[0]], beam_type=BeamType.STATIC)

    def test_mlc_positions_match_meter_sets(self):
        with self.assertRaises(ValueError):
            create_beam(mlc_positions=[[0] * 5], meter_sets=[1, 2, 3])

    def test_gantry_must_be_dynamic_for_multiple_angles(self):
        with self.assertRaises(ValueError):
            create_beam(gantry_angles=[0, 90], beam_type=BeamType.STATIC)

        # valid
        create_beam(gantry_angles=[0, 90], beam_type=BeamType.DYNAMIC)

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
