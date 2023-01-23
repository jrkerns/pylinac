import copy
import io
import math
import tempfile
from unittest import TestCase

import matplotlib.pyplot as plt

import pylinac
from pylinac import WinstonLutz, WinstonLutzMultiTarget
from pylinac.core.geometry import Vector, vector_is_close, Line, Point
from pylinac.core.io import TemporaryZipDirectory
from pylinac.core.scale import MachineScale
from pylinac.winston_lutz import Axis, WinstonLutzResult, WinstonLutz2D, bb_projection_long, bb_projection_gantry_plane, \
    BB, bb_ray_line, BBArrangement
from tests_basic.utils import (
    save_file,
    CloudFileMixin,
    get_folder_from_cloud_test_repo,
    get_file_from_cloud_test_repo,
    FromDemoImageTesterMixin,
    FromURLTesterMixin,
)

TEST_DIR = "Winston-Lutz"


class TestProjection(TestCase):
    """Test the BB isoplane projections"""

    def test_longitudinal_projection(self):
        # in coordinate space, positive is in, but in plotting space, positive is out
        # thus, we return the opposite sign than the coordinate space
        # dead center
        assert bb_projection_long(offset_in=0, offset_up=0, offset_left=0, sad=1000, gantry=0) == 0
        # up-only won't change it
        assert bb_projection_long(offset_in=0, offset_up=30, offset_left=0, sad=1000, gantry=0) == 0
        # long-only won't change it
        assert bb_projection_long(offset_in=20, offset_up=0, offset_left=0, sad=1000, gantry=0) == 20
        # left-only won't change it
        assert bb_projection_long(offset_in=0, offset_up=0, offset_left=15, sad=1000, gantry=0) == 0
        # in and up will make it look further away at gantry 0
        assert math.isclose(bb_projection_long(offset_in=10, offset_up=10, offset_left=0, sad=1000, gantry=0), 10.1, abs_tol=0.005)
        # in and down will make it closer at gantry 0
        assert math.isclose(bb_projection_long(offset_in=10, offset_up=-10, offset_left=0, sad=1000, gantry=0), 9.9, abs_tol=0.005)
        # in and up will make it look closer at gantry 180
        assert math.isclose(bb_projection_long(offset_in=10, offset_up=10, offset_left=0, sad=1000, gantry=180), 9.9, abs_tol=0.005)
        # in and down will make it further away at gantry 180
        assert math.isclose(bb_projection_long(offset_in=10, offset_up=-10, offset_left=0, sad=1000, gantry=180), 10.1, abs_tol=0.005)
        # in and left will make it closer at gantry 90
        assert math.isclose(
            bb_projection_long(offset_in=10, offset_up=0, offset_left=10, sad=1000, gantry=90), 9.9,
            abs_tol=0.005)
        # in and right will make it further away at gantry 90
        assert math.isclose(
            bb_projection_long(offset_in=10, offset_up=0, offset_left=-10, sad=1000, gantry=90), 10.1,
            abs_tol=0.005)
        # in and right will make it closer at gantry 270
        assert math.isclose(
            bb_projection_long(offset_in=10, offset_up=0, offset_left=-10, sad=1000, gantry=270), 9.9,
            abs_tol=0.005)
        # in and left won't change at gantry 0
        assert math.isclose(
            bb_projection_long(offset_in=10, offset_up=0, offset_left=10, sad=1000, gantry=0), 10,
            abs_tol=0.005)
        # double the sad will half the effect:
        # in and up will make it look further away at gantry 0
        assert math.isclose(bb_projection_long(offset_in=10, offset_up=10, offset_left=0, sad=1000, gantry=0), 10.1, abs_tol=0.005)
        # out and up will make it look further away at gantry 0
        assert math.isclose(bb_projection_long(offset_in=-10, offset_up=10, offset_left=0, sad=1000, gantry=0), -10.1, abs_tol=0.005)
        # out and up will make it look closer at gantry 180
        assert math.isclose(bb_projection_long(offset_in=-10, offset_up=10, offset_left=0, sad=1000, gantry=180), -9.9, abs_tol=0.005)
        # out and down will make it look closer at gantry 0
        assert math.isclose(bb_projection_long(offset_in=-10, offset_up=-10, offset_left=0, sad=1000, gantry=0), -9.9, abs_tol=0.005)
        # out and down will make it look further out at gantry 180
        assert math.isclose(bb_projection_long(offset_in=-10, offset_up=-10, offset_left=0, sad=1000, gantry=180), -10.1, abs_tol=0.005)

    def test_gantry_plane_projection(self):
        # left is negative, right is positive
        # dead center
        assert bb_projection_gantry_plane(offset_up=0, offset_left=0, sad=1000, gantry=0) == 0
        # up-only at gantry 0 is still 0
        assert bb_projection_gantry_plane(offset_up=10, offset_left=0, sad=1000, gantry=0) == 0
        # up-only at gantry 90 is exactly negative the offset
        assert bb_projection_gantry_plane(offset_up=10, offset_left=0, sad=1000, gantry=90) == -10
        # down-only at gantry 90 is exactly the offset
        assert bb_projection_gantry_plane(offset_up=-10, offset_left=0, sad=1000, gantry=90) == 10
        # left-only at gantry 0 is exactly negative the offset
        assert bb_projection_gantry_plane(offset_up=0, offset_left=10, sad=1000, gantry=0) == -10
        # right-only at gantry 0 is exactly negative the offset
        assert bb_projection_gantry_plane(offset_up=0, offset_left=-10, sad=1000, gantry=0) == 10
        # left-only at gantry 180 is exactly the offset
        assert bb_projection_gantry_plane(offset_up=0, offset_left=10, sad=1000, gantry=180) == 10
        # left and up at gantry 0 makes the bb appear away from CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=10, offset_left=20, sad=1000, gantry=0), -20.2, abs_tol=0.005)
        # left and down at gantry 0 makes the bb appear closer to the CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=-10, offset_left=20, sad=1000, gantry=0), -19.8, abs_tol=0.005)
        # left and up at gantry 180 makes the bb appear closer to CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=10, offset_left=20, sad=1000, gantry=180), 19.8, abs_tol=0.005)
        # left and up at gantry 90 makes the bb appear closer to CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=10, offset_left=20, sad=1000, gantry=90), -9.8, abs_tol=0.005)
        # left and down at gantry 90 makes the bb appear closer to CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=-10, offset_left=20, sad=1000, gantry=90), 9.8, abs_tol=0.005)
        # left and down at gantry 270 makes the bb appear further from the CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=-10, offset_left=20, sad=1000, gantry=270), -10.2, abs_tol=0.005)
        # right and down at gantry 270 makes the bb appear closer the CAX
        assert math.isclose(bb_projection_gantry_plane(offset_up=-10, offset_left=-20, sad=1000, gantry=270), -9.8, abs_tol=0.005)


class TestBB(TestCase):

    @property
    def nominal_0(self):
        return {'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0, 'bb_size_mm': 5}

    def test_measured_0_at_0(self):
        """The determined BB position when the rays are about the origin should be the origin"""
        line_g0 = Line(point1=Point(x=1000, y=0, z=0), point2=Point(x=-1000, y=0, z=0))
        line_g90 = Line(Point(x=0, y=0, z=1000), Point(x=0, y=0, z=-1000))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert bb.measured_position.x == bb.measured_position.y == bb.measured_position.z == 0

    def test_measured_in_1mm(self):
        """The determined BB position when the rays are shifted in 1mm should be 1mm"""
        line_g0 = Line(point1=Point(x=1000, y=1, z=0), point2=Point(x=-1000, y=1, z=0))
        line_g90 = Line(Point(x=0, y=1, z=1000), Point(x=0, y=1, z=-1000))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, 1, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, 0, abs_tol=0.001)

    def test_measured_out_1mm(self):
        """The determined BB position when the rays are shifted out"""
        line_g0 = Line(point1=Point(x=1000, y=-1, z=0), point2=Point(x=-1000, y=-1, z=0))
        line_g90 = Line(Point(x=0, y=-1, z=1000), Point(x=0, y=-1, z=-1000))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, -1, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, 0, abs_tol=0.001)

    def test_measured_up_1mm(self):
        """The determined BB position when the rays are shifted in 1mm should be 1mm"""
        line_g0 = Line(point1=Point(x=0, y=0, z=1000), point2=Point(x=0, y=0, z=-1000))
        line_g90 = Line(Point(x=-1000, y=0, z=1), Point(x=1000, y=0, z=1))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, 1, abs_tol=0.001)

    def test_measured_down_1mm(self):
        """The determined BB position when the rays are shifted in 1mm should be 1mm"""
        line_g0 = Line(point1=Point(x=0, y=0, z=1000), point2=Point(x=0, y=0, z=-1000))
        line_g90 = Line(Point(x=-1000, y=0, z=-1), Point(x=1000, y=0, z=-1))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, -1, abs_tol=0.001)

    def test_measured_left_1mm(self):
        """The determined BB position when the rays are shifted in 1mm should be 1mm"""
        line_g0 = Line(point1=Point(x=-1, y=0, z=1000), point2=Point(x=-1, y=0, z=-1000))
        line_g90 = Line(Point(x=-1000, y=0, z=0), Point(x=1000, y=0, z=0))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, -1, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, 0, abs_tol=0.001)

    def test_measured_right_1mm(self):
        """The determined BB position when the rays are shifted in 1mm should be 1mm"""
        line_g0 = Line(point1=Point(x=1, y=0, z=1000), point2=Point(x=1, y=0, z=-1000))
        line_g90 = Line(Point(x=-1000, y=0, z=0), Point(x=1000, y=0, z=0))
        bb = BB(self.nominal_0, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.measured_position.x, 1, abs_tol=0.001)
        assert math.isclose(bb.measured_position.y, 0, abs_tol=0.001)
        assert math.isclose(bb.measured_position.z, 0, abs_tol=0.001)

    def test_delta_0_at_0(self):
        """The delta vector should be 0 when the measured bb and nominal position are the same"""
        line_g0 = Line(point1=Point(x=1000, y=0, z=0), point2=Point(x=-1000, y=0, z=0))
        line_g90 = Line(Point(x=0, y=0, z=1000), Point(x=0, y=0, z=-1000))
        bb = BB(nominal_bb={'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0}, ray_lines=[line_g0, line_g90])
        assert bb.delta_vector.x == 0
        assert bb.delta_vector.y == 0
        assert bb.delta_vector.z == 0

    def test_delta_vector_1_at_shift_1(self):
        """When the nominal BB is at 0 and the measured BB position is 1 it should be 1"""
        line_g0 = Line(point1=Point(x=1000, y=1, z=0), point2=Point(x=-1000, y=1, z=0))
        line_g90 = Line(Point(x=0, y=1, z=1000), Point(x=0, y=1, z=-1000))
        bb = BB(nominal_bb={'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0}, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.delta_vector.x, 0, abs_tol=0.001)
        assert math.isclose(bb.delta_vector.y, 1, abs_tol=0.001)
        assert math.isclose(bb.delta_vector.z, 0, abs_tol=0.001)

    def test_delta_scalar_1_at_shift_1(self):
        """When the nominal BB is at 0 and the measured BB position is 1 it should be 1"""
        line_g0 = Line(point1=Point(x=1000, y=1, z=0), point2=Point(x=-1000, y=1, z=0))
        line_g90 = Line(Point(x=0, y=1, z=1000), Point(x=0, y=1, z=-1000))
        bb = BB(nominal_bb={'offset_left_mm': 0, 'offset_up_mm': 0, 'offset_in_mm': 0}, ray_lines=[line_g0, line_g90])
        assert math.isclose(bb.delta_distance, 1, abs_tol=0.001)


class TestRays(TestCase):
    sad = 1000
    dpmm = 1
    image_center = Point(0, 0)

    def test_g0_centered(self):
        """A BB at iso should produce a perfect line"""
        bb = Point(0, 0)
        gantry = 0
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert line.point1.x == 0
        assert line.point2.x == 0
        assert math.isclose(line.point1.y, 0)
        assert math.isclose(line.point2.y, 0)
        assert math.isclose(line.point1.z, 1000)
        assert math.isclose(line.point2.z, -1000)

    def test_g0_out_1mm(self):
        bb = Point(0, 1)
        gantry = 0
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert line.point1.y == 0
        assert line.point2.y == 2

    def test_g0_in_1mm(self):
        bb = Point(0, -1)
        gantry = 0
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert line.point1.y == 0
        assert line.point2.y == -2

    def test_g0_left_1mm(self):
        bb = Point(-1, 0)
        gantry = 0
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert line.point1.x == 0
        assert line.point2.x == -2

    def test_g0_right_1mm(self):
        bb = Point(1, 0)
        gantry = 0
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert line.point1.x == 0
        assert line.point2.x == 2

    def test_g90_centered(self):
        bb = Point(0, 0)
        gantry = 90
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert math.isclose(line.point1.z, 0, abs_tol=0.005)
        assert math.isclose(line.point2.z, 0, abs_tol=0.005)
        assert math.isclose(line.point1.x, 1000, abs_tol=0.005)
        assert math.isclose(line.point2.x, -1000, abs_tol=0.005)

    def test_g90_up_1mm(self):
        bb = Point(-1, 0)
        gantry = 90
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert math.isclose(line.point1.z, 0, abs_tol=0.005)
        assert math.isclose(line.point2.z, 2, abs_tol=0.005)

    def test_g90_down_1mm(self):
        bb = Point(1, 0)
        gantry = 90
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert math.isclose(line.point1.z, 0, abs_tol=0.005)
        assert math.isclose(line.point2.z, -2, abs_tol=0.005)

    def test_g270_up_1mm(self):
        bb = Point(1, 0)
        gantry = 270
        line = bb_ray_line(bb=bb, gantry_angle=gantry, sad=self.sad, image_center=self.image_center, dpmm=self.dpmm)
        assert math.isclose(line.point1.z, 0, abs_tol=0.005)
        assert math.isclose(line.point2.z, 2, abs_tol=0.005)


class TestWLMultiImage(TestCase):

    def test_demo_images(self):
        wl = WinstonLutzMultiTarget.from_demo_images()
        # shouldn't raise
        wl.analyze(BBArrangement.SNC_MULTIMET)

    def test_demo(self):
        # shouldn't raise
        WinstonLutzMultiTarget.run_demo()

    def test_plot_locations(self):
        wl = WinstonLutzMultiTarget.from_demo_images()
        wl.analyze(BBArrangement.SNC_MULTIMET)
        wl.plot_locations()

    def test_publish_pdf(self):
        wl = WinstonLutzMultiTarget.from_demo_images()
        wl.analyze(BBArrangement.SNC_MULTIMET)
        wl.publish_pdf('output.pdf')


class WinstonLutzMultiTargetMixin(CloudFileMixin):
    dir_path = ["Winston-Lutz"]
    num_images = 0
    zip = True
    bb_size = 5
    print_results = False
    wl: WinstonLutzMultiTarget

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutzMultiTarget.from_zip(filename, use_filenames=cls.use_filenames)
        else:
            cls.wl = WinstonLutzMultiTarget(filename, use_filenames=cls.use_filenames)
        cls.wl.analyze(bb_size_mm=cls.bb_size, machine_scale=cls.machine_scale)
        if cls.print_results:
            print(cls.wl.results())

    def test_number_of_images(self):
        self.assertEqual(self.num_images, len(self.wl.images))

    def test_bb_max_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="max"), self.cax2bb_max_distance, delta=0.15
        )

    def test_bb_median_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="median"),
            self.cax2bb_median_distance,
            delta=0.1,
        )

    def test_bb_mean_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="mean"), self.cax2bb_mean_distance, delta=0.1
        )


class TestWLLoading(TestCase, FromDemoImageTesterMixin, FromURLTesterMixin):
    klass = WinstonLutz
    demo_load_method = "from_demo_images"
    url = "winston_lutz.zip"

    def test_loading_from_config_mapping(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "noisy_WL_30x5.zip"])
        with TemporaryZipDirectory(path) as z:
            config = {
                "WL G=0, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                    11,
                    12,
                    13,
                ),
                "WL G=90, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                    21,
                    22,
                    23,
                ),
                "WL G=180, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                    31,
                    32,
                    33,
                ),
                "WL G=270, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                    41,
                    42,
                    43,
                ),
            }
            wl = WinstonLutz(z, axis_mapping=config)
        wl.analyze()
        self.assertEqual(wl.images[0].gantry_angle, 11)
        self.assertEqual(wl.images[2].collimator_angle, 32)
        self.assertEqual(wl.images[3].couch_angle, 43)

    def test_loading_from_config_mapping_from_zip(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "noisy_WL_30x5.zip"])
        config = {
            "WL G=0, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                11,
                12,
                13,
            ),
            "WL G=90, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                21,
                22,
                23,
            ),
            "WL G=180, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                31,
                32,
                33,
            ),
            "WL G=270, C=0, P=0; Field=(30, 30)mm; BB=5mm @ left=0, in=0, up=0; Gantry tilt=0, Gantry sag=0.dcm": (
                41,
                42,
                43,
            ),
        }
        wl = WinstonLutz.from_zip(path, axis_mapping=config)
        wl.analyze()
        self.assertEqual(wl.images[0].gantry_angle, 11)
        self.assertEqual(wl.images[2].collimator_angle, 32)
        self.assertEqual(wl.images[3].couch_angle, 43)

    def test_loading_1_image_fails(self):
        with self.assertRaises(ValueError):
            folder = get_folder_from_cloud_test_repo(
                ["Winston-Lutz", "lutz", "1_image"]
            )
            WinstonLutz(folder)

    def test_invalid_dir(self):
        with self.assertRaises(ValueError):
            WinstonLutz(r"nonexistant/dir")

    def test_load_from_file_object(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "noisy_WL_30x5.zip"])
        ref_w = WinstonLutz.from_zip(path)
        ref_w.analyze()
        with open(path, "rb") as f:
            w = WinstonLutz.from_zip(f)
            w.analyze()
        self.assertIsInstance(w, WinstonLutz)
        self.assertEqual(w.gantry_iso_size, ref_w.gantry_iso_size)

    def test_load_from_stream(self):
        path = get_file_from_cloud_test_repo([TEST_DIR, "noisy_WL_30x5.zip"])
        ref_w = WinstonLutz.from_zip(path)
        ref_w.analyze()
        with open(path, "rb") as f:
            s = io.BytesIO(f.read())
            w = WinstonLutz.from_zip(s)
            w.analyze()
        self.assertIsInstance(w, WinstonLutz)
        self.assertEqual(w.gantry_iso_size, ref_w.gantry_iso_size)

    def test_load_2d_from_stream(self):
        path = get_file_from_cloud_test_repo(
            ["Winston-Lutz", "lutz", "1_image", "gantry0.dcm"]
        )
        ref_w = WinstonLutz2D(path)
        ref_w.analyze()
        with open(path, "rb") as f:
            s = io.BytesIO(f.read())
            w = WinstonLutz2D(s)
            w.analyze()
        self.assertIsInstance(w, WinstonLutz2D)
        self.assertEqual(w.bb, ref_w.bb)


class GeneralTests(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze(machine_scale=MachineScale.VARIAN_IEC)

    def test_run_demo(self):
        WinstonLutz.run_demo()  # shouldn't raise

    def test_results(self):
        print(self.wl.results())  # shouldn't raise

    def test_not_yet_analyzed(self):
        wl = WinstonLutz.from_demo_images()
        with self.assertRaises(ValueError):
            wl.results()  # not yet analyzed

        with self.assertRaises(ValueError):
            wl.plot_images()

        with self.assertRaises(ValueError):
            wl.plot_summary()

    def test_str_or_enum(self):
        # shouldn't raise
        self.wl.plot_images("Gantry")
        self.wl.plot_images(Axis.GANTRY)

        self.wl.plot_axis_images("Gantry")
        self.wl.plot_axis_images(Axis.GANTRY)

    def test_bb_override(self):
        with self.assertRaises(ValueError):
            wl = pylinac.WinstonLutz.from_demo_images()
            wl.analyze(bb_size_mm=8)

    def test_bb_shift_instructions(self):
        move = self.wl.bb_shift_instructions()
        self.assertTrue("RIGHT" in move)

        move = self.wl.bb_shift_instructions(couch_vrt=-2, couch_lat=1, couch_lng=100)
        self.assertTrue("RIGHT" in move)
        self.assertTrue("VRT" in move)

    def test_results_data(self):
        data = self.wl.results_data()
        self.assertIsInstance(data, WinstonLutzResult)
        self.assertEqual(
            data.num_couch_images,
            self.wl._get_images(axis=(Axis.COUCH, Axis.REFERENCE))[0],
        )
        self.assertEqual(data.max_2d_cax_to_epid_mm, self.wl.cax2epid_distance("max"))
        self.assertEqual(
            data.median_2d_cax_to_epid_mm, self.wl.cax2epid_distance("median")
        )

        data_dict = self.wl.results_data(as_dict=True)
        self.assertIn("pylinac_version", data_dict)
        self.assertEqual(
            data_dict["gantry_3d_iso_diameter_mm"], self.wl.gantry_iso_size
        )

    def test_bb_too_far_away_fails(self):
        """BB is >20mm from CAX"""
        file = get_file_from_cloud_test_repo([TEST_DIR, "bb_too_far_away.zip"])
        wl = WinstonLutz.from_zip(file)
        with self.assertRaises(ValueError):
            wl.analyze()


class TestPublishPDF(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze()

    def test_publish_pdf(self):
        # normal publish; shouldn't raise
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t)

    def test_publish_w_metadata_and_notes(self):
        with tempfile.TemporaryFile() as t:
            self.wl.publish_pdf(t, notes="stuff", metadata={"Unit": "TB1"})


class TestPlottingSaving(TestCase):
    def setUp(self):
        self.wl = WinstonLutz.from_demo_images()
        self.wl.analyze()

    @classmethod
    def tearDownClass(cls):
        plt.close("all")

    def test_plot(self):
        self.wl.plot_images()  # shouldn't raise
        self.wl.plot_images(axis=Axis.GANTRY)
        self.wl.plot_images(axis=Axis.COLLIMATOR)
        self.wl.plot_images(axis=Axis.COUCH)
        self.wl.plot_images(axis=Axis.GB_COMBO)
        self.wl.plot_images(axis=Axis.GBP_COMBO)

    def test_save_to_stream(self):
        items = self.wl.save_images_to_stream()
        assert isinstance(items, dict)
        assert str(self.wl.images[0]) in items.keys()
        assert len(items) == 15

    def test_plot_split_plots(self):
        figs, names = self.wl.plot_images(show=False, split=True)
        assert isinstance(figs[0], plt.Figure)
        assert isinstance(names[0], str)
        assert len(figs) == 9

    def test_save(self):
        save_file(self.wl.save_summary)
        save_file(self.wl.save_images)

    def test_plot_wo_all_axes(self):
        # test that analyzing images w/o gantry images doesn't fail
        wl_zip = get_file_from_cloud_test_repo([TEST_DIR, "Naming.zip"])
        wl = WinstonLutz.from_zip(wl_zip, use_filenames=True)
        wl.analyze()
        wl.plot_summary()  # shouldn't raise


class WinstonLutzMixin(CloudFileMixin):
    dir_path = ["Winston-Lutz"]
    num_images = 0
    zip = True
    bb_size = 5
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    cax2bb_mean_distance = 0
    epid_deviation = None
    bb_shift_vector = Vector()  # vector to place BB at iso
    machine_scale = MachineScale.IEC61217
    axis_of_rotation = {
        0: Axis.REFERENCE
    }  # fill with as many {image#: known_axis_of_rotation} pairs as desired
    print_results = False
    use_filenames = False

    @classmethod
    def setUpClass(cls):
        filename = cls.get_filename()
        if cls.zip:
            cls.wl = WinstonLutz.from_zip(filename, use_filenames=cls.use_filenames)
        else:
            cls.wl = WinstonLutz(filename, use_filenames=cls.use_filenames)
        cls.wl.analyze(bb_size_mm=cls.bb_size, machine_scale=cls.machine_scale)
        if cls.print_results:
            print(cls.wl.results())
            print(cls.wl.bb_shift_vector)

    def test_number_of_images(self):
        self.assertEqual(self.num_images, len(self.wl.images))

    def test_gantry_iso(self):
        # test iso size
        self.assertAlmostEqual(
            self.wl.gantry_iso_size, self.gantry_iso_size, delta=0.15
        )

    def test_collimator_iso(self):
        # test iso size
        if self.collimator_iso_size is not None:
            self.assertAlmostEqual(
                self.wl.collimator_iso_size, self.collimator_iso_size, delta=0.15
            )

    def test_couch_iso(self):
        # test iso size
        if self.couch_iso_size is not None:
            self.assertAlmostEqual(
                self.wl.couch_iso_size, self.couch_iso_size, delta=0.15
            )

    def test_epid_deviation(self):
        if self.epid_deviation is not None:
            self.assertAlmostEqual(
                max(self.wl.axis_rms_deviation(Axis.EPID)),
                self.epid_deviation,
                delta=0.15,
            )

    def test_bb_max_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="max"), self.cax2bb_max_distance, delta=0.15
        )

    def test_bb_median_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="median"),
            self.cax2bb_median_distance,
            delta=0.1,
        )

    def test_bb_mean_distance(self):
        self.assertAlmostEqual(
            self.wl.cax2bb_distance(metric="mean"), self.cax2bb_mean_distance, delta=0.1
        )

    def test_bb_shift_vector(self):
        self.assertTrue(
            vector_is_close(self.wl.bb_shift_vector, self.bb_shift_vector, delta=0.15),
            msg="The vector {} is not sufficiently close to vector {}".format(
                self.wl.bb_shift_vector, self.bb_shift_vector
            ),
        )

    def test_known_axis_of_rotation(self):
        for idx, axis in self.axis_of_rotation.items():
            self.assertEqual(axis, self.wl.images[idx].variable_axis)


class WLDemo(WinstonLutzMixin, TestCase):
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1.2
    couch_iso_size = 2.3
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.7
    cax2bb_mean_distance = 0.6
    machine_scale = MachineScale.VARIAN_IEC
    epid_deviation = 1.3
    axis_of_rotation = {0: Axis.REFERENCE}
    bb_shift_vector = Vector(x=0.4, y=-0.4, z=-0.2)
    delete_file = False

    @classmethod
    def setUpClass(cls):
        cls.wl = WinstonLutz.from_demo_images()
        cls.wl.analyze(machine_scale=cls.machine_scale)

    def test_different_scale_has_different_shift(self):
        assert "RIGHT" in self.wl.bb_shift_instructions()
        assert self.wl.bb_shift_vector.x > 0.1
        new_wl = WinstonLutz.from_demo_images()
        new_wl.analyze(machine_scale=MachineScale.IEC61217)
        assert new_wl.bb_shift_vector.x < 0.1
        assert "LEFT" in new_wl.bb_shift_instructions()

    def test_multiple_analyses_gives_same_result(self):
        original_vector = copy.copy(self.wl.bb_shift_vector)
        # re-analyze w/ same settings
        self.wl.analyze(machine_scale=self.machine_scale)
        new_vector = self.wl.bb_shift_vector
        assert vector_is_close(original_vector, new_vector, delta=0.05)


class WLPerfect30x8(WinstonLutzMixin, TestCase):
    """30x30mm field, 8mm BB"""

    file_name = "perfect_WL_30x8.zip"
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLPerfect30x2(WinstonLutzMixin, TestCase):
    """30x30mm field, 2mm BB"""

    file_name = "perfect_WL_30x2mm.zip"
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()
    bb_size = 2


class WLPerfect10x4(WinstonLutzMixin, TestCase):
    """10x10mm field, 4mm BB"""

    file_name = "perfect_WL_10x4.zip"
    num_images = 4
    gantry_iso_size = 0
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLNoisy30x5(WinstonLutzMixin, TestCase):
    """30x30mm field, 5mm BB. S&P noise added"""

    file_name = "noisy_WL_30x5.zip"
    num_images = 4
    gantry_iso_size = 0.08
    collimator_iso_size = 0
    couch_iso_size = 0
    cax2bb_max_distance = 0
    cax2bb_median_distance = 0
    epid_deviation = 0
    bb_shift_vector = Vector()


class WLLateral3mm(WinstonLutzMixin, TestCase):
    # verified independently
    file_name = "lat3mm.zip"
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    cax2bb_mean_distance = 2.3
    bb_shift_vector = Vector(x=-3.6, y=0.5, z=0.6)


class WLLongitudinal3mm(WinstonLutzMixin, TestCase):
    # verified independently
    file_name = "lng3mm.zip"
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.9
    cax2bb_median_distance = 3.7
    cax2bb_mean_distance = 3.7
    bb_shift_vector = Vector(x=-0.63, y=3.6, z=0.6)


class WLVertical3mm(WinstonLutzMixin, TestCase):
    file_name = "vrt3mm.zip"
    num_images = 4
    gantry_iso_size = 0.5
    cax2bb_max_distance = 3.8
    cax2bb_median_distance = 2.3
    cax2bb_mean_distance = 2.3
    bb_shift_vector = Vector(x=-0.5, y=0.5, z=3.6)
    print_results = True


class WLDontUseFileNames(WinstonLutzMixin, TestCase):
    file_name = "Naming.zip"
    num_images = 4
    gantry_iso_size = 0.3
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    cax2bb_mean_distance = 0.8
    bb_shift_vector = Vector(x=-0.4, y=0.6, z=0.6)
    axis_of_rotation = {
        0: Axis.REFERENCE,
        1: Axis.GANTRY,
        2: Axis.GANTRY,
        3: Axis.GANTRY,
    }


class WLUseFileNames(WinstonLutzMixin, TestCase):
    file_name = "Naming.zip"
    use_filenames = True
    num_images = 4
    collimator_iso_size = 1.2
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.8
    cax2bb_mean_distance = 0.8
    bb_shift_vector = Vector(y=0.6)
    axis_of_rotation = {
        0: Axis.COLLIMATOR,
        1: Axis.COLLIMATOR,
        2: Axis.COLLIMATOR,
        3: Axis.COLLIMATOR,
    }


class WLBadFilenames(TestCase):
    def test_bad_filenames(self):
        # tests_basic that using filenames with incorrect syntax will fail
        wl_dir = get_file_from_cloud_test_repo([TEST_DIR, "Bad-Names.zip"])
        with self.assertRaises(ValueError):
            wl = WinstonLutz.from_zip(wl_dir, use_filenames=True)
            wl.analyze()


class KatyiX0(WinstonLutzMixin, TestCase):
    # independently verified
    file_name = ["Katy iX", "0.zip"]
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.8
    cax2bb_mean_distance = 0.7
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=-0.5, y=0.4, z=-0.5)
    print_results = True


class KatyiX1(WinstonLutzMixin, TestCase):
    file_name = ["Katy iX", "1.zip"]
    num_images = 17
    gantry_iso_size = 1.1
    collimator_iso_size = 0.7
    couch_iso_size = 0.6
    cax2bb_max_distance = 1.2
    cax2bb_median_distance = 0.3
    cax2bb_mean_distance = 0.4
    bb_shift_vector = Vector(x=0.3, y=-0.2, z=0.3)


class KatyiX2(WinstonLutzMixin, TestCase):
    file_name = ["Katy iX", "2.zip"]
    num_images = 17
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.5
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.5
    cax2bb_mean_distance = 0.6
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=0.4, y=-0.1, z=0.1)


class KatyiX3(WinstonLutzMixin, TestCase):
    file_name = ["Katy iX", "3 (with crosshair).zip"]
    num_images = 17
    gantry_iso_size = 1.1
    collimator_iso_size = 1.3
    couch_iso_size = 1.8
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.8
    cax2bb_mean_distance = 0.75
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=-0.3, y=0.4, z=-0.5)


class KatyTB0(WinstonLutzMixin, TestCase):
    file_name = ["Katy TB", "0.zip"]
    num_images = 17
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.07
    cax2bb_median_distance = 0.8
    cax2bb_mean_distance = 0.8
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=-0.7, y=-0.1, z=-0.2)


class KatyTB1(WinstonLutzMixin, TestCase):
    file_name = ["Katy TB", "1.zip"]
    num_images = 16
    gantry_iso_size = 0.9
    collimator_iso_size = 0.8
    couch_iso_size = 1.1
    cax2bb_max_distance = 1
    cax2bb_median_distance = 0.7
    cax2bb_mean_distance = 0.6
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=-0.6, y=-0.2)


class KatyTB2(WinstonLutzMixin, TestCase):
    file_name = ["Katy TB", "2.zip"]
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 0.7
    couch_iso_size = 0.7
    cax2bb_max_distance = 1.1
    cax2bb_median_distance = 0.4
    cax2bb_mean_distance = 0.5
    bb_shift_vector = Vector(x=0.0, y=-0.2, z=-0.6)


class ChicagoTBFinal(WinstonLutzMixin, TestCase):
    # verified independently
    file_name = ["Chicago", "WL-Final_C&G&C_Final.zip"]
    num_images = 17
    gantry_iso_size = 0.91
    collimator_iso_size = 0.1
    couch_iso_size = 0.3
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    cax2bb_mean_distance = 0.3
    bb_shift_vector = Vector(y=0.1)


class ChicagoTB52915(WinstonLutzMixin, TestCase):
    file_name = ["Chicago", "WL_05-29-15_Final.zip"]
    num_images = 16
    gantry_iso_size = 0.6
    collimator_iso_size = 0.3
    couch_iso_size = 0.3
    cax2bb_max_distance = 0.5
    cax2bb_median_distance = 0.3
    cax2bb_mean_distance = 0.3
    bb_shift_vector = Vector(z=0.2)


class TrueBeam3120213(WinstonLutzMixin, TestCase):
    file_name = ["TrueBeam 3", "120213.zip"]
    num_images = 26
    cax2bb_max_distance = 0.9
    cax2bb_median_distance = 0.35
    cax2bb_mean_distance = 0.4
    gantry_iso_size = 1.1
    collimator_iso_size = 0.7
    couch_iso_size = 0.7
    bb_shift_vector = Vector(x=-0.1, y=-0.2, z=0.2)


class SugarLandiX2015(WinstonLutzMixin, TestCase):
    file_name = ["Sugarland iX", "2015", "Lutz2.zip"]
    num_images = 17
    gantry_iso_size = 1.3
    collimator_iso_size = 0.5
    couch_iso_size = 1.3
    cax2bb_max_distance = 1.67
    cax2bb_median_distance = 1.05
    cax2bb_mean_distance = 1.1
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=0.4, y=-0.7, z=0.1)


class BayAreaiX0(WinstonLutzMixin, TestCase):
    file_name = ["Bay Area iX", "0.zip"]
    num_images = 17
    gantry_iso_size = 1
    collimator_iso_size = 1.1
    couch_iso_size = 2.3
    cax2bb_max_distance = 1.25
    cax2bb_median_distance = 0.6
    cax2bb_mean_distance = 0.6
    machine_scale = MachineScale.VARIAN_IEC
    bb_shift_vector = Vector(x=0.3, y=-0.4, z=-0.2)


class DAmoursElektaOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""

    file_name = ["Michel DAmours - WLGantry_Offset_x=-1cm,y=+1cm,z=-1cm.zip"]
    num_images = 8
    gantry_iso_size = 1.1
    cax2bb_max_distance = 17.5
    cax2bb_median_distance = 14.3
    cax2bb_mean_distance = 13.8
    bb_shift_vector = Vector(x=10.2, y=-9.2, z=-11.1)  # independently verified


class DAmoursElektaXOffset(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""

    file_name = ["Michel D'Amours - WL_Shift_x=+1cm.zip"]
    num_images = 8
    gantry_iso_size = 1.1
    cax2bb_max_distance = 9.5
    cax2bb_median_distance = 6.9
    cax2bb_mean_distance = 6
    bb_shift_vector = Vector(x=-9.5, y=0.3, z=0.1)  # independently verified


class DAmoursElektaCentered(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""

    file_name = ["Michel D'Amours - GantryWL_BBCentered.zip"]
    num_images = 8
    gantry_iso_size = 1.1
    collimator_iso_size = None
    couch_iso_size = None
    cax2bb_max_distance = 0.8
    cax2bb_median_distance = 0.5
    cax2bb_mean_distance = 0.6
    bb_shift_vector = Vector(y=0.4)


class DeBr6XElekta(WinstonLutzMixin, TestCase):
    """An Elekta dataset, with the BB centered."""

    file_name = ["DeBr", "6X_Elekta_Ball_Bearing.zip"]
    num_images = 8
    gantry_iso_size = 1.8
    collimator_iso_size = 1.8
    couch_iso_size = None
    cax2bb_max_distance = 1.0
    cax2bb_median_distance = 0.7
    cax2bb_mean_distance = 0.6
    bb_shift_vector = Vector(x=0.4, y=-0.2)


class LargeFieldCouchPresent(WinstonLutzMixin, TestCase):
    """A very large field where the couch is present"""

    file_name = ["large_field_couch_present.zip"]
    num_images = 4
    gantry_iso_size = 0.8
    collimator_iso_size = None
    couch_iso_size = None
    cax2bb_max_distance = 1.3
    cax2bb_median_distance = 1
    cax2bb_mean_distance = 1
    bb_shift_vector = Vector(x=0.5, y=-0.7, z=0.8)
