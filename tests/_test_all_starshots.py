"""Travis CI memory can't handle all the starshots; thus only test them when explicitly asked to."""
import unittest
import os.path as osp

from tests.test_starshot import StarMixin, test_file_dir, Point


class Starshot2(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#2.tif')
    wobble_center = Point(566, 590)
    wobble_diameter_mm = 0.2
    num_rad_lines = 4


class Starshot3(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#3.tif')
    wobble_center = Point(466, 595)
    wobble_diameter_mm = 0.32
    num_rad_lines = 6


class Starshot4(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#4.tif')
    wobble_center = Point(446, 565)
    wobble_diameter_mm = 0.38
    num_rad_lines = 6


class Starshot5(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#5.tif')
    wobble_center = Point(557, 580)
    wobble_diameter_mm = 0.15
    num_rad_lines = 4


class Starshot6(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#6.tif')
    wobble_center = Point(528, 607)
    wobble_diameter_mm = 0.3
    num_rad_lines = 7


class Starshot7(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#7.tif')
    wobble_center = Point(469, 646)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot8(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#8.tiff')
    wobble_center = Point(686, 669)
    wobble_diameter_mm = 0.35
    num_rad_lines = 5


class Starshot9(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#9.tiff')
    wobble_center = Point(714, 611)
    wobble_diameter_mm = 0.3
    num_rad_lines = 5


class Starshot10(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#10.tiff')
    wobble_center = Point(725, 802)
    wobble_diameter_mm = 0.65
    num_rad_lines = 5


class Starshot11(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#11.tiff')
    wobble_center = Point(760, 650)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot12(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#12.tiff')
    wobble_center = Point(315, 292)
    wobble_diameter_mm = 0.88
    num_rad_lines = 4


class Starshot13(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#13.tiff')
    wobble_center = Point(376, 303)
    wobble_diameter_mm = 0.25
    num_rad_lines = 4


class Starshot14(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#14.tiff')
    wobble_center = Point(334, 282)
    wobble_diameter_mm = 0.55
    num_rad_lines = 4


class Starshot15(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#15.tiff')
    wobble_center = Point(346, 309)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot16(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#16.tiff')
    wobble_center = Point(1444, 1452)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot17(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#17.tiff')
    wobble_center = Point(1475, 1361)
    wobble_diameter_mm = 0.44
    num_rad_lines = 6


class Starshot18(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#18.tiff')
    wobble_center = Point(1516, 1214)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot19(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#19.tiff')
    wobble_center = Point(1475, 1276)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot20(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#20.tiff')
    wobble_center = Point(347, 328)
    wobble_diameter_mm = 0.71
    num_rad_lines = 4


class Starshot21(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#21.tiff')
    wobble_center = Point(354, 294)
    wobble_diameter_mm = 1.4
    wobble_tolerance = 0.5
    num_rad_lines = 4
    passes = False


class Starshot22(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#22.tiff')
    wobble_center = Point(1305, 1513)
    wobble_diameter_mm = 0.8
    wobble_tolerance = 0.3
    num_rad_lines = 9

    def test_bad_input_no_recursion_fails(self):
        """Test that without recursion, a bad setup fails."""
        with self.assertRaises(RuntimeError):
            self.star.analyze(radius=0.3, min_peak_height=0.95, recursive=False)

        # but will pass with recursion
        self.star.analyze()
        self.test_passed()


class Starshot23(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#23.tiff')
    wobble_center = Point(1297, 1699)
    wobble_diameter_mm = 0.38
    num_rad_lines = 9


class Starshot24(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#24.tiff')
    wobble_center = Point(1370, 1454)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot25(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#25.tiff')
    wobble_center = Point(286, 279)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot26(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#26.tiff')
    wobble_center = Point(1511, 1452)
    wobble_diameter_mm = 0.5
    num_rad_lines = 4


class Starshot27(StarMixin, unittest.TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#27.tiff')
    wobble_center = Point(1105, 1306)
    wobble_diameter_mm = 0.4
    num_rad_lines = 6

