"""Travis CI memory can't handle all the starshots; thus only test them when explicitly asked to."""
import concurrent.futures
import os
import os.path as osp
from unittest import TestCase
import time

from pylinac.core.image import Image
from tests.test_starshot import StarMixin, test_file_dir, Point, Starshot


class Starshot2(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#2.tif')
    wobble_center = Point(566, 590)
    wobble_diameter_mm = 0.2
    num_rad_lines = 4


class Starshot3(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#3.tif')
    wobble_center = Point(466, 595)
    wobble_diameter_mm = 0.32
    num_rad_lines = 6


class Starshot4(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#4.tif')
    wobble_center = Point(446, 565)
    wobble_diameter_mm = 0.38
    num_rad_lines = 6


class Starshot5(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#5.tif')
    wobble_center = Point(557, 580)
    wobble_diameter_mm = 0.15
    num_rad_lines = 4


class Starshot6(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#6.tif')
    wobble_center = Point(528, 607)
    wobble_diameter_mm = 0.3
    num_rad_lines = 7


class Starshot7(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#7.tif')
    wobble_center = Point(469, 646)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot8(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#8.tiff')
    wobble_center = Point(686, 669)
    wobble_diameter_mm = 0.35
    num_rad_lines = 5


class Starshot9(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#9.tiff')
    wobble_center = Point(714, 611)
    wobble_diameter_mm = 0.3
    num_rad_lines = 5


class Starshot10(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#10.tiff')
    wobble_center = Point(725, 802)
    wobble_diameter_mm = 0.65
    num_rad_lines = 5


class Starshot11(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#11.tiff')
    wobble_center = Point(760, 650)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot12(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#12.tiff')
    wobble_center = Point(315, 292)
    wobble_diameter_mm = 0.88
    num_rad_lines = 4


class Starshot13(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#13.tiff')
    wobble_center = Point(376, 303)
    wobble_diameter_mm = 0.25
    num_rad_lines = 4


class Starshot14(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#14.tiff')
    wobble_center = Point(334, 282)
    wobble_diameter_mm = 0.55
    num_rad_lines = 4


class Starshot15(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#15.tiff')
    wobble_center = Point(346, 309)
    wobble_diameter_mm = 0.6
    num_rad_lines = 4


class Starshot16(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#16.tiff')
    wobble_center = Point(1444, 1452)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot17(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#17.tiff')
    wobble_center = Point(1475, 1361)
    wobble_diameter_mm = 0.44
    num_rad_lines = 6


class Starshot18(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#18.tiff')
    wobble_center = Point(1516, 1214)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot19(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#19.tiff')
    wobble_center = Point(1475, 1276)
    wobble_diameter_mm = 0.6
    num_rad_lines = 6


class Starshot20(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#20.tiff')
    wobble_center = Point(347, 328)
    wobble_diameter_mm = 0.71
    num_rad_lines = 4


class Starshot21(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#21.tiff')
    wobble_center = Point(354, 294)
    wobble_diameter_mm = 1.4
    wobble_tolerance = 0.5
    num_rad_lines = 4
    passes = False


class Starshot22(StarMixin, TestCase):
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


class Starshot23(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#23.tiff')
    wobble_center = Point(1297, 1699)
    wobble_diameter_mm = 0.38
    num_rad_lines = 9


class Starshot24(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#24.tiff')
    wobble_center = Point(1370, 1454)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot25(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#25.tiff')
    wobble_center = Point(286, 279)
    wobble_diameter_mm = 0.3
    num_rad_lines = 4


class Starshot26(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#26.tiff')
    wobble_center = Point(1511, 1452)
    wobble_diameter_mm = 0.5
    num_rad_lines = 4


class Starshot27(StarMixin, TestCase):
    star_file = osp.join(test_file_dir, 'Starshot#27.tiff')
    wobble_center = Point(1105, 1306)
    wobble_diameter_mm = 0.4
    num_rad_lines = 6


def run_star(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        mystar = Starshot(path)
        mystar.analyze()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


class TestImageBank(TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""
    image_bank_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Starshots'))

    def test_all(self):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.image_bank_dir):
                for file in files:
                    filepath = osp.join(pdir, file)
                    try:
                        Image.load(filepath)
                    except:
                        pass
                    else:
                        future = exec.submit(run_star, filepath)
                        futures.append(future)
            for future in concurrent.futures.as_completed(futures):
                print(future.result())
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))
