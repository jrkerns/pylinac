from unittest import TestCase

from pydicom import Dataset

from pylinac.core.geometry import Point
from pylinac.core.image import XIM, DicomImage
from pylinac.core.image_generator import AS1000Image, PerfectBBLayer, RandomNoiseLayer
from pylinac.core.image_generator.layers import (
    FilteredFieldLayer,
    GaussianFilterLayer,
    Layer,
)
from pylinac.core.metrics import (
    DiskLocator,
    GlobalDiskLocator,
    GlobalSizedFieldLocator,
    deduplicate_points_and_boundaries,
)
from tests_basic.utils import get_file_from_cloud_test_repo


def create_bb_image(
    field_size=(50, 50),
    bb_size=5,
    offset=(0, 0),
    field: type[Layer] = FilteredFieldLayer,
) -> Dataset:
    as1000 = AS1000Image(
        sid=1000
    )  # this will set the pixel size and shape automatically
    as1000.add_layer(
        field(field_size_mm=field_size, cax_offset_mm=(0, 0))
    )  # create a 50x50mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=2)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    as1000.add_layer(PerfectBBLayer(bb_size_mm=bb_size, cax_offset_mm=offset))
    as1000.add_layer(RandomNoiseLayer())
    return as1000.as_dicom()


def create_open_field_image(
    field_size=(50, 50),
    offset=(0, 0),
    field: type[Layer] = FilteredFieldLayer,
) -> Dataset:
    as1000 = AS1000Image(
        sid=1000
    )  # this will set the pixel size and shape automatically
    as1000.add_layer(
        field(field_size_mm=field_size, cax_offset_mm=offset)
    )  # create a 50x50mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=2)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    as1000.add_layer(RandomNoiseLayer())
    return as1000.as_dicom()


def create_multi_open_field(
    field_sizes=((50, 50), (50, 50)),
    offsets=((0, 0), (40, 40)),
    field: type[Layer] = FilteredFieldLayer,
) -> Dataset:
    as1000 = AS1000Image(
        sid=1000
    )  # this will set the pixel size and shape automatically
    for field_size, offset in zip(field_sizes, offsets):
        as1000.add_layer(
            field(field_size_mm=field_size, cax_offset_mm=offset)
        )  # create a 50x50mm square field
    as1000.add_layer(
        GaussianFilterLayer(sigma_mm=2)
    )  # add an image-wide gaussian to simulate penumbra/scatter
    as1000.add_layer(RandomNoiseLayer())
    return as1000.as_dicom()


class TestGeneralMetric(TestCase):
    def test_modifying_image_fails(self):
        # user should not modify the image under the hood
        class DiskModifier(DiskLocator):
            def calculate(self) -> Point:
                calc = super().calculate()
                self.image.crop(10)  # modifying the image
                return calc

        bb_diameter_mm = 5
        ds = create_bb_image(bb_size=bb_diameter_mm)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(RuntimeError):
            img.compute(
                metrics=[
                    DiskModifier(
                        expected_position=(511.5, 383.5),
                        search_window=(50, 50),
                        radius=6,
                        radius_tolerance=1,
                    )
                ]
            )


class TestGlobalDiskLocator(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        path = get_file_from_cloud_test_repo(["IsoCal-kV-08.xim"])
        xim = XIM(path)
        cls.xim = xim

    def test_isocal_image_default(self):
        bbs = self.xim.compute(
            metrics=GlobalDiskLocator(radius_mm=3.5, radius_tolerance_mm=1)
        )
        self.assertEqual(len(bbs), 15)

    def test_min_number_fails(self):
        # we should get a failure if the min number isn't found
        with self.assertRaises(ValueError):
            self.xim.compute(
                metrics=GlobalDiskLocator(
                    radius_mm=3.5, radius_tolerance_mm=1, min_number=16
                )
            )

    def test_max_number_is_close(self):
        # the max number is soft, but we can stop early if we find at least the max number
        bbs = self.xim.compute(
            metrics=GlobalDiskLocator(
                radius_mm=3.5, radius_tolerance_mm=1, max_number=5
            )
        )
        self.assertEqual(len(bbs), 5)

    def test_none_found(self):
        # set the radius too high to test that no bbs are found
        with self.assertRaises(ValueError):
            self.xim.compute(
                metrics=GlobalDiskLocator(radius_mm=10, radius_tolerance_mm=1)
            )

    def test_min_separation(self):
        # set the separation to 35mm. this will reduce the number found by 2 compared to no separation
        bbs = self.xim.compute(
            metrics=GlobalDiskLocator(
                radius_mm=3.5, radius_tolerance_mm=1, min_separation_mm=35
            )
        )
        self.assertEqual(len(bbs), 13)

    def test_min_separation_massive(self):
        # with a massive separation, we should get 1 bb
        bbs = self.xim.compute(
            metrics=GlobalDiskLocator(
                radius_mm=3.5, radius_tolerance_mm=1, min_separation_mm=1000
            )
        )
        self.assertEqual(len(bbs), 1)


class TestDeduplicatePoints(TestCase):
    def test_normal(self):
        # normal case: no duplicates
        original_points = [Point(1, 1), Point(5, 5)]
        new_points = [Point(10, 10)]
        deduped = deduplicate_points_and_boundaries(
            original_points, new_points, min_separation_px=3
        )
        self.assertEqual(len(deduped), 3)

    def test_no_recursion(self):
        # tests we don't recurse infinitely
        original_points = [Point(1, 1), Point(5, 5)]
        deduped = deduplicate_points_and_boundaries(
            original_points, original_points, min_separation_px=3
        )
        self.assertEqual(len(deduped), 2)

    def test_exclude_from_original_list(self):
        # tests that points in the original list are not added to the new list
        original_points = [Point(1, 1), Point(5, 5)]
        deduped = deduplicate_points_and_boundaries(
            [], original_points, min_separation_px=30
        )
        self.assertEqual(len(deduped), 1)


class TestDiskLocatorPixels(TestCase):
    def test_perfect_image(self):
        bb_diameter_mm = 5
        ds = create_bb_image(bb_size=bb_diameter_mm)
        img = DicomImage.from_dataset(ds)
        position = img.compute(
            metrics=[
                DiskLocator(
                    expected_position=(511.5, 383.5),
                    search_window=(50, 50),
                    radius=6,
                    radius_tolerance=1,
                )
            ]
        )
        self.assertAlmostEqual(position.x, 511.5, delta=1)
        self.assertAlmostEqual(position.y, 383.5, delta=1)

    def test_wrong_area(self):
        # if the user selects a position without a bb it should raise an error
        ds = create_bb_image(offset=(20, 20))
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator(
                        expected_position=(511.5, 383.5),
                        search_window=(10, 10),
                        radius=20,
                        radius_tolerance=1,
                    )
                ]
            )

    def test_bb_too_small(self):
        ds = create_bb_image(bb_size=1)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator(
                        expected_position=(511.5, 383.5),
                        search_window=(10, 10),
                        radius=20,
                        radius_tolerance=2,
                    )
                ]
            )


class TestDiskLocatorPhysical(TestCase):
    def test_perfect_image(self):
        bb_diameter = 5
        ds = create_bb_image(bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        position = img.compute(
            metrics=[
                DiskLocator.from_physical(
                    expected_position_mm=(200, 150),
                    search_window_mm=(10, 10),
                    radius_mm=bb_diameter / 2,
                    radius_tolerance_mm=1,
                )
            ]
        )
        self.assertAlmostEqual(position.x, 511.5, delta=1)
        self.assertAlmostEqual(position.y, 383.5, delta=1)

    def test_wrong_area(self):
        # if the user selects a position without a bb it should raise an error
        bb_diameter = 5
        ds = create_bb_image(offset=(20, 20), bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator.from_physical(
                        expected_position_mm=(200, 150),
                        search_window_mm=(10, 10),
                        radius_mm=bb_diameter / 2,
                        radius_tolerance_mm=1,
                    )
                ]
            )

    def test_bb_too_small(self):
        bb_diameter = 5
        ds = create_bb_image(bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator.from_physical(
                        expected_position_mm=(200, 150),
                        search_window_mm=(10, 10),
                        radius_mm=bb_diameter,
                        radius_tolerance_mm=1,
                    )
                ]
            )

    def test_barely_too_small(self):
        bb_diameter = 10
        ds = create_bb_image(bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            # expected radius is 2mm too big w/ tolerance of only 1mm
            img.compute(
                metrics=[
                    DiskLocator.from_physical(
                        expected_position_mm=(200, 150),
                        search_window_mm=(20, 20),
                        radius_mm=bb_diameter / 2 + 2,
                        radius_tolerance_mm=1,
                    )
                ]
            )

    def test_barely_too_big(self):
        bb_diameter = 10
        ds = create_bb_image(bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            # expected radius is 2mm too small w/ tolerance of only 1mm
            img.compute(
                metrics=[
                    DiskLocator.from_physical(
                        expected_position_mm=(200, 150),
                        search_window_mm=(20, 20),
                        radius_mm=bb_diameter / 2 - 2,
                        radius_tolerance_mm=1,
                    )
                ]
            )


class TestDiskLocatorCenterPixels(TestCase):
    def test_perfect_image(self):
        bb_diameter_mm = 5
        ds = create_bb_image(bb_size=bb_diameter_mm)
        img = DicomImage.from_dataset(ds)
        position = img.compute(
            metrics=[
                DiskLocator.from_center(
                    expected_position=(0, 0),
                    search_window=(50, 50),
                    radius=6,
                    radius_tolerance=1,
                )
            ]
        )
        self.assertAlmostEqual(position.x, 511.5, delta=1)
        self.assertAlmostEqual(position.y, 383.5, delta=1)

    def test_wrong_area(self):
        # if the user selects a position without a bb it should raise an error
        ds = create_bb_image(offset=(20, 20))
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator.from_center(
                        expected_position=(100, 100),
                        search_window=(10, 10),
                        radius=20,
                        radius_tolerance=1,
                    )
                ]
            )

    def test_bb_too_small(self):
        ds = create_bb_image(bb_size=2)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator.from_center(
                        expected_position=(0, 0),
                        search_window=(10, 10),
                        radius=15,
                        radius_tolerance=2,
                    )
                ]
            )


class TestDiskLocatorCenterPhysical(TestCase):
    def test_perfect_image(self):
        bb_diameter = 5
        ds = create_bb_image(bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        position = img.compute(
            metrics=[
                DiskLocator.from_center_physical(
                    expected_position_mm=(0, 0),
                    search_window_mm=(10, 10),
                    radius_mm=bb_diameter / 2,
                    radius_tolerance_mm=1,
                )
            ]
        )
        self.assertAlmostEqual(position.x, 511.5, delta=1)
        self.assertAlmostEqual(position.y, 383.5, delta=1)

    def test_wrong_area(self):
        # if the user selects a position without a bb it should raise an error
        bb_diameter = 5
        ds = create_bb_image(offset=(20, 20), bb_size=bb_diameter)
        img = DicomImage.from_dataset(ds)
        with self.assertRaises(ValueError):
            img.compute(
                metrics=[
                    DiskLocator.from_center_physical(
                        expected_position_mm=(0, 0),
                        search_window_mm=(10, 10),
                        radius_mm=bb_diameter / 2,
                        radius_tolerance_mm=1,
                    )
                ]
            )

    def test_shifted_bb(self):
        # both bb and expected position are shifted
        bb_diameter = 5
        ds = create_bb_image(bb_size=bb_diameter, offset=(20, 20), field_size=(75, 75))
        img = DicomImage.from_dataset(ds)
        img.compute(
            metrics=[
                DiskLocator.from_center_physical(
                    expected_position_mm=(20, 20),
                    search_window_mm=(10, 10),
                    radius_mm=bb_diameter / 2,
                    radius_tolerance_mm=1,
                )
            ]
        )


class TestGlobalSizedFieldLocator(TestCase):
    def test_perfect_image(self):
        ds = create_open_field_image(field_size=(60, 60))
        img = DicomImage.from_dataset(ds)
        fields = img.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                field_width_mm=60,
                field_height_mm=60,
                field_tolerance_mm=2,
                max_number=1,
            )
        )
        self.assertAlmostEqual(fields[0].x, 511.5, delta=1)
        self.assertAlmostEqual(fields[0].y, 383.5, delta=1)

    def test_small_image(self):
        ds = create_open_field_image(field_size=(10, 10))
        img = DicomImage.from_dataset(ds)
        fields = img.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                field_width_mm=10,
                field_height_mm=10,
                field_tolerance_mm=2,
                max_number=1,
            )
        )
        self.assertAlmostEqual(fields[0].x, 511.5, delta=1)
        self.assertAlmostEqual(fields[0].y, 383.5, delta=1)

    def test_large_image(self):
        ds = create_open_field_image(field_size=(250, 250))
        img = DicomImage.from_dataset(ds)
        fields = img.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                field_width_mm=250,
                field_height_mm=250,
                field_tolerance_mm=5,
                max_number=1,
            )
        )
        self.assertAlmostEqual(fields[0].x, 511.5, delta=1)
        self.assertAlmostEqual(fields[0].y, 383.5, delta=1)

    def test_multiple_medium_fields(self):
        ds = create_multi_open_field(
            field_sizes=((30, 30), (30, 30)), offsets=((0, 0), (40, 40))
        )
        img = DicomImage.from_dataset(ds)
        fields = img.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                field_width_mm=30,
                field_height_mm=30,
                field_tolerance_mm=5,
                max_number=2,
            )
        )
        self.assertAlmostEqual(fields[0].x, 511.5, delta=1)
        self.assertAlmostEqual(fields[0].y, 383.5, delta=1)
        self.assertAlmostEqual(fields[1].x, 613.6, delta=1)
        self.assertAlmostEqual(fields[1].y, 485.6, delta=1)

    def test_lots_of_fields(self):
        ds = create_multi_open_field(
            field_sizes=((30, 30), (30, 30), (30, 30), (30, 30)),
            offsets=((0, 0), (40, 40), (-40, 0), (-60, -60)),
        )
        img = DicomImage.from_dataset(ds)
        fields = img.compute(
            metrics=GlobalSizedFieldLocator.from_physical(
                field_width_mm=30,
                field_height_mm=30,
                field_tolerance_mm=5,
                max_number=4,
            )
        )
        self.assertEqual(len(fields), 4)

    def test_not_from_physical(self):
        ds = create_open_field_image(field_size=(10, 10))
        img = DicomImage.from_dataset(ds)
        dpmm = img.dpmm
        fields = img.compute(
            metrics=GlobalSizedFieldLocator(
                field_width_px=10 * dpmm,
                field_height_px=10 * dpmm,
                field_tolerance_px=2 * dpmm,
                max_number=1,
            )
        )
        self.assertAlmostEqual(fields[0].x, 511.5, delta=1)
        self.assertAlmostEqual(fields[0].y, 383.5, delta=1)
