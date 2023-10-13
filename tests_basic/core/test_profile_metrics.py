from unittest import TestCase

from pydicom import Dataset

from pylinac.core.image import DicomImage
from pylinac.core.image_generator import AS1000Image, PerfectBBLayer, RandomNoiseLayer
from pylinac.core.image_generator.layers import (
    FilteredFieldLayer,
    GaussianFilterLayer,
    Layer,
)
from pylinac.core.metrics import DiskLocator


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
