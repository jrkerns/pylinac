import copy
import io
import tempfile
from unittest import TestCase

import numpy as np
import numpy.testing
import pydicom

import tests_basic  # noqa; load env settings and filters
from pylinac import Interpolation, Normalization
from pylinac.core.image import DicomImage, load
from pylinac.core.image_generator import (
    AS500Image,
    AS1000Image,
    AS1200Image,
    ConstantLayer,
    FilteredFieldLayer,
    FilterFreeFieldLayer,
    GaussianFilterLayer,
    PerfectBBLayer,
    PerfectConeLayer,
    PerfectFieldLayer,
    RandomNoiseLayer,
)
from pylinac.core.image_generator.layers import Layer, SlopeLayer, clip_add, even_round
from pylinac.core.image_generator.simulators import Simulator
from pylinac.core.profile import SingleProfile
from pylinac.metrics.image import GlobalFieldLocator

np.random.seed(1234)  # noqa


class TestClipAdd(TestCase):
    def test_clip_add_normal(self):
        image1 = np.zeros((10, 10), dtype=np.uint16)
        image2 = np.ones((10, 10), dtype=np.uint16)
        output = clip_add(image1, image2, dtype=np.uint16)
        self.assertEqual(output.dtype, np.uint16)
        self.assertEqual(output.shape, image1.shape)
        numpy.testing.assert_array_equal(
            output, image2
        )  # arrays are equal because image1 is zeros

    def test_clip_doesnt_flip_bit(self):
        image1 = np.zeros((10, 10), dtype=np.uint16)
        image1.fill(np.iinfo(np.uint16).max)  # set image to max value
        image2 = np.ones((10, 10), dtype=np.uint16)
        output = clip_add(image1, image2, dtype=np.uint16)
        # adding 1 to an array at the max would normally flip bits; ensure it doesn't
        self.assertEqual(output.dtype, np.uint16)
        numpy.testing.assert_array_equal(
            output, image1
        )  # output is same as image1 because image2 didn't actually add anything


class TestEvenRound(TestCase):
    def test_even_round(self):
        self.assertEqual(even_round(3), 4)
        self.assertEqual(even_round(2), 2)
        self.assertEqual(even_round(15), 16)


def profiles_from_simulator(
    simulator: Simulator,
    interpolation: Interpolation = Interpolation.LINEAR,
    y_position: float = 0.5,
    x_position: float = 0.5,
) -> (SingleProfile, SingleProfile):
    stream = io.BytesIO()
    simulator.generate_dicom(stream)
    stream.seek(0)
    img = load(stream)
    y_pixel = int(round(simulator.shape[0] * y_position))
    x_pixel = int(round(simulator.shape[1] * x_position))
    inplane_profile = SingleProfile(
        img[:, x_pixel],
        dpmm=img.dpmm,
        interpolation=interpolation,
        normalization_method=Normalization.NONE,
    )
    cross_profile = SingleProfile(
        img[y_pixel, :],
        dpmm=img.dpmm,
        interpolation=interpolation,
        normalization_method=Normalization.NONE,
    )
    return inplane_profile, cross_profile


class TestFilteredFieldLayer(TestCase):
    def test_50x50_1000sid_centered(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(FilteredFieldLayer(field_size_mm=(50, 50)))
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the center
            self.assertAlmostEqual(centers[0].x, img.center.x, delta=1)
            self.assertAlmostEqual(centers[0].y, img.center.y, delta=1)

    def test_50x50_1000sid_offset(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(
                FilteredFieldLayer(field_size_mm=(50, 50), cax_offset_mm=(30, 50))
            )
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the offset
            self.assertAlmostEqual(
                centers[0].x, img.center.x + 50 / sim.pixel_size, delta=1
            )
            self.assertAlmostEqual(
                centers[0].y, img.center.y + 30 / sim.pixel_size, delta=1
            )

    def test_50x50_1500sid_centered(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(FilteredFieldLayer(field_size_mm=(50, 50)))
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the center
            self.assertAlmostEqual(centers[0].x, img.center.x, delta=1)
            self.assertAlmostEqual(centers[0].y, img.center.y, delta=1)

    def test_50x50_1500sid_offset(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(
                FilteredFieldLayer(field_size_mm=(50, 50), cax_offset_mm=(30, 50))
            )
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the offset
            self.assertAlmostEqual(
                centers[0].x, img.center.x + 50 * 1.5 / sim.pixel_size, delta=1
            )
            self.assertAlmostEqual(
                centers[0].y, img.center.y + 30 * 1.5 / sim.pixel_size, delta=1
            )


class TestPerfectFieldLayer(TestCase):
    def test_50x50_1000sid_centered(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(50, 50)))
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the center
            self.assertAlmostEqual(centers[0].x, img.center.x, delta=1)
            self.assertAlmostEqual(centers[0].y, img.center.y, delta=1)

    def test_50x50_1000sid_offset(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(
                PerfectFieldLayer(field_size_mm=(50, 50), cax_offset_mm=(30, 50))
            )
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the offset
            self.assertAlmostEqual(
                centers[0].x, img.center.x + 50 / sim.pixel_size, delta=1
            )
            self.assertAlmostEqual(
                centers[0].y, img.center.y + 30 / sim.pixel_size, delta=1
            )

    def test_50x50_1500sid_centered(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(50, 50)))
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the center
            self.assertAlmostEqual(centers[0].x, img.center.x, delta=1)
            self.assertAlmostEqual(centers[0].y, img.center.y, delta=1)

    def test_50x50_1500sid_offset(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(
                PerfectFieldLayer(field_size_mm=(50, 50), cax_offset_mm=(30, 50))
            )
            as1200_ds = as1200.as_dicom()
            img = DicomImage.from_dataset(as1200_ds)
            centers = img.compute(GlobalFieldLocator(max_number=1))
            # test we're at the offset
            self.assertAlmostEqual(
                centers[0].x, img.center.x + 50 * 1.5 / sim.pixel_size, delta=1
            )
            self.assertAlmostEqual(
                centers[0].y, img.center.y + 30 * 1.5 / sim.pixel_size, delta=1
            )

    def test_10x10_100sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(100, 100)))
            as1200.add_layer(GaussianFilterLayer(sigma_mm=0.5))
            # no interpolation
            inplane_profile, cross_profile = profiles_from_simulator(
                as1200, interpolation=Interpolation.NONE
            )
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            # linear interp
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            # spline interp
            inplane_profile, cross_profile = profiles_from_simulator(
                as1200, interpolation=Interpolation.SPLINE
            )
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )

    def test_10x10_150sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(100, 100)))
            as1200.add_layer(GaussianFilterLayer(sigma_mm=0.5))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                150,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                150,
                delta=sim.pixel_size * 0.6,
            )

    def test_offset_150sid(self):
        as1200 = AS1200Image(sid=1500)
        as1200.add_layer(
            PerfectFieldLayer(field_size_mm=(100, 100), cax_offset_mm=(20, 30))
        )
        as1200.add_layer(GaussianFilterLayer(sigma_mm=0.5))
        y_position = 0.5 + (10 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[0]
        )
        x_position = 0.5 + (20 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[1]
        )
        inplane_profile, cross_profile = profiles_from_simulator(
            as1200, y_position=y_position, x_position=x_position
        )
        self.assertAlmostEqual(
            inplane_profile.fwxm_data()["width (exact) mm"], 150, delta=1
        )
        self.assertAlmostEqual(
            cross_profile.fwxm_data()["width (exact) mm"], 150, delta=1
        )


class TestSlopeLayer(TestCase):
    def test_slope_x(self):
        as1200 = AS1200Image(sid=1000)
        # create field that will cover the whole image; want a uniform field
        as1200.add_layer(PerfectFieldLayer(field_size_mm=(400, 400), alpha=0.6))
        as1200.add_layer(SlopeLayer(slope_x=0.1, slope_y=0))
        # test that the left is less than the right edge
        left = as1200.image[:, 100].max()
        right = as1200.image[:, -100].max()
        self.assertLess(left, right)

    def test_negative_slope_x(self):
        as1200 = AS1200Image(sid=1000)
        # create field that will cover the whole image; want a uniform field
        as1200.add_layer(PerfectFieldLayer(field_size_mm=(400, 400), alpha=0.6))
        as1200.add_layer(SlopeLayer(slope_x=-0.1, slope_y=0))
        # test that the left is greater than the right edge
        left = as1200.image[:, 100].max()
        right = as1200.image[:, -100].max()
        self.assertGreater(left, right)

    def test_slope_y(self):
        as1200 = AS1200Image(sid=1000)
        # create field that will cover the whole image; want a uniform field
        as1200.add_layer(PerfectFieldLayer(field_size_mm=(400, 400), alpha=0.6))
        as1200.add_layer(SlopeLayer(slope_x=0, slope_y=0.1))
        # test that the top is less than the bottom edge
        top = as1200.image[100, :].max()
        bottom = as1200.image[-100, :].max()
        self.assertLess(top, bottom)


class TestPerfectConeLayer(TestCase):
    def test_alpha(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(PerfectConeLayer(cone_size_mm=15, alpha=0.5))
        self.assertAlmostEqual(
            as1200.image.max(), np.iinfo(np.uint16).max * 0.5, delta=1
        )

    def test_15mm_1000sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(PerfectConeLayer(cone_size_mm=15))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"], 15, delta=1
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"], 15, delta=1
            )
            self.assertAlmostEqual(
                as1200.image.max(), np.iinfo(np.uint16).max
            )  # default alpha is 1, thus max value

    def test_15mm_1500sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(PerfectConeLayer(cone_size_mm=15))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"], 15 * 1.5, delta=1
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"], 15 * 1.5, delta=1
            )
            self.assertAlmostEqual(
                as1200.image.max(), np.iinfo(np.uint16).max
            )  # default alpha is 1, thus max value

    def test_offset_1000sid(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(PerfectConeLayer(cone_size_mm=15, cax_offset_mm=(10, 20)))
        y_position = 0.5 + (10 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[0]
        )
        x_position = 0.5 + (20 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[1]
        )
        inplane_profile, cross_profile = profiles_from_simulator(
            as1200, y_position=y_position, x_position=x_position
        )
        self.assertAlmostEqual(
            inplane_profile.fwxm_data()["width (exact) mm"], 15, delta=0.3
        )
        self.assertAlmostEqual(
            cross_profile.fwxm_data()["width (exact) mm"], 15, delta=0.3
        )

    def test_out_20mm(self):
        as1200 = AS1200Image(sid=1000)
        # 10mm out, meaning down
        as1200.add_layer(PerfectConeLayer(cone_size_mm=15, cax_offset_mm=(10, 0)))
        ds = as1200.as_dicom()
        img = DicomImage.from_dataset(ds)
        centers = img.compute(GlobalFieldLocator(max_number=1))
        # y will be 10mm down from center (positive)
        y_position = img.center.y + 10 / as1200.pixel_size
        self.assertAlmostEqual(centers[0].y, y_position, delta=0.1)

    def test_offset_1500sid(self):
        as1200 = AS1200Image(sid=1500)
        as1200.add_layer(PerfectConeLayer(cone_size_mm=15, cax_offset_mm=(10, 20)))
        y_position = 0.5 + (10 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[0]
        )
        x_position = 0.5 + (20 * as1200.mag_factor) / (
            as1200.pixel_size * as1200.shape[1]
        )
        inplane_profile, cross_profile = profiles_from_simulator(
            as1200, y_position=y_position, x_position=x_position
        )
        self.assertAlmostEqual(
            inplane_profile.fwxm_data()["width (exact) mm"], 15 * 1.5, delta=0.3
        )
        self.assertAlmostEqual(
            cross_profile.fwxm_data()["width (exact) mm"], 15 * 1.5, delta=0.3
        )


class TestPerfectBBLayer(TestCase):
    def test_10mm_100sid(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(
            ConstantLayer(constant=1.0)
        )  # simulate huge field for easier analysis later on
        as1200.add_layer(PerfectBBLayer(bb_size_mm=10))
        stream = io.BytesIO()
        as1200.generate_dicom(stream)
        stream.seek(0)
        img = load(stream)
        img.invert()  # we invert so the BB looks like a profile, not a dip
        inplane_profile = SingleProfile(img[:, int(as1200.shape[1] / 2)], dpmm=img.dpmm)
        cross_profile = SingleProfile(img[int(as1200.shape[0] / 2), :], dpmm=img.dpmm)
        self.assertAlmostEqual(
            inplane_profile.fwxm_data()["width (exact) mm"],
            10,
            delta=as1200.pixel_size * 0.6,
        )
        self.assertAlmostEqual(
            cross_profile.fwxm_data()["width (exact) mm"],
            10,
            delta=as1200.pixel_size * 0.6,
        )

    def test_10mm_150sid(self):
        as1200 = AS1200Image(sid=1500)
        as1200.add_layer(
            ConstantLayer(constant=1.0)
        )  # simulate huge field for easier analysis later on
        as1200.add_layer(PerfectBBLayer(bb_size_mm=10))
        stream = io.BytesIO()
        as1200.generate_dicom(stream)
        stream.seek(0)
        img = load(stream)
        img.invert()  # we invert so the BB looks like a profile, not a dip
        inplane_profile = SingleProfile(img[:, int(as1200.shape[1] / 2)], dpmm=img.dpmm)
        cross_profile = SingleProfile(img[int(as1200.shape[0] / 2), :], dpmm=img.dpmm)
        self.assertAlmostEqual(
            inplane_profile.fwxm_data()["width (exact) mm"], 15, delta=1
        )
        self.assertAlmostEqual(
            cross_profile.fwxm_data()["width (exact) mm"], 15, delta=1
        )


class TestFFFLayer(TestCase):
    def test_10x10_100sid(self):
        for sim in (AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(FilterFreeFieldLayer(field_size_mm=(100, 100)))
            # no interpolation
            inplane_profile, cross_profile = profiles_from_simulator(
                as1200, interpolation=Interpolation.NONE
            )
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            # linear interp
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            # spline interp
            as1200.add_layer(
                GaussianFilterLayer(sigma_mm=0.2)
            )  # spline causes ringing artifacts for ultra-sharp gradients, this is also more realistic anyway
            inplane_profile, cross_profile = profiles_from_simulator(
                as1200, interpolation=Interpolation.SPLINE
            )
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                100,
                delta=sim.pixel_size * 0.6,
            )

    def test_10x10_150sid(self):
        for sim in (AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(FilterFreeFieldLayer(field_size_mm=(150, 150)))
            inplane_profile, cross_profile = profiles_from_simulator(
                as1200, interpolation=Interpolation.NONE
            )
            self.assertAlmostEqual(
                inplane_profile.fwxm_data()["width (exact) mm"],
                150,
                delta=sim.pixel_size * 0.6,
            )
            self.assertAlmostEqual(
                cross_profile.fwxm_data()["width (exact) mm"],
                150,
                delta=sim.pixel_size * 0.6,
            )


class TestRandomNoise(TestCase):
    def test_mean_doesnt_change(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35000))
        as1200.add_layer(RandomNoiseLayer(mean=0, sigma=0.001))
        self.assertAlmostEqual(as1200.image.mean(), 35000, delta=1)

    def test_std(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35000))
        as1200.add_layer(RandomNoiseLayer(mean=0, sigma=0.003))
        std = np.iinfo(np.uint16).max * 0.003
        self.assertAlmostEqual(as1200.image.std(), std, delta=1)


class TestConstantLayer(TestCase):
    def test_constant(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35))
        self.assertAlmostEqual(as1200.image.max(), 35)
        self.assertAlmostEqual(as1200.image.min(), 35)

    def test_two_constants(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35))
        as1200.add_layer(ConstantLayer(constant=11))
        self.assertAlmostEqual(as1200.image.max(), 46)
        self.assertAlmostEqual(as1200.image.min(), 46)

    def test_constant_wont_flip_bits_over(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35))
        as1200.add_layer(ConstantLayer(constant=777777777))
        self.assertAlmostEqual(as1200.image.max(), np.iinfo(np.uint16).max)
        self.assertAlmostEqual(as1200.image.min(), np.iinfo(np.uint16).max)

    def test_constant_wont_flip_bits_under(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=35))
        as1200.add_layer(ConstantLayer(constant=-777777777))
        self.assertAlmostEqual(as1200.image.max(), np.iinfo(np.uint16).min)
        self.assertAlmostEqual(as1200.image.min(), np.iinfo(np.uint16).min)


class NOOPLayer(Layer):
    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        return image


class SimulatorTestMixin:
    simulator: Simulator
    pixel_size: float
    shape: (int, int)
    mag_factor = 1.5

    def test_pixel_size(self):
        self.assertEqual(self.simulator().pixel_size, self.pixel_size)

    def test_shape(self):
        self.assertEqual(self.simulator().shape, self.shape)

    def test_image(self):
        sim = self.simulator()
        self.assertEqual(sim.image.shape, self.shape)
        self.assertEqual(sim.image.dtype, np.uint16)

    def test_mag_factor(self):
        self.assertEqual(self.simulator().mag_factor, self.mag_factor)
        ssd1000 = self.simulator(sid=1000)
        self.assertEqual(ssd1000.mag_factor, 1)

    def test_noop_layer_doesnt_change_image(self):
        sim = self.simulator()
        orig_img = copy.deepcopy(sim.image)
        sim.add_layer(NOOPLayer())
        numpy.testing.assert_array_equal(sim.image, orig_img)

    def test_save_dicom(self):
        sim = self.simulator()
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            sim.generate_dicom(tf.name, gantry_angle=12, coll_angle=33, table_angle=5)
            # shouldn't raise
            ds = pydicom.dcmread(tf.name)
        self.assertEqual(ds.pixel_array.shape, self.shape)
        self.assertEqual(ds.GantryAngle, 12)
        self.assertEqual(ds.BeamLimitingDeviceAngle, 33)
        self.assertEqual(ds.PatientSupportAngle, 5)


class TestAS500(SimulatorTestMixin, TestCase):
    simulator = AS500Image
    pixel_size = 0.78125
    shape = (384, 512)


class TestAS1000(SimulatorTestMixin, TestCase):
    simulator = AS1000Image
    pixel_size = 0.390625
    shape = (768, 1024)


class TestAS1200(SimulatorTestMixin, TestCase):
    simulator = AS1200Image
    pixel_size = 0.336
    shape = (1280, 1280)


class CustomSim(Simulator):
    pixel_size = 0.123
    shape = (500, 700)


class TestCustomSimulator(SimulatorTestMixin, TestCase):
    simulator = CustomSim
    pixel_size = 0.123
    shape = (500, 700)

    def test_save_dicom(self):
        sim = self.simulator()
        with self.assertRaises(NotImplementedError):
            with tempfile.NamedTemporaryFile(delete=False) as tf:
                sim.generate_dicom(
                    tf.name, gantry_angle=12, coll_angle=33, table_angle=5
                )
