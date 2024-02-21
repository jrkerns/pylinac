import math
from unittest import TestCase

import numpy as np

from pylinac.core.nps import (
    average_power,
    max_frequency,
    noise_power_spectrum_1d,
    noise_power_spectrum_2d,
)

np.random.seed(123)  # noqa: NPY002


def generate_gaussian_noise_map(
    shape: (int, int), scale: int = 10, intensity: float = 0.5
):
    """Generate a Gaussian noise map with varying intensities (clumps).

    Parameters:
    - shape: Shape of the noise map (height, width).
    - scale: Scale of the Gaussian clumps.
    - intensity: Intensity of the noise.

    Returns:
    - Gaussian noise map as a NumPy array.
    """
    # Create low-resolution noise
    low_res_shape = (shape[0] // scale, shape[1] // scale)
    low_res_noise = np.random.normal(  # noqa: NPY002
        loc=0, scale=intensity, size=low_res_shape
    )

    # Upscale the noise to the original resolution
    noise_map = np.kron(low_res_noise, np.ones((scale, scale)))

    # Ensure the noise map matches the exact original shape, trimming if necessary
    noise_map = noise_map[: shape[0], : shape[1]]

    return noise_map


def apply_noise_clumps(image: np.ndarray, noise_map: np.ndarray) -> np.ndarray:
    """
    Apply a Gaussian noise map to an image.

    Parameters:
    - image: Original image as a NumPy array.
    - noise_map: Gaussian noise map.

    Returns:
    - Noisy image as a NumPy array.
    """
    noisy_image = image + noise_map
    dtype_min, dtype_max = np.iinfo(image.dtype).min, np.iinfo(image.dtype).max
    noisy_image = np.clip(noisy_image, dtype_min, dtype_max)

    return noisy_image


def generate_noisy_image(
    shape: (int, int),
    scale: int = 10,
    intensity: float = 0.5,
    dtype: np.dtype = np.uint16,
) -> np.ndarray:
    """
    Generate a noisy image with Gaussian clumps.

    Parameters:
    - image: Original image as a NumPy array.
    - scale: Scale of the Gaussian clumps.
    - intensity: Intensity of the noise.

    Returns:
    - Noisy image as a NumPy array.
    """
    latent = np.zeros(shape, dtype=dtype)
    noise_map = generate_gaussian_noise_map(latent.shape, scale, intensity)
    noisy_image = apply_noise_clumps(latent, noise_map)
    return noisy_image


class Test2DSpectrum(TestCase):
    def test_single_roi(self):
        roi = generate_noisy_image((300, 300), scale=30, intensity=500, dtype=np.uint16)
        nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[roi])
        self.assertEqual(nps2d.shape, roi.shape)

    def test_multiple_rois(self):
        roi1 = generate_noisy_image(
            (300, 300), scale=30, intensity=500, dtype=np.uint16
        )
        roi2 = generate_noisy_image(
            (300, 300), scale=10, intensity=100, dtype=np.uint16
        )
        nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[roi1, roi2])
        self.assertEqual(nps2d.shape, roi1.shape)

    def test_take_smallest_shape(self):
        roi1 = generate_noisy_image(
            (300, 300), scale=30, intensity=500, dtype=np.uint16
        )
        roi2 = generate_noisy_image(
            (200, 200), scale=10, intensity=100, dtype=np.uint16
        )
        nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[roi1, roi2])
        self.assertEqual(nps2d.shape, roi2.shape)


class Test1DSpectrum(TestCase):
    def test_1d_spectrum_uniform(self):
        nps2d = np.ones((300, 300))
        nps1d = noise_power_spectrum_1d(nps2d)
        self.assertAlmostEqual(nps1d[0], 1, delta=0.0001)

    def test_1d_spectrum(self):
        roi = generate_noisy_image((300, 300), scale=30, intensity=500, dtype=np.uint16)
        nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[roi])
        nps1d = noise_power_spectrum_1d(nps2d)
        # shape is same as diagonal distance from center to corner
        self.assertEqual(len(nps1d), math.ceil(300 * math.sqrt(2) / 2))


class TestAvgPower(TestCase):
    def setUp(self) -> None:
        self.roi = generate_noisy_image(
            (300, 300), scale=30, intensity=500, dtype=np.uint16
        )
        self.nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[self.roi])

    def test_avg_power(self):
        nps1d = noise_power_spectrum_1d(self.nps2d)
        avg_power = average_power(nps1d)
        self.assertAlmostEqual(avg_power, 0.0207, delta=0.005)

    def test_odd_roi_size_same_as_even(self):
        nps1d = noise_power_spectrum_1d(self.nps2d)
        avg_power_even = average_power(nps1d)
        nps2d_odd = noise_power_spectrum_2d(pixel_size=1, rois=[self.roi[:-1, :-1]])
        avg_power_odd = average_power(noise_power_spectrum_1d(nps2d_odd))
        self.assertAlmostEqual(avg_power_even, avg_power_odd, delta=0.0005)


class TestFrequency(TestCase):
    def test_frequency(self):
        roi = generate_noisy_image((300, 300), scale=30, intensity=500, dtype=np.uint16)
        nps2d = noise_power_spectrum_2d(pixel_size=1, rois=[roi])
        nps1d = noise_power_spectrum_1d(nps2d)
        f = max_frequency(nps1d)
        self.assertAlmostEqual(f, 0.0094, delta=0.0001)
