from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from skimage import draw, filters

from ..array_utils import geometric_center_idx


def clip_add(
    image1: np.ndarray, image2: np.ndarray, dtype: type[np.dtype] = np.uint16
) -> np.ndarray:
    """Clip the image to the dtype extrema. Otherwise, the bits will flip."""
    # convert to float first so we don't flip bits initially
    combined_img = image1.astype(float) + image2.astype(float)
    return np.clip(combined_img, np.iinfo(dtype).min, np.iinfo(dtype).max).astype(dtype)


def even_round(num: float) -> int:
    """Return an even number"""
    num = int(round(num))
    return num + num % 2


def gaussian2d(
    mx: float,
    my: float,
    height: float,
    center_x: float,
    center_y: float,
    width_x: float,
    width_y: float,
    constant: float = 0,
) -> np.ndarray:
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return (
        height
        * np.exp(
            -(((center_x - mx) / width_x) ** 2 + ((center_y - my) / width_y) ** 2) / 2
        )
        + constant
    )


class Layer(ABC):
    """Abstract base for layers"""

    @abstractmethod
    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.ndarray:
        """Apply the layer. Takes a 2D array and pixel size value in and returns a modified array."""
        pass


class PerfectConeLayer(Layer):
    """A cone without flattening filter effects"""

    def __init__(
        self,
        cone_size_mm: float = 10,
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = 1.0,
    ):
        """
        Parameters
        ----------

        cone_size_mm
            Cone size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        """
        self.cone_size_mm = cone_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.ndarray:
        image, _, _ = self._create_perfect_field(image, pixel_size, mag_factor)
        return image

    def _create_perfect_field(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> (np.ndarray, ...):
        cone_size_pix = ((self.cone_size_mm / 2) / pixel_size) * mag_factor**2
        cax_offset_pix = tuple(
            x * mag_factor / pixel_size + (shape / 2 - 0.5)
            for x, shape in zip(self.cax_offset_mm, image.shape)
        )
        rr, cc = draw.disk(cax_offset_pix, cone_size_pix, shape=image.shape)
        rr = np.round(rr).astype(int)
        cc = np.round(cc).astype(int)
        temp_array = np.zeros(image.shape)
        temp_array[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        image = clip_add(image, temp_array)
        return image, rr, cc


class FilterFreeConeLayer(PerfectConeLayer):
    """A cone with flattening filter effects."""

    def __init__(
        self,
        cone_size_mm: float = 10,
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = 1.0,
        filter_magnitude: float = 0.4,
        filter_sigma_mm: float = 80,
    ):
        """
        Parameters
        ----------

        cone_size_mm
            Cone size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        filter_magnitude
            The magnitude of the CAX peak. Larger values result in "pointier" fields.
        filter_sigma_mm
            Proportional to the width of the CAX peak. Larger values produce wider curves.
        """
        super().__init__(cone_size_mm, cax_offset_mm, alpha)
        self.filter_magnitude = filter_magnitude
        self.filter_sigma_mm = filter_sigma_mm

    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.ndarray:
        image, rr, cc = self._create_perfect_field(image, pixel_size, mag_factor)
        # add filter effect
        center_x = geometric_center_idx(image[:, 0])
        center_y = geometric_center_idx(image[0, :])
        n = gaussian2d(
            rr,
            cc,
            self.filter_magnitude * np.iinfo(image.dtype).max,
            center_x,
            center_y,
            self.filter_sigma_mm / pixel_size,
            self.filter_sigma_mm / pixel_size,
            constant=-self.filter_magnitude * np.iinfo(image.dtype).max,
        )
        image[rr, cc] += n.astype(image.dtype)
        return image


class PerfectFieldLayer(Layer):
    """A square field without flattening filter effects"""

    def __init__(
        self,
        field_size_mm: (float, float) = (10, 10),
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = 1.0,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        """
        self.field_size_mm = field_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def _create_perfect_field(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> (np.ndarray, ...):
        field_size_pix = [
            even_round(f * mag_factor**2 / pixel_size) for f in self.field_size_mm
        ]
        cax_offset_mm_mag = [v * mag_factor for v in self.cax_offset_mm]
        field_start = [
            x / pixel_size + (shape / 2) - field_size / 2
            for x, shape, field_size in zip(
                cax_offset_mm_mag, image.shape, field_size_pix
            )
        ]
        field_end = [
            x / pixel_size + (shape / 2) + field_size / 2 - 1
            for x, shape, field_size in zip(
                cax_offset_mm_mag, image.shape, field_size_pix
            )
        ]
        # -1 due to skimage implementation of [start:(end+1)]
        rr, cc = draw.rectangle(field_start, end=field_end, shape=image.shape)
        rr = np.round(rr).astype(int)
        cc = np.round(cc).astype(int)
        temp_array = np.zeros(image.shape)
        temp_array[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        image = clip_add(image, temp_array)
        return image, rr, cc

    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.array:
        image, _, _ = self._create_perfect_field(image, pixel_size, mag_factor)
        return image


class FilteredFieldLayer(PerfectFieldLayer):
    """A square field with flattening filter effects"""

    def __init__(
        self,
        field_size_mm: (float, float) = (10, 10),
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = 1.0,
        gaussian_height: float = 0.03,
        gaussian_sigma_mm: float = 32,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        gaussian_height
            The intensity of the "horns", or more accurately, the CAX dip. Proportional to the max value allowed for the data type.
            Increase to make the horns more prominent.
        gaussian_sigma_mm
            The width of the "horns". A.k.a. the CAX dip width. Increase to create a wider
            horn effect.
        """
        super().__init__(field_size_mm, cax_offset_mm, alpha)
        self.gaussian_height = gaussian_height
        self.gaussian_sigma_mm = gaussian_sigma_mm

    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        image, rr, cc = self._create_perfect_field(image, pixel_size, mag_factor)
        # add filter effect
        height = -self.gaussian_height * np.iinfo(image.dtype).max
        width = self.gaussian_sigma_mm / pixel_size
        center_x = geometric_center_idx(image[:, 0])
        center_y = geometric_center_idx(image[0, :])
        horns = gaussian2d(
            rr,
            cc,
            height=height,
            center_x=center_x,
            center_y=center_y,
            width_x=width,
            width_y=width,
        )
        image[rr, cc] += horns.astype(image.dtype)
        return image


class FilterFreeFieldLayer(FilteredFieldLayer):
    """A square field with flattening filter free (FFF) effects"""

    def __init__(
        self,
        field_size_mm: (float, float) = (10, 10),
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = 1.0,
        gaussian_height: float = 0.4,
        gaussian_sigma_mm: float = 80,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        gaussian_height
            The magnitude of the CAX peak. Larger values result in "pointier" fields.
        gaussian_sigma_mm
            Proportional to the width of the CAX peak. Larger values produce wider curves.
        """
        super().__init__(
            field_size_mm, cax_offset_mm, alpha, gaussian_height, gaussian_sigma_mm
        )

    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        image, rr, cc = self._create_perfect_field(image, pixel_size, mag_factor)
        # add filter effect
        center_x = geometric_center_idx(image[:, 0])
        center_y = geometric_center_idx(image[0, :])
        n = gaussian2d(
            rr,
            cc,
            self.gaussian_height * np.iinfo(image.dtype).max,
            center_x,
            center_y,
            self.gaussian_sigma_mm / pixel_size,
            self.gaussian_sigma_mm / pixel_size,
            constant=-self.gaussian_height * np.iinfo(image.dtype).max,
        )
        image[rr, cc] += n.astype(image.dtype)
        return image


class PerfectBBLayer(PerfectConeLayer):
    """A BB-like layer. Like a cone, but with lower alpha (i.e. higher opacity)"""

    def __init__(
        self,
        bb_size_mm: float = 5,
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = -0.5,
    ):
        super().__init__(
            cone_size_mm=bb_size_mm, cax_offset_mm=cax_offset_mm, alpha=alpha
        )


class GaussianFilterLayer(Layer):
    """A Gaussian filter. Simulates the effects of scatter on the field"""

    def __init__(self, sigma_mm: float = 2):
        self.sigma_mm = sigma_mm

    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        sigma_pix = self.sigma_mm / pixel_size
        return filters.gaussian(image, sigma_pix, preserve_range=True).astype(
            image.dtype
        )


class RandomNoiseLayer(Layer):
    """A salt and pepper noise, simulating dark current"""

    def __init__(self, mean: float = 0.0, sigma: float = 0.001):
        self.mean = mean
        self.sigma = sigma

    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        normalized_sigma = self.sigma * np.iinfo(image.dtype).max
        rng = np.random.default_rng()
        noise = rng.normal(self.mean, normalized_sigma, size=image.shape)
        return clip_add(image, noise, dtype=image.dtype)


class ConstantLayer(Layer):
    """A constant layer. Can be used to simulate scatter or background."""

    def __init__(self, constant: float):
        self.constant = constant

    def apply(self, image: np.array, pixel_size: float, mag_factor: float) -> np.array:
        constant_img = np.full(image.shape, fill_value=self.constant)
        return clip_add(image, constant_img, dtype=image.dtype)
