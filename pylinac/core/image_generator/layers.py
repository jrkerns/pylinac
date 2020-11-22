from abc import abstractmethod, ABC

import matplotlib.pyplot as plt
import numpy as np
from skimage import draw, filters


def clip_add(image1: np.ndarray, image2: np.ndarray, dtype: np.dtype = np.uint16):
    """Clip the image to the dtype extrema. Otherwise the bits will flip."""
    return np.clip(image1 + image2, np.iinfo(dtype).min, np.iinfo(dtype).max).astype(dtype)


def even_round(num: float) -> int:
    """Return an even number"""
    num = int(round(num))
    return num + num % 2


def gaussian2d(mx, my, height, center_x, center_y, width_x, width_y, constant=0):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return height * np.exp(-(((center_x - mx) / width_x) ** 2 + ((center_y - my) / width_y) ** 2) / 2) + constant


class Layer(ABC):

    @abstractmethod
    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        """Apply the layer. Takes a 2D array and pixel size value in and returns a modified array."""
        pass


class PerfectConeLayer(Layer):
    """A cone without flattening filter effects"""

    def __init__(self, cone_size_mm=10, cax_offset_mm=(0, 0), alpha=1.0):
        self.cone_size_mm = cone_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        image, _, _ = self._create_perfect_field(image, pixel_size)
        return image

    def _create_perfect_field(self, image, pixel_size):
        cone_size_pix = (self.cone_size_mm / 2) / pixel_size
        cax_offset_pix = [x / pixel_size + (shape / 2 - 0.5) for x, shape in zip(self.cax_offset_mm, image.shape)]
        rr, cc = draw.disk(cax_offset_pix, cone_size_pix, shape=image.shape)
        rr = np.round(rr).astype(np.int)
        cc = np.round(cc).astype(np.int)
        temp_array = np.zeros(image.shape)
        temp_array[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        image = clip_add(image, temp_array)
        return image, rr, cc


class FilterFreeConeLayer(PerfectConeLayer):
    """A cone with flattening filter effects."""

    def __init__(self, cone_size_mm=10, cax_offset_mm=(0, 0), alpha=1.0, filter_magnitude=0.4, filter_sigma_mm=80,
                 penumbra_mm=2):
        super().__init__(cone_size_mm, cax_offset_mm, alpha)
        self.filter_magnitude = filter_magnitude
        self.filter_sigma_mm = filter_sigma_mm
        self.penumbra_mm = penumbra_mm

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        image, rr, cc = self._create_perfect_field(image, pixel_size)
        # add filter effect
        n = gaussian2d(rr, cc, self.filter_magnitude * np.iinfo(image.dtype).max, image.shape[0] / 2,
                       image.shape[1] / 2,
                       self.filter_sigma_mm / pixel_size, self.filter_sigma_mm / pixel_size,
                       constant=-self.filter_magnitude * np.iinfo(image.dtype).max)
        image[rr, cc] += n.astype(image.dtype)
        return image


class PerfectFieldLayer(Layer):
    """A square field without flattening filter effects"""

    def __init__(self, field_size_mm=(10, 10), cax_offset_mm=(0, 0), alpha=1.0):
        self.field_size_mm = field_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def _create_perfect_field(self, image, pixel_size):
        field_size_pix = [even_round(f / pixel_size) for f in self.field_size_mm]
        field_start = [x / pixel_size + (shape / 2) - field_size / 2 for x, shape, field_size in
                       zip(self.cax_offset_mm, image.shape, field_size_pix)]
        field_end = [x / pixel_size + (shape / 2) + field_size / 2 - 1 for x, shape, field_size in
                     zip(self.cax_offset_mm, image.shape, field_size_pix)]
        # -1 due to skimage implementation of [start:(end+1)]
        rr, cc = draw.rectangle(field_start, end=field_end, shape=image.shape)
        rr = np.round(rr).astype(np.int)
        cc = np.round(cc).astype(np.int)
        temp_array = np.zeros(image.shape)
        temp_array[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        image = clip_add(image, temp_array)
        return image, rr, cc

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        image, _, _ = self._create_perfect_field(image, pixel_size)
        return image


class FilteredFieldLayer(PerfectFieldLayer):
    """A square field with flattening filter effects"""

    def __init__(self, field_size_mm=(10, 10), cax_offset_mm=(0, 0), alpha=1.0, gaussian_height=0.03,
                 gaussian_sigma_mm=32):
        super().__init__(field_size_mm, cax_offset_mm, alpha)
        self.gaussian_height = gaussian_height
        self.gaussian_sigma_mm = gaussian_sigma_mm

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        image, rr, cc = self._create_perfect_field(image, pixel_size)
        # add filter effect
        height = -self.gaussian_height * np.iinfo(image.dtype).max
        width = self.gaussian_sigma_mm / pixel_size
        horns = gaussian2d(rr, cc, height=height, center_x=image.shape[0] / 2,
                           center_y=image.shape[1] / 2,
                           width_x=width, width_y=width)
        image[rr, cc] += horns.astype(image.dtype)
        return image


class FilterFreeFieldLayer(FilteredFieldLayer):
    """A square field with flattening filter free (FFF) effects"""

    def __init__(self, field_size_mm=(10, 10), cax_offset_mm=(0, 0), alpha=1.0, gaussian_height=0.4,
                 gaussian_sigma_mm=80):
        super().__init__(field_size_mm, cax_offset_mm, alpha, gaussian_height, gaussian_sigma_mm)

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        image, rr, cc = self._create_perfect_field(image, pixel_size)
        # add filter effect
        n = gaussian2d(rr, cc, self.gaussian_height * np.iinfo(image.dtype).max, image.shape[0] / 2, image.shape[1] / 2,
                       self.gaussian_sigma_mm / pixel_size, self.gaussian_sigma_mm / pixel_size,
                       constant=-self.gaussian_height * np.iinfo(image.dtype).max)
        image[rr, cc] += n.astype(image.dtype)
        return image


class PerfectBBLayer(PerfectConeLayer):
    """A BB-like layer. Like a cone, but with lower alpha (i.e. higher opacity)"""

    def __init__(self, bb_size_mm=5, cax_offset_mm=(0, 0), alpha=-0.5):
        super().__init__(cone_size_mm=bb_size_mm, cax_offset_mm=cax_offset_mm, alpha=alpha)


class GaussianFilterLayer(Layer):
    """A Gaussian filter. Simulates the effects of scatter on the field"""

    def __init__(self, sigma_mm=2):
        self.sigma_mm = sigma_mm

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        sigma_pix = self.sigma_mm / pixel_size
        return filters.gaussian(image, sigma_pix, preserve_range=True).astype(image.dtype)


class RandomNoiseLayer(Layer):
    """A salt and pepper noise, simulating dark current"""

    def __init__(self, mean=0.0, sigma=0.01):
        self.mean = mean
        self.sigma = sigma

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        noise = np.random.normal(self.mean, self.sigma, size=image.shape)
        return clip_add(image, noise, dtype=image.dtype)


class ConstantLayer(Layer):
    """A constant layer. Can be used to simulate scatter or background."""

    def __init__(self, constant):
        self.constant = constant

    def apply(self, image: np.ndarray, pixel_size: float) -> np.ndarray:
        return clip_add(image, self.constant, dtype=image.dtype)
