from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from skimage import draw, filters
from skimage.draw import polygon

from ..array_utils import geometric_center_idx


def clip_add(
    image1: np.ndarray, image2: np.ndarray, dtype: type[np.dtype] = np.uint16
) -> np.ndarray:
    """Clip the image to the dtype extrema. Otherwise, the bits will flip."""
    # convert to float first so we don't flip bits initially
    combined_img = image1.astype(float) + image2.astype(float)
    return np.clip(combined_img, np.iinfo(dtype).min, np.iinfo(dtype).max).astype(dtype)


def clip_multiply(
    image1: np.ndarray, image2: np.ndarray, dtype: type[np.dtype] = np.uint16
) -> np.ndarray:
    """Clip the image to the dtype extrema. Otherwise, the bits will flip."""
    # convert to float first so we don't flip bits initially
    combined_img = image1.astype(float) * image2.astype(float)
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
        rotation: float = 0,
    ):
        """
        Parameters
        ----------

        cone_size_mm
            Cone size in mm at the iso plane
        cax_offset_mm
            The offset in mm. (down, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        rotation: float
            The amount of rotation in degrees. When there is an offset, this acts like a couch kick.
        """
        self.cone_size_mm = cone_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha
        self.rotation = rotation

    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.ndarray:
        image, _, _ = self._create_perfect_field(image, pixel_size, mag_factor)
        return image

    def _create_perfect_field(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> (np.ndarray, ...):
        cone_size_pix = ((self.cone_size_mm / 2) / pixel_size) * mag_factor**2
        # we rotate the point around the center of the image
        offset_pix_y, offset_pix_x = rotate_point(
            x=self.cax_offset_mm[0] * mag_factor / pixel_size,
            y=self.cax_offset_mm[1] * mag_factor / pixel_size,
            angle=self.rotation,
        )
        # convert to pixels and shift to center
        cax_offset_pix = (
            offset_pix_y + (image.shape[0] / 2 - 0.5),
            offset_pix_x + (image.shape[1] / 2 - 0.5),
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
        rotation: float = 0,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane as (height, width)
        cax_offset_mm
            The offset in mm. (down, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        rotation: float
            The amount of rotation in degrees. This acts like a collimator rotation.
        """
        self.field_size_mm = field_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha
        self.rotation = rotation

    def _create_perfect_field(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> (np.ndarray, ...):
        field_size_pix = [
            even_round(f * mag_factor**2 / pixel_size) for f in self.field_size_mm
        ]
        cax_offset_pix_mag = [v * mag_factor / pixel_size for v in self.cax_offset_mm]
        field_center = [
            offset + (shape / 2) - 0.5
            for offset, shape in zip(cax_offset_pix_mag, image.shape)
        ]
        rr, cc = draw_rotated_rectangle(
            image.shape, center=field_center, extent=field_size_pix, angle=self.rotation
        )
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
        rotation: float = 0,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane (height, width)
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
        rotation: float
            The amount of rotation in degrees. This acts like a collimator rotation.
        """
        super().__init__(
            field_size_mm=field_size_mm,
            cax_offset_mm=cax_offset_mm,
            alpha=alpha,
            rotation=rotation,
        )
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
        # the horns are negative, so we don't have
        # to worry about clipping here
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
        rotation: float = 0,
    ):
        """
        Parameters
        ----------

        field_size_mm
            Field size in mm at the iso plane (height, width).
        cax_offset_mm
            The offset in mm. (out, right)
        alpha
            The intensity of the layer. 1 is full saturation/radiation. 0 is none.
        gaussian_height
            The magnitude of the CAX peak. Larger values result in "pointier" fields.
        gaussian_sigma_mm
            Proportional to the width of the CAX peak. Larger values produce wider curves.
        rotation: float
            The amount of rotation in degrees. This acts like a collimator rotation.
        """
        super().__init__(
            field_size_mm,
            cax_offset_mm,
            alpha,
            gaussian_height,
            gaussian_sigma_mm,
            rotation=rotation,
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
        # the horns are negative, so we don't have
        # to worry about clipping here
        image[rr, cc] += n.astype(image.dtype)
        return image


class PerfectBBLayer(PerfectConeLayer):
    """A BB-like layer. Like a cone, but with lower alpha (i.e. higher opacity)"""

    def __init__(
        self,
        bb_size_mm: float = 5,
        cax_offset_mm: (float, float) = (0, 0),
        alpha: float = -0.5,
        rotation: float = 0,
    ):
        super().__init__(
            cone_size_mm=bb_size_mm,
            cax_offset_mm=cax_offset_mm,
            alpha=alpha,
            rotation=rotation,
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


class SlopeLayer(Layer):
    """Adds a slope in both directions of the image. Usually used for simulating asymmetry or a-flatness.

    Parameters
    ----------
    slope_x : float
        The slope in the x-direction (left/right). If positive, will increase the right side.
        The value is multiplicative to the current state of the image. E.g. a value of 0.1 will increase the right side by 10% and 0% on the left.
    slope_y : float
        The slope in the y-direction (up/down). If positive, will increase the bottom side.
    """

    def __init__(self, slope_x: float, slope_y: float):
        self.slope_x = slope_x
        self.slope_y = slope_y

    def apply(
        self, image: np.ndarray, pixel_size: float, mag_factor: float
    ) -> np.ndarray:
        nrows, ncols = image.shape
        # Create the row and column scaling factors
        y_scaling = (1 + self.slope_y * np.arange(nrows) / nrows).reshape(-1, 1)
        x_scaling = (1 + self.slope_x * np.arange(ncols) / ncols).reshape(1, -1)

        y_scaled = clip_multiply(image, y_scaling)
        xy_scaled = clip_multiply(y_scaled, x_scaling)
        return xy_scaled


def rotate_point(x: float, y: float, angle: float) -> (float, float):
    """
    Rotate a point (px, py) about the origin by a given angle in degrees.

    Parameters
    ----------
    x : float
        The x-coordinate of the point to rotate.
    y : float
        The y-coordinate of the point to rotate.
    angle : float
        The angle of rotation in degrees.

    Returns
    -------
    px_rotated, py_rotated : tuple
        The rotated point.

    Notes
    -----
    ChatGPT-generated 🥰
    """
    # Convert angle from degrees to radians
    theta = np.radians(angle)

    # Apply rotation
    px_rotated = x * np.cos(theta) - y * np.sin(theta)
    py_rotated = x * np.sin(theta) + y * np.cos(theta)

    return px_rotated, py_rotated


def draw_rotated_rectangle(
    shape: tuple[int, int],
    center: list[float, float],
    extent: list[int, int],
    angle: float,
):
    """Generates the coordinate points for a rectangle rotated about the image center.

    Parameters
    ----------
    shape : tuple
        The shape of the image.
    center : list
        The center of the rectangle.
    extent : list
        The height and width of the rectangle.
    angle : float
        The angle of rotation in degrees.

    Returns
    -------
    rr, cc : tuple
        The row and column coordinates of the rectangle.

    Notes
    -----
    ChatGPT-generated 🥰
    """
    # Calculate rectangle coordinates before rotation
    # follows row,col convention of numpy; y = row, x = col
    x0 = center[1] - extent[1] / 2
    x1 = center[1] + extent[1] / 2
    y0 = center[0] - extent[0] / 2
    y1 = center[0] + extent[0] / 2

    rect_coords = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])

    # Calculate rotation matrix
    theta = np.radians(angle)
    c, s = np.cos(theta), np.sin(theta)
    rotation = np.array([[c, -s], [s, c]])

    # Rotate rectangle coordinates
    rotated_coords = np.dot(rect_coords - np.array(center), rotation) + np.array(center)

    # Draw rotated rectangle
    rr, cc = polygon(rotated_coords[:, 1], rotated_coords[:, 0], shape)
    return rr, cc
