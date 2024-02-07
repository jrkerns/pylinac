"""This module holds classes for image loading and manipulation."""
from __future__ import annotations

import copy
import io
import json
import math
import os
import os.path as osp
import re
import warnings
from collections import Counter
from datetime import datetime
from functools import cached_property
from io import BufferedReader, BytesIO
from pathlib import Path
from typing import Any, BinaryIO, Iterable, Sequence, Union

import argue
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import scipy.ndimage as spf
from PIL import Image as pImage
from PIL.PngImagePlugin import PngInfo
from PIL.TiffTags import TAGS
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.errors import InvalidDicomError
from pydicom.uid import UID, generate_uid
from scipy import ndimage
from skimage.draw import disk
from skimage.transform import rotate

from ..metrics.image import MetricBase
from ..settings import PATH_TRUNCATION_LENGTH, get_dicom_cmap
from .array_utils import (
    bit_invert,
    convert_to_dtype,
    filter,
    get_dtype_info,
    ground,
    invert,
    normalize,
)
from .geometry import Point
from .io import (
    TemporaryZipDirectory,
    get_url,
    is_dicom_image,
    retrieve_dicom_file,
    retrieve_filenames,
)
from .profile import stretch as stretcharray
from .scale import wrap360
from .utilities import decode_binary, is_close, simple_round

ARRAY = "Array"
DICOM = "DICOM"
IMAGE = "Image"

FILE_TYPE = "file"
STREAM_TYPE = "stream"

XIM_PROP_INT = 0
XIM_PROP_DOUBLE = 1
XIM_PROP_STRING = 2
XIM_PROP_DOUBLE_ARRAY = 4
XIM_PROP_INT_ARRAY = 5

MM_PER_INCH = 25.4

ImageLike = Union["DicomImage", "ArrayImage", "FileImage", "LinacDicomImage"]


def equate_images(image1: ImageLike, image2: ImageLike) -> tuple[ImageLike, ImageLike]:
    """Crop and resize two images to make them:
      * The same pixel dimensions
      * The same DPI

    The usefulness of the function comes when trying to compare images from different sources.
    The best example is calculating gamma on a machine log fluence and EPID image. The physical
    and pixel dimensions must be normalized, the SID normalized

    Parameters
    ----------
    image1 : {:class:`~pylinac.core.image.ArrayImage`, :class:`~pylinac.core.image.DicomImage`, :class:`~pylinac.core.image.FileImage`}
        Must have DPI and SID.
    image2 : {:class:`~pylinac.core.image.ArrayImage`, :class:`~pylinac.core.image.DicomImage`, :class:`~pylinac.core.image.FileImage`}
        Must have DPI and SID.

    Returns
    -------
    image1 : :class:`~pylinac.core.image.ArrayImage`
        The first image equated.
    image2 : :class:`~pylinac.core.image.ArrayImage`
        The second image equated.
    """
    image1 = copy.deepcopy(image1)
    image2 = copy.deepcopy(image2)
    # crop images to be the same physical size
    # ...crop height
    physical_height_diff = image1.physical_shape[0] - image2.physical_shape[0]
    if physical_height_diff < 0:  # image2 is bigger
        img = image2
    else:
        img = image1
    pixel_height_diff = abs(int(round(-physical_height_diff * img.dpmm / 2)))
    if pixel_height_diff > 0:
        img.crop(pixel_height_diff, edges=("top", "bottom"))

    # ...crop width
    physical_width_diff = image1.physical_shape[1] - image2.physical_shape[1]
    if physical_width_diff > 0:
        img = image1
    else:
        img = image2
    pixel_width_diff = abs(int(round(physical_width_diff * img.dpmm / 2)))
    if pixel_width_diff > 0:
        img.crop(pixel_width_diff, edges=("left", "right"))

    # resize images to be of the same shape
    zoom_factor = image1.shape[1] / image2.shape[1]
    image2_array = ndimage.interpolation.zoom(image2.as_type(float), zoom_factor)
    image2 = load(image2_array, dpi=image2.dpi * zoom_factor)

    return image1, image2


def is_image(path: str | io.BytesIO | ImageLike | np.ndarray) -> bool:
    """Determine whether the path is a valid image file.

    Returns
    -------
    bool
    """
    return any((_is_array(path), _is_dicom(path), _is_image_file(path)))


def retrieve_image_files(path: str) -> list[str]:
    """Retrieve the file names of all the valid image files in the path.

    Returns
    -------
    list
        Contains strings pointing to valid image paths.
    """
    return retrieve_filenames(directory=path, func=is_image)


def load(path: str | Path | ImageLike | np.ndarray | BinaryIO, **kwargs) -> ImageLike:
    r"""Load a DICOM image, JPG/TIF/BMP image, or numpy 2D array.

    Parameters
    ----------
    path : str, file-object
        The path to the image file or data stream or array.
    kwargs
        See :class:`~pylinac.core.image.FileImage`, :class:`~pylinac.core.image.DicomImage`,
        or :class:`~pylinac.core.image.ArrayImage` for keyword arguments.

    Returns
    -------
    ::class:`~pylinac.core.image.FileImage`, :class:`~pylinac.core.image.ArrayImage`, or :class:`~pylinac.core.image.DicomImage`
        Return type depends on input image.

    Examples
    --------
    Load an image from a file and then apply a filter::

        >>> from pylinac.core.image import load
        >>> my_image = r"C:\QA\image.tif"
        >>> img = load(my_image)  # returns a FileImage
        >>> img.filter(5)

    Loading from an array is just like loading from a file::

        >>> arr = np.arange(36).reshape(6, 6)
        >>> img = load(arr)  # returns an ArrayImage
    """
    if isinstance(path, BaseImage):
        return path

    if _is_array(path):
        return ArrayImage(path, **kwargs)
    elif _is_dicom(path):
        return DicomImage(path, **kwargs)
    elif _is_image_file(path):
        return FileImage(path, **kwargs)
    else:
        raise TypeError(
            f"The argument `{path}` was not found to be a valid DICOM file, Image file, or array"
        )


def load_url(url: str, progress_bar: bool = True, **kwargs) -> ImageLike:
    """Load an image from a URL.

    Parameters
    ----------
    url : str
        A string pointing to a valid URL that points to a file.

        .. note:: For some images (e.g. Github), the raw binary URL must be used, not simply the basic link.

    progress_bar: bool
        Whether to display a progress bar of download status.
    """
    filename = get_url(url, progress_bar=progress_bar)
    return load(filename, **kwargs)


def load_multiples(
    image_file_list: Sequence,
    method: str = "mean",
    stretch_each: bool = True,
    loader: callable = load,
    **kwargs,
) -> ImageLike:
    """Combine multiple image files into one superimposed image.

    Parameters
    ----------
    image_file_list : list
        A list of the files to be superimposed.
    method : {'mean', 'max', 'sum'}
        A string specifying how the image values should be combined.
    stretch_each : bool
        Whether to normalize the images being combined by stretching their high/low values to the same values across images.
    loader: callable
        The function to use to load the images. If a special image subclass is used, this is how it can be passed.
    kwargs :
        Further keyword arguments are passed to the load function and stretch function.

    Examples
    --------
    Load multiple images::

        >>> from pylinac.core.image import load_multiples
        >>> paths = ['starshot1.tif', 'starshot2.tif']
        >>> superimposed_img = load_multiples(paths)
    """
    # load images
    img_list = [loader(path, **kwargs) for path in image_file_list]
    first_img = img_list[0]

    # check that all images are the same size and stretch if need be
    for img in img_list:
        if img.shape != first_img.shape:
            raise ValueError("Images were not the same shape")
        if stretch_each:
            img.array = stretcharray(img.array, fill_dtype=kwargs.get("dtype"))

    # stack and combine arrays
    new_array = np.dstack(tuple(img.array for img in img_list))
    if method == "mean":
        combined_arr = np.mean(new_array, axis=2)
    elif method == "max":
        combined_arr = np.max(new_array, axis=2)
    elif method == "sum":
        combined_arr = np.sum(new_array, axis=2)

    # replace array of first object and return
    first_img.array = combined_arr
    # set the raw pixels flag; this will mark the image to be converted if we save out
    first_img._raw_pixels = True
    return first_img


def _rescale_dicom_values(
    unscaled_array: np.ndarray, metadata: Dataset, raw_pixels: bool
) -> np.ndarray:
    """Rescale the DICOM pixel values depending on the tags available.

    See Also
    --------
    https://pylinac.readthedocs.io/en/latest/topics/images.html#pixel-data-inversion
    """
    has_all_rescale_tags = (
        hasattr(metadata, "RescaleSlope")
        and hasattr(metadata, "RescaleIntercept")
        and hasattr(metadata, "PixelIntensityRelationshipSign")
    )
    has_some_rescale_tags = hasattr(metadata, "RescaleSlope") and hasattr(
        metadata, "RescaleIntercept"
    )
    is_ct_storage = metadata.SOPClassUID.name == "CT Image Storage"
    is_mr_storage = metadata.SOPClassUID.name == "MR Image Storage"
    if raw_pixels:
        return unscaled_array
    elif has_all_rescale_tags:
        scaled_array = (
            (metadata.RescaleSlope * unscaled_array) + metadata.RescaleIntercept
        ) * metadata.PixelIntensityRelationshipSign
    elif is_ct_storage or has_some_rescale_tags:
        scaled_array = (
            metadata.RescaleSlope * unscaled_array
        ) + metadata.RescaleIntercept
    elif is_mr_storage:
        # signal is usually correct as-is, no inversion needed
        scaled_array = unscaled_array
    else:
        # invert it
        orig_array = unscaled_array
        scaled_array = -orig_array + orig_array.max() + orig_array.min()
    return scaled_array


def _unscale_dicom_values(
    scaled_array: np.ndarray, metadata: Dataset, raw_pixels: bool
) -> np.ndarray:
    """Unscale the DICOM pixel values depending on the tags available.

    This is the inverse of _rescale_dicom_values; specifically, when we
    want to save the DICOM image we want to save the raw values
    back to such that when re-importing and rescaling we will get the same array.
    """
    has_all_rescale_tags = (
        hasattr(metadata, "RescaleSlope")
        and hasattr(metadata, "RescaleIntercept")
        and hasattr(metadata, "PixelIntensityRelationshipSign")
    )
    has_some_rescale_tags = hasattr(metadata, "RescaleSlope") and hasattr(
        metadata, "RescaleIntercept"
    )
    is_ct_storage = metadata.SOPClassUID.name == "CT Image Storage"
    is_mr_storage = metadata.SOPClassUID.name == "MR Image Storage"
    if raw_pixels:
        return scaled_array
    elif has_all_rescale_tags:
        unscaled_array = scaled_array * metadata.PixelIntensityRelationshipSign
        unscaled_array = (
            unscaled_array - metadata.RescaleIntercept
        ) / metadata.RescaleSlope
    elif is_ct_storage or has_some_rescale_tags:
        unscaled_array = (
            scaled_array - metadata.RescaleIntercept
        ) / metadata.RescaleSlope
    elif is_mr_storage:
        # signal is usually correct as-is, no inversion needed
        unscaled_array = scaled_array
    else:
        # invert it
        orig_array = scaled_array
        unscaled_array = -orig_array + orig_array.max() + orig_array.min()
    return unscaled_array


def _is_dicom(path: str | Path | io.BytesIO | ImageLike | np.ndarray) -> bool:
    """Whether the file is a readable DICOM file via pydicom."""
    return is_dicom_image(file=path)


def _is_image_file(path: str | Path) -> bool:
    """Whether the file is a readable image file via Pillow."""
    try:
        with pImage.open(path):
            return True
    except Exception:
        return False


def _is_array(obj: Any) -> bool:
    """Whether the object is a numpy array."""
    return isinstance(obj, np.ndarray)


class BaseImage:
    """Base class for the Image classes.

    Attributes
    ----------
    path : str
        The path to the image file.
    array : numpy.ndarray
        The actual image pixel array.
    """

    array: np.ndarray
    path: str | Path
    metrics: list[MetricBase]
    metric_values: dict[str, Any]
    source: FILE_TYPE | STREAM_TYPE

    def __init__(
        self, path: str | Path | BytesIO | ImageLike | np.ndarray | BufferedReader
    ):
        """
        Parameters
        ----------
        path : str
            The path to the image.
        """
        self.metrics = []
        self.metric_values = {}
        if isinstance(path, (str, Path)) and not osp.isfile(path):
            raise FileExistsError(
                f"File `{path}` does not exist. Verify the file path name."
            )
        elif isinstance(path, (str, Path)) and osp.isfile(path):
            self.path = path
            self.base_path = osp.basename(path)
            self.source = FILE_TYPE
        else:
            self.source = STREAM_TYPE
            path.seek(0)
            try:
                self.path = str(Path(path.name))
            except AttributeError:
                self.path = ""

    @property
    def truncated_path(
        self,
    ) -> str:  # TODO: Use textwrap or pull out into util function
        if self.source == FILE_TYPE:
            path = str(self.path)
            if len(path) > PATH_TRUNCATION_LENGTH:
                return (
                    path[: PATH_TRUNCATION_LENGTH // 2]
                    + "..."
                    + path[-PATH_TRUNCATION_LENGTH // 2 :]
                )
            else:
                return path
        else:
            return ""  # was from stream, no path

    @classmethod
    def from_multiples(
        cls,
        filelist: list[str],
        method: str = "mean",
        stretch: bool = True,
        **kwargs,
    ) -> ImageLike:
        """Load an instance from multiple image items. See :func:`~pylinac.core.image.load_multiples`."""
        return load_multiples(filelist, method, stretch, **kwargs)

    @property
    def center(self) -> Point:
        """Return the center position of the image array as a Point.
        Even-length arrays will return the midpoint between central two indices. Odd will return the central index.
        """
        x_center = (self.shape[1] / 2) - 0.5
        y_center = (self.shape[0] / 2) - 0.5
        return Point(x_center, y_center)

    @property
    def physical_shape(self) -> (float, float):
        """The physical size of the image in mm."""
        return self.shape[0] / self.dpmm, self.shape[1] / self.dpmm

    def date_created(self, format: str = "%A, %B %d, %Y") -> str:
        """The date the file was created. Tries DICOM data before falling back on OS timestamp.
        The method use one or more inputs of formatted code, where % means a placeholder and the letter the time unit of interest.
        For a full description of the several formatting codes see `strftime() documentation. <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`_

        Parameters
        ----------
        format : str
            %A means weekday full name, %B month full name, %d day of the month as a zero-padded decimal number and %Y year with century as a decimal number.

        Returns
        -------
        str
            The date the file was created.
        """
        date = None
        try:
            date = datetime.strptime(
                self.metadata.InstanceCreationDate
                + str(round(float(self.metadata.InstanceCreationTime))),
                "%Y%m%d%H%M%S",
            )
            date = date.strftime(format)
        except (AttributeError, ValueError):
            try:
                date = datetime.strptime(self.metadata.StudyDate, "%Y%m%d")
                date = date.strftime(format)
            except Exception:
                pass
        if date is None:
            try:
                date = datetime.fromtimestamp(osp.getctime(self.path)).strftime(format)
            except AttributeError:
                date = "Unknown"
        return date

    def plot(
        self,
        ax: plt.Axes = None,
        show: bool = True,
        clear_fig: bool = False,
        metric_kwargs: dict | None = None,
        **kwargs,
    ) -> plt.Axes:
        """Plot the image.

        Parameters
        ----------
        ax : matplotlib.Axes instance
            The axis to plot the image to. If None, creates a new figure.
        show : bool
            Whether to actually show the image. Set to false when plotting multiple items.
        clear_fig : bool
            Whether to clear the prior items on the figure before plotting.
        metric_kwargs : dict
            kwargs passed to the metric plot method.
        kwargs
            kwargs passed to plt.imshow()
        """
        if metric_kwargs is None:
            metric_kwargs = {}
        if ax is None:
            fig, ax = plt.subplots()
        if clear_fig:
            plt.clf()
        ax.imshow(self.array, cmap=get_dicom_cmap(), **kwargs)
        # plot the metrics
        for metric in self.metrics:
            metric.plot(axis=ax, **metric_kwargs)
        if show:
            plt.show()
        return ax

    def plot_metrics(self, show: bool = True) -> list[plt.figure]:
        """Plot any additional figures from the metrics.

        Returns a list of figures of the metrics. These metrics are not
        drawn on the original image but rather are something complete separate.
        E.g. a profile plot or a histogram of the metric."""
        figs = []
        for metric in self.metrics:
            figs.append(metric.additional_plots())
        if show:
            plt.show()
        return figs

    def filter(
        self,
        size: float | int = 0.05,
        kind: str = "median",
    ) -> None:
        """Filter the profile in place.

        Parameters
        ----------
        size : int, float
            Size of the median filter to apply.
            If a float, the size is the ratio of the length. Must be in the range 0-1.
            E.g. if size=0.1 for a 1000-element array, the filter will be 100 elements.
            If an int, the filter is the size passed.
        kind : {'median', 'gaussian'}
            The kind of filter to apply. If gaussian, *size* is the sigma value.
        """
        self.array = filter(self.array, size=size, kind=kind)

    def crop(
        self,
        pixels: int = 15,
        edges: tuple[str, ...] = ("top", "bottom", "left", "right"),
    ) -> None:
        """Removes pixels on all edges of the image in-place.

        Parameters
        ----------
        pixels : int
            Number of pixels to cut off all sides of the image.
        edges : tuple
            Which edges to remove from. Can be any combination of the four edges.
        """
        if pixels <= 0:
            raise ValueError("Pixels to remove must be a positive number")
        if "top" in edges:
            self.array = self.array[pixels:, :]
        if "bottom" in edges:
            self.array = self.array[:-pixels, :]
        if "left" in edges:
            self.array = self.array[:, pixels:]
        if "right" in edges:
            self.array = self.array[:, :-pixels]

    def flipud(self) -> None:
        """Flip the image array upside down in-place. Wrapper for np.flipud()"""
        self.array = np.flipud(self.array)

    def fliplr(self) -> None:
        """Flip the image array upside down in-place. Wrapper for np.fliplr()"""
        self.array = np.fliplr(self.array)

    def invert(self) -> None:
        """Invert (imcomplement) the image."""
        self.array = invert(self.array)

    def bit_invert(self) -> None:
        """Invert the image bit-wise"""
        self.array = bit_invert(self.array)

    def roll(self, direction: str = "x", amount: int = 1) -> None:
        """Roll the image array around in-place. Wrapper for np.roll().

        Parameters
        ----------
        direction : {'x', 'y'}
            The axis to roll over.
        amount : int
            The amount of elements to roll over.
        """
        axis = 1 if direction == "x" else 0
        self.array = np.roll(self.array, amount, axis=axis)

    def rot90(self, n: int = 1) -> None:
        """Wrapper for numpy.rot90; rotate the array by 90 degrees CCW n times."""
        self.array = np.rot90(self.array, n)

    def rotate(self, angle: float, mode: str = "edge", *args, **kwargs):
        """Rotate the image counter-clockwise. Simple wrapper for scikit-image. See https://scikit-image.org/docs/stable/api/skimage.transform.html#skimage.transform.rotate.
        All parameters are passed to that function."""
        self.array = rotate(self.array, angle, mode=mode, *args, **kwargs)

    def threshold(self, threshold: float, kind: str = "high") -> None:
        """Apply a high- or low-pass threshold filter.

        Parameters
        ----------
        threshold : int
            The cutoff value.
        kind : str
            If ``high`` (default), will apply a high-pass threshold. All values above the cutoff are left as-is.
            Remaining points are set to 0.
            If ``low``, will apply a low-pass threshold.
        """
        if kind == "high":
            self.array = np.where(self.array >= threshold, self, 0)
        else:
            self.array = np.where(self.array <= threshold, self, 0)

    def as_binary(self, threshold: int) -> ImageLike:
        """Return a binary (black & white) image based on the given threshold.

        Parameters
        ----------
        threshold : int, float
            The threshold value. If the value is above or equal to the threshold it is set to 1, otherwise to 0.

        Returns
        -------
        ArrayImage
        """
        array = np.where(self.array >= threshold, 1, 0)
        return ArrayImage(array)

    def dist2edge_min(self, point: Point | tuple) -> float:
        """Calculates distance from given point to the closest edge.

        Parameters
        ----------
        point : geometry.Point, tuple

        Returns
        -------
        float
        """
        if isinstance(point, tuple):
            point = Point(point)
        rows = self.shape[0]
        cols = self.shape[1]
        disttoedge = np.zeros(4)
        disttoedge[0] = rows - point.y
        disttoedge[1] = cols - point.x
        disttoedge[2] = point.y
        disttoedge[3] = point.x
        return min(disttoedge)

    def ground(self) -> float:
        """Ground the profile in place such that the lowest value is 0.

        .. note::
            This will also "ground" profiles that are negative or partially-negative.
            For such profiles, be careful that this is the behavior you desire.

        Returns
        -------
        float
            The amount subtracted from the image.
        """
        min_val = self.array.min()
        self.array = ground(self.array)
        return min_val

    def normalize(self, norm_val: str | float | None = None) -> None:
        """Normalize the profile to the given value.

        Parameters
        ----------
        value : number or None
            If a number, normalize the array to that number. If None, normalizes to the maximum value.
        """
        # backwards compatibility
        if norm_val == "max":
            norm_val = None
        self.array = normalize(self.array, value=norm_val)

    def check_inversion(
        self, box_size: int = 20, position: (float, float) = (0.0, 0.0)
    ) -> None:
        """Check the image for inversion by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.

        Parameters
        ----------
        box_size : int
            The size in pixels of the corner box to detect inversion.
        position : 2-element sequence
            The location of the sampling boxes.
        """
        row_pos = max(int(position[0] * self.array.shape[0]), 1)
        col_pos = max(int(position[1] * self.array.shape[1]), 1)
        lt_upper = self.array[
            row_pos : row_pos + box_size, col_pos : col_pos + box_size
        ]
        rt_upper = self.array[
            row_pos : row_pos + box_size, -col_pos - box_size : -col_pos
        ]
        lt_lower = self.array[
            -row_pos - box_size : -row_pos, col_pos : col_pos + box_size
        ]
        rt_lower = self.array[
            -row_pos - box_size : -row_pos, -col_pos - box_size : -col_pos
        ]
        avg = np.mean((lt_upper, lt_lower, rt_upper, rt_lower))
        if avg > np.mean(self.array.flatten()):
            self.invert()

    def check_inversion_by_histogram(
        self, percentiles: (float, float, float) = (5, 50, 95)
    ) -> bool:
        """Check the inversion of the image using histogram analysis. The assumption is that the image
        is mostly background-like values and that there is a relatively small amount of dose getting to the image
        (e.g. a picket fence image). This function looks at the distance from one percentile to another to determine
        if the image should be inverted.

        Parameters
        ----------
        percentiles : 3-element tuple
            The 3 percentiles to compare. Default is (5, 50, 95). Recommend using (x, 50, y). To invert the other way
            (where pixel value is *decreasing* with dose, reverse the percentiles, e.g. (95, 50, 5).

        Returns
        -------
        bool: Whether an inversion was performed.
        """
        was_inverted = False
        p_low = np.percentile(self.array, percentiles[0])
        p_mid = np.percentile(self.array, percentiles[1])
        p_high = np.percentile(self.array, percentiles[2])
        mid_to_low = abs(p_mid - p_low)
        mid_to_high = abs(p_mid - p_high)
        if mid_to_low > mid_to_high:
            was_inverted = True
            self.invert()
        return was_inverted

    @argue.bounds(threshold=(0.0, 1.0))
    def gamma(
        self,
        comparison_image: ImageLike,
        doseTA: float = 1,
        distTA: float = 1,
        threshold: float = 0.1,
        ground: bool = True,
        normalize: bool = True,
    ) -> np.ndarray:
        """Calculate the gamma between the current image (reference) and a comparison image.

        .. versionadded:: 1.2

        The gamma calculation is based on `Bakai et al
        <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
        which is a quicker alternative to the standard Low gamma equation.

        Parameters
        ----------
        comparison_image : {:class:`~pylinac.core.image.ArrayImage`, :class:`~pylinac.core.image.DicomImage`, or :class:`~pylinac.core.image.FileImage`}
            The comparison image. The image must have the same DPI/DPMM to be comparable.
            The size of the images must also be the same.
        doseTA : int, float
            Dose-to-agreement in percent; e.g. 2 is 2%.
        distTA : int, float
            Distance-to-agreement in mm.
        threshold : float
            The dose threshold percentage of the maximum dose, below which is not analyzed.
            Must be between 0 and 1.
        ground : bool
            Whether to "ground" the image values. If true, this sets both datasets to have the minimum value at 0.
            This can fix offset errors in the data.
        normalize : bool
            Whether to normalize the images. This sets the max value of each image to the same value.

        Returns
        -------
        gamma_map : numpy.ndarray
            The calculated gamma map.

        See Also
        --------
        :func:`~pylinac.core.image.equate_images`
        """
        # error checking
        if not is_close(self.dpi, comparison_image.dpi, delta=0.1):
            raise AttributeError(
                f"The image DPIs to not match: {self.dpi:.2f} vs. {comparison_image.dpi:.2f}"
            )
        same_x = is_close(self.shape[1], comparison_image.shape[1], delta=1.1)
        same_y = is_close(self.shape[0], comparison_image.shape[0], delta=1.1)
        if not (same_x and same_y):
            raise AttributeError(
                f"The images are not the same size: {self.shape} vs. {comparison_image.shape}"
            )

        # set up reference and comparison images
        ref_img = ArrayImage(copy.copy(self.array))
        ref_img.check_inversion_by_histogram()
        if ground:
            ref_img.ground()
        if normalize:
            ref_img.normalize()
        comp_img = ArrayImage(copy.copy(comparison_image.array))
        comp_img.check_inversion_by_histogram()
        if ground:
            comp_img.ground()
        if normalize:
            comp_img.normalize()

        # invalidate dose values below threshold so gamma doesn't calculate over it
        ref_img.array[ref_img < threshold * np.max(ref_img)] = np.NaN

        # convert distance value from mm to pixels
        distTA_pixels = self.dpmm * distTA

        # construct image gradient using sobel filter
        img_x = spf.sobel(ref_img.as_type(np.float32), 1)
        img_y = spf.sobel(ref_img.as_type(np.float32), 0)
        grad_img = np.hypot(img_x, img_y)

        # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        subtracted_img = np.abs(comp_img - ref_img)
        denominator = np.sqrt(
            ((doseTA / 100.0) ** 2) + ((distTA_pixels**2) * (grad_img**2))
        )
        gamma_map = subtracted_img / denominator

        return gamma_map

    def as_type(self, dtype: np.dtype) -> np.ndarray:
        return self.array.astype(dtype)

    def compute(self, metrics: list[MetricBase] | MetricBase) -> Any | dict[str, Any]:
        """Compute the given metrics on the image.

        This can be called multiple times to compute different metrics.
        Metrics are appended on each call. This allows for modification
        of the image between metric calls as well as the ability to compute
        different metrics on the same image that might depend on
        earlier metrics.

        Metrics are both returned and stored in the ``metrics`` attribute.
        The ``metrics`` attribute will store all metrics every calculated.
        The metrics returned are only those passed in the ``metrics`` argument.

        Parameters
        ----------
        metrics : list[MetricBase] | MetricBase
            The metric(s) to compute.
        """
        metric_data = {}
        if isinstance(metrics, MetricBase):
            metrics = [metrics]
        for metric in metrics:
            metric.inject_image(self)
            self.metrics.append(metric)
            value = metric.context_calculate()
            metric_data[metric.name] = value
        self.metric_values |= metric_data
        if len(metrics) == 1:
            return metric_data[metrics[0].name]
        return metric_data

    @property
    def shape(self) -> (int, int):
        return self.array.shape

    @property
    def size(self) -> int:
        return self.array.size

    @property
    def ndim(self) -> int:
        return self.array.ndim

    @property
    def dtype(self) -> np.dtype:
        return self.array.dtype

    def sum(self) -> float:
        return self.array.sum()

    def ravel(self) -> np.ndarray:
        return self.array.ravel()

    @property
    def flat(self) -> np.ndarray:
        return self.array.flat

    def __len__(self):
        return len(self.array)

    def __getitem__(self, item):
        return self.array[item]


class XIM(BaseImage):
    """A class to open, read, and/or export an .xim image, Varian's custom image format which is 99.999% PNG

    This had inspiration from a number of places:
    - https://gist.github.com/1328/7da697c71f9c4ef12e1e
    - https://medium.com/@duhroach/how-png-works-f1174e3cc7b7
    - https://www.mathworks.com/matlabcentral/answers/419228-how-to-write-for-loop-and-execute-data
    - https://www.w3.org/TR/PNG-Filters.html
    - https://bitbucket.org/dmoderesearchtools/ximreader/src/master/
    """

    array: np.ndarray  #:
    properties: dict  #:

    def __init__(self, file_path: str | Path, read_pixels: bool = True):
        """
        Parameters
        ----------
        file_path
            The path to the file of interest.
        read_pixels
            Whether to read and parse the pixel information. Doing so is quite slow.
            Set this to false if, e.g., you are searching for images only via tags or doing
            a pre-filtering of image selection.
        """
        super().__init__(path=file_path)
        with open(self.path, "rb") as xim:
            self.format_id = decode_binary(xim, str, 8)
            self.format_version = decode_binary(xim, int)
            self.img_width_px = decode_binary(xim, int)
            self.img_height_px = decode_binary(xim, int)
            self.bits_per_pixel = decode_binary(xim, int)
            self.bytes_per_pixel = decode_binary(xim, int)
            self.compression = decode_binary(xim, int)
            if not self.compression:
                pixel_buffer_size = decode_binary(xim, int)
                self.pixel_buffer = decode_binary(
                    xim, str, num_values=pixel_buffer_size
                )
            else:
                lookup_table_size = decode_binary(xim, int)
                self.lookup_table = decode_binary(
                    xim, "B", num_values=lookup_table_size
                )
                comp_pixel_buffer_size = decode_binary(xim, int)
                if read_pixels:
                    lookup_keys = self._parse_lookup_table(self.lookup_table)
                    self.array = self._parse_compressed_bytes(
                        xim, lookup_table=lookup_keys
                    )
                else:
                    _ = decode_binary(xim, "c", num_values=comp_pixel_buffer_size)
                decode_binary(xim, int)
            self.num_hist_bins = decode_binary(xim, int)
            self.histogram = decode_binary(xim, int, num_values=self.num_hist_bins)
            self.num_properties = decode_binary(xim, int)
            self.properties = {}
            for prop in range(self.num_properties):
                name_length = decode_binary(xim, int)
                name = decode_binary(xim, str, num_values=name_length)
                tipe = decode_binary(xim, int)
                if tipe == XIM_PROP_INT:
                    value = decode_binary(xim, int)
                elif tipe == XIM_PROP_DOUBLE:
                    value = decode_binary(xim, "d")
                elif tipe == XIM_PROP_STRING:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(xim, str, num_values=num_bytes)
                elif tipe == XIM_PROP_DOUBLE_ARRAY:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(
                        xim, "d", num_values=int(num_bytes // 8)
                    )  # doubles are 8 bytes
                elif tipe == XIM_PROP_INT_ARRAY:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(
                        xim, int, num_values=int(num_bytes // 4)
                    )  # ints are 4 bytes
                self.properties[name] = value

    @staticmethod
    def _parse_lookup_table(lookup_table_bytes: np.ndarray) -> np.ndarray:
        """The lookup table doesn't follow normal structure conventions like 1, 2, or 4 byte values. They
        got smart and said each value is 2 bits. Yes, bits. This means each byte is actually 4 values.
        Python only reads things as granular as bytes. To get around this the general logic is:

        1) interpret the data as integers at the single byte level
        2) convert those integers back into bit representation; e.g. 115 => 01110011. Note the representation must contain the full byte. I.e. 3 => 11 does not work.
        3) split the binary representation into the 2-bit representations; generates 4x the number of elements. 01110011 => (01, 11, 00, 11)
        4) Convert the 2-bit representation back into integers (01, 11, 00, 11) => (1, 3, 0, 3)

        .. note::

            This is ripe for optimization, but brevity and clarity won out. Options include bit-shifting (fastest)
            and numpy.packbits/unpackbits.
        """
        table = []
        extend = table.extend  # prevent python having to do a lookup on each iteration
        for byte in lookup_table_bytes:
            byte_repr = f"{byte:08b}"
            # didn't actually check these indexes but I think they're right.
            extend(
                [
                    int(byte_repr[6:8], 2),
                    int(byte_repr[4:6], 2),
                    int(byte_repr[2:4], 2),
                    int(byte_repr[0:2], 2),
                ]
            )
        return np.asarray(table, dtype=np.int8)

    def _parse_compressed_bytes(
        self, xim: BinaryIO, lookup_table: np.ndarray
    ) -> np.ndarray:
        """Parse the compressed pixels. We have to do this pixel-by-pixel because each
        pixel can have a different number of bytes representing it

        Per the readme:

        1) The first row is uncompressed
        2) The first element of the second row is uncompressed
        3) all other elements are represented by 1, 2, or 4 bytes of data (the annoying part)
        4) The byte size of the element is given in the lookup table

        So, we have to read in 1, 2, or 4 bytes and convert to an integer depending on
        the lookup table, which tells us how many bytes to read in

        .. note::

            Optimization can help here. A few ideas:

            - reading in groups of data of the same byte size. I already tried this, and I think it will work, but I couldn't get it going.
            - reading in rows of data where no byte change occurred in that row. Similar to above.
            - Using joblib or a processpool
        """
        img_height = self.img_height_px
        img_width = self.img_width_px
        dtype = np.int8 if self.bytes_per_pixel == 1 else np.int16
        compressed_array = a = np.zeros((img_height * img_width), dtype=dtype)
        # first row and 1st element, 2nd row is uncompressed
        # this SHOULD work by reading the # of bytes specified in the header but AFAICT this is just a standard int (4 bytes)
        compressed_array[: img_width + 1] = decode_binary(
            xim, int, num_values=img_width + 1
        )
        diffs = self._get_diffs(lookup_table, xim)
        for diff, idx in zip(diffs, range(img_width + 1, img_width * img_height)):
            # intermediate math can cause overflow errors. Use float for intermediate, then back to int
            left = float(a[idx - 1])
            above = float(a[idx - img_width])
            upper_left = float(a[idx - img_width - 1])
            a[idx] = np.asarray(diff + left + above - upper_left, dtype=dtype)
        return a.reshape((img_height, img_width))

    @staticmethod
    def _get_diffs(lookup_table: np.ndarray, xim: BinaryIO):
        """Read in all the pixel value 'diffs'. These can be 1, 2, or 4 bytes in size,
        so instead of just reading N pixels of M bytes which would be SOOOO easy, we have to read dynamically

        We optimize here by reading bytes in clumps, which is way faster than reading one at a time.
        Knowing that most values are single bytes with an occasional 2-byte element
        we read chunks that all look like (n 1-bytes and 1 2-byte)
        """
        byte_changes = lookup_table.nonzero()
        byte_changes = np.insert(byte_changes, 0, -1)
        byte_changes = np.append(byte_changes, len(lookup_table) - 1)
        diffs = [5000] * (
            len(lookup_table) - 1
        )  # pre-allocate for speed; 5000 is just for debugging
        LOOKUP_CONVERSION = {0: "b", 1: "h", 2: "i"}
        for start, stop in zip(byte_changes[:-1], byte_changes[1:]):
            if stop - start > 1:
                vals = decode_binary(xim, "b", num_values=stop - start - 1)
                if not isinstance(vals, Iterable):
                    vals = [
                        vals,
                    ]
                diffs[start + 1 : stop] = vals
            if stop != byte_changes[-1]:
                diffs[stop] = decode_binary(xim, LOOKUP_CONVERSION[lookup_table[stop]])
        return np.asarray(diffs, dtype=float)

    @property
    def dpmm(self) -> float:
        """The dots/mm value of the XIM images. The value appears to be in cm in the file."""
        if self.properties["PixelWidth"] != self.properties["PixelHeight"]:
            raise ValueError(
                "The XIM image does not have the same pixel height and width"
            )
        return 1 / (10 * self.properties["PixelHeight"])

    def save_as(self, file: str, format: str | None = None) -> None:
        """Save the image to a NORMAL format. PNG is highly suggested. Accepts any format supported by Pillow.
        Ironically, an equivalent PNG image (w/ metadata) is ~50% smaller than an .xim image.

        .. warning::

            Any format other than PNG will not include the properties included in the .xim image!

        Parameters
        ----------
        file
            The file to save the image to. E.g. my_xim.png
        format
            The format to save the image as. Uses the Pillow logic, which will infer the format if the file name has one.
        """
        img = pImage.fromarray(self.array)
        # we construct the custom PNG tags; it won't be included for tiff or jpeg, etc but it won't error it either.
        metadata = PngInfo()
        for prop, value in self.properties.items():
            if isinstance(value, np.ndarray):
                value = value.tolist()
            if not isinstance(value, str):
                value = json.dumps(value)
            metadata.add_text(prop, value)
        img.save(file, format=format, pnginfo=metadata)


class DicomImage(BaseImage):
    """An image from a DICOM RTImage file.

    Attributes
    ----------
    metadata : pydicom Dataset
        The dataset of the file as returned by pydicom without pixel data.
    """

    metadata: pydicom.FileDataset
    _sid: float
    _dpi: float
    _sad: float

    def __init__(
        self,
        path: str | Path | BytesIO | BufferedReader,
        *,
        dtype: np.dtype | None = None,
        dpi: float = None,
        sid: float = None,
        sad: float = 1000,
        raw_pixels: bool = False,
    ):
        """
        Parameters
        ----------
        path : str, file-object
            The path to the file or the data stream.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        dpi : int, float
            The dots-per-inch of the image, defined at isocenter.

            .. note:: If a DPI tag is found in the image, that value will override the parameter, otherwise this one
                will be used.

        sid : int, float
            The Source-to-Image distance in mm.
        sad : float
            The Source-to-Axis distance in mm.
        raw_pixels : bool
            Whether to apply pixel intensity correction to the DICOM data.
            Typically, Rescale Slope, Rescale Intercept, and other tags
            are included and meant to be applied to the raw pixel data, which
            is potentially compressed.
            If True, no correction will be applied. This is typically used
            for scenarios when you want to match behavior to older or different
            software.
        """
        super().__init__(path)
        self._sid = sid
        self._dpi = dpi
        self._sad = sad
        # read the file once to get just the DICOM metadata
        self.metadata = retrieve_dicom_file(path)
        self._original_dtype = self.metadata.pixel_array.dtype
        self._raw_pixels = raw_pixels
        # read a second time to get pixel data
        try:
            path.seek(0)
        except AttributeError:
            pass
        ds = retrieve_dicom_file(path)
        if dtype is not None:
            self.array = ds.pixel_array.astype(dtype)
        else:
            self.array = ds.pixel_array.copy()
        # convert values to HU or CU
        self.array = _rescale_dicom_values(self.array, ds, raw_pixels=raw_pixels)

    @classmethod
    def from_dataset(cls, dataset: Dataset):
        """Create a DICOM image instance from a pydicom Dataset."""
        stream = io.BytesIO()
        dataset.save_as(stream)
        return cls(path=stream)

    def save(self, filename: str | Path) -> str | Path:
        """Save the image instance back out to a .dcm file.

        Parameters
        ----------
        filename : str, Path
            The filename to save the DICOM file as.

        Returns
        -------
        A string pointing to the new filename.
        """
        unscaled_array = _unscale_dicom_values(
            self.array, self.metadata, self._raw_pixels
        )
        # if we will have bit overflows, stretch instead
        max_is_too_high = (
            unscaled_array.max() > get_dtype_info(self._original_dtype).max
        )
        min_is_too_low = unscaled_array.min() < get_dtype_info(self._original_dtype).min
        if min_is_too_low or max_is_too_high:
            warnings.warn(
                "The pixel values of image were detected to be outside"
                f"the range of {self._original_dtype} values and will be normalized to fit the original dtype. "
                f"The maximum value will be the maximum value of the original datatype: ({get_dtype_info(self._original_dtype).max})."
            )
            unscaled_array = convert_to_dtype(unscaled_array, self._original_dtype)
        if self._raw_pixels:
            # for raw pixels, often the values are wacky; convert them to a normal range
            # e.g. the range might be 0-1 for a float array. This will convert to the original dtype
            # This prevents such arrays from being converted to only 0s and 1s for an int dtype.
            unscaled_array = convert_to_dtype(unscaled_array, self._original_dtype)
        self.metadata.PixelData = unscaled_array.astype(self._original_dtype).tobytes()
        self.metadata.Columns = unscaled_array.shape[1]
        self.metadata.Rows = unscaled_array.shape[0]
        self.metadata.save_as(filename)
        return filename

    @property
    def z_position(self) -> float:
        """The z-position of the slice. Relevant for CT and MR images."""
        return z_position(self.metadata)

    @property
    def slice_spacing(self) -> float:
        """Determine the distance between slices. The spacing can be greater than the slice thickness (i.e. gaps).
        Uses the absolute version as it can apparently be negative: https://dicom.innolitics.com/ciods/nm-image/nm-reconstruction/00180088

        This attempts to use the slice spacing attr and if it doesn't exist, use the slice thickness attr
        """

        try:
            return abs(self.metadata.SpacingBetweenSlices)
        except AttributeError:
            return self.metadata.SliceThickness

    @property
    def sid(self) -> float:
        """The Source-to-Image in mm."""
        try:
            return float(self.metadata.RTImageSID)
        except (AttributeError, ValueError, TypeError):
            return self._sid

    @property
    def sad(self) -> float:
        """The source to axis (iso) in mm"""
        try:
            return float(self.metadata.RadiationMachineSAD)
        except (AttributeError, ValueError, TypeError):
            return self._sad

    @property
    def dpi(self) -> float:
        """The dots-per-inch of the image, defined at isocenter."""
        try:
            return self.dpmm * MM_PER_INCH
        except Exception:
            return self._dpi

    @property
    def dpmm(self) -> float:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        dpmm = None
        for tag in ("PixelSpacing", "ImagePlanePixelSpacing"):
            mmpd = self.metadata.get(tag)
            if mmpd is not None:
                dpmm = 1 / mmpd[0]
                break
        if dpmm is not None and self.sid is not None:
            dpmm *= self.sid / self.sad
        elif dpmm is None and self._dpi is not None:
            dpmm = self._dpi / MM_PER_INCH
        return dpmm

    @property
    def cax(self) -> Point:
        """The position of the beam central axis. If no DICOM translation tags are found then the center is returned.
        Uses this tag: https://dicom.innolitics.com/ciods/rt-beams-delivery-instruction/rt-beams-delivery-instruction/00741020/00741030/3002000d
        """
        try:
            x = self.center.x - self.metadata.XRayImageReceptorTranslation[0]
            y = self.center.y - self.metadata.XRayImageReceptorTranslation[1]
        except (AttributeError, ValueError, TypeError):
            return self.center
        else:
            return Point(x, y)


class LinacDicomImage(DicomImage):
    """DICOM image taken on a linac. Also allows passing of gantry/coll/couch values via the filename."""

    gantry_keyword = "Gantry"
    collimator_keyword = "Coll"
    couch_keyword = "Couch"

    _use_filenames: bool

    def __init__(
        self,
        path: str | Path | BinaryIO,
        use_filenames: bool = False,
        axes_precision: int | None = None,
        **kwargs,
    ):
        self._gantry = kwargs.pop("gantry", None)
        self._coll = kwargs.pop("coll", None)
        self._couch = kwargs.pop("couch", None)
        self._axes_precision = axes_precision
        super().__init__(path, **kwargs)
        self._use_filenames = use_filenames

    @property
    def gantry_angle(self) -> float:
        """Gantry angle of the irradiation."""
        if self._gantry is not None:
            g = self._gantry
        else:
            g = self._get_axis_value(self.gantry_keyword, "GantryAngle")
        return wrap360(simple_round(g, self._axes_precision))

    @property
    def collimator_angle(self) -> float:
        """Collimator angle of the irradiation."""
        if self._coll is not None:
            c = self._coll
        else:
            c = self._get_axis_value(self.collimator_keyword, "BeamLimitingDeviceAngle")
        return wrap360(simple_round(c, self._axes_precision))

    @property
    def couch_angle(self) -> float:
        """Couch angle of the irradiation."""
        if self._couch is not None:
            c = self._couch
        else:
            c = self._get_axis_value(self.couch_keyword, "PatientSupportAngle")
        return wrap360(simple_round(c, self._axes_precision))

    def _get_axis_value(self, axis_str: str, axis_dcm_attr: str) -> float:
        """Retrieve the value of the axis. This will first look in the file name for the value.
        If not in the filename then it will look in the DICOM metadata. If the value can be found in neither
        then a value of 0 is assumed.

        Parameters
        ----------
        axis_str : str
            The string to look for in the filename.
        axis_dcm_attr : str
            The DICOM attribute that should contain the axis value.

        Returns
        -------
        float
        """
        axis_found = False
        if self._use_filenames:
            filename = osp.basename(self.path)
            # see if the keyword is in the filename
            keyword_in_filename = axis_str.lower() in filename.lower()
            # if it's not there, then assume it's zero
            if not keyword_in_filename:
                axis = 0
                axis_found = True
            # if it is, then make sure it follows the naming convention of <axis###>
            else:
                match = re.search(rf"(?<={axis_str.lower()})\d+", filename.lower())
                if match is None:
                    raise ValueError(
                        f"The filename contains '{axis_str}' but could not read a number following it. Use the format '...{axis_str}<#>...'"
                    )
                else:
                    axis = float(match.group())
                    axis_found = True
        # try to interpret from DICOM data
        if not axis_found:
            try:
                axis = float(getattr(self.metadata, axis_dcm_attr))
            except AttributeError:
                axis = 0
        return axis


class FileImage(BaseImage):
    """An image from a "regular" file (.tif, .jpg, .bmp).

    Attributes
    ----------
    info : dict
        The info dictionary as generated by Pillow.
    sid : float
        The SID value as passed in upon construction.
    """

    def __init__(
        self,
        path: str | Path | BinaryIO,
        *,
        dpi: float | None = None,
        sid: float | None = None,
        dtype: np.dtype | None = None,
    ):
        """
        Parameters
        ----------
        path : str, file-object
            The path to the file or a data stream.
        dpi : int, float
            The dots-per-inch of the image, defined at isocenter.

            .. note:: If a DPI tag is found in the image, that value will override the parameter, otherwise this one
                will be used.
        sid : int, float
            The Source-to-Image distance in mm.
        dtype : numpy.dtype
            The data type to cast the array as. If None, will use the datatype stored in the file.
            If the file is multi-channel (e.g. RGB), it will be converted to int32
        """
        super().__init__(path)
        pil_image = pImage.open(path)
        # convert from multichannel if need be
        if len(pil_image.getbands()) > 1:
            pil_image = pil_image.convert(
                "I"
            )  # int32; uint16 preferred but not reliable using PIL
        self.info = pil_image.info
        try:  # tiff tags
            self.tags = {TAGS[key]: pil_image.tag_v2[key] for key in pil_image.tag_v2}
        except AttributeError:
            pass
        self.array = np.array(pil_image, dtype=dtype)
        self._dpi = dpi
        self.sid = sid

    @property
    def dpi(self) -> float | None:
        """The dots-per-inch of the image, defined at isocenter."""
        dpi = None
        for key in ("dpi", "resolution"):
            dpi = self.info.get(key)
            if dpi is not None:
                dpi = float(dpi[0])
                if dpi < 3 and not self._dpi:
                    raise ValueError(
                        f"The DPI setting is abnormal or nonsensical. Got resolution of {dpi}. Pass in the dpi manually."
                    )
                if dpi < 3:
                    dpi = None
                break
        if dpi is None:
            dpi = self._dpi
        if self.sid is not None and dpi is not None:
            dpi *= self.sid / 1000
        return dpi

    @property
    def dpmm(self) -> float | None:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        try:
            return self.dpi / MM_PER_INCH
        except TypeError:
            return


class ArrayImage(BaseImage):
    """An image constructed solely from a numpy array."""

    def __init__(
        self,
        array: np.ndarray,
        *,
        dpi: float = None,
        sid: float = None,
        dtype=None,
    ):
        """
        Parameters
        ----------
        array : numpy.ndarray
            The image array.
        dpi : int, float
            The dots-per-inch of the image, defined at isocenter.

            .. note:: If a DPI tag is found in the image, that value will override the parameter, otherwise this one
                will be used.
        sid : int, float
            The Source-to-Image distance in mm.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        if dtype is not None:
            self.array = np.array(array, dtype=dtype)
        else:
            self.array = array
        self._dpi = dpi
        self.sid = sid
        self.metrics = []
        self.metric_values = {}

    @property
    def dpmm(self) -> float | None:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        try:
            return self.dpi / MM_PER_INCH
        except Exception:
            return

    @property
    def dpi(self) -> float | None:
        """The dots-per-inch of the image, defined at isocenter."""
        dpi = None
        if self._dpi is not None:
            dpi = self._dpi
            if self.sid is not None:
                dpi *= self.sid / 1000
        return dpi

    def __sub__(self, other):
        return ArrayImage(self.array - other.array)


class LazyDicomImageStack:
    _image_path_keys: list[Path]
    metadatas: list[pydicom.Dataset]

    def __init__(
        self,
        folder: str | Path | Sequence[str | Path],
        dtype: np.dtype | None = None,
        min_number: int = 39,
        check_uid: bool = True,
    ):
        """Load a folder with DICOM CT images. This variant is more memory efficient than the standard DicomImageStack.

        This is done by loading images from disk on the fly. This assumes all images remain on disk for the lifetime of the instance. This does not
        need to be true for the original implementation.

        See the documentation for DicomImageStack for parameter descriptions.
        """
        self.dtype = dtype
        paths = []
        # load in images in their received order
        if isinstance(folder, (list, tuple)):
            paths = folder
        elif osp.isdir(folder):
            for pdir, sdir, files in os.walk(folder):
                for file in files:
                    paths.append(osp.join(pdir, file))
        # we only want to read the metadata once
        # so we read it here and then filter and sort
        metadatas, paths = self._get_path_metadatas(paths)

        # check that at least 1 image was loaded
        if len(paths) < 1:
            raise FileNotFoundError(
                f"No files were found in the specified location: {folder}"
            )

        # error checking
        if check_uid:
            most_common_uid = self._get_common_uid_imgs(metadatas, min_number)
            metadatas = [m for m in metadatas if m.SeriesInstanceUID == most_common_uid]
        # sort according to physical order
        order = np.argsort([m.ImagePositionPatient[-1] for m in metadatas])
        self.metadatas = [metadatas[i] for i in order]
        self._image_path_keys = [paths[i] for i in order]

    @classmethod
    def from_zip(cls, zip_path: str | Path, dtype: np.dtype | None = None):
        """Load a DICOM ZIP archive.

        Parameters
        ----------
        zip_path : str
            Path to the ZIP archive.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        with TemporaryZipDirectory(zip_path, delete=False) as tmpzip:
            obj = cls(tmpzip, dtype)
        return obj

    def _get_common_uid_imgs(
        self, metadata: list[pydicom.Dataset], min_number: int
    ) -> pydicom.DataElement:
        """Check that all the images are from the same study."""
        most_common_uid = Counter(i.SeriesInstanceUID for i in metadata).most_common(1)[
            0
        ]
        if most_common_uid[1] < min_number:
            raise ValueError(
                "The minimum number images from the same study were not found"
            )
        return most_common_uid[0]

    def _get_path_metadatas(
        self, paths: list[Path]
    ) -> (list[pydicom.Dataset], list[Path]):
        """Get the metadata for the images. This also filters out non-image files."""
        metadata = []
        matched_paths = []
        for path in paths:
            try:
                ds = pydicom.dcmread(path, force=True, stop_before_pixels=True)
                if "Image Storage" in ds.SOPClassUID.name:
                    metadata.append(ds)
                    matched_paths.append(path)
            except (InvalidDicomError, AttributeError, MemoryError):
                pass
        return metadata, matched_paths

    def side_view(self, axis: int) -> np.ndarray:
        """Return the side view of the stack. E.g. if axis=0, return the maximum value along the 0th axis."""
        return np.stack([i for i in self], axis=-1).max(axis=axis)

    @cached_property
    def metadata(self) -> pydicom.FileDataset:
        """The metadata of the first image; shortcut attribute. Only attributes that are common throughout the stack should be used,
        otherwise the individual image metadata should be used."""
        return self[0].metadata

    def __getitem__(self, item: int) -> DicomImage:
        return DicomImage(self._image_path_keys[item], dtype=self.dtype)

    def __setitem__(self, key: int, value: DicomImage):
        """Save the passed image to disk in place of the current image."""
        current_path = self._image_path_keys[key]
        value.save(current_path)

    def __delitem__(self, key):
        """Delete the image from the stack and OS."""
        current_path = self._image_path_keys[key]
        try:
            os.remove(current_path)
        except Exception:
            pass
        del self._image_path_keys[key]

    def __len__(self):
        return len(self._image_path_keys)


class DicomImageStack(LazyDicomImageStack):
    """A class that loads and holds a stack of DICOM images (e.g. a CT dataset). The class can take
    a folder or zip file and will read CT images. The images must all be the same size. Supports
    indexing to individual images.

    Attributes
    ----------
    images : list
        Holds instances of :class:`~pylinac.core.image.DicomImage`. Can be accessed via index;
        i.e. self[0] == self.images[0].

    Examples
    --------
    Load a folder of Dicom images
    >>> from pylinac import image
    >>> img_folder = r"folder/qa/cbct/june"
    >>> dcm_stack = image.DicomImageStack(img_folder)  # loads and sorts the images
    >>> dcm_stack.plot(3)  # plot the 3rd image

    Load a zip archive
    >>> img_folder_zip = r"archive/qa/cbct/june.zip"  # save space and zip your CBCTs
    >>> dcm_stack = image.DicomImageStack.from_zip(img_folder_zip)

    Load as a certain data type
    >>> dcm_stack_uint32 = image.DicomImageStack(img_folder, dtype=np.uint32)
    """

    images: list[ImageLike]

    def __init__(
        self,
        folder: str | Path,
        dtype: np.dtype | None = None,
        min_number: int = 39,
        check_uid: bool = True,
        raw_pixels: bool = False,
    ):
        """Load a folder with DICOM CT images.

        Parameters
        ----------
        folder : str
            Path to the folder.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        super().__init__(folder, dtype, min_number, check_uid)
        self.images = [
            DicomImage(path, dtype=dtype, raw_pixels=raw_pixels)
            for path in self._image_path_keys
        ]

    @classmethod
    def from_zip(cls, zip_path: str | Path, dtype: np.dtype | None = None):
        """Load a DICOM ZIP archive.

        Parameters
        ----------
        zip_path : str
            Path to the ZIP archive.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        with TemporaryZipDirectory(zip_path) as tmpzip:
            obj = cls(tmpzip, dtype)
        return obj

    def plot_3view(self):
        """Plot the stack in 3 views: axial, coronal, and sagittal."""
        fig, axes = plt.subplots(1, 3)
        names = ("Coronal", "Sagittal", "Axial")
        for idx, (ax, name) in enumerate(zip(axes, names)):
            arry = self.side_view(idx)
            ax.imshow(arry, cmap="gray", aspect="equal")
            ax.set_title(name)
        plt.show()

    def roll(self, direction: str, amount: int):
        for img in self.images:
            img.roll(direction, amount)

    def __getitem__(self, item) -> DicomImage:
        return self.images[item]

    def __setitem__(self, key, value: DicomImage):
        self.images[key] = value

    def __delitem__(self, key):
        """Delete the image from the stack."""
        del self.images[key]

    def __len__(self):
        return len(self.images)


class NMImageStack:
    """A class of frames of a nuclear medicine image.
    A single image can have N frames. For our purposes, we
    can treat this as a stack of images."""

    metadata: Dataset
    path: str | Path
    frames: list[DicomImage]

    def __init__(self, path: str | Path):
        """Load a single NM image with N frames."""
        self.path = path
        self.frames = []
        ds = pydicom.dcmread(path, force=True, stop_before_pixels=True)
        if ds.Modality != "NM":
            raise TypeError("The file is not a NM image")
        self.metadata = ds
        full_ds = pydicom.dcmread(path, force=True)
        for i in range(full_ds.NumberOfFrames):
            full_array = full_ds.pixel_array
            if full_array.ndim == 2:
                array = full_array
            else:
                array = full_array[i]
            img = DicomImage(self.path)
            img.array = array
            self.frames.append(img)

    def as_3d_array(self) -> np.ndarray:
        """Return the frames as a 3D array."""
        return np.stack([i.array for i in self.frames], axis=0)

    def __len__(self):
        return len(self.frames)


def tiff_to_dicom(
    tiff_file: str | Path | BytesIO,
    sid: float,
    gantry: float,
    coll: float,
    couch: float,
    dpi: float | None = None,
) -> Dataset:
    """Converts a TIFF file into a **simplistic** DICOM file. Not meant to be a full-fledged tool. Used for conversion so that tools that are traditionally oriented
    towards DICOM have a path to accept TIFF. Currently used to convert files for WL.

    .. note::

        This will convert the image into an uint16 datatype to match the native EPID datatype.

    Parameters
    ----------
    tiff_file
        The TIFF file to be converted.
    sid
        The Source-to-Image distance in mm.
    dpi
        The dots-per-inch value of the TIFF image.
    gantry
        The gantry value that the image was taken at.
    coll
        The collimator value that the image was taken at.
    couch
        The couch value that the image was taken at.
    """
    tiff_img = FileImage(tiff_file, dpi=dpi, sid=sid)
    if not tiff_img.dpmm:
        raise ValueError(
            "Automatic detection of `dpi` failed. A `dpi` value must be passed to the constructor."
        )
    uint_array = convert_to_dtype(tiff_img.array, np.uint16)
    mm_pixel = 25.4 / tiff_img.dpi
    file_meta = FileMetaDataset()
    # Main data elements
    ds = Dataset()
    ds.SOPClassUID = UID("1.2.840.10008.5.1.4.1.1.481.1")
    ds.SOPInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()
    ds.Modality = "RTIMAGE"
    ds.ConversionType = "WSD"
    ds.PatientName = "Lutz^Test Tool"
    ds.PatientID = "Someone Important"
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.Rows = tiff_img.shape[0]
    ds.Columns = tiff_img.shape[1]
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.ImagePlanePixelSpacing = [mm_pixel, mm_pixel]
    ds.RadiationMachineSAD = "1000.0"
    ds.RTImageSID = sid
    ds.PrimaryDosimeterUnit = "MU"
    ds.GantryAngle = str(gantry)
    ds.BeamLimitingDeviceAngle = str(coll)
    ds.PatientSupportAngle = str(couch)
    ds.PixelData = uint_array

    ds.file_meta = file_meta
    ds.is_implicit_VR = True
    ds.is_little_endian = True
    return ds


def gamma_2d(
    reference: np.ndarray,
    evaluation: np.ndarray,
    dose_to_agreement: float = 1,
    distance_to_agreement: int = 1,
    gamma_cap_value: float = 2,
    global_dose: bool = True,
    dose_threshold: float = 5,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Compute a 2D gamma of two 2D numpy arrays. This does NOT do size or spatial resolution checking.
    It performs an element-by-element evaluation. It is the responsibility
    of the caller to ensure the reference and evaluation have comparable spatial resolution.

    The algorithm follows Table I of D. Low's 2004 paper: Evaluation of the gamma dose distribution comparison method: https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.1598711

    This is similar to the gamma_1d function for profiles, except we must search a 2D grid around the reference point.

    Parameters
    ----------
    reference
        The reference 2D array.
    evaluation
        The evaluation 2D array.
    dose_to_agreement
        The dose to agreement in %. E.g. 1 is 1% of global reference max dose.
    distance_to_agreement
        The distance to agreement in **elements**. E.g. if the value is 4 this means 4 elements from the reference point under calculation.
        Must be >0
    gamma_cap_value
        The value to cap the gamma at. E.g. a gamma of 5.3 will get capped to 2. Useful for displaying data with a consistent range.
    global_dose
        Whether to evaluate the dose to agreement threshold based on the global max or the dose point under evaluation.
    dose_threshold
        The dose threshold as a number between 0 and 100 of the % of max dose under which a gamma is not calculated.
        This is not affected by the global/local dose normalization and the threshold value is evaluated against the global max dose, period.
    fill_value
        The value to give pixels that were not calculated because they were under the dose threshold. Default
        is NaN, but another option would be 0. If NaN, allows the user to calculate mean/median gamma over just the
        evaluated portion and not be skewed by 0's that should not be considered.
    """
    if reference.ndim != 2 or evaluation.ndim != 2:
        raise ValueError(
            f"Reference and evaluation arrays must be 2D. Got reference: {reference.ndim} and evaluation: {evaluation.ndim}"
        )
    threshold = reference.max() / 100 * dose_threshold
    # convert dose to agreement to % of global max; ignored later if local dose
    dose_ta = dose_to_agreement / 100 * reference.max()
    # pad eval array on both edges so our search does not go out of bounds
    eval_padded = np.pad(evaluation, distance_to_agreement, mode="edge")
    # iterate over each reference element, computing distance value and dose value
    gamma = np.zeros(reference.shape)
    for row_idx, row in enumerate(reference):
        for col_idx, ref_point in enumerate(row):
            # skip if below dose threshold
            if ref_point < threshold:
                gamma[row_idx, col_idx] = fill_value
                continue
            # use scikit-image to compute the indices of a disk around the reference point
            # we can then compute gamma over the eval points at these indices
            # unlike the 1D computation, we have to search at an index offset by the distance to agreement
            # we use DTA+1 in disk because it looks like the results are exclusive of edges.
            # https://scikit-image.org/docs/stable/api/skimage.draw.html#disk
            rs, cs = disk(
                (row_idx + distance_to_agreement, col_idx + distance_to_agreement),
                distance_to_agreement + 1,
            )

            capital_gammas = []
            for r, c in zip(rs, cs):
                eval_point = eval_padded[r, c]
                # for the distance, we compare the ref row/col to the eval padded matrix
                # but remember the padded array is padded by DTA, so to compare distances, we
                # have to cancel the offset we used for dose purposes.
                dist = math.dist(
                    (row_idx, col_idx),
                    (r - distance_to_agreement, c - distance_to_agreement),
                )
                dose = eval_point - ref_point
                if not global_dose:
                    dose_ta = dose_to_agreement / 100 * ref_point
                capital_gamma = math.sqrt(
                    dist**2 / distance_to_agreement**2 + dose**2 / dose_ta**2
                )
                capital_gammas.append(capital_gamma)
            gamma[row_idx, col_idx] = min(np.nanmin(capital_gammas), gamma_cap_value)
    return np.asarray(gamma)


def z_position(metadata: pydicom.Dataset) -> float:
    """The 'z-position' of the image. Relevant for CT and MR images."""
    try:
        return metadata.ImagePositionPatient[-1]
    except AttributeError:
        return metadata.SliceLocation
