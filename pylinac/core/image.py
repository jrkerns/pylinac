"""This module holds classes for image loading and manipulation."""
import copy
from collections import Counter
from datetime import datetime
from io import BytesIO
import re
import os.path as osp
import os
from typing import Union, Sequence, List, Any, Tuple, Optional

import pydicom
from pydicom.errors import InvalidDicomError
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as pImage
from scipy import ndimage
import scipy.ndimage.filters as spf
import argue

from .utilities import is_close
from .geometry import Point
from .io import get_url, TemporaryZipDirectory, retrieve_filenames, is_dicom_image, retrieve_dicom_file
from .profile import stretch as stretcharray
from .typing import NumberLike
from ..settings import get_dicom_cmap, PATH_TRUNCATION_LENGTH

ARRAY = 'Array'
DICOM = 'DICOM'
IMAGE = 'Image'

FILE_TYPE = 'file'
STREAM_TYPE = 'stream'

MM_PER_INCH = 25.4

ImageLike = Union['DicomImage', 'ArrayImage', 'FileImage', 'LinacDicomImage']


def equate_images(image1: ImageLike, image2: ImageLike) -> Tuple[ImageLike, ImageLike]:
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
    image2 : :class:`~pylinac.core.image.ArrayImage`
        The returns are new instances of Images.
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
    img.remove_edges(pixel_height_diff, edges=('top', 'bottom'))

    # ...crop width
    physical_width_diff = image1.physical_shape[1] - image2.physical_shape[1]
    if physical_width_diff > 0:
        img = image1
    else:
        img = image2
    pixel_width_diff = abs(int(round(physical_width_diff*img.dpmm/2)))
    img.remove_edges(pixel_width_diff, edges=('left', 'right'))

    # resize images to be of the same shape
    zoom_factor = image1.shape[1] / image2.shape[1]
    image2_array = ndimage.interpolation.zoom(image2.as_type(np.float), zoom_factor)
    image2 = load(image2_array, dpi=image2.dpi * zoom_factor)

    return image1, image2


def is_image(path: str) -> bool:
    """Determine whether the path is a valid image file.

    Returns
    -------
    bool
    """
    return any((_is_array(path), _is_dicom(path), _is_image_file(path)))


def retrieve_image_files(path: str) -> List:
    """Retrieve the file names of all the valid image files in the path.

    Returns
    -------
    list
        Contains strings pointing to valid image paths.
    """
    return retrieve_filenames(directory=path, func=is_image)


def load(path: str, **kwargs) -> ImageLike:
    """Load a DICOM image, JPG/TIF/BMP image, or numpy 2D array.

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
        >>> my_image = "C:\QA\image.tif"
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
        raise TypeError(f"The argument `{path}` was not found to be a valid DICOM file, Image file, or array")


def load_url(url: str, progress_bar: bool=True, **kwargs) -> ImageLike:
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


@argue.options(method=('mean', 'max', 'sum'))
def load_multiples(image_file_list: List, method: str='mean', stretch_each: bool=True, **kwargs) -> ImageLike:
    """Combine multiple image files into one superimposed image.

    Parameters
    ----------
    image_file_list : list
        A list of the files to be superimposed.
    method : {'mean', 'max', 'sum'}
        A string specifying how the image values should be combined.
    stretch_each : bool
        Whether to normalize the images being combined by stretching their high/low values to the same values across images.
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
    img_list = [load(path, **kwargs) for path in image_file_list]
    first_img = img_list[0]

    # check that all images are the same size and stretch if need be
    for img in img_list:
        if img.shape != first_img.shape:
            raise ValueError("Images were not the same shape")
        if stretch_each:
            img.array = stretcharray(img.array, fill_dtype=kwargs.get('dtype'))

    # stack and combine arrays
    new_array = np.dstack(tuple(img.array for img in img_list))
    if method == 'mean':
        combined_arr = np.mean(new_array, axis=2)
    elif method == 'max':
        combined_arr = np.max(new_array, axis=2)
    elif method == 'sum':
        combined_arr = np.sum(new_array, axis=2)

    # replace array of first object and return
    first_img.array = combined_arr
    return first_img


def _is_dicom(path: str) -> bool:
    """Whether the file is a readable DICOM file via pydicom."""
    return is_dicom_image(file=path)


def _is_image_file(path: str) -> bool:
    """Whether the file is a readable image file via Pillow."""
    try:
        pImage.open(path)
        return True
    except:
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

    def __init__(self, path: str):
        """
        Parameters
        ----------
        path : str
            The path to the image.
        """
        path: str
        source: Union[FILE_TYPE, STREAM_TYPE]

        if isinstance(path, str) and not osp.isfile(path):
            raise FileExistsError(f"File `{path}` does not exist. Verify the file path name.")
        elif isinstance(path, str) and osp.isfile(path):
            self.path = path
            self.base_path = osp.basename(path)
            self.source = FILE_TYPE
        else:
            self.source = STREAM_TYPE

    @property
    def truncated_path(self) -> str:
        if len(self.path) > PATH_TRUNCATION_LENGTH:
            return self.path[:PATH_TRUNCATION_LENGTH // 2] + '...' + self.path[-PATH_TRUNCATION_LENGTH // 2:]
        else:
            return self.path

    @classmethod
    def from_multiples(cls, filelist: List[str], method: str='mean', stretch: bool=True, **kwargs) -> ImageLike:
        """Load an instance from multiple image items. See :func:`~pylinac.core.image.load_multiples`."""
        return load_multiples(filelist, method, stretch, **kwargs)

    @property
    def center(self) -> Point:
        """Return the center position of the image array as a Point."""
        x_center = (self.shape[1] / 2) - 0.5
        y_center = (self.shape[0] / 2) - 0.5
        return Point(x_center, y_center)

    @property
    def physical_shape(self) -> Tuple[float, float]:
        """The physical size of the image in mm."""
        return self.shape[0] / self.dpmm, self.shape[1] / self.dpmm

    def date_created(self, format: str="%A, %B %d, %Y") -> str:
        date = None
        try:
            date = datetime.strptime(self.metadata.InstanceCreationDate+str(round(float(self.metadata.InstanceCreationTime))), "%Y%m%d%H%M%S")
            date = date.strftime(format)
        except (AttributeError, ValueError):
            try:
                date = datetime.strptime(self.metadata.StudyDate, "%Y%m%d")
                date = date.strftime(format)
            except:
                pass
        if date is None:
            try:
                date = datetime.fromtimestamp(osp.getctime(self.path)).strftime(format)
            except AttributeError:
                date = 'Unknown'
        return date

    def plot(self, ax: plt.Axes=None, show: bool=True, clear_fig: bool=False, **kwargs) -> plt.Axes:
        """Plot the image.

        Parameters
        ----------
        ax : matplotlib.Axes instance
            The axis to plot the image to. If None, creates a new figure.
        show : bool
            Whether to actually show the image. Set to false when plotting multiple items.
        clear_fig : bool
            Whether to clear the prior items on the figure before plotting.
        """
        if ax is None:
            fig, ax = plt.subplots()
        if clear_fig:
            plt.clf()
        ax.imshow(self.array, cmap=get_dicom_cmap(), **kwargs)
        if show:
            plt.show()
        return ax

    @argue.options(kind=('median', 'gaussian'))
    def filter(self, size: Union[float, int]=0.05, kind: str='median') -> None:
        """Filter the profile.

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
        if isinstance(size, float):
            if 0 < size < 1:
                size *= len(self.array)
                size = max(size, 1)
            else:
                raise TypeError("Float was passed but was not between 0 and 1")

        if kind == 'median':
            self.array = ndimage.median_filter(self.array, size=size)
        elif kind == 'gaussian':
            self.array = ndimage.gaussian_filter(self.array, sigma=size)

    def crop(self, pixels: int=15, edges: Tuple[str, ...]=('top', 'bottom', 'left', 'right')) -> None:
        """Removes pixels on all edges of the image in-place.

        Parameters
        ----------
        pixels : int
            Number of pixels to cut off all sides of the image.
        edges : tuple
            Which edges to remove from. Can be any combination of the four edges.
        """
        if pixels < 0:
            raise ValueError("Pixels to remove must be a positive number")
        if 'top' in edges:
            self.array = self.array[pixels:, :]
        if 'bottom' in edges:
            self.array = self.array[:-pixels, :]
        if 'left' in edges:
            self.array = self.array[:, pixels:]
        if 'right' in edges:
            self.array = self.array[:, :-pixels]

    def remove_edges(self, pixels: int=15, edges: Tuple[str, ...]=('top', 'bottom', 'left', 'right')) -> None:
        """Removes pixels on all edges of the image in-place.

        Parameters
        ----------
        pixels : int
            Number of pixels to cut off all sides of the image.
        edges : tuple
            Which edges to remove from. Can be any combination of the four edges.
        """
        DeprecationWarning("`remove_edges` is deprecated and will be removed in a future version. Use `crop` instead")
        self.crop(pixels=pixels, edges=edges)

    def flipud(self) -> None:
        """ Flip the image array upside down in-place. Wrapper for np.flipud()"""
        self.array = np.flipud(self.array)

    def invert(self) -> None:
        """Invert (imcomplement) the image."""
        orig_array = self.array
        self.array = -orig_array + orig_array.max() + orig_array.min()

    def roll(self, direction: str='x', amount: int=1) -> None:
        """Roll the image array around in-place. Wrapper for np.roll().

        Parameters
        ----------
        direction : {'x', 'y'}
            The axis to roll over.
        amount : int
            The amount of elements to roll over.
        """
        axis = 1 if direction == 'x' else 0
        self.array = np.roll(self.array, amount, axis=axis)

    def rot90(self, n: int=1) -> None:
        """Wrapper for numpy.rot90; rotate the array by 90 degrees CCW."""
        self.array = np.rot90(self.array, n)

    @argue.options(kind=('high', 'low'))
    def threshold(self, threshold: int, kind: str='high') -> None:
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
        if kind == 'high':
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

    def dist2edge_min(self, point: Union[Point, Tuple]) -> float:
        """Calculates minimum distance from given point to image edges.

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
        """Ground the profile such that the lowest value is 0.

        .. note::
            This will also "ground" profiles that are negative or partially-negative.
            For such profiles, be careful that this is the behavior you desire.

        Returns
        -------
        float
            The amount subtracted from the image.
        """
        min_val = self.array.min()
        self.array -= min_val
        return min_val

    def normalize(self, norm_val: Union[str, NumberLike]='max') -> None:
        """Normalize the image values to the given value.

        Parameters
        ----------
        norm_val : str, number
            If a string, must be 'max', which normalizes the values to the maximum value.
            If a number, normalizes all values to that number.
        """
        if norm_val == 'max':
            val = self.array.max()
        else:
            val = norm_val
        self.array = self.array / val

    def check_inversion(self, box_size: int=20, position: Sequence=(0.0, 0.0)) -> None:
        """Check the image for inversion by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.

        Parameters
        ----------
        box_size : int
            The size in pixels of the corner box to detect inversion.
        position : 2-element sequence
            The location of the sampling boxes.
        """
        row_pos = max(int(position[0]*self.array.shape[0]), 1)
        col_pos = max(int(position[1]*self.array.shape[1]), 1)
        lt_upper = self.array[row_pos: row_pos+box_size, col_pos: col_pos+box_size]
        rt_upper = self.array[row_pos: row_pos+box_size, -col_pos-box_size: -col_pos]
        lt_lower = self.array[-row_pos-box_size:-row_pos, col_pos: col_pos+box_size]
        rt_lower = self.array[-row_pos-box_size:-row_pos, -col_pos-box_size:-col_pos]
        avg = np.mean((lt_upper, lt_lower, rt_upper, rt_lower))
        if avg > np.mean(self.array.flatten()):
            self.invert()

    def check_inversion_by_histogram(self, percentiles=(5, 50, 95)) -> bool:
        """Check the inversion of the image using histogram analysis. The assumption is that the image
        is mostly background-like values and that there is a relatively small amount of dose getting to the image
        (e.g. a picket fence image). This function looks at the distance from one percentile to another to determine
        if the image should be inverted.

        Parameters
        ----------
        percentiles : 3-element tuple
            The 3 percentiles to compare. Default is (5, 50, 95). Recommend using (x, 50, y). To invert the other way
            (where pixel value is *decreasing* with dose, reverse the percentiles, e.g. (95, 50, 5).
        """
        was_inverted = False
        p5 = np.percentile(self.array, percentiles[0])
        p50 = np.percentile(self.array, percentiles[1])
        p95 = np.percentile(self.array, percentiles[2])
        dist_to_5 = abs(p50 - p5)
        dist_to_95 = abs(p50 - p95)
        if dist_to_5 > dist_to_95:
            was_inverted = True
            self.invert()
        return was_inverted

    @argue.bounds(threshold=(0.0, 1.0))
    def gamma(self, comparison_image: ImageLike, doseTA: NumberLike=1, distTA: NumberLike=1,
              threshold: NumberLike=0.1, ground: bool=True, normalize: bool=True) -> np.ndarray:
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
            raise AttributeError(f"The image DPIs to not match: {self.dpi:.2f} vs. {comparison_image.dpi:.2f}")
        same_x = is_close(self.shape[1], comparison_image.shape[1], delta=1.1)
        same_y = is_close(self.shape[0], comparison_image.shape[0], delta=1.1)
        if not (same_x and same_y):
            raise AttributeError(f"The images are not the same size: {self.shape} vs. {comparison_image.shape}")

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
        denominator = np.sqrt(((doseTA / 100.0) ** 2) + ((distTA_pixels ** 2) * (grad_img ** 2)))
        gamma_map = subtracted_img / denominator

        return gamma_map

    def as_type(self, dtype) -> np.ndarray:
        return self.array.astype(dtype)

    @property
    def shape(self) -> Tuple[int, int]:
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


class DicomImage(BaseImage):
    """An image from a DICOM RTImage file.

    Attributes
    ----------
    metadata : pydicom Dataset
        The dataset of the file as returned by pydicom without pixel data.
    """
    metadata: pydicom.FileDataset
    _sid = NumberLike
    _dpi = NumberLike

    def __init__(self, path: str, *, dtype=None, dpi: NumberLike=None, sid: NumberLike=None):
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
        """
        super().__init__(path)
        self._sid = sid
        self._dpi = dpi
        # read the file once to get just the DICOM metadata
        self.metadata = retrieve_dicom_file(path)
        self._original_dtype = self.metadata.pixel_array.dtype
        # read a second time to get pixel data
        if isinstance(path, BytesIO):
            path.seek(0)
        ds = retrieve_dicom_file(path)
        if dtype is not None:
            self.array = ds.pixel_array.astype(dtype)
        else:
            self.array = ds.pixel_array.copy()
        # convert values to HU or CU: real_values = slope * raw + intercept
        is_ct_storage = self.metadata.SOPClassUID.name == 'CT Image Storage'
        has_rescale_tags = hasattr(self.metadata, 'RescaleSlope') and hasattr(self.metadata, 'RescaleIntercept')
        if is_ct_storage or has_rescale_tags:
            self.array = (self.metadata.RescaleSlope*self.array) + self.metadata.RescaleIntercept

    def save(self, filename: str) -> str:
        """Save the image instance back out to a .dcm file.

        Returns
        -------
        A string pointing to the new filename.
        """
        if self.metadata.SOPClassUID.name == 'CT Image Storage':
            self.array = (self.array - int(self.metadata.RescaleIntercept)) / int(self.metadata.RescaleSlope)
        self.metadata.PixelData = self.array.astype(self._original_dtype).tobytes()
        self.metadata.save_as(filename)
        return filename

    @property
    def sid(self) -> NumberLike:
        """The Source-to-Image in mm."""
        try:
            return float(self.metadata.RTImageSID)
        except:
            return self._sid

    @property
    def dpi(self) -> NumberLike:
        """The dots-per-inch of the image, defined at isocenter."""
        try:
            return self.dpmm * MM_PER_INCH
        except:
            return self._dpi

    @property
    def dpmm(self) -> NumberLike:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        dpmm = None
        for tag in ('PixelSpacing', 'ImagePlanePixelSpacing'):
            mmpd = self.metadata.get(tag)
            if mmpd is not None:
                dpmm = 1 / mmpd[0]
                break
        if dpmm is not None and self.sid is not None:
            dpmm *= self.sid / 1000
        elif dpmm is None and self._dpi is not None:
            dpmm = self._dpi / MM_PER_INCH
        return dpmm

    @property
    def cax(self) -> Point:
        """The position of the beam central axis. If no DICOM translation tags are found then the center is returned."""
        try:
            x = self.center.x - self.metadata.XRayImageReceptorTranslation[0]
            y = self.center.y - self.metadata.XRayImageReceptorTranslation[1]
        except AttributeError:
            return self.center
        else:
            return Point(x, y)


class LinacDicomImage(DicomImage):
    """DICOM image taken on a linac. Also allows passing of gantry/coll/couch values via the filename."""
    gantry_keyword = 'Gantry'
    collimator_keyword = 'Coll'
    couch_keyword = 'Couch'

    _use_filenames: bool

    def __init__(self, path: str, use_filenames: bool=False):
        super().__init__(path)
        self._use_filenames = use_filenames

    @property
    def gantry_angle(self) -> float:
        """Gantry angle of the irradiation."""
        return self._get_axis_value(self.gantry_keyword.lower(), 'GantryAngle')

    @property
    def collimator_angle(self) -> float:
        """Collimator angle of the irradiation."""
        return self._get_axis_value(self.collimator_keyword.lower(), 'BeamLimitingDeviceAngle')

    @property
    def couch_angle(self) -> float:
        """Couch angle of the irradiation."""
        return self._get_axis_value(self.couch_keyword.lower(), 'PatientSupportAngle')

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
                match = re.search('(?<={})\d+'.format(axis_str.lower()), filename.lower())
                if match is None:
                    raise ValueError(
                            f"The filename contains '{axis_str}' but could not read a number following it. Use the format '...{axis_str}<#>...'")
                else:
                    axis = float(match.group())
                    axis_found = True
        # try to interpret from DICOM data
        if not axis_found:
            try:
                axis = float(getattr(self.metadata, axis_dcm_attr))
            except AttributeError:
                axis = 0
        # if the value is close to 0 or 360 then peg at 0
        if is_close(axis, [0, 360], delta=1):
            return 0
        else:
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

    def __init__(self, path: str, *, dpi: NumberLike=None, sid: NumberLike=None, dtype=None):
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
            The data type to cast the array as.
        """
        super().__init__(path)
        pil_image = pImage.open(path)
        # convert to gray if need be
        if pil_image.mode not in ('F', 'L', '1'):
            pil_image = pil_image.convert('F')
        self.info = pil_image.info
        if dtype is not None:
            self.array = np.array(pil_image, dtype=dtype)
        else:
            self.array = np.array(pil_image)
        self._dpi = dpi
        self.sid = sid

    @property
    def dpi(self) -> float:
        """The dots-per-inch of the image, defined at isocenter."""
        dpi = None
        for key in ('dpi', 'resolution'):
            dpi = self.info.get(key)
            if dpi is not None:
                dpi = float(dpi[0])
                break
        if dpi is None:
            dpi = self._dpi
        if self.sid is not None and dpi is not None:
            dpi *= self.sid / 1000
        return dpi

    @property
    def dpmm(self) -> float:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        try:
            return self.dpi / MM_PER_INCH
        except TypeError:
            return


class ArrayImage(BaseImage):
    """An image constructed solely from a numpy array."""

    def __init__(self, array: np.array, *, dpi: NumberLike=None, sid: NumberLike=None, dtype=None):
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

    @property
    def dpmm(self) -> Optional[float]:
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        try:
            return self.dpi / MM_PER_INCH
        except:
            return

    @property
    def dpi(self) -> Optional[float]:
        """The dots-per-inch of the image, defined at isocenter."""
        dpi = None
        if self._dpi is not None:
            dpi = self._dpi
            if self.sid is not None:
                dpi *= self.sid / 1000
        return dpi

    def __sub__(self, other):
        return ArrayImage(self.array - other.array)


class DicomImageStack:
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
    images: List

    def __init__(self, folder: str, dtype=None, min_number: int=39, check_uid: bool=True):
        """Load a folder with DICOM CT images.

        Parameters
        ----------
        folder : str
            Path to the folder.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        self.images = []
        # load in images in their received order
        for pdir, sdir, files in os.walk(folder):
            for file in files:
                path = osp.join(pdir, file)
                if self.is_CT_slice(path):
                    img = DicomImage(path, dtype=dtype)
                    self.images.append(img)

        # check that at least 1 image was loaded
        if len(self.images) < 1:
            raise FileNotFoundError(f"No files were found in the specified location: {folder}")

        # error checking
        if check_uid:
            self.images = self._check_number_and_get_common_uid_imgs(min_number)
        # sort according to physical order
        self.images.sort(key=lambda x: x.metadata.ImagePositionPatient[-1])

    @classmethod
    def from_zip(cls, zip_path: str, dtype=None):
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

    @staticmethod
    def is_CT_slice(file: str) -> bool:
        """Test if the file is a CT Image storage DICOM file."""
        try:
            ds = pydicom.dcmread(file, force=True, stop_before_pixels=True)
            return ds.SOPClassUID.name == 'CT Image Storage'
        except (InvalidDicomError, AttributeError, MemoryError):
            return False

    def _check_number_and_get_common_uid_imgs(self, min_number: int) -> List:
        """Check that all the images are from the same study."""
        most_common_uid = Counter(i.metadata.SeriesInstanceUID for i in self.images).most_common(1)[0]
        if most_common_uid[1] < min_number:
            raise ValueError("The minimum number images from the same study were not found")
        return [i for i in self.images if i.metadata.SeriesInstanceUID == most_common_uid[0]]

    def plot(self, slice: int=0) -> None:
        """Plot a slice of the DICOM dataset.

        Parameters
        ----------
        slice : int
            The slice to plot.
        """
        self.images[slice].plot()

    @property
    def metadata(self) -> pydicom.FileDataset:
        """The metadata of the first image; shortcut attribute. Only attributes that are common throughout the stack should be used,
        otherwise the individual image metadata should be used."""
        return self.images[0].metadata

    def __getitem__(self, item) -> DicomImage:
        return self.images[item]

    def __setitem__(self, key, value: DicomImage):
        self.images[key] = value

    def __len__(self):
        return len(self.images)
