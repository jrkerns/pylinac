"""This module holds classes for image loading and manipulation."""
import copy
from io import BytesIO
import os.path as osp
import os

import dicom
from dicom.errors import InvalidDicomError
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as pImage
from scipy import ndimage
from scipy.misc import imresize
import scipy.ndimage.filters as spf

from .decorators import type_accept, value_accept
from .geometry import Point
from .io import get_url, TemporaryZipDirectory
from .profile import stretch as stretcharray

ARRAY = 'Array'
DICOM = 'DICOM'
IMAGE = 'Image'
MM_PER_INCH = 25.4


class Image:
    """A swiss-army knife, delegate class for loading in images and image-like things.

    The class should not be instantiated directly, but through its class methods. These methods
    return not an `Image` class but one of three specialized image classes:

    * :class:`~pylinac.core.image.DicomImage` : Handles all DICOM images; utilizes pydicom.
    * :class:`~pylinac.core.image.FileImage` : Handles JPEG, BMP, TIF, and other "regular" image files; utilizes Pillow.
    * :class:`~pylinac.core.image.ArrayImage` : Handles 2D numpy arrays; convenient for doing processing of arrays
      that represent an image.

    There are two methods to construct these classes:

    * :meth:`~pylinac.core.image.Image.load` : For loading single images/arrays.
    * :meth:`~pylinac.core.image.Image.load_multiples` : For loading and superimposing multiple images/arrays. All the images must be the same size.

    Examples
    --------
    Load an image from a file::

        >>> my_image = "C:\QA\image.tif"
        >>> img = Image.load(my_image)  # returns a FileImage
        >>> img.filter(5)

    Loading from an array is just like loading from a file::

        >>> arr = np.arange(36).reshape(6, 6)
        >>> img = Image.load(arr)  # returns an ArrayImage

    Load multiple images::

        >>> paths = ['starshot1.tif', 'starshot2.tif']
        >>> superimposed_img = Image.load_multiples(paths)
    """

    @classmethod
    def load(cls, path, **kwargs):
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
        """
        if cls._is_array(path):
            return ArrayImage(path, **kwargs)
        elif cls._is_dicom(path):
            return DicomImage(path, **kwargs)
        elif cls._is_image_file(path):
            return FileImage(path, **kwargs)
        else:
            raise TypeError("The argument `{0}` was not found to be a valid DICOM file, Image file, or array".format(path))

    @classmethod
    def load_url(cls, url, **kwargs):
        """Load an image from a URL.

        Parameters
        ----------
        url : str
            A string pointing to a valid URL that points to a file.

            .. note:: For some images (e.g. Github), the raw binary URL must be used, not simply the basic link.
        """
        filename = get_url(url)
        return cls.load(filename, **kwargs)

    @classmethod
    @value_accept(method=('mean', 'max', 'sum'))
    def load_multiples(cls, image_file_list, method='mean', stretch=True, **kwargs):
        """Combine multiple image files into one superimposed image.

        .. versionadded:: 0.5.1
        """
        # load images
        img_list = [cls.load(path, **kwargs) for path in image_file_list]
        first_img = img_list[0]

        # check that all images are the same size and stretch if need be
        for img in img_list:
            if img.shape != first_img.shape:
                raise ValueError("Images were not the same shape")
            if stretch:
                img.array = stretcharray(img.array)

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
        first_img.check_inversion()
        return first_img

    @staticmethod
    def _is_dicom(path):
        """Whether the file is a readable DICOM file via pydicom."""
        try:
            ds = dicom.read_file(path, stop_before_pixels=True, force=True)
            ds.SOPClassUID
            return True
        except:
            return False

    @staticmethod
    def _is_image_file(path):
        """Whether the file is a readable image file via Pillow."""
        try:
            pImage.open(path)
            return True
        except:
            return False

    @staticmethod
    def _is_array(obj):
        """Whether the object is a numpy array."""
        return isinstance(obj, np.ndarray)


class BaseImage:
    """Base class for the Image classes.

    Attributes
    ----------
    array : numpy.ndarray
        The actual image pixel array.
    """

    def __init__(self, path):
        if isinstance(path, BytesIO):
            path.seek(0)
        elif not osp.isfile(path):
            raise FileExistsError("File `{0}` does not exist".format(path))
        else:
            self.filename = path

    @classmethod
    def from_multiples(cls, filelist, method='mean', stretch=True, **kwargs):
        img_list = [cls(path, **kwargs) for path in filelist]
        first_img = img_list[0]

        # check that all images are the same size and stretch if need be
        for img in img_list:
            if img.shape != first_img.shape:
                raise ValueError("Images were not the same shape")
            if stretch:
                img.array = stretcharray(img.array)

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
        first_img.check_inversion()
        return first_img

    @property
    def center(self):
        """Return the center position of the image array as a Point."""
        x_center = self.shape[1] / 2
        y_center = self.shape[0] / 2
        return Point(x_center, y_center)

    def plot(self, ax=None, show=True, clear_fig=False):
        """Plot the image."""
        if ax is None:
            fig, ax = plt.subplots()
        if clear_fig:
            plt.clf()
        ax.imshow(self.array, cmap=plt.cm.Greys)
        if show:
            plt.show()
        return ax

    @value_accept(kind=('median', 'gaussian'))
    def filter(self, size=0.05, kind='median'):
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

    @type_accept(pixels=int)
    def remove_edges(self, pixels=15):
        """Removes pixels on all edges of the image.

        Parameters
        ----------
        pixels : int
            Number of pixels to cut off all sides of the image.
        """
        self.array = self.array[pixels - 1:-pixels - 1, pixels - 1:-pixels - 1]

    def invert(self):
        """Invert (imcomplement) the image."""
        orig_array = self.array
        self.array = -orig_array + orig_array.max() + orig_array.min()

    def roll(self, direction='x', amount=1):
        axis = 1 if direction == 'x' else 0
        self.array = np.roll(self.array, amount, axis=axis)

    def rot90(self, n=1):
        """Wrapper for numpy.rot90."""
        self.array = np.rot90(self.array, n)

    def resize(self, size, interp='bilinear'):
        """Resize/scale the image.

        Wrapper for scipy's `imresize <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.misc.imresize.html>`_:
        """
        self.array = imresize(self.array, size=size, interp=interp, mode='F')

    @value_accept(kind=('high', 'low'))
    def threshold(self, threshold, kind='high'):
        """Apply a threshold filter.

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

    def as_binary(self, threshold):
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

    @type_accept(point=(Point, tuple))
    def dist2edge_min(self, point):
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

    def ground(self):
        """Ground the profile such that the lowest value is 0.

        .. note::
            This will also "ground" profiles that are negative or partially-negative.
            For such profiles, be careful that this is the behavior you desire.
        """
        min_val = self.array.min()
        self.array -= min_val
        return min_val

    def normalize(self, norm_val='max'):
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

    def check_inversion(self, box_size=20, offset=10):
        """Check the image for inversion by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.
        """
        outer_edge = offset
        inner_edge = offset + box_size
        TL_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BL_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        TR_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BR_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        corner_avg = np.mean((TL_corner, BL_corner, TR_corner, BR_corner))
        if corner_avg > np.mean(self.array.flatten()):
            self.invert()

    def gamma(self, comparison_image, doseTA=1, distTA=1, threshold=0.1):
        """Calculate the gamma between the current image (reference) and a comparison image.

        .. versionadded:: 1.2

        The gamma calculation is based on `Bakai et al
        <http://iopscience.iop.org/0031-9155/48/21/006/>`_ eq.6,
        which is a quicker alternative to the standard Low gamma equation.

        Parameters
        ----------
        comparison_image :
            The comparison image. Must be a :class:`~pylinac.core.image.ArrayImage`,
            :class:`~pylinac.core.image.DicomImage`, or :class:`~pylinac.core.image.FileImage`.
            The image must have the same DPI/DPMM to be comparable.
            The size of the images must also be the same.
        doseTA : int, float
            Dose-to-agreement in percent; e.g. 2 is 2%.
        distTA : int, float
            Distance-to-agreement in mm.
        threshold : float
            The dose threshold percentage of the maximum dose, below which is not analyzed.
            Must be between 0 and 1.

        Returns
        -------
        gamma_map : numpy.ndarray
            The calculated gamma map.
        """
        # error checking
        if self.dpi != comparison_image.dpi:
            raise AttributeError("The image DPIs to not match; cannot perform gamma calculation.")
        if self.shape != comparison_image.shape:
            raise AttributeError("The images are not the same size; cannot perform gamma calculation.")

        # set up reference and comparison images
        ref_img = ArrayImage(copy.copy(self.array))
        ref_img.check_inversion()
        ref_img.ground()
        ref_img.normalize()
        comp_img = ArrayImage(copy.copy(comparison_image.array))
        comp_img.check_inversion()
        comp_img.ground()
        comp_img.normalize()

        # invalidate dose values below threshold so gamma doesn't calculate over it
        ref_img.array[ref_img < threshold * np.max(ref_img)] = np.NaN

        # convert distance value from mm to pixels
        distTA_pixels = self.dpmm * distTA

        # construct image gradient using sobel filter
        img_x = spf.sobel(ref_img, 1)
        img_y = spf.sobel(ref_img, 0)
        grad_img = img_x + img_y

        # equation: (measurement - reference) / sqrt ( doseTA^2 + distTA^2 * image_gradient^2 )
        subtracted_img = np.abs(comp_img - ref_img)
        denominator = np.sqrt((doseTA / 100.0 ** 2) + ((distTA_pixels ** 2) * (grad_img ** 2)))
        gamma_map = subtracted_img / denominator

        return gamma_map

    @property
    def shape(self):
        return self.array.shape

    @property
    def size(self):
        return self.array.size

    @property
    def ndim(self):
        return self.array.ndim

    @property
    def dtype(self):
        return self.array.dtype

    def sum(self):
        return self.array.sum()

    def ravel(self):
        return self.array.ravel()

    @property
    def flat(self):
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

    def __init__(self, path, *, dtype=None):
        """
        Parameters
        ----------
        path : str, file-object
            The path to the file or the data stream.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        super().__init__(path)
        # read the file once to get just the DICOM metadata
        self.metadata = dicom.read_file(path, force=True, stop_before_pixels=True)
        # read a second time to get pixel data
        if isinstance(path, BytesIO):
            path.seek(0)
        ds = dicom.read_file(path, force=True)
        if dtype is not None:
            self.array = ds.pixel_array.astype(dtype)
        else:
            self.array = ds.pixel_array
        # convert values to proper HU: real_values = slope * raw + intercept
        if self.metadata.SOPClassUID.name == 'CT Image Storage':
            self.array = int(self.metadata.RescaleSlope)*self.array + int(self.metadata.RescaleIntercept)

    @property
    def sid(self):
        """The Source-to-Image in mm."""
        return float(self.metadata.RTImageSID)

    @property
    def dpi(self):
        """The dots-per-inch of the image, defined at isocenter."""
        return self.dpmm * MM_PER_INCH

    @property
    def dpmm(self):
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        try:
            # most dicom files have this tag
            dpmm = 1 / self.metadata.PixelSpacing[0]
        except AttributeError:
            try:
                # EPID images sometimes have this tag
                dpmm = 1 / self.metadata.ImagePlanePixelSpacing[0]
            except AttributeError:
                raise ("No pixel/distance conversion tag found")
        dpmm *= self.sid / 1000
        return dpmm

    @property
    def cax(self):
        """The position of the beam central axis. If no DICOM translation tags are found then the center is returned."""
        try:
            x = self.center.x - self.metadata.XRayImageReceptorTranslation[0]
            y = self.center.y - self.metadata.XRayImageReceptorTranslation[1]
        except AttributeError:
            return self.center
        else:
            return Point(x, y)


class FileImage(BaseImage):
    """An image from a "regular" file (.tif, .jpg, .bmp).

    Attributes
    ----------
    info : dict
        The info dictionary as generated by Pillow.
    sid : float
        The SID value as passed in upon construction.
    """

    def __init__(self, path, *, dpi=None, sid=1000):
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
        """
        super().__init__(path)
        pil_image = pImage.open(path)
        # convert to gray if need be
        if pil_image.mode not in ('F', 'L', '1'):
            pil_image = pil_image.convert('F')
        self.info = pil_image.info
        self.array = np.array(pil_image)
        self._dpi = dpi
        self.sid = sid

    @property
    def dpi(self):
        """The dots-per-inch of the image, defined at isocenter."""
        try:
            dpi = self.info['dpi'][0]
        except (IndexError, KeyError):
            try:
                dpi = self.info['resolution'][0]
            except:
                dpi = None
        if dpi is None:
            if self._dpi is None:
                raise AttributeError("No pixel/distance conversion tag found. If you know the DPI, pass it in during construction.")
            else:
                dpi = self._dpi
        dpi *= self.sid / 1000
        return dpi

    @property
    def dpmm(self):
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        return self.dpi / MM_PER_INCH


class ArrayImage(BaseImage):
    """An image constructed solely from a numpy array."""

    def __init__(self, array, *, dpi=None, sid=1000, dtype=None):
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
    def dpmm(self):
        """The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm."""
        return self.dpi / MM_PER_INCH

    @property
    def dpi(self):
        """The dots-per-inch of the image, defined at isocenter."""
        if self._dpi is not None:
            return self._dpi
        else:
            raise AttributeError("No pixel/distance conversion value; pass one in during construction")

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
    """

    def __init__(self, folder=None, dtype=None):
        """Load a folder with DICOM CT images.

        Parameters
        ----------
        folder : str
            Path to the folder.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        self.images = []
        if folder is not None:
            self._get_datasets(folder, dtype=dtype)

    @classmethod
    def from_zip(cls, zip_path, dtype=None):
        """Load a DICOM ZIP archive.

        Parameters
        ----------
        zip_path : str
            Path to the ZIP archive.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        obj = cls()
        obj._get_datasets(zip_path, dtype=dtype, zip=True)
        return obj

    def _get_datasets(self, folder, zip=False, dtype=None):
        """Read and load DICOM files from a folder or zip archive.

        Parameters
        ----------
        folder : str
            Path to folder with CT images or path to zip archive of CT images.
        zip : bool
            Whether ``folder`` is a zip archive or not.
        dtype : dtype, None, optional
            The data type to cast the image data as. If None, will use whatever raw image format is.
        """
        # unpack zip archive if passed one
        if zip:
            if osp.isfile(folder):
                temp_folder = TemporaryZipDirectory(folder)
                folder = temp_folder.name
            else:
                raise FileExistsError("Zip archive `{0}` not found".format(folder))

        # load in images in their received order
        for pdir, sdir, files in os.walk(folder):
            for file in files:
                path = osp.join(pdir, file)
                if self._is_CT_slice(path):
                    img = DicomImage(path)
                    self.images.append(img)

        # cleanup zip archive if passed one
        if zip:
            temp_folder.cleanup()

        # check that at least 1 image was loaded
        if len(self.images) < 1:
            raise FileNotFoundError("No files were found in the specified location: {0}".format(folder))

        # error checking
        self._check_all_from_same_study()
        # get the original image order
        original_img_order = [int(round(image.metadata.ImagePositionPatient[-1])) for image in self.images]
        # correctly reorder the images
        self.images = [self.images[i] for i in np.argsort(original_img_order)]

    @staticmethod
    def _is_CT_slice(file):
        """Test if the file is a CT Image storage DICOM file."""
        try:
            ds = dicom.read_file(file, force=True, stop_before_pixels=True)
            return ds.SOPClassUID.name == 'CT Image Storage'
        except (InvalidDicomError, AttributeError, MemoryError):
            return False

    def _check_all_from_same_study(self):
        """Check that all the images are from the same study."""
        initial_uid = self.images[0].metadata.SeriesInstanceUID
        if not all(image.metadata.SeriesInstanceUID == initial_uid for image in self.images):
            raise ValueError("The images were not all from the same study")

    def plot(self, slice=0):
        """Plot a slice of the DICOM dataset."""
        self.images[slice].plot()

    @property
    def metadata(self):
        """The metadata of the first image; shortcut attribute. Only attributes that are common throughout the stack should be used,
        otherwise the individual image metadata should be used."""
        return self.images[0].metadata

    def __getitem__(self, item):
        return self.images[item]

    def __len__(self):
        return len(self.images)
