
"""This module holds classes for image loading and manipulation."""
from io import BytesIO
import os.path as osp
import os
import zipfile

import dicom
from dicom.errors import InvalidDicomError
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as pImage
from scipy import ndimage
from scipy.misc import imresize

from pylinac.core.decorators import type_accept, value_accept
from pylinac.core.geometry import Point
from pylinac.core.profile import stretch as stretcharray
from pylinac.core.utilities import typed_property


ARRAY = 'Array'
DICOM = 'DICOM'
IMAGE = 'Image'
MM_PER_INCH = 25.4


class DICOMStack:
    """A class that loads and holds a stack of DICOM images (e.g. a CT dataset). The class can take
    a folder or zip file and read any and all DICOM images. The images must all be the same size."""
    metadata = None
    array = None

    def __init__(self, folder=None, dtype=int):
        if folder is not None:
            self._instantiate(data=folder, dtype=dtype)

    def _instantiate(self, data, dtype, is_zip=False):
        if is_zip:
            metadatalist = self._get_ct_images_metadata_from_zip(data)
        else:
            metadatalist = self._get_ct_images_metadata_list(data)
        self.metadata = metadatalist[-1]
        self._check_all_from_same_study(metadatalist)
        self._create_array(metadatalist, dtype)

    def _get_ct_images_metadata_list(self, folder):
        # read each file to see if it's a CT Image file
        filelist = []
        for par_dir, sub_dir, files in os.walk(folder):
            for name in files:
                try:
                    ds = dicom.read_file(osp.join(par_dir, name), force=True)
                    if ds.SOPClassUID.name == 'CT Image Storage':
                        filelist.append(ds)
                except (InvalidDicomError, AttributeError, MemoryError):
                    pass
            if filelist:
                return filelist
        raise FileNotFoundError("No CT images were found within the specified folder.")

    def _get_ct_images_metadata_from_zip(self, zfile):
        """Get the CT image file names from a zip file."""
        allfiles = zfile.namelist()
        filelist = []
        for name in allfiles:
            try:
                ds = dicom.read_file(BytesIO(zfile.read(name)), force=True)
                if ds.SOPClassUID.name == 'CT Image Storage':
                    filelist.append(ds)
            except (InvalidDicomError, AttributeError):
                pass
        if filelist:
            return filelist
        raise FileNotFoundError("No CT images were found within the specified folder.")

    def _check_all_from_same_study(self, metadatalist):
        initial_uid = metadatalist[0].SeriesInstanceUID
        if not all(ds.SeriesInstanceUID == initial_uid for ds in metadatalist):
            raise ValueError("The images were not all from the same study")

    def _create_array(self, metadatalist, dtype):
        # create empty array the size of the images
        image_shape = metadatalist[0].pixel_array.shape
        self.array = np.zeros((image_shape[0], image_shape[1], len(metadatalist)), dtype=dtype)

        # get the original image order
        original_img_order = [ds.ImagePositionPatient[-1] for ds in metadatalist]

        # set the image array according to the sorted order
        for new, old in enumerate(np.argsort(original_img_order)):
            self.array[:, :, new] = metadatalist[old].pixel_array

        # convert values to proper HU
        self.array *= int(self.metadata.RescaleSlope)
        self.array += int(self.metadata.RescaleIntercept)

    @classmethod
    def from_zip(cls, zip_path, dtype=int):
        if isinstance(zip_path, zipfile.ZipFile):
            zfiles = zip_path
        elif zipfile.is_zipfile(zip_path):
            zfiles = zipfile.ZipFile(zip_path)
        else:
            raise FileExistsError("File given was not a valid zip file")

        obj = cls()
        obj._instantiate(data=zfiles, dtype=dtype, is_zip=True)
        return obj

    def plot(self, slice=0):
        """Plot a slice of the DICOM dataset."""
        plt.imshow(self.slice(slice))

    def slice(self, slice=0):
        """Return an array for the given slice."""
        return self.array[:, :, slice]

    @property
    def shape(self):
        return self.array.shape


class Image:
    """A swiss-army knife, delegate class for loading in images and image-like things.

    The class should not be instantiated directly, but through its class methods. These methods
    return not an `Image` class but one of three specialized image classes:

    * `~pylinac.core.image.DicomImage` : Handles all DICOM images; utilizes pydicom.
    * `~pylinac.core.image.FileImage` : Handles JPEG, BMP, TIF, and other "regular" image files; utilizes Pillow.
    * `~pylinac.core.image.ArrayImage` : Handles 2D numpy arrays; convenient for doing processing of arrays
      that represent an image.

    There are two methods to construct these classes:

    * `.load()` : For loading single images/arrays.
    * `.load_multiples()` : For loading and superimposing multiple images/arrays. All the images must be the same size.

    Examples
    --------
    Load an image from a file::

        >>> my_image = "C:\QA\image.tif"
        >>> img = Image.load(my_image)  # returns a FileImage
        >>> img.median_filter(5)

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
            See `~pylinac.core.image.FileImage` or `~pylinac.core.image.ArrayImage` for keyword arguments.

        Returns
        -------
        `~pylinac.core.image.FileImage`, `~pylinac.core.image.ArrayImage`,
        or `~pylinac.core.image.DicomImage` instance, depending on file type.
        """
        if cls._is_array(path):
            return ArrayImage(path, **kwargs)
        elif cls._is_dicom(path):
            return DicomImage(path)
        elif cls._is_image_file(path):
            return FileImage(path, **kwargs)
        else:
            raise TypeError("The argument `{}` was not found to be a valid DICOM file, Image file, or array".format(path))

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

    @classmethod
    @value_accept(method=('mean', 'max', 'sum'))
    def load_multiples(cls, image_file_list, method='mean', stretch=True):
        """Combine multiple image files into one superimposed image.

        .. versionadded:: 0.5.1
        """
        # load images
        img_list = [cls.load(path) for path in image_file_list]
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


class BaseImage:
    """Base class for the Image classes.

    Attributes
    ----------
    array : numpy.ndarray
        The actual image pixel array.
    center : geometry.Point
        The center pixel of the image as a Point.
    """
    array = typed_property('array', np.ndarray)

    def __init__(self, path):
        if isinstance(path, BytesIO):
            path.seek(0)
        elif not osp.isfile(path):
            raise FileExistsError("File `{}` does not exist".format(path))
        else:
            self.filename = path

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

    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image.

        Wrapper for scipy's `median <http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html>`_ filter function:
        """
        if isinstance(size, float):
            if size < 1:
                size = max(int(self.array.shape[0] * size), 1)
            else:
                raise ValueError("If size is a float, it must be <1.0")
        self.array = ndimage.median_filter(self.array, size=size, mode=mode)

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

        # def rotate(self, angle, order=3):
        #     raise NotImplementedError()
        # self.array = ndimage.interpolation.rotate(self.array, angle, order=order, mode='wrap', reshape=False)

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

    def threshold(self, threshold):
        """Use a high-pass threshold on the array.

        Parameters
        ----------
        threshold : int
            If the value is less than the threshold it is set to 0, otherwise the original value is left as-is.
        """
        self.array = np.where(self.array >= threshold, self, 0)

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

    def check_inversion(self):
        """Check the image for inversion by sampling the 4 image corners.
        If the average value of the four corners is above the average pixel value, then it is very likely inverted.
        """
        outer_edge = 10
        inner_edge = 30
        TL_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BL_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        TR_corner = self.array[outer_edge:inner_edge, outer_edge:inner_edge]
        BR_corner = self.array[-inner_edge:-outer_edge, -inner_edge:-outer_edge]
        corner_avg = np.mean((TL_corner, BL_corner, TR_corner, BR_corner))
        if corner_avg > np.mean(self.array.flatten()):
            self.invert()

    @property
    def shape(self):
        return self.array.shape

    @property
    def size(self):
        return self.array.size

    @property
    def ndim(self):
        return self.array.ndim

    def sum(self):
        return self.array.sum()

    def __len__(self):
        return len(self.array)

    def __getitem__(self, item):
        return self.array[item]


class DicomImage(BaseImage):
    """An image from a DICOM RTImage file.

    Attributes
    ----------
    dicom_dataset : pydicom Dataset
        The dataset of the file as returned by pydicom.
    """

    def __init__(self, path):
        """
        Parameters
        ----------
        path : str, file-object
            The path to the file or the data stream.
        """
        super().__init__(path)
        self.dicom_dataset = dicom.read_file(path, force=True)
        self.array = self.dicom_dataset.pixel_array

    @property
    def sid(self):
        """The Source-to-Image in mm."""
        return float(self.dicom_dataset.RTImageSID)

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
            dpmm = 1 / self.dicom_dataset.PixelSpacing[0]
        except AttributeError:
            try:
                # EPID images sometimes have this tag
                dpmm = 1 / self.dicom_dataset.ImagePlanePixelSpacing[0]
            except AttributeError:
                raise ("No pixel/distance conversion tag found")
        dpmm *= self.sid / 1000
        return dpmm

    @property
    def cax(self):
        """The position of the beam central axis. If no DICOM translation tags are found then the center is returned."""
        try:
            x = self.center.x - self.dicom_dataset.XRayImageReceptorTranslation[0]
            y = self.center.y - self.dicom_dataset.XRayImageReceptorTranslation[1]
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
        return self.dpi / MM_PER_INCH


class ArrayImage(BaseImage):
    """An image constructed solely from a numpy array."""

    def __init__(self, array, *, dpi=None, sid=1000):
        self.array = array
        self._dpi = dpi
        self.sid = sid

    @property
    def dpmm(self):
        return self.dpi / MM_PER_INCH

    @property
    def dpi(self):
        if self._dpi is not None:
            return self._dpi
        else:
            raise AttributeError("No pixel/distance conversion value; pass one in during construction")
