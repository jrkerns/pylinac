
"""This module holds classes for image loading and manipulation."""
import warnings

from dicom.errors import InvalidDicomError
import numpy as np
from PIL import Image as pImage
import dicom
from scipy import ndimage
from scipy.misc import imresize

from pylinac.core.decorators import type_accept
from pylinac.core.geometry import Point
from pylinac.core.io import get_filepath_UI
from pylinac.core.utilities import array2logical, typed_property


DICOM = 'DICOM'
IMAGE = 'Image'
ARRAY = 'Array'
MM_per_INCH = 25.4


class Image:
    """A class that holds an image as a numpy array, relevant image metadata (dpi, SID, etc),
    and methods for image (array) manipulation (resize, rotate, etc).

    There are 3 types of Images: DICOM, Image, and Array. For DICOM and Image types,
    relevant metadata (see Attributes) is extracted where possible. Metadata can also be set directly.

    * A "DICOM" image is just what it sounds like. It could be an EPID capture, or a scanned film/CR cassette.
    * An "Image" image is what you think normally of an image. It could be a JPEG, BMP, TIF, etc.
    * An "Array" image is one created directly from an existing array.

    Attributes
    ----------
    array : numpy.ndarray
        The actual image pixel array.
    dpi : int, float
        The Dots-per-inch of the image.
    dpmm : int, float
        The Dots-per-mm of the image.
    SID : int, float
        The Source-to-Image distance in cm.
    im_type : {'DICOM', 'Image', 'Array'}
        Image type.
    center : geometry.Point
        The center pixel of the image as a Point.

    Examples
    --------
    Load an image from a file::

        >>> my_image = "C:/QA/image.tif"
        >>> img = Image(my_image)
        >>> img.median_filter(5)

    Additionally, load from a UI dialog box::

        >>> img = Image.open_UI()

    Or, load from an existing array::

        >>> arr = np.arange(36).reshape(6,6)
        >>> img = Image(arr)
    """
    SID = typed_property('SID', (int, float, np.number))
    im_type = typed_property('im_type', str)
    array = typed_property('array', np.ndarray)

    def __init__(self, file_or_array, to_gray=True):
        """
        Parameters
        ----------
        file_path : str
            Path to the image file.
        to_gray : bool
            If True (default), will convert RGB and HSV type images to greyscale.
            If False, will not do any conversion.
        """
        if isinstance(file_or_array, str):
            self._load_file(file_or_array, to_gray)
        elif isinstance(file_or_array, np.ndarray):
            self._load_array(file_or_array)
        else:
            raise TypeError("Image input type not understood")

    @property
    def dpi(self):
        return getattr(self, '_dpi', None)

    @dpi.setter
    @type_accept(dpi=(int, float, np.number))
    def dpi(self, dpi):
        self._dpi = dpi
        # update the DPmm attr if need be
        updated_dpmm = dpi / MM_per_INCH
        if self.dpmm != updated_dpmm:
            self.dpmm = updated_dpmm

    @property
    def dpmm(self):
        return getattr(self, '_dpmm', None)

    @dpmm.setter
    @type_accept(dpmm=(int, float, np.number))
    def dpmm(self, dpmm):
        self._dpmm = dpmm
        # update the DPI attr if need be
        updated_dpi = dpmm * MM_per_INCH
        if self.dpi != updated_dpi:
            self.dpi = updated_dpi

    @property
    def center(self):
        """Return the center position of the image array as a Point."""
        x_center = self.shape[1] / 2
        y_center = self.shape[0] / 2
        return Point(x_center, y_center)

    @classmethod
    def open_UI(cls, caption='', to_gray=True):
        """Load an image using a UI dialog."""
        file_path = get_filepath_UI()
        cls(file_path, to_gray)
        return cls

    def _load_array(self, array):
        """Load an array."""
        self.array = array
        self.im_type = ARRAY

    def _load_file(self, file_path, to_gray):
        """Load a file."""
        try:
            self._construct_dicom(file_path)
        except InvalidDicomError:
            try:
                self._construct_image(file_path, to_gray)
            except OSError:
                raise IOError("Image type not supported")

    def _construct_image(self, file_path, to_gray):
        """Construct an object from an image file (TIF, JPEG, etc)."""
        img = pImage.open(file_path)
        if to_gray:
            if img.mode == 'RGB' or img.mode == 'HSV':
                img = img.convert('F')
        self.array = np.array(img)
        self.im_type = IMAGE

        # explicitly set dpi, which also sets dpmm
        if 'dpi' in img.info:
            if len(img.info['dpi']) > 1:
                if img.info['dpi'][0] != img.info['dpi'][1]:
                    raise TypeError("Input image has different DPI ratios for horizontal and vertical")
                self.dpi = img.info['dpi'][0]
            else:
                self.dpi = img.info['dpi']
        else:
            warnings.warn("Pixel distance information was not available. You can set the attributed"
                          "explicitly if you know it.")

    def _construct_dicom(self, file_path):
        """Construct an object from a DICOM file (.dcm)."""
        dcm = dicom.read_file(file_path)
        self.array = dcm.pixel_array
        self.im_type = DICOM
        # try to set the pixel spacing
        try:
            # most dicom files have this tag
            dpmm = dcm.PixelSpacing[0]
        except AttributeError:
            # EPID images sometimes have this tag
            dpmm = dcm.ImagePlanePixelSpacing[0]
        except:
            warnings.warn("Pixel distance information for the DICOM file was not determined. You can set the "
                          "attribute (.dpi or .dpmm) explicitly if you know it.")
        finally:
            self.dpmm = float(dpmm)

        # try to set the SID
        try:
            sid = float(dcm.RTImageSID) / 10
            self.SID = sid
        except AttributeError:
            pass  # just don't set the SID

    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image.

        Wrapper for scipy's `median <http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html>`_ filter function:
        """
        self.array = ndimage.median_filter(self.array, size=size, mode=mode)

    @type_accept(pixels=int)
    def remove_edges(self, pixels=15):
        """Removes pixels on all edges of the image.

        Parameters
        ----------
        pixels : int
            Number of pixels to cut off all sides of the image.
        """
        self.array = self.array[pixels - 1:-pixels-1, pixels - 1:-pixels-1]

    def invert(self):
        """Invert (imcomplement) the image."""
        orig_array = self.array
        self.array = -orig_array + orig_array.max() + orig_array.min()

    def rotate(self, angle, order=3):
        raise NotImplementedError()
        # self.array = ndimage.interpolation.rotate(self.array, angle, order=order, mode='wrap', reshape=False)

    def resize(self, size, interp='bilinear'):
        """Resize/scale the image.

        Wrapper for scipy's `imresize <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.misc.imresize.html>`_:
        """
        self.array = imresize(self.array, size=size, interp=interp, mode='F')

    def threshold(self, threshold):
        """Convert the pixel array to a black & white array based on the threshold.

        Parameters
        ----------
        threshold : int
            If the value is less than the threshold it is set to 0, otherwise to 1.

        Returns
        -------
        A numpy array the same size as the original image.
        """
        arr = array2logical(self.array, threshold)
        return Image(arr)

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

    def __getattr__(self, item):
        """Set the Attribute getter to grab from the array if possible (for things like .shape, .size, etc)."""
        return getattr(self.array, item)

