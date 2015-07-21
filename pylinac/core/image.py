
"""This module holds classes for image loading and manipulation."""
import warnings

from dicom.errors import InvalidDicomError
from dicom.dataset import Dataset
import numpy as np
from PIL import Image as pImage
import dicom
from scipy import ndimage
from scipy.misc import imresize
import matplotlib.pyplot as plt

from pylinac.core.decorators import type_accept
from pylinac.core.geometry import Point
from pylinac.core.io import get_filepath_UI, get_filenames_UI
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
        The Dots-per-mm of the image, defined at isocenter. E.g. if an EPID image is taken at 150cm SID,
        the dpmm will scale back to 100cm.
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

        >>> img = Image.from_UI()

    Or, load from an existing array::

        >>> arr = np.arange(36).reshape(6,6)
        >>> img = Image(arr)
    """
    # SID = typed_property('SID', (int, float, np.number))
    im_type = typed_property('im_type', str)
    array = typed_property('array', np.ndarray)

    def __init__(self, filename=None):
        """
        Parameters
        ----------
        filename : str
            Path to the image file.
        """
        self._SID = 1000
        if filename is not None:
            try:
                self._load_file(filename)
            except (IOError, AttributeError):
                raise TypeError("Image input type not understood")

    @classmethod
    def from_array(cls, array):
        obj = cls()
        obj.array = array
        obj.im_type = ARRAY
        return obj

    @property
    def dpi(self):
        if self.im_type == DICOM:
            dpi = getattr(self, 'dpmm', None)
            if dpi is not None:
                dpi *= MM_per_INCH
            return dpi
        elif self.im_type == IMAGE:
            try:
                dpi = self._img_meta['dpi'][0]
            except (IndexError, KeyError):
                dpi = self._img_meta.get('dpi', None)
            if dpi is not None:
                dpi *= self.SID / 1000
            return dpi
        else:
            return None

    @property
    def dpmm(self):
        if self.im_type == DICOM:
            try:
                # most dicom files have this tag
                dpmm = 1/self._dcm_meta.PixelSpacing[0]
            except AttributeError:
                try:
                    # EPID images sometimes have this tag
                    dpmm = 1/self._dcm_meta.ImagePlanePixelSpacing[0]
                except AttributeError:
                    dpmm = None
            if dpmm is not None:
                dpmm *= self.SID / 1000
            return dpmm
        elif self.im_type == IMAGE:
            dpmm = self.dpi
            if dpmm is not None:
                dpmm /= MM_per_INCH
            return dpmm
        else:  # Array type
            return None

    @property
    def center(self):
        """Return the center position of the image array as a Point."""
        x_center = self.shape[1] / 2
        y_center = self.shape[0] / 2
        return Point(x_center, y_center)

    @property
    def cax(self):
        """Return the position of the CAX."""
        try:
            x = self.center.x - self._dcm_meta.XRayImageReceptorTranslation[0]
            y = self.center.y - self._dcm_meta.XRayImageReceptorTranslation[1]
        except AttributeError:
            return self.center
        else:
            return Point(x, y)

    @property
    def SID(self):
        """Return the SID."""
        if self.im_type == DICOM:
            sid = getattr(self._dcm_meta, 'RTImageSID', self._SID)
        else:
            sid = self._SID
        return float(sid)

    @SID.setter
    @type_accept(value=(int, float, np.number))
    def SID(self, value):
        """Set the SID value; must be in mm."""
        # if not isinstance(value, (int, float, np.number)):
        #     raise ValueError("SID must be a number")
        if self.im_type == DICOM:
            raise AttributeError("Cannot set the SID for DICOM Images")
        else:
            self._SID = value

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

    @classmethod
    def from_UI(cls, caption='', to_gray=True):
        """Load an image using a UI dialog."""
        file_path = get_filepath_UI()
        if file_path:
            obj = cls(file_path, to_gray)
            return obj

    @classmethod
    def from_multiple_UI(cls, caption='', to_gray=True):
        """Load multiple images using a UI dialog.

        .. versionadded:: 0.5.1

        All files must be images, and must be the same size and shape.
        Image metadata, e.g. DPI, is all based on the first image selected.
        """
        file_list = get_filenames_UI()
        if file_list:
            obj = cls.from_multiples(file_list)
            return obj

    def _load_file(self, file_path):
        """Load a file."""
        try:
            self._construct_dicom(file_path)
        except InvalidDicomError:
            try:
                self._construct_image(file_path)
            except OSError:
                raise IOError("Image type not supported")

    def _construct_image(self, file_path):
        """Construct an object from an image file (TIF, JPEG, etc)."""
        try:
            file_path.seek(0)
        except AttributeError:
            pass
        img = pImage.open(file_path)

        # convert to gray if need be
        if img.mode == 'RGB' or img.mode == 'HSV':
            img = img.convert('F')

        self._img_meta = img.info

        self.array = np.array(img)
        self.im_type = IMAGE

    def _construct_dicom(self, file_path):
        """Construct an object from a DICOM file (.dcm)."""
        try:
            file_path.seek(0)
        except AttributeError:
            pass
        dcm = dicom.read_file(file_path)
        self.array = dcm.pixel_array
        self.im_type = DICOM

        # attach the metadata
        self._dcm_meta = dicom.read_file(file_path, stop_before_pixels=True)

    def plot(self):
        """Plot the image."""
        plt.clf()
        plt.imshow(self.array, cmap=plt.cm.Greys)
        plt.show()

    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image.

        Wrapper for scipy's `median <http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html>`_ filter function:
        """
        if isinstance(size, float):
            if size < 1:
                size = max(int(self.array.shape[0]*size), 1)
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
        self.array = self.array[pixels - 1:-pixels-1, pixels - 1:-pixels-1]

    def invert(self):
        """Invert (imcomplement) the image."""
        orig_array = self.array
        self.array = -orig_array + orig_array.max() + orig_array.min()

    # def rotate(self, angle, order=3):
    #     raise NotImplementedError()
        # self.array = ndimage.interpolation.rotate(self.array, angle, order=order, mode='wrap', reshape=False)

    def rot90(self, n=1):
        """Wrapper for numpy.rot90."""
        self.array = np.rot90(self.array, n)

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
        return Image.from_array(arr)

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

    @classmethod
    def from_multiples(cls, image_file_list):
        """Combine multiple image files into one superimposed image.

        .. versionadded:: 0.5.1
        """
        # open first one to get initial settings
        init_obj = cls(image_file_list[0])
        concat_arr = init_obj.array
        initial_shape = init_obj.shape
        for img_file in image_file_list[1:]:
            obj = cls(img_file)
            if obj.shape != initial_shape:
                 raise AttributeError("Images must be the same size when combining.")
            concat_arr = np.dstack((concat_arr, obj.array))

        # create new mean array
        combined_arr = np.mean(concat_arr, axis=2)
        # use the initial Image object and replace its array (thus keeping all the other properties)
        init_obj.array = combined_arr
        return init_obj

    @property
    def shape(self):
        return self.array.shape

    @property
    def size(self):
        return self.array.size

    # def __getattr__(self, item):
    #     """Set the Attribute getter to grab from the array if possible (for things like .shape, .size, etc)."""
    #     return getattr(self.array, item)

if __name__ == '__main__':
    img = Image.from_multiple_UI()

