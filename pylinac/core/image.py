
"""Several classes used in pylinac derive from base classes which provide lower-level functionality like image loading, basic image
manipulation (invert, rotate, etc), and other misc things. These are
the API docs for those classes. """
import copy
import warnings

import numpy as np
from scipy import ndimage
from scipy.misc import imresize
from PIL import Image
import dicom

from pylinac.core.decorators import type_accept
from pylinac.core.geometry import Point
from pylinac.core.io import get_filepath_UI, is_valid_file
from pylinac.core.utilities import array2logical


MM_per_INCH = 25.4

class ImageObj:
    """An analysis module component that utilizes a single image in its analysis.
        Contains methods to load and manipulate the image and its properties.
        """
    # TODO: see about inheriting from PIL's Image class.

    def __init__(self, image_array=None):
        """Initialize some attributes."""
        if image_array is None:
            self.pixel_array = np.array([])
        else:
            self.pixel_array = image_array

        self.dpi = 0
        self.dpmm = 0
        self.SID = 0

    @property
    def center(self):
        x_center = self.pixel_array.shape[1] / 2
        y_center = self.pixel_array.shape[0] / 2
        return Point(x_center, y_center)

    @property
    def pixel_array(self):
        return self._pixel_array

    @pixel_array.setter
    def pixel_array(self, array):
        self.set_pixel_array(array)

    @type_accept(image_array=np.ndarray)
    def set_pixel_array(self, image_array=None):
        """Set the image from a pre-existing numpy array"""
        self._pixel_array = image_array

    def load_image_UI(self, caption='', to_gray=True):
        """Load an image using a UI dialog."""
        file_path = get_filepath_UI()
        self.load_image(file_path)

    def load_image(self, file_path, to_gray=True, apply_filter=False, return_it=False):
        """Load an image using PIL or pydicom as a numpy array from a filestring.

        :param file_path: Specifies the file location
        :type file_path: str
        :param to_gray: Indicates whether to convert the image to black & white or not if an RGB image
        :type to_gray: bool
        :param return_it: Will *return* the image and improps if True, otherwise, will save as attr.
            Useful for loading images that won't be used (e.g. for use with combine_images).
        :type return_it: bool
        """

        # Check that filestring points to valid file
        is_valid_file(file_path, raise_error=True)

        # Read image depending on file type
        im_file, image = self._return_img_file(file_path, to_gray)
        # Read in image properties
        self._set_im_props(im_file)

        if apply_filter:
            image = self.median_filter()

        self.pixel_array = image

    def _return_img_file(self, filestring, to_gray=True):
        """Return the file and image, depending on if it's a normal image type (JPG, PNG, etc) or DICOM."""
        # TODO: try incorporating the DICOM SOP; http://www.dicomlibrary.com/dicom/sop/
        try: # try loading dicom first
            im_file = dicom.read_file(filestring)
            image = im_file.pixel_array
            self.set_im_type('DICOM')
        except:  # load as a normal image
            im_file = Image.open(filestring)
            if to_gray:
                im_file = im_file.convert("L")
            image = np.array(im_file)
            self.set_im_type('IMAGE')
        return im_file, image

    def _set_im_props(self, image_file):
        """Return the properties of an image file."""
        if self.im_type == 'DICOM':
            try:
                self.SID = float(image_file.RTImageSID) / 10
            except:
                self.SID = 100
                warnings.warn("Source-to-Image not determined; assuming 100 cm")
            try:
                pixel_spacing = float(image_file.ImagePlanePixelSpacing[0])
                self.dpmm = 1 / (pixel_spacing * self.SID / 1000)
            except AttributeError:
                pixel_spacing = float(image_file.PixelSpacing[0])
                self.dpmm = 1 / pixel_spacing
            except:
                warnings.warn("DICOM Image pixel spacing not set")
        elif self.im_type == 'IMAGE':
            try:
                dpi = image_file.info['dpi']
                if len(dpi) > 1:
                    # ensure all values are the same, i.e. x-resolution is the same as y-resolution
                    if dpi[0] != dpi[1]:
                        raise ValueError("Image DPI is not equal in both directions")
                    dpi = dpi[0]
                    self.dpi = dpi
            except:  # DPI unable to be determined
                self.dpi = 0
                warnings.warn("Image DPI unable to be determined")

    def _combine_images(self, normalize_maximums=True, *images):
        """Combine multiple images together into one image.

        :param normalize_maximums: Specifies whether to normalize the images so that they have the same maximum value. Good for
            images of much different magnitudes.
        :type normalize_maximums: bool
        """
        #TODO: work on this
        raise NotImplementedError("Combine images has not yet been implemented")

    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image.

        Wrapper for scipy's median filter function:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html
        """
        self.pixel_array = ndimage.median_filter(self.pixel_array, size=size, mode=mode)

    def remove_edges(self, pixels=15):
        """Removes the edge pixels on all sides of the image.

        :param pixels: Number of pixels to cut off all sides of the image.
        :type pixels: int
        """
        self.pixel_array = self.pixel_array[pixels - 1:-pixels, pixels - 1:-pixels]

    def invert_array(self):
        """Return the imcomplement of the image.

        Equivalent to Matlab's imcomplement function.
        """
        self.pixel_array = -self.pixel_array + np.max(self.pixel_array) + np.min(self.pixel_array)

    def rotate_array_ccw90(self, n=1):
        """Rotate the image counter-clockwise by 90 degrees n times.

        :param n: Number of times to rotate the image.
        :type n: int
        """
        self.pixel_array = np.rot90(self.pixel_array, n)

    def resize_image(self, size, interp='bilinear'):
        """Resize/scale the image.

        Wrapper for scipy.misc.imresize:
        http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.misc.imresize.html

        :param size: If int, returns % of current size. If float, fraction of current size. If tuple, output size.
        :type size: float, int, tuple
        """
        self.pixel_array = imresize(self.pixel_array, size=size, interp=interp, mode='F')

    def convert2BW(self, threshold, return_it=False):
        """Convert the pixel array to a black & white array based on the threshold."""
        array = array2logical(self.pixel_array, threshold)
        if return_it:
            new_self = copy.copy(self)
            new_self.pixel_array = array
            return new_self
        else:
            self.pixel_array = array

    def dist2edge_min(self, point):
        """Calculates minimum distance from user point to image edges

        :param point: The point to calculate from
        :type point: Point
        """
        rows = np.size(self.pixel_array, 0)
        cols = np.size(self.pixel_array, 1)
        disttoedge = np.zeros(4)
        disttoedge[0] = rows - point.y
        disttoedge[1] = cols - point.x
        disttoedge[2] = point.y
        disttoedge[3] = point.x
        return min(disttoedge)

    @property
    def dpi(self):
        return self._dpi

    @dpi.setter
    def dpi(self, dpi):
        self.set_dpi(dpi)

    def set_dpi(self, dpi):
        """Set the dots-per-inch attribute directly."""
        self._dpi = dpi
        # update the DPmm attr if need be
        updated_dpmm = dpi / MM_per_INCH
        if hasattr(self, 'dpmm') and self.dpmm != updated_dpmm:
            self.dpmm = updated_dpmm

    @property
    def dpmm(self):
        return self._dpmm

    @dpmm.setter
    def dpmm(self, dpmm):
        self.set_dpmm(dpmm)

    def set_dpmm(self, dpmm):
        """Set the dots-per-mm attr directly."""
        self._dpmm = dpmm
        # update the DPI attr if need be
        updated_dpi = dpmm * MM_per_INCH
        if hasattr(self, 'dpi') and self.dpi != updated_dpi:
            self.dpi = updated_dpi

    @property
    def SID(self):
        return self._SID

    @SID.setter
    def SID(self, SID):
        self.set_SID(SID)

    def set_SID(self, SID):
        """Set the SID; units are in mm."""
        self._SID = SID

    @property
    def im_type(self):
        return self._im_type

    @im_type.setter
    def im_type(self, im_type):
        self.set_im_type(im_type)

    @type_accept(im_type=str)
    def set_im_type(self, im_type):
        """Set the image type"""
        self._im_type = im_type


class MultiImageObject:
    """A class to be inherited for multiple image analysis (e.g. CBCT QA)"""
    pass
