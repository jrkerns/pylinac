
"""Several classes used in pylinac derive from base classes which provide lower-level functionality like image loading, basic image
manipulation (invert, rotate, etc), and other misc things. These are
the API docs for those classes. """

import os.path as osp
from tkinter.filedialog import askdirectory

import numpy as np
from scipy import ndimage
from scipy.misc import imresize
from PIL import Image
import dicom

from pylinac.core.geometry import Point
from pylinac.core.utilities import withdraw_tkinter


class ImageObj(object):
    """An analysis module component that utilizes a single image in its analysis.
    Contains methods to load and manipulate the image and its properties.
    """

    def __init__(self, image_array=None):
        """Initialize some attributes."""
        self.pixel_array = image_array
        self.properties = {'DPI': None,  # Dots (pixels) per inch
                           'DPmm': None,  # Dots (pixels) per mm
                           'SID mm': None,  # Source (linac target) to Image distance in mm
                           'Image Type': '',  # Image type; either 'DICOM' or 'IMAGE'
        }

    @property
    def dpi(self):
        pass



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

    def set_pixel_array(self, image_array=None):
        """Set the image from a pre-existing numpy array"""
        if image_array is None:
            self._pixel_array = np.array([])
        elif not isinstance(image_array, np.ndarray):
            raise TypeError("Image array must be a numpy array")
        else:
            self._pixel_array = image_array

    def load_image(self, file_path, to_gray=True, return_it=False, apply_filter=False):
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
        if not osp.isfile(file_path):
            raise FileExistsError("{} did not point to a valid file".format(file_path))
        # Read image depending on file type
        im_file, image = self._return_img_file(file_path)
        # Read in image properties
        im_props = self._return_im_props(im_file)

        if apply_filter:
            image = self.median_filter()

        if return_it:
            return image, im_props
        else:
            self.pixel_array = image
            self.properties = im_props

    def _return_img_file(self, filestring):
        """Return the file and image, depending on if it's a normal image type (JPG, PNG, etc) or DICOM."""
        # TODO: try incorporating the DICOM SOP; http://www.dicomlibrary.com/dicom/sop/
        try: # try loading dicom first
            im_file = dicom.read_file(filestring)
            image = im_file.pixel_array
            self.properties['Image Type'] = 'DICOM'
        except:  # load as a normal image
            im_file = Image.open(filestring)
            image = np.array(im_file)
            self.properties['Image Type'] = 'IMAGE'
        return im_file, image

    def _return_im_props(self, image_file):
        """Return the properties of an image file."""
        im_props = self.properties
        if self.properties['Image Type'] == 'DICOM':
            try:
                im_props['SID mm'] = float(image_file.RTImageSID)
            except:
                # assume 1000 mm otherwise
                im_props['SID mm'] = 1000
            try:
                pixel_spacing = float(image_file.ImagePlanePixelSpacing[0])
                im_props['DPmm'] = 1 / (pixel_spacing * im_props['SID mm'] / 1000)
                im_props['DPI'] = im_props['DPmm'] * 25.4  # 1 inch = 25.4 mm
            except:
                pass
        else:
            try:
                dpi = image_file.info['dpi']
                if len(dpi) > 1:
                    # ensure all values are the same, i.e. x-resolution is the same as y-resolution
                    if dpi[0] != dpi[1]:
                        raise ValueError("Image DPI is not equal in both directions")
                    dpi = dpi[0]
                    im_props['DPI'] = dpi
                    im_props['DPmm'] = dpi / 25.4
            except:  # DPI unable to be determined
                pass

        return im_props

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

    def set_dpi(self, dpi):
        """Set the dots-per-inch attribute directly."""
        self.properties['DPI'] = dpi

    def set_dpmm(self, dpmm):
        """Set the dots-per-mm attr directly."""
        self.properties['DPmm'] = dpmm


class MultiImageObject(object):
    """A class to be inherited for multiple image analysis (e.g. CBCT QA)"""

    def __init__(self):
        self.images = []  # the image attr is a LIST, not numpy array (allows different size images; more flexible)
        # self._using_pyqa = using_pyqa  # boolean describing whether tool is being used in the PyQA GUI.

    def get_folder_UI(self):
        """Return the string of the location of the folder using a UI."""

        withdraw_tkinter()
        folderstring = askdirectory()
        return folderstring

    # def get_filenames_UI(UIdir=None, UIcaption='', UIfilters=''):
    #     """
    #     Custom function that is equivalent to Matlab's uigetfile command.
    #     :return: str; filenames
    #     """
        # if a QApplication isn't running turn one on; necessary to have one running to use QFileDialog()
        # app = check_app_running()
        # #TODO: update return to give list, or something other than a string: http://stackoverflow.com/questions/16790328/open-multiple-filenames-in-tkinter-and-add-the-filesnames-to-a-list
        # filenames = askopenfilenames()
        # return filenames

    def load_folder(self, filestring, append=False):
        """Load images from a folder.
        :param append: boolean; specifies whether to append images in folder to images attr.
        """
        pass

    def load_folder_UI(self):
        """Load the images from a folder using a UI."""
        pass
