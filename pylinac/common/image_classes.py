
"""Several classes used in pylinac derive from base classes which provide lower-level functionality like image loading, basic image
manipulation (invert, rotate, etc), and other misc things. These are
the API docs for those classes. """

from __future__ import print_function, division, absolute_import, unicode_literals
import os.path as osp

try:
    # Python 3.x
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
except ImportError:
    # Python 2.x
    from Tkinter import Tk
    from tkFileDialog import askopenfilename, askopenfilenames, askdirectory

from future.builtins import object

import numpy as np
from scipy import ndimage
from scipy.misc import imresize
from PIL import Image
import dicom

from pylinac.common.decorators import type_accept

class SingleImageObject(object):
    """A class to be inherited by classes that utilize a single image in its analysis. Contains simple methods for such a class."""

    @type_accept(image=(None, np.ndarray))
    def __init__(self, image=None):
        """Initialize some attributes."""
        self.image = image  # if None, will eventually be a numpy array
        self.im_props = {'DPI': 0,  # Dots (pixels) per inch
                         'DPmm': 0,  # Dots (pixels) per mm
                         'SID mm': None,  # Source (linac target) to Image distance in mm
                         'Image Type': '',  # Image type; either 'DICOM' or 'IMAGE'
        }

    @type_accept(image=np.ndarray)
    def set_image(self, image):
        """Set the image from a pre-existing numpy array"""
        self.image = image

    def load_demo_image(self):
        """To be overloaded by each specific tool. Loads a demo image for the given class."""
        pass

    @type_accept(filestring=str, to_gray=bool, return_it=bool, apply_filter=bool)
    def load_image(self, filestring, to_gray=True, return_it=False, apply_filter=False):
        """Load an image using PIL or pydicom as a numpy array from a filestring.

        :param filestring: Specifies the file location
        :type filestring: str
        :param to_gray: Indicates whether to convert the image to black & white or not if an RGB image
        :type to_gray: bool
        :param return_it: Will *return* the image and improps if True, otherwise, will save as attr.
            Useful for loading images that won't be used (e.g. for use with combine_images).
        :type return_it: bool
        """

        # Check that filestring points to valid file
        if not osp.isfile(filestring):
            raise FileExistsError("{} did not point to a valid file".format(filestring))
        # Read image depending on file type
        im_file, image = self._return_img_file(filestring)
        # Read in image properties
        im_props = self._return_im_props(im_file)

        if apply_filter:
            image = self.median_filter()

        if return_it:
            return image, im_props
        else:
            self.image = image
            self.im_props = im_props

    @type_accept(to_gray=bool, return_it=bool)
    def load_image_UI(self, dir='', caption='', filters='', to_gray=True, return_it=False):
        """Load an image using a UI Dialog.

        :param to_gray: Indicates whether to convert the image to black & white or not if an RGB image
        :type to_gray: bool
        :param return_it: Will *return* the image and improps if True, otherwise, will save as attr.
            Useful for loading images that won't be used (e.g. for use with combine_images).
        :type return_it: bool
        """

        fs = self._get_imagepath_UI(dir=dir,caption=caption, filters=filters)

        if fs:  # if user didn't hit cancel
            if return_it:
                img, props = self.load_image(fs, to_gray=to_gray, return_it=return_it)
                return img, props
            else:
                self.load_image(fs, to_gray=to_gray, return_it=return_it)

    def _get_imagepath_UI(self, dir=None, caption=None, filters=None):
        """Return the path of the file chosen with the UI as a string."""

        app = check_app_running()
        filestring = askopenfilename()
        return filestring

    def _return_img_file(self, filestring):
        """Return the file and image, depending on if it's a normal image type (JPG, PNG, etc) or DICOM."""
        # TODO: try incorporating the DICOM SOP; http://www.dicomlibrary.com/dicom/sop/
        try: # try loading dicom first
            im_file = dicom.read_file(filestring)
            image = im_file.pixel_array
            self.im_props['Image Type'] = 'DICOM'
        except:  # load as a normal image
            im_file = Image.open(filestring)
            image = np.array(im_file)
            self.im_props['Image Type'] = 'IMAGE'
        return im_file, image

    def _return_im_props(self, image_file):
        """Return the properties of an image file."""
        im_props = self.im_props
        if self.im_props['Image Type'] == 'DICOM':
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
        pass

    def analyze(self):
        """To be overloaded by subclass."""
        pass

    @type_accept(size=int, mode=str)
    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image. Wrapper for scipy's median filter function:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html
        """
        self.image = ndimage.median_filter(self.image, size=size, mode=mode)

    @type_accept(pixels=int)
    def remove_edges(self, pixels=15):
        """
        Removes the edge pixels on all sides of the image.

        :param pixels: Number of pixels to cut off all sides of the image
        :type pixels: int
        """
        self.image = self.image[pixels - 1:-pixels, pixels - 1:-pixels]

    def invert_image(self):
        """
        Return the imcomplement of the image. Equivalent to Matlab's imcomplement function.
        """
        self.image = -self.image + np.max(self.image) + np.min(self.image)

    @type_accept(n=int)
    def rotate_image_ccw90(self, n=1):
        """Rotate the image counter-clockwise by 90 degrees n times.

        :param n: Number of times to rotate the image.
        :type n: int
        """
        self.image = np.rot90(self.image, n)

    @type_accept(size=(float, int, tuple), interp=str)
    def resize_image(self, size, interp='bilinear'):
        """Scale the image. Wrapper for scipy.misc.imresize: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.misc.imresize.html

        :param size: If int, returns % of current size. If float, fraction of current size. If tuple, output size.
        :type size: float, int, tuple
        """
        self.image = imresize(self.image, size=size, interp=interp, mode='F')
        # self.image = spint.zoom(self.image, zoom=size, mode='nearest')

    def set_dpi(self, dpi):
        """Set the dots-per-inch attribute directly."""
        self.im_props['DPI'] = dpi

    def set_dpmm(self, dpmm):
        """Set the dots-per-mm attr directly."""
        self.im_props['DPmm'] = dpmm


class MultiImageObject(object):
    """A class to be inherited for multiple image analysis (e.g. CBCT QA)"""

    def __init__(self):
        self.images = []  # the image attr is a LIST, not numpy array (allows different size images; more flexible)
        # self._using_pyqa = using_pyqa  # boolean describing whether tool is being used in the PyQA GUI.

    def get_folder_UI(self):
        """Return the string of the location of the folder using a UI."""

        app = check_app_running()
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

    def load_demo_images(self):
        """To be overloaded by subclass."""
        pass

    def load_folder_UI(self):
        """Load the images from a folder using a UI."""
        pass


def check_app_running():
    """
    Opens a base Tkinter window if need be
    """
    Tk().withdraw()