
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


class AnalysisModule(object):
    """An abstract class for distinct analysis modules (Starshot, VMAT, etc). Its purpose is to define basic method
    names to be consistent across analysis modules. E.g. the starshot and VMAT modules will both have an "analyze" and
    "load_image" method.
    """
    def __init__(self):
        self.test_passed = False  # set initial test result to fail

    def load_demo_image(self):
        """To be overloaded by each specific tool. Loads a demo image for the given class."""
        raise NotImplementedError("Loading the demo image(s) for this module has not been implemented yet.")

    def analyze(self):
        """To be overloaded by subclass. Main analysis method for module"""
        raise NotImplementedError("Analysis has not been implemented for this module.")

    def run_demo(self):
        """Demo of module's abilities."""
        raise NotImplementedError("The demo for this module has not been built yet.")

class ImageObj(object):
    """An analysis module component that utilizes a single image in its analysis.
    Contains methods to load and manipulate the image and its properties.
    """

    def __init__(self, image_array=None):
        """Initialize some attributes."""
        self.pixel_array = image_array  # if None, will eventually be a numpy array
        self.properties = {'DPI': None,  # Dots (pixels) per inch
                         'DPmm': None,  # Dots (pixels) per mm
                         'SID mm': None,  # Source (linac target) to Image distance in mm
                         'Image Type': '',  # Image type; either 'DICOM' or 'IMAGE'
        }

    def set_pixel_array(self, image_array):
        """Set the image from a pre-existing numpy array"""
        self.pixel_array = image_array

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
            # Python 2 doesn't have FileExistsError
            try:
                raise FileExistsError("{} did not point to a valid file".format(file_path))
            except:
                raise IOError("{} did not point to a valid file".format(file_path))
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
        """Apply a median filter to the image. Wrapper for scipy's median filter function:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html
        """
        self.pixel_array = ndimage.median_filter(self.pixel_array, size=size, mode=mode)

    def remove_edges(self, pixels=15):
        """Removes the edge pixels on all sides of the image.

        :param pixels: Number of pixels to cut off all sides of the image
        :type pixels: int
        """
        self.pixel_array = self.pixel_array[pixels - 1:-pixels, pixels - 1:-pixels]

    def invert_array(self):
        """
        Return the imcomplement of the image. Equivalent to Matlab's imcomplement function.
        """
        self.pixel_array = -self.pixel_array + np.max(self.pixel_array) + np.min(self.pixel_array)

    def rotate_array_ccw90(self, n=1):
        """Rotate the image counter-clockwise by 90 degrees n times.

        :param n: Number of times to rotate the image.
        :type n: int
        """
        self.pixel_array = np.rot90(self.pixel_array, n)

    def resize_image(self, size, interp='bilinear'):
        """Scale the image. Wrapper for scipy.misc.imresize: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.misc.imresize.html

        :param size: If int, returns % of current size. If float, fraction of current size. If tuple, output size.
        :type size: float, int, tuple
        """
        self.pixel_array = imresize(self.pixel_array, size=size, interp=interp, mode='F')
        # self.image = spint.zoom(self.image, zoom=size, mode='nearest')

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