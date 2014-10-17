
from __future__ import print_function, division
from future.builtins import str
from future import standard_library
standard_library.install_hooks()
from future.builtins import object

import numpy as np
import scipy.ndimage.filters as spfilt
from scipy.misc import imresize
from PIL import Image
import dicom

from pylinac import running_py3, has_pyqt

if has_pyqt:
    from PyQt4 import QtGui
else:
    if running_py3:
        from tkinter import Tk
        from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
    else:
        from Tkinter import Tk
        from tkFileDialog import askopenfilename, askopenfilenames, askdirectory



"""This holds the basic classes to be inhereted by single-image tools (e.g. starshot, striptest, etc) for simple methods like loading,
rotating the image, and other things."""


class SingleImageObject(object):
    """A class to be inherited by classes that utilize a single image in its analysis. Contains simple methods for such a class."""

    def __init__(self, image=None):
        """Initialize some attributes."""
        self.image = image  # if None, will eventually be a numpy array
        self.im_props = {'DPI': 0,  # Dots (pixels) per inch
                         'DPmm': 0,  # Dots (pixels) per mm
                         'SID mm': None,  # Source (linac target) to Image distance in mm
                         'Image Type': '',  # Image type; either 'DICOM' or 'IMAGE'
        }

    def set_image(self, image):
        """Set the image from a pre-existing numpy array"""
        self.image = image

    def load_demo_image(self):
        """To be overloaded by each specific tool. Loads a demo image for the given class."""
        pass

    def load_image(self, filestring, to_gray=True, return_it=False, apply_filter=False):
        """Load an image using PIL or pydicom as a numpy array from a filestring.
        :param filestring: str; string specifying the file location
        :param to_gray: boolean; indicates whether to convert the image to black & white or not if an RGB image
        :param return_it: boolean; will *return* the image and improps if True, otherwise, will save as attr.
            Useful for loading images that won't be used (e.g. for use with combine images).
        """

        # Read image depending on file type
        im_file, image = self.load_img_file(filestring)
        # Read in image properties
        im_props = self.load_im_props(im_file)

        if apply_filter:
            image = self.median_filter()

        if return_it:
            return image, im_props
        else:
            self.image = image
            self.im_props = im_props

    def load_image_UI(self, dir=None, caption=None, filters=None, to_gray=True, return_it=False):
        """Load an image using a UI"""

        fs = self.get_imagepath_UI(dir=dir,caption=caption, filters=filters)

        if fs:  # if user didn't hit cancel
            if return_it:
                img, props = self.load_image(fs, to_gray=to_gray, return_it=return_it)
                return img, props
            else:
                self.load_image(fs, to_gray=to_gray, return_it=return_it)

    def get_imagepath_UI(self, dir=None, caption=None, filters=None):
        """Return the path of the file chosen with the UI as a string."""

        app = check_app_running()
        if has_pyqt:
            filestring = str(QtGui.QFileDialog.getOpenFileName())
        else:
            filestring = askopenfilename()
        return filestring

    def load_img_file(self, filestring):
        """Return the file and image, depending on if it's a normal image type (JPG, PNG, etc) or DICOM."""
        try: # try loading dicom first
            im_file = dicom.read_file(filestring)
            image = im_file.pixel_array
            self.im_props['Image Type'] = 'DICOM'
            pass
        except:  # load as a normal image
            im_file = Image.open(filestring)
            image = np.array(im_file)
            self.im_props['Image Type'] = 'IMAGE'
        return im_file, image

    def load_im_props(self, image_file):
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

    def combine_images(self, normalize_maximums=True, *images):
        """Combine multiple images together into one image.
        :param normalize_maximums: boolean; specifies whether to normalize the images so that they have the same maximum value. Good for
            images of much different magnitudes.
        """
        #TODO: work on this
        pass

    def analyze(self):
        """To be overloaded by subclass."""
        pass

    def median_filter(self, size=3, mode='reflect'):
        """Apply a median filter to the image. See scipy function for more input information.
        (http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.median_filter.html)
        """
        self.image = spfilt.median_filter(self.image, size=size, mode=mode)

    def remove_edges(self, pixels=15):
        """
        Removes the edge pixels on all sides of the image by the pixel value amount.
        """
        self.image = self.image[pixels - 1:-pixels, pixels - 1:-pixels]

    def invert_image(self):
        """
        Return the imcomplement of the image. Equivalent to Matlab's imcomplement function.
        """
        self.image = -self.image + np.max(self.image) + np.min(self.image)

    def rotate_image_ccw90(self, n=1):
        """Rotate the image counter-clockwise by 90 degrees n times.
        :param n: int; how many times to rotate the image
        """
        self.image = np.rot90(self.image, n)

    def resize_image(self, size, interp='bilinear'):
        """Scale the image. See scipy.misc.pilutil.imresize for further parameter options.
        :param size: int, float, or tuple; if int, returns % of current size. If float, fraction of current size. If tuple, output size.
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
        if has_pyqt:
            folderstring = str(QtGui.QFileDialog.getExistingDirectory())
        else:
            folderstring = askdirectory()
        return folderstring

    def get_filenames_UI(UIdir=None, UIcaption='', UIfilters=''):
        """
        Custom function that is equivalent to Matlab's uigetfile command.
        :return: str; filenames
        """
        # if a QApplication isn't running turn one on; necessary to have one running to use QFileDialog()
        app = check_app_running()
        if has_pyqt:
            filenamesqt = (QtGui.QFileDialog.getOpenFileNames(caption=UIcaption))
            filenames = [str(filename) for filename in filenamesqt]  # convert the PyQt string list to a list of standard strings
        else:
            #TODO: update return to give list, or something other than a string: http://stackoverflow.com/questions/16790328/open-multiple-filenames-in-tkinter-and-add-the-filesnames-to-a-list
            filenames = askopenfilenames()
        return filenames

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
    Opens a QtGui Application or a base Tkinter window if need be; necessary for using things like QFileDialog/askopen*
    """
    if has_pyqt:
        if QtGui.QApplication.instance() is None:
            return QtGui.QApplication([])
    else:
        Tk().withdraw()