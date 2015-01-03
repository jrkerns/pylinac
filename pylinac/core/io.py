from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import os.path as osp

from pylinac.core.utilities import withdraw_tkinter


def is_valid_file(file_path, raise_error=False):
    """Check if path points to a valid file."""
    if osp.isfile(file_path):
        return True
    elif not raise_error:
        return False
    else:
        raise FileExistsError("Path does not point to valid file")

# def get_image_file(file_path):
#     # Check that file_path points to valid file
#     if not osp.isfile(file_path):
#         raise FileExistsError("{} did not point to a valid file".format(file_path))
#     # try to read image
#     # ...try to read as DICOM (most common type in Pylinac)
#     try:
#         imfile = dicom.read_file(file_path)
#     # ...otherwise try as a regular image file
#     except:
#         imfile = Image.open(file_path)
#     return imfile

# def get_image(file_path, to_gray=True):
#     """Return an image (DICOM or anything that PIL can handle) as a numpy array.
#
#     :param to_gray: Specifies whether to return a simple numpy array of image or image object (with properties)
#     :type to_gray: boolean
#     """
#     raw_image = get_image_file(file_path)
#     try:
#         image = raw_image.pixel_array
#     except:
#         image = np.asarray(raw_image)
#
#     return image

# def get_image_properties(image_object):
#     pass



# def get_image(filestring, to_gray=True):
#     """Return a valid image (DICOM or anything that PIL can handle).
#
#     :param to_gray: Specifies whether to return a simple numpy array of image or image object (with properties)
#     :type to_gray: boolean
#     """
#
#     # Check that filestring points to valid file
#     if not osp.isfile(file_path):
#         raise FileExistsError("{} did not point to a valid file".format(file_path))
#
#     # create empty dictionary to put image info into
#     improps = {'DPI': None, 'SID': None, 'Image Type': ''}
#
#     # read image depending on file type
#     if not filestring.endswith('dcm'):
#         imfile = Image.open(filestring)
#         image = np.array(imfile)
#         try:
#             dpi = imfile.info['dpi']
#             if len(dpi) > 1:
#                 # ensure all values are the same, i.e. x-resolution is the same as y-resolution
#                 if dpi[0] != dpi[1]:
#                     raise ValueError("Image DPI is not equal in both directions")
#                 dpi = dpi[0]
#         except:
#             # DPI unable to be determined
#             pass
#         finally:
#             improps['DPI'] = dpi
#     else:  # if dicom file, use pydicom to import
#         imfile = dicom.read_file(filestring)
#         image = imfile.pixel_array
#         improps['Image Type'] = 'DICOM'
#         improps['SID'] = float(imfile.RTImageSID)
#         try:
#             pixel_spacing = float(imfile.ImagePlanePixelSpacing[0])
#             improps['DPI'] = pixel_spacing * improps['SID'] / 1000
#         except:
#             pass
#
#     return image, improps


# def get_image_UI(UIdir=None, UIcaption='', UIfilters='', to_gray=True):  # TODO: update args to fit getopenfilename's args
#     """Return an image using a UI Dialog."""
#
#     filestring = get_filepath_UI(UIdir=UIdir, UIcaption=UIcaption, UIfilters=UIfilters)
#
#     if filestring:
#         filestring = str(filestring)
#         image = get_image(filestring, to_gray)
#         return image

def get_filepath_UI(dir=None, caption='', filters=''):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    # if a QApplication isn't running turn one on; necessary to have one running to use QFileDialog()
    withdraw_tkinter()

    # get user-defined image file
    filename = askopenfilename()
    return filename

def get_filenames_UI(UIdir=None, UIcaption='', UIfilters=''):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    withdraw_tkinter()
    # get user-defined files
    filenames = askopenfilenames(caption=UIcaption)
    filenames = [str(f) for f in filenames]  # convert the PyQt string list to a list of standard strings

    if filenames:
        return filenames

def get_folder_UI(UIdir=None, UIcaption=''):
    """
    :returns: string of Folder the user selected
    """
    withdraw_tkinter()
    folderstring = askdirectory()
    return folderstring
