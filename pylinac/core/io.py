"""I/O helper functions for pylinac."""

from io import BytesIO
import os.path as osp
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import zipfile


def load_zipfile(zfilename, read=False):
    """Read a .zip file.

    Parameters
    ----------
    zfilename : str
        Path to the zip file.
    read : bool
        Whether to read the zip files.

    Returns
    -------
    files
        If read is False (default), returns a python zipfile.ZipFile object. If read is True, returns
        a list of BytesIO objects.
    """
    if isinstance(zfilename, zipfile.ZipFile):
        zfiles = zfilename
    elif zipfile.is_zipfile(zfilename):
        zfiles = zipfile.ZipFile(zfilename)
    else:
        raise FileExistsError("File '{}' given was not a valid zip file".format(zfilename))
    if read:
        return [BytesIO(zfiles.read(name)) for name in zfiles.namelist()]
    else:
        return zfiles


def get_url(url):
    """Get the response from the URL."""
    try:
        import requests
    except ImportError:
        raise ImportError("Requests is not installed; cannot get the log from a URL")
    response = requests.get(url, timeout=10)
    if response.status_code != 200:
        raise ConnectionError("Could not connect to the URL")
    return response


def is_valid_file(file_path, raise_error=True):
    """Check if path points to a valid file.

    Parameters
    ----------
    file_path : str
        Path to the file.
    raise_error : boolean
        If False (default), will simply return false.
        If True, will raise an error if path is not valid.

    Raises
    ------
    FileExistsError
        If file_path does not point to a valid file.
    """
    try:
        file_path.seek(0)
        return True
    except ValueError:
        return True
    except AttributeError:
        if osp.isfile(file_path):
            return True
        elif not raise_error:
            return False
        else:
            raise FileExistsError("{} is not a valid file".format(file_path))

def open_file(file, mode='rb'):
    """Open a file if a file is a string, or open the object if file is a file object.
        If an already-opened file object, reset the pointer to 0."""
    if isinstance(file, str):
        openfile = open(file, mode)
    else:
        try:
            file.open()
            openfile = file
        except (ValueError, AttributeError):
            file.seek(0)
            openfile = file
        except:
            raise TypeError("Could not open file: {}".format(file))
    return openfile

def is_valid_dir(dir_path, raise_error=True):
    """Check if path points to a valid directory."""
    if osp.isdir(dir_path):
        return True
    elif not raise_error:
        return False
    else:
        raise NotADirectoryError("{} does not point to a valid directory".format(dir_path))

def get_filepath_UI(dir=None, caption='', filters=''):
    """Display a UI dialog box to select a file.

    Returns
    -------
    str
        Path to the file chosen.
    """
    withdraw_tkinter()
    filename = askopenfilename()
    return filename

def get_filenames_UI(UIdir=None, UIcaption='', UIfilters=''):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    withdraw_tkinter()
    filenames = askopenfilenames()
    return filenames

def get_folder_UI(UIdir=None, UIcaption=''):
    """Display a UI dialog box to select a folder.

    Returns
    -------
    str
        Path to the folder chosen.
    """
    withdraw_tkinter()
    folderstring = askdirectory()
    return folderstring

def withdraw_tkinter():
    """Opens and withdraws a Tk window. Necessary so a base window doesn't open."""
    Tk().withdraw()