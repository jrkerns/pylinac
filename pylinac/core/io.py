"""I/O helper functions for pylinac."""
import os.path as osp
from tempfile import TemporaryDirectory
from urllib.error import HTTPError
from urllib.request import urlretrieve
import zipfile


class TemporaryZipDirectory(TemporaryDirectory):
    """Creates a temporary directory that unpacks a ZIP archive."""
    def __init__(self, zfile):
        super().__init__()
        zfiles = zipfile.ZipFile(zfile)
        zfiles.extractall(path=self.name)


def get_url(url, destination=None):
    """Download a URL to a local file.

    Parameters
    ----------
    destination : str, None
        The destination of the file. If None is given the file is saved to a temporary directory.
    """
    try:
        filename, _ = urlretrieve(url, filename=destination)
    except HTTPError as e:
        raise e
    return filename


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
            raise FileExistsError("{0} is not a valid file".format(file_path))

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
            raise TypeError("Could not open file: {0}".format(file))
    return openfile

def is_valid_dir(dir_path, raise_error=True):
    """Check if path points to a valid directory."""
    if osp.isdir(dir_path):
        return True
    elif not raise_error:
        return False
    else:
        raise NotADirectoryError("{0} does not point to a valid directory".format(dir_path))
