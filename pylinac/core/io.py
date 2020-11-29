"""I/O helper functions for pylinac."""
import os
import os.path as osp
import struct
from tempfile import TemporaryDirectory
from typing import Callable, List
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve, urlopen
import zipfile

import pydicom
from tqdm import tqdm


def is_dicom(file: str) -> bool:
    """Boolean specifying if file is a proper DICOM file.

    This function is a pared down version of read_preamble meant for a fast return.
    The file is read for a proper preamble ('DICM'), returning True if so,
    and False otherwise. This is a conservative approach.

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    with open(file, 'rb') as fp:
        fp.read(0x80)
        prefix = fp.read(4)
        return prefix == b"DICM"


def is_dicom_image(file: str) -> bool:
    """Boolean specifying if file is a proper DICOM file with a image

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    result = False
    try:
        img = pydicom.dcmread(file, force=True)
        if 'TransferSyntaxUID' not in img.file_meta:
            img.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        img.pixel_array
        result = True
    except (AttributeError, TypeError, KeyError, struct.error):
        pass
    return result


def retrieve_dicom_file(file: str) -> pydicom.FileDataset:
    """Read and return the DICOM dataset.

    Parameters
    ----------
    file : str
        The path to the file.
    """
    img = pydicom.dcmread(file, force=True)
    if 'TransferSyntaxUID' not in img.file_meta:
        img.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
    return img


def is_zipfile(file: str) -> bool:
    """Wrapper function for detecting if file is a true ZIP archive"""
    return zipfile.is_zipfile(file)


class TemporaryZipDirectory(TemporaryDirectory):
    """Creates a temporary directory that unpacks a ZIP archive."""
    def __init__(self, zfile):
        """
        Parameters
        ----------
        zfile : str
            String that points to a ZIP archive.
        """
        super().__init__()
        zfiles = zipfile.ZipFile(zfile)
        zfiles.extractall(path=self.name)


def retrieve_filenames(directory: str, func: Callable=None, recursive: bool=True, **kwargs) -> List:
    """Retrieve file names in a directory.

    Parameters
    ----------
    directory : str
        The directory to walk over recursively.
    func : function, None
        The function that validates if the file name should be kept.
        If None, no validation will be performed and all file names will be returned.
    recursive : bool
        Whether to search only the root directory.
    kwargs
        Additional arguments passed to the func parameter.
    """
    filenames = []
    if func is None:
        func = lambda x: True
    for pdir, _, files in os.walk(directory):
        for file in files:
            filename = osp.join(pdir, file)
            if func(filename, **kwargs):
                filenames.append(filename)
        if not recursive:
            break
    return filenames


def retrieve_demo_file(url: str) -> str:
    """Retrieve the demo file either by getting it from file or from a URL.

    If the file is already on disk it returns the file name. If the file isn't
    on disk, get the file from the URL and put it at the expected demo file location
    on disk for lazy loading next time.

    Parameters
    ----------
    url : str
        The suffix to the url (location within the S3 bucket) pointing to the demo file.
    """
    true_url = 'https://s3.amazonaws.com/pylinac/' + url
    demo_file = osp.join(osp.dirname(osp.dirname(__file__)), 'demo_files', url)
    if not osp.isfile(demo_file):
        demo_dir = osp.dirname(demo_file)
        if not osp.exists(demo_dir):
            os.makedirs(demo_dir)
        get_url(true_url, destination=demo_file)
    return demo_file


def is_url(url: str) -> bool:
    """Determine whether a given string is a valid URL.

    Parameters
    ----------
    url : str

    Returns
    -------
    bool
    """
    try:
        with urlopen(url) as r:
            return r.status == 200
    except:
        return False


def get_url(url: str, destination: str=None, progress_bar: bool=True) -> str:
    """Download a URL to a local file.

    Parameters
    ----------
    url : str
        The URL to download.
    destination : str, None
        The destination of the file. If None is given the file is saved to a temporary directory.
    progress_bar : bool
        Whether to show a command-line progress bar while downloading.

    Returns
    -------
    filename : str
        The location of the downloaded file.

    Notes
    -----
    Progress bar use/example adapted from tqdm documentation: https://github.com/tqdm/tqdm
    """

    def my_hook(t):
        last_b = [0]

        def inner(b=1, bsize=1, tsize=None):
            if tsize is not None:
                t.total = tsize
            if b > 0:
                t.update((b - last_b[0]) * bsize)
            last_b[0] = b
        return inner

    try:
        if progress_bar:
            with tqdm(unit='B', unit_scale=True, miniters=1, desc=url.split('/')[-1]) as t:
                filename, _ = urlretrieve(url, filename=destination, reporthook=my_hook(t))
        else:
            filename, _ = urlretrieve(url, filename=destination)
    except (HTTPError, URLError, ValueError) as e:
        raise e
    return filename
