"""I/O helper functions for pylinac."""
from tempfile import TemporaryDirectory
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve
import zipfile

from tqdm import tqdm


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


def get_url(url, destination=None, progress_bar=True):
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
