"""Pylinac settings"""
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import cm

# use a string or colormap option. See options here: http://matplotlib.org/examples/color/colormaps_reference.html
DICOM_COLORMAP = cm.gray
ARRAY_COLORMAP = cm.gray
PATH_TRUNCATION_LENGTH = 80


def get_dicom_cmap() -> LinearSegmentedColormap:
    """Return the DICOM colormap. Passed to cmap parameter in matplotlib calls."""
    return DICOM_COLORMAP


def get_array_cmap() -> LinearSegmentedColormap:
    """Return the array colormap. Passed to cmap parameter in matplotlib calls."""
    return ARRAY_COLORMAP
