"""Pylinac settings"""
from matplotlib.pyplot import cm

# use a string or colormap option. See options here: http://matplotlib.org/examples/color/colormaps_reference.html
DICOM_COLORMAP = cm.gray
ARRAY_COLORMAP = cm.viridis
PATH_TRUNCATION_LENGTH = 80


def get_dicom_cmap():
    """Return the DICOM colormap. Passed to cmap parameter in matplotlib calls."""
    return DICOM_COLORMAP


def get_array_cmap():
    """Return the array colormap. Passed to cmap parameter in matplotlib calls."""
    return ARRAY_COLORMAP
