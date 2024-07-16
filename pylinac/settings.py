"""Pylinac settings"""

# use a string or colormap option. See options here: http://matplotlib.org/examples/color/colormaps_reference.html
DICOM_COLORMAP = "gray"
ARRAY_COLORMAP = "viridis"
PATH_TRUNCATION_LENGTH = 80


def get_dicom_cmap() -> str:
    """Return the DICOM colormap. Passed to cmap parameter in matplotlib calls."""
    return DICOM_COLORMAP


def get_array_cmap() -> str:
    """Return the array colormap. Passed to cmap parameter in matplotlib calls."""
    return ARRAY_COLORMAP
