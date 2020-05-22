"""The CatPhan module automatically analyzes DICOM images of a CatPhan 504, 503, or 600 acquired when doing CBCT or CT quality assurance.
It can load a folder or zip file that the images are in and automatically correct for translational and rotational errors.
It can analyze the HU regions and image scaling (CTP404), the high-contrast line pairs (CTP528) to calculate the modulation transfer function (MTF),
the HU uniformity (CTP486), and Low Contrast (CTP515) on the corresponding slices.

Features:

* **Automatic phantom registration** - Your phantom can be tilted, rotated, or translated--pylinac will automatically register the phantom.
* **Automatic testing of all major modules** - Major modules are automatically registered and analyzed.
* **Any scan protocol** - Scan your CatPhan with any protocol; even scan it in a regular CT scanner.
  Any field size or field extent is allowed.
"""
from collections import OrderedDict
from datetime import datetime
from functools import lru_cache
import io
from os import path as osp
import os
import webbrowser
import zipfile
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from skimage import filters, measure, segmentation

from .core import image
from .core.io import TemporaryZipDirectory
from .core.decorators import value_accept
from .core.geometry import Point, Line
from .core.io import get_url, retrieve_demo_file
from .core import pdf
from .core.mtf import MTF
from .core.profile import CollapsedCircleProfile, SingleProfile
from .core.roi import DiskROI, RectangleROI, LowContrastDiskROI
from .settings import get_dicom_cmap


AIR = -1000
PMP = -196
LDPE = -104
POLY = -47
ACRYLIC = 115
DELRIN = 365
TEFLON = 1000
BONE_20 = 237
BONE_50 = 725


class HUDiskROI(DiskROI):
    """An HU ROI object. Represents a circular area measuring either HU sample (Air, Poly, ...)
    or HU uniformity (bottom, left, ...).
    """
    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, nominal_value=None, tolerance=None,
                 background_median=None, background_std=None):
        """
        Parameters
        ----------
        nominal_value : int
            The nominal pixel value of the HU ROI.
        tolerance : int
            The roi pixel value tolerance.
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.nominal_val = nominal_value
        self.tolerance = tolerance
        self.background_median = background_median
        self.background_std = background_std

    @property
    def cnr(self):
        """The contrast-to-noise value of the HU disk"""
        return 2*abs(self.pixel_value - self.background_median) / (self.std + self.background_std)

    @property
    def value_diff(self):
        """The difference in HU between measured and nominal."""
        return self.pixel_value - self.nominal_val

    @property
    def passed(self):
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return abs(self.value_diff) <= self.tolerance

    @property
    def plot_color(self):
        """Return one of two colors depending on if ROI passed."""
        return 'green' if self.passed else 'red'


class ThicknessROI(RectangleROI):
    """A rectangular ROI that measures the angled wire rod in the HU linearity slice which determines slice thickness."""

    @property
    @lru_cache(maxsize=1)
    def long_profile(self):
        """The profile along the axis perpendicular to ramped wire."""
        img = image.load(self.pixel_array)
        img.filter(size=1, kind='gaussian')
        prof = SingleProfile(img.array.max(axis=np.argmin(img.shape)))
        return prof

    @property
    @lru_cache(maxsize=1)
    def wire_fwhm(self):
        """The FWHM of the wire in pixels."""
        return self.long_profile.fwxm(x=50, interpolate=True)

    @property
    def plot_color(self):
        """The plot color."""
        return 'blue'


class Slice:
    """Base class for analyzing specific slices of a CBCT dicom set."""
    def __init__(self, catphan, slice_num=None, combine=True, combine_method='mean', num_slices=0):
        """
        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        slice_num : int
            The slice number of the DICOM array desired. If None, will use the ``slice_num`` property of subclass.
        combine : bool
            If True, combines the slices +/- ``num_slices`` around the slice of interest to improve signal/noise.
        combine_method : {'mean', 'max'}
            How to combine the slices if ``combine`` is True.
        num_slices : int
            The number of slices on either side of the nominal slice to combine to improve signal/noise; only
            applicable if ``combine`` is True.
        """
        if slice_num is not None:
            self.slice_num = slice_num
        if combine:
            array = combine_surrounding_slices(catphan.dicom_stack, self.slice_num, mode=combine_method, slices_plusminus=num_slices)
        else:
            array = catphan.dicom_stack[self.slice_num].array
        self.image = image.load(array)
        self.catphan_size = catphan.catphan_size
        self.mm_per_pixel = catphan.mm_per_pixel

    @property
    def __getitem__(self, item):
        return self.image.array[item]

    @property
    @lru_cache(maxsize=1)
    def phan_center(self):
        """Determine the location of the center of the phantom.

        The image is analyzed to see if:
        1) the CatPhan is even in the image (if there were any ROIs detected)
        2) an ROI is within the size criteria of the catphan
        3) the ROI area that is filled compared to the bounding box area is close to that of a circle

        Raises
        ------
        ValueError
            If any of the above conditions are not met.
        """
        # convert the slice to binary and label ROIs
        edges = filters.scharr(self.image.as_type(np.float))
        if np.max(edges) < 0.1:
            raise ValueError("Unable to locate Catphan")
        larr, regionprops, num_roi = get_regions(self, fill_holes=True, threshold='mean')
        # check that there is at least 1 ROI
        if num_roi < 1 or num_roi is None:
            raise ValueError("Unable to locate the CatPhan")
        catphan_region = sorted(regionprops, key=lambda x: np.abs(x.filled_area - self.catphan_size))[0]
        if (self.catphan_size * 1.2 < catphan_region.filled_area) or (catphan_region.filled_area < self.catphan_size / 1.2):
            raise ValueError("Unable to locate Catphan")
        center_pixel = catphan_region.centroid
        return Point(center_pixel[1], center_pixel[0])


class CatPhanModule(Slice):
    """Base class for a CTP module.
    """
    combine_method = 'mean'
    num_slices = 0
    roi_settings = {}
    background_roi_settings = {}
    roi_dist_mm = float
    roi_radius_mm = float
    rois = {}  # dicts of HUDiskROIs
    background_rois = {}  # dict of HUDiskROIs; possibly empty

    def __init__(self, catphan, tolerance, offset=0):
        """
        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        tolerance : float
        offset : int, float
        """
        self.model = ''
        self._offset = offset
        self.origin_slice = catphan.origin_slice
        self.tolerance = tolerance
        self.slice_thickness = catphan.dicom_stack.metadata.SliceThickness
        self.catphan_roll = catphan.catphan_roll
        self.mm_per_pixel = catphan.mm_per_pixel
        self.rois = {}  # dicts of HUDiskROIs
        self.background_rois = {}  # dict of HUDiskROIs; possibly empty
        Slice.__init__(self, catphan, combine_method=self.combine_method, num_slices=self.num_slices)
        self._convert_units_in_settings()
        self.preprocess(catphan)
        self._setup_rois()

    def _convert_units_in_settings(self):
        for roi, settings in self.roi_settings.items():
            self.roi_settings[roi]['distance_pixels'] = settings['distance'] / self.mm_per_pixel
            self.roi_settings[roi]['angle_corrected'] = settings['angle'] + self.catphan_roll
            self.roi_settings[roi]['radius_pixels'] = settings['radius'] / self.mm_per_pixel
        for roi, settings in self.background_roi_settings.items():
            self.background_roi_settings[roi]['distance_pixels'] = settings['distance'] / self.mm_per_pixel
            self.background_roi_settings[roi]['angle_corrected'] = settings['angle'] + self.catphan_roll
            self.background_roi_settings[roi]['radius_pixels'] = settings['radius'] / self.mm_per_pixel

    def preprocess(self, catphan):
        """A preprocessing step before analyzing the CTP module.

        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        """
        pass

    @property
    def slice_num(self):
        """The slice number of the spatial resolution module.

        Returns
        -------
        float
        """
        return int(self.origin_slice+round(self._offset/self.slice_thickness))

    def _setup_rois(self):
        for name, setting in self.background_roi_settings.items():
            self.background_rois[name] = HUDiskROI(self.image, setting['angle_corrected'], setting['radius_pixels'], setting['distance_pixels'],
                                        self.phan_center)
        background_median = np.mean([roi.pixel_value for roi in self.background_rois.values()])
        background_std = np.std([roi.pixel_value for roi in self.background_rois.values()])

        for name, setting in self.roi_settings.items():
            nominal_value = setting.get('value', 0)
            self.rois[name] = HUDiskROI(self.image, setting['angle_corrected'], setting['radius_pixels'], setting['distance_pixels'],
                                        self.phan_center, nominal_value, self.tolerance,
                                        background_median=background_median, background_std=background_std)

    def plot_rois(self, axis, threshold=None):
        """Plot the ROIs to the axis."""
        for roi in self.rois.values():
            if not threshold:
                roi.plot2axes(axis, edgecolor=roi.plot_color)
            else:
                roi.plot2axes(axis, edgecolor=roi.plot_color_cnr)
        for roi in self.background_rois.values():
            roi.plot2axes(axis, edgecolor='blue')

    @property
    def roi_vals_as_str(self):
        return ', '.join(f'{name}: {roi.pixel_value}' for name, roi in self.rois.items())


class CTP404CP504(CatPhanModule):
    """Class for analysis of the HU linearity, geometry, and slice thickness regions of the CTP404.
    """
    attr_name = 'ctp404'
    common_name = 'HU Linearity'
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        'Air': {'value': AIR, 'angle': -90, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'PMP': {'value': PMP, 'angle': -120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'LDPE': {'value': LDPE, 'angle': 180, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Poly': {'value': POLY, 'angle': 120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Acrylic': {'value': ACRYLIC, 'angle': 60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Delrin': {'value': DELRIN, 'angle': 0, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Teflon': {'value': TEFLON, 'angle': -60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }
    background_roi_settings = {
        '1': {'angle': -30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '2': {'angle': -150, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '3': {'angle': -210, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '4': {'angle': 30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }
    # thickness
    thickness_roi_height = 40
    thickness_roi_width = 10
    thickness_roi_distance_mm = 38
    thickness_roi_settings = {
        'Left': {'angle': 180, 'width': thickness_roi_width, 'height': thickness_roi_height, 'distance': thickness_roi_distance_mm},
        'Bottom': {'angle': 90, 'width': thickness_roi_height, 'height': thickness_roi_width, 'distance': thickness_roi_distance_mm},
        'Right': {'angle': 0, 'width': thickness_roi_width, 'height': thickness_roi_height, 'distance': thickness_roi_distance_mm},
        'Top': {'angle': -90, 'width': thickness_roi_height, 'height': thickness_roi_width, 'distance': thickness_roi_distance_mm},
    }
    # geometry
    geometry_roi_size_mm = 35
    geometry_roi_settings = {
        'Top-Horizontal': (0, 1),
        'Bottom-Horizontal': (2, 3),
        'Left-Vertical': (0, 2),
        'Right-Vertical': (1, 3),
    }

    def __init__(self, catphan, offset, hu_tolerance, thickness_tolerance, scaling_tolerance):
        """
        Parameters
        ----------
        catphan : `~pylinac.cbct.CatPhanBase` instance.
        offset : float
        hu_tolerance : float
        thickness_tolerance : float
        scaling_tolerance : float
        """
        self.mm_per_pixel = catphan.mm_per_pixel
        self.hu_tolerance = hu_tolerance
        self.thickness_tolerance = thickness_tolerance
        self.scaling_tolerance = scaling_tolerance
        self.thickness_rois = {}
        self.lines = {}
        super().__init__(catphan, tolerance=hu_tolerance, offset=offset)

    def _convert_units_in_settings(self):
        super()._convert_units_in_settings()
        for roi, settings in self.thickness_roi_settings.items():
            self.thickness_roi_settings[roi]['width_pixels'] = settings['width'] / self.mm_per_pixel
            self.thickness_roi_settings[roi]['height_pixels'] = settings['height'] / self.mm_per_pixel
            self.thickness_roi_settings[roi]['angle_corrected'] = settings['angle'] + self.catphan_roll
            self.thickness_roi_settings[roi]['distance_pixels'] = settings['distance'] / self.mm_per_pixel

    def preprocess(self, catphan):
        # for the thickness analysis image, combine thin slices or just use one slice if slices are thick
        if float(catphan.dicom_stack.metadata.SliceThickness) < 3.5:
            self.pad = 1
        else:
            self.pad = 0
        self.thickness_image = Slice(catphan, combine_method='mean', num_slices=self.num_slices+self.pad, slice_num=self.slice_num).image

    def _setup_rois(self):
        super()._setup_rois()
        self._setup_thickness_rois()
        self._setup_geometry_rois()

    def _setup_thickness_rois(self):
        for name, setting in self.thickness_roi_settings.items():
            self.thickness_rois[name] = ThicknessROI(self.thickness_image, setting['width_pixels'],
                                                     setting['height_pixels'], setting['angle_corrected'],
                                                     setting['distance_pixels'], self.phan_center)

    def _setup_geometry_rois(self):
        boxsize = self.geometry_roi_size_mm / self.mm_per_pixel
        xbounds = (int(self.phan_center.x-boxsize), int(self.phan_center.x+boxsize))
        ybounds = (int(self.phan_center.y-boxsize), int(self.phan_center.y+boxsize))
        geo_img = self.image[ybounds[0]:ybounds[1], xbounds[0]:xbounds[1]]
        larr, regionprops, num_roi = get_regions(geo_img, fill_holes=True, clear_borders=False)
        # check that there is at least 1 ROI
        if num_roi < 4:
            raise ValueError("Unable to locate the Geometric nodes")
        elif num_roi > 4:
            regionprops = sorted(regionprops, key=lambda x: x.filled_area, reverse=True)[:4]
        sorted_regions = sorted(regionprops, key=lambda x: (2*x.centroid[0]+x.centroid[1]))
        centers = [Point(r.weighted_centroid[1]+xbounds[0], r.weighted_centroid[0]+ybounds[0]) for r in sorted_regions]
        #  setup the geometric lines
        for name, order in self.geometry_roi_settings.items():
            self.lines[name] = GeometricLine(centers[order[0]], centers[order[1]], self.mm_per_pixel, self.scaling_tolerance)

    @property
    def lcv(self):
        """The low-contrast visibility"""
        return 2 * abs(self.rois['LDPE'].pixel_value - self.rois['Poly'].pixel_value) / (self.rois['LDPE'].std + self.rois['Poly'].std)

    def plot_linearity(self, axis=None, plot_delta=True):
        """Plot the HU linearity values to an axis.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the values on. If None, will create a new figure.
        plot_delta : bool
            Whether to plot the actual measured HU values (False), or the difference from nominal (True).
        """
        nominal_x_values = [roi.nominal_val for roi in self.rois.values()]
        if axis is None:
            fig, axis = plt.subplots()
        if plot_delta:
            values = [roi.value_diff for roi in self.rois.values()]
            nominal_measurements = [0]*len(values)
            ylabel = 'HU Delta'
        else:
            values = [roi.pixel_value for roi in self.rois.values()]
            nominal_measurements = nominal_x_values
            ylabel = 'Measured Values'
        points = axis.plot(nominal_x_values, values, 'g+', markersize=15, mew=2)
        axis.plot(nominal_x_values, nominal_measurements)
        axis.plot(nominal_x_values, np.array(nominal_measurements) + self.hu_tolerance, 'r--')
        axis.plot(nominal_x_values, np.array(nominal_measurements) - self.hu_tolerance, 'r--')
        axis.margins(0.05)
        axis.grid(True)
        axis.set_xlabel("Nominal Values")
        axis.set_ylabel(ylabel)
        axis.set_title("HU linearity")
        return points

    @property
    def passed_hu(self):
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())

    def plot_rois(self, axis):
        """Plot the ROIs onto the image, as well as the background ROIs"""
        # plot HU linearity ROIs
        super().plot_rois(axis)
        # plot thickness ROIs
        for roi in self.thickness_rois.values():
            roi.plot2axes(axis, edgecolor='blue')
        # plot geometry lines
        for line in self.lines.values():
            line.plot2axes(axis, color=line.pass_fail_color)

    @property
    def passed_thickness(self):
        """Whether the slice thickness was within tolerance from nominal."""
        return self.slice_thickness-self.thickness_tolerance<self.meas_slice_thickness<self.slice_thickness+self.thickness_tolerance

    @property
    def meas_slice_thickness(self):
        """The average slice thickness for the 4 wire measurements in mm."""
        return np.mean(sorted(roi.wire_fwhm*self.mm_per_pixel*0.42 for roi in self.thickness_rois.values()))/(1+2*self.pad)

    @property
    def avg_line_length(self):
        return np.mean([line.length_mm for line in self.lines.values()])

    @property
    def passed_geometry(self):
        """Returns whether all the line lengths were within tolerance."""
        return all(line.passed for line in self.lines.values())


class CTP404CP503(CTP404CP504):
    """Alias for namespace consistency"""
    pass


class CTP404CP600(CTP404CP504):
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        'Air': {'value': AIR, 'angle': 90, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'PMP': {'value': PMP, 'angle': 60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'LDPE': {'value': LDPE, 'angle': 0, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Poly': {'value': POLY, 'angle': -60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Acrylic': {'value': ACRYLIC, 'angle': -120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Delrin': {'value': DELRIN, 'angle': -180, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Teflon': {'value': TEFLON, 'angle': 120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }


class CTP404CP604(CTP404CP504):
    roi_dist_mm = 58.7
    roi_radius_mm = 5
    roi_settings = {
        'Air': {'value': AIR, 'angle': -90, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'PMP': {'value': PMP, 'angle': -120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '50% Bone': {'value': BONE_50, 'angle': -150, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'LDPE': {'value': LDPE, 'angle': 180, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Poly': {'value': POLY, 'angle': 120, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Acrylic': {'value': ACRYLIC, 'angle': 60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '20% Bone': {'value': BONE_20, 'angle': 30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Delrin': {'value': DELRIN, 'angle': 0, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Teflon': {'value': TEFLON, 'angle': -60, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }
    background_roi_settings = {
        '1': {'angle': -30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '2': {'angle': -210, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }


class CTP486(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """
    attr_name = 'ctp486'
    common_name = 'HU Uniformity'
    roi_dist_mm = 53
    roi_radius_mm = 10
    nominal_value = 0
    roi_settings = {
        'Top': {'value': nominal_value, 'angle': 90, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Right': {'value': nominal_value, 'angle': 0, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Bottom': {'value': nominal_value, 'angle': -90, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Left': {'value': nominal_value, 'angle': 180, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        'Center': {'value': nominal_value, 'angle': 0, 'distance': 0, 'radius': roi_radius_mm},
    }

    def plot_profiles(self, axis=None):
        """Plot the horizontal and vertical profiles of the Uniformity slice.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot on; if None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        horiz_data = self.image[int(self.phan_center.y), :]
        vert_data = self.image[:, int(self.phan_center.x)]
        axis.plot(horiz_data, 'g', label='Horizontal')
        axis.plot(vert_data, 'b', label='Vertical')
        axis.autoscale(tight=True)
        axis.axhline(self.tolerance, color='r', linewidth=3)
        axis.axhline(-self.tolerance, color='r', linewidth=3)
        axis.grid(True)
        axis.set_ylabel("HU")
        axis.legend(loc=8, fontsize='small', title="")
        axis.set_title("Uniformity Profiles")

    @property
    def overall_passed(self):
        """Boolean specifying whether all the ROIs passed within tolerance."""
        return all(roi.passed for roi in self.rois.values())

    @property
    def uniformity_index(self):
        """The Uniformity Index"""
        center = self.rois['Center']
        uis = [100*((roi.pixel_value-center.pixel_value)/(center.pixel_value+1000)) for roi in self.rois.values()]
        abs_uis = np.abs(uis)
        return uis[np.argmax(abs_uis)]

    @property
    def integral_non_uniformity(self):
        """The Integral Non-Uniformity"""
        maxhu = max(roi.pixel_value for roi in self.rois.values())
        minhu = min(roi.pixel_value for roi in self.rois.values())
        return (maxhu - minhu)/(maxhu + minhu + 2000)


class CTP528CP504(CatPhanModule):
    """Class for analysis of the Spatial Resolution slice of the CBCT dicom data set.

    A collapsed circle profile is taken of the line-pair region. This profile is search for
    peaks and valleys. The MTF is calculated from those peaks & valleys.

    Attributes
    ----------
    radius2linepairs_mm : float
        The radius in mm to the line pairs.
    """
    attr_name = 'ctp528'
    common_name = 'Spatial Resolution'
    radius2linepairs_mm = 47
    combine_method = 'max'
    num_slices = 3
    boundaries = (0, 0.107, 0.173, 0.236, 0.286, 0.335, 0.387, 0.434, 0.479)
    start_angle = np.pi
    ccw = True
    roi_settings = {
        'region 1': {'start': boundaries[0], 'end': boundaries[1], 'num peaks': 2, 'num valleys': 1,
                     'peak spacing': 0.021, 'gap size (cm)': 0.5, 'lp/mm': 0.1},
        'region 2': {'start': boundaries[1], 'end': boundaries[2], 'num peaks': 3, 'num valleys': 2,
                     'peak spacing': 0.01, 'gap size (cm)': 0.25, 'lp/mm': 0.2},
        'region 3': {'start': boundaries[2], 'end': boundaries[3], 'num peaks': 4, 'num valleys': 3,
                     'peak spacing': 0.006, 'gap size (cm)': 0.167, 'lp/mm': 0.3},
        'region 4': {'start': boundaries[3], 'end': boundaries[4], 'num peaks': 4, 'num valleys': 3,
                     'peak spacing': 0.00557, 'gap size (cm)': 0.125, 'lp/mm': 0.4},
        'region 5': {'start': boundaries[4], 'end': boundaries[5], 'num peaks': 4, 'num valleys': 3,
                     'peak spacing': 0.004777, 'gap size (cm)': 0.1, 'lp/mm': 0.5},
        'region 6': {'start': boundaries[5], 'end': boundaries[6], 'num peaks': 5, 'num valleys': 4,
                     'peak spacing': 0.00398, 'gap size (cm)': 0.083, 'lp/mm': 0.6},
        'region 7': {'start': boundaries[6], 'end': boundaries[7], 'num peaks': 5, 'num valleys': 4,
                     'peak spacing': 0.00358, 'gap size (cm)': 0.071, 'lp/mm': 0.7},
        'region 8': {'start': boundaries[7], 'end': boundaries[8], 'num peaks': 5, 'num valleys': 4,
                     'peak spacing': 0.0027866, 'gap size (cm)': 0.063, 'lp/mm': 0.8},
    }

    def _setup_rois(self):
        pass

    def _convert_units_in_settings(self):
        pass

    @property
    @lru_cache(maxsize=1)
    def mtf(self):
        """The Relative MTF of the line pairs, normalized to the first region.

        Returns
        -------
        dict
        """
        maxs = list()
        mins = list()
        for key, value in self.roi_settings.items():
            max_values = self.circle_profile.find_peaks(min_distance=value['peak spacing'], max_number=value['num peaks'],
                                                        search_region=(value['start'], value['end']), kind='value')
            # check that the right number of peaks were found before continuing, otherwise stop searching for regions
            if len(max_values) != value['num peaks']:
                break
            maxs.append(max_values.mean())
            max_indices = self.circle_profile.find_peaks(min_distance=value['peak spacing'], max_number=value['num peaks'],
                                                         search_region=(value['start'], value['end']), kind='index')
            lower_mean = self.circle_profile.find_valleys(min_distance=value['peak spacing'], max_number=value['num valleys'],
                                                          search_region=(min(max_indices), max(max_indices)), kind='value').mean()
            mins.append(lower_mean)
        if not maxs:
            raise ValueError("Did not find any spatial resolution pairs to analyze. File an issue on github (https://github.com/jrkerns/pylinac/issues) if this is a valid dataset.")

        spacings = [roi['lp/mm'] for roi in self.roi_settings.values()]
        mtf = MTF(lp_spacings=spacings, lp_maximums=maxs, lp_minimums=mins)
        return mtf

    @property
    def radius2linepairs(self):
        """Radius from the phantom center to the line-pair region, corrected for pixel spacing."""
        return self.radius2linepairs_mm / self.mm_per_pixel

    def plot_rois(self, axis):
        """Plot the circles where the profile was taken within."""
        self.circle_profile.plot2axes(axis, edgecolor='blue', plot_peaks=False)

    @property
    @lru_cache(maxsize=1)
    def circle_profile(self):
        """Calculate the median profile of the Line Pair region.

        Returns
        -------
        :class:`pylinac.core.profile.CollapsedCircleProfile` : A 1D profile of the Line Pair region.
        """
        circle_profile = CollapsedCircleProfile(self.phan_center, self.radius2linepairs, image_array=self.image,
                                                start_angle=self.start_angle + np.deg2rad(self.catphan_roll),
                                                width_ratio=0.04, sampling_ratio=2, ccw=self.ccw)
        circle_profile.filter(0.001, kind='gaussian')
        circle_profile.ground()
        return circle_profile

    def plot_mtf(self, axis=None):
        """Plot the Relative MTF.

        Parameters
        ----------
        axis : None, matplotlib.Axes
            The axis to plot the MTF on. If None, will create a new figure.
        """
        if axis is None:
            fig, axis = plt.subplots()
        points = axis.plot(list(self.mtf.norm_mtfs.keys()), list(self.mtf.norm_mtfs.values()), marker='o')
        axis.margins(0.05)
        axis.grid(True)
        axis.set_xlabel('Line pairs / mm')
        axis.set_ylabel("Relative MTF")
        axis.set_title('RMTF')
        return points


class CTP528CP604(CTP528CP504):
    """Alias for namespace consistency."""
    pass


class CTP528CP600(CTP528CP504):
    start_angle = np.pi - 0.1
    ccw = False
    boundaries = (0, 0.116, 0.182, 0.244, 0.294, 0.344, 0.396, 0.443, 0.488)


class CTP528CP503(CTP528CP504):
    start_angle = 0
    ccw = False
    boundaries = (0, 0.111, 0.176, 0.240, 0.289, 0.339, 0.390, 0.436, 0.481)


class GeometricLine(Line):
    """Represents a line connecting two nodes/ROIs on the Geometry Slice.

    Attributes
    ----------
    nominal_length_mm : int, float
        The nominal distance between the geometric nodes, in mm.
    """
    nominal_length_mm = 50

    def __init__(self, geo_roi1, geo_roi2, mm_per_pixel, tolerance):
        """
        Parameters
        ----------
        geo_roi1 : GEO_ROI
            One of two ROIs representing one end of the line.
        geo_roi2 : GEO_ROI
            The other ROI which is the other end of the line.
        mm_per_pixel : float
            The mm/pixel value.
        tolerance : int, float
            The tolerance of the geometric line, in mm.
        """
        super().__init__(geo_roi1, geo_roi2)
        self.mm_per_pixel = mm_per_pixel
        self.tolerance = tolerance

    @property
    def passed(self):
        """Whether the line passed tolerance."""
        return self.nominal_length_mm - self.tolerance < self.length_mm < self.nominal_length_mm + self.tolerance

    @property
    def pass_fail_color(self):
        """Plot color for the line, based on pass/fail status."""
        return 'blue' if self.passed else 'red'

    @property
    def length_mm(self):
        """Return the length of the line in mm."""
        return self.length*self.mm_per_pixel


class CTP515(CatPhanModule):
    """Class for analysis of the low contrast slice of the CTP module. Low contrast is measured by obtaining
    the average pixel value of the contrast ROIs and comparing that value to the average background value. To obtain
    a more "human" detection level, the contrast (which is largely the same across different-sized ROIs) is multiplied
    by the diameter. This value is compared to the contrast threshold to decide if it can be "seen".
    """
    attr_name = 'ctp515'
    common_name = 'Low Contrast'
    num_slices = 1
    roi_dist_mm = 50
    roi_radius_mm = [6, 3.5, 3, 2.5, 2, 1.5]
    roi_angles = [-87.4, -69.1, -52.7, -38.5, -25.1, -12.9]
    roi_settings = {
        '15': {'angle': roi_angles[0], 'distance': roi_dist_mm, 'radius': roi_radius_mm[0]},
        '9': {'angle': roi_angles[1], 'distance': roi_dist_mm, 'radius': roi_radius_mm[1]},
        '8': {'angle': roi_angles[2], 'distance': roi_dist_mm, 'radius': roi_radius_mm[2]},
        '7': {'angle': roi_angles[3], 'distance': roi_dist_mm, 'radius': roi_radius_mm[3]},
        '6': {'angle': roi_angles[4], 'distance': roi_dist_mm, 'radius': roi_radius_mm[4]},
        '5': {'angle': roi_angles[5], 'distance': roi_dist_mm, 'radius': roi_radius_mm[5]},
    }
    background_roi_dist_ratio = 0.75
    background_roi_radius_mm = 4

    def __init__(self, catphan, tolerance, cnr_threshold, offset):
        self.cnr_threshold = cnr_threshold
        super().__init__(catphan, tolerance=tolerance, offset=offset)

    def _setup_rois(self):
        # create both background rois dynamically, then create the actual sample ROI as normal
        for name, setting in self.roi_settings.items():
            self.background_rois[name+'-outer'] = LowContrastDiskROI(self.image, setting['angle_corrected'],
                                                                     self.background_roi_radius_mm / self.mm_per_pixel,
                                                                     setting['distance_pixels'] * (2-self.background_roi_dist_ratio),
                                                                     self.phan_center)
            self.background_rois[name+'-inner'] = LowContrastDiskROI(self.image, setting['angle_corrected'],
                                                                     self.background_roi_radius_mm / self.mm_per_pixel,
                                                                     setting['distance_pixels'] * self.background_roi_dist_ratio,
                                                                     self.phan_center)
            background_val = np.mean([self.background_rois[name+'-outer'].pixel_value, self.background_rois[name+'-inner'].pixel_value])

            self.rois[name] = LowContrastDiskROI(self.image, setting['angle_corrected'], setting['radius_pixels'], setting['distance_pixels'],
                                                 self.phan_center, background=background_val, cnr_threshold=self.cnr_threshold)

    def plot_rois(self, axis, threshold='constant'):
        super().plot_rois(axis, threshold)

    @property
    def rois_visible(self):
        """The number of ROIs "visible"."""
        return sum(roi.passed_cnr_constant for roi in self.rois.values())


class CatPhanBase:
    """A class for loading and analyzing CT DICOM files of a CatPhan 504 & CatPhan 503. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528), Image Scaling & HU Linearity (CTP404).
    """
    _demo_url = ''
    _model = ''
    air_bubble_radius_mm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__(self, folderpath, check_uid=True):
        """
        Parameters
        ----------
        folderpath : str
            String that points to the CBCT image folder location.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.

        Raises
        ------
        NotADirectoryError : If folder str passed is not a valid directory.
        FileNotFoundError : If no CT images are found in the folder
        """
        self.origin_slice = 0
        self.catphan_roll = 0
        if not osp.isdir(folderpath):
            raise NotADirectoryError("Path given was not a Directory/Folder")
        self.dicom_stack = image.DicomImageStack(folderpath, check_uid=check_uid)
        self.localize()

    @classmethod
    def from_demo_images(cls):
        """Construct a CBCT object from the demo images."""
        demo_file = retrieve_demo_file(url=cls._demo_url)
        return cls.from_zip(demo_file)

    @classmethod
    def from_url(cls, url, check_uid=True):
        """Instantiate a CBCT object from a URL pointing to a .zip object.

        Parameters
        ----------
        url : str
            URL pointing to a zip archive of CBCT images.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.
        """
        filename = get_url(url)
        return cls.from_zip(filename, check_uid=check_uid)

    @classmethod
    def from_zip(cls, zip_file, check_uid=True):
        """Construct a CBCT object and pass the zip file.

        Parameters
        ----------
        zip_file : str, ZipFile
            Path to the zip file or a ZipFile object.
        check_uid : bool
            Whether to enforce raising an error if more than one UID is found in the dataset.

        Raises
        ------
        FileExistsError : If zip_file passed was not a legitimate zip file.
        FileNotFoundError : If no CT images are found in the folder
        """
        with TemporaryZipDirectory(zip_file) as temp_zip:
            obj = cls(temp_zip, check_uid=check_uid)
        obj.was_from_zip = True
        return obj

    def plot_analyzed_image(self, show=True):
        """Plot the images used in the calculate and summary data.

        Parameters
        ----------
        show : bool
            Whether to plot the image or not.
        """
        def plot(ctp_module, axis):
            axis.imshow(ctp_module.image.array, cmap=get_dicom_cmap())
            ctp_module.plot_rois(axis)
            axis.autoscale(tight=True)
            axis.set_title(ctp_module.common_name)
            axis.axis('off')

        # set up grid and axes
        grid_size = (2, 4)
        hu_ax = plt.subplot2grid(grid_size, (0, 1))
        plot(self.ctp404, hu_ax)
        hu_lin_ax = plt.subplot2grid(grid_size, (0, 2))
        self.ctp404.plot_linearity(hu_lin_ax)
        if self._has_module(CTP486):
            unif_ax = plt.subplot2grid(grid_size, (0, 0))
            plot(self.ctp486, unif_ax)
            unif_prof_ax = plt.subplot2grid(grid_size, (1, 2), colspan=2)
            self.ctp486.plot_profiles(unif_prof_ax)
        if self._has_module(CTP528CP504):
            sr_ax = plt.subplot2grid(grid_size, (1, 0))
            plot(self.ctp528, sr_ax)
            mtf_ax = plt.subplot2grid(grid_size, (0, 3))
            self.ctp528.plot_mtf(mtf_ax)
        if self._has_module(CTP515):
            locon_ax = plt.subplot2grid(grid_size, (1, 1))
            plot(self.ctp515, locon_ax)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()

    def save_analyzed_image(self, filename, **kwargs):
        """Save the analyzed summary plot.

        Parameters
        ----------
        filename : str, file object
            The name of the file to save the image to.
        kwargs :
            Any valid matplotlib kwargs.
        """
        self.plot_analyzed_image(show=False)
        plt.savefig(filename, **kwargs)

    def plot_analyzed_subimage(self, subimage='hu', delta=True, show=True):
        """Plot a specific component of the CBCT analysis.

        Parameters
        ----------
        subimage : {'hu', 'un', 'sp', 'lc', 'mtf', 'lin', 'prof'}
            The subcomponent to plot. Values must contain one of the following letter combinations.
            E.g. ``linearity``, ``linear``, and ``lin`` will all draw the HU linearity values.

            * ``hu`` draws the HU linearity image.
            * ``un`` draws the HU uniformity image.
            * ``sp`` draws the Spatial Resolution image.
            * ``mtf`` draws the RMTF plot.
            * ``lin`` draws the HU linearity values. Used with ``delta``.
            * ``prof`` draws the HU uniformity profiles.
        delta : bool
            Only for use with ``lin``. Whether to plot the HU delta or actual values.
        show : bool
            Whether to actually show the plot.
        """
        subimage = subimage.lower()
        plt.clf()
        plt.axis('off')

        if 'hu' in subimage:  # HU, GEO & thickness objects
            plt.imshow(self.ctp404.image.array, cmap=get_dicom_cmap())
            self.ctp404.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'un' in subimage:  # uniformity
            plt.imshow(self.ctp486.image.array, cmap=get_dicom_cmap())
            self.ctp486.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'sp' in subimage:  # SR objects
            plt.imshow(self.ctp528.image.array, cmap=get_dicom_cmap())
            self.ctp528.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'mtf' in subimage:
            plt.axis('on')
            self.ctp528.plot_mtf(plt.gca())
        elif 'lc' in subimage:
            plt.imshow(self.ctp515.image.array, cmap=get_dicom_cmap())
            self.ctp515.plot_rois(plt.gca())
            plt.autoscale(tight=True)
        elif 'lin' in subimage:
            plt.axis('on')
            self.ctp404.plot_linearity(plt.gca(), delta)
        elif 'prof' in subimage:
            plt.axis('on')
            self.ctp486.plot_profiles(plt.gca())
        else:
            raise ValueError(f"Subimage parameter {subimage} not understood")

        if show:
            plt.show()

    def save_analyzed_subimage(self, filename, subimage='hu', **kwargs):
        """Save a component image to file.

        Parameters
        ----------
        filename : str, file object
            The file to write the image to.
        subimage : str
            See :meth:`~pylinac.cbct.CBCT.plot_analyzed_subimage` for parameter info.
        """
        self.plot_analyzed_subimage(subimage, show=False)
        plt.savefig(filename, **kwargs)
        if isinstance(filename, str):
            print(f"CatPhan subimage figure saved to {osp.abspath(filename)}")

    def _results(self):
        """Helper function to spit out values that will be tested."""
        print(self.results())
        print(f"Phantom roll: {self.find_phantom_roll()}")
        print(f"Origin slice: {self.find_origin_slice()}")
        mtfs = {}
        for mtf in (95, 90, 80, 50, 30):
            mtfval = self.ctp528.mtf.relative_resolution(mtf)
            mtfs[mtf] = mtfval
        print(f'MTFs: {mtfs}')

    def localize(self):
        """Find the slice number of the catphan's HU linearity module and roll angle"""
        self.origin_slice = self.find_origin_slice()
        self.catphan_roll = self.find_phantom_roll()

    @property
    def mm_per_pixel(self):
        """The millimeters per pixel of the DICOM images."""
        return self.dicom_stack.metadata.PixelSpacing[0]

    @lru_cache(maxsize=1)
    def find_origin_slice(self):
        """Using a brute force search of the images, find the median HU linearity slice.

        This method walks through all the images and takes a collapsed circle profile where the HU
        linearity ROIs are. If the profile contains both low (<800) and high (>800) HU values and most values are the same
        (i.e. it's not an artifact), then
        it can be assumed it is an HU linearity slice. The median of all applicable slices is the
        center of the HU slice.

        Returns
        -------
        int
            The middle slice of the HU linearity module.
        """
        hu_slices = []
        for image_number in range(0, self.num_images, 2):
            slice = Slice(self, image_number, combine=False)
            #print(image_number)
            # slice.image.plot()
            try:
                center = slice.phan_center
            except ValueError:  # a slice without the phantom in view
                pass
            else:
                circle_prof = CollapsedCircleProfile(center, radius=self.localization_radius/self.mm_per_pixel, image_array=slice.image, width_ratio=0.05, num_profiles=5)
                prof = circle_prof.values
                # determine if the profile contains both low and high values and that most values are the same
                low_end, high_end = np.percentile(prof, [2, 98])
                median = np.median(prof)
                middle_variation = np.percentile(prof, 80) - np.percentile(prof, 20)
                variation_limit = max(100, self.dicom_stack.metadata.SliceThickness*-100+300)
                if (low_end < median - 400) and (high_end > median + 400) and (middle_variation < variation_limit):
                    hu_slices.append(image_number)

        if not hu_slices:
            raise ValueError("No slices were found that resembled the HU linearity module")
        hu_slices = np.array(hu_slices)
        c = int(round(np.median(hu_slices)))
        ln = len(hu_slices)
        # drop slices that are way far from median
        hu_slices = hu_slices[((c + ln/2) >= hu_slices) & (hu_slices >= (c - ln/2))]
        center_hu_slice = int(round(np.median(hu_slices)))
        if self._is_within_image_extent(center_hu_slice):
            #print(center_hu_slice)
            return center_hu_slice

    @lru_cache(maxsize=1)
    def find_phantom_roll(self):
        """Determine the "roll" of the phantom.

        This algorithm uses the two air bubbles in the HU slice and the resulting angle between them.

        Returns
        -------
        float : the angle of the phantom in **degrees**.
        """
        def is_right_area(region):
            thresh = np.pi * ((self.air_bubble_radius_mm / self.mm_per_pixel) ** 2)
            return thresh * 2 > region.filled_area > thresh / 2

        def is_right_eccentricity(region):
            return region.eccentricity < 0.5

        # get edges and make ROIs from it
        slice = Slice(self, self.origin_slice)
        larr, regions, _ = get_regions(slice)
        # find appropriate ROIs and grab the two most centrally positioned ones
        hu_bubbles = [r for r in regions if (is_right_area(r) and is_right_eccentricity(r))]
        central_bubbles = sorted(hu_bubbles, key=lambda x: abs(x.centroid[1] - slice.phan_center.x))[:2]
        sorted_bubbles = sorted(central_bubbles, key=lambda x: x.centroid[0])  # top, bottom
        y_dist = sorted_bubbles[1].centroid[0] - sorted_bubbles[0].centroid[0]
        x_dist = sorted_bubbles[1].centroid[1] - sorted_bubbles[0].centroid[1]
        phan_roll = np.arctan2(y_dist, x_dist)
        anglroll = np.rad2deg(phan_roll) - 90
        return anglroll

    @property
    def num_images(self):
        """Return the number of images loaded.

        Returns
        -------
        float
        """
        return len(self.dicom_stack)

    def _is_within_image_extent(self, image_num):
        """Determine if the image number is beyond the edges of the images (negative or past last image)."""
        if self.num_images - 1 > image_num > 1:
            return True
        else:
            raise ValueError("The determined image number is beyond the image extent. Either the entire dataset "
                             "wasn't loaded or the entire phantom wasn't scanned.")

    @property
    def catphan_size(self):
        """The expected size of the phantom in pixels, based on a 20cm wide phantom.

        Returns
        -------
        float
        """
        phan_area = np.pi*(self.catphan_radius_mm**2)
        return phan_area/(self.mm_per_pixel**2)

    def publish_pdf(self, filename: str, notes: str=None, open_file: bool=False, metadata: Optional[dict]=None):
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        """
        analysis_title = f'CatPhan {self._model} Analysis'
        module_texts = [
            [' - CTP404 Results - ',
             f'HU Linearity tolerance: {self.ctp404.hu_tolerance}',
             f'HU Linearity ROIs: {self.ctp404.roi_vals_as_str}',
             f'Geometric node spacing (mm): {self.ctp404.avg_line_length:2.2f}',
             f'Slice thickness (mm): {self.ctp404.meas_slice_thickness:2.2f}',
             f'Low contrast visibility: {self.ctp404.lcv:2.2f}',
            ],
        ]
        module_images = [('hu', 'lin')]
        if self._has_module(CTP528CP504):
            add = [' - CTP528 Results - ',
             f'MTF 80% (lp/mm): {self.ctp528.mtf.relative_resolution(80):2.2f}',
             f'MTF 50% (lp/mm): {self.ctp528.mtf.relative_resolution(50):2.2f}',
             f'MTF 30% (lp/mm): {self.ctp528.mtf.relative_resolution(30):2.2f}',
            ]
            module_texts.append(add)
            module_images.append(('sp', 'mtf'))
        if self._has_module(CTP486):
            add = [' - CTP486 Results - ',
             f'Uniformity tolerance: {self.ctp486.tolerance}',
             f'Uniformity ROIs: {self.ctp486.roi_vals_as_str}',
             f'Uniformity Index: {self.ctp486.uniformity_index:2.2f}',
             f'Integral non-uniformity: {self.ctp486.integral_non_uniformity:2.4f}',
            ]
            module_texts.append(add)
            module_images.append(('un', 'prof'))
        if self._has_module(CTP515):
            add = [' - CTP515 Results - ',
             f'CNR threshold: {self.ctp515.cnr_threshold}',
             f'Low contrast ROIs "seen": {self.ctp515.rois_visible}'
            ]
            module_texts.append(add)
            module_images.append(('lc', None))

        self._publish_pdf(filename, metadata, notes, analysis_title,
                          module_texts, module_images)
        if open_file:
            webbrowser.open(filename)

    def _publish_pdf(self, filename, metadata, notes, analysis_title, texts, imgs):
        try:
            date = datetime.strptime(self.dicom_stack[0].metadata.InstanceCreationDate, "%Y%m%d").strftime("%A, %B %d, %Y")
        except:
            date = "Unknown"
        canvas = pdf.PylinacCanvas(filename, page_title=analysis_title, metadata=metadata)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        for page, ((img1, img2), text) in enumerate(zip(imgs, texts)):
            for img, offset in zip((img1, img2), (12, 2)):
                if img is not None:
                    data = io.BytesIO()
                    self.save_analyzed_subimage(data, img)
                    canvas.add_image(data, location=(4, offset), dimensions=(15, 10))
            canvas.add_text(text=text, location=(1.5, 23))
            canvas.add_new_page()
        canvas.finish()

    def _zip_images(self):
        """Compress the raw images into a ZIP archive and remove the uncompressed images."""
        zip_name = f'{osp.dirname(self.dicom_stack[0].path)}\CBCT - {self.dicom_stack[0].date_created(format="%A, %I-%M-%S, %B %d, %Y")}.zip'
        with zipfile.ZipFile(zip_name, 'w', compression=zipfile.ZIP_DEFLATED) as zfile:
            for image in self.dicom_stack:
                zfile.write(image.path, arcname=osp.basename(image.path))
        for image in self.dicom_stack:
            try:
                os.remove(image.path)
            except:
                pass

    def analyze(self, hu_tolerance=40, scaling_tolerance=1, thickness_tolerance=0.2,
                low_contrast_tolerance=1, cnr_threshold=15, zip_after=False):
        """Single-method full analysis of CBCT DICOM files.

        Parameters
        ----------
        hu_tolerance : int
            The HU tolerance value for both HU uniformity and linearity.
        scaling_tolerance : float, int
            The scaling tolerance in mm of the geometric nodes on the HU linearity slice (CTP404 module).
        thickness_tolerance : float, int
            The tolerance of the thickness calculation in mm, based on the wire ramps in the CTP404 module.

            .. warning:: Thickness accuracy degrades with image noise; i.e. low mAs images are less accurate.

        low_contrast_tolerance : int
            The number of low-contrast bubbles needed to be "seen" to pass.
        cnr_threshold : float, int
            The threshold for "detecting" low-contrast image. See RTD for calculation info.
        zip_after : bool
            If the CT images were not compressed before analysis and this is set to true, pylinac will compress
            the analyzed images into a ZIP archive.
        """
        ctp404, offset = self._get_module(CTP404CP504, raise_empty=True)
        self.ctp404 = ctp404(self, offset=offset, hu_tolerance=hu_tolerance, thickness_tolerance=thickness_tolerance,
                             scaling_tolerance=scaling_tolerance)
        if self._has_module(CTP486):
            ctp486, offset = self._get_module(CTP486)
            self.ctp486 = ctp486(self, offset=offset, tolerance=hu_tolerance)
        if self._has_module(CTP528CP504):
            ctp528, offset = self._get_module(CTP528CP504)
            self.ctp528 = ctp528(self, offset=offset, tolerance=None)
        if self._has_module(CTP515):
            ctp515, offset = self._get_module(CTP515)
            self.ctp515 = ctp515(self, tolerance=low_contrast_tolerance, cnr_threshold=cnr_threshold,
                                 offset=offset)
        if zip_after and not self.was_from_zip:
            self._zip_images()

    def _has_module(self, module_of_interest):
        return any(issubclass(module, module_of_interest) for module in self.modules.keys())

    def _get_module(self, module_of_interest, raise_empty=False):
        """Grab the module that is, or is a subclass of, the module of interest. This allows users to subclass a CTP module and pass that in."""
        for module, values in self.modules.items():
            if issubclass(module, module_of_interest):
                return module, values['offset']
        if raise_empty:
            raise ValueError(f"Tried to find the {module_of_interest} or a subclass of it. Did you override `modules` and not pass this module in?")

    def results(self):
        """Return the results of the analysis as a string. Use with print()."""
        string = (f'\n - CatPhan {self._model} QA Test - \n'
                  f'HU Linearity ROIs: {self.ctp404.roi_vals_as_str}\n'
                  f'HU Passed?: {self.ctp404.passed_hu}\n'
                  f'Low contrast visibility: {self.ctp404.lcv:2.2f}\n'
                  f'Geometric Line Average (mm): {self.ctp404.avg_line_length:2.2f}\n'
                  f'Geometry Passed?: {self.ctp404.passed_geometry}\n'
                  f'Measured Slice Thickness (mm): {self.ctp404.meas_slice_thickness:2.3f}\n'
                  f'Slice Thickness Passed? {self.ctp404.passed_thickness}\n')
        if self._has_module(CTP486):
            add = (f'Uniformity ROIs: {self.ctp486.roi_vals_as_str}\n'
                  f'Uniformity index: {self.ctp486.uniformity_index:2.3f}\n'
                  f'Integral non-uniformity: {self.ctp486.integral_non_uniformity:2.4f}\n'
                  f'Uniformity Passed?: {self.ctp486.overall_passed}\n')
            string += add
        if self._has_module(CTP528CP504):
            add = (f'MTF 50% (lp/mm): {self.ctp528.mtf.relative_resolution(50):2.2f}\n')
            string += add
        if self._has_module(CTP515):
            add = (f'Low contrast ROIs "seen": {self.ctp515.rois_visible}\n')
            string += add
        return string


class CatPhan503(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 503.
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528), Image Scaling & HU Linearity (CTP404).
    """
    _demo_url = 'CatPhan503.zip'
    _model = '503'
    catphan_radius_mm = 97
    modules = {
        CTP404CP503: {'offset': 0},
        CTP486: {'offset': -110},
        CTP528CP503: {'offset': -30},
    }

    @staticmethod
    def run_demo(show=True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan503.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


class CatPhan504(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 504. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """
    _demo_url = 'CatPhan504.zip'
    _model = '504'
    catphan_radius_mm = 101
    modules = {
        CTP404CP504: {'offset': 0},
        CTP486: {'offset': -65},
        CTP528CP504: {'offset': 30},
        CTP515: {'offset': -30}
    }

    @staticmethod
    def run_demo(show=True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan504.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


class CatPhan604(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 604. Can be from a CBCT or CT scanner
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """
    _demo_url = 'CatPhan604.zip'
    _model = '604'
    catphan_radius_mm = 101
    modules = {
        CTP404CP604: {'offset': 0},
        CTP486: {'offset': -80},
        CTP528CP604: {'offset': 42},
        CTP515: {'offset': -40}
    }

    @staticmethod
    def run_demo(show=True):
        """Run the CBCT demo using high-quality head protocol images."""
        cbct = CatPhan604.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


class CatPhan600(CatPhanBase):
    """A class for loading and analyzing CT DICOM files of a CatPhan 600.
    Analyzes: Uniformity (CTP486), High-Contrast Spatial Resolution (CTP528),
    Image Scaling & HU Linearity (CTP404), and Low contrast (CTP515).
    """
    _demo_url = 'CatPhan600.zip'
    _model = '600'
    catphan_radius_mm = 101
    modules = {
        CTP404CP600: {'offset': 0},
        CTP486: {'offset': -160},
        CTP515: {'offset': -110},
        CTP528CP600: {'offset': -70},
    }

    @staticmethod
    def run_demo(show=True):
        """Run the CatPhan 600 demo."""
        cbct = CatPhan600.from_demo_images()
        cbct.analyze()
        print(cbct.results())
        cbct.plot_analyzed_image(show)


def get_regions(slice_or_arr, fill_holes=False, clear_borders=True, threshold='otsu'):
    """Get the skimage regions of a black & white image."""
    if threshold == 'otsu':
        thresmeth = filters.threshold_otsu
    elif threshold == 'mean':
        thresmeth = np.mean
    if isinstance(slice_or_arr, Slice):
        edges = filters.scharr(slice_or_arr.image.array.astype(np.float))
        center = slice_or_arr.image.center
    elif isinstance(slice_or_arr, np.ndarray):
        edges = filters.scharr(slice_or_arr.astype(np.float))
        center = (int(edges.shape[1]/2), int(edges.shape[0]/2))
    edges = filters.gaussian(edges, sigma=1)
    if isinstance(slice_or_arr, Slice):
        box_size = 100/slice_or_arr.mm_per_pixel
        thres_img = edges[int(center.y-box_size):int(center.y+box_size),
                          int(center.x-box_size):int(center.x+box_size)]
        thres = thresmeth(thres_img)
    else:
        thres = thresmeth(edges)
    bw = edges > thres
    if clear_borders:
        segmentation.clear_border(bw, buffer_size=int(max(bw.shape)/50), in_place=True)
    if fill_holes:
        bw = ndimage.binary_fill_holes(bw)
    labeled_arr, num_roi = measure.label(bw, return_num=True)
    regionprops = measure.regionprops(labeled_arr, edges)
    return labeled_arr, regionprops, num_roi


@value_accept(mode=('mean', 'median', 'max'))
def combine_surrounding_slices(dicomstack, nominal_slice_num, slices_plusminus=1, mode='mean'):
    """Return an array that is the combination of a given slice and a number of slices surrounding it.

    Parameters
    ----------
    dicomstack : `~pylinac.core.image.DicomImageStack`
        The CBCT DICOM stack.
    nominal_slice_num : int
        The slice of interest (along 3rd dim).
    slices_plusminus: int
        How many slices plus and minus to combine (also along 3rd dim).
    mode : {'mean', 'median', 'max}
        Specifies the method of combination.

    Returns
    -------
    combined_array : numpy.array
        The combined array of the DICOM stack slices.
    """
    slices = range(nominal_slice_num - slices_plusminus, nominal_slice_num + slices_plusminus + 1)
    arrays = tuple(dicomstack[s].array for s in slices)
    array_stack = np.dstack(arrays)
    if mode == 'mean':
        combined_array = np.mean(array_stack, 2)
    elif mode == 'median':
        combined_array = np.median(array_stack, 2)
    else:
        combined_array = np.max(array_stack, 2)
    return combined_array
