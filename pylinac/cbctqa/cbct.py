
from __future__ import print_function, division, absolute_import, unicode_literals
import os
import os.path as osp
from time import sleep
import zipfile as zp
import shutil

from future.builtins import range
import numpy as np
from scipy.ndimage import binary_fill_holes
from scipy import ndimage
import scipy.ndimage.measurements as meas
from scipy.misc import imresize
import dicom
from matplotlib.patches import Circle
import matplotlib.pyplot as plt

from pylinac.common.decorators import type_accept, value_accept
from pylinac.common.image_classes import MultiImageObject
from pylinac.common.common_functions import invert, dist_2points, sector_mask, peak_detect


"""Default constants"""
nominal_GEO = 50.0
tolerances = {'HU': 40, 'UN': 40, 'SR': None, 'GEO': 1}
slices_of_interest = {'HU': 32, 'UN': 9, 'SR': 44, 'LOCON': 23}
known_manufacturers = ('Varian Medical Systems', 'ELEKTA')
FOV_thresh = 300  # the field of view size threshold (in mm) between "small" and "large"
protocols = ('hi_head', 'thorax', 'pelvis')

class Slice(object):
    """A subclass for analyzing slices of a CBCT dicom set."""
    def __init__(self, images, slice_num, mode):
        self.image = combine_surrounding_slices(images, slice_num, mode=mode)

    def set_roi_names(self, names):
        """Set the names of the ROIs in the Slice.
        :param names: Names of ROIs
        :type names: tuple, list
        """
        self.roi_names = names

    def find_phan_center(self):
        """Determine the location of the center of the phantom."""
        SOI_bw = array2logical(self.image, self._algo_settings['BW threshold'])  # convert slice to binary based on threshold
        SOI_bw = binary_fill_holes(SOI_bw)  # fill in air pockets to make one solid ROI
        SOI_labeled, no_roi = meas.label(SOI_bw)  # identify the ROIs
        if no_roi < 1 or no_roi is None:
            raise ValueError("Unable to locate the CatPhan")
        hist, bins = np.histogram(SOI_labeled, bins=no_roi)  # hist will give the size of each label
        SOI_bw_clean = np.where(SOI_labeled == np.argmax(hist), 1, 0)  # remove all ROIs except the largest one (the CatPhan)
        center_pixel = meas.center_of_mass(SOI_bw_clean)


class CBCT(MultiImageObject):
    """
    A class for loading and analyzing Cone-Beam CT DICOM files of a CatPhan 503 (Elekta) & 504 (Varian). Analyzes: Uniformity,
    Spatial Resolution, Contrast, & HU Linearity.
    """
    # TODO: Low Contrast
    def __init__(self):
        MultiImageObject.__init__(self)
        """Image properties and algorithm attrs."""
        self._dcm_props = {}  # dict of relevant properties of the DICOM files
        self._algo_settings = {}  # dict of settings the algo will use to compute the HU, SR, etc.
        self._slice_centers = {}  # dict y,x point indicating the center of the phantom for that slice.
        self._phan_roll = 0  # the angle in degrees of the roll of the phantom from a perfect setup.
        self._roll_found = False  # boolean specifying whether the algo successfully found the roll.

        """HU settings and attrs."""
        # HU ROI angles
        self.HU_roi_angles = {'Air': -90, 'PMP': -120, 'LDPE': 180,
                              'Poly': 120, 'Acrylic': 60, 'Delrin': 0, 'Teflon': -60}
        # Nominal HU values, i.e. what they should be.
        self.HU_nominal_vals = {'Air': -1000, 'PMP': -200, 'LDPE': -100,
                                'Poly': -35, 'Acrylic': 120, 'Delrin': 340, 'Teflon': 990}
        # The measured HU values
        self.HU_actual_vals = {'Air': None, 'PMP': None, 'LDPE': None,
                               'Poly': None, 'Acrylic': None, 'Delrin': None, 'Teflon': None}
        # The difference in value from nominal to actual
        self.HU_diff_vals = {'Air': None, 'PMP': None, 'LDPE': None,
                               'Poly': None, 'Acrylic': None, 'Delrin': None, 'Teflon': None}
        # The box bounds are the pixel positions for determining the location of the ROIs
        self._HU_roi_centers = {'Air': None, 'PMP': None, 'LDPE': None,
                             'Poly': None, 'Acrylic': None, 'Delrin': None, 'Teflon': None}
        self.HU_tolerance = 40  # tolerance the measured HU must be within to be considered passing
        self.HU_passed = True

        """Spatial Resolution (SR) and Modulation Transfer Function (MTF) settings.
        The MTF is a measurement mechanism (and thus, subsidiary) of the spatial resolution."""
        self.SR_lpmm = np.array((0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.41, 1.59))  # the line pair/mm values for the CatPhan SR slice
        self.MTF_vals = np.zeros(8)  # The values of the MTF at the 8 SR ROIs
        # The y, x points of each corner of each SR ROI
        self.SR_ROI_yxpoints = [{'Top-Left': (0, 0), 'Top-Right': (0, 0),
                                'Bottom-Left': (0, 0), 'Bottom-Right': (0, 0)}] * 8

        """Geometrical Distance and Distortion."""
        self.GEO_dist = {'Top-Horiz': None, 'Bottom-Horiz': None,
                         'Left-Vert': None, 'Right-Vert': None}  # dictionary of distances of geometric nodes; will be 4 values
        self.GEO_diff = {'Top-Horiz': None, 'Bottom-Horiz': None,
                         'Left-Vert': None, 'Right-Vert': None}  # dictionary of differences of geometric distances from nominal
        self.GEO_nominal = 50.0
        self.GEO_tolerance = 1.0
        self._GEO_ROI_centers = {'Top-Left': (0, 0), 'Top-Right': (0, 0),
                                 'Bottom-Left': (0, 0), 'Bottom-Right': (0, 0)}
        self.GEO_passed = True

        """Uniformity"""
        self.UN_vals = {'Center': None, 'Top': None, 'Right': None,
                       'Bottom': None, 'Left': None}
        self.UN_diff = {'Center': None, 'Top': None, 'Right': None,
                        'Bottom': None, 'Left': None}
        self._UN_roi_bounds = {'Center': None, 'Top': None, 'Right': None,
                               'Bottom': None, 'Left': None}
        self._UN_angles = {'Bottom': 90, 'Top': -90, 'Left': 180, 'Right': 0, 'Center': 0}
        self._UN_roi_centers = {'Center': None, 'Top': None, 'Right': None,
                                'Bottom': None, 'Left': None}
        self.UNIF_passed = True

    @type_accept(protocol=str)
    @value_accept(protocol=protocols)
    def load_demo_images(self, protocol='hi_head'):
        """Load demo images based from a given protocol.

        :param protocol: Options include 'hi_head', 'thorax', 'pelvis'.
        :type protocol: str
        """

        demos_folder = osp.join(osp.split(osp.abspath(__file__))[0], 'demo_files')
        if protocol.lower() == protocols[0]:  # high quality head
            demo_folder = osp.join(demos_folder, 'High quality head')
        elif protocol.lower() == protocols[1]:  # thorax
            demo_folder = osp.join(demos_folder, 'Low dose thorax')
        elif protocol.lower() == protocols[2]:  # pelvis
            demo_folder = osp.join(demos_folder, 'Pelvis')

        # extract files from the zip file
        zp.ZipFile(demo_folder + '.zip').extractall(demos_folder)

        # walk through the extracted demo folder and pull the names of the CT dicom files
        for par_dir, sub_dir, demo_folder in os.walk(demo_folder):
            file_list = [osp.join(par_dir, item) for item in demo_folder if item.endswith('.dcm') and item.startswith('CT')]

        self._load_files(file_list, demo_files=True)

    def load_folder_UI(self):
        """Load the CT DICOM files from a folder using a UI."""
        folder = self.get_folder_UI()
        self.load_folder(folder)

    @type_accept(folder=str, append=bool)
    def load_folder(self, folder, append=False):
        """Load the CT DICOM files from a folder string input

        :param folder: Path to the folder
        :type folder: str
        :param append: Whether or not the images in the folder should be appended to those already loaded. Not yet implemented.
        :type append: bool
        """
        # check that folder is valid
        if not osp.isdir(folder):
            raise NotADirectoryError("Path given was not a Directory/Folder")

        for par_dir, sub_dir, files in os.walk(folder):
            filelist = [osp.join(par_dir, item) for item in files if item.endswith('.dcm') and item.startswith('CT')]
        self._load_files(filelist)

    def _sort_images(self, im_order, images):
        """Sort the images according to the image order."""
        #TODO: convert this to some sort of in-place method
        nimages = np.zeros(images.shape, dtype=int)
        for idx, orig in enumerate(im_order):
            nimages[:, :, orig] = images[:, :, idx]
        self.images = nimages

    def _load_files(self, file_list, demo_files=False, im_size=512):
        """Load CT DICOM files given a list of image paths.

        :param file_list: list of file paths to load
        :type file_list: list
        :param demo_files: Specifies whether the method is using the demo files; if so, it will delete the extracted folder to
            keep things clean.
        :type demo_files: bool
        :param im_size: int specifying the size of images the DICOM images should be resized to; used for simplicity of settings.
        """
        # initialize image array
        images = np.zeros([im_size, im_size, len(file_list)], dtype=int)


        im_order = np.zeros(len(file_list), dtype=int)
        # check that enough slices were loaded
        #TODO: select better method for determining if enough slices were chosen
        if len(file_list) < slices_of_interest['SR']:
            raise ValueError("Not enough slices were selected. Select all the files or change the SR slice value.")

        # load dicom files from list names and get the image slice position
        #TODO: figure out more memory-efficient way to sort images; maybe based on number of slices and slice thickness
        for idx, item in enumerate(file_list):
            dcm = dicom.read_file(item)
            im_order[idx] = np.round(dcm.ImagePositionPatient[-1]/dcm.SliceThickness)
            # resize image if need be
            if dcm.pixel_array.shape != (im_size, im_size):
                image = imresize(dcm.pixel_array,(im_size, im_size))
            else:
                image = dcm.pixel_array
            # place image into images array
            images[:, :, idx] = image
        # shift order list from -xxx:+xxx to 0:2xxx
        im_order += np.round(abs(min(im_order)))
        # resort the images to be in the correct order
        self._sort_images(im_order, images)

        # if using the demo files (which are originally a .zip file), delete the temporary directory they were unzipped to.
        if demo_files:
            try:
                for root, _, files in os.walk(osp.split(file_list[0])[0]):
                    for file in files:
                        os.remove(osp.join(root, file))
                sleep(0.05)  # see if adding a pause prevents rmtree errors
                os.rmdir(osp.split(file_list[0])[0])
            except OSError:
                shutil.rmtree(osp.split(file_list[0])[0])
            except:
                print("Demo files were extracted from zip but were not able to be deleted after loading. This must be done manually.")

        # determine settings needed for given CBCT.
        self._get_settings(dcm)
        # convert images from CT# to HU using dicom tags
        self._convert2HU()

    def _convert2HU(self):
        """Convert the images from CT# to HU."""
        self.images *= self._dcm_props['Rescale Slope']
        self.images += self._dcm_props['Rescale Intercept']

    def _get_settings(self, dicom_file):
        """Based on the images input, determine the starting conditions for the analysis. There are various manufacturers, protocols,
        and field-of-views for CBCTs. Based on this, the starting conditions can be determined.

        :param dicom_file: dicom file as read by pydicom
        """
        dcm = dicom_file
        # get general image information
        self._dcm_props = {'FOV': dcm.DataCollectionDiameter,
                        'Rescale Slope': dcm.RescaleSlope,
                        'Rescale Intercept': dcm.RescaleIntercept,
                        'mm/Pixel': dcm.PixelSpacing[0],
                        'Manufacturer': dcm.Manufacturer}
        """determine algorithm conditions"""
        # determine scaling ratio; this is basically a magnification factor for the following thresholds and radii and so forth.
        scaling_ratio = 0.488/self._dcm_props['mm/Pixel']
        self._algo_settings = {'Catphan size thresh': 100000 * scaling_ratio,
                              'pix thresh': 500,
                              'BW threshold': -800,

                              'GEO min size': 25 * scaling_ratio**2,
                              'GEO box width': 65 / self._dcm_props['mm/Pixel'],

                              'HU radius': 120 * scaling_ratio,
                              'HU air bubble size': 450 * scaling_ratio**2,
                              'HU ROI radius': 9 * scaling_ratio,

                              'CON angles': [0.4, -0.0 - 5, 2.03],
                              'CON radius': 100 * scaling_ratio,
                              'CON ROI size': 20 * scaling_ratio**2,

                              'SR radii': np.arange(95, 100, 1) * scaling_ratio,

                              'UN ROI radius': 20 * scaling_ratio,

                              'UN radius': 110 * scaling_ratio}

    def _find_phan_centers(self):
        """
        Find the center of the phantom for each slice that is analyzed. All algorithms are relative to the center point.
        The slice is first converted to a binary image based on a threshold, then all other ROIs but the Catphan are removed.
        Then the center of ROI is found.
        """
        for (slice_name, slice_num) in list(slices_of_interest.items()):
            SOI = self.images[:, :, slice_num]
            SOI_bw = array2logical(SOI, self._algo_settings['BW threshold'])  # convert slice to binary based on threshold
            SOI_bw = binary_fill_holes(SOI_bw)  # fill in air pockets to make one solid ROI
            SOI_labeled, no_roi = meas.label(SOI_bw)  # identify the ROIs
            if no_roi < 1 or no_roi is None:
                raise ValueError("Unable to locate the CatPhan")
            hist, bins = np.histogram(SOI_labeled, bins=no_roi) # hist will give the size of each label
            SOI_bw_clean = np.where(SOI_labeled == np.argmax(hist), 1, 0)  # remove all ROIs except the largest one (the CatPhan)
            self._slice_centers[slice_name] = meas.center_of_mass(SOI_bw_clean)

            # if this has in any way failed, set the center to the center pixel
            if np.isnan(self._slice_centers[slice_name][0]):
                self._slice_centers[slice_name] = (np.size(SOI,0), np.size(SOI, 1))

    def _find_roll(self):
        """
        Determine the roll of the phantom by calculating the angle of the two air bubbles on the HU slice.
        """
        SOI = combine_surrounding_slices(self.images, slices_of_interest['HU'])
        # convert SOI to logical
        SOI = array2logical(SOI, self._algo_settings['BW threshold'])
        # invert the SOI; this makes the Air == 1 and phantom == 0
        SOI_inv = invert(SOI)
        # determine labels and number of rois of inverted SOI
        labels, no_roi = meas.label(SOI_inv)
        # calculate ROI sizes of each label TODO: simplify the air bubble-finding
        roi_sizes = [meas.sum(SOI_inv, labels, index=item) for item in range(1,no_roi+1)]
        # extract air bubble ROIs (based on size threshold
        bubble_thresh = self._algo_settings['HU air bubble size']
        air_bubbles = [idx+1 for idx, item in enumerate(roi_sizes) if item < bubble_thresh*1.5 and item > bubble_thresh/1.5]
        # if the algo has worked correctly, it has found 2 and only 2 ROIs (the air bubbles)
        if len(air_bubbles) != 2:
            self._phan_roll = 0
            # raise RuntimeWarning("Roll unable to be determined; assuming 0")
            print("Warning: CBCT phantom roll unable to be determined; assuming 0")
        else:
            air_bubble_CofM = meas.center_of_mass(SOI_inv, labels, air_bubbles)
            y_dist = air_bubble_CofM[0][0] - air_bubble_CofM[1][0]
            x_dist = air_bubble_CofM[0][1] - air_bubble_CofM[1][1]
            angle = np.arctan2(y_dist, x_dist) * 180/np.pi
            if angle < 0:
                roll = abs(angle) - 90
            else:
                roll = angle - 90
            self._phan_roll = roll

            self._roll_found = True

    def find_HU(self):
        """Determine HU values from HU slice. Averages 3 slices."""

        # average 3 slices around the nominal slice
        HU_slice = combine_surrounding_slices(self.images, slices_of_interest['HU'])

        # For each HU ROI...
        for idx, (roi_key, val) in enumerate(list(self.HU_nominal_vals.items())):
            # rename some things for convenience
            phan_cen_y = self._slice_centers['HU'][0]  # HU slice center y-coordinate
            phan_cen_x = self._slice_centers['HU'][1]  # ditto for x-coordinate
            radius = self._algo_settings['HU radius']  # radius from slice center to extract HU rois

            # Calculate the additional shift of the ROI coords, given the phantom roll.
            corrected_roll_x = np.cos(np.radians(self._phan_roll + self.HU_roi_angles[roi_key]))
            corrected_roll_y = np.sin(np.radians(self._phan_roll + self.HU_roi_angles[roi_key]))

            # Calculate the center pixel of the ROI
            roi_center = [phan_cen_x + radius*corrected_roll_x, phan_cen_y + radius*corrected_roll_y]
            self._HU_roi_centers[roi_key] = roi_center

            mask = sector_mask(HU_slice.shape, roi_center, self._algo_settings['HU ROI radius'], (0, 360))
            masked_img = np.where(mask == True, HU_slice, np.NaN)

            # Calculate the mean HU value and the difference from nominal
            self.HU_actual_vals[roi_key] = np.nanmean(masked_img)
            self.HU_diff_vals[roi_key] = self.HU_actual_vals[roi_key] - self.HU_nominal_vals[roi_key]

        # Check if all ROIs passed tolerance
        for val in self.HU_diff_vals.values():
            if val < -self.HU_tolerance or val > self.HU_tolerance:
                self.HU_passed = False

    def _get_LP_profile(self, SR_slice):
        """Extract circular profiles of the Line-Pairs on the SR slice. Will extract 5 profiles, then take the median to make 1 profile.

        :param SR_slice: Image/Slice to extract profiles from.
        :returns: 1-D profile of all Line Pairs. Plot this for a nice view of all line pairs.
        """
        # rename a few things for convenience
        phan_cent_x = self._slice_centers['SR'][0]  # phantom center x-coord
        phan_cent_y = self._slice_centers['SR'][1]  # ditto for y-coord
        # create index and cos, sin points which will be the circle's rectilinear coordinates
        deg = np.arange(180 + self._phan_roll, 360 - 0.01 + 180 + self._phan_roll, 0.01)[::-1]
        profs = np.zeros((len(self._algo_settings['SR radii']), len(deg)))
        for row, radius in enumerate(self._algo_settings['SR radii']):
            x = np.cos(np.deg2rad(deg)) * radius + phan_cent_x
            y = np.sin(np.deg2rad(deg)) * radius + phan_cent_y
            # this scipy function pulls the values of the image along the y,x points defined above
            profs[row, :] = ndimage.map_coordinates(SR_slice, [y, x], order=0)

        m_prof = np.median(profs - profs.min(), 0)  # take the median of the grounded profile
        return m_prof

    def _find_LP_peaks(self, profile):
        """Find the peaks along the line pair profile extracted. Because of the varying width of lead/no lead, 3 searches are done
        with varying widths of peak spacing. This is to ensure that only 1 peak is found for the larger LPs, but does not pass over the
        thinner LPs further down the profile."""

        region_1_bound = 4000  # approximate index between 1st and 2nd LP regions
        region_2_bound = 10500  # approximate index between 4th and 5th LP regions
        # region_3_bound = 17500  # approximate index between 8th and 9th LP regions; after this, regions become very hard to distinguish
        region_3_bound = 12300  # for head this can be 17500, but for thorax (low dose => low quality), we can only sample the first
            # 5 line pairs accurately
        # Find 1st LP max vals and idxs
        max_vals_1, max_idx_1 = peak_detect(profile[:region_1_bound],
                                            threshold=0.2, min_peak_width=600, max_num_peaks=2)
        max_vals_2, max_idx_2 = peak_detect(profile[region_1_bound:region_2_bound],
                                            x=np.arange(region_1_bound, region_2_bound),
                                            threshold=0.2, min_peak_width=250)
        max_vals_3, max_idx_3 = peak_detect(profile[region_2_bound: region_3_bound],
                                            x=np.arange(region_2_bound, region_3_bound),
                                            threshold=0.2, min_peak_width=100)
        max_vals = np.concatenate((max_vals_1, max_vals_2, max_vals_3))
        max_idxs = np.concatenate((max_idx_1, max_idx_2, max_idx_3))
        if len(max_idxs) != 17:
            raise ArithmeticError("Did not find the correct number of line pairs")
        return max_vals, max_idxs

    def _find_LP_valleys(self, profile, max_idxs):
        """Find the line pair valleys. This has to be done carefully to avoid catching valleys between LP regions; this is done by
        passing the indices of the peaks. The valleys are searched only between these peaks."""

        idx2del = np.array((1, 4, 7, 11))
        min_vals = np.zeros(16)
        min_idxs = np.zeros(16)
        for idx in range(len(max_idxs)-1):
            min_val, min_idx = peak_detect(profile[max_idxs[idx]:max_idxs[idx+1]], x=np.arange(max_idxs[idx], max_idxs[idx + 1]),
                        max_num_peaks=1, threshold=0.05, find_min_instead=True)
            min_vals[idx] = min_val
            min_idxs[idx] = min_idx
        # now delete the valleys *in between* the LP regions
        min_vals = np.delete(min_vals, idx2del)
        min_idxs = np.delete(min_idxs, idx2del)
        return min_vals, min_idxs

    def _calc_MTF(self, max_vals, min_vals):
        """Calculate the Modulation Transfer Function of the Line-Pair profile.
        See http://en.wikipedia.org/wiki/Transfer_function#Optics for calculation. Maximum and minimum values are calculated by
        averaging the pixel values of the peaks/valleys.
        """
        num_peaks = np.array((0, 2, 3, 3, 4, 4, 4, 5, 5)).cumsum()
        num_valleys = np.array((0, 1, 2, 2, 3, 3, 3, 4, 4)).cumsum()
        for region in range(len(num_peaks)-1):
            region_max = max_vals[num_peaks[region]:num_peaks[region+1]].mean()
            region_min = min_vals[num_valleys[region]:num_valleys[region + 1]].mean()
            self.MTF_vals[region] = (region_max - region_min) / (region_max + region_min)


    def find_SR(self):
        """Determine the Spatial Resolution from the Line-Pair slice."""
        # average 3 slices centered around SR slice
        SR_slice = combine_surrounding_slices(self.images, slices_of_interest['SR'], mode='max')

        prof = self._get_LP_profile(SR_slice)

        max_vals, max_idxs = self._find_LP_peaks(prof)
        min_vals, min_idxs = self._find_LP_valleys(prof, max_idxs)
        self._calc_MTF(max_vals, min_vals)

    @type_accept(percent=int)
    @value_accept(percent=(40, 80))
    def get_MTF(self, percent=50):
        """Return the MTF value for the percent passed in.

        :param percent: The line-pair/mm value for the given MTF percentage.
        :type percent: int
        """

        # calculate x and y interpolations from Line Pair values and from the MTF measured
        x_vals_intrp = np.arange(self.SR_lpmm[0], self.SR_lpmm[-1], 0.01)
        y_vals_intrp = np.interp(x_vals_intrp, self.SR_lpmm, self.MTF_vals)

        mtf_val = x_vals_intrp[np.argmin(np.abs(y_vals_intrp - percent/100))]
        return mtf_val

    def find_GEO(self):
        """Determine the geometric distortion of the images.
         The distortion is found by finding the 4 small points on the inner region of the HU slice
         and calculating their distances to each other. Nominal distance is 5cm apart.
        """
        # median 3 slices around GEO slice, which is the HU slice)
        GEO_slice = combine_surrounding_slices(self.images, slices_of_interest['HU'], mode='median')

        # Rename for convenience
        left_bound = self._slice_centers['HU'][1] - self._algo_settings['GEO box width'] / 2
        right_bound = self._slice_centers['HU'][1] + self._algo_settings['GEO box width'] / 2
        top_bound = self._slice_centers['HU'][0] - self._algo_settings['GEO box width'] / 2
        bottom_bound = self._slice_centers['HU'][0] + self._algo_settings['GEO box width'] / 2

        # extract zoomed-in section of GEO slice (around the 4 geometric markers).
        GEO_slice_zoomed = GEO_slice[left_bound:right_bound, top_bound:bottom_bound]

        # construct black & white image
        GEO_bw1 = np.where(GEO_slice_zoomed > np.median(GEO_slice_zoomed) * 1.4, 1, 0)
        GEO_bw2 = np.where(GEO_slice_zoomed < np.median(GEO_slice_zoomed) * 0.6, 1, 0)
        GEO_bw = GEO_bw1 + GEO_bw2

        # find centers of geometric points
        labels, no_roi = meas.label(GEO_bw)
        # determine the size of the ROIs found
        roi_sizes = np.array([meas.sum(GEO_bw, labels, index=item) for item in range(1, no_roi + 1)])
        if len(roi_sizes) <  4:
            raise ValueError("Unable to locate the geometric nodes. May need to adjust the geometric ROI slice.")
        # grab indices of largest 4
        geo_node_indx = np.argsort(roi_sizes)[-4:]
        geo_CofM = meas.center_of_mass(GEO_bw, labels, index=geo_node_indx+1)
        #TODO: add calc of geo point centers
        geo_CofM2 = [(x+left_bound, y+top_bound) for (x, y) in geo_CofM]
        self.geo_CofM = geo_CofM2

        # distance calculations (result is in mm; nominal distance is 50mm (5cm)
        # for key, value in list(self.GEO_dist.items()):
        self.GEO_dist['Top-Horiz'] = dist_2points(geo_CofM[0], geo_CofM[1]) * self._dcm_props['mm/Pixel']
        self.GEO_dist['Left-Vert'] = dist_2points(geo_CofM[0], geo_CofM[2]) * self._dcm_props['mm/Pixel']
        self.GEO_dist['Right-Vert'] = dist_2points(geo_CofM[3], geo_CofM[1]) * self._dcm_props['mm/Pixel']
        self.GEO_dist['Bottom-Horiz'] = dist_2points(geo_CofM[3], geo_CofM[2]) * self._dcm_props['mm/Pixel']
        self.GEO_diff['Top-Horiz'] = self.GEO_dist['Top-Horiz'] - self.GEO_nominal
        self.GEO_diff['Left-Vert'] = self.GEO_dist['Left-Vert'] - self.GEO_nominal
        self.GEO_diff['Right-Vert'] = self.GEO_dist['Right-Vert'] - self.GEO_nominal
        self.GEO_diff['Bottom-Horiz'] = self.GEO_dist['Bottom-Horiz'] - self.GEO_nominal

        # check if passed
        for val in self.GEO_dist.values():
            if (val > self.GEO_nominal + self.GEO_tolerance) or (val < self.GEO_nominal - self.GEO_tolerance):
                self.GEO_passed = False

    def find_UNIF(self):
        """
        Determine Uniformity from Uniformity slice by analyzing 5 ROIs: center, top, right, bottom, left
        """
        # mean from 3 slices around nominal slice
        UNIF_slice = combine_surrounding_slices(self.images, slices_of_interest['UN'])

        # calculate ROI centers & HU values thereof
        for key in self.UN_vals:
            if key == 'Center':
                self._UN_roi_centers[key] = (self._slice_centers['UN'][1], self._slice_centers['UN'][0])
            else:
                x_shift = self._algo_settings['UN radius'] * np.cos(np.radians(self._phan_roll + self._UN_angles[key]))
                y_shift = self._algo_settings['UN radius'] * np.sin(np.radians(self._phan_roll + self._UN_angles[key]))
                self._UN_roi_centers[key] = (self._slice_centers['UN'][1] + x_shift,
                                             self._slice_centers['UN'][0] + y_shift)

            mask = sector_mask(UNIF_slice.shape, self._UN_roi_centers[key], self._algo_settings['UN ROI radius'], (0, 360))
            masked_img = np.where(mask == True, UNIF_slice, np.NaN)

            # Calculate the mean UNIF value and the difference from nominal
            self.UN_vals[key] = np.nanmean(masked_img)

        # check if values passed
        for value in self.UN_vals.values():
            if value < -self.HU_tolerance or value > self.HU_tolerance:
                self.UNIF_passed = False

    def plot_analyzed_image(self):
        """Draw the ROIs and lines the calculations were done on or based on."""
        # create figure
        fig, ((UN_ax, HU_ax), (SR_ax, LOCON_ax)) = plt.subplots(2,2)
        # Uniformity objects
        UN_ax.imshow(self.images[:,:,slices_of_interest['UN']])
        for (key, values) in self._UN_roi_centers.items():
            UN_ax.add_patch(Circle((values[0], values[1]), radius=self._algo_settings['UN ROI radius'],
                                   fill=False, edgecolor='blue'))
        UN_ax.autoscale(tight=True)
        UN_ax.set_title('Uniformity Slice')

        # HU objects
        HU_ax.imshow(self.images[:, :, slices_of_interest['HU']])
        for (key, value) in self._HU_roi_centers.items():
            HU_ax.add_patch(Circle((value[0], value[1]), radius=self._algo_settings['HU ROI radius'],
                                   fill=False, edgecolor='blue'))
        HU_ax.autoscale(tight=True)
        HU_ax.set_title('HU & Geometric Slice')

        # GEO objects
        for idx in [0,3]: # for 2 corner geo points...
            for idx2 in [1,2]:  # connect to the other 2 geo points
                HU_ax.plot([self.geo_CofM[idx][1],
                            self.geo_CofM[idx2][1]],
                           [self.geo_CofM[idx][0],
                            self.geo_CofM[idx2][0]],
                           'black')

        # SR objects
        SR_ax.imshow(self.images[:,:,slices_of_interest['SR']])
        phan_cent_y = self._slice_centers['SR'][0]  # phantom center x-coord
        phan_cent_x = self._slice_centers['SR'][1]  # ditto for y-coord
        for radius in [self._algo_settings['SR radii'][0], self._algo_settings['SR radii'][-1]]:
            SR_ax.add_patch(Circle((phan_cent_x, phan_cent_y), radius=radius,
                                   fill=False, edgecolor='blue'))
        SR_ax.autoscale(tight=True)
        SR_ax.set_title('Spatial Resolution Slice')

        #TODO: Low contrast
        LOCON_ax.imshow(self.images[:,:,slices_of_interest['LOCON']])
        LOCON_ax.set_title('Low Contrast (In Development)')

        plt.show()

    def return_results(self):
        """Return and print the results of the analysis as a string."""

        #TODO: make a bit prettier
        print('HU Regions: ', self.HU_actual_vals)
        print('HU Passed?: ', self.HU_passed)
        print('Uniformity: ', self.UN_vals)
        print('Uniformity Passed?: ', self.UNIF_passed)
        print('MTF 50% (lp/mm): ', self.get_MTF())
        print('Geometric distances: ', self.GEO_dist)
        print('Geometry Passed?: ', self.GEO_passed)

    def analyze(self):
        """One-method for full analysis of CBCT DICOM files."""
        self._find_phan_centers()
        self._find_roll()

        self.find_HU()
        self.find_SR()
        self.find_GEO()
        self.find_UNIF()

def array2logical(array, threshold_value):
    """Return a 'logical' (binary) version of the input array.
    :param array: numerical array to be analyzed
    :param threshold_value: int or float specifying the threshold value. If an array value is below the
        threshold value, it is converted to 0, otherwise to 1.
    """
    return np.where(array >= threshold_value, 1, 0)

@type_accept(mode=str)
def combine_surrounding_slices(im_array, nominal_slice_num, slices_plusminus=1, mode='mean'):
    """Return an array that is the combination of a given slice and a number of slices surrounding it.
    :param im_array: numpy array, assumed 3 dims
    :param nominal_slice_num: int specifying the slice of interest (along 3rd dim)
    :param slices_plusminus: int specifying how many slices plus and minus to combine (also along 3rd dim)
    :param mode: string specifying the method of combination. Options are: 'mean', 'median', or 'max'.
    """
    slices = im_array[:,:,nominal_slice_num-slices_plusminus:nominal_slice_num+slices_plusminus]
    if mode == 'mean':
        comb_slice = np.mean(slices, 2)
    elif mode == 'median':
        comb_slice = np.median(slices, 2)
    elif mode == 'max':
        comb_slice = np.max(slices, 2)
    else:
        raise KeyError("Slice combination key invalid")
    return comb_slice



# ----------------------------------------
# CBCT Demo
# ----------------------------------------
if __name__ == '__main__':
    cbct = CBCT()
    cbct.load_demo_images('pelvis')
    cbct.analyze()
    cbct.return_results()
    cbct.plot_analyzed_image()
