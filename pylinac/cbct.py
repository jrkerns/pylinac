
"""The CBCT module automatically analyzes DICOM images of a CatPhan acquired when doing CBCT quality assurance. It can load a folder
the images are in and automatically correct for phantom setup by determining the yaw, pitch, and roll of the phantom.
It can analyze the HU regions and image scaling (CTP404), the high-contrast line pairs (CTP528) to calculate the modulation transfer function (MTF), and the HU
uniformity (CTP486) on the corresponding slice.

Currently only Varian (CatPhan 504) is supported, but Elekta (CatPhan 503) support is being worked on.
"""

import os
import os.path as osp
import zipfile as zp
import shutil

import numpy as np
from scipy import ndimage
from scipy.misc import imresize
import dicom
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept
from pylinac.core.image import ImageObj
from pylinac.core.common_functions import invert, dist_2points, peak_detect
from pylinac.core.geometry import Point, Circle, sector_mask


known_manufacturers = ('Varian Medical Systems', 'ELEKTA')

class ROI(object):

    def __init__(self, name, angle, radius_from_center):
        self.name = name
        self.angle = angle
        self.radius_from_center = radius_from_center


class HU_ROI(Circle, ROI):
    """An HU ROI object. Represents a circular area measuring either HU sample (Air, Poly, ...) or uniformity (bottom, left, ...)."""
    def __init__(self, name, angle, nominal_value, radius, radius_from_center):
        Circle.__init__(self, radius=radius)
        ROI.__init__(self, name, angle, radius_from_center)
        self.nominal_val = nominal_value

    def get_pixel_value(self):
        pass

    def compare_to_nominal(self):
        pass

    @property
    def value_diff(self):
        pass

class GEO_ROI(Circle, ROI):
    """A circle ROI, much like the HU ROI, but with methods to find the center of the geometric "node"."""
    def __init__(self, name, angle, radius, radius_from_center):
        Circle.__init__(self, radius=radius)
        ROI.__init__(self, name, angle, radius_from_center)

    def _threshold(self, threshold):
        pass

    def find_node_center(self):
        pass


class Slice(object):
    """An abstract subclass for analyzing specific slices of a CBCT dicom set."""
    def __init__(self, images, settings, slice, mode):
        """
        :param images: The CBCT image set as a 3D numpy array.
        :type images: ndarray
        :param slice: Slice number of interest
        :type slice: int
        :param mode: Mode of combining surrounding images to lower noise. Options: 'mean', 'max', 'median'
        :type mode: str
        """
        self.image = ImageObj(combine_surrounding_slices(images, slice, mode=mode))
        self.settings = settings
        self.phan_center = Point()  # center pixel of the phantom
        self.object_names = []
        self.object_angles = {}
        self.object_values = {}  # The pixel values of the ROIs
        self.object_valdiff = {}
        self.object_tolerance = 0
        self.object_nominal_values = {}  # The nominal values of the ROIs, if applicable
        self.object_center_pixel = {}
        self.radius2objs = 0
        self.object_radius = 0
        self.passed_test = True

    def set_obj_names(self, names):
        """Set the names of the Objects of Interest in the Slice.

        :param names: Names of ROIs
        :type names: tuple, list
        """
        self.object_names = names

    def set_obj_angles(self, angles):
        """Set the angles of the ROIs calculated from the center pixel in degrees. Directions: 0 is right, -90 is top, -90 is bottom."""
        for name, angle in zip(self.object_names, angles):
            self.object_angles[name] = angle

    def set_obj_values(self, values):
        """Set the values of the ROIs."""
        for name, value in zip(self.object_names, values):
            self.object_values[name] = value

    def set_obj_nominal_values(self, values):
        """Set the nominal values of the ROIs."""
        for name, value in zip(self.object_names, values):
            self.object_nominal_values[name] = value

    def set_radius2objs(self, radius):
        """Set the radius from center pixel to the ROIs."""
        self.radius2objs = radius

    def set_obj_radius(self, radius):
        """Set the radius of the ROIs."""
        self.object_radius = radius

    def set_obj_tolerance(self, tolerance):
        """Set the ROI passing tolerance."""
        self.object_tolerance = tolerance

    def find_phan_center(self, threshold):
        """Determine the location of the center of the phantom."""
        SOI_bw = array2logical(self.image, threshold)  # convert slice to binary based on threshold
        SOI_bw = ndimage.binary_fill_holes(SOI_bw)  # fill in air pockets to make one solid ROI
        SOI_labeled, no_roi = ndimage.label(SOI_bw)  # identify the ROIs
        if no_roi < 1 or no_roi is None:
            raise ValueError("Unable to locate the CatPhan")
        hist, bins = np.histogram(SOI_labeled, bins=no_roi)  # hist will give the size of each label
        SOI_bw_clean = np.where(SOI_labeled == np.argmax(hist), 1, 0)  # remove all ROIs except the largest one (the CatPhan)
        center_pixel = ndimage.measurements.center_of_mass(SOI_bw_clean)
        self.phan_center = Point(center_pixel)

    def _return_masked_ROI(self, roi_center, nan_outside_mask=True):
        # Create mask of ROI
        mask = sector_mask(self.image.pixel_array.shape, roi_center, self.object_radius, (0, 360))
        # Apply mask to image and replace all other values with NaNs.
        if nan_outside_mask:
            masked_img = np.where(mask == True, self.image, np.NaN)
        else:
            masked_img = np.where(mask == True, self.image, 0)
        return masked_img

    def _return_ROI_center(self, angle, name):
        # rename some things for convenience
        phan_cen_y = self.phan_center.y  # HU slice center y-coordinate
        phan_cen_x = self.phan_center.x  # ditto for x-coordinate
        radius = self.radius2objs  # radius from slice center to extract HU rois
        roll = self.settings['phantom roll']
        # Calculate the additional shift of the ROI coords, given the phantom roll. Exception is the Center ROI
        # of the Uniformity slice, which has a radius of 0.
        if name != 'Center':
            corrected_roll_x = np.cos(np.radians(roll + angle))
            corrected_roll_y = np.sin(np.radians(roll + angle))
        else:
            corrected_roll_x = 0
            corrected_roll_y = 0
        # Calculate the center pixel of the ROI
        roi_center = [phan_cen_x + radius * corrected_roll_x, phan_cen_y + radius * corrected_roll_y]
        return roi_center

    def calc_OOI(self, calc_diff=True):
        """Calculate the mean HU value of the 7 ROI areas."""
        for (name, angle) in self.object_angles.items():
            self.object_center_pixel[name] = self._return_ROI_center(angle, name)
            masked_img = self._return_masked_ROI(self.object_center_pixel[name])
            # Calculate the mean HU value and the difference from nominal
            self.object_values[name] = np.nanmean(masked_img)

            if calc_diff:
                # Calculate different from the nominal value if applicable
                self.object_valdiff[name] = self.object_values[name] - self.object_nominal_values[name]

        if calc_diff:
            self.check_objs_passed()

    def check_objs_passed(self):
        """Check whether the ROIs passed within tolerance."""
        for val in self.object_valdiff.values():
            if val < -self.object_tolerance or val > self.object_tolerance:
                self.passed_test = False


class HU(Slice):
    """Class for analysis of the HU slice of the CBCT dicom data set."""
    object_names = ('Air', 'PMP', 'LDPE', 'Poly', 'Acrylic', 'Delrin', 'Teflon')
    object_nominal_values = (-1000, -200, -100, -35, 120, 340, 990)
    object_angles = (-90, -120, 180, 120, 60, 0, -60)
    radius2objs = 120
    object_radius = 9
    tolerance = 40
    air_bubble_size = 450

    def __init__(self, images, settings, mode):
        Slice.__init__(self, images, settings, settings['HU slice'], mode)
        self.set_obj_names(HU.object_names)
        self.set_obj_angles(HU.object_angles)
        self.set_obj_nominal_values(HU.object_nominal_values)
        self.set_radius2objs(HU.radius2objs * settings['scaling ratio'])
        self.set_obj_radius(HU.object_radius * settings['scaling ratio'])
        self.set_obj_tolerance(HU.tolerance)

    def return_phan_roll(self):
        """Determine the "roll" of the phantom based on the two air bubbles in the HU slice."""
        # convert slice to logical
        SOI = array2logical(self.image, self.settings['BW threshold'])
        # invert the SOI; this makes the Air == 1 and phantom == 0
        SOI_inv = invert(SOI)
        # determine labels and number of rois of inverted SOI
        labels, no_roi = ndimage.measurements.label(SOI_inv)
        # calculate ROI sizes of each label TODO: simplify the air bubble-finding
        roi_sizes = [ndimage.measurements.sum(SOI_inv, labels, index=item) for item in range(1, no_roi + 1)]
        # extract air bubble ROIs (based on size threshold)
        bubble_thresh = self.air_bubble_size * self.settings['scaling ratio'] ** 2
        air_bubbles = [idx + 1 for idx, item in enumerate(roi_sizes) if
                       item < bubble_thresh * 1.5 and item > bubble_thresh / 1.5]
        # if the algo has worked correctly, it has found 2 and only 2 ROIs (the air bubbles)
        if len(air_bubbles) == 2:
            air_bubble_CofM = ndimage.measurements.center_of_mass(SOI_inv, labels, air_bubbles)
            y_dist = air_bubble_CofM[0][0] - air_bubble_CofM[1][0]
            x_dist = air_bubble_CofM[0][1] - air_bubble_CofM[1][1]
            angle = np.arctan2(y_dist, x_dist) * 180 / np.pi
            if angle < 0:
                roll = abs(angle) - 90
            else:
                roll = angle - 90
            phan_roll = roll
        else:
            phan_roll = 0
            print("Warning: CBCT phantom roll unable to be determined; assuming 0")

        return phan_roll


class Locon(Slice):
    """Class for analysis of the low contrast slice of the CBCT dicom data set."""
    # TODO: work on this
    def __init__(self, images, settings, mode):
        Slice.__init__(self, images, settings, settings['LC slice'], mode)


class SR(Slice):
    """Class for analysis of the HU slice of the CBCT dicom data set."""
    object_names = (0.2, 0.4, 0.6, 0.8, 1, 1.2)
    radius2objs = np.arange(95, 100)

    def __init__(self, images, settings, mode):
        Slice.__init__(self, images, settings, settings['SR slice'], mode)
        self.set_obj_names(SR.object_names)
        self.set_radius2objs(SR.radius2objs * settings['scaling ratio'])

    def _return_LP_profile(self):
        """Extract circular profiles of the Line-Pairs on the SR slice. Will extract 5 profiles, then take the median to make 1 profile.

        :returns: 1-D profile of all Line Pairs. Plot this for a nice view of all line pairs.
        """
        # rename some things for convenience
        phan_cent_y = self.phan_center.y  # HU slice center y-coordinate
        phan_cent_x = self.phan_center.x  # ditto for x-coordinate
        roll = self.settings['phantom roll']
        # create index and cos, sin points which will be the circle's rectilinear coordinates
        circle_idx = np.radians(np.arange(180 + roll, 360 - 0.01 + 180 + roll, 0.01)[::-1])
        profs = np.zeros((len(SR.radius2objs), len(circle_idx)))
        for row, radius in enumerate(self.radius2objs):
            x = np.cos(circle_idx) * radius + phan_cent_x
            y = np.sin(circle_idx) * radius + phan_cent_y
            # this scipy function pulls the values of the image along the y,x points defined above
            profs[row, :] = ndimage.map_coordinates(self.image, [y, x], order=0)

        median_profile = np.median(profs - profs.min(), 0)  # take the median of the grounded profile
        return median_profile

    def _find_LP_peaks(self, profile):
        """Find the peaks along the line pair profile extracted. Because of the varying width of lead/no lead, 3 searches are done
        with varying widths of peak spacing. This is to ensure that only 1 peak is found for the larger LPs, but does not pass over the
        thinner LPs further down the profile.

        :param profile: 1-D profile of the Line Pairs (normally from what is returned by return_LP_profile).
        :returns: values of peaks, indices of peaks. Will be 1-D numpy arrays.
        """

        region_1_bound = 4000  # approximate index between 1st and 2nd LP regions
        region_2_bound = 10500  # approximate index between 4th and 5th LP regions
        # region_3_bound = 17500  # approximate index between 8th and 9th LP regions; after this, regions become very hard to distinguish
        region_3_bound = 12300  # for head this can be 17500, but for thorax (low dose => low quality), we can only sample the first
        # 5 line pairs accurately
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
        """Find the line pair valleys. This is done by passing the indices of the peaks.
        The valleys are searched only between these peaks.

        :param profile: 1-D profile of the Line Pairs (normally from what is returned by return_LP_profile).
        :param max_idxs: 1-D array containing the indices of peaks
        """
        idx2del = np.array((1, 4, 7, 11))
        min_vals = np.zeros(16)
        min_idxs = np.zeros(16)
        for idx in range(len(max_idxs) - 1):
            min_val, min_idx = peak_detect(profile[max_idxs[idx]:max_idxs[idx + 1]],
                                           x=np.arange(max_idxs[idx], max_idxs[idx + 1]),
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
        averaging the pixel values of the peaks/valleys found.
        """
        num_peaks = np.array((0, 2, 3, 3, 4, 4, 4)).cumsum()
        num_valleys = np.array((0, 1, 2, 2, 3, 3, 3)).cumsum()
        for name, LP_pair in zip(self.object_names, range(len(num_peaks) - 1)):
            region_max = max_vals[num_peaks[LP_pair]:num_peaks[LP_pair + 1]].max()
            region_min = min_vals[num_valleys[LP_pair]:num_valleys[LP_pair + 1]].min()
            self.object_values[name] = (region_max - region_min) / (region_max + region_min)
        # normalize the values by the first LP
        max_mtf = np.array(list(self.object_values.values())).max()
        for name, value in self.object_values.items():
            self.object_values[name] /= max_mtf

    def calc_OOI(self):
        """Calculate the line pairs of the SR slice."""
        LP_profile = self._return_LP_profile()
        max_vals, max_idxs = self._find_LP_peaks(LP_profile)
        min_vals, min_idxs = self._find_LP_valleys(LP_profile, max_idxs)
        self._calc_MTF(max_vals, min_vals)

    @value_accept(percent=(60, 95))
    def get_MTF(self, percent=80):
        """Return the MTF value for the percent passed in.

        :param percent: The line-pair/mm value for the given MTF percentage.
        :type percent: int
        """
        # calculate x and y interpolations from Line Pair values and from the MTF measured
        x_vals_intrp = np.arange(self.object_names[0], self.object_names[-1], 0.01)
        x_vals = np.array(sorted(self.object_values.keys()))
        y_vals = np.array(sorted(self.object_values.values())[::-1])
        y_vals_intrp = np.interp(x_vals_intrp, x_vals, y_vals)

        mtf_percent = x_vals_intrp[np.argmin(np.abs(y_vals_intrp - (percent / 100)))]
        return mtf_percent

class UNIF(Slice):
    """Class for analysis of the Uniformity slice of the CBCT dicom data set."""
    object_names = ('Center', 'Right', 'Top', 'Left', 'Bottom')
    object_angles = (0, 0, -90, 180, 90)  # The center has an angle of 0, but it doesn't matter; the radius for that ROI is 0
    object_nominal_values = (0, 0, 0, 0, 0)
    radius2objs = 110
    object_radius = 20
    tolerance = 40

    def __init__(self, images, settings, mode):
        Slice.__init__(self, images, settings, settings['UN slice'], mode)
        self.set_obj_names(UNIF.object_names)
        self.set_obj_angles(UNIF.object_angles)
        self.set_obj_nominal_values(UNIF.object_nominal_values)
        self.set_radius2objs(UNIF.radius2objs * settings['scaling ratio'])
        self.set_obj_radius(UNIF.object_radius * settings['scaling ratio'])
        self.set_obj_tolerance(UNIF.tolerance)


class GEO(Slice):
    """Class for analysis of the Geometry slice of the CBCT dicom data set."""
    object_names = ('Top-Left', 'Top-Right', 'Bottom-Right', 'Bottom-Left')
    dist_names = ('Top-Horiz', 'Right-Vert','Bottom-Horiz', 'Left-Vert')
    object_angles = (-135, -45, 45, 135)  # The center has an angle of 0, but it doesn't matter; the radius for that ROI is 0
    roi_nominal_value = 50
    radius2objs = 72
    object_radius = 20
    tolerance = 1

    def __init__(self, images, settings, mode):
        Slice.__init__(self, images, settings, settings['HU slice'], mode)
        self.set_obj_names(GEO.object_names)
        self.set_obj_angles(GEO.object_angles)
        self.set_radius2objs(GEO.radius2objs * settings['scaling ratio'])
        self.set_obj_radius(GEO.object_radius * settings['scaling ratio'])
        self.set_obj_tolerance(GEO.tolerance)
        # map the lines to the nodes that they connect to
        self.node_map = {}
        for idx, dname in enumerate(GEO.dist_names):
            try:
                self.node_map[dname] = (GEO.object_names[idx], GEO.object_names[idx+1])
            except IndexError:
                self.node_map[dname] = (GEO.object_names[idx], GEO.object_names[0])

    def calc_node_CoM(self):
        """Calculate the center-of-mass of the geometric nodes."""
        for (name, angle) in self.object_angles.items():
            self.object_center_pixel[name] = self._return_ROI_center(angle, name)
            masked_img = self._return_masked_ROI(self.object_center_pixel[name], nan_outside_mask=False)
            # threshold image
            GEO_bw1 = np.where(masked_img > np.nanmedian(masked_img) * 1.4, 1, 0)
            GEO_bw2 = np.where(masked_img < np.nanmedian(masked_img) * 0.6, 1, 0)
            GEO_bw = GEO_bw1 + GEO_bw2
            # find center of the geometric node and number of ROIs found
            label, no_roi = ndimage.measurements.label(GEO_bw)
            if no_roi != 1:
                raise ValueError("Did not find the geometric node.")
            # determine the center of mass of the geometric node
            geo_center = ndimage.measurements.center_of_mass(GEO_bw, label)
            self.object_center_pixel[name] = geo_center

    def calc_node_distances(self):
        """Calculate the distances from node to node."""
        for dname, (node1, node2) in self.node_map.items():
            self.object_values[dname] = dist_2points(self.object_center_pixel[node1], self.object_center_pixel[node2]) * self.settings['mm/Pixel']
            self.object_valdiff[dname] = self.object_values[dname] - GEO.roi_nominal_value


class CBCT(object):
    """A class for loading and analyzing Cone-Beam CT DICOM files of a CatPhan 504 (Varian; Elekta 503 is being developed.
    Analyzes: Uniformity, Spatial Resolution, Image Scaling & HU Linearity.
    """
    BW_threshold = -800  # The threshold to apply to images to convert to B&W, which will highlight the CatPhan.

    def __init__(self):
        """Image properties and algorithm attrs."""
        self.settings = {}  # dict of settings the algo will use to compute the HU, SR, etc.
        self._phan_roll = 0  # the angle in degrees of the roll of the phantom from a perfect setup.
        self._roll_found = False  # boolean specifying whether the algo successfully found the roll.

    def load_demo_images(self):
        """Load the CBCT demo images."""
        demos_folder = osp.join(osp.split(osp.abspath(__file__))[0], 'demo_files')
        demo_folder = osp.join(demos_folder, 'High quality head')


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
        if len(file_list) < 20:
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
        self.images *= self.settings['Rescale Slope']
        self.images += self.settings['Rescale Intercept']

    def _get_settings(self, dicom_file):
        """Based on the images input, determine the starting conditions for the analysis. There are various manufacturers, protocols,
        and field-of-views for CBCTs. Based on this, the starting conditions can be determined.

        :param dicom_file: dicom file as read by pydicom
        """
        dcm = dicom_file
        """determine algorithm conditions"""
        # determine scaling ratio; this is basically a magnification factor for the following thresholds and radii and so forth.
        scaling_ratio = 0.488 / dcm.PixelSpacing[0]
        self.settings = {'BW threshold': CBCT.BW_threshold,
                         'scaling ratio': scaling_ratio,
                         'Rescale Slope': dcm.RescaleSlope,
                         'Rescale Intercept': dcm.RescaleIntercept,
                         'mm/Pixel': dcm.PixelSpacing[0],
                         'Manufacturer': dcm.Manufacturer,
                         'phantom roll': 0}
        # Set slice locations based on manufacturer
        if dcm.Manufacturer in known_manufacturers:
            if dcm.Manufacturer == known_manufacturers[0]:
                self.settings['HU slice'] = 32
                self.settings['UN slice'] = 9
                self.settings['SR slice'] = 44
                self.settings['LC slice'] = 23
            else:
                raise NotImplementedError("Elekta not yet implemented")
        else:
            raise ValueError("Unknown Manufacturer")

    def find_roll(self):
        """Determine the roll of the phantom by calculating the angle of the two air bubbles on the HU slice."""
        self.HU = HU(self.images, self.settings, 'mean')
        self.settings['phantom roll'] = self.HU.return_phan_roll()

    def find_HU(self):
        """Determine HU values from HU slice. Averages 3 slices."""
        if not hasattr(self, 'HU'):
            self.HU = HU(self.images, self.settings, 'mean')
        self.HU.find_phan_center(self.settings['BW threshold'])
        self.HU.calc_OOI()

    def find_SR(self):
        """Determine the Spatial Resolution from the Line-Pair slice."""
        self.SR = SR(self.images, self.settings, 'max')
        self.SR.find_phan_center(self.settings['BW threshold'])
        self.SR.calc_OOI()

    def find_GEO(self):
        """Determine the geometric distortion of the images.
         The distortion is found by finding the 4 small points on the inner region of the HU slice
         and calculating their distances to each other. Nominal distance is 5cm apart.
        """
        self.GEO = GEO(self.images, self.settings, 'median')
        self.GEO.find_phan_center(self.settings['BW threshold'])
        self.GEO.calc_node_CoM()
        self.GEO.calc_node_distances()
        self.GEO.check_objs_passed()

    def find_UNIF(self):
        """Determine the HU uniformity from Uniformity slice by analyzing 5 ROIs: center, top, right, bottom, left."""
        self.UN = UNIF(self.images, self.settings, 'mean')
        self.UN.find_phan_center(self.settings['BW threshold'])
        self.UN.calc_OOI()

    def find_Locon(self):
        """Determine Low Contrast values."""
        self.Locon = Locon(self.images, self.settings, 'mean')

    def plot_analyzed_image(self):
        """Draw the ROIs and lines the calculations were done on or based on."""
        # create figure
        fig, ((UN_ax, HU_ax), (SR_ax, LOCON_ax)) = plt.subplots(2,2)

        # Uniformity objects
        UN_ax.imshow(self.UN.image)
        for (key, values) in self.UN.object_center_pixel.items():
            UN_ax.add_patch(Circle((values[0], values[1]), radius=self.UN.object_radius,
                                   fill=False, edgecolor='blue'))
        UN_ax.autoscale(tight=True)
        UN_ax.set_title('Uniformity Slice')

        # HU objects
        HU_ax.imshow(self.HU.image)
        for (key, value) in self.HU.object_center_pixel.items():
            HU_ax.add_patch(Circle((value[0], value[1]), radius=self.HU.object_radius,
                                   fill=False, edgecolor='blue'))
        HU_ax.autoscale(tight=True)
        HU_ax.set_title('HU & Geometric Slice')

        # GEO objects
        for node1, node2 in self.GEO.node_map.values(): # for 2 corner geo points...
            HU_ax.plot([self.GEO.object_center_pixel[node1][1],
                        self.GEO.object_center_pixel[node2][1]],
                       [self.GEO.object_center_pixel[node1][0],
                        self.GEO.object_center_pixel[node2][0]],
                       'black')

        # SR objects
        SR_ax.imshow(self.SR.image)
        phan_cent_y = self.SR.phan_center.x  # phantom center x-coord
        phan_cent_x = self.SR.phan_center.y  # ditto for y-coord
        for radius in [self.SR.radius2objs[0], self.SR.radius2objs[-1]]:
            SR_ax.add_patch(Circle((phan_cent_x, phan_cent_y), radius=radius,
                                   fill=False, edgecolor='blue'))
        SR_ax.autoscale(tight=True)
        SR_ax.set_title('Spatial Resolution Slice')

        # Locon objects
        LOCON_ax.imshow(self.Locon.image)
        LOCON_ax.set_title('Low Contrast (In Development)')

        # show it all
        plt.show()

    def return_results(self):
        """Return and print the results of the analysis as a string."""

        #TODO: make prettier
        print(' - CBCT QA Test - ')
        print('HU Regions: ', self.HU.object_values)
        print('HU Passed?: ', self.HU.passed_test)
        print('Uniformity: ', self.UN.object_values)
        print('Uniformity Passed?: ', self.UN.passed_test)
        print('MTF 80% (lp/mm): ', self.SR.get_MTF(80))
        print('Geometric distances: ', self.GEO.object_values)
        print('Geometry Passed?: ', self.GEO.passed_test)

    def analyze(self):
        """Single-method full analysis of CBCT DICOM files."""
        self.find_roll()
        self.find_HU()
        self.find_UNIF()
        self.find_GEO()
        self.find_SR()
        self.find_Locon()

    def run_demo(self):
        """Run the CBCT demo using the high-quality head protocol."""
        cbct = CBCT()
        cbct.load_demo_images()
        cbct.analyze()
        cbct.return_results()
        cbct.plot_analyzed_image()

def array2logical(array, threshold_value):
    """Return a 'logical' (binary) version of the input array based on a threshold.

    :param array: numerical array to be analyzed
    :param threshold_value: int or float specifying the threshold value. If an array value is below the
        threshold value, it is converted to 0, otherwise to 1.
    """
    return np.where(array >= threshold_value, 1, 0)

@value_accept(mode=('mean','median','max'))
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
    else:
        comb_slice = np.max(slices, 2)
    return comb_slice



# ----------------------------------------
# CBCT Demo
# ----------------------------------------
if __name__ == '__main__':
    CBCT().run_demo()
    # cb = CBCT()
    # cb.load_folder_UI()
    # cb.load_folder(r"C:\Users\JRKerns\Dropbox\Programming\MATLAB\Projects\UVC Suite\Data\Full Data Set\CBCT Images\ACB1\High quality head")
    # cb.analyze()
    # cb.return_results()
    # cb.plot_analyzed_image()