
"""The CBCT module automatically analyzes DICOM images of a CatPhan acquired when doing CBCT quality assurance. It can load a folder
the images are in and automatically correct for phantom setup by determining the yaw, pitch, and roll of the phantom.
It can analyze the HU regions and image scaling (CTP404), the high-contrast line pairs (CTP528) to calculate the modulation transfer function (MTF), and the HU
uniformity (CTP486) on the corresponding slice.

Currently only Varian (CatPhan 504) is supported, but Elekta (CatPhan 503) support is being worked on.
"""
from abc import ABCMeta
from collections import OrderedDict
import os
import os.path as osp
import shutil

import numpy as np
from scipy import ndimage
from scipy.misc import imresize
import dicom
import matplotlib.pyplot as plt

from pylinac.core.decorators import value_accept, lazyproperty, type_accept
from pylinac.core.image import ImageObj
from pylinac.core.geometry import Point, Circle, sector_mask, Line
from pylinac.core.profile import CircleProfile, Profile
from pylinac.core.io import get_folder_UI


np.seterr(invalid='ignore')  # ignore warnings for invalid numpy operations. Used for np.where() operations on partially-NaN arrays.

known_manufacturers = ('Varian Medical Systems', 'ELEKTA')


class ROI(metaclass=ABCMeta):
    """Abstract base class for CBCT regions of interest."""
    def __init__(self, name, img_array):
        """
        :param name: Name of the ROI
        :type name: str
        :param img_array: 2D array the ROI is on.
        :type img_array: numpy.ndarray
        """
        self.name = name
        self.img_array = img_array


class ROI_Disk(Circle, ROI):
    """An abstract class representing a circular/disk Region of Interest on a CBCT slice."""
    def __init__(self, name, img_array, angle, radius=None, dist_from_center=None):
        ROI.__init__(self, name, img_array)
        Circle.__init__(self, radius=radius)
        self.angle = angle
        self.dist_from_center = dist_from_center

    @type_accept(phan_cent_point=Point)
    def set_center_via_phan_center(self, phan_cent_point):
        """Set the center of the ROI based on phantom center."""
        y_shift = -np.sin(np.deg2rad(self.angle))*self.dist_from_center
        x_shift = np.cos(np.deg2rad(self.angle))*self.dist_from_center
        self.center = Point(phan_cent_point.x+x_shift, phan_cent_point.y+y_shift)

    def get_roi_mask(self, outside=np.NaN):
        """Return a masked array of the ROI, with outside values being NaNs."""
        # create mask
        mask = sector_mask(self.img_array.shape, self.center, self.radius)
        # Apply mask to image
        masked_img = np.where(mask == True, self.img_array, outside)
        return masked_img


class HU_ROI(ROI_Disk):
    """An HU ROI object. Represents a circular area measuring either HU sample (Air, Poly, ...) or uniformity (bottom, left, ...)."""
    def __init__(self, name, angle, nominal_value, img_array=None, radius=None, dist_from_center=None, tolerance=None):
        super().__init__(name, img_array, angle, radius, dist_from_center)
        self.nominal_val = nominal_value
        self.tolerance = tolerance

    def get_pixel_value(self, mode='mean'):
        """Return the pixel value calculation within the ROI."""
        masked_img = self.get_roi_mask()
        if mode == 'mean':
            pix_val = np.nanmean(masked_img)
        if mode == 'median':
            pix_val = np.nanmedian(masked_img)
        return pix_val

    @property
    def value_diff(self):
        """The difference in HU between measured and nominal."""
        return self.pixel_value - self.nominal_val

    @property
    def pixel_value(self):
        """The mean pixel value of the ROI."""
        return self.get_pixel_value()

    @property
    def passed(self):
        """Boolean specifying if ROI pixel value was within tolerance of the nominal value."""
        return self.value_diff <= self.tolerance

    def get_pass_fail_color(self, passed='blue', failed='red'):
        if self.passed:
            return passed
        else:
            return failed

class Slice:
    """A base class for analyzing specific slices of a CBCT dicom set."""
    def __init__(self, images, slice_num, mode, phan_roll=None):
        """
        :param images: The CBCT image set as a 3D numpy array.
        :type images: ndarray
        :param slice_num: Slice number of interest
        :type slice: int
        :param mode: Mode of combining surrounding images to lower noise. Options: 'mean', 'max', 'median'
        :type mode: str
        """
        self.image = ImageObj(combine_surrounding_slices(images, slice_num, mode=mode))
        self.ROIs = OrderedDict()
        self.phan_roll = phan_roll

    def add_ROI(self, *ROIs):
        """Add ROIs to the slice.

        The ROI is added to a dictionary using the ROI.name attr as the key, with the
        ROI itself as the value.

        :type ROI: pylinac.cbct.ROI, subclass of ROI
        """
        for roi in ROIs:
            if roi.name in self.ROIs.keys():
                print("ROI name already instantiated. Use a different name")
                return
            self.ROIs[roi.name] = roi

    def find_phan_center(self, threshold):
        """Determine the location of the center of the phantom."""
        SOI_bw = self.image.convert2BW(threshold, return_it=True)  # convert slice to binary based on threshold
        SOI_bw = ndimage.binary_fill_holes(SOI_bw.pixel_array)  # fill in air pockets to make one solid ROI
        SOI_labeled, no_roi = ndimage.label(SOI_bw)  # identify the ROIs
        if no_roi < 1 or no_roi is None:
            raise ValueError("Unable to locate the CatPhan")
        hist, bins = np.histogram(SOI_labeled, bins=no_roi)  # hist will give the size of each label
        SOI_bw_clean = np.where(SOI_labeled == np.argmax(hist), 1, 0)  # remove all ROIs except the largest one (the CatPhan)
        center_pixel = ndimage.measurements.center_of_mass(SOI_bw_clean)  # returns (y,x)
        self.phan_center = Point(center_pixel[1], center_pixel[0])

        # Propagate the phantom center out to the ROIs (they initially don't have a center because it's relative
        # to the phantom center)
        for roi in self.ROIs.values():
            roi.set_center_via_phan_center(self.phan_center)


class Base_HU_Slice(Slice):
    """Base class for the HU and Uniformity Slices."""

    def get_ROI_vals(self):
        """Return the HU values of the HU ROIs."""
        return {key: val.pixel_value for key, val in self.ROIs.items()}

    def get_ROI_passing(self):
        """Return the pass/fails for the ROIs."""
        return {key: val.passed for key, val in self.ROIs.items()}

    @property
    def overall_passed(self):
        if all(self.get_ROI_passing().values()):
            return True
        else:
            return False

class HU_Slice(Base_HU_Slice):
    """Class for analysis of the HU slice of the CBCT dicom data set."""
    radius2objs = 120
    object_radius = 9
    tolerance = 40
    air_bubble_size = 450


    def __init__(self, images, slice_num, threshold, mode, phan_roll=None):
        super().__init__(images, slice_num, mode, phan_roll)

        air = HU_ROI('Air', 90, -1000)
        pmp = HU_ROI('PMP', 120, -200)
        ldpe = HU_ROI('LDPE', 180, -100)
        poly = HU_ROI('Poly', -120, -35)
        acrylic = HU_ROI('Acrylic', -60, 120)
        delrin = HU_ROI('Delrin', 0, 340)
        teflon = HU_ROI('Teflon', 60, 990)
        self.add_ROI(air, pmp, ldpe, poly, acrylic, delrin, teflon)

        for roi in self.ROIs.values():
            roi.img_array = self.image.pixel_array
            roi.radius = self.object_radius
            roi.tolerance = self.tolerance
            roi.dist_from_center = self.radius2objs

        super().find_phan_center(threshold)

    def determine_phantom_roll(self, threshold, scaling_ratio):
        """Determine the "roll" of the phantom based on the two air bubbles in the HU slice."""
        # convert slice to logical
        SOI = self.image.convert2BW(threshold, return_it=True)
        # invert the SOI; this makes the Air == 1 and phantom == 0
        SOI.invert_array()
        # determine labels and number of rois of inverted SOI
        labels, no_roi = ndimage.measurements.label(SOI.pixel_array)
        # calculate ROI sizes of each label TODO: simplify the air bubble-finding
        roi_sizes = [ndimage.measurements.sum(SOI.pixel_array, labels, index=item) for item in range(1, no_roi + 1)]
        # extract air bubble ROIs (based on size threshold)
        bubble_thresh = self.air_bubble_size * (scaling_ratio ** 2)
        air_bubbles = [idx + 1 for idx, item in enumerate(roi_sizes) if
                       item < bubble_thresh * 1.5 and item > bubble_thresh / 1.5]
        # if the algo has worked correctly, it has found 2 and only 2 ROIs (the air bubbles)
        if len(air_bubbles) == 2:
            air_bubble_CofM = ndimage.measurements.center_of_mass(SOI.pixel_array, labels, air_bubbles)
            y_dist = air_bubble_CofM[0][0] - air_bubble_CofM[1][0]
            x_dist = air_bubble_CofM[0][1] - air_bubble_CofM[1][1]
            angle = np.arctan2(y_dist, x_dist)
            if angle < 0:
                roll = abs(angle) - np.pi/2
            else:
                roll = angle - np.pi/2
            phan_roll = roll
        else:
            phan_roll = 0
            print("Warning: CBCT phantom roll unable to be determined; assuming 0")

        return phan_roll


class UNIF_Slice(Base_HU_Slice):
    """Class for analysis of the Uniformity slice of the CBCT dicom data set."""
    radius2objs = 110
    obj_radius = 20
    tolerance = 40

    def __init__(self, images, slice_num, threshold, mode, phan_roll):
        super().__init__(images, slice_num, mode, phan_roll)
        center = HU_ROI('Center', 0, 0, dist_from_center=0)  # angle is irrelevant for center ROI
        right = HU_ROI('Right', 0, 0)
        top = HU_ROI('Top', -90, 0)
        left = HU_ROI('Left', 180, 0)
        bottom = HU_ROI('Bottom', 90, 0)
        self.add_ROI(center, right, top, left, bottom)

        for roi in self.ROIs.values():
            roi.img_array = self.image.pixel_array
            roi.tolerance = self.tolerance
            roi.radius = self.obj_radius
            if roi.name != center.name:
                roi.dist_from_center = self.radius2objs

        super().find_phan_center(threshold)

class Locon_Slice(Slice):
    """Class for analysis of the low contrast slice of the CBCT dicom data set."""
    # TODO: work on this
    def __init__(self, images, slice_num, mode, phan_roll):
        super().__init__(images, slice_num, mode, phan_roll)
        raise NotImplementedError


class SR_Circle_ROI(CircleProfile, ROI):
    def __init__(self, name, img_array, radius):
        CircleProfile.__init__(self, radius=radius)
        ROI.__init__(self, name, img_array)

    @type_accept(phan_cent_point=Point)
    def set_center_via_phan_center(self, phan_cent_point):
        """For the SR ROIs, the phantom center is also the SR ROI center."""
        self.center = Point(phan_cent_point.x, phan_cent_point.y)

class SR_Slice(Slice):
    """Class for analysis of the Spatial Resolution slice of the CBCT dicom data set.

    This slice is quite different from the other CBCT slices. Rather than having ROIs like
    the HU and UNIF slices, this one calculates the resolution using several CircleProfiles.
    It computes 5 profiles, each one pixel smaller than the other, averages them, and then
    computes the spatial resolution from that.

    """
    LP_freq = (0.2, 0.4, 0.6, 0.8, 1, 1.2)
    radius2profs = np.arange(95, 100)

    def __init__(self, images, slice_num, threshold, mode, phan_roll):
        super().__init__(images, slice_num, mode, phan_roll)
        self.LP_MTF = {}
        for idx, radius in enumerate(self.radius2profs):
            c = SR_Circle_ROI(idx, self.image.pixel_array, radius=radius)
            self.add_ROI(c)

        super().find_phan_center(threshold)

    def calc_median_profile(self, roll_offset=0):
        # extract the profile for each ROI (5 adjacent profiles)
        for roi in self.ROIs.values():
            roi.get_profile(self.image.pixel_array, size=2*np.pi*1000, start=np.pi+roll_offset)
        # average profiles together
        prof = np.zeros(len(roi.y_values))
        for idx, roi in enumerate(self.ROIs.values()):
            prof += roi.y_values
        prof /= len(self.ROIs)
        # ground profile
        new_prof = Profile(prof)
        new_prof.filter_profile(6)
        new_prof.ground_profile()
        return new_prof

    def _find_LP_peaks(self, profile):
        """Find the peaks along the line pair profile extracted. Because of the varying width of lead/no lead, 3 searches are done
        with varying widths of peak spacing. This is to ensure that only 1 peak is found for the larger LPs, but does not pass over the
        thinner LPs further down the profile.

        :param profile: 1-D profile of the Line Pairs (normally from what is returned by return_LP_profile).
        :type profile: pylinac.core.profile.Profile
        :returns: values of peaks, indices of peaks. Will be 1-D numpy arrays.
        """

        region_1_bound = 4000  # approximate index between 1st and 2nd LP regions
        region_2_bound = 10500  # approximate index between 4th and 5th LP regions
        # region_3_bound = 17500  # approximate index between 8th and 9th LP regions; after this, regions become very hard to distinguish
        region_3_bound = 12300  # for head this can be 17500, but for thorax (low dose => low quality), we can only sample the first
        # 5 line pairs accurately
        max_vals_1, max_idx_1 = profile.find_peaks(min_peak_distance=150, max_num_peaks=2, exclude_rt_edge=0.9, return_it=True)
        max_vals_2, max_idx_2 = profile.find_peaks(min_peak_distance=48, exclude_lt_edge=0.1, exclude_rt_edge=0.7, return_it=True)
        max_vals_3, max_idx_3 = profile.find_peaks(min_peak_distance=25, exclude_lt_edge=0.3, exclude_rt_edge=0.65, return_it=True)
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
            min_val, min_idx = profile.find_valleys(exclude_lt_edge=max_idxs[idx], exclude_rt_edge=max_idxs[idx+1], max_num_peaks=1, \
                                                                                                                                 return_it=True)
            min_vals[idx] = min_val
            min_idxs[idx] = min_idx
        # now delete the valleys *in between* the LP regions
        min_vals = np.delete(min_vals, idx2del)
        min_idxs = np.delete(min_idxs, idx2del)
        return min_vals, min_idxs

    def _calc_MTF(self, max_vals, min_vals):
        """Calculate the Modulation Transfer Function of the Line-Pair profile.
        See http://en.wikipedia.org/wiki/Transfer_function#Optics for calculation.
        Maximum and minimum values are calculated by averaging the pixel
        values of the peaks/valleys found.
        """
        num_peaks = np.array((0, 2, 3, 3, 4, 4, 4)).cumsum()
        num_valleys = np.array((0, 1, 2, 2, 3, 3, 3)).cumsum()
        for key, LP_pair in zip(self.LP_freq, range(len(num_peaks) - 1)):
            region_max = max_vals[num_peaks[LP_pair]:num_peaks[LP_pair + 1]].max()
            region_min = min_vals[num_valleys[LP_pair]:num_valleys[LP_pair + 1]].min()
            self.LP_MTF[key] = (region_max - region_min) / (region_max + region_min)
        # normalize the values by the first LP
        max_mtf = np.array(list(self.LP_MTF.values())).max()
        for name, value in self.LP_MTF.items():
            self.LP_MTF[name] /= max_mtf

    def calc_MTF(self):
        """Calculate the line pairs of the SR slice."""
        profile = self.calc_median_profile(roll_offset=self.phan_roll)
        max_vals, max_idxs = self._find_LP_peaks(profile)
        min_vals, min_idxs = self._find_LP_valleys(profile, max_idxs)
        self._calc_MTF(max_vals, min_vals)

    @value_accept(percent=(60, 95))
    def get_MTF(self, percent=80):
        """Return the MTF value for the percent passed in.

        :param percent: The line-pair/mm value for the given MTF percentage.
        :type percent: int
        """
        # calculate x and y interpolations from Line Pair values and from the MTF measured
        x_vals_intrp = np.arange(self.LP_freq[0], self.LP_freq[-1], 0.01)
        x_vals = np.array(sorted(self.LP_MTF.keys()))
        y_vals = np.array(sorted(self.LP_MTF.values())[::-1])
        y_vals_intrp = np.interp(x_vals_intrp, x_vals, y_vals)

        mtf_percent = x_vals_intrp[np.argmin(np.abs(y_vals_intrp - (percent / 100)))]
        return mtf_percent


class GEO_ROI(ROI_Disk):
    """A circle ROI, much like the HU ROI, but with methods to find the center of the geometric "node"."""
    def __init__(self, name, img_array, angle, radius, dist_from_center):
        super().__init__(name, img_array, angle, radius, dist_from_center)
        self.node_CoM = None

    def _threshold_node(self):
        # create mask
        masked_img = self.get_roi_mask()
        # threshold image
        upper_band_pass = np.where(masked_img > np.nanmedian(masked_img) * 1.4, 1, 0)
        lower_band_pass = np.where(masked_img < np.nanmedian(masked_img) * 0.6, 1, 0)
        bw_node = upper_band_pass + lower_band_pass
        return bw_node

    def find_node_center(self):
        bw_node = self._threshold_node()
        # label ROIs found
        label, no_roi = ndimage.measurements.label(bw_node)
        hist, bins = np.histogram(label, bins=no_roi)  # hist will give the size of each label
        bw_node_cleaned = np.where(label == np.argsort(hist)[-2], 1, 0)  # remove all ROIs except the largest one
        label, no_roi = ndimage.measurements.label(bw_node_cleaned)
        if no_roi != 1:
            raise ValueError("Did not find the geometric node.")
        # determine the center of mass of the geometric node
        node_CoM = ndimage.measurements.center_of_mass(bw_node_cleaned, label)
        self.node_CoM = Point(node_CoM[1], node_CoM[0])  # the scipy com function returns (y, x), thus inversion


class GEO_Line(Line):
    """Represents a line connecting two nodes on the Geometry Slice."""
    def __init__(self, name, node1, node2):
        """
        :param name:
        :param node1:
        :type node1: GEO_ROI
        :param node2:
        :type node2: GEO_ROI
        :return:
        """
        super().__init__()
        self.name = name
        self.node1 = node1
        self.node2 = node2

    @property
    def point1(self):
        return self.node1.node_CoM

    @property
    def point2(self):
        return self.node2.node_CoM


class GEO_Slice(Slice):
    """Class for analysis of the Geometry slice of the CBCT dicom data set.

    The Geometry class is slightly more complex than the HU and Uniformity classes.
    Four ROIs are set, which correspond to the locations of the 1 air and 3 acrylic "nodes".
    Within these ROIs the center of the nodes must be found.

    Once the nodes centers are found four lines are constructed by linking the node centers,
    which should be 50mm apart.
    """
    roi_nominal_value = 50
    radius2objs = 72
    obj_radius = 20
    tolerance = 1

    def __init__(self, images, slice_num, threshold, mode, phan_roll, mm_per_pixel):
        super().__init__(images, slice_num, mode, phan_roll)
        self.mm_per_pixel = mm_per_pixel

        arr = self.image.pixel_array
        tl = GEO_ROI('Top-Left', arr, -135, self.obj_radius, self.radius2objs)
        tr = GEO_ROI('Top-Right', arr, -45, self.obj_radius, self.radius2objs)
        br = GEO_ROI('Bottom-Right', arr, 45, self.obj_radius, self.radius2objs)
        bl = GEO_ROI('Bottom-Left', arr, 135, self.obj_radius, self.radius2objs)
        self.add_ROI(tl, tr, br, bl)

        # Construct the Lines, mapping to the nodes they connect to
        lv = GEO_Line('Left-Vert', tl, bl)
        bh = GEO_Line('Bottom-Horiz', bl, br)
        rv = GEO_Line('Right-Vert', tr, br)
        th = GEO_Line('Top-Horiz', tl, tr)
        self.add_line(lv, bh, rv, th)

        super().find_phan_center(threshold)

    def add_line(self, *lines):
        """Add GEO_Lines; similar to add_ROI of Slice class."""
        if not hasattr(self, 'lines'):
            self.lines = {}
        for line in lines:
            self.lines[line.name] = line

    def calc_node_centers(self):
        """Calculate the center-of-mass of the geometric nodes."""
        for roi in self.ROIs.values():
            roi.find_node_center()

    def get_line_lengths(self):
        return {line_key: line.length*self.mm_per_pixel for line_key, line in self.lines.items()}

    @property
    def overall_passed(self):
        #TODO: add check against tolerance
        if all(self.get_line_lengths().values()):
            return True
        else:
            return False

class CBCT:
    """A class for loading and analyzing Cone-Beam CT DICOM files of a CatPhan 504 (Varian; Elekta 503 is being developed.
    Analyzes: Uniformity, Spatial Resolution, Image Scaling & HU Linearity.
    """
    BW_threshold = -800  # The threshold to apply to images to convert to B&W, which will highlight the CatPhan.
    cbct_demo_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'demo_files', 'cbct')
    demo_zip = osp.join(cbct_demo_dir, 'High quality head.zip')
    demo_folder = osp.join(cbct_demo_dir, 'High quality head')

    def __init__(self):
        """Image properties and algorithm attrs."""
        self.settings = {}  # dict of settings the algo will use to compute the HU, SR, etc.

    def load_demo_images(self):
        """Load the CBCT demo images."""

        # unpack zip folder
        shutil.unpack_archive(self.demo_zip, self.cbct_demo_dir)

        # walk through the extracted demo folder and pull the names of the CT dicom files
        for par_dir, sub_dir, folder in os.walk(self.demo_folder):
            file_list = [osp.join(par_dir, item) for item in folder if item.endswith('.dcm') and item.startswith('CT')]

        self._load_files(file_list)

        # delete the unpacked demo folder
        try:
            shutil.rmtree(self.demo_folder)
        except:
            print("Extracted demo images were not able to be deleted. You can manually delete them if you "
                  "like from %s" % self.demo_folder)

    def load_folder_UI(self):
        """Load the CT DICOM files from a folder using a UI."""
        folder = get_folder_UI()
        if folder:
            self.load_folder(folder)

    def load_folder(self, folder):
        """Load the CT DICOM files from a folder string input

        :param folder: Path to the folder
        :type folder: str
        """
        # check that folder is valid
        if not osp.isdir(folder):
            raise NotADirectoryError("Path given was not a Directory/Folder")

        for par_dir, sub_dir, files in os.walk(folder):
            filelist = [osp.join(par_dir, item) for item in files if item.endswith('.dcm') and item.startswith('CT')]
        self._load_files(filelist)

    def _load_files(self, file_list, im_size=512):
        """Load CT DICOM files given a list of image paths.

        :param file_list: list of file paths to load
        :type file_list: list
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

        self.dicom_data = dicom.read_file(item, stop_before_pixels=True)
        # determine settings needed for given CBCT.
        self._get_settings(dcm)
        # convert images from CT# to HU using dicom tags
        self._convert2HU()


    def _sort_images(self, im_order, images):
        """Sort the images according to the image order."""
        #TODO: convert this to some sort of in-place method
        sorted_images = np.zeros(images.shape, dtype=int)
        for idx, orig in enumerate(im_order):
            sorted_images[:, :, orig] = images[:, :, idx]
        self.images = sorted_images

    def _get_settings(self, dicom_file):
        """Based on the images input, determine the starting conditions for the analysis. There are various manufacturers, protocols,
        and field-of-views for CBCTs. Based on this, the starting conditions can be determined.

        :param dicom_file: dicom file as read by pydicom
        """
        dcm = dicom_file
        """determine algorithm conditions"""
        # determine scaling ratio; this is basically a magnification factor for the following thresholds and radii and so forth.
        self.scaling_ratio = 0.488 / dcm.PixelSpacing[0]
        self.mm_per_pixel = dcm.PixelSpacing[0]
        # Set slice locations based on manufacturer
        if dcm.Manufacturer in known_manufacturers:
            if dcm.Manufacturer == known_manufacturers[0]:
                self.HU_slice_num = 32
                self.UN_slice_num = 9
                self.SR_slice_num = 44
                self.LC_slice_num = 23
            else:
                raise NotImplementedError("Elekta not yet implemented")
        else:
            raise ValueError("Unknown Manufacturer")

    def _convert2HU(self):
        """Convert the images from CT# to HU."""
        self.images *= self.dicom_data.RescaleSlope
        self.images += self.dicom_data.RescaleIntercept

    @lazyproperty
    def phantom_roll(self):
        return self.calc_phantom_roll()

    def calc_phantom_roll(self):
        """Determine the roll of the phantom by calculating the angle of the two
        air bubbles on the HU slice. Delegation method."""
        HU = HU_Slice(self.images, self.HU_slice_num, self.BW_threshold, 'mean')
        return HU.determine_phantom_roll(self.BW_threshold, self.scaling_ratio)

    def construct_HU(self):
        """Determine HU ROI values from HU Slice."""
        self.HU = HU_Slice(self.images, self.HU_slice_num, self.BW_threshold, 'mean', self.phantom_roll)

    def construct_SR(self):
        """Determine the Spatial Resolution from the Line-Pair slice."""
        self.SR = SR_Slice(self.images, self.SR_slice_num, self.BW_threshold, 'max', self.phantom_roll)
        self.SR.calc_MTF()

    def construct_GEO(self):
        """Determine the geometric distortion of the images.
         The distortion is found by finding the 4 small points on the inner region of the HU slice
         and calculating their distances to each other. Nominal distance is 5cm apart.
        """
        self.GEO = GEO_Slice(self.images, self.HU_slice_num, self.BW_threshold, 'median', self.phantom_roll, self.mm_per_pixel)
        self.GEO.calc_node_centers()

    def construct_UNIF(self):
        """Determine the HU uniformity from Uniformity slice by analyzing 5 ROIs: center, top, right, bottom, left."""
        self.UN = UNIF_Slice(self.images, self.UN_slice_num, self.BW_threshold, 'mean', self.phantom_roll)

    def construct_Locon(self):
        """Determine Low Contrast values."""
        self.Locon = Locon_Slice(self.images, self.LC_slice_num, 'mean', self.phantom_roll)

    def plot_analyzed_image(self):
        """Draw the ROIs and lines the calculations were done on or based on."""
        # create figure
        fig, ((UN_ax, HU_ax), (SR_ax, LOCON_ax)) = plt.subplots(2,2)

        # Uniformity objects
        UN_ax.imshow(self.UN.image.pixel_array)
        for roi in self.UN.ROIs.values():
            color = roi.get_pass_fail_color()
            roi.add_to_axes(UN_ax, edgecolor=color)
        UN_ax.autoscale(tight=True)
        UN_ax.set_title('Uniformity Slice')

        # HU objects
        HU_ax.imshow(self.HU.image.pixel_array)
        for roi in self.HU.ROIs.values():
            color = roi.get_pass_fail_color()
            roi.add_to_axes(HU_ax, edgecolor=color)
        HU_ax.autoscale(tight=True)
        HU_ax.set_title('HU & Geometric Slice')

        # GEO objects
        for line in self.GEO.lines.values():
            line.add_to_axes(HU_ax, color='blue')

        # SR objects
        SR_ax.imshow(self.SR.image.pixel_array)
        for roi in [self.SR.ROIs[0], self.SR.ROIs[4]]:
            roi.add_to_axes(SR_ax, edgecolor='blue')
        SR_ax.autoscale(tight=True)
        SR_ax.set_title('Spatial Resolution Slice')

        # Locon objects
        # LOCON_ax.imshow(self.Locon.image.pixel_array)
        # LOCON_ax.set_title('Low Contrast (In Development)')

        # show it all
        plt.show()

    def return_results(self):
        """Return and print the results of the analysis as a string."""

        #TODO: make prettier
        print(' - CBCT QA Test - ')
        print('HU Regions: ', self.HU.get_ROI_vals())
        print('HU Passed?: ', self.HU.overall_passed)
        print('Uniformity: ', self.UN.get_ROI_vals())
        print('Uniformity Passed?: ', self.UN.overall_passed)
        print('MTF 80% (lp/mm): ', self.SR.get_MTF(80))
        print('Geometric distances: ', self.GEO.get_line_lengths())
        print('Geometry Passed?: ', self.GEO.overall_passed)

    def analyze(self):
        """Single-method full analysis of CBCT DICOM files."""
        # self.calc_phantom_roll()
        self.construct_HU()
        self.construct_UNIF()
        self.construct_GEO()
        self.construct_SR()
        # self.construct_Locon()

    def run_demo(self):
        """Run the CBCT demo using the high-quality head protocol."""
        cbct = CBCT()
        cbct.load_demo_images()
        cbct.analyze()
        cbct.return_results()
        cbct.plot_analyzed_image()


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