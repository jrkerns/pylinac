"""Module for processing EPID images that have acquired Winston-Lutz type images."""
import math
import os
import os.path as osp

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.optimize import minimize

from pylinac.core.image import DicomImage, Image
from pylinac.core.profile import SingleProfile
from pylinac.core.geometry import Point, Line
from pylinac.core.utilities import is_close
from pylinac.cbct import get_filled_area_ratio, get_bounding_box


class WinstonLutz:

    def __init__(self, directory):
        # load directory
        self.images = ImageHandler(directory)

    @property
    def iso_size_3d(self):

        def distance(p, lines):
            """Calculate the maximum distance to any line from the given point."""
            return max(line.distance_to(Point(p[0], p[1], p[2])) for line in lines)

        # res = differential_evolution(distance, bounds=[(-2, 2), (-2, 2), (-2, 2)], args=(self.images.cax_projections,),
        #                              mutation=0.0001)
        bounds = [(-2, 2), (-2, 2), (-2, 2)]
        res = minimize(distance, np.array([0, 0, 0]), args=self.images.cax_projections, bounds=bounds)

        self.wobble.radius = res.fun
        self.wobble.center = Point(res.x[0], res.x[1])

    @property
    def iso_bb_distance(self):
        pass

    def plot_gantry_sag(self):
        gantry_imgs = [image for image in self.images if image.variable_axis in ('Gantry', 'Reference')]
        gantry_angles = [image.gantry_angle for image in gantry_imgs]
        gantry_sag = [image.z_shift for image in gantry_imgs]
        plt.plot(gantry_angles, gantry_sag, 'bo')
        plt.grid('on')
        plt.show()


class ImageHandler(list):

    def __init__(self, directory):
        super().__init__()
        for basefile in os.listdir(directory):
            file = osp.join(directory, basefile)
            if file.endswith('.dcm'):
                image = WLImage(file)
                self.append(image)

    @property
    def cax_projections(self):
        return [proj.cax_line_projection for proj in self]


class WLImage(DicomImage):
    """Holds individual Winston-Lutz images and performs various processing for analysis."""

    def __init__(self, file):
        super().__init__(file)
        self.filename = osp.basename(file)
        self.check_inversion()
        self.field_cax = self._find_field_centroid()
        self.bb = self._find_bb()

    def __repr__(self):
        return "WLImage(G={}, B={}, P={})".format(self.gantry_angle, self.collimator_angle, self.couch_angle)

    def _find_field_centroid(self):
        threshold_img = Image.load(self.array)
        min, max = np.percentile(self.array, [5, 99.5])
        threshold_img.threshold((max - min)/2 + min)
        coords = ndimage.measurements.center_of_mass(threshold_img)
        p = Point(x=coords[-1], y=coords[0])
        return p

    def _find_bb(self):

        def is_boxlike(array):
            ymin, ymax, xmin, xmax = get_bounding_box(array)
            y = abs(ymax - ymin)
            x = abs(xmax - xmin)
            if x > y * 1.05 or x < y * 0.95:
                return False
            return True

        hmin, hmax = np.percentile(self.array, [5, 100])
        spread = hmax - hmin
        max_thresh = hmax
        found = False
        while not found:
            try:
                lower_thresh = hmax - spread / 2
                t = np.where((max_thresh > self) & (self >= lower_thresh), 1, 0)
                labeled_arr, num_roi = ndimage.measurements.label(t)
                roi_sizes, bin_edges = np.histogram(labeled_arr, bins=num_roi + 1)
                roi_sizes = [roi for roi in roi_sizes if roi > 20]
                bw_node_cleaned = np.where(labeled_arr == np.argsort(roi_sizes)[-3], 1, 0)
                expected_fill_ratio = np.pi / 4
                actual_fill_ratio = get_filled_area_ratio(bw_node_cleaned)
                if (expected_fill_ratio * 1.05 < actual_fill_ratio) or (actual_fill_ratio < expected_fill_ratio * 0.95):
                    raise ValueError
                if not is_boxlike(bw_node_cleaned):
                    raise ValueError
            except (IndexError, ValueError):
                max_thresh -= 0.03 * spread
                if max_thresh < hmin:
                    raise ValueError("Unable to locate the BB")
            else:
                found = True

        # determine the center of mass of the geometric node
        inv_img = Image.load(self.array)
        inv_img.invert()
        x_arr = np.abs(np.average(bw_node_cleaned, weights=inv_img, axis=0))
        x_com = SingleProfile(x_arr).fwxm_center(interpolate=True)
        y_arr = np.abs(np.average(bw_node_cleaned, weights=inv_img, axis=1))
        y_com = SingleProfile(y_arr).fwxm_center(interpolate=True)
        return Point(x_com, y_com)

    @property
    def gantry_angle(self):
        return self.dicom_dataset.GantryAngle

    @property
    def collimator_angle(self):
        return self.dicom_dataset.BeamLimitingDeviceAngle

    @property
    def couch_angle(self):
        return self.dicom_dataset.PatientSupportAngle

    @property
    def y_shift(self):
        return math.sin(math.radians(self.gantry_angle)) * self.cax_dist2bb.x

    @property
    def x_shift(self):
        return math.cos(math.radians(self.gantry_angle)) * self.cax_dist2bb.x

    @property
    def z_shift(self):
        if is_close(self.couch_angle, [0, 360], delta=2):
            return self.cax_dist2bb.y

    @property
    def cax_line_projection(self):
        p1 = Point()
        p2 = Point()
        UPPER_QUADRANT = self.gantry_angle <= 45 or self.gantry_angle >= 315 or 225 >= self.gantry_angle > 135
        LR_QUADRANT = 45 < self.gantry_angle <= 135 or 225 < self.gantry_angle <= 315
        if UPPER_QUADRANT:
            p1.y = 2
            p2.y = -2
            p1.z = self.z_shift
            p2.z = self.z_shift
            p1.x = 2 * math.tan(math.radians(self.gantry_angle)) + self.x_shift / math.cos(math.radians(self.gantry_angle))
            p2.x = 2 * math.tan(math.radians(self.gantry_angle + 180)) + self.x_shift / math.cos(math.radians(self.gantry_angle))
        elif LR_QUADRANT:
            p1.x = 2
            p2.x = -2
            p1.z = self.z_shift
            p2.z = self.z_shift
            p1.y = 2 / math.tan(math.radians(self.gantry_angle)) + self.y_shift / math.cos(math.radians(self.gantry_angle - 90))
            p2.y = 2 / math.tan(math.radians(self.gantry_angle + 180)) + self.y_shift / math.cos(math.radians(self.gantry_angle - 90))
        l = Line(p1, p2)
        return l

    @property
    def cax_dist2bb(self):
        dist = self.field_cax - self.bb
        return dist / self.dpmm

    @property
    def cax_dist2iso(self):
        pass

    @property
    def epid_line_projection(self):
        pass

    @property
    def epid_dist2bb(self):
        pass

    @property
    def epid_dist2iso(self):
        pass

    def plot(self):
        super().plot()
        plt.plot(self.bb.x, self.bb.y, 'r+')
        plt.plot(self.field_cax.x, self.field_cax.y, 'mx')

    @property
    def variable_axis(self):
        """The axis that is varying; i.e. the axis that isn't at 0."""
        G0 = is_close(self.gantry_angle, [0, 360])
        B0 = is_close(self.collimator_angle, [0, 360])
        P0 = is_close(self.couch_angle, [0, 360])
        if G0 and B0 and not P0:
            return 'Couch'
        elif G0 and P0 and not B0:
            return 'Collimator'
        elif P0 and B0 and not G0:
            return 'Gantry'
        elif P0 and B0 and G0:
            return 'Reference'
        else:
            return 'Combo'



if __name__ == '__main__':
    directory = r'C:\Users\JRKerns\Dropbox\Programming\Python\Projects\unorganized linac data\Winston-Lutz\Katy TB\with crosshairs\\Only gantry'
    wl = WinstonLutz(directory)
    wl.plot_gantry_sag()
    ttt = 1