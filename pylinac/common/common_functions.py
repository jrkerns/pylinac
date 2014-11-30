from __future__ import division, print_function, unicode_literals, absolute_import

from builtins import str

from future import standard_library

standard_library.install_aliases()

import os
import os.path as osp
import inspect
try:
    # Python 3.x
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
except ImportError:
    # Python 2.x
    from Tkinter import Tk
    from tkFileDialog import askopenfilename, askopenfilenames, askdirectory

import numpy as np
from numpy import sqrt
from PIL import Image
import dicom

"""Common Functions used in Pylinac"""

def load_image(filestring, to_gray=True):
    """
    Load image using PIL or pydicom via the filepath/string. Also attempt to get some image properties (DPI, SID, etc).
    :param to_gray: boolean specifying whether to return a simple numpy array of image or image object (with properties)

    """
    # create empty dictionary to put image info into
    improps = {'DPI':None, 'SID':None, 'Image Type':''}

    # read image depending on file type
    if not filestring.endswith('dcm'):
        imfile = Image.open(filestring)
        image = np.array(imfile)
        try:
            dpi = imfile.info['dpi']
            if len(dpi) > 1:
                # ensure all values are the same, i.e. x-resolution is the same as y-resolution
                if dpi[0] != dpi[1]:
                    raise ValueError("Image DPI is not equal in both directions")
                dpi = dpi[0]
        except:
            #DPI unable to be determined
            pass
        finally:
            improps['DPI'] = dpi
    else:  # if dicom file, use pydicom to import
        imfile = dicom.read_file(filestring)
        image = imfile.pixel_array
        improps['Image Type'] = 'DICOM'
        improps['SID'] = float(imfile.RTImageSID)
        try:
            pixel_spacing = float(imfile.ImagePlanePixelSpacing[0])
            improps['DPI'] = pixel_spacing * improps['SID']/1000
        except:
            pass

    return image, improps

def load_image_UI(UIdir=None, UIcaption='', UIfilters='', togray=True,
                  multiselect=False):  #TODO: update args to fit getopenfilename's args
    """Load an image using a UI Dialog."""

    filestring = get_filename(UIdir=UIdir, UIcaption=UIcaption, UIfilters=UIfilters, multiselect=multiselect)

    if filestring:
        # convert the filename (natively a QString) to a Python string
        filestring = str(filestring)
        image, improps = load_image(filestring, togray)
        return image, improps
    else:  # user hit cancel...
        return None, None

def get_filename(UIdir=None, UIcaption='', UIfilters='', multiselect=False):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    # if a QApplication isn't running turn one on; necessary to have one running to use QFileDialog()
    withdraw_tkinter()

    # get user-defined image file
    filename = askopenfilename(caption=UIcaption)
    return filename


def get_filenames(UIdir=None, UIcaption='', UIfilters=''):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    withdraw_tkinter()
    # get user-defined files
    filenames = askopenfilenames(caption=UIcaption)
    filenames = [str(f) for f in filenames]  #convert the PyQt string list to a list of standard strings

    if filenames:
        return filenames

def get_folder(UIdir=None, UIcaption=''):
    """
    :returns: string of Folder the user selected
    """
    withdraw_tkinter()
    folderstring = askdirectory()
    return folderstring

def open_PDF_file(pdf):
    """
    Open a PDF file from the Files folder in PyTG142
    """
    filepath = osp.join(osp.split(osp.split(osp.abspath(__file__))[0])[0], 'Files', pdf)
    os.startfile(filepath, 'open')  #TODO: see if this is most appropriate


def remove_edges(image, edge=15):
    """Removes the edge pixels on all sides of an image by the edge value amount."""
    im_out = image[edge - 1:-edge, edge - 1:-edge]
    return im_out


class Prof_Penum(object):
    def __init__(self, ydata, xdata=None, normalize_sides=True):
        #TODO: add error checking for array size, etc
        #        self.normsides = normsides
        ydata = ydata.astype(int)
        self.ydata = ydata
        ymaxind = np.argmax(ydata)
        #        self.ymax = np.max(ydata)
        if normalize_sides:
            self.ydata_left = ydata - np.min(ydata[0:ymaxind])
            self.ydata_right = ydata - np.min(ydata[ymaxind:-1])
            self.ymax_left = np.max(self.ydata_left)
            self.ymax_right = np.max(self.ydata_right)
        else:
            self.ydata_left = ydata
            self.ydata_right = ydata
            self.ymax_left = np.max(self.ydata_left)
            self.ymax_right = np.max(self.ydata_right)

        if xdata is None:
            self.xdata = np.arange(len(ydata))
        else:
            self.xdata = xdata

    def get_X_penum_idx(self, side='left', penum_point=50):
        """
        the main method of ProfPenum

        returns: the index of the point and the value of the point

        side indicates which side the method is looking from: left or right

        penumpoint indicates what value percentage is being sought. E.g. if
        penumpoint=20 the method looks for the index & value of the point along
        the 1D array that closest matches the 20% of maximum point on the given
        side
        """
        found = False
        if side == 'left':
            i, thresh = 0, self.ymax_left * penum_point / 100
            while found == False:
                if self.ydata_left[i] > thresh:
                    found = True
                    i -= 1
                elif i == len(self.ydata_left):
                    return None
                i += 1
        elif side == 'right':
            i, thresh = len(self.ydata) - 1, self.ymax_right * penum_point / 100
            while found == False:
                if self.ydata_right[i] > thresh:
                    found = True
                    i += 1
                elif i == 0:
                    return None
                i -= 1
        else:
            raise TypeError("side was not correctly specified; use 'left' or 'right'")

        return self.xdata[i]  #, self.ydata[i]

    def get_FWXM(self, X=50):
        """
        get the Full-Width X-Max, where X is the percentage height.
        E.g. X = 50 is 50% height, a.k.a half max, thus X=50 is the FWHM
        returns the width in number of elements of the FWXM
        """
        li = self.get_X_penum_idx('left', X)
        ri = self.get_X_penum_idx('right', X)
        fwxm = np.abs(ri - li)
        return fwxm

    def get_FWXM_center(self, X=50):
        """
        returns the center point (index) of FWXM
        """
        fwxm = self.get_FWXM(X)
        li = self.get_X_penum_idx('left', X)
        fwxmcen = np.abs(li + fwxm / 2)
        return fwxmcen

    def get_penum_width(self, side='left', lower_penum=20, upper_penum=80):
        """
        return the actual penumbral width of the profile. This is the
        standard "penumbra width" that med. phys. talks about in
        radiation profiles. Standard is the 80/20 width, although 90/10
        is sometimes used.

        side options include: left, right, & both (average)
        """
        if side == 'left':
            li = self.get_X_penum_idx('left', lower_penum)
            ui = self.get_X_penum_idx('left', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'right':
            li = self.get_X_penum_idx('right', lower_penum)
            ui = self.get_X_penum_idx('right', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'both':
            li = self.get_X_penum_idx('left', lower_penum)
            ui = self.get_X_penum_idx('left', upper_penum)
            lpen = np.abs(ui - li)
            li = self.get_X_penum_idx('right', lower_penum)
            ui = self.get_X_penum_idx('right', upper_penum)
            rpen = np.abs(ui - li)
            pen = np.mean([lpen, rpen])
            return pen
        else:
            raise NameError("getpenumwidth input parameter not acceptable")

    def get_field_value(self, field_width_percent=80, value='mean'):
        """Get the value of the field in the profile, either mean, median, or max, within the given field width.

        :param field_width_percent: The percent width relative to the FWHM to sample.
        :type field_width_percent: int < 100
        :param value: Value type to extract. Either 'mean', 'median', or 'max' of field values.
        :type value: str
        """

        fwhmc = self.get_FWXM_center()
        fwhm = self.get_FWXM() * 0.8
        left = fwhmc - fwhm/2
        right = fwhmc + fwhm/2

        field_values = self.ydata[left:right]

        if value == 'mean':
            return field_values.mean()
        elif value == 'median':
            return np.median(field_values)
        elif value == 'max':
            return field_values.max()

def point_to_2point_line_dist(point, line_point):
    """
    Determine the minimum distance from a point to a line defined by two points.
    Based on Wikipedia article: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line, paragraph: Line defined by two points

    :param point: (y,x)
    :type point: Point
    :param line_point1: (y,x)
    :return:
    """
    x0 = point.x
    y0 = point.y
    x1 = line_point[0,1]
    y1 = line_point[0,0]
    x2 = line_point[1,1]
    y2 = line_point[1,0]

    Dx = x2 - x1
    Dy = y2 - y1
    numerator = np.abs(Dy*x0 - Dx*y0 - x1*y2 + x2*y1)
    denom = np.sqrt(Dx**2 + Dy**2)
    distance = numerator/denom
    return distance

def point_line_dist(p, seg, testSegmentEnds=False):
    """
    Minimum Distance between a Point and a Line
    Written by Paul Bourke,    October 1988
    http://astronomy.swin.edu.au/~pbourke/geometry/pointline/

    input:
    p: point, y,x
    seg: y1,x1,y2,x2

    """

    y3, x3 = p
    y1, x1, y2, x2 = seg[0, 0], seg[0, 1], seg[1, 0], seg[1, 1]

    dx21 = (x2 - x1)
    dy21 = (y2 - y1)

    lensq21 = dx21 * dx21 + dy21 * dy21
    if lensq21 == 0:
        #20080821 raise ValueError, "zero length line segment"
        dy = y3 - y1
        dx = x3 - x1
        return sqrt(dx * dx + dy * dy)  # return point to point distance

    u = (x3 - x1) * dx21 + (y3 - y1) * dy21
    u = u / float(lensq21)

    x = x1 + u * dx21
    y = y1 + u * dy21

    if testSegmentEnds:
        if u < 0:
            x, y = x1, y1
        elif u > 1:
            x, y = x2, y2

    dx30 = x3 - x
    dy30 = y3 - y

    return sqrt(dx30 * dx30 + dy30 * dy30)

def point_line_dist_multiline(p, segs, minormax='max'):
    """The smallest/biggest distance of a point to a sequence of line segments. Used in Starshot module."""
    if minormax == 'min':
        return min([point_line_dist(p, seg) for seg in segs])
    elif minormax == 'max':
        return max([point_line_dist(p, seg) for seg in segs])

def point2edge_min(image, point):
    """Calculates minimum distance from user point to image edges

    :param image: numpy image array
    :type image: numpy.ndarray
    :param point: The point to calculate from
    :type point: Point

    """
    rows = np.size(image, 0)
    cols = np.size(image, 1)
    disttoedge = np.zeros(4)
    disttoedge[0] = rows - point.y
    disttoedge[1] = cols - point.x
    disttoedge[2] = point.y
    disttoedge[3] = point.x
    return min(disttoedge)

def invert(matrix):
    """Return the imcomplement of the matrix/image. Equivalent to Matlab's imcomplement function."""
    newmatrix = -matrix + np.max(matrix) + np.min(matrix)
    return newmatrix

def dist_2points(point1, point2):
    """
    Find the distance from point1 to point2
    """
    #TODO: make this multi-dimensional
    dist = np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)
    return dist

def withdraw_tkinter():
    """Opens and withdraws a Tk window. Necessary so a base window doesn't open."""
    Tk.withdraw()

def go_up_dirlevel(levels=0):
    """Go up directory levels from where the caller file is located.

    :param levels: Specifies how many levels to go up. 0 goes to the current directory.
    :type levels: int
    """
    calling_file = inspect.stack()[1][1]
    calling_dir = osp.dirname(calling_file)
    new_dir = calling_dir
    while levels != 0:
        old_dir = new_dir
        new_dir = osp.dirname(old_dir)
        levels -= 1
    return new_dir

def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = list(range(len(y_axis)))

    if len(y_axis) != len(x_axis):
        raise ValueError

    # needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis

def peak_detect(y, x=None, threshold=0, min_peak_width=10, delta=0, max_num_peaks=None, find_min_instead=False):
    """Find the peaks or valleys of a 1-D signal. Uses the difference (np.diff) in signal to find peaks. Current limitations
    include:
    1) Only for use in 1-D data; 2-D may be possible with the gradient function. 2) Will not detect peaks at the very edge of array
    (i.e. 0 or -1 index)

    :param y: 1-D y-data of signal.
    :type y: numpy array
    :param x: 1-D x-data of signal. If left as None, will create a uniform range the length of the y-data.
    :type x: numpy array
    :param threshold: The value the peak must be above to be considered a peak.
        This removes "peaks" that are in a low-value region. If passed an int, the actual value is the threshold;
        if a float <1.0 is passed, it will threshold that percent. E.g. when passed 15, any peak less than 15 is
        removed. When passed 0.4, any peak less than 40% of the maximum value will be removed.
    :type threshold: int, float
    :param min_peak_width: The number of elements apart a peak must be from neighboring peaks.
    :type min_peak_width: int
    :param delta: The value a peak candidate must be higher than its neighbors by to be considered a true peak. Not yet implemented.
    :type delta: int
    :param max_num_peaks: Specify up to how many peaks will be returned. E.g. if 3 is passed in and 5 peaks are found, only the 3 largest
        peaks will be returned.
    :param find_min_instead: If True, algorithm will find minimums of y instead of maximums.
    :type find_min_instead: bool
    :returns: two 1-D numpy arrays: max_vals, max_idxs; max_vals contains the y-values of the peaks, max_idxs contains the x-index of the
        peaks.
    """
    peak_vals = []  # a list to hold the y-values of the peaks. Will be converted to a numpy array
    peak_idxs = []  # ditto for x-values (index) of y data.

    x, y = _datacheck_peakdetect(x, y)  # check length and convert to numpy arrays

    if find_min_instead:
        y = -y

    y_diff = np.diff(y.astype(float))  # Had problems with uint input. y_diff *must* be converted to signed type.

    if isinstance(threshold, float):
        if threshold >= 1:
            raise ValueError("When threshold is passed a float, value must be less than 1")
        else:
            data_range = y.max() - y.min()
            threshold = threshold * data_range + y.min()

    """Find all potential peaks"""
    for idx in range(len(y_diff)-1):
        # For each item of the diff array, check if:
        # 1) The y-value is above the threshold.
        # 2) The value of y_diff is positive (negative for valley search), it means the y-value changed upward.
        # 3) The next y_diff value is zero or negative (or positive for valley search); a positive-then-negative diff value means the value
        # is a peak of some kind. If the diff is zero it could be a flat peak, which still counts.

        # 1)
        if y[idx + 1] < threshold:
            continue

        y1_gradient = y_diff[idx] > 0
        y2_gradient = y_diff[idx + 1] <= 0

        # 2) & 3)
        if y1_gradient and y2_gradient:
            # If the next value isn't zero it's a single-pixel peak. Easy enough.
            if y_diff[idx+1] != 0:
                peak_vals.append(y[idx + 1])
                peak_idxs.append(x[idx + 1])
            # Else if the diff value is zero, it could be a flat peak, or it could keep going up; we don't know yet.
            else:
                # Continue on until we find the next nonzero diff value.
                shift = 1
                while y_diff[(idx + 1) + shift] == 0:
                    shift += 1
                # If the next diff is negative (or positive for min), we've found a peak. Also put the peak at the center of the flat
                # region.
                is_a_peak = y_diff[(idx + 1) + shift] < 0
                if is_a_peak:
                    peak_vals.append(y[(idx + 1) + np.round(shift / 2)])
                    peak_idxs.append(x[(idx + 1) + np.round(shift / 2)])

    # convert to numpy arrays
    peak_vals = np.array(peak_vals)
    peak_idxs = np.array(peak_idxs)

    """Enforce the min_peak_distance by removing smaller peaks."""
    index = 0
    while index < len(peak_idxs) - 1:
        # For each peak, determine if the next peak is within the look_ahead range.

        # If the second peak is closer than min_peak_distance to the first peak, find the larger peak and remove the other one.
        if peak_idxs[index] > peak_idxs[index+1] - min_peak_width:
            if peak_vals[index] > peak_vals[index+1]:
                idx2del = index + 1
            else:
                idx2del = index
            peak_vals = np.delete(peak_vals, idx2del)
            peak_idxs = np.delete(peak_idxs, idx2del)
        else:
            index += 1

    """If Maximum Number passed, return only up to number given based on a sort of peak values."""
    if max_num_peaks is not None and len(peak_idxs) > max_num_peaks:
        sorted_peak_vals = peak_vals.argsort()  # sorts low to high
        peak_vals = peak_vals[sorted_peak_vals[-max_num_peaks:]]
        peak_idxs = peak_idxs[sorted_peak_vals[-max_num_peaks:]]

    # If we were looking for minimums, convert the values back to the original sign
    if find_min_instead:
        peak_vals = -peak_vals

    return peak_vals, peak_idxs


def sector_mask(shape, centre, radius, angle_range):
    """
    Return a boolean mask for a circular sector. The start/stop angles in
    `angle_range` should be given in clockwise order.

    Found here: https://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array/18354475#18354475
    """

    x, y = np.ogrid[:shape[0], :shape[1]]
    cy, cx = centre
    # tmin, tmax = np.deg2rad(angle_range)
    tmin, tmax = angle_range

    # ensure stop angle > start angle
    if tmax < tmin:
        tmax += 2 * np.pi

    # convert cartesian --> polar coordinates
    r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
    theta = np.arctan2(x - cx, y - cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2 * np.pi)

    # circular mask
    circmask = r2 <= radius * radius

    # angular mask
    anglemask = theta <= (tmax - tmin)

    return circmask * anglemask