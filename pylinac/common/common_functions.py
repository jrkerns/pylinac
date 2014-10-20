from __future__ import division, print_function
import os
import os.path as osp
import inspect

from future.builtins import str
from future.builtins import object
from PySide import QtGui
import numpy as np
from numpy import sqrt
from PIL import Image
import dicom


"""Common Functions used in PyTG142"""

def open_TG142():
    """
    Open AAPM TG-142 in PDF
    """
    open_PDF_file('TG-142 Linac QA.pdf')

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
    """
    Custom function that is equivalent to Matlab's uigetfile command + return the image/matrix. Returns numpy array.

    numpy_image = OpenLoadImage(UIdir=None,UIcaption='',UIfilters='')
    """

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
    app = check_QApp_running()

    # get user-defined image file
    filename = str(QtGui.QFileDialog.getOpenFileName(caption=UIcaption))
    return filename


def get_filenames(UIdir=None, UIcaption='', UIfilters=''):
    """
    Custom function that is equivalent to Matlab's uigetfile command. Returns filename as a string.

    filenamestring = GetFile(UIdir=None,UIcaption='',UIfilters='')
    """
    # if a QApplication isn't running turn one on; necessary to have one running to use QFileDialog()
    app = check_QApp_running()

    # get user-defined files
    filenames = (QtGui.QFileDialog.getOpenFileNames(caption=UIcaption))
    filenames = [str(f) for f in filenames]  #convert the PyQt string list to a list of standard strings

    if filenames:
        return filenames

def get_folder(UIdir=None, UIcaption=''):
    """
    :return string of Folder the user selected
    """
    app = check_QApp_running()

    folderstring = str(QtGui.QFileDialog.getExistingDirectory())
    return folderstring

def open_PDF_file(pdf):
    """
    Open a PDF file from the Files folder in PyTG142
    """
    filepath = osp.join(osp.split(osp.split(osp.abspath(__file__))[0])[0], 'Files', pdf)
    os.startfile(filepath, 'open')  #TODO: see if this is most appropriate


def remove_edges(image, edge=15):
    """
    Removes the edge pixels on all sides of an image by the edge value amount
    """
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

    def get_penum(self, side='left', penum_point=50):
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
        li = self.get_penum('left', X)
        ri = self.get_penum('right', X)
        fwxm = np.abs(ri - li)
        return fwxm

    def get_FWXM_center(self, X=50):
        """
        returns the center point (index) of FWXM
        """
        fwxm = self.get_FWXM(X)
        li = self.get_penum('left', X)
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
            li = self.get_penum('left', lower_penum)
            ui = self.get_penum('left', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'right':
            li = self.get_penum('right', lower_penum)
            ui = self.get_penum('right', upper_penum)
            pen = np.abs(ui - li)
            return pen
        elif side == 'both':
            li = self.get_penum('left', lower_penum)
            ui = self.get_penum('left', upper_penum)
            lpen = np.abs(ui - li)
            li = self.get_penum('right', lower_penum)
            ui = self.get_penum('right', upper_penum)
            rpen = np.abs(ui - li)
            pen = np.mean([lpen, rpen])
            return pen
        else:
            raise NameError("getpenumwidth input parameter not acceptable")

def point_to_2point_line_dist(point, line_point):
    """
    Determine the minimum distance from a point to a line defined by two points.
    Based on Wikipedia article: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line, paragraph: Line defined by two points

    :param point: (y,x)
    :param line_point1: (y,x)
    :return:
    """
    x0 = point[1]
    y0 = point[0]
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
    """
    smallest/biggest distance of a point to a sequence of line segments. Used in Starshot module.
    """
    if minormax == 'min':
        return min([point_line_dist(p, seg) for seg in segs])
    elif minormax == 'max':
        return max([point_line_dist(p, seg) for seg in segs])

def point2edge_min(image, point):
    """
    Calculates minimum distance from user point to image edges
    point = (y,x)
    """
    rows, cols = size(image)
    disttoedge = np.zeros(4)
    disttoedge[0] = rows - point[0]
    disttoedge[1] = cols - point[1]
    disttoedge[2] = point[0]
    disttoedge[3] = point[1]
    return min(disttoedge)

def size(matrix):
    """
    Matlab equivalent of size; returns the size of the matrix in [rows, columns]
    """
    rows = np.size(matrix, 0)
    cols = np.size(matrix, 1)
    return rows, cols

def invert(matrix):
    """
    Return the imcomplement of the matrix/image. Equivalent to Matlab's imcomplement function.
    """
    newmatrix = -matrix + np.max(matrix) + np.min(matrix)
    return newmatrix

def dist_2points(point1, point2):
    """
    Find the distance from point1 to point2
    """
    #TODO: make this multi-dimensional
    dist = np.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)
    return dist

def check_QApp_running():
    """
    opens a QtGui Application if need be; necessary for using things like QFileDialog
    """
    if QtGui.QApplication.instance() is None:
        return QtGui.QApplication([])

def go_up_dirlevel(levels=0):
    """Go up directory levels from where the caller file is located.
        :param levels: int specifying how many levels to go up. 0 goes to the current directory.
    """
    calling_file = inspect.stack()[1][1]
    calling_dir = osp.dirname(calling_file)
    new_dir = calling_dir
    while levels != 0:
        old_dir = new_dir
        new_dir = osp.dirname(old_dir)
        levels -= 1
    return new_dir