# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 10:33:57 2013

PyQA is a GUI that utilizes the tool set of pylinac. It is meant both as a demo of pylinac as well as a standalone tool for physicists
who need the power of software to do good physics quality assurance.

@author: James
"""

from __future__ import print_function, division, absolute_import
from functools import partial

from PySide import QtGui, QtCore
from PySide.QtCore import QObject as QO
import matplotlib.cm as cm

# Would do relative imports (..) but when running as __main__ then imports break down. See PEP 366:
# http://legacy.python.org/dev/peps/pep-0366/
from pyqa_gui.PyQAui import Ui_PyQA_MainWindow
from pylinac.starshot.starshot import Starshot
from pylinac.vmatqa.vmat import VMAT

DEBUG = True
analysis_types = ('star', 'vmat')
im_types = ('vmat_open', 'vmat_dmlc')
sub_analysis_types = ('drgs', 'drmlc')

class PyQA(Ui_PyQA_MainWindow):
    def __init__(self, MainWindow):
        # call setupUi from parent class
        self.setupUi(MainWindow)

        """ Initialize pylinac classes"""
        # Ideally, these classes would inherit the pyqa main window, but because the attrs do not
        # actually construct until setupUi is called (as opposed to __init__) then the attrs cannot be inherited, but must be passed.
        self.star = Starshot()
        self.vmat = VMAT()

        # Turn buttons/widgets on/off
        self.qual_mplw_left.setVisible(False)
        self.qual_mplw_main.setVisible(False)
        self.qual_mplw_right.setVisible(False)
        self.quan_mplw.setVisible(False)
        self.View_Data_mplw_left.setVisible(False)
        self.View_Data_mplw_right.setVisible(False)

        self.settings_window = None
        self._readyforpoint = False
        # self.SID_image = None

        """connect slots & signals"""
        """Load buttons"""
        # Load starshot Image button
        QO.connect(self.star_load_img, QtCore.SIGNAL("clicked()"), partial(self.load_image_button, analysis_types[0]))

        """Analyze button"""
        QO.connect(self.analyze_button, QtCore.SIGNAL("clicked()"), self.analyze)

        """Debug Items"""
        # Turn off the debugging buttons if in Debug mode. Else connect signals and slots of debug buttons.
        if not DEBUG:
            self.autoloadStar.setVisible(False)
            self.autoloadvmat_drgs.setVisible(False)
            self.autoloadvmat_drmlc.setVisible(False)
        else:
            QO.connect(self.autoloadStar, QtCore.SIGNAL("clicked()"), partial(self.load_image_button, analysis_types[0], demo=True))
            QO.connect(self.autoloadvmat_drgs, QtCore.SIGNAL("clicked()"), partial(self.load_image_button, analysis_types[1],
                                                                                   sub_analysis_types[0], demo=True))
            QO.connect(self.autoloadvmat_drmlc, QtCore.SIGNAL("clicked()"), partial(self.load_image_button, analysis_types[1],
                                                                                    sub_analysis_types[1], demo=True))

    def load_image_button(self, analysis_type, sub_analysis_type=None, im_type=None, demo=False):
        """Have user get an image file using standard Open File-type dialog and plot it.
        :param type: string specifying type of image to load if a demo is desired.
        """
        if analysis_type == analysis_types[0]:  # starshot
            # set analysis type to self
            self.analysis_type = analysis_types[0]
            # either load demo or user's file
            if demo:
                Starshot.load_demo_image(self.star)
            else:
                Starshot.load_image_UI(self.star)
            # Name data
            img_data, improps = self.star.image, self.star.im_props
            # Set up MPLW views
            self.View_Data_mplw_right.setVisible(False)
            self.View_Data_mplw_left.setVisible(True)
            # Draw image to Data View plot
            self.View_Data_mplw_left.axes.imshow(img_data, cmap=cm.Greys_r)
            # Switch to starshot tab in settings
            self.Settings_tabWidget.setCurrentIndex(0)

        elif analysis_type == analysis_types[1]:  # VMAT
            self.analysis_type = analysis_types[1]
            if demo:
                if sub_analysis_type == sub_analysis_types[0]:  # DRGS
                    VMAT.load_demo_image(self.vmat, 'drgs')
                else:
                    VMAT.load_demo_image(self.vmat, 'drmlc')
            else:
                VMAT.load_image_UI(self.vmat)
            # Set up MPLW views
            self.View_Data_mplw_right.setVisible(True)
            self.View_Data_mplw_left.setVisible(True)
            # Draw image to Data View plot
            self.View_Data_mplw_left.axes.imshow(self.vmat.image_open, cmap=cm.Greys_r)
            self.View_Data_mplw_right.axes.imshow(self.vmat.image_mlc, cmap=cm.Greys_r)
            # Switch to VMAT tab in settings
            self.Settings_tabWidget.setCurrentIndex(1)

        #TODO: elif's for vmat, cbct, etc, when built
        else:
            raise NameError("image type not properly specified")

        # self.data = img_data
        # self.improps = improps

        # Switch to Data View tab
        self.main_widget.setCurrentIndex(1)

    def analyze(self):
        """Analyze image/data."""

        if self.analysis_type == analysis_types[0]:  # starshot
            # grab settings
            allow_inv = self.star_img_inversion.isChecked()
            radius = self.star_radius.value()
            min_peak_height = self.star_min_height.value()

            # analyze
            self.progressBar.setValue(25)
            self.statusbar.showMessage("Analyzing...")
            Starshot.analyze(self.star, allow_inv, radius, min_peak_height)
            self.progressBar.setValue(60)
            self.statusbar.showMessage("Plotting Results...")

            # draw analyzed images to PyQA
            Starshot.plot_analyzed_image(self.star, self.qual_mplw_left)
            Starshot.plot_analyzed_image(self.star, self.qual_mplw_right)
            Starshot.plot_analyzed_image(self.star, self.quan_mplw)
            # Tighten zoom on right qualitative widget
            self.qual_mplw_right.axes.set_xlim([self.star._wobble_center[1] - 10 * self.star._wobble_radius_pix,
                                                self.star._wobble_center[1] + 10 * self.star._wobble_radius_pix])
            self.qual_mplw_right.axes.set_ylim([self.star._wobble_center[0] - 10 * self.star._wobble_radius_pix,
                                                self.star._wobble_center[0] + 10 * self.star._wobble_radius_pix])
            # Tighten zoom on quantitative widget
            self.quan_mplw.axes.set_xlim([self.star._wobble_center[1] - 1.2 * self.star._profile_radius,
                                          self.star._wobble_center[1] + 1.2 * self.star._profile_radius])
            self.quan_mplw.axes.set_ylim([self.star._wobble_center[0] - 1.2 * self.star._profile_radius,
                                          self.star._wobble_center[0] + 1.2 * self.star._profile_radius])
            self.statusbar.showMessage("Found! See Analysis tabs", 3000)
            self.progressBar.setValue(0)

            # Post results
            self.results_browser.setText(self.star.return_string_results())

        elif self.analysis_type == analysis_types[1]:  # VMAT
            # grab settings
            test_idx = self.vmat_test_box.currentIndex()
            if test_idx == 0:  # DRGS
                test = 'drgs'
            else:
                test = 'drmlc'
            tolerance = self.vmat_tolerance_box.value()

            # analyze
            self.progressBar.setValue(25)
            self.statusbar.showMessage("Analyzing...")
            self.vmat.analyze(test,tolerance)
            self.progressBar.setValue(60)
            self.statusbar.showMessage("Plotting Results...")

            # draw analyzed images to PyQA
            VMAT.show_img_results(self.vmat, self.qual_mplw_left, self.qual_mplw_right)

            # draw on quantitative widget
            VMAT.show_img_results(self.vmat, self.quan_mplw)

            self.statusbar.showMessage("Found! See Analysis tabs", 3000)
            self.progressBar.setValue(0)

            # Post results
            self.results_browser.setText(self.vmat.get_string_results())
        else:
            raise NameError("Analysis type not properly set")

        # show applicable widgets
        self.qual_mplw_left.setVisible(True)
        self.qual_mplw_right.setVisible(True)
        self.quan_mplw.setVisible(True)
        # Switch to Qualitative View tab
        self.main_widget.setCurrentIndex(3)


def set_bground_on_flag(QtWidg, flag):
    """Set the background color of the given widget to green if the flag is true and red if the flag is false."""
    palette = QtGui.QPalette()
    if flag:
        brush = QtGui.QBrush(QtGui.QColor(0, 255, 0))  # Green
    else:
        brush = QtGui.QBrush(QtGui.QColor(255, 0, 0))  # Red
    brush.setStyle(QtCore.Qt.SolidPattern)
    palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
    QtWidg.setPalette(palette)






if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = PyQA(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

    