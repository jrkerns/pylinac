import dataclasses
import io
import webbrowser
from pathlib import Path
from typing import Callable, Dict, Optional, Union

import numpy as np
from matplotlib import pyplot as plt
from pylinac.core.profile import CollapsedCircleProfile
from pylinac.core.roi import DiskROI
from pylinac.core.utilities import ResultBase
from pylinac.ct import CatPhanBase, CatPhanModule, Slice

from pylinac.core import pdf


@dataclasses.dataclass
class TomoCheeseResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    origin_slice: int  #:
    num_images: int  #:
    phantom_roll: float  #:
    roi_1: dict  #:
    roi_2: dict  #:
    roi_3: dict  #:
    roi_4: dict  #:
    roi_5: dict  #:
    roi_6: dict  #:
    roi_7: dict  #:
    roi_8: dict  #:
    roi_9: dict  #:
    roi_10: dict  #:
    roi_11: dict  #:
    roi_12: dict  #:
    roi_13: dict  #:
    roi_14: dict  #:
    roi_15: dict  #:
    roi_16: dict  #:
    roi_17: dict  #:
    roi_18: dict  #:
    roi_19: dict  #:
    roi_20: dict  #:


class TomoCheeseModule(CatPhanModule):
    """The pluggable module with user-accessible holes.

    The ROIs of the inner circle are ~45 degrees apart. The ROIs of the outer circle are ~30 degrees apart
    """

    attr_name = "cheese"
    common_name = "Tomo Cheese"
    inner_roi_dist_mm = 65
    outer_roi_dist_mm = 110
    roi_radius_mm = 12
    rois: Dict[str, DiskROI]
    roi_settings = {
        "1": {
            "angle": -75,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "2": {
            "angle": -67.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "3": {
            "angle": -45,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "4": {
            "angle": -22.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "5": {
            "angle": -15,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "6": {
            "angle": 15,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "7": {
            "angle": 22.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "8": {
            "angle": 45,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "9": {
            "angle": 67.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "10": {
            "angle": 75,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "11": {
            "angle": 105,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "12": {
            "angle": 112.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "13": {
            "angle": 135,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "14": {
            "angle": 157.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "15": {
            "angle": 165,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "16": {
            "angle": -165,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "17": {
            "angle": -157.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "18": {
            "angle": -135,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "19": {
            "angle": -112.5,
            "distance": inner_roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "20": {
            "angle": -105,
            "distance": outer_roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }

    def _setup_rois(self) -> None:
        # unlike its super, we use simple disk ROIs as we're not doing complicated things.
        for name, setting in self.roi_settings.items():
            self.rois[name] = DiskROI(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
            )

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis. We add the ROI # to help the user differentiate"""
        for name, roi in self.rois.items():
            roi.plot2axes(axis, edgecolor="blue", text=name)


class TomoCheese(CatPhanBase):
    """A class for analyzing the TomoTherapy 'Cheese' Phantom containing insert holes and plugs for HU analysis."""

    _model = "Tomo Cheese"
    _demo_url = "TomoCheese.zip"
    air_bubble_radius_mm = 14
    localization_radius = 110
    min_num_images = 10
    catphan_radius_mm = 150
    roi_config: dict
    hu: TomoCheeseModule
    modules = {
        TomoCheeseModule: {"offset": 0},
    }

    @staticmethod
    def run_demo(show: bool = True):
        """Run the Tomotherapy Cheese demo"""
        cheese = TomoCheese.from_demo_images()
        cheese.analyze()
        print(cheese.results())
        cheese.plot_analyzed_image(show)

    def analyze(self, roi_config: Optional[dict] = None) -> None:
        """Analyze the Tomo Cheese phantom.

        Parameters
        ----------
        roi_config : dict
            The configuration of the ROIs, specifically the known densities.
        """
        self.localize()
        self.hu = TomoCheeseModule(self)
        self.roi_config = roi_config

    def find_phantom_roll(self, func: Optional[Callable] = None) -> float:
        """Examine the phantom for the maximum HU delta insert position. Roll the phantom by the
        measured angle to the nearest nominal angle if nearby. If not nearby, default to 0
        """
        # get edges and make ROIs from it
        slice = Slice(self, self.origin_slice)
        circle = CollapsedCircleProfile(
            slice.phan_center,
            self.localization_radius / self.mm_per_pixel,
            slice.image.array,
            ccw=False,
            width_ratio=0.05,
            num_profiles=5,
        )
        peak_idxs, _ = circle.find_fwxm_peaks(max_number=1)
        if peak_idxs:
            angle = peak_idxs[0] / len(circle) * 360
            # see if angle is near an ROI node
            shifts = [angle - a for a in np.arange(15, 360, 30)]
            min_shift = shifts[np.argmin([abs(shift) for shift in shifts])]
            if -5 < min_shift < 5:
                return min_shift
            else:
                print(
                    f"Detected shift of {min_shift} was >5 degrees; automatic roll compensation aborted. Setting roll to 0."
                )
                return 0
        else:
            print(
                "No low-HU regions found in the outer ROI circle; automatic roll compensation aborted. Setting roll to 0."
            )
            return 0

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs: dict) -> None:
        """Plot the images used in the calculation and summary data.

        Parameters
        ----------
        show : bool
            Whether to plot the image or not.
        plt_kwargs : dict
            Keyword args passed to the plt.figure() method. Allows one to set things like figure size.
        """
        fig, ax = plt.subplots(**plt_kwargs)
        self.hu.plot(ax)
        plt.tight_layout()
        if show:
            plt.show()

    def plot_density_curve(self, show: bool = True, **plt_kwargs: dict):
        """Plot the densities of the ROIs vs the measured HU. This will sort the ROIs by measured HU before plotting.

        Parameters
        ----------
        show : bool
            Whether to plot the image or not.
        plt_kwargs : dict
            Keyword args passed to the plt.figure() method. Allows one to set things like figure size.
        """
        if not self.roi_config:
            raise ValueError(
                "No ROI density configuration was passed to the analyze method. Re-analyze with densities first."
            )
        xs = []
        ys = []
        for roi_num, roi_data in self.roi_config.items():
            xs.append(self.hu.rois[roi_num].pixel_value)
            ys.append(roi_data["density"])
        # sort by HU so it looks like a normal curve; ROI densities can be out of order
        sorted_args = np.argsort(xs)
        xs = np.array(xs)[sorted_args]
        ys = np.array(ys)[sorted_args]
        # plot
        fig, ax = plt.subplots(**plt_kwargs)
        ax.plot(xs, ys)
        ax.set_title("Density vs HU curve")
        ax.set_ylabel("Density")
        ax.set_xlabel("HU")
        ax.grid("on")
        plt.tight_layout()
        if show:
            plt.show()

    def results(self, as_list: bool = False) -> Union[str, list[str]]:
        """Return the results of the analysis as a string. Use with print().

        Parameters
        ----------
        as_list : bool
            Whether to return as a list of strings vs single string. Pretty much for internal usage.
        """
        results = [
            " - TomoTherapy Cheese Phantom Analysis - ",
            " - HU Module - ",
        ]
        results += [
            f"ROI {name} median: {roi.pixel_value:.1f}, stdev: {roi.std:.1f}"
            for name, roi in self.hu.rois.items()
        ]
        if as_list:
            return results
        else:
            return "\n".join(results)

    def results_data(self, as_dict: bool = False) -> Union[TomoCheeseResult, dict]:
        """Return the results of the analysis as a structure dataclass"""
        data = TomoCheeseResult(
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            phantom_roll=self.catphan_roll,
            roi_1=self.hu.rois["1"].as_dict(),
            roi_2=self.hu.rois["2"].as_dict(),
            roi_3=self.hu.rois["3"].as_dict(),
            roi_4=self.hu.rois["4"].as_dict(),
            roi_5=self.hu.rois["5"].as_dict(),
            roi_6=self.hu.rois["6"].as_dict(),
            roi_7=self.hu.rois["7"].as_dict(),
            roi_8=self.hu.rois["8"].as_dict(),
            roi_9=self.hu.rois["9"].as_dict(),
            roi_10=self.hu.rois["10"].as_dict(),
            roi_11=self.hu.rois["11"].as_dict(),
            roi_12=self.hu.rois["12"].as_dict(),
            roi_13=self.hu.rois["13"].as_dict(),
            roi_14=self.hu.rois["14"].as_dict(),
            roi_15=self.hu.rois["15"].as_dict(),
            roi_16=self.hu.rois["16"].as_dict(),
            roi_17=self.hu.rois["17"].as_dict(),
            roi_18=self.hu.rois["18"].as_dict(),
            roi_19=self.hu.rois["19"].as_dict(),
            roi_20=self.hu.rois["20"].as_dict(),
        )

        if as_dict:
            return dataclasses.asdict(data)
        return data

    def publish_pdf(
        self,
        filename: Union[str, Path],
        notes: Optional[str] = None,
        open_file: bool = False,
        metadata: Optional[dict] = None,
        logo: Optional[Union[Path, str]] = None,
    ) -> None:
        """Publish (print) a PDF containing the analysis and quantitative results.

        Parameters
        ----------
        filename : (str, file-like object}
            The file to write the results to.
        notes : str, list of strings
            Text; if str, prints single line.
            If list of strings, each list item is printed on its own line.
        open_file : bool
            Whether to open the file using the default program after creation.
        metadata : dict
            Extra data to be passed and shown in the PDF. The key and value will be shown with a colon.
            E.g. passing {'Author': 'James', 'Unit': 'TrueBeam'} would result in text in the PDF like:
            --------------
            Author: James
            Unit: TrueBeam
            --------------
        logo: Path, str
            A custom logo to use in the PDF report. If nothing is passed, the default pylinac logo is used.
        """
        canvas = pdf.PylinacCanvas(
            filename,
            page_title="TomoTherapy Cheese Phantom",
            metadata=metadata,
            logo=logo,
        )
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        canvas.add_text(text=self.results(as_list=True), location=(3, 23), font_size=16)
        data = io.BytesIO()
        self.save_analyzed_image(data)
        canvas.add_new_page()
        canvas.add_image(data, location=(0, 4), dimensions=(22, 22))
        canvas.finish()
        if open_file:
            webbrowser.open(filename)

    def save_analyzed_subimage(self) -> None:
        raise NotImplementedError("There are no sub-images for the Tomo Cheese phantom")

    def plot_analyzed_subimage(self) -> None:
        raise NotImplementedError("There are no sub-images for the Tomo Cheese phantom")
