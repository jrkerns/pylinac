import dataclasses
import io
import webbrowser
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Union, List, Dict, Optional

from matplotlib import pyplot as plt

from pylinac.core import pdf
from pylinac.core.mtf import MTF
from pylinac.core.roi import HighContrastDiskROI
from pylinac.core.utilities import ResultBase
from pylinac.ct import CatPhanBase, CatPhanModule

UNIFORMITY_MODULE_OFFSET_MM = 70
SPATIAL_RESOLUTION_MODULE_OFFSET_MM = 100
LOW_CONTRAST_MODULE_OFFSET_MM = 30


class CTModule(CatPhanModule):
    attr_name = "ct_calibration_module"
    roi_dist_mm = 63
    roi_radius_mm = 10
    roi_settings = {
        "Air": {"angle": 45, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Poly": {"angle": 225, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Acrylic": {"angle": 135, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Bone": {"angle": -45, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Water": {"angle": 180, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    window_min = -200
    window_max = 200


@dataclass
class CTModuleOutput:
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    offset: int
    roi_distance_from_center_mm: int
    roi_radius_mm: int
    roi_settings: dict
    rois: dict


class UniformityModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "uniformity_module"
    common_name = "HU Uniformity"
    roi_dist_mm = 66
    roi_radius_mm = 11
    roi_settings = {
        "Top": {"angle": -90, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Right": {"angle": 0, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Bottom": {"angle": 90, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Left": {"angle": 180, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "Center": {"angle": 0, "distance": 0, "radius": roi_radius_mm},
    }
    window_min = -50
    window_max = 50


@dataclass
class UniformityModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    center_roi_stdev: float


class SpatialResolutionModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "spatial_resolution_module"
    common_name = "Spatial Resolution"
    rois: Dict[str, HighContrastDiskROI]
    roi_dist_mm = 70
    roi_radius_mm = 6
    roi_settings = {
        "10oclock": {
            "angle": -135,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.4,
        },
        "9oclock": {
            "angle": -180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.5,
        },
        "7oclock": {
            "angle": 135,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.6,
        },
        "6oclock": {
            "angle": 90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.7,
        },
        "4oclock": {
            "angle": 45,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.8,
        },
        "3oclock": {
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 0.9,
        },
        "2oclock": {
            "angle": -45,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 1.0,
        },
        "12oclock": {
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
            "lp/mm": 1.2,
        },
    }

    def _setup_rois(self) -> None:
        for name, setting in self.roi_settings.items():
            self.rois[name] = HighContrastDiskROI(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_threshold=1.0,  # fixed to 1 so everything passes. We aren't evaluating pass/fail here
            )

    @property
    def mtf(self) -> MTF:
        spacings = [roi["lp/mm"] for roi in self.roi_settings.values()]
        return MTF.from_high_contrast_diskset(
            spacings=spacings, diskset=list(self.rois.values())
        )

    def plot_rois(self, axis: plt.Axes) -> None:
        """Plot the ROIs to the axis. Override to set the color"""
        for roi, mtf in zip(self.rois.values(), self.mtf.norm_mtfs.values()):
            roi.plot2axes(axis, edgecolor='g')


@dataclass
class SpatialResolutionModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    lpmm_to_rmtf: dict


class LowContrastModule(CatPhanModule):
    """Class for analysis of the Uniformity slice of the CTP module. Measures 5 ROIs around the slice that
    should all be close to the same value.
    """

    attr_name = "low_contrast_module"
    common_name = "Low Contrast"
    roi_dist_mm = 60
    roi_radius_mm = 6
    nominal_value = 0
    roi_settings = {
        "ROI": {"angle": -90, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    background_roi_settings = {
        "ROI": {"angle": -115, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }
    window_min = 50
    window_max = 150

    def cnr(self) -> float:
        """Given in the guidance doc as |A-B|/SD where A is the contrast ROI, B is the background, and SD is stdev of B"""
        return (
            abs(self.rois["ROI"].pixel_value - self.background_rois["ROI"].pixel_value)
            / self.background_rois["ROI"].std
        )


@dataclass
class LowContrastModuleOutput(CTModuleOutput):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    cnr: float


@dataclass
class ACRCTResult(ResultBase):
    """This class should not be called directly. It is returned by the ``results_data()`` method.
    It is a dataclass under the hood and thus comes with all the dunder magic.

    Use the following attributes as normal class attributes."""

    phantom_model: str  #:
    phantom_roll_deg: float  #:
    origin_slice: int  #:
    num_images: int  #:
    ct_module: CTModuleOutput  #:
    uniformity_module: UniformityModuleOutput  #:
    low_contrast_module: LowContrastModuleOutput  #:
    spatial_resolution_module: SpatialResolutionModuleOutput  #:


class CT464(CatPhanBase):
    _model = "ACR CT 464"
    catphan_radius_mm = 100
    air_bubble_radius_mm = 18
    min_num_images = 4
    localization_radius = 70
    results: dict
    ct_calibration_module = CTModule
    low_contrast_module = LowContrastModule
    spatial_resolution_module = SpatialResolutionModule
    uniformity_module = UniformityModule
    clear_borders = False

    def plot_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `plot_images`")

    def save_analyzed_subimage(self, *args, **kwargs):
        raise NotImplementedError("Use `save_images`")

    def analyze(self) -> None:
        """Analyze the ACR CT phantom"""
        self.ct_calibration_module = self.ct_calibration_module(self, offset=0, clear_borders=self.clear_borders)
        self.uniformity_module = self.uniformity_module(
            self, offset=UNIFORMITY_MODULE_OFFSET_MM, clear_borders=self.clear_borders
        )
        self.spatial_resolution_module = self.spatial_resolution_module(
            self,
            offset=SPATIAL_RESOLUTION_MODULE_OFFSET_MM, clear_borders=self.clear_borders
        )
        self.low_contrast_module = self.low_contrast_module(
            self, offset=LOW_CONTRAST_MODULE_OFFSET_MM, clear_borders=self.clear_borders
        )

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs) -> plt.Figure:
        """Plot the analyzed image

        Parameters
        ----------
        show
            Whether to show the image.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        # set up grid and axes
        fig = plt.figure(**plt_kwargs)
        grid_size = (2, 3)
        hu_ax = plt.subplot2grid(grid_size, (0, 0))
        self.ct_calibration_module.plot(hu_ax)
        unif_ax = plt.subplot2grid(grid_size, (0, 1))
        self.uniformity_module.plot(unif_ax)
        sr_ax = plt.subplot2grid(grid_size, (0, 2))
        self.spatial_resolution_module.plot(sr_ax)
        locon_ax = plt.subplot2grid(grid_size, (1, 0))
        self.low_contrast_module.plot(locon_ax)
        spatial_res_graph = plt.subplot2grid(grid_size, (1, 1), colspan=2)
        self.spatial_resolution_module.mtf.plot(spatial_res_graph)

        # finish up
        plt.tight_layout()
        if show:
            plt.show()
        return fig

    def save_analyzed_image(self, filename: Union[str, Path, BytesIO], **plt_kwargs) -> None:
        """Save the analyzed image to disk or stream

        Parameters
        ----------
        filename
            Where to save the image to
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        fig = self.plot_analyzed_image(show=False, **plt_kwargs)
        fig.savefig(filename)

    def plot_images(self, show: bool = True, **plt_kwargs) -> Dict[str, plt.Figure]:
        """Plot all the individual images separately

        Parameters
        ----------
        show
            Whether to show the images.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        figs = {}
        # plot the images
        modules = {'hu': self.ct_calibration_module, 'uniformity': self.uniformity_module,
                   'spatial resolution': self.spatial_resolution_module,
                   'low contrast': self.low_contrast_module}
        for key, module in modules.items():
            fig, ax = plt.subplots(**plt_kwargs)
            module.plot(ax)
            figs[key] = fig
        # plot the one-off MTF image
        fig, ax = plt.subplots(**plt_kwargs)
        figs['mtf'] = fig
        self.spatial_resolution_module.mtf.plot(ax)
        plt.tight_layout()

        if show:
            plt.show()
        return figs

    def save_images(self, directory: Optional[Union[Path, str]] = None, to_stream: bool = False, **plt_kwargs) -> List[Union[Path, BytesIO]]:
        """Save separate images to disk or stream.

        Parameters
        ----------
        directory
            The directory to write the images to. If None, will use current working directory
        to_stream
            Whether to write to stream or disk. If True, will return streams. Directory is ignored in that scenario.
        plt_kwargs
            Keywords to pass to matplotlib for figure customization.
        """
        figs = self.plot_images(show=False, **plt_kwargs)
        paths = []
        for name, fig in figs.items():
            if to_stream:
                path = io.BytesIO()
            else:
                destination = Path(directory) or Path.cwd()
                path = (destination / name).with_suffix('.png').absolute()
            fig.savefig(path)
            paths.append(path)
        return paths

    def find_phantom_roll(self, func=lambda roi: roi.bbox_area) -> float:
        """Determine the "roll" of the phantom.

        Only difference of base method is that we sort the ROIs by size,
        not by being in the center since the two we're looking for are both right-sided.
        """
        return super().find_phantom_roll(func)

    def results(self) -> str:
        """Return the results of the analysis as a string. Use with print()."""
        string = (
            f"\n - ACR CT 464 QA Test - \n"
            f"HU ROIs: {self.ct_calibration_module.roi_vals_as_str}\n"
            f"Contrast to Noise Ratio: {self.low_contrast_module.cnr():2.2f}\n"
            f"Uniformity ROIs: {self.uniformity_module.roi_vals_as_str}\n"
            f'Uniformity Center ROI standard deviation: {self.uniformity_module.rois["Center"].std:2.2f}\n'
            f"MTF 50% (lp/mm): {self.spatial_resolution_module.mtf.relative_resolution(50):2.2f}\n"
        )
        return string

    def results_data(self, as_dict=False) -> Union[ACRCTResult, dict]:
        data = ACRCTResult(
            phantom_model="ACR CT 464",
            phantom_roll_deg=self.catphan_roll,
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            ct_module=CTModuleOutput(
                offset=0,
                roi_distance_from_center_mm=self.ct_calibration_module.roi_dist_mm,
                roi_radius_mm=self.ct_calibration_module.roi_radius_mm,
                roi_settings=self.ct_calibration_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.ct_calibration_module.rois.items()
                },
            ),
            uniformity_module=UniformityModuleOutput(
                offset=70,
                roi_distance_from_center_mm=self.uniformity_module.roi_dist_mm,
                roi_radius_mm=self.uniformity_module.roi_radius_mm,
                roi_settings=self.uniformity_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.uniformity_module.rois.items()
                },
                center_roi_stdev=self.uniformity_module.rois["Center"].std,
            ),
            spatial_resolution_module=SpatialResolutionModuleOutput(
                offset=0,
                roi_distance_from_center_mm=self.spatial_resolution_module.roi_dist_mm,
                roi_radius_mm=self.spatial_resolution_module.roi_radius_mm,
                roi_settings=self.spatial_resolution_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.spatial_resolution_module.rois.items()
                },
                lpmm_to_rmtf=self.spatial_resolution_module.mtf.norm_mtfs,
            ),
            low_contrast_module=LowContrastModuleOutput(
                offset=0,
                roi_distance_from_center_mm=self.low_contrast_module.roi_dist_mm,
                roi_radius_mm=self.low_contrast_module.roi_radius_mm,
                roi_settings=self.low_contrast_module.roi_settings,
                rois={
                    name: roi.pixel_value
                    for name, roi in self.low_contrast_module.rois.items()
                },
                cnr=self.low_contrast_module.cnr(),
            ),
        )
        if as_dict:
            return dataclasses.asdict(data)
        return data

    def publish_pdf(self, filename: Union[str, Path], notes: Optional[str] = None, open_file: bool = False,
                    metadata: Optional[dict] = None) -> None:
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
        """
        analysis_title = f'{self._model} Analysis'
        texts = [
            ' - ACR CT 464 Results - ',
             f'HU Linearity ROIs: {self.ct_calibration_module.roi_vals_as_str}',
             f'Low contrast visibility: {self.low_contrast_module.cnr():2.2f}',
             f'Uniformity ROIs: {self.uniformity_module.roi_vals_as_str}',
        ]
        analysis_images = self.save_images(to_stream=True)

        canvas = pdf.PylinacCanvas(filename, page_title=analysis_title, metadata=metadata)
        if notes is not None:
            canvas.add_text(text="Notes:", location=(1, 4.5), font_size=14)
            canvas.add_text(text=notes, location=(1, 4))

        for idx, text in enumerate(texts):
            canvas.add_text(text=text, location=(1.5, 23-idx*0.5))
        for page, img in enumerate(analysis_images):
            canvas.add_new_page()
            canvas.add_image(img, location=(1, 5), dimensions=(18, 18))
        canvas.finish()

        if open_file:
            webbrowser.open(filename)
