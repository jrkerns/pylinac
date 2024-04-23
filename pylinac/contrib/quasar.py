from typing import List

from .. import StandardImagingFC2
from ..core.geometry import Point
from ..metrics.image import SizedDiskLocator


class QuasarLightRadScaling(StandardImagingFC2):
    """A light/rad and also scaling analysis for the Quasar phantom. The user
    also uses custom edge blocks with offset BBs to detect the light corners.
    It's a mix of scaling validation and light/rad."""

    common_name = "Quasar Light/Rad Scaling"
    bb_sampling_box_size_mm = 10
    bb_size_mm = 5
    field_strip_width_mm = 20
    light_rad_bb_offset_mm = 11
    scaling_centers: List[Point]

    def analyze(
        self, invert: bool = False, fwxm: int = 50, bb_edge_threshold_mm: float = 10
    ) -> None:
        """Analyze the image for the light/rad and scaling"""
        super().analyze(
            invert=invert, fwxm=fwxm, bb_edge_threshold_mm=bb_edge_threshold_mm
        )
        self.scaling_centers = self._detect_scaling_centers()

    def _determine_bb_set(self, fwxm: int) -> dict:
        """We determine the BB set to use for the CAX (the light/rad part is separate).
        We do this by first finding the field edges and then offsetting inward by 10 mm.
        """
        fs_y = self.field_width_y / 2
        fs_x = self.field_width_x / 2
        positions_offsets = {
            "TL": (
                -fs_x + self.light_rad_bb_offset_mm,
                -fs_y + self.light_rad_bb_offset_mm,
            ),
            "BL": (
                -fs_x + self.light_rad_bb_offset_mm,
                fs_y - self.light_rad_bb_offset_mm,
            ),
            "TR": (
                fs_x - self.light_rad_bb_offset_mm,
                fs_y - self.light_rad_bb_offset_mm,
            ),
            "BR": (
                fs_x - self.light_rad_bb_offset_mm,
                -fs_y + self.light_rad_bb_offset_mm,
            ),
        }
        return positions_offsets

    def _detect_scaling_centers(self) -> List[Point]:
        """Sample a 10x10mm square about each BB to detect it. Adjustable using self.bb_sampling_box_size_mm"""
        scaling_centers = self.image.compute(
            SizedDiskLocator.from_center_physical(
                expected_position_mm=Point(0, 0),
                search_window_mm=(35, 35),
                radius_mm=self.bb_size_mm / 2,
                radius_tolerance_mm=self.bb_size_mm / 2,
                min_number=5,
                max_number=5,
                min_separation_mm=4,
            )
        )
        return scaling_centers
