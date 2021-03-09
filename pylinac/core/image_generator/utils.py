import copy
from typing import Union, List, Tuple
import os.path as osp

from .layers import *
from .simulators import Simulator
from ..geometry import cos, sin


def generate_picketfence(simulator: Simulator, field_layer: Union[FilterFreeFieldLayer, FilteredFieldLayer],
                         file_out: str,
                         final_layers: List[Layer] = None,
                         pickets: int = 11, picket_spacing_mm: float = 20, picket_width_mm: int = 2,
                         picket_height_mm=300, gantry_angle=0) -> None:
    """Create a mock picket fence image"""
    picket_pos_mm = range(-int((pickets - 1) * picket_spacing_mm / 2),
                          int((pickets - 1) * picket_spacing_mm / 2) + 1,
                          picket_spacing_mm)
    for pos in picket_pos_mm:
        simulator.add_layer(field_layer((picket_height_mm, picket_width_mm), cax_offset_mm=(0, pos)))
    if final_layers is not None:
        for layer in final_layers:
            simulator.add_layer(layer)
    simulator.generate_dicom(file_out, gantry_angle=gantry_angle)


def generate_winstonlutz(simulator: Simulator, field_layer: Layer, dir_out: str, field_size_mm=(30, 30),
                         final_layers: List[Layer]=None, bb_size_mm=5, offset_mm_left=0, offset_mm_up=0, offset_mm_in=0,
                         image_axes: List[Tuple[int, int, int]]=[(0, 0, 0), (90, 0, 0), (180, 0, 0), (270, 0, 0)],
                         gantry_tilt=0, gantry_sag=0) -> List[str]:
    """Create a mock set of WL images, simulating gantry sag effects. Produces one image for each item in image_axes.

    Parameters
    ----------
    image_axes
        List of axes for the images. Sequence is (Gantry, Coll, Couch).
    """
    file_names = []
    for (gantry, coll, couch) in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(field_layer(field_size_mm=field_size_mm, cax_offset_mm=(gantry_tilt*cos(gantry), gantry_sag*sin(gantry))))
        sim_single.add_layer(PerfectBBLayer(cax_offset_mm=(offset_mm_in,
                                                           offset_mm_left*cos(gantry)+offset_mm_up*sin(gantry))))
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; BB={bb_size_mm}mm @ left={offset_mm_left}, in={offset_mm_in}, up={offset_mm_up}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(osp.join(dir_out, file_name), gantry_angle=gantry, coll_angle=coll, table_angle=couch)
        file_names.append(file_name)
    return file_names


