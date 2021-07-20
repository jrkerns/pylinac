import copy
from typing import Union, List, Tuple, Type, Optional, Sequence
import os.path as osp

from .layers import FilteredFieldLayer, FilterFreeFieldLayer, Layer, PerfectBBLayer, PerfectFieldLayer
from .simulators import Simulator
from ..geometry import cos, sin
from ...picketfence import Orientation


def generate_picketfence(simulator: Simulator, field_layer: Type[Union[FilterFreeFieldLayer, FilteredFieldLayer, PerfectFieldLayer]],
                         file_out: str,
                         final_layers: List[Layer] = None,
                         pickets: int = 11, picket_spacing_mm: float = 20, picket_width_mm: int = 2,
                         picket_height_mm: int = 300, gantry_angle: int = 0, orientation=Orientation.UP_DOWN,
                         picket_offset_error: Optional[Sequence] = None) -> None:
    """Create a mock picket fence image"""
    picket_pos_mm = range(-int((pickets - 1) * picket_spacing_mm / 2),
                          int((pickets - 1) * picket_spacing_mm / 2) + 1,
                          picket_spacing_mm)
    for idx, pos in enumerate(picket_pos_mm):
        if picket_offset_error is not None:
            if len(picket_offset_error) != pickets:
                raise ValueError("The length of the error array must be the same as the number of pickets.")
            pos += picket_offset_error[idx]
        if orientation == orientation.UP_DOWN:
            position = (0, pos)
            layout = (picket_height_mm, picket_width_mm)
        else:
            position = (pos, 0)
            layout = (picket_width_mm, picket_height_mm)
        simulator.add_layer(field_layer(layout, cax_offset_mm=position))
    if final_layers is not None:
        for layer in final_layers:
            simulator.add_layer(layer)
    simulator.generate_dicom(file_out, gantry_angle=gantry_angle)


def generate_winstonlutz(simulator: Simulator, field_layer: Type[Layer], dir_out: str, field_size_mm=(30, 30),
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
        sim_single.add_layer(PerfectBBLayer(cax_offset_mm=(-offset_mm_in,
                                                           -offset_mm_left*cos(gantry)-offset_mm_up*sin(gantry)), bb_size_mm=bb_size_mm))
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Field={field_size_mm}mm; BB={bb_size_mm}mm @ left={offset_mm_left}, in={offset_mm_in}, up={offset_mm_up}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(osp.join(dir_out, file_name), gantry_angle=gantry, coll_angle=coll, table_angle=couch)
        file_names.append(file_name)
    return file_names


def generate_winstonlutz_cone(simulator: Simulator, cone_layer: Layer, dir_out: str, cone_size_mm=17.5,
                         final_layers: List[Layer]=None, bb_size_mm=5, offset_mm_left=0, offset_mm_up=0, offset_mm_in=0,
                         image_axes: List[Tuple[int, int, int]]=[(0, 0, 0), (90, 0, 0), (180, 0, 0), (270, 0, 0)],
                         gantry_tilt=0, gantry_sag=0) -> List[str]:
    """Create a mock set of WL images with a cone field, simulating gantry sag effects. Produces one image for each item in image_axes.

    Parameters
    ----------
    image_axes
        List of axes for the images. Sequence is (Gantry, Coll, Couch).
    """
    file_names = []
    for (gantry, coll, couch) in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(cone_layer(cone_size_mm=cone_size_mm, cax_offset_mm=(gantry_tilt*cos(gantry), gantry_sag*sin(gantry))))
        sim_single.add_layer(PerfectBBLayer(cax_offset_mm=(-offset_mm_in,
                                                           -offset_mm_left*cos(gantry)-offset_mm_up*sin(gantry))))
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Cone={cone_size_mm}mm; BB={bb_size_mm}mm @ left={offset_mm_left}, in={offset_mm_in}, up={offset_mm_up}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(osp.join(dir_out, file_name), gantry_angle=gantry, coll_angle=coll, table_angle=couch)
        file_names.append(file_name)
    return file_names

