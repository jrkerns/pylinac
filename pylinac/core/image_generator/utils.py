import copy
import os
import os.path as osp
from typing import Sequence, Union, List, Tuple, Type, Optional

from . import GaussianFilterLayer
from .layers import (
    FilteredFieldLayer,
    FilterFreeFieldLayer,
    Layer,
    PerfectBBLayer,
    PerfectFieldLayer,
    FilterFreeConeLayer,
    PerfectConeLayer,
)
from .simulators import Simulator
from ..geometry import cos, sin
from ...picketfence import Orientation


def generate_lightrad(
    file_out: str,
    simulator: Simulator,
    field_layer: Type[
        Union[FilterFreeFieldLayer, FilteredFieldLayer, PerfectFieldLayer]
    ],
    field_size_mm: (float, float) = (150, 150),
    cax_offset_mm: (float, float) = (0, 0),
    final_layers: List[Layer] = [
        GaussianFilterLayer(),
    ],
    bb_size_mm: float = 3,
    bb_positions: Tuple[Tuple[float, float], ...] = (
        (-40, -40),
        (-40, 40),
        (40, -40),
        (40, 40),
        (-65, -65),
        (-65, 65),
        (65, -65),
        (65, 65),
    ),
) -> None:
    """Create a mock light/rad image with BBs.

    Parameters
    ----------
    simulator
        The image simulator
    field_layer
        The primary field layer
    file_out
        The name of the file to save the DICOM file to.
    final_layers
        Optional layers to apply at the end of the procedure. Useful for noise or blurring.
    bb_size_mm
        The size of the phantom BBs
    bb_positions
        The position of the BBs relative to the CAX.
    """
    # open field layer
    simulator.add_layer(
        field_layer(field_size_mm=field_size_mm, cax_offset_mm=cax_offset_mm)
    )
    # bbs
    for bb in bb_positions:
        simulator.add_layer(PerfectBBLayer(bb_size_mm=bb_size_mm, cax_offset_mm=bb))

    if final_layers is not None:
        for layer in final_layers:
            simulator.add_layer(layer)
    simulator.generate_dicom(file_out)


def generate_picketfence(
    simulator: Simulator,
    field_layer: Type[
        Union[FilterFreeFieldLayer, FilteredFieldLayer, PerfectFieldLayer]
    ],
    file_out: str,
    final_layers: List[Layer] = None,
    pickets: int = 11,
    picket_spacing_mm: float = 20,
    picket_width_mm: int = 2,
    picket_height_mm: int = 300,
    gantry_angle: int = 0,
    orientation=Orientation.UP_DOWN,
    picket_offset_error: Optional[Sequence] = None,
) -> None:
    """Create a mock picket fence image. Will always be up-down.

    Parameters
    ----------
    simulator
        The image simulator
    field_layer
        The primary field layer
    file_out
        The name of the file to save the DICOM file to.
    final_layers
        Optional layers to apply at the end of the procedure. Useful for noise or blurring.
    pickets
        The number of pickets
    picket_spacing_mm
        The space between pickets
    picket_width_mm
        Picket width parallel to leaf motion
    picket_height_mm
        Picket height parallel to leaf motion
    gantry_angle
        Gantry angle; sets the DICOM tag.
    """
    picket_pos_mm = range(
        -int((pickets - 1) * picket_spacing_mm / 2),
        int((pickets - 1) * picket_spacing_mm / 2) + 1,
        picket_spacing_mm,
    )
    for idx, pos in enumerate(picket_pos_mm):
        if picket_offset_error is not None:
            if len(picket_offset_error) != pickets:
                raise ValueError(
                    "The length of the error array must be the same as the number of pickets."
                )
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


def generate_winstonlutz(
    simulator: Simulator,
    field_layer: Type[Layer],
    dir_out: str,
    field_size_mm: Tuple[float, float] = (30, 30),
    final_layers: Optional[List[Layer]] = None,
    bb_size_mm: float = 5,
    offset_mm_left: float = 0,
    offset_mm_up: float = 0,
    offset_mm_in: float = 0,
    image_axes: List[Tuple[int, int, int]] = [
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ],
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
) -> List[str]:
    """Create a mock set of WL images, simulating gantry sag effects. Produces one image for each item in image_axes.

    Parameters
    ----------
    simulator
        The image simulator
    field_layer
        The primary field layer simulating radiation
    dir_out
        The directory to save the images to.
    field_size_mm
        The field size of the radiation field in mm
    final_layers
        Layers to apply after generating the primary field and BB layer. Useful for blurring or adding noise.
    bb_size_mm
        The size of the BB. Must be positive.
    offset_mm_left
        How far left (lat) to set the BB. Can be positive or negative.
    offset_mm_up
        How far up (vert) to set the BB. Can be positive or negative.
    offset_mm_in
        How far in (long) to set the BB. Can be positive or negative.
    image_axes
        List of axis values for the images. Sequence is (Gantry, Coll, Couch).
    gantry_tilt
        The tilt of the gantry that affects the position at 0 and 180. Simulates a simple cosine function.
    gantry_sag
        The sag of the gantry that affects the position at gantry=90 and 270. Simulates a simple sine function.
    clean_dir
        Whether to clean out the output directory. Useful when iterating.
    """
    if not osp.isdir(dir_out):
        os.mkdir(dir_out)
    if clean_dir:
        for pdir, _, files in os.walk(dir_out):
            [os.remove(osp.join(pdir, f)) for f in files]
    file_names = []
    for (gantry, coll, couch) in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(
            field_layer(
                field_size_mm=field_size_mm,
                cax_offset_mm=(gantry_tilt * cos(gantry), gantry_sag * sin(gantry)),
            )
        )
        sim_single.add_layer(
            PerfectBBLayer(
                cax_offset_mm=(
                    -offset_mm_in,
                    -offset_mm_left * cos(gantry) - offset_mm_up * sin(gantry),
                ),
                bb_size_mm=bb_size_mm,
            )
        )
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Field={field_size_mm}mm; BB={bb_size_mm}mm @ left={offset_mm_left}, in={offset_mm_in}, up={offset_mm_up}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(
            osp.join(dir_out, file_name),
            gantry_angle=gantry,
            coll_angle=coll,
            table_angle=couch,
        )
        file_names.append(file_name)
    return file_names


def generate_winstonlutz_cone(
    simulator: Simulator,
    cone_layer: Union[Type[FilterFreeConeLayer], Type[PerfectConeLayer]],
    dir_out: str,
    cone_size_mm: float = 17.5,
    final_layers: Optional[List[Layer]] = None,
    bb_size_mm: float = 5,
    offset_mm_left: float = 0,
    offset_mm_up: float = 0,
    offset_mm_in: float = 0,
    image_axes: List[Tuple[int, int, int]] = [
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ],
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
) -> List[str]:
    """Create a mock set of WL images with a cone field, simulating gantry sag effects. Produces one image for each item in image_axes.

    Parameters
    ----------
    simulator
        The image simulator
    cone_layer
        The primary field layer simulating radiation
    dir_out
        The directory to save the images to.
    cone_size_mm
        The field size of the radiation field in mm
    final_layers
        Layers to apply after generating the primary field and BB layer. Useful for blurring or adding noise.
    bb_size_mm
        The size of the BB. Must be positive.
    offset_mm_left
        How far left (lat) to set the BB. Can be positive or negative.
    offset_mm_up
        How far up (vert) to set the BB. Can be positive or negative.
    offset_mm_in
        How far in (long) to set the BB. Can be positive or negative.
    image_axes
        List of axis values for the images. Sequence is (Gantry, Coll, Couch).
    gantry_tilt
        The tilt of the gantry that affects the position at 0 and 180. Simulates a simple cosine function.
    gantry_sag
        The sag of the gantry that affects the position at gantry=90 and 270. Simulates a simple sine function.
    clean_dir
        Whether to clean out the output directory. Useful when iterating.
    """
    if not osp.isdir(dir_out):
        os.mkdir(dir_out)
    if clean_dir:
        for pdir, _, files in os.walk(dir_out):
            [os.remove(osp.join(pdir, f)) for f in files]
    file_names = []
    for (gantry, coll, couch) in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(
            cone_layer(
                cone_size_mm=cone_size_mm,
                cax_offset_mm=(gantry_tilt * cos(gantry), gantry_sag * sin(gantry)),
            )
        )
        sim_single.add_layer(
            PerfectBBLayer(
                cax_offset_mm=(
                    -offset_mm_in,
                    -offset_mm_left * cos(gantry) - offset_mm_up * sin(gantry),
                )
            )
        )
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Cone={cone_size_mm}mm; BB={bb_size_mm}mm @ left={offset_mm_left}, in={offset_mm_in}, up={offset_mm_up}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(
            osp.join(dir_out, file_name),
            gantry_angle=gantry,
            coll_angle=coll,
            table_angle=couch,
        )
        file_names.append(file_name)
    return file_names
