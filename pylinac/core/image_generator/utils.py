from __future__ import annotations

import copy
import os
import os.path as osp
import random
from typing import Sequence

from ...picketfence import Orientation
from ...winston_lutz import bb_projection_gantry_plane, bb_projection_long
from ..geometry import cos, sin
from . import GaussianFilterLayer
from .layers import (
    FilteredFieldLayer,
    FilterFreeConeLayer,
    FilterFreeFieldLayer,
    Layer,
    PerfectBBLayer,
    PerfectConeLayer,
    PerfectFieldLayer,
)
from .simulators import Simulator


def generate_lightrad(
    file_out: str,
    simulator: Simulator,
    field_layer: type[FilterFreeFieldLayer | FilteredFieldLayer | PerfectFieldLayer],
    field_size_mm: (float, float) = (150, 150),
    cax_offset_mm: (float, float) = (0, 0),
    final_layers: list[Layer] = [
        GaussianFilterLayer(),
    ],
    bb_size_mm: float = 3,
    bb_positions: ((float, float), ...) = (
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
    field_layer: type[FilterFreeFieldLayer | FilteredFieldLayer | PerfectFieldLayer],
    file_out: str,
    final_layers: list[Layer] = None,
    pickets: int = 11,
    picket_spacing_mm: float = 20,
    picket_width_mm: int = 2,
    picket_height_mm: int = 300,
    gantry_angle: int = 0,
    orientation: Orientation = Orientation.UP_DOWN,
    picket_offset_error: Sequence | None = None,
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
    field_layer: type[Layer],
    dir_out: str,
    field_size_mm: tuple[float, float] = (30, 30),
    final_layers: list[Layer] | None = None,
    bb_size_mm: float = 5,
    offset_mm_left: float = 0,
    offset_mm_up: float = 0,
    offset_mm_in: float = 0,
    image_axes: ((int, int, int), ...) = (
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ),
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
    field_alpha: float = 1.0,
    bb_alpha: float = -0.5,
) -> list[str]:
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
    field_alpha
        The normalized alpha (i.e. signal) of the radiation field. Use in combination
        with bb_alpha such that the sum of the two is always <= 1.
    bb_alpha
        The normalized alpha (in the case of the BB think of it as attenuation) of the BB against the radiation field. More negative values
        attenuate (remove signal) more.
    """
    if field_alpha + bb_alpha > 1:
        raise ValueError("field_alpha and bb_alpha must sum to <=1")
    if field_alpha - bb_alpha < 0:
        raise ValueError("field_alpha and bb_alpha must have a sum >=0")
    if not osp.isdir(dir_out):
        os.mkdir(dir_out)
    if clean_dir:
        for pdir, _, files in os.walk(dir_out):
            [os.remove(osp.join(pdir, f)) for f in files]
    file_names = []
    for gantry, coll, couch in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(
            field_layer(
                field_size_mm=field_size_mm,
                cax_offset_mm=(gantry_tilt * cos(gantry), gantry_sag * sin(gantry)),
                alpha=field_alpha,
            )
        )
        # we return the negative because this function
        # will return the offset in PLOTTING space, not coordinate space
        # which is inverted in the long direction
        long_offset = -bb_projection_long(
            offset_in=offset_mm_in,
            offset_up=offset_mm_up,
            offset_left=offset_mm_left,
            sad=1000,
            gantry=gantry,
        )
        gplane_offset = bb_projection_gantry_plane(
            offset_left=offset_mm_left, offset_up=offset_mm_up, sad=1000, gantry=gantry
        )
        sim_single.add_layer(
            PerfectBBLayer(
                cax_offset_mm=(long_offset, gplane_offset),
                bb_size_mm=bb_size_mm,
                alpha=bb_alpha,
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


def generate_winstonlutz_multi_bb_single_field(
    simulator: Simulator,
    field_layer: type[Layer],
    dir_out: str,
    offsets: list[list[float]] | list[dict[str, float]],
    field_size_mm: tuple[float, float] = (30, 30),
    final_layers: list[Layer] | None = None,
    bb_size_mm: float = 5,
    image_axes: ((int, int, int), ...) = (
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ),
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
    jitter_mm: float = 0,
) -> list[str]:
    """Create a mock set of WL images, simulating gantry sag effects. Produces one image for each item in image_axes.
    This will also generate multiple BBs on the image, one per item in `offsets`. Each offset should be a list of
    the shifts of the BB relative to isocenter like so: [<left>, <up>, <in>] OR an arrangement from the WL module.

    Parameters
    ----------
    simulator
        The image simulator
    field_layer
        The primary field layer simulating radiation
    dir_out
        The directory to save the images to.
    offsets
        A list of lists containing the shift of the BBs from iso; each sublist should be a 3-item list/tuple of left, up, in.
        Negative values are acceptable and will go the opposite direction.
    field_size_mm
        The field size of the radiation field in mm
    final_layers
        Layers to apply after generating the primary field and BB layer. Useful for blurring or adding noise.
    bb_size_mm
        The size of the BB. Must be positive.
    image_axes
        List of axis values for the images. Sequence is (Gantry, Coll, Couch).
    gantry_tilt
        The tilt of the gantry in degrees that affects the position at 0 and 180. Simulates a simple cosine function.
    gantry_sag
        The sag of the gantry that affects the position at gantry=90 and 270. Simulates a simple sine function.
    clean_dir
        Whether to clean out the output directory. Useful when iterating.
    jitter_mm
        The amount of jitter to add to the in/left/up location of the BB in MM.
    """
    if not osp.isdir(dir_out):
        os.mkdir(dir_out)
    if clean_dir:
        for pdir, _, files in os.walk(dir_out):
            [os.remove(osp.join(pdir, f)) for f in files]
    file_names = []
    for gantry, coll, couch in image_axes:
        sim_single = copy.copy(simulator)
        sim_single.add_layer(
            field_layer(
                field_size_mm=field_size_mm,
                cax_offset_mm=(gantry_tilt * cos(gantry), gantry_sag * sin(gantry)),
            )
        )
        for offset in offsets:
            if isinstance(offset, dict):
                offset_mm_left = offset["offset_left_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
                offset_mm_up = offset["offset_up_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
                offset_mm_in = -offset["offset_in_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
            else:
                offset_mm_left = offset[0] + random.uniform(-jitter_mm, jitter_mm)
                offset_mm_up = offset[1] + random.uniform(-jitter_mm, jitter_mm)
                offset_mm_in = -offset[2] + random.uniform(
                    -jitter_mm, jitter_mm
                )  # negative because pixels increase as we go out, so to go in we subtract

            long_offset = bb_projection_long(
                offset_in=offset_mm_in,
                offset_up=offset_mm_up,
                offset_left=offset_mm_left,
                sad=1000,
                gantry=gantry,
            )
            gplane_offset = bb_projection_gantry_plane(
                offset_left=offset_mm_left,
                offset_up=offset_mm_up,
                sad=1000,
                gantry=gantry,
            )
            sim_single.add_layer(
                PerfectBBLayer(
                    cax_offset_mm=(
                        long_offset,
                        gplane_offset,
                    ),
                    bb_size_mm=bb_size_mm,
                )
            )
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Field={field_size_mm}mm; BB={bb_size_mm}mm @ left={offset_mm_left:.2f}, in={offset_mm_in:.2f}, up={offset_mm_up:.2f}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
        sim_single.generate_dicom(
            osp.join(dir_out, file_name),
            gantry_angle=gantry,
            coll_angle=coll,
            table_angle=couch,
        )
        file_names.append(file_name)
    return file_names


def generate_winstonlutz_multi_bb_multi_field(
    simulator: Simulator,
    field_layer: type[Layer],
    dir_out: str,
    field_offsets: list[list[float]],
    bb_offsets: list[list[float]] | list[dict[str, float]],
    field_size_mm: tuple[float, float] = (20, 20),
    final_layers: list[Layer] | None = None,
    bb_size_mm: float = 5,
    image_axes: ((int, int, int), ...) = (
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ),
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
    jitter_mm: float = 0,
    align_to_pixels: bool = True,
) -> list[str]:
    """Create a mock set of WL images, simulating gantry sag effects. Produces one image for each item in image_axes.
    This will also generate multiple BBs on the image, one per item in `offsets`. Each offset should be a list of
    the shifts of the BB relative to isocenter like so: [<left>, <up>, <in>] OR an arrangement from the WL module.

    Parameters
    ----------
    simulator
        The image simulator
    field_layer
        The primary field layer simulating radiation
    dir_out
        The directory to save the images to.
    field_offsets
        A list of lists containing the shift of the fields. Format is the same as bb_offsets.
    bb_offsets
        A list of lists containing the shift of the BBs from iso; each sublist should be a 3-item list/tuple of left, up, in.
        Negative values are acceptable and will go the opposite direction.
    field_size_mm
        The field size of the radiation field in mm
    final_layers
        Layers to apply after generating the primary field and BB layer. Useful for blurring or adding noise.
    bb_size_mm
        The size of the BB. Must be positive.
    image_axes
        List of axis values for the images. Sequence is (Gantry, Coll, Couch).
    gantry_tilt
        The tilt of the gantry in degrees that affects the position at 0 and 180. Simulates a simple cosine function.
    gantry_sag
        The sag of the gantry that affects the position at gantry=90 and 270. Simulates a simple sine function.
    clean_dir
        Whether to clean out the output directory. Useful when iterating.
    jitter_mm
        The amount of jitter to add to the in/left/up location of the BB in MM.
    """
    if not osp.isdir(dir_out):
        os.mkdir(dir_out)
    if clean_dir:
        for pdir, _, files in os.walk(dir_out):
            [os.remove(osp.join(pdir, f)) for f in files]
    file_names = []
    for gantry, coll, couch in image_axes:
        sim_single = copy.copy(simulator)
        for field_offset in field_offsets:
            offset_mm_left = field_offset[0] + random.uniform(-jitter_mm, jitter_mm)
            offset_mm_up = field_offset[1] + random.uniform(-jitter_mm, jitter_mm)
            offset_mm_in = -field_offset[2] + random.uniform(
                -jitter_mm, jitter_mm
            )  # negative because pixels increase as we go out, so to go in we subtract
            long_offset = bb_projection_long(
                offset_in=offset_mm_in,
                offset_up=offset_mm_up,
                offset_left=offset_mm_left,
                sad=1000,
                gantry=gantry,
            )
            gplane_offset = bb_projection_gantry_plane(
                offset_left=offset_mm_left,
                offset_up=offset_mm_up,
                sad=1000,
                gantry=gantry,
            )
            long_offset += gantry_tilt * cos(gantry)
            gplane_offset += gantry_sag * sin(gantry)
            if align_to_pixels:
                long_offset = pixel_align(sim_single.pixel_size, long_offset)
                gplane_offset = pixel_align(sim_single.pixel_size, gplane_offset)
            sim_single.add_layer(
                field_layer(
                    field_size_mm=field_size_mm,
                    cax_offset_mm=(long_offset, gplane_offset),
                )
            )
        for offset in bb_offsets:
            if isinstance(offset, dict):
                offset_mm_left = offset["offset_left_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
                offset_mm_up = offset["offset_up_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
                offset_mm_in = -offset["offset_in_mm"] + random.uniform(
                    -jitter_mm, jitter_mm
                )
            else:
                offset_mm_left = offset[0] + random.uniform(-jitter_mm, jitter_mm)
                offset_mm_up = offset[1] + random.uniform(-jitter_mm, jitter_mm)
                offset_mm_in = -offset[2] + random.uniform(
                    -jitter_mm, jitter_mm
                )  # negative because pixels increase as we go out, so to go in we subtract

            long_offset = bb_projection_long(
                offset_in=offset_mm_in,
                offset_up=offset_mm_up,
                offset_left=offset_mm_left,
                sad=1000,
                gantry=gantry,
            )
            gplane_offset = bb_projection_gantry_plane(
                offset_left=offset_mm_left,
                offset_up=offset_mm_up,
                sad=1000,
                gantry=gantry,
            )
            if align_to_pixels:
                long_offset = pixel_align(sim_single.pixel_size, long_offset)
                gplane_offset = pixel_align(sim_single.pixel_size, gplane_offset)
            sim_single.add_layer(
                PerfectBBLayer(
                    cax_offset_mm=(
                        long_offset,
                        gplane_offset,
                    ),
                    bb_size_mm=bb_size_mm,
                )
            )
        if final_layers is not None:
            for layer in final_layers:
                sim_single.add_layer(layer)
        file_name = f"WL G={gantry}, C={coll}, P={couch}; Field={field_size_mm}mm (shifts={field_offsets}); BB={bb_size_mm}mm @ left={offset_mm_left:.2f}, in={offset_mm_in:.2f}, up={offset_mm_up:.2f}; Gantry tilt={gantry_tilt}, Gantry sag={gantry_sag}.dcm"
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
    cone_layer: type[FilterFreeConeLayer] | type[PerfectConeLayer],
    dir_out: str,
    cone_size_mm: float = 17.5,
    final_layers: list[Layer] | None = None,
    bb_size_mm: float = 5,
    offset_mm_left: float = 0,
    offset_mm_up: float = 0,
    offset_mm_in: float = 0,
    image_axes: ((int, int, int), ...) = (
        (0, 0, 0),
        (90, 0, 0),
        (180, 0, 0),
        (270, 0, 0),
    ),
    gantry_tilt: float = 0,
    gantry_sag: float = 0,
    clean_dir: bool = True,
) -> list[str]:
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
        The tilt of the gantry in degrees that affects the position at 0 and 180. Simulates a simple cosine function.
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
    for gantry, coll, couch in image_axes:
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


def pixel_align(pixel_size: float, length_mm: float) -> float:
    """Due to the finite pixel size, a desired shift may not be possible
    at the exact distance desired. This may cause benchmarking issues when
    the user thinks the field or bb is at X when really it's at X +/- 1/2 * pixel size

    This corrects the value to be aligned to the nearest pixel"""
    return round(length_mm / pixel_size) * pixel_size
