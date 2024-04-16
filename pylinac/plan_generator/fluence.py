import math
from itertools import zip_longest

import numpy as np
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from pydicom import Dataset


def generate_fluence(
    rt_plan: Dataset, width_mm: int, resolution_mm: float = 0.1, dtype=np.uint16
) -> np.array:
    """Generate the fluence map from the RT Plan. This will create an
    MxN array where M is the width in mm * resolution and N is the height set by the bounds of the MLCs * resolution.
    """

    # zoom factor
    z = 1 / resolution_mm
    # Get the MLC extent
    # if no MLCs, abort
    mlc_bounds_um = (
        rt_plan.BeamSequence[0]
        .BeamLimitingDeviceSequence[-1]
        .get("LeafPositionBoundaries")
    )
    if not mlc_bounds_um:
        return np.zeros((1, 10, 10))
    mlc_bounds_um = np.asarray(mlc_bounds_um, dtype=float) * z
    y_offset_um = np.abs(mlc_bounds_um.min())
    mlc_bounds_um += np.abs(y_offset_um)
    num_beams = len(rt_plan.BeamSequence)
    fluences = np.zeros(
        (num_beams, (int(mlc_bounds_um[-1] - mlc_bounds_um[0])), (int(width_mm * z))),
        dtype=dtype,
    )
    # iterate through each beam
    for beam_idx, beam in enumerate(rt_plan.BeamSequence):
        # if setup field, skip beam
        if beam.TreatmentDeliveryType == "SETUP":
            continue
        # if no MLCs defined, skip beam (I.e. jaw-based open field)
        if (
            len(
                beam.ControlPointSequence[0].get(
                    "BeamLimitingDevicePositionSequence", []
                )
            )
            < 3
        ):
            continue
        running_meterset = 0
        # Fill the fluence matrix
        for cp_idx, cp in enumerate(beam.ControlPointSequence):
            if (
                not cp.get("BeamLimitingDevicePositionSequence")
                and cp.CumulativeMetersetWeight == 1
            ):
                # shortcut for a simple 2-element fluence. the end position is the same as the start
                # position, so we can set this control point to the same as the last one and set the
                # meterset to 1
                cp = beam.ControlPointSequence[cp_idx - 1]
                cp.CumulativeMetersetWeight = 1
            if len(cp.BeamLimitingDevicePositionSequence) == 3:
                pass
                # left_x_jaw = (
                #     cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[0]
                # ) * z
                # right_x_jaw = (
                #     cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[1]
                # ) * z
                # top_y_jaw = (
                #     cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[0] * z
                #     + y_offset_um
                # )
                # bottom_y_jaw = (
                #     cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[1] * z
                #     + y_offset_um
                # )
            else:
                pass
                # left_x_jaw = 0
                # right_x_jaw = 0
                # top_y_jaw = 0
                # bottom_y_jaw = fluences[0].shape[0]
            left_leaves = (
                np.asarray(
                    cp.BeamLimitingDevicePositionSequence[-1].LeafJawPositions[
                        : int(len(mlc_bounds_um) - 1)
                    ]
                )
                * z
                + width_mm * z / 2
            ).astype(int)
            left_leaves[left_leaves < 0] = 0
            right_leaves = (
                np.asarray(
                    cp.BeamLimitingDevicePositionSequence[-1].LeafJawPositions[
                        int(len(mlc_bounds_um) - 1) :
                    ]
                )
                * z
                + width_mm * z / 2
            ).astype(int)
            right_leaves[right_leaves < 0] = 0
            mu_delta = int((cp.CumulativeMetersetWeight - running_meterset) * 1000)
            for idx in range(len(mlc_bounds_um) - 1):
                leaf_bounds_um = mlc_bounds_um[idx : idx + 2].astype(int)
                for px in np.arange(leaf_bounds_um[0], leaf_bounds_um[1]):
                    s = slice(left_leaves[idx], right_leaves[idx])
                    fluences[beam_idx, px, s] += mu_delta
            running_meterset = cp.CumulativeMetersetWeight
    return fluences


def generate_control_point_fluence(
    rt_plan: Dataset,
    beam_idx: int,
    width_mm: int,
    resolution_mm: float = 0.1,
    dtype=np.uint16,
) -> (np.array, list[dict[str, dict[str, int]]]):
    """Generate the fluence map from the RT Plan. This will create an
    MxN array where M is the width in mm * resolution and N is the height set by the bounds of the MLCs * resolution.
    """
    z = 1 / resolution_mm
    # Get the MLC extent
    # if no MLCs, abort
    mlc_bounds_um = (
        rt_plan.BeamSequence[0]
        .BeamLimitingDeviceSequence[-1]
        .get("LeafPositionBoundaries")
    )
    if not mlc_bounds_um:
        return np.zeros((1, 10, 10))
    mlc_bounds_um = np.asarray(mlc_bounds_um, dtype=float) * z
    y_offset_um = np.abs(mlc_bounds_um.min())
    mlc_bounds_um += np.abs(y_offset_um)
    beam = rt_plan.BeamSequence[beam_idx]
    num_cps = len(beam.ControlPointSequence)
    fluences = np.zeros(
        (num_cps, (int(mlc_bounds_um[-1] - mlc_bounds_um[0])), (int(width_mm * z))),
        dtype=dtype,
    )

    mlc_leaves = []
    # iterate through each beam
    # if setup field, skip beam
    if beam.TreatmentDeliveryType == "SETUP":
        return np.zeros((1, 10, 10))
    # if no MLCs defined, skip beam (I.e. jaw-based open field)
    if (
        len(beam.ControlPointSequence[0].get("BeamLimitingDevicePositionSequence", []))
        < 3
    ):
        return np.zeros((1, 10, 10))
    running_meterset = 0
    # Fill the fluence matrix
    for cp_idx, cp in enumerate(beam.ControlPointSequence):
        if (
            not cp.get("BeamLimitingDevicePositionSequence")
            and cp.CumulativeMetersetWeight == 1
        ):
            # shortcut for a simple 2-element fluence. the end position is the same as the start
            # position, so we can set this control point to the same as the last one and set the
            # meterset to 1
            cp = beam.ControlPointSequence[cp_idx - 1]
            cp.CumulativeMetersetWeight = 1
        if len(cp.BeamLimitingDevicePositionSequence) == 3:
            pass
            # left_x_jaw = (
            #     cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[0]
            # ) * z
            # right_x_jaw = (
            #     cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[1]
            # ) * z
            # top_y_jaw = (
            #     cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[0] * z
            #     + y_offset_um
            # )
            # bottom_y_jaw = (
            #     cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[1] * z
            #     + y_offset_um
            # )
        else:
            pass
            # left_x_jaw = 0
            # right_x_jaw = 0
            # top_y_jaw = 0
            # bottom_y_jaw = fluences[0].shape[0]
        left_leaves = (
            np.asarray(
                cp.BeamLimitingDevicePositionSequence[-1].LeafJawPositions[
                    : int(len(mlc_bounds_um) - 1)
                ]
            )
            * z
            + width_mm * z / 2
        ).astype(int)
        left_leaves[left_leaves < 0] = 0
        right_leaves = (
            np.asarray(
                cp.BeamLimitingDevicePositionSequence[-1].LeafJawPositions[
                    int(len(mlc_bounds_um) - 1) :
                ]
            )
            * z
            + width_mm * z / 2
        ).astype(int)
        right_leaves[right_leaves < 0] = 0
        mu_delta = int((cp.CumulativeMetersetWeight - running_meterset) * 1000)
        fluences[cp_idx] = fluences[cp_idx - 1]
        mlc_positions = {}
        for idx in range(len(mlc_bounds_um) - 1):
            # fluence of the very first control point is always 0
            if cp_idx == 0:
                continue
            else:
                leaf_bounds_um = mlc_bounds_um[idx : idx + 2].astype(int)
                for px in np.arange(leaf_bounds_um[0], leaf_bounds_um[1]):
                    s = slice(left_leaves[idx], right_leaves[idx])
                    fluences[cp_idx, px, s] += mu_delta
            mlc_positions[f"A{idx+1}"] = {
                "x0": 0,
                "x1": left_leaves[idx],
                "y0": leaf_bounds_um[0],
                "y1": leaf_bounds_um[1],
            }
            mlc_positions[f"B{idx+1}"] = {
                "x0": right_leaves[idx],
                "x1": fluences.shape[2],
                "y0": leaf_bounds_um[0],
                "y1": leaf_bounds_um[1],
            }
        mlc_leaves.append(mlc_positions)
        running_meterset = cp.CumulativeMetersetWeight
    return fluences, mlc_leaves


def plot_fluences(
    fluences: np.array, plan: Dataset, width_mm: int, resolution_mm: float
):
    """Plot the fluences of the dataset"""
    # figure out the right rows/cols
    if fluences.shape[0] < 3:
        fig, axes = plt.subplots(1, fluences.shape[0])
        if fluences.shape[0] == 1:
            axes = (axes,)
    else:
        fig, axes = plt.subplots(math.ceil(fluences.shape[0] / 3), 3)
        axes = axes.flat
    fig.suptitle(plan.PatientID)
    m = fluences.max()
    for i, (fluence, ax) in enumerate(zip_longest(fluences, axes)):
        if fluence is None or np.isclose(fluence.max(), 0):
            ax.axis("off")
        else:
            ax.imshow(fluence, vmin=0, vmax=m)
            ax.set_title(f"{plan.BeamSequence[i].BeamName}")
            ax.set_xticks([])
            ax.set_yticks([])
        # now plot the jaw positions as vertical/horizontal lines
        beam = plan.BeamSequence[i]
        cp = beam.ControlPointSequence[0]
        scale = 1 / resolution_mm
        x_offset = width_mm * scale / 2
        y_offset = fluence.shape[0] / 2
        left_x_jaw = (
            cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[0] * scale
            + x_offset
        )
        right_x_jaw = (
            cp.BeamLimitingDevicePositionSequence[0].LeafJawPositions[1] * scale
            + x_offset
        )
        top_y_jaw = (
            cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[0] * scale
            + y_offset
        )
        bottom_y_jaw = (
            cp.BeamLimitingDevicePositionSequence[1].LeafJawPositions[1] * scale
            + y_offset
        )
        rect = Rectangle(
            xy=(left_x_jaw, bottom_y_jaw),
            width=right_x_jaw - left_x_jaw,
            height=top_y_jaw - bottom_y_jaw,
            fill=False,
            color="r",
        )
        ax.add_patch(rect)
    plt.show()


def plot_control_point_fluence(
    fluences: np.ndarray, mlc_positions: list[dict[str, dict[str, int]]]
):
    frames = [fluence for fluence in fluences]  # Replace with your actual data
    overall_max = max([frame.max() for frame in frames])
    # Example rectangle definitions for each frame (x0, y0, x1, y1)
    # You would replace these with your actual rectangle coordinates
    # rectangles = [((1+i, 1+i), (3+i, 3+i)) for i in range(len(frames))]
    rectangles = mlc_positions
    # create rectangles that represent the MLCs
    # rectangles = []
    # for frame in frames:

    # Create initial figure
    fig = go.Figure()

    # Add initial image
    fig.add_trace(
        go.Heatmap(
            z=frames[0], colorscale="Viridis", showscale=False, zmin=0, zmax=overall_max
        )
    )

    # Add initial rectangles (shapes) to the layout
    initial_mlc_positions = rectangles[0]
    initial_shapes = []
    for leaf in initial_mlc_positions.values():
        shape = dict(
            type="rect",
            x0=leaf["x0"],
            y0=leaf["y0"],
            x1=leaf["x1"],
            y1=leaf["y1"],
            line=dict(color="Black"),
        )
        initial_shapes.append(shape)
    fig.update_layout(
        shapes=initial_shapes,
        title="Frame 0",
    )

    # Animation
    frames_data = []
    for i, frame in enumerate(frames[1:], start=1):
        # generate all leaf rectangles
        mlcs = rectangles[i]
        leaves = []
        for leaf in mlcs.values():
            shape = dict(
                type="rect",
                x0=leaf["x0"],
                y0=leaf["y0"],
                x1=leaf["x1"],
                y1=leaf["y1"],
                line=dict(color="Black"),
            )
            leaves.append(shape)
        frame_data = go.Frame(
            data=[
                go.Heatmap(
                    z=frame,
                    colorscale="Viridis",
                    showscale=False,
                    zmin=0,
                    zmax=overall_max,
                )
            ],
            name=f"Frame {i}",
            layout=dict(
                shapes=leaves,
                title_text=f"Frame {i}",
            ),
        )
        frames_data.append(frame_data)

    fig.frames = frames_data

    # Animation layout options
    fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "buttons": [
                    {
                        "label": "Play",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "frame": {"duration": 500, "redraw": True},
                                "fromcurrent": True,
                            },
                        ],
                    },
                    {
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": False},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                        "label": "Pause",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ]
    )

    fig.show()
