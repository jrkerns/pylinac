import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from pydicom import Dataset


def generate_fluences(
    rt_plan: Dataset,
    width_mm: float,
    resolution_mm: float = 0.1,
    dtype: np.dtype = np.uint16,
) -> np.ndarray:
    """Generate the fluence map from the RT Plan. This will create an
    MxN array where M is the width in mm * resolution and N is the height set by the bounds of the MLCs * resolution.

    Parameters
    ----------
    rt_plan : pydicom.Dataset
        The RT Plan dataset. Must contain BeamSequence.
    width_mm : int
        The width of the fluence map in mm. Use smaller values for faster calculation.
    resolution_mm : float, optional
        The resolution of the fluence map in mm. Smaller values will take longer to calculate.
    dtype : type, optional
        The data type of the fluence map. Default is uint16.

    Returns
    -------
    np.array
        The fluence map. Will be of shape (num_beams, height, width).
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


def plot_fluences(
    plan: Dataset,
    width_mm: float,
    resolution_mm: float,
    dtype: np.dtype = np.uint16,
    show: bool = True,
) -> list[Figure]:
    """Plot the fluences of the dataset. Generates N figures where N is the number of Beams in the plan BeamSequence.

    Parameters
    ----------
    plan : pydicom.Dataset
        The RT Plan dataset. Must contain BeamSequence.
    width_mm : int
        The width of the fluence map in mm. Use smaller values for faster calculation.
    resolution_mm : float, optional
        The resolution of the fluence map in mm. Smaller values will take longer to calculate.
    dtype : type, optional
        The data type of the fluence map. Default is uint16.
    show : bool, optional
        Whether to show the plots. Default is True.

    Returns
    -------
    list[Figure]
        A list of matplotlib figures, one for each beam in the plan.
    """
    fluences = generate_fluences(plan, width_mm, resolution_mm, dtype)
    m = fluences.max()
    figs = []
    for i, fluence in enumerate(fluences):
        fig, ax = plt.subplots()
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
        figs.append(fig)
    if show:
        plt.show()
    return figs
