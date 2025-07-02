import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from py_linq import Enumerable
from pydicom import Dataset


def get_scaled_leaf_boundaries(boundaries: list[float], scale: float) -> np.ndarray:
    """Get the MLC boundaries scaled to the resolution.

    Parameters
    ----------
    boundaries: list[float]
        The MLC boundaries.
    scale : float
        The resolution scale.

    Returns
    -------
    np.ndarray
        The MLC boundaries scaled to the resolution.
    """
    return np.asarray(boundaries, dtype=float) * scale


def get_mlc_y_extent(ds: Dataset, scale: float) -> int:
    """Get the **fluence** extent of the MLC boundaries in the y direction. I.e. the height of the fluence map.
    Can account for multiple MLC stacks. Assumes the first beam contains the MLCs.
    E.g. if the MLCs have a range of -100 to 100 and scale of 2, the extent will be 200*2=400.
    """
    extent = 0
    for device in ds.BeamSequence[0].BeamLimitingDeviceSequence:
        if "MLC" in device.get("RTBeamLimitingDeviceType"):
            mlc_bounds = device.get("LeafPositionBoundaries")
            if mlc_bounds:
                extent = max(
                    extent, np.ptp(get_scaled_leaf_boundaries(mlc_bounds, scale))
                )
    return int(extent)


def get_absolute_scaled_leaf_boundaries(
    ds: Dataset, boundaries: list[float], scale: float
) -> np.ndarray:
    """Get the MLC boundaries scaled to the resolution and offset to start at 0. Used to
    get the boundaries of the MLCs in the fluence map (vs physical which might be negative).
    E.g. if the MLC boundaries go (-100, -90, ... 100), the result will be (0, 20, ... 400) if the scale is 2.

    Parameters
    ----------
    ds : pydicom.Dataset
        The RT Plan dataset.
    scale : float
        The resolution scale.

    Returns
    -------
    np.ndarray
        The MLC boundaries scaled to the resolution and offset by the given amount.
    """
    extent = get_mlc_y_extent(ds, scale)
    return (get_scaled_leaf_boundaries(boundaries, scale) + extent / 2).astype(int)


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
    num_beams = len(rt_plan.BeamSequence)
    extent = get_mlc_y_extent(rt_plan, z)
    fluences = np.zeros(
        (num_beams, extent, int(width_mm * z)),
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
        # For each MLC stack stated in the Beam limiting sequence, construct the fluence
        mlc_stacks = (
            Enumerable(beam.BeamLimitingDeviceSequence)
            .where(lambda c: "MLC" in c.RTBeamLimitingDeviceType)
            .to_list()
        )
        stack_fluences = []
        for stack in mlc_stacks:
            running_meterset = 0
            stack_fluence = np.zeros_like(fluences[0])
            for cp_idx, cp in enumerate(beam.ControlPointSequence):
                # for each MLC stack in the CP
                mlc_cp = (
                    Enumerable(cp.BeamLimitingDevicePositionSequence)
                    .where(
                        lambda cp: cp.RTBeamLimitingDeviceType
                        == stack.RTBeamLimitingDeviceType
                    )
                    .single()
                )
                mid_leaf_index = len(mlc_cp.LeafJawPositions) // 2
                left_leaves = (
                    np.asarray(mlc_cp.LeafJawPositions[:mid_leaf_index]) * z
                    + width_mm * z / 2
                ).astype(int)
                left_leaves[left_leaves < 0] = 0
                right_leaves = (
                    np.asarray(mlc_cp.LeafJawPositions[mid_leaf_index:]) * z
                    + width_mm * z / 2
                ).astype(int)
                right_leaves[right_leaves < 0] = 0
                mu_delta = int((cp.CumulativeMetersetWeight - running_meterset) * 1000)
                boundaries = get_absolute_scaled_leaf_boundaries(
                    rt_plan, stack.LeafPositionBoundaries, z
                )
                for leaf_num in range(len(left_leaves)):
                    width_slice = slice(left_leaves[leaf_num], right_leaves[leaf_num])
                    height_slice = slice(boundaries[leaf_num], boundaries[leaf_num + 1])
                    stack_fluence[height_slice, width_slice] += mu_delta
                running_meterset = cp.CumulativeMetersetWeight
            stack_fluences.append(stack_fluence)
        if len(stack_fluences) == 1:
            fluences[beam_idx] += stack_fluences[0]
        else:
            # for multiple MLC stacks, take the minimum fluence across each one.
            fluences[beam_idx] += np.minimum(*stack_fluences)
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
