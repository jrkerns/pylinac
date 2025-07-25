import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
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

    # number of beams
    num_beams = len(rt_plan.BeamSequence)

    # number of rows (height)
    leaf_boundaries_elem = [
        elem for elem in rt_plan.iterall() if elem.tag == "LeafPositionBoundaries"
    ]
    leaf_boundaries = np.array([(x[0], x[-1]) for x in leaf_boundaries_elem])
    y = np.arange(
        np.min(leaf_boundaries), np.max(leaf_boundaries) + resolution_mm, resolution_mm
    )

    # number of cols (width)
    x = np.arange(-width_mm / 2, width_mm / 2 + resolution_mm, resolution_mm)

    fluences = np.zeros((num_beams, len(y), len(x)), dtype=dtype)
    # iterate through each beam
    for beam_idx, beam in enumerate(rt_plan.BeamSequence):
        # if setup field, skip beam
        if beam.TreatmentDeliveryType == "SETUP":
            continue

        cumulative_meterset = 1000 * np.array(
            [float(cps.CumulativeMetersetWeight) for cps in beam.ControlPointSequence]
        )
        cumulative_meterset_per_cp = np.diff(cumulative_meterset, prepend=0)

        mlc_stacks = [
            (blds.RTBeamLimitingDeviceType, blds.NumberOfLeafJawPairs)
            for blds in beam.BeamLimitingDeviceSequence
            if "MLC" in blds.RTBeamLimitingDeviceType
        ]
        stack_fluences = np.zeros((len(mlc_stacks), len(y), len(x)), dtype=dtype)
        for stack_idx, stack in enumerate(mlc_stacks):
            mlc_id = stack[0]
            number_of_leaf_pairs = stack[1]
            stack_fluence_compact = np.zeros((number_of_leaf_pairs, len(x)))

            leaf_positions = np.array(
                [
                    bld.LeafJawPositions
                    for cps in beam.ControlPointSequence
                    for bld in cps.BeamLimitingDevicePositionSequence
                    if bld.RTBeamLimitingDeviceType == mlc_id
                ]
            )
            leaf_positions_b = leaf_positions[:, 0:number_of_leaf_pairs]
            leaf_positions_a = leaf_positions[:, number_of_leaf_pairs:]

            for cp_idx, cp in enumerate(beam.ControlPointSequence):
                leaves_b = leaf_positions_b[cp_idx, :]
                leaves_a = leaf_positions_a[cp_idx, :]
                mu = cumulative_meterset_per_cp[cp_idx]
                mask = (x >= leaves_b[None].T) & (x <= leaves_a[None].T)
                stack_fluence_compact[mask] += mu

            boundaries = [
                bld.LeafPositionBoundaries
                for bld in beam.BeamLimitingDeviceSequence
                if bld.RTBeamLimitingDeviceType == mlc_id
            ]
            row_to_leaf_map = np.argmax(np.array(boundaries).T - y > 0, axis=0) - 1
            for row in range(len(y)):
                stack_fluences[stack_idx, row, :] = stack_fluence_compact[
                    row_to_leaf_map[row], :
                ]

        if len(stack_fluences) == 1:
            fluences[beam_idx] = stack_fluences
        else:
            # for multiple MLC stacks, take the minimum fluence across each one.
            fluences[beam_idx] = np.min(stack_fluences, axis=0)
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
