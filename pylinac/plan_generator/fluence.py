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

    # number of beams
    num_beams = len(rt_plan.BeamSequence)
    if num_beams == 0:
        return np.empty(0)

    # y-axis (height)
    leaf_boundaries_elem = [
        elem for elem in rt_plan.iterall() if elem.tag == "LeafPositionBoundaries"
    ]
    leaf_boundaries = np.array([(x[0], x[-1]) for x in leaf_boundaries_elem])
    y = np.arange(
        np.min(leaf_boundaries), np.max(leaf_boundaries) + resolution_mm, resolution_mm
    )

    # x-axis (width)
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
        meterset_per_cp = np.diff(cumulative_meterset, prepend=0)

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

            # this pre-allocation is unnecessary since DICOM standard mandates the leaves to be defined on the first control point, but pycharm will issue a warning otherwise.
            leaves_b = leaves_a = np.zeros(number_of_leaf_pairs)
            for cp_idx, cp in enumerate(beam.ControlPointSequence):
                beam_limiting_device_position_sequence = cp.get(
                    "BeamLimitingDevicePositionSequence"
                )
                if cp_idx == 0 or beam_limiting_device_position_sequence is not None:
                    # update mlc positions, otherwise the previous positions will be used
                    leaf_positions = [
                        bld.LeafJawPositions
                        for bld in beam_limiting_device_position_sequence
                        if bld.RTBeamLimitingDeviceType == mlc_id
                    ]
                    leaves_b = np.array(leaf_positions)[0, 0:number_of_leaf_pairs]
                    leaves_a = np.array(leaf_positions)[0, number_of_leaf_pairs:]
                mu = meterset_per_cp[cp_idx]
                mask = (x > leaves_b[np.newaxis].T) & (x <= leaves_a[np.newaxis].T)
                stack_fluence_compact[mask] += mu

            boundaries = [
                bld.LeafPositionBoundaries
                for bld in beam.BeamLimitingDeviceSequence
                if bld.RTBeamLimitingDeviceType == mlc_id
            ]
            row_to_leaf_map = np.argmax(np.array(boundaries).T - y > 0, axis=0) - 1
            for row in range(len(y)):
                leaf = row_to_leaf_map[row]
                if leaf < 0:
                    continue
                stack_fluences[stack_idx, row, :] = stack_fluence_compact[leaf, :]

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
    if len(fluences) == 0:
        return []
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
