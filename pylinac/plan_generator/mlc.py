from __future__ import annotations

import numpy as np

from pylinac.core import validators


class MLCShaper:
    """The MLC Shaper is a tool for generating MLC positions and sequences to create a given pattern.
    Meterset values can be given and set.
    It can also create 'sacrifices' of MLC leaves to set the MLC speed/dose rate to a certain value.
    """

    control_points: list[list[float]]
    metersets: list[float]

    def __init__(
        self,
        leaf_y_positions: list[float],
        max_x_mm: float,
        sacrifice_gap_mm: float = 5,
        sacrifice_max_move_mm: float = 50,
        max_overtravel_mm: float = 140,
    ):
        """
        Parameters
        ----------
        leaf_y_positions
            The y-positions of the MLC leaves. This is the same as the LeafJawPositions in the DICOM RT Plan.
        max_x_mm
            The maximum x-position of the MLC leaves. E.g. 200mm away from the isocenter is 200.
        sacrifice_gap_mm
            The gap between the sacrificial leaves. This is used to separate the leaves that are being moved out of the way.
        sacrifice_max_move_mm
            The maximum distance a sacrificial leaf can move in one control point.
        max_overtravel_mm
            The maximum distance a leaf can move beyond another MLC leaf and also the limit of exposure of the MLC tail.
        """
        self.leaf_y_positions = leaf_y_positions
        self.max_x = max_x_mm  # mm
        self.sacrifice_gap = sacrifice_gap_mm  # mm gap
        self.sacrifice_max_move_mm = sacrifice_max_move_mm  # mm
        self.max_overtravel_mm = max_overtravel_mm  # mm
        self.control_points = []
        self.metersets = []

    @property
    def centers(self) -> list[float]:
        """The center positions of the MLC leaves"""
        centers = []
        for leaf_start, leaf_end in zip(
            self.leaf_y_positions[:-1], self.leaf_y_positions[1:]
        ):
            centers.append(float(np.mean([leaf_start, leaf_end])))
        return centers

    @property
    def num_leaves(self) -> int:
        """The number of leaves in the MLC"""
        return int((len(self.leaf_y_positions) - 1) * 2)

    @property
    def num_pairs(self) -> int:
        """The number of leaf pairs in the MLC"""
        return int(self.num_leaves / 2)

    def as_control_points(self) -> list[list[float]]:
        """Return the MLC positions in DICOM format as a list of positions for each control point"""
        return self.control_points

    def as_metersets(self) -> list[float]:
        """Return the MLC metersets in DICOM format as a list for each control point"""
        return self.metersets

    def add_rectangle(
        self,
        left_position: float,
        right_position: float,
        x_outfield_position: float,
        top_position: float,
        bottom_position: float,
        outer_strip_width: float,
        meterset_at_target: float,
        meterset_transition: float = 0,
        sacrificial_distance: float = 0,
        initial_sacrificial_gap: float | None = None,
    ) -> None:
        """Create a rectangle using the MLCs.

        Parameters
        ----------
        left_position
            The left positions the MLCs should have when they are on the "infield" area in mm
        right_position
            The right side of the rectangle in mm.
        x_outfield_position
            The position the MLCs should have when they are on the "outfield" area. Typically, these are unused/outside leaves and just need to be put somewhere.
        top_position
            The upper y-bound that defines the out/in boundary in mm.
        bottom_position
            The lower y-bound that defines the out/in boundary in mm.
        outer_strip_width
            The separation width in mm of the leaves for the leaves outside the rectangle.
        meterset_at_target
            The ratio of MU that should be delivered AT the target rectangle position.
            This is for delivered rectangles of dose. E.g. for a picket fence the
            meterset at target might be 0.25 to deliver 25% of the dose at each picket.
        meterset_transition
            The ratio of MU that should be delivered between the MLC state before the rectangle and the MLC state AT the target rectangle position.
            Set to 0 to transition immediately without dose. E.g. if delivering a picket fence
            the dose between pickets is low or 0. For an MLC speed test or an ROI where
            the goal is to deliver through a transition region, this might be high compared
            to the meterset at target.
        sacrificial_distance
            The distance to move the sacrificial leaves. This is used to module the dose rate or MLC speed.
            If this is set, the meterset_transition must be set to a non-zero value.
        initial_sacrificial_gap
            The initial gap between the sacrificial leaves. This is only used for the first control point.
        """
        positions: list = [0] * self.num_leaves
        left_pos, right_pos, x_out, y_up, y_lo = (
            left_position,
            right_position,
            x_outfield_position,
            top_position,
            bottom_position,
        )
        for idx, leaf_center in enumerate(self.centers):
            positions[idx] = left_pos if y_lo < leaf_center < y_up else x_out
            positions[idx + self.num_pairs] = (
                right_pos if y_lo < leaf_center < y_up else x_out
            )
            if (leaf_center > y_up) or (leaf_center < y_lo):
                positions[idx] -= outer_strip_width / 2
                positions[idx + self.num_pairs] += outer_strip_width / 2
        if initial_sacrificial_gap:
            # add the initial sacrificial gap
            positions[0] -= initial_sacrificial_gap / 2
            positions[int(self.num_leaves / 2) - 1] -= initial_sacrificial_gap / 2
            positions[int(self.num_leaves / 2)] += initial_sacrificial_gap / 2
            positions[-1] += initial_sacrificial_gap / 2
        start_meterset = self.metersets[-1] if self.metersets else 0
        end_meterset = start_meterset + meterset_at_target + meterset_transition
        if end_meterset > 1.0:
            raise ValueError("Meterset exceeds 1.0")
        if sacrificial_distance > 0 and meterset_transition == 0:
            raise ValueError(
                "Sacrificial distance > 0 but transition meterset was 0. Sacrifices are only used in transitions."
            )
        if sacrificial_distance > 0 and initial_sacrificial_gap is not None:
            raise ValueError(
                "Cannot specify both a sacrificial distance and an initial sacrificial gap."
            )
        if initial_sacrificial_gap and len(self.control_points) > 0:
            raise ValueError(
                "Cannot specify an initial sacrificial gap if there are already control points."
            )
        if initial_sacrificial_gap and meterset_transition:
            raise ValueError(
                "Cannot specify an initial sacrificial gap if there is a transition dose."
            )
        # if we have transition doses and sacrificial moves,
        # we might need to interpolate the control points
        if meterset_transition > 0:
            if len(self.control_points) == 0:
                raise ValueError(
                    "Cannot have a transition without a starting control point. Add a control point first."
                )
            if sacrificial_distance > 0:
                # split the sacrifices into chunks in case
                # the distance is too large
                sacrifice_chunks = split_sacrifice_travel(
                    sacrificial_distance, self.sacrifice_max_move_mm
                )
                # calculate the number of interpolation points
                interpolation_ratios = np.cumsum(
                    [m / sum(sacrifice_chunks) for m in sacrifice_chunks]
                )
                interpolated_control_points = interpolate_control_points(
                    control_point_start=self.control_points[-1],
                    control_point_end=positions,
                    interpolation_ratios=interpolation_ratios,
                    sacrifice_chunks=sacrifice_chunks,
                    max_overtravel=self.max_overtravel_mm,
                )
                self.control_points.extend(interpolated_control_points)
                self.metersets.extend(
                    [
                        start_meterset + meterset_transition * ratio
                        for ratio in interpolation_ratios
                    ]
                )
            else:
                # we have transition doses but no sacrifices
                # this just adds a control point
                self.control_points.append(positions)
                self.metersets.append(start_meterset + meterset_transition)
        else:
            # add starting control point; no transition dose
            self.control_points.append(positions)
            self.metersets.append(start_meterset)
            # if there is no dose delivered at the target, we can skip
            # adding another control point
            if end_meterset != start_meterset:
                self.control_points.append(positions)
                self.metersets.append(end_meterset)

    def add_strip(
        self,
        position_mm: float,
        strip_width_mm: float,
        meterset_at_target: float,
        meterset_transition: float = 0,
        sacrificial_distance_mm: float = 0,
        initial_sacrificial_gap_mm: float | None = None,
    ) -> None:
        """Create a single strip composed of MLCs.
        This is a subset of the `add_rectangle` method, but centers the strip about the x_infield_position and uses
        all the leaves.

        Parameters
        ----------
        position_mm
            The central x-position of the leaves for the leaves on the 'infield' in mm.
        strip_width_mm
            The width of the strip in mm, centered about the x_infield_position.
        meterset_at_target
            The ratio of MU that should be delivered within this control point.
            Set to 0 for a "transition" control point, such as at the beginning or end of a beam
            or when moving from one ROI to another.
        meterset_transition
            The ratio of MU that should be delivered between the MLC state before the strip and the MLC state AT the target strip position.
            Set to 0 to transition immediately without dose. E.g. if delivering a picket fence
            the dose between pickets is low or 0. For an MLC speed test or an ROI where
            the goal is to deliver through a transition region, this might be high compared
            to the meterset at target.
        sacrificial_distance_mm
            The distance to move the sacrificial leaves. This is used to module the dose rate.
            If this is set, the meterset_transition must be set to a non-zero value.
        initial_sacrificial_gap_mm
            The initial gap between the sacrificial leaves. This is only used for the first control point.
        """
        self.add_rectangle(
            left_position=position_mm - strip_width_mm / 2,
            right_position=position_mm + strip_width_mm / 2,
            x_outfield_position=-200,  # not relevant/used since we will use all the leaves
            top_position=max(self.leaf_y_positions),
            bottom_position=min(self.leaf_y_positions),
            outer_strip_width=1,  # not relevant/used since we will use all the leaves
            meterset_at_target=meterset_at_target,
            meterset_transition=meterset_transition,
            sacrificial_distance=sacrificial_distance_mm,
            initial_sacrificial_gap=initial_sacrificial_gap_mm,
        )


def next_sacrifice_shift(
    current_position_mm: float,
    travel_mm: float,
    x_width_mm: float,
    other_mlc_position: float,
    max_overtravel_mm: float,
) -> float:
    """Calculate the next position of a sacrificial leaf.
    This will calculate the next position, accounting for the MLC movement range. It
    will try to oscillate the target position.

    Parameters
    ----------
    current_position_mm
        The current position of the leaf.
    x_width_mm
        The width of the MLCs in the x-direction.
    other_mlc_position
        The position of the other leaves. I.e. where are the rest of the leaves generally at.
        This provides a target for the sacrifice to move toward so the leaves stay in generally
        the same area and don't reach overtravel.
    max_overtravel_mm
        The max overtravel allowed by the MLCs.
    """
    largest_travel_allowed = max_overtravel_mm + abs(
        other_mlc_position - current_position_mm
    )
    if travel_mm > largest_travel_allowed:
        raise ValueError("Travel distance exceeds allowed range")
    if x_width_mm < max_overtravel_mm:
        raise ValueError("Max overtravel exceeds MLC width")
    movement_direction = 1 if current_position_mm < other_mlc_position else -1
    target_shift = movement_direction * travel_mm
    # if we go beyond the MLC width limit, go the other way
    if (target_shift + current_position_mm < -x_width_mm / 2) or (
        target_shift + current_position_mm > x_width_mm / 2
    ):
        target_shift = -movement_direction * travel_mm
    return target_shift


def interpolate_control_points(
    control_point_start: list[float],
    control_point_end: list[float],
    interpolation_ratios: list[float],
    sacrifice_chunks: list[float],
    max_overtravel: float,
) -> list[list[float]]:
    """Interpolate between two control points, including the sacrifices needed.
    This interpolates the start and end positions of everything except the first and last pair.
    For those, the sacrifices are injected into the control points per the length of the sacrifice chunks.

    Parameters
    ----------
    control_point_start
        The starting control point.
    control_point_end
        The ending control point.
    interpolation_ratios
        The ratios at which to interpolate between the two control points.
    sacrifice_chunks
        The distances to move the sacrificial leaves for each interpolation ratio.
    max_overtravel
        The maximum overtravel allowed by the MLCs in mm.
    """
    if len(control_point_start) != len(control_point_end):
        raise ValueError("Control points must be the same length")
    if any(
        r < 0 or r > 1.001 for r in interpolation_ratios
    ):  # 1.001 for floating point error
        raise ValueError("Interpolation ratios must be between 0 and 1")
    if len(interpolation_ratios) == 0:
        raise ValueError("Interpolation ratios must be provided")
    if len(interpolation_ratios) != len(sacrifice_chunks):
        raise ValueError(
            "Interpolation ratios must be the same length as the sacrifice chunks"
        )
    num_leaves = int(len(control_point_start) / 2)
    all_cps = [control_point_start]
    for idx, (ratio, sacrifice) in enumerate(
        zip(interpolation_ratios, sacrifice_chunks)
    ):
        last_cp = all_cps[-1]
        sacrificial_shift = next_sacrifice_shift(
            current_position_mm=last_cp[0],
            travel_mm=sacrifice,
            x_width_mm=400,
            other_mlc_position=last_cp[1],
            max_overtravel_mm=max_overtravel,
        )
        new_cp = [
            start + (end - start) * ratio
            for start, end in zip(control_point_start, control_point_end)
        ]
        # set the sacrifical leaf positions
        new_cp[0] = last_cp[0] + sacrificial_shift
        new_cp[num_leaves - 1] = last_cp[num_leaves - 1] + sacrificial_shift
        new_cp[num_leaves] = last_cp[num_leaves] + sacrificial_shift
        new_cp[-1] = last_cp[-1] + sacrificial_shift
        all_cps.append(new_cp)
    return all_cps[1:]


def split_sacrifice_travel(distance: float, max_travel: float) -> list[float]:
    """
    Split a number into a list of multiples of a max travel distance and a remainder. E.g. 66 => [50, 16].

    Parameters
    ----------
    distance : float
        The distance to split.
    max_travel : float
        The maximum travel distance allowed.

    Returns
    -------
    list
        A list containing multiples of the max travel and the remainder.
    """
    validators.is_positive(distance)
    validators.is_positive(max_travel)
    result = []
    while distance >= max_travel:
        result.append(max_travel)
        distance -= max_travel
    if distance > 0:
        result.append(distance)
    return result
