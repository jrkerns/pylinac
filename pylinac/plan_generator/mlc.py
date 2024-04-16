import numpy as np


class MLCShaper:
    """The MLC Shaper is a tool for generating MLC positions and sequences to create a given pattern.
    It can also create 'sacrifices' of MLC leaves to set the MLC or gantry speed to a certain value.
    """

    control_points: list[list[float]]

    def __init__(
        self, leaf_y_positions: list[float], max_x: float, sacrifice_gap: int = 5
    ):
        """
        Parameters
        ----------
        leaf_y_positions
            The y-positions of the MLC leaves. This is the same as the LeafJawPositions in the DICOM RT Plan.
        max_x
            The maximum x-position of the MLC leaves.
        sacrifice_gap
            The gap between the sacrificial leaves. This is used to separate the leaves that are being moved out of the way.
        """
        self.leaf_y_positions = leaf_y_positions
        self.max_x = max_x  # mm
        self.sacrifice_gap = sacrifice_gap  # mm gap
        self.control_points = []

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
        return len(self.leaf_y_positions) - 1

    def as_control_points(self) -> list[list[float]]:
        """Return the MLC positions in DICOM format as a list of positions for each control point"""
        return self.control_points

    def add_rectangle(
        self,
        left_position: float,
        right_position: float,
        x_outfield_position: float,
        top_position: float,
        bottom_position: float,
        outer_strip_width: float,
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
        """
        positions: list = [0] * self.num_leaves * 2
        left_pos, right_pos, x_out, y_up, y_lo = (
            left_position,
            right_position,
            x_outfield_position,
            top_position,
            bottom_position,
        )
        for idx, leaf_center in enumerate(self.centers):
            positions[idx] = left_pos if y_lo < leaf_center < y_up else x_out
            positions[idx + self.num_leaves] = (
                right_pos if y_lo < leaf_center < y_up else x_out
            )
            if (leaf_center > y_up) or (leaf_center < y_lo):
                positions[idx] -= outer_strip_width / 2
                positions[idx + self.num_leaves] += outer_strip_width / 2
        self.control_points.append(positions)

    def add_strip(
        self,
        position: float,
        strip_width: float,
    ) -> None:
        """Create a single strip composed of MLCs.
        This is a subset of the `add_rectangle` method, but centers the strip about the x_infield_position and uses
        all the leaves.

        Parameters
        ----------
        position
            The central x-position of the leaves for the leaves on the 'infield' in mm.
        strip_width
            The width of the strip in mm, centered about the x_infield_position
        """
        self.add_rectangle(
            left_position=position - strip_width / 2,
            right_position=position + strip_width / 2,
            x_outfield_position=-200,  # not relevant/used since we will use all the leaves
            top_position=max(self.leaf_y_positions),
            bottom_position=min(self.leaf_y_positions),
            outer_strip_width=1,  # not relevant/used since we will use all the leaves
        )
