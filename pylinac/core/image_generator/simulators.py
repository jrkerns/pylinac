from __future__ import annotations

from abc import ABC

import numpy as np
from plotly import graph_objects as go
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.uid import UID

from ..array_utils import array_to_dicom
from ..plotly_utils import add_title
from .layers import Layer


def generate_file_metadata() -> Dataset:
    file_meta = FileMetaDataset()
    file_meta.TransferSyntaxUID = UID(
        "1.2.840.10008.1.2"
    )  # default DICOM transfer syntax
    return file_meta


class Simulator(ABC):
    """Abstract class for an image simulator"""

    pixel_size: float
    shape: (int, int)
    image: np.ndarray

    def __init__(self, sid: float = 1500):
        """

        Parameters
        ----------
        sid
            Source to image distance in mm.
        """
        self.image = np.zeros(self.shape, np.uint16)
        self.sid = sid
        self.mag_factor = sid / 1000

    def add_layer(self, layer: Layer) -> None:
        """Add a layer to the image"""
        self.image = layer.apply(self.image, self.pixel_size, self.mag_factor)

    def as_dicom(
        self,
        gantry_angle: float = 0.0,
        coll_angle: float = 0.0,
        table_angle: float = 0.0,
        invert_array: bool = False,
        tags: dict | None = None,
    ) -> Dataset:
        """Create and return a pydicom Dataset. I.e. create a pseudo-DICOM image."""
        if invert_array:
            array = -self.image + self.image.max() + self.image.min()
        else:
            array = self.image
        return array_to_dicom(
            array=array,
            sid=self.sid,
            gantry=gantry_angle,
            coll=coll_angle,
            couch=table_angle,
            dpi=25.4 / self.pixel_size,
            extra_tags=tags or {},
        )

    def generate_dicom(self, file_out_name: str, *args, **kwargs) -> None:
        """Save the simulated image to a DICOM file.

        See Also
        --------
        as_dicom
        """
        ds = self.as_dicom(*args, **kwargs)
        ds.save_as(file_out_name, write_like_original=False)

    def plot(self, show: bool = True) -> go.Figure:
        """Plot the simulated image."""
        fig = go.Figure()
        fig.add_heatmap(
            z=self.image,
            colorscale="gray",
            x0=-self.image.shape[1] / 2 * self.pixel_size,
            dx=self.pixel_size,
            y0=-self.image.shape[0] / 2 * self.pixel_size,
            dy=self.pixel_size,
        )
        fig.update_layout(
            yaxis_constrain="domain",
            xaxis_scaleanchor="y",
            xaxis_constrain="domain",
            xaxis_title="Crossplane (mm)",
            yaxis_title="Inplane (mm)",
        )
        add_title(fig, f"Simulated {self.__class__.__name__} @{self.sid}mm SID")
        if show:
            fig.show()
        return fig


class AS500Image(Simulator):
    """Simulates an AS500 EPID image."""

    pixel_size: float = 0.78125
    shape: (int, int) = (384, 512)


class AS1000Image(Simulator):
    """Simulates an AS1000 EPID image."""

    pixel_size: float = 0.390625
    shape: (int, int) = (768, 1024)


class AS1200Image(Simulator):
    """Simulates an AS1200 EPID image."""

    pixel_size: float = 0.336
    shape: (int, int) = (1280, 1280)
