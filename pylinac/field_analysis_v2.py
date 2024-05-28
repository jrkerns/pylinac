from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from .core import image
from .core.image import ImageLike
from .core.profile import (
    FWXMProfilePhysical,
    HillProfilePhysical,
    InflectionDerivativeProfilePhysical,
)
from .core.utilities import ResultBase, ResultsDataMixin
from .metrics.profile import PenumbraLeftMetric, PenumbraRightMetric, TopDistanceMetric


class FieldResult(ResultBase):
    data: str


DEFAULT_METRICS = (
    PenumbraLeftMetric(),
    PenumbraRightMetric(),
    TopDistanceMetric(),
)


class FieldAnalysis(ResultsDataMixin[FieldResult]):
    x_profile: FWXMProfilePhysical | HillProfilePhysical | InflectionDerivativeProfilePhysical
    y_profile: FWXMProfilePhysical | HillProfilePhysical | InflectionDerivativeProfilePhysical

    def __init__(self, path: str | Path, **kwargs):
        self.image: ImageLike = image.load(path, **kwargs)
        self.image.check_inversion_by_histogram()

    def analyze(
        self,
        position: tuple | str = (0.5, 0.5),
        width: float = 0.05,
        edge_type: str = "FWXM",
        invert: bool = False,
        metrics=DEFAULT_METRICS,
    ):
        if invert:
            self.image.invert()
        # determine the position to extract

        top = self.image.shape[0] * position[0] - self.image.shape[0] * width / 2
        bottom = self.image.shape[0] * position[0] + self.image.shape[0] * width / 2
        x_values = self.image[round(top) : round(bottom), :].mean(axis=0)

        if edge_type == "FWXM":
            self.x_profile = FWXMProfilePhysical(values=x_values, dpmm=self.image.dpmm)
        self.x_profile.compute(metrics=metrics)

    def plot_analyzed_images(self):
        fig, ax = plt.subplots()
        self.x_profile.plot(axis=ax, show=False)
        fig, ax = plt.subplots()
        self.image.plot(ax=ax, show=False)
        plt.show()
