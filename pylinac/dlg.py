from math import floor
from typing import Iterable

from pylinac.core import image
from pylinac.picketfence import MLCs, MLCArrangement

path = r'E:\OneDrive - Falcon Ridge Software\Programming\Python\Projects\pylinac\tests_basic\test_files\DLG 1.5+0.2.dcm'


class DLG:

    def __init__(self, path):
        self.image = image.LinacDicomImage(path)

    def analyze(self, gaps: Iterable, mlc: MLCArrangement, y_field_size: float = 100):
        profiles = []
        padding = 20
        padding_px = self.image.dpmm * padding
        mid_point = self.image.shape[0]/2
        for idx, center in enumerate(mlc.centers):
            if -y_field_size/2 < center < y_field_size/2:
                top = floor(center + mlc.widths[idx]/4)
                bottom = floor(center - mlc.widths[idx]/4)
                window = self.image[bottom:top, mid_point - padding_px:mid_point+padding_px]
                profiles.append(window.mean(axis=0))



# dlg = DLG(path)
# dlg.analyze(gaps=(-0.9, -1.1, -1.3, -1.5, -1.7, -1.9), mlc=MLCs['Millennium'])
