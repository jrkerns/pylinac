"""Per-model CatPhan phantom module definitions.

Each sub-module contains the CTP module classes and the CatPhan model class
for a specific phantom model. The classes import shared infrastructure from
``pylinac.ct`` and build on each other following the physical inheritance of
hardware:

    CP504 (503/504 base)
      └── CP600 (600, introduces CTP763/CTP764)
      └── CP604 (604, introduces CTP732/CTP730/BeadMTF)
      └── CP700 (700, introduces CTP682/CTP714/CTP712)
"""

from . import CP504, CP600, CP604, CP700

__all__ = ["CP504", "CP600", "CP604", "CP700"]
