"""
Tests the pre-fabricated plans generated by the prefab_plan_generator.py script
and listed on the Pylinac docs site for use by the public.
"""

import os
import re
from pathlib import Path
from unittest import TestCase

import pydicom

from tests_basic.utils import get_folder_from_cloud_repo


class TestGeneratedPlans(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # if we're in CI/CD, download from the GCP bucket
        if os.environ.get("CI"):
            cls.cloud_folder = get_folder_from_cloud_repo(
                folder=["rtplans"], cloud_repo="pylinac_demo_files"
            )
        # if we're local, use the locally generated plans
        else:
            cls.cloud_folder = Path(__file__).parent.parent / "scripts"

    def test_num_rev2_plans(self):
        self.assertEqual(len(list(Path(self.cloud_folder).glob("R2*"))), 12)

    def test_all_energies_of_plan_are_the_same(self):
        for plan in Path(self.cloud_folder).glob("R2*"):
            energy_str = plan.name.split("_")[1]
            energy = int(re.match(r"(\d+)(MV|MVFFF)?", energy_str).group(1))
            ds = pydicom.dcmread(plan)
            for beam in ds.BeamSequence:
                self.assertEqual(
                    int(beam.ControlPointSequence[0].NominalBeamEnergy),
                    energy,
                    msg=f"Plan {plan} has a beam {beam.BeamName} with a different energy than the plan",
                )

    def test_fluence_mode_of_plan_are_the_same(self):
        FLUENCE_MAP = {"NON_STANDARD": True, "STANDARD": False}
        for plan in Path(self.cloud_folder).glob("R2*"):
            is_fff_beam = "FFF" in plan.name
            ds = pydicom.dcmread(plan)
            for beam in ds.BeamSequence:
                fluence = beam.PrimaryFluenceModeSequence[0].FluenceMode
                self.assertEqual(
                    FLUENCE_MAP[fluence],
                    is_fff_beam,
                    msg=f"Plan {plan} has a beam {beam.BeamName} with a different MLC than the plan",
                )

    def test_nominal_dose_rate_is_the_same_per_plan(self):
        for plan in Path(self.cloud_folder).glob("R2*"):
            ds = pydicom.dcmread(plan)
            first_dose_rate = int(
                ds.BeamSequence[0].ControlPointSequence[0].DoseRateSet
            )
            for beam in ds.BeamSequence[1:]:
                self.assertEqual(
                    int(beam.ControlPointSequence[0].DoseRateSet),
                    first_dose_rate,
                    msg=f"Plan {plan} has a beam {beam.BeamName} with a different dose rate than the plan",
                )
