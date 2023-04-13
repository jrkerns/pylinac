import os.path
import tempfile
import unittest
from pathlib import Path

from pylinac import (
    DRGS,
    CatPhan604,
    Dynalog,
    FieldAnalysis,
    LeedsTOR,
    PicketFence,
    QuartDVT,
    StandardImagingFC2,
    Starshot,
    WinstonLutz,
)

CUSTOM_LOGO = Path(__file__).parent / "rm-logo.png"  # stuffcicles


class TestPDF(unittest.TestCase):
    def test_normal_pdf(self):
        pf = PicketFence.from_demo_image()
        pf.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            pf.publish_pdf(tf)

    def test_custom_logo_str(self):
        pf = PicketFence.from_demo_image()
        pf.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            pf.publish_pdf(
                tf, logo=os.path.join(os.path.dirname(__file__), "rm-logo.png")
            )

    def test_custom_logo_path(self):
        pf = PicketFence.from_demo_image()
        pf.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            pf.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_catphan(self):
        ct = CatPhan604.from_demo_images()
        ct.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            ct.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_field_analysis(self):
        fa = FieldAnalysis.from_demo_image()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_dynalog(self):
        dyn = Dynalog.from_demo()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            dyn.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_planar_imaging(self):
        fa = LeedsTOR.from_demo_image()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_fc2(self):
        fa = StandardImagingFC2.from_demo_image()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_quart(self):
        fa = QuartDVT.from_demo_images()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_starshot(self):
        fa = Starshot.from_demo_image()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_vmat(self):
        fa = DRGS.from_demo_images()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)

    def test_winston_lutz(self):
        fa = WinstonLutz.from_demo_images()
        fa.analyze()
        with tempfile.NamedTemporaryFile() as tf:
            # shouldn't raise.
            fa.publish_pdf(tf, logo=CUSTOM_LOGO)
