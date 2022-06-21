import io
from unittest import TestCase

from pylinac.core.image import load
from pylinac.core.image_generator import AS1200Image, PerfectFieldLayer, PerfectBBLayer, PerfectConeLayer, AS500Image, \
    AS1000Image, ConstantLayer
from pylinac.core.image_generator.simulators import Simulator
from pylinac.core.profile import SingleProfile


def profiles_from_simulator(simulator: Simulator) -> (SingleProfile, SingleProfile):
    stream = io.BytesIO()
    simulator.generate_dicom(stream)
    stream.seek(0)
    img = load(stream)
    inplane_profile = SingleProfile(img[:, int(simulator.shape[1]/2)], dpmm=img.dpmm)
    cross_profile = SingleProfile(img[int(simulator.shape[0]/2), :], dpmm=img.dpmm)
    return inplane_profile, cross_profile


class TestPerfectFieldLayer(TestCase):

    def test_10x10_100sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(100, 100)))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 100, delta=0.1)
            self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 100, delta=0.1)

    def test_10x10_150sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(PerfectFieldLayer(field_size_mm=(100, 100)))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 150, delta=0.2)
            self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 150, delta=0.2)


class TestPerfectConeLayer(TestCase):

    def test_15mm_100sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1000)
            as1200.add_layer(PerfectConeLayer(cone_size_mm=15))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 15, delta=1)
            self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 15, delta=1)

    def test_15mm_150sid(self):
        for sim in (AS500Image, AS1000Image, AS1200Image):
            as1200 = sim(sid=1500)
            as1200.add_layer(PerfectConeLayer(cone_size_mm=15))
            inplane_profile, cross_profile = profiles_from_simulator(as1200)
            self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 15*1.5, delta=1)
            self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 15*1.5, delta=1)


class TestPerfectBBLayer(TestCase):

    def test_10mm_100sid(self):
        as1200 = AS1200Image(sid=1000)
        as1200.add_layer(ConstantLayer(constant=1.0))  # simulate huge field for easier analysis later on
        as1200.add_layer(PerfectBBLayer(bb_size_mm=10))
        stream = io.BytesIO()
        as1200.generate_dicom(stream)
        stream.seek(0)
        img = load(stream)
        img.invert()  # we invert so the BB looks like a profile, not a dip
        inplane_profile = SingleProfile(img[:, int(as1200.shape[1] / 2)], dpmm=img.dpmm)
        cross_profile = SingleProfile(img[int(as1200.shape[0] / 2), :], dpmm=img.dpmm)
        self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 10, delta=1)
        self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 10, delta=1)

    def test_10mm_150sid(self):
        as1200 = AS1200Image(sid=1500)
        as1200.add_layer(ConstantLayer(constant=1.0))  # simulate huge field for easier analysis later on
        as1200.add_layer(PerfectBBLayer(bb_size_mm=10))
        stream = io.BytesIO()
        as1200.generate_dicom(stream)
        stream.seek(0)
        img = load(stream)
        img.invert()  # we invert so the BB looks like a profile, not a dip
        inplane_profile = SingleProfile(img[:, int(as1200.shape[1] / 2)], dpmm=img.dpmm)
        cross_profile = SingleProfile(img[int(as1200.shape[0] / 2), :], dpmm=img.dpmm)
        self.assertAlmostEqual(inplane_profile.fwxm_data()['width (exact) mm'], 15, delta=1)
        self.assertAlmostEqual(cross_profile.fwxm_data()['width (exact) mm'], 15, delta=1)
