"""Tests for the flatsym module of pylinac."""

import unittest
import os.path as osp

import matplotlib.pyplot as plt

from pylinac.core.profile import SingleProfile
from pylinac.flatsym import BeamImage, _is_both_planes, _is_crossplane, _is_inplane


class UserTests(unittest.TestCase):
    """Tests for the user-facing methods."""

    def test_init(self):
        """Test BeamImage creation."""
        # test without file passed
        b = BeamImage()
        # test with file passed
        demo_file = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'flatsym', 'flatsym_demo.dcm')
        b = BeamImage(demo_file)

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        b = BeamImage()
        b.load_demo_image()

    def test_run_demo(self):
        """Test that running the demo works without error."""
        BeamImage().run_demo()

    def test_flatness(self):
        """Test flatness values"""
        img = BeamImage()
        img.load_demo_image()

        self.assertAlmostEqual(img.flatness('x'), 1.91, delta=0.03)
        self.assertAlmostEqual(img.flatness('in'), 1.60, delta=0.03)
        self.assertAlmostEqual(img.flatness('x', method='elekta'), 103.9, delta=0.1)
        self.assertAlmostEqual(img.flatness('both')[0], 1.91, delta=0.03)
        self.assertAlmostEqual(img.flatness('both')[1], 1.60, delta=0.03)
        self.assertRaises(ValueError, img.flatness, method='varen')

    def test_symmetry(self):
        """Test symmetry values"""
        img = BeamImage()
        img.load_demo_image()

        self.assertAlmostEqual(img.symmetry('x'), 3.08, delta=0.05)
        self.assertAlmostEqual(img.symmetry('x', method='elekta'), 103.1, delta=0.5)
        self.assertAlmostEqual(img.symmetry('in'), 2.63, delta=0.02)
        self.assertAlmostEqual(img.symmetry('both')[0], 3.08, delta=0.05)
        self.assertAlmostEqual(img.symmetry('both')[1], 2.63, delta=0.02)
        self.assertRaises(ValueError, img.symmetry, method='varen')

    def test_plot_flatness(self):
        """Test that plots work as expected."""
        img = BeamImage()
        img.load_demo_image()

        img.plot_flatness()
        img.plot_flatness(plane='x')
        img.plot_flatness(plane='in')

    def test_plot_symmeetry(self):
        """Test that plots work as expected."""
        img = BeamImage()
        img.load_demo_image()

        img.plot_symmetry()
        img.plot_symmetry(plane='x')
        img.plot_symmetry(plane='in')

    def test_plot_flatsym(self):
        """Test that plots work as expected."""
        img = BeamImage()
        img.load_demo_image()

        img.plot_flatsym()
        img.plot_flatsym(plane='x')
        img.plot_flatsym(plane='in')

    def test_analyze_wo_img(self):
        b = BeamImage()
        self.assertRaises(AttributeError, b.flatness)


class UnitTests(unittest.TestCase):

    def setUp(self):
        self.img = BeamImage()
        self.img.load_demo_image()

    def test_get_symmetry(self):
        prof = self.img._get_profile('x', 'auto')

        # varian
        for method in ('varian', 'point difference'):
            symmetry, lt_edge, rt_edge, max_idx = self.img._get_symmetry(prof, method)
            self.assertAlmostEqual(symmetry, 3.08, delta=0.01)
            self.assertAlmostEqual(lt_edge, 393, delta=3)
            self.assertAlmostEqual(rt_edge, 681, delta=3)
            self.assertAlmostEqual(max_idx, 429, delta=3)

        # elekta
        for method in ('elekta', 'pdq-iec'):
            symmetry, lt_edge, rt_edge, max_idx = self.img._get_symmetry(prof, method)
            self.assertAlmostEqual(symmetry, 103.11, delta=0.05)

    def test_convert_position(self):

        input_pos = 'auto'

        out_pos = self.img._convert_position(input_pos, 'both')
        self.assertEqual(out_pos[0], 366)
        self.assertEqual(out_pos[1], 537)

        out_pos = self.img._convert_position(input_pos, 'x')
        self.assertEqual(out_pos[0], 366)

        input_pos = (300, 0.6)
        out_pos = self.img._convert_position(input_pos, 'both')
        self.assertEqual(out_pos[0], 300)
        self.assertEqual(out_pos[1], 614)

        input_pos = 0.6
        out_pos = self.img._convert_position(input_pos, 'y')
        self.assertEqual(out_pos[0], 614)

        self.assertRaises(ValueError, self.img._convert_position, (300, 400, 0.5), 'both')
        self.assertRaises(ValueError, self.img._convert_position, (300, 400), 'x')
        self.assertRaises(ValueError, self.img._convert_position, (300, 400), 'notaplane')

    def test_parse_position(self):

        position = 30
        int_input = self.img._parse_position(position, 'x')
        self.assertEqual(int_input, position)

        float_pos = 0.4
        fl_input = self.img._parse_position(float_pos, 'i')
        self.assertAlmostEqual(fl_input, float_pos*self.img.array.shape[1], delta=1)
        fl_input = self.img._parse_position(float_pos, 'x')
        self.assertAlmostEqual(fl_input, float_pos * self.img.array.shape[0], delta=1)

        self.assertRaises(ValueError, self.img._parse_position, position, 'both')
        self.assertRaises(ValueError, self.img._parse_position, 1.1, 'x')

    def test_get_profile(self):
        """Test various inputs for _get_profile()."""
        b = BeamImage()
        self.assertRaises(AttributeError, b._get_profile, 'x', 'auto')

        # test output type
        self.assertIsInstance(self.img._get_profile('x', 'auto'), SingleProfile)

        # test axes
        for axis in ('crossplane', 'cross', 'x', 'inplane', 'in'):
            self.img._get_profile(axis, 'auto')

        # self.assertRaises(ValueError, self.img._get_profile, 'notaplane', 'auto')

        # make sure extraction is along correct axis
        size_xplane = 1024
        prof = self.img._get_profile('x', 'auto')
        self.assertEqual(prof.y_values.size, size_xplane)

        size_yplane = 768
        prof = self.img._get_profile('in', 'auto')
        self.assertEqual(prof.y_values.size, size_yplane)

    def test_check_position_bounds(self):
        # test positions
        for pos in (100, 500, 0.5):
            self.img._check_position_inbounds(pos, 'x')
        self.assertRaises(IndexError, self.img._check_position_inbounds, 1100, 'x')
        self.assertRaises(IndexError, self.img._check_position_inbounds, 800, 'i')

    def test_get_flatness(self):
        prof = self.img._get_profile('x', 'auto')

        # varian
        for method in ('varian', 'VoM80', 'siemens'):
            flatness, dmax, dmin, lt_edge, rt_edge = self.img._get_flatness(prof, method)
            self.assertAlmostEqual(flatness, 1.91, delta=0.01)
            self.assertAlmostEqual(dmax, 58533, delta=20)
            self.assertAlmostEqual(dmin, 56335, delta=20)
            self.assertAlmostEqual(lt_edge, 393, delta=3)
            self.assertAlmostEqual(rt_edge, 681, delta=3)

        # elekta
        for method in ('elekta', 'IEC'):
            flatness, dmax, dmin, lt_edge, rt_edge = self.img._get_flatness(prof, method)
            self.assertAlmostEqual(flatness, 103.9, delta=0.05)

    def test_plot_annotation(self):
        fig, ax = plt.subplots()

        # axis = self.img._plot_annotation(ax, )

    def test_plot_title(self):

        planes = ('x', 'in')
        prefixs = ('Crossplane', 'Inplane')

        flatsym = ('flat', 'sym')
        suffixs = (' Flatness', ' Symmetry')

        for plane, prefix in zip(planes, prefixs):
            for fs, suffix in zip(flatsym, suffixs):
                fig, ax = plt.subplots()
                title = prefix + suffix
                axis = self.img._plot_title(ax, plane, fs)
                self.assertEqual(title, axis.get_title())

    def test_determine_center(self):
        """Test the determined center is correct"""
        y, x = self.img._determine_center('both')
        self.assertAlmostEqual(y, 366, delta=4)
        self.assertAlmostEqual(x, 537, delta=4)

    def test_check_inversion(self):
        """Test that inverting the image doesn't change the flat/sym values."""
        sym = self.img.symmetry('x')
        flat = self.img.flatness('in')
        self.img.invert()
        i_sym = self.img.symmetry('x')
        i_flat = self.img.flatness('in')
        self.assertAlmostEqual(sym, i_sym, delta=0.01)
        self.assertAlmostEqual(flat, i_flat, delta=0.01)

    def test_img_is_loaded(self):
        """Test image is loaded property."""
        b = BeamImage()
        self.assertFalse(b._img_is_loaded)
        b.load_demo_image()
        self.assertTrue(b._img_is_loaded)


class TestMisc(unittest.TestCase):

    def test_is_crossplane(self):
        self.assertTrue(_is_crossplane('crossplane'))
        self.assertTrue(_is_crossplane('x'))
        self.assertFalse(_is_crossplane('inplane'))
        self.assertFalse(_is_crossplane('in'))

    def test_is_inplane(self):
        self.assertTrue(_is_inplane('inplane'))
        self.assertTrue(_is_inplane('in'))
        self.assertTrue(_is_inplane('y'))
        self.assertFalse(_is_inplane('x'))
        self.assertFalse(_is_inplane('crossplane'))

    def test_is_both_planes(self):
        self.assertTrue(_is_both_planes('both'))
        self.assertTrue(_is_both_planes('b'))
        self.assertFalse(_is_both_planes('x'))
        self.assertFalse(_is_both_planes('in'))


