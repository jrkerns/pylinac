"""Tests for the flatsym module of pylinac."""

import unittest
import os.path as osp

from pylinac.flatsym import BeamImage


class IOTests(unittest.TestCase):

    def test_demo_loads_properly(self):
        """Loading the demo shouldn't raise an error"""
        b = BeamImage()
        b.load_demo_image()

    def test_load(self):
        demo_file = osp.join(osp.dirname(osp.dirname(__file__)), 'pylinac', 'demo_files', 'flatsym', 'flatsym_demo.dcm')
        b = BeamImage(demo_file)


class BeamDemoTest(unittest.TestCase):

    def setUp(self):
        self.img = BeamImage()
        self.img.load_demo_image()

    def test_analyze_wo_img(self):
        b = BeamImage()
        self.assertRaises(AttributeError, b.flatness)

    def test_demo(self):
        self.img.run_demo()
        self.img.run_demo('x')
        self.img.run_demo('i')

    def test_inverted_image(self):
        sym = self.img.symmetry('x')
        self.img.invert()
        i_sym = self.img.symmetry('x')
        self.assertAlmostEqual(sym, i_sym, delta=0.01)

    def test_get_profile(self):
        """Test various inputs for _get_profile()."""
        # test axes
        for axis in ('crossplane', 'cross', 'x', 'inplane', 'in'):
            self.img._get_profile(axis, 'auto')
        self.assertRaises(ValueError, self.img._get_profile, 'crossplain', 'auto')

        # test positions
        for pos in (100, 500, 0.5):
            self.img._get_profile('x', pos)
        self.assertRaises(IndexError, self.img._get_profile, 'x', 2000)
        self.assertRaises(IndexError, self.img._get_profile, 'i', 2000)
        self.assertRaises(ValueError, self.img.flatness, position=(0.5, 1.5))

    def test_flatness(self):
        """Test flatness values"""
        self.assertAlmostEqual(self.img.flatness('x'), 1.91, delta=0.03)
        self.assertAlmostEqual(self.img.flatness('in'), 1.60, delta=0.03)
        self.assertAlmostEqual(self.img.flatness('x', method='elekta'), 103.9, delta=0.1)
        self.assertAlmostEqual(self.img.flatness('both')[0], 1.91, delta=0.03)
        self.assertAlmostEqual(self.img.flatness('both')[1], 1.60, delta=0.03)
        self.assertRaises(ValueError, self.img.flatness, method='varen')

    def test_symmetry(self):
        """Test symmetry values"""
        self.assertAlmostEqual(self.img.symmetry('x'), 3.08, delta=0.05)
        self.assertAlmostEqual(self.img.symmetry('x', method='elekta'), 103.1, delta=0.5)
        self.assertAlmostEqual(self.img.symmetry('in'), 2.63, delta=0.02)
        self.assertAlmostEqual(self.img.symmetry('both')[0], 3.08, delta=0.05)
        self.assertAlmostEqual(self.img.symmetry('both')[1], 2.63, delta=0.02)
        self.assertRaises(ValueError, self.img.symmetry, method='varen')

    def test_determine_center(self):
        """Test the determined center is correct"""
        y, x = self.img._determine_center()
        self.assertAlmostEqual(y, 366, delta=4)
        self.assertAlmostEqual(x, 537, delta=4)

    def test_plot_flatness(self):
        """Test that plots work as expected."""
        self.img.plot_flatness()
        self.img.plot_flatness(plane='x')
        self.img.plot_flatness(plane='in')

    def test_plot_symmeetry(self):
        """Test that plots work as expected."""
        self.img.plot_symmetry()
        self.img.plot_symmetry(plane='x')
        self.img.plot_symmetry(plane='in')

    def test_plot_flatsym(self):
        """Test that plots work as expected."""
        self.img.plot_flatsym()
        self.img.plot_flatsym(plane='x')
        self.img.plot_flatsym(plane='in')



