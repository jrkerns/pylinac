import unittest

from pylinac.core.utilities import *

class Test_Utilities(unittest.TestCase):

    def test_isnumeric(self):
        # test numerics
        numerics = (5, 3.2, np.array((5,6))[0])
        for number in numerics:
            self.assertTrue(isnumeric(number))

        # test non-numerics
        notnumerics = ('5', np.array((5,6)))
        for notnumeric in notnumerics:
            self.assertFalse(isnumeric(notnumeric))

    def test_go_up_dirlevel(self):
        path2here = osp.abspath(__file__)
        path_dir = osp.dirname(path2here)
        path_dirdir = osp.dirname(path_dir)

        dir = go_up_dirlevel()
        self.assertEqual(dir, path_dir)
        dir = go_up_dirlevel(1)
        self.assertEqual(dir, path_dirdir)

    def test_is_iterable(self):
        # test iterables
        iters = ((1,2,'t'), [4, 8, 'r'], np.array((5,6,7)))
        for iter in iters:
            self.assertTrue(is_iterable(iter))
        # test non-iterables
        noniters = (5,)
        for iter in noniters:
            self.assertFalse(is_iterable(iter))