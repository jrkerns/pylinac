import concurrent.futures
import os
import os.path as osp
import time
from unittest import TestCase

import matplotlib.pyplot as plt

from pylinac.core.image import Image
from pylinac import LeedsTOR, PipsProQC3

IMAGE_BANK_DIR = osp.abspath(osp.join('..', '..', 'unorganized linac data', '2D Image quality phantoms'))


def run_leeds(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        # print('Processing: {}'.format(path))
        leeds = LeedsTOR(path)
        leeds.analyze()
        leeds.plot_analyzed_image()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


def run_pipspro(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        pp = PipsProQC3(path)
        pp.analyze()
        pp.plot_analyzed_image()
    except:
        return 'Failure at {}'.format(path)
    else:
        plt.close('all')
        return 'Success'


class ImageBankMixin:
    func = object
    IMAGE_DIR = ''

    @classmethod
    def setUpClass(cls):
        cls.IMAGE_DIR = osp.join(IMAGE_BANK_DIR, cls.IMAGE_DIR)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_all(self, func):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.IMAGE_DIR):
                for file in files:
                    filepath = osp.join(pdir, file)
                    try:
                        Image.load(filepath)
                    except:
                        pass
                    else:
                        future = exec.submit(func, filepath)
                        futures.append(future)
            for future in concurrent.futures.as_completed(futures):
                if future.result() != 'Success':
                    print(future.result())
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))


class TestLeedsImageBank(TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""
    IMAGE_DIR = 'Leeds'

    @classmethod
    def setUpClass(cls):
        cls.IMAGE_DIR = osp.join(IMAGE_BANK_DIR, cls.IMAGE_DIR)

    def test_all(self):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.IMAGE_DIR):
                for file in files:
                    filepath = osp.join(pdir, file)
                    try:
                        Image.load(filepath)
                    except:
                        pass
                    else:
                        future = exec.submit(run_leeds, filepath)
                        futures.append(future)
            for future in concurrent.futures.as_completed(futures):
                if future.result() != 'Success':
                    print(future.result())
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))


class TestPipsProImageBank(ImageBankMixin, TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""
    IMAGE_DIR = 'PipsPro'

    @classmethod
    def setUpClass(cls):
        cls.IMAGE_DIR = osp.join(IMAGE_BANK_DIR, cls.IMAGE_DIR)

    def test_all(self):
        # super().test_all(run_pipspro)
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.IMAGE_DIR):
                for file in files:
                    filepath = osp.join(pdir, file)
                    try:
                        Image.load(filepath)
                    except:
                        pass
                    else:
                        future = exec.submit(run_pipspro, filepath)
                        futures.append(future)
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                # if future.result() != 'Success':
                print(future.result(), idx)
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))
