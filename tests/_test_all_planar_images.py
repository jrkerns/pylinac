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
    """Function to pass to the process pool executor to process leeds phantom images."""
    try:
        # print('Processing: {}'.format(path))
        leeds = LeedsTOR(path)
        leeds.analyze()
        leeds.plot_analyzed_image()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


def run_pipspro(path):
    """Function to pass to the process pool executor to process pipspro QC3 images."""
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
    IMAGE_DIR = ''

    @classmethod
    def setUpClass(cls):
        cls.IMAGE_DIR = osp.join(IMAGE_BANK_DIR, cls.IMAGE_DIR)

    def test_all(self, func):
        futures = []
        passes = 0
        fails = 0
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
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                if future.result() == 'Success':
                    passes += 1
                else:
                    fails += 1
                print(future.result(), idx)
        end = time.time() - start
        print('Processing of {} files took {:3.1f}s. {} passed; {} failed.'.format(len(futures), end, passes, fails))


class TestLeedsImageBank(ImageBankMixin, TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""
    IMAGE_DIR = 'Leeds'

    def test_all(self):
        super().test_all(run_leeds)


class TestPipsProImageBank(ImageBankMixin, TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""
    IMAGE_DIR = 'PipsPro'

    def test_all(self):
        super().test_all(run_pipspro)
