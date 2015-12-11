import concurrent.futures
import os
import os.path as osp
import time
from unittest import TestCase

from pylinac.core.image import Image
from pylinac import LeedsTOR

IMAGE_BANK_DIR = osp.abspath(osp.join('..', '..', 'unorganized linac data', '2D Image quality phantoms', 'Leeds'))


def run_leeds(path):
    """Function to pass to the process pool executor to process picket fence images."""
    try:
        # print('Processing: {}'.format(path))
        leeds = LeedsTOR(path)
        leeds.analyze()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


class TestImageBank(TestCase):
    """Test the picket fences in the large image bank. Only tests the analysis runs; no details are tested."""

    def test_all(self):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(IMAGE_BANK_DIR):
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