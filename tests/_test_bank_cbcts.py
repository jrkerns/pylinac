"""Run through the CBCT image bank. Since the criteria and analysis are different than the default DataBankMixin, we must create our own."""
import concurrent.futures
import time
import os
import os.path as osp
from unittest import TestCase

from pylinac import CBCT
from tests import TEST_BANK_DIR


def run_cbct(path):
    """Function to pass to the process pool executor to process cbct images."""
    try:
        mypf = CBCT(path)
        mypf.analyze()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


class TestCBCTBank(TestCase):
    image_bank_dir = osp.join(TEST_BANK_DIR, 'CBCTs')

    def test_all(self):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.image_bank_dir):
                if files and files[0].endswith('.dcm'):
                    future = exec.submit(run_cbct, pdir)
                    futures.append(future)
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                print(future.result(), idx)
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))
