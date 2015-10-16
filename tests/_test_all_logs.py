import concurrent.futures
import os
import os.path as osp
import time
from unittest import TestCase

from pylinac import MachineLog


def run_log(path):
    """Function to pass to the process pool executor to process machine logs."""
    try:
        mypf = MachineLog(path)
        mypf.fluence.actual.calc_map()
        return 'Success'
    except:
        return 'Failure at {}'.format(path)


class TestImageBank(TestCase):
    """Test the Machine Logs in the large image bank. No details are tested other than that it runs."""
    log_bank_dir = osp.abspath(osp.join('..', '..', 'unorganized linac data', 'Machine logs'))

    def test_all(self):
        futures = []
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.log_bank_dir):
                for file in files:
                    filepath = osp.join(pdir, file)
                    if filepath.endswith('.dlg') or filepath.endswith('.bin'):
                        future = exec.submit(run_log, filepath)
                        futures.append(future)
            for future in concurrent.futures.as_completed(futures):
                print(future.result())
        end = time.time() - start
        print('Processing of {} files took {}s'.format(len(futures), end))
