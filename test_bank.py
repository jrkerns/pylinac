"""Test the entire package; an underscore precedes this file name
so it does not include itself in the test discovery."""
import os.path as osp
from tests_basic import run_tests


test_dir = osp.join(osp.dirname(__file__), 'tests_bank')
run_tests(test_dir, pattern='_test*.py')
