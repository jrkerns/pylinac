"""Test the entire package; an underscore precedes this file name
so it does not include itself in the test discovery."""
import os.path as osp
from . import run_tests

test_dir = osp.dirname(__file__)
run_tests(test_dir)
