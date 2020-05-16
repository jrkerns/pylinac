from pathlib import Path
import warnings

from tests_basic import run_tests

warnings.filterwarnings("ignore")


run_tests(Path(__file__).parent)
