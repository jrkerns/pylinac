import os
import os.path as osp

import matplotlib

TEST_FILES_DIR = osp.abspath(osp.join(osp.dirname(__file__), "test_files"))
TEST_BANK_DIR = osp.abspath(
    osp.join(osp.dirname(__file__), "..", "..", "pylinac test files")
)
HIDE_PLOTS = True

if os.environ.get("CI") or HIDE_PLOTS:
    matplotlib.use("Agg")

DELETE_FILES = bool(os.environ.get("DELETE_FILES", default=False))
