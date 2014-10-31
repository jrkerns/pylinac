from __future__ import division, absolute_import, print_function

from unittest import TestCase
import os.path as osp

from pylinac.loganalyzer.log_analyzer import MachineLog, log_types
from pylinac.common.common_functions import go_up_dirlevel


class Test_loading_schemes(TestCase):
    """A test class for ensuring that MachineLog will open log files and handle appropriately if it can't."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_load_demo_dyn(self):
        """Load the dynalog demo files properly."""
        self.mlc.load_demo_file_dynalog()
        self.assertTrue(self.mlc.log_type == log_types[0])  # ensure dynalog type

    def test_load_demo_traj(self):
        """Ensure demo trajectory log loads properly."""
        self.mlc.load_demo_file_trajectorylog()
        self.assertTrue(self.mlc.log_type == log_types[1])  # traj log type

    def test_bad_dyn_ext(self):
        """Test that loading the badly-named dynalog still loads properly"""
        bad_log_path = osp.join(go_up_dirlevel(1), 'demo files', 'Bbadext')
        self.mlc.load_logfile(bad_log_path)
        self.assertTrue(self.mlc.log_type == log_types[0])

    def test_bad_traj_ext(self):
        """Test that loading the badly-named trajectory log still loads properly"""
        bad_log_path = osp.join(go_up_dirlevel(1), 'demo files', 'tlogbadext')
        self.mlc.load_logfile(bad_log_path)
        self.assertTrue(self.mlc.log_type == log_types[1])


class Test_Dynalog_analysis(TestCase):
    """Tests having to do with dynalog files."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_Bfile_exists(self):
        """Test that a dynalog finds the associated B-file, and responds appropriately if it doesn't find it."""
        pass

class Test_Trajectory_analysis(TestCase):
    """Tests having to do with trajectory log files."""
    pass

class Test_end2end(TestCase):
    """Assess an end-to-end test of both demo types."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_dynalog_95th(self):
        """Test dynalog 95th percentile error."""
        pass
