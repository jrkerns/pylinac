from __future__ import division, absolute_import, print_function

from unittest import TestCase

from pylinac.log_analyzer import MachineLog, log_types


class Test_dynalog(TestCase):
    """Tests of dynalog files, mostly using the demo file."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_Bfile_exists(self):
        """Test that a dynalog finds the associated B-file, and responds appropriately if it doesn't find it."""
        pass

    def test_log_type(self):
        """Ensure correct log type identification."""
        self.mlc.load_demo_dynalog()
        self.assertTrue(self.mlc.log_type == log_types[0])  # dynalog

    # def test_bad_ext(self):
    #     """Test that loading the badly-named dynalog still loads and identifies properly."""
    #     bad_log_path = osp.join(go_up_dirlevel(1), 'demo files', 'Bbadext')
    #     self.mlc.load_logfile(bad_log_path)
    #     self.assertTrue(self.mlc.log_type == log_types[0])

class Test_trajectory(TestCase):
    """Tests having to do with trajectory log files."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_log_type(self):
        """Ensure demo trajectory log is identified properly."""
        self.mlc.load_demo_trajectorylog()
        self.assertTrue(self.mlc.log_type == log_types[1])  # traj log type

    # def test_bad_ext(self):
    #     """Test that loading the badly-named trajectory log still loads and identifies properly"""
    #     bad_log_path = osp.join(go_up_dirlevel(1), 'demo files', 'tlogbadext')
    #     self.mlc.load_logfile(bad_log_path)
    #     self.assertTrue(self.mlc.log_type == log_types[1])

class Test_end2end(TestCase):
    """Assess an end-to-end test of both demo types."""
    def setUp(self):
        self.mlc = MachineLog()

    def test_dynalog_95th(self):
        """Test dynalog 95th percentile error."""
        pass
