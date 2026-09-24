"""Test suite for the pylinac.io module."""

import os
import os.path as osp
import unittest

from pylinac import Interpolation
from pylinac.core.io import (
    LoadableMixin,
    Path,
    SNCProfiler,
    TemporaryZipDirectory,
    URLError,
    get_url,
    is_dicom,
)
from pylinac.core.profile import SingleProfile
from tests_basic.utils import get_file_from_cloud_test_repo


class TestIO(unittest.TestCase):
    def test_temp_zip_dir(self):
        """Test the TemporaryZipDirectory."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])

        # test context manager use; shouldn't raise
        with TemporaryZipDirectory(zfile) as tmpzip:
            files = [osp.join(tmpzip, file) for file in os.listdir(tmpzip)]
            # test that they are real files
            self.assertTrue(osp.isfile(files[0]))
            # test that both images were unpacked
            self.assertEqual(len(files), 2, msg="There were not 2 files found")

    def test_get_url(self):
        """Test the URL retreiver."""
        # test webpage
        webpage_url = "http://google.com"
        get_url(webpage_url)  # shouldn't raise
        # test file
        file_url = "https://storage.googleapis.com/pylinac_demo_files/winston_lutz.zip"
        local_file = get_url(file_url)
        osp.isfile(local_file)
        # bad URL
        with self.assertRaises(URLError):
            get_url("http://asdfasdfasdfasdfasdfasdfasdfasdf.org")

    def test_is_dicom(self):
        """Test the is_dicom function."""

        test_file = get_file_from_cloud_test_repo(["VMAT", "DRGSdmlc-105-example.dcm"])
        invalid_file = test_file.replace("DR", "DR_")
        notdicom_file = osp.abspath(__file__)

        # valid file returns True
        self.assertTrue(is_dicom(test_file))

        # return false for real file but not dicom
        self.assertFalse(is_dicom(notdicom_file))

        # test invalid path
        self.assertRaises(IOError, is_dicom, invalid_file)


class TestTempZipDir(unittest.TestCase):
    def test_dir_is_deleted_normally(self):
        """Test that the directory is deleted normally."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])
        with TemporaryZipDirectory(zfile) as tmpzip:
            self.assertTrue(osp.isdir(tmpzip))
        self.assertFalse(osp.exists(tmpzip))

    def test_dir_remains_when_delete_false(self):
        """Test that the directory is not deleted when delete=False."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])
        with TemporaryZipDirectory(zfile, delete=False) as tmpzip:
            self.assertTrue(osp.isdir(tmpzip))
        self.assertTrue(osp.exists(tmpzip))


class TestSNCProfiler(unittest.TestCase):
    def test_loading(self):
        path = get_file_from_cloud_test_repo(["9E-GA0.prs"])
        prof = SNCProfiler(path)
        self.assertEqual(len(prof.detectors), 254)

    def test_to_profiles(self):
        path = get_file_from_cloud_test_repo(["9E-GA0.prs"])
        prof = SNCProfiler(path)
        profs = prof.to_profiles()
        self.assertEqual(len(profs), 4)
        self.assertIsInstance(profs[0], SingleProfile)

    def test_detectors(self):
        path = get_file_from_cloud_test_repo(["6XFFF.prs"])
        prof = SNCProfiler(path)
        profs = prof.to_profiles(interpolation=Interpolation.NONE)
        self.assertEqual(len(profs), 4)
        crossplane_prof = profs[0]
        self.assertTrue(crossplane_prof[30] < crossplane_prof[31])
        self.assertTrue(crossplane_prof[32] < crossplane_prof[31])
        self.assertEqual(len(crossplane_prof), 63)
        self.assertEqual(len(profs[1]), 65)


class TestLoadableMixin(unittest.TestCase):
    """Unit tests for the LoadableMixin behaviors."""

    def _make_zip(self, files: dict[str, bytes]) -> str:
        """Create a temporary zip file containing given files. Returns path."""
        import tempfile
        import zipfile

        fd, zpath = tempfile.mkstemp(suffix=".zip")
        os.close(fd)
        with zipfile.ZipFile(zpath, "w") as zf:
            for name, data in files.items():
                zf.writestr(name, data)
        return zpath

    def test_from_zip_uses_tmpdir_when_no_glob(self):
        class Dummy(LoadableMixin):
            def __init__(self, path, **kwargs):
                # record the received argument
                self.init_arg = path

        zpath = self._make_zip({"a.txt": b"1"})
        inst = Dummy.from_zip(zpath)
        # Should have received a path-like string to the extracted folder
        self.assertTrue(isinstance(inst.init_arg, str))
        self.assertTrue(len(inst.init_arg) > 0)

    def test_from_zip_collects_files_when_glob_set(self):
        class DummyList(LoadableMixin):
            ZIP_FILE_GLOB = "*.txt"

            def __init__(self, files, **kwargs):
                # files should be an iterable of Path objects
                self.files = list(files)

        zpath = self._make_zip({"b.txt": b"b", "a.txt": b"a"})
        inst = DummyList.from_zip(zpath)
        # Should be list of Path objects sorted alphabetically
        self.assertEqual(len(inst.files), 2)
        self.assertEqual([p.name for p in inst.files], ["a.txt", "b.txt"])

    def test_from_url_dispatches_to_from_zip_or_constructor(self):
        # Prepare files
        import tempfile

        fd, nonzip = tempfile.mkstemp(suffix=".bin")
        os.close(fd)
        with open(nonzip, "wb") as f:
            f.write(b"data")
        zpath = self._make_zip({"x.txt": b"x"})

        # Monkeypatch get_url in the module
        import pylinac.core.io as io_mod

        original_get_url = io_mod.get_url

        try:
            io_mod.get_url = lambda url, progress_bar=True: zpath

            class FromZipDummy(LoadableMixin):
                ZIP_FILE_GLOB = "*.txt"

                def __init__(self, files, **kwargs):
                    self.files = list(files)

            inst = FromZipDummy.from_url("http://example.com/fake.zip")
            self.assertEqual(len(inst.files), 1)

            # Now non-zip
            io_mod.get_url = lambda url, progress_bar=True: nonzip

            class FromFileDummy(LoadableMixin):
                def __init__(self, path, **kwargs):
                    self.path = path

            inst2 = FromFileDummy.from_url("http://example.com/file.bin")
            self.assertTrue(inst2.path.endswith(".bin"))
        finally:
            io_mod.get_url = original_get_url
            try:
                os.remove(nonzip)
            except Exception:
                pass

    def test_from_demo_uses_retrieve_demo_file_and_errors(self):
        """Cover demo error cases and delegate behavior for demo/image helpers."""

        import pylinac.core.io as io_mod

        # Case: no DEMO_FILES configured
        class NoDemo(LoadableMixin):
            DEMO_FILES = None

        with self.assertRaises(AttributeError):
            NoDemo.from_demo()

        # Case: empty DEMO_FILES list
        class EmptyDemo(LoadableMixin):
            DEMO_FILES = []

        with self.assertRaises(AttributeError):
            EmptyDemo.from_demo()

        # Case: empty but truthy DEMO_FILES (force the len==0 branch)
        class TruthyEmpty(list):
            def __bool__(self):
                return True

        class WeirdDemo(LoadableMixin):
            DEMO_FILES = TruthyEmpty()

        with self.assertRaises(AttributeError):
            # this should hit the len(cls.DEMO_FILES) == 0 branch
            WeirdDemo.from_demo()

        # Case: delegate to from_zip or constructor via retrieve_demo_file
        import tempfile

        fd, nonzip = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        with open(nonzip, "wb") as f:
            f.write(b"x")
        zpath = self._make_zip({"y.txt": b"y"})

        original_retrieve = io_mod.retrieve_demo_file
        try:
            # demo returns zip
            io_mod.retrieve_demo_file = lambda name: Path(zpath)

            class DemoZip(LoadableMixin):
                DEMO_FILES = ["demo.zip"]
                ZIP_FILE_GLOB = "*.txt"

                def __init__(self, files, **kwargs):
                    self.files = list(files)

            inst = DemoZip.from_demo()
            self.assertEqual(len(inst.files), 1)

            # demo returns non-zip
            io_mod.retrieve_demo_file = lambda name: Path(nonzip)

            class DemoFile(LoadableMixin):
                DEMO_FILES = ["demo.dat"]

                def __init__(self, path, **kwargs):
                    self.path = path

            inst2 = DemoFile.from_demo()
            self.assertTrue(str(inst2.path).endswith(".dat"))

            # Test from_demo_image/from_demo_images delegate to from_demo
            called = {}

            io_mod.retrieve_demo_file = lambda name: Path(nonzip)

            class DemoImage(LoadableMixin):
                DEMO_FILES = ["demo.dat"]

                def __init__(self, path, **kwargs):
                    called["init"] = path

            # call the convenience wrappers
            DemoImage.from_demo_image()
            self.assertIn("init", called)
            called.clear()
            DemoImage.from_demo_images()
            self.assertIn("init", called)

            # Test from_image directly
            class FromImageDummy(LoadableMixin):
                def __init__(self, image, **kwargs):
                    self.image = image

            inst_img = FromImageDummy.from_image("rawimage")
            self.assertEqual(inst_img.image, "rawimage")

        finally:
            io_mod.retrieve_demo_file = original_retrieve
            try:
                os.remove(nonzip)
            except Exception:
                pass
