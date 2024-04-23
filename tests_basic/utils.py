import base64
import concurrent.futures
import contextlib
import hashlib
import multiprocessing
import os
import os.path as osp
import pprint
import shutil
import time
from functools import lru_cache
from io import BytesIO, StringIO
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from typing import Callable, List, Sequence, Union
from urllib.request import urlopen

from google.cloud import storage
from py_linq import Enumerable

from pylinac.core import image
from tests_basic import DELETE_FILES

GCP_BUCKET_NAME = "pylinac_test_files"
LOCAL_TEST_DIR = "test_files"

# make the local test dir if it doesn't exist
if not osp.isdir(osp.join(osp.dirname(__file__), LOCAL_TEST_DIR)):
    os.mkdir(osp.join(osp.dirname(__file__), LOCAL_TEST_DIR))


@contextlib.contextmanager
def access_gcp() -> storage.Client:
    # access GCP
    credentials_file = Path(__file__).parent.parent / "GCP_creds.json"
    # check if the credentials file is available (local dev)
    # if not, load from the env var (test pipeline)
    if not credentials_file.is_file():
        with open(credentials_file, "wb") as f:
            creds = base64.b64decode(os.environ.get("GOOGLE_CREDENTIALS", ""))
            f.write(creds)
    client = storage.Client.from_service_account_json(str(credentials_file))
    try:
        yield client
    finally:
        del client


@lru_cache
def gcp_bucket_object_list(bucket_name: str) -> list:
    with access_gcp() as storage_client:
        return list(storage_client.list_blobs(bucket_name))


def get_folder_from_cloud_test_repo(folder: List[str], skip_exists: bool = True) -> str:
    """Get a folder from GCP.

    Parameters
    ----------
    skip_exists
        If True, only checks that the destination folder exists and isn't empty.
        This is helpful for avoiding network calls since querying GCP can cost significant time.
    """
    dest_folder = Path(LOCAL_TEST_DIR, *folder)
    if skip_exists and dest_folder.exists() and len(list(dest_folder.iterdir())) > 0:
        return str(dest_folder)
    # get the folder data
    all_blobs = gcp_bucket_object_list(GCP_BUCKET_NAME)
    blobs = (
        Enumerable(all_blobs)
        .where(lambda b: len(b.name.split("/")) > len(folder))
        .where(lambda b: b.name.split("/")[1] != "")
        .where(
            lambda b: all(f in b.name.split("/")[idx] for idx, f in enumerate(folder))
        )
        .to_list()
    )

    # make root folder if need be
    dest_folder = osp.join(osp.dirname(__file__), LOCAL_TEST_DIR, *folder)
    if not osp.isdir(dest_folder):
        os.makedirs(dest_folder)

    # make subfolders if need be
    subdirs = [b.name.split("/")[1:-2] for b in blobs if len(b.name.split("/")) > 2]
    dest_sub_folders = [
        osp.join(osp.dirname(__file__), LOCAL_TEST_DIR, *folder, *f) for f in subdirs
    ]
    for sdir in dest_sub_folders:
        if not osp.isdir(sdir):
            os.makedirs(sdir)

    for blob in blobs:
        # download file
        path = osp.join(dest_folder, blob.name.split("/")[-1])
        if not os.path.exists(path):
            blob.download_to_filename(path)

    return osp.join(osp.dirname(__file__), LOCAL_TEST_DIR, *folder)


def get_file_from_cloud_test_repo(path: List[str], force: bool = False) -> str:
    """Get a single file from GCP storage. Returns the path to disk it was downloaded to"""
    local_filename = osp.join(osp.dirname(__file__), LOCAL_TEST_DIR, *path)
    if osp.isfile(local_filename) and not force:
        with open(local_filename, "rb") as f:
            file_hash = hashlib.md5(f.read()).hexdigest()
            print(f"Local file found: {local_filename}@{file_hash}")
        return local_filename
    with access_gcp() as client:
        bucket = client.bucket(GCP_BUCKET_NAME)
        blob = bucket.blob(
            str(PurePosixPath(*path))
        )  # posix because google storage is on unix and won't find path w/ windows path
        # make any necessary subdirs leading up to the file
        if len(path) > 1:
            for idx in range(1, len(path)):
                local_dir = osp.join(osp.dirname(__file__), LOCAL_TEST_DIR, *path[:idx])
                if not osp.isdir(local_dir):
                    os.mkdir(local_dir)

        blob.download_to_filename(local_filename)
        time.sleep(2)
        print(
            f"Downloaded from GCP: {local_filename}@{hashlib.md5(open(local_filename, 'rb').read()).hexdigest()}"
        )
        return local_filename


def has_www_connection():
    try:
        with urlopen("http://www.google.com") as r:
            return r.status == 200
    except Exception:
        return False


def save_file(method, *args, as_file_object=None, to_single_file=True, **kwargs):
    """Save a file using the passed method and assert it exists after saving.
    Also deletes the file after checking for existence."""
    if as_file_object is None:  # regular file
        with TemporaryDirectory() as tmpdir:
            if to_single_file:
                tmpfile = osp.join(tmpdir, "myfile")
                method(tmpfile, *args, **kwargs)
            else:
                method(tmpdir, *args, **kwargs)
    else:
        if "b" in as_file_object:
            temp = BytesIO
        elif "s" in as_file_object:
            temp = StringIO
        with temp() as t:
            method(t, *args, **kwargs)


def point_equality_validation(point1, point2):
    if point1.x != point2.x:
        raise ValueError(f"{point1.x} does not equal {point2.x}")
    if point1.y != point2.y:
        raise ValueError(f"{point1.y} does not equal {point2.y}")


class CloudFileMixin:
    """A mixin that provides attrs and a method to get the absolute location of the
    file using syntax.

    Usage
    -----
    0. Override ``dir_location`` with the directory of the destination.
    1. Override ``file_path`` with a list that contains the subfolder(s) and file name.
    """

    file_name: Union[str, Sequence[str]]
    dir_path: Sequence[str]
    delete_file = True

    @classmethod
    def get_filename(cls) -> str:
        """Return the canonical path to the file on disk. Download if it doesn't exist."""
        if isinstance(cls.file_name, (list, tuple)):
            full_path = Path(LOCAL_TEST_DIR, *cls.dir_path, *cls.file_name).absolute()
        else:
            full_path = Path(LOCAL_TEST_DIR, *cls.dir_path, cls.file_name).absolute()
        if full_path.is_file():
            return str(full_path)
        else:
            if isinstance(cls.file_name, (list, tuple)):
                return get_file_from_cloud_test_repo([*cls.dir_path, *cls.file_name])
            else:
                return get_file_from_cloud_test_repo([*cls.dir_path, cls.file_name])

    @classmethod
    def tearDownClass(cls):
        if cls.delete_file and DELETE_FILES:
            file = cls.get_filename()
            if osp.isfile(file):
                os.remove(file)
            elif osp.isdir(file):
                shutil.rmtree(file)


class MixinTesterBase:
    klass = Callable


class FromZipTesterMixin(MixinTesterBase):
    zip: List[str]
    zip_kwargs: dict = {}

    def test_from_zip(self):
        file = get_file_from_cloud_test_repo(self.zip)
        inst = self.klass.from_zip(file, **self.zip_kwargs)
        self.assertIsInstance(inst, self.klass)


class FromDemoImageTesterMixin(MixinTesterBase):
    demo_load_method: str = "from_demo_image"

    def test_from_demo(self):
        inst = getattr(self.klass, self.demo_load_method)()
        self.assertIsInstance(inst, self.klass)


class InitTesterMixin(MixinTesterBase):
    init_file: List[str]
    init_kwargs: dict = {}
    is_folder = False

    @property
    def full_init_file(self) -> str:
        if self.is_folder:
            return get_folder_from_cloud_test_repo(self.init_file)
        else:
            return get_file_from_cloud_test_repo(self.init_file)

    def test_init(self):
        inst = self.klass(self.full_init_file, **self.init_kwargs)
        self.assertIsInstance(inst, self.klass)


class FromURLTesterMixin(MixinTesterBase):
    url: str
    url_kwargs: dict = {}

    @property
    def full_url(self):
        return r"https://storage.googleapis.com/pylinac_demo_files/" + self.url

    def test_from_url(self):
        self.klass.from_url(self.full_url, **self.url_kwargs)


class DataBankMixin:
    """Mixin class for running through a bank of images/datasets. No details are tested; only pass/fail results.

    DATA_DIR must be set and is a relative location.
    ``file_should_be_processed()`` should be overloaded if the criteria is not an Image.
    ``test_all()`` should be overloaded to pass an analysis function. Because of the pool executor, methods nor attrs can be passed.
    ``executor`` can be either 'ProcessPoolExecutor' or 'ThreadPoolExecutor'.
    ``workers`` specifies how many workers (threads or processes) will be started.

    This class walks through a directory and runs through the files, analyzing them if the file should be analyzed.
    Note that an executor is started **per directory** rather than one executor for all the files.
    This is on purpose.
    Some test runs, seemingly due to the executor, bog down when running a very large number of files. By opening a new executor at
    every directory, memory leaks are minimized.
    """

    DATA_BANK_DIR = osp.abspath(
        osp.join(osp.abspath(__file__), "..", "..", "..", "pylinac test files")
    )
    DATA_DIR = []
    executor = "ProcessPoolExecutor"
    workers = multiprocessing.cpu_count() - 1
    print_success_path = False
    write_failures_to_file = False

    @classmethod
    def setUpClass(cls):
        cls.DATA_DIR = osp.join(cls.DATA_BANK_DIR, *cls.DATA_DIR)
        if not osp.isdir(cls.DATA_DIR):
            raise NotADirectoryError(f"Directory {cls.DATA_DIR} is not valid")

    def file_should_be_processed(self, filepath):
        """Decision of whether file should be run. Returns boolean."""
        try:
            image.load(filepath)
            return True
        except Exception:
            return False

    def test_all(self, func):
        """Test method that runs through the bank data in an executor pool.

        Parameters
        ----------
        func : A function that processes the filepath and determines if it passes. Must return ``Success`` if passed.
        """
        passes = 0
        fails = []
        start = time.time()
        futures = {}
        # open an executor
        with getattr(concurrent.futures, self.executor)(
            max_workers=self.workers
        ) as exec:
            # walk through datasets
            for pdir, sdir, files in os.walk(self.DATA_DIR):
                for file in files:
                    # if the file needs processing, submit it into the queue
                    filepath = osp.join(pdir, file)
                    if self.file_should_be_processed(filepath):
                        future = exec.submit(func, filepath)
                        futures[future] = filepath

            # return results
            for test_num, future in enumerate(concurrent.futures.as_completed(futures)):
                stuff_to_print = [test_num, future.result()]
                if future.result() == "Success":
                    passes += 1
                    if self.print_success_path:
                        stuff_to_print.append(futures[future])
                else:
                    fails += [futures[future]]
                print(*stuff_to_print)

        end = time.time() - start
        print(
            "Processing of {} files took {:3.1f}s ({:3.2f}s/item). {} passed; {} failed.".format(
                test_num, end, end / test_num, passes, len(fails)
            )
        )
        if len(fails) > 0:
            pprint.pprint(f"Failures: {fails}")
            if self.write_failures_to_file:
                with open(f"failures_{osp.basename(self.DATA_DIR)}.txt", mode="w") as f:
                    for file in fails:
                        f.write(file + "\n")
                print("Failures written to file")
