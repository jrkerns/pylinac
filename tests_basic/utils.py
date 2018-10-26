import concurrent.futures
from io import BytesIO, StringIO
import multiprocessing
import os
import os.path as osp
import pprint
import time
from tempfile import TemporaryDirectory
from urllib.request import urlopen

from pylinac.core import image


def has_www_connection():
    try:
        with urlopen('http://www.google.com') as r:
            return r.status == 200
    except:
        return False


def save_file(method, *args, as_file_object=None, **kwargs):
    """Save a file using the passed method and assert it exists after saving.
    Also deletes the file after checking for existence."""
    if as_file_object is None:  # regular file
        with TemporaryDirectory() as tmpdir:
            tmpfile = osp.join(tmpdir, 'myfile')
            method(tmpfile, *args, **kwargs)
    else:
        if 'b' in as_file_object:
            temp = BytesIO
        elif 's' in as_file_object:
            temp = StringIO
        with temp() as t:
            method(t, *args, **kwargs)


def test_point_equality(point1, point2):
    if point1.x != point2.x:
        raise ValueError("{} does not equal {}".format(point1.x, point2.x))
    if point1.y != point2.y:
        raise ValueError("{} does not equal {}".format(point1.y, point2.y))


class LocationMixin:
    """A mixin that provides attrs and a method to get the absolute location of the
    file using syntax.

    Usage
    -----
    0. Override ``dir_location`` with the directory of the destination.
    1. Override ``file_path`` with a list that contains the subfolder(s) and file name.
    """
    dir_location = ''
    file_path = []

    @classmethod
    def get_filename(cls):
        """Return the canonical path to the file."""
        return osp.join(cls.dir_location, *cls.file_path)


class LoadingTestBase:
    """This class can be used as a base for a module's loading test class."""
    klass = object
    constructor_input = None
    demo_load_method = 'from_demo_image'
    url = None
    zip = None
    kwargs = {}

    @property
    def full_url(self):
        return 'https://s3.amazonaws.com/pylinac/' + self.url

    def test_consructor(self):
        if self.constructor_input is not None:
            self.klass(self.constructor_input, **self.kwargs)

    def test_from_demo(self):
        if self.demo_load_method is not None:
            getattr(self.klass, self.demo_load_method)()

    def test_from_url(self):
        if self.url is not None:
            self.klass.from_url(self.full_url, **self.kwargs)

    def test_from_zip(self):
        if self.zip is not None:
            self.klass.from_zip(self.zip, **self.kwargs)


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
    DATA_BANK_DIR = osp.abspath(osp.join(osp.abspath(__file__), '..', '..', '..', 'pylinac test files'))
    DATA_DIR = []
    executor = 'ProcessPoolExecutor'
    workers = multiprocessing.cpu_count() - 1
    print_success_path = False
    write_failures_to_file = False

    @classmethod
    def setUpClass(cls):
        cls.DATA_DIR = osp.join(cls.DATA_BANK_DIR, *cls.DATA_DIR)
        if not osp.isdir(cls.DATA_DIR):
            raise NotADirectoryError("Directory {} is not valid".format(cls.DATA_DIR))

    def file_should_be_processed(self, filepath):
        """Decision of whether file should be run. Returns boolean."""
        try:
            image.load(filepath)
            return True
        except:
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
        with getattr(concurrent.futures, self.executor)(max_workers=self.workers) as exec:
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
                if future.result() == 'Success':
                    passes += 1
                    if self.print_success_path:
                        stuff_to_print.append(futures[future])
                else:
                    fails += [futures[future]]
                print(*stuff_to_print)

        end = time.time() - start
        print('Processing of {} files took {:3.1f}s ({:3.2f}s/item). {} passed; {} failed.'.format(test_num, end, end/test_num, passes, len(fails)))
        if len(fails) > 0:
            pprint.pprint("Failures: {}".format(fails))
            if self.write_failures_to_file:
                with open('failures_{}.txt'.format(osp.basename(self.DATA_DIR)), mode='w') as f:
                    for file in fails:
                        f.write(file + '\n')
                print("Failures written to file")
