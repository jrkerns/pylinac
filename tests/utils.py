import concurrent.futures
from io import BytesIO, StringIO
import os
import os.path as osp
import time
from tempfile import TemporaryFile

from pylinac.core.image import Image

DATA_BANK_DIR = osp.abspath(osp.join('..', '..', 'pylinac test files'))


def save_file(method, *args, as_file_object=None, **kwargs):
    """Save a file using the passed method and assert it exists after saving.
    Also deletes the file after checking for existence."""
    if as_file_object is None:  # regular file
        with TemporaryFile(mode='w') as tfile:
            method(tfile, *args, **kwargs)
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


class LoadingTestBase:
    """This class can be used as a base for a module's loading test class."""
    klass = object
    constructor_input = None
    demo_method = 'from_demo_image'
    url = None
    zip = None

    @property
    def real_url(self):
        return 'https://s3.amazonaws.com/assuranceqa-staging/uploads/imgs/' + self.url

    def test_consructor(self):
        if self.constructor_input is not None:
            self.klass(self.constructor_input)

    def test_from_demo(self):
        if self.demo_method is not None:
            getattr(self.klass, self.demo_method)()

    def test_from_url(self):
        if self.url is not None:
            self.klass.from_url(self.real_url)

    def test_from_zip(self):
        if self.zip is not None:
            self.klass.from_zip(self.zip)


class DataBankMixin:
    """Mixin class for running through a bank of images/datasets. No details are tested; only pass/fail results.

    DATA_DIR must be set and is a relative location.
    ``file_should_be_processed()`` should be overloaded if the criteria is not an Image.
    ``test_all()`` should be overloaded to pass an analysis function. Because of the pool executor, methods nor attrs can be passed.
    """
    DATA_DIR = ''

    @classmethod
    def setUpClass(cls):
        cls.DATA_DIR = osp.join(DATA_BANK_DIR, cls.DATA_DIR)

    def file_should_be_processed(self, filepath):
        """Decision of whether file should be run. Returns boolean."""
        try:
            Image.load(filepath)
            return True
        except:
            return False

    def test_all(self, func):
        """Test method that runs through the bank data in a process pool.

        Parameters
        ----------
        func : A function that processes the filepath and determines if it passes. Must return ``Success`` if passed.
        """
        futures = []
        passes = 0
        fails = 0
        start = time.time()
        # create a process pool and walk through datasets
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, sdir, files in os.walk(self.DATA_DIR):
                for file in files:
                    filepath = osp.join(pdir, file)
                    # if data needs processing, submit it into the queue
                    if self.file_should_be_processed(filepath):
                        future = exec.submit(func, filepath)
                        futures.append(future)
            # return results
            for idx, future in enumerate(concurrent.futures.as_completed(futures)):
                if future.result() == 'Success':
                    passes += 1
                else:
                    fails += 1
                print(future.result(), idx)
        end = time.time() - start
        print('Processing of {} files took {:3.1f}s. {} passed; {} failed.'.format(len(futures), end, passes, fails))
