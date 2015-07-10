import os.path as osp
import time
import os


def save_file(filename, method):

    method(filename)
    time.sleep(0.1)  # sleep just to let OS work
    assert osp.isfile(filename), "Save file did not successfully save the image"
    # cleanup
    os.remove(filename)
    assert not osp.isfile(filename), "Save file test did not clean up saved image"
