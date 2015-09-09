import os.path as osp
import time
import os


def save_file(filename, method):

    method(filename)
    time.sleep(0.15)  # sleep just to let OS work
    assert osp.isfile(filename), "Save file did not successfully save the image"
    # cleanup
    os.remove(filename)
    assert not osp.isfile(filename), "Save file test did not clean up saved image"


def test_point_equality(point1, point2):
    if point1.x != point2.x:
        raise ValueError("{} does not equal {}".format(point1.x, point2.x))
    if point1.y != point2.y:
        raise ValueError("{} does not equal {}".format(point1.y, point2.y))