from io import BytesIO, StringIO
from tempfile import TemporaryFile


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