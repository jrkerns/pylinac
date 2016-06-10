"""A script to generate thumbnails of pylinac images for machine learning"""
import os.path as osp
import concurrent.futures
import time

import numpy as np
from pylinac import image
from sklearn import preprocessing

from machinelearning.tools import get_files, is_dicom, process_image


def build_images():
    """Completely load, resize, and save the images for training. Main function."""
    # get image file paths for each image type
    path_stub = r'D:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files'
    pf_files = get_files(osp.join(path_stub, 'Picket Fences'), is_dicom)
    pipspro_files = get_files(osp.join(path_stub, '2D Image quality phantoms', 'PipsPro'), is_dicom)
    leeds_files = get_files(osp.join(path_stub, '2D Image quality phantoms', 'Leeds'), is_dicom)
    # wl_files = get_files(osp.join(path_stub, 'Winston-Lutz'), is_dicom)
    filepaths = pf_files + pipspro_files + leeds_files # + wl_files
    print("{} files found".format(len(filepaths)))

    # preallocate
    total_array = np.zeros((len(filepaths), 10000), dtype=np.float32)
    print("Training array preallocated")

    # resize each image and add to a training array
    start = time.time()
    futures = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as exec:
        for idx, path in enumerate(filepaths):
            future = exec.submit(process_image, path)
            futures[future] = idx
    for idx, future in enumerate(concurrent.futures.as_completed(futures)):
        total_array[futures[future], :] = future.result()
    print("Training array set in {:.2f}s".format(time.time() - start))

    # feature scale the images
    scaled_array = preprocessing.minmax_scale(total_array, feature_range=(0, 1), axis=1)
    print("Training array scaled")

    # save arrays to disk for future use
    np.save(osp.join(osp.dirname(osp.abspath(__file__)), 'images'), scaled_array)
    np.save(osp.join(osp.dirname(osp.abspath(__file__)), 'labels'), np.concatenate(
        (np.repeat(0, len(pf_files)), np.repeat(1, len(pipspro_files)), np.repeat(2, len(leeds_files)))))
    print("Images build")


if __name__ == '__main__':
    build_images()

