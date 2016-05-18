"""A script to generate thumbnails of pylinac images for machine learning"""
import os.path as osp
import os
import concurrent.futures
import time

import matplotlib.pyplot as plt
import numpy as np
from pylinac import image
from scipy.misc import imresize
from sklearn import preprocessing


def get_image_files(folder):
    """Get a list of files that are valid images from the folder."""
    futures = {}
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as exec:
        for pdir, _, files in os.walk(folder):
            for file in files:
                filepath = osp.join(pdir, file)
                future = exec.submit(image.is_image, filepath)
                futures[future] = filepath
    filepaths = []
    for idx, future in enumerate(concurrent.futures.as_completed(futures)):
        if future.result():
            filepaths.append(futures[future])
    print("Done with {} in {:.2f}s".format(osp.basename(folder), time.time() - start))
    return filepaths


def process_image(path):
    """Load and resize the images and return as flattened numpy array"""
    img = image.load(path, dtype=np.float32)
    return imresize(img.array, size=(100, 100), mode='F').flatten()


def build_images():
    """Completely load, resize, and save the images for training. Main function."""
    # get image file paths for each image type
    path_stub = r'D:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files'
    pf_files = get_image_files(osp.join(path_stub, 'Picket Fences'))
    pipspro_files = get_image_files(osp.join(path_stub, '2D Image quality phantoms', 'PipsPro'))
    leeds_files = get_image_files(osp.join(path_stub, '2D Image quality phantoms', 'Leeds'))
    wl_files = get_image_files(osp.join(path_stub, 'Winston-Lutz'))
    # cbct_files = get_image_files(osp.join(path_stub, 'CBCTs'))
    filepaths = pf_files + pipspro_files + leeds_files + wl_files
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
        (np.repeat(0, len(pf_files)), np.repeat(1, len(pipspro_files)), np.repeat(2, len(leeds_files)), np.repeat(3, len(wl_files)))))
    print("Images build")


def load_images():
    """Load the built images for training."""
    return np.load('images.npy'), np.load('labels.npy')

