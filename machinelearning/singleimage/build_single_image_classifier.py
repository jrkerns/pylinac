"""A script to generate thumbnails of pylinac images for machine learning"""
import os.path as osp
import concurrent.futures
import time

import numpy as np
from pylinac import image

from machinelearning.tools import get_files, is_dicom, process_image


def build_images(use_pool=True):
    """Completely load, resize, and save the images for training. Main function."""
    # get image file paths for each image type
    path_stub = r'D:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files'
    # fetch path filenames
    pf_files = get_files(osp.join(path_stub, 'Picket Fences'), is_dicom, use_pool=True)
    pipspro_files = get_files(osp.join(path_stub, '2D Image quality phantoms', 'QC-3'), is_dicom)
    leeds_files = get_files(osp.join(path_stub, '2D Image quality phantoms', 'Leeds'), is_dicom)
    star_files = get_files(osp.join(path_stub, 'Starshots'), image.is_image, use_pool=True)
    wl_files = get_files(osp.join(path_stub, 'Winston-Lutz'), is_dicom, use_pool=True)
    # vmat_files = get_files(osp.join(path_stub, 'VMATs'), is_dicom, use_pool=True)
    lv_files = get_files(osp.join(path_stub, '2D Image quality phantoms', 'Las Vegas'), is_dicom)
    filepaths = pf_files + pipspro_files + leeds_files + star_files + wl_files + lv_files
    print("{} total training files found".format(len(filepaths)))

    # generate label data
    pf_labels = np.repeat(1, len(pf_files))
    pp_labels = np.repeat(2, len(pipspro_files))
    leeds_labels = np.repeat(3, len(leeds_files))
    star_labels = np.repeat(4, len(star_files))
    wl_labels = np.repeat(5, len(wl_files))
    # vmat_labels = np.repeat(6, len(vmat_files))
    lv_labels = np.repeat(6, len(lv_files))
    all_labels = np.concatenate((pf_labels, pp_labels, leeds_labels, star_labels, wl_labels, lv_labels))

    # preallocate
    total_array = np.zeros((len(filepaths), 10000), dtype=np.float32)
    print("Training array preallocated")

    # resize each image and add to a training array
    start = time.time()
    if use_pool:
        futures = {}
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for idx, path in enumerate(filepaths):
                future = exec.submit(process_image, path)
                futures[future] = idx
        for idx, future in enumerate(concurrent.futures.as_completed(futures)):
            total_array[futures[future], :] = future.result()
    else:
        for idx, path in enumerate(filepaths):
            future = process_image(path)
            total_array[idx, :] = future.result()
    print("Training array scaled/processed in {:.2f}s".format(time.time() - start))

    # save arrays to disk for future use
    np.save(osp.join(osp.dirname(osp.abspath(__file__)), 'data', 'images'), total_array)
    np.save(osp.join(osp.dirname(osp.abspath(__file__)), 'data', 'labels'), all_labels)
    print("Images/labels written to disk")


if __name__ == '__main__':
    build_images(use_pool=True)

