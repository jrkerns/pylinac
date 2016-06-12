import math
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from pylinac.core.io import TemporaryZipDirectory

from machinelearning.tools import get_files, is_dicom, process_image


def identify_images(zip_file):
    """Interactively identify images from a folder, writing the labels to an array for later training"""
    with TemporaryZipDirectory(zip_file) as zfiles:
        filepaths = get_files(zfiles, is_dicom)
        feature_array = np.zeros((len(filepaths), 10000), dtype=np.float32)
        labels = np.zeros(len(filepaths))
        split_val = 25
        length = len(filepaths)
        rounds = int(math.ceil(length / split_val))
        for n in range(rounds):
            fig, axes = plt.subplots(5, 5)
            for axis, (idx, fp) in zip(axes.flatten(), enumerate(filepaths[split_val*n:split_val*(n+1)])):
                img = process_image(fp)
                plt.sca(axis)
                plt.imshow(img.array, cmap=plt.cm.Greys)
                plt.axis('off')
                plt.title(idx+split_val*n)
            plt.show()
        not_done = True
        while not_done:
            label = input("Input the HU indices sequentially, one at a time. Type 'done' when finished:")
            if label == 'done':
                not_done = False
            else:
                labels[int(label)] = 1
        # labels = np.array(labels)
        for idx, fp in enumerate(filepaths):
            feature_array[idx, :] = process_image(fp)
        # scaled_features = preprocessing.minmax_scale(feature_array, axis=1)
        scaled_features = feature_array
        dir2write = osp.dirname(zip_file)
        np.save(osp.join(dir2write, 'images_' + osp.splitext(osp.basename(zip_file))[0]), scaled_features)
        np.save(osp.join(dir2write, 'labels_' + osp.splitext(osp.basename(zip_file))[0]), labels)


if __name__ == '__main__':
    data_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'data')
    zsets = (
        # osp.join(data_dir, 'Elekta_7.zip'),
        # osp.join(data_dir, 'Elekta_8.zip'),
        # osp.join(data_dir, 'Elekta_11.zip'),
        # osp.join(data_dir, 'Elekta_12.zip'),
        osp.join(data_dir, 'CBCT_17.zip'),
        # osp.join(data_dir, 'Standard head.zip'),
    )
    for zset in zsets:
        # path = osp.join(osp.dirname(osp.abspath(__file__)), 'data', 'thorax.zip')
        identify_images(zset)
