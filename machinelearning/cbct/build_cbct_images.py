import math
import os
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from pylinac.core.io import TemporaryZipDirectory
from pylinac.core.image import load

from machinelearning.tools import get_files, is_dicom, process_image


def identify_images(zip_file):
    """Interactively identify images from a folder, writing the labels to an array for later training"""
    with TemporaryZipDirectory(zip_file) as zfiles:
        filepaths = get_files(zfiles, is_dicom)
        labels = np.zeros(len(filepaths))
        split_val = 25
        length = len(filepaths)
        rounds = int(math.ceil(length / split_val))
        for n in range(rounds):
            fig, axes = plt.subplots(5, 5, figsize=(10, 10))
            for axis, (idx, fp) in zip(axes.flatten(), enumerate(filepaths[split_val*n:split_val*(n+1)])):
                img = load(fp)
                plt.sca(axis)
                plt.imshow(img, cmap=plt.cm.Greys)
                plt.axis('off')
                plt.title(idx+split_val*n)
            plt.show()
        not_done = True
        while not_done:
            label = input("Input the HU indices as a number or range. E.g. '66' or '25-47'. Type 'done' when finished:")
            if label == 'done':
                not_done = False
            else:
                items = label.split('-')
                if len(items) > 1:
                    labels[int(items[0]):int(items[1])] = 1
                else:
                    labels[int(items[0])] = 1
        scaled_features = np.zeros((len(filepaths), 10000), dtype=np.float32)
        for idx, fp in enumerate(filepaths):
            scaled_features[idx, :] = process_image(fp)
    dir2write = osp.dirname(zip_file)
    np.save(osp.join(dir2write, 'images_' + osp.splitext(osp.basename(zip_file))[0]), scaled_features)
    np.save(osp.join(dir2write, 'labels_' + osp.splitext(osp.basename(zip_file))[0]), labels)
    os.remove(zip_file)


if __name__ == '__main__':
    data_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'data', 'CatPhan 600')
    zsets = get_files(data_dir, func=lambda x: x.endswith('.zip'))
    for zset in zsets:
        identify_images(zset)
