import os.path as osp
import os
import math

import dicom
import matplotlib.pyplot as plt
import numpy as np
from pylinac import image
from pylinac.core.io import TemporaryZipDirectory
from scipy.misc import imresize
from sklearn import preprocessing


def is_dicom(path):
    """Whether the file is a readable DICOM file via pydicom."""
    try:
        ds = dicom.read_file(path, force=True)
        ds.pixel_array
        return True
    except:
        return False


def get_files(folder, func):
    """Get a list of files that are valid images from the folder."""
    paths = []
    for pdir, _, files in os.walk(folder):
        for file in files:
            filepath = osp.join(pdir, file)
            if func(filepath):
                paths.append(filepath)
    return paths


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
                img = image.load(fp)
                plt.sca(axis)
                plt.imshow(img.array, cmap=plt.cm.Greys)
                plt.axis('off')
                plt.title(idx+split_val*n)
            plt.show()
        # for idx, fp in enumerate(filepaths):
        #     img = image.load(fp)
        #     img.plot()
        #     label = input("Input 0 or nothing if not an HU slice, 1 if it is:")
        #     if label == '':
        #         label = 0
        #     else:
        #         label = 1
        #     labels.append(label)
        #     feature_array[idx, :] = process_image(fp)
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
        scaled_features = preprocessing.minmax_scale(feature_array, axis=1)
        dir2write = osp.dirname(zip_file)
        np.save(osp.join(dir2write, 'images_' + osp.splitext(osp.basename(zip_file))[0]), scaled_features)
        np.save(osp.join(dir2write, 'labels_' + osp.splitext(osp.basename(zip_file))[0]), labels)


def process_image(path):
    """Load and resize the images and return as flattened numpy array"""
    img = image.load(path, dtype=np.float32)
    return imresize(img.array, size=(100, 100), mode='F').flatten()


def load_images():
    """Load the built images for training."""
    path = osp.join(osp.dirname(osp.abspath(__file__)), 'data')
    imgs = get_files(path, lambda x: 'images' in x)
    img_arr = np.vstack([np.load(f) for f in imgs])
    labels = get_files(path, lambda x: 'labels' in x)
    labels_arr = np.concatenate([np.load(f) for f in labels])
    return img_arr, labels_arr


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
