import os.path as osp

import matplotlib.pyplot as plt
import numpy as np

from machinelearning.cbct.build_cbct_images import get_files, is_dicom, process_image

data_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'data')


def identify_images(folder, name=None):
    print("Obtaining file paths...")
    filepaths = get_files(folder, is_dicom)
    print(len(filepaths), "found")
    feature_array = np.zeros((len(filepaths), 10000), dtype=np.float32)
    labels = np.zeros(len(filepaths))
    for idx, path in enumerate(filepaths):
        img = process_image(path)
        plt.imshow(img, cmap=plt.cm.viridis)
        plt.axis('off')
        plt.show()
        input_invalid = True
        break_out = False
        while input_invalid:
            label_str = input("Input the classification: 0 if not a VMAT field, 1 if open field, 2 if drgs, 3 if mlcs. 'done' to finish: ")
            if label_str == 'done':
                break_out = True
                input_invalid = False
            elif int(label_str) in [0, 1, 2, 3]:
                input_invalid = False
        if break_out:
            break
        labels[idx] = int(label_str)
        feature_array[idx, :] = process_image(path)
    # scaled_features = preprocessing.minmax_scale(feature_array, axis=1)
    scaled_features = feature_array
    if name is None:
        name = osp.basename(folder)
    np.save(osp.join(data_dir, 'images_' + name + '_vmat'), scaled_features)
    np.save(osp.join(data_dir, 'labels_' + name + '_vmat'), labels)


if __name__ == '__main__':
    base_dir = r'C:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files\VMATs'
    img_dir = osp.join(base_dir, 'ACB 1')
    identify_images(img_dir)
