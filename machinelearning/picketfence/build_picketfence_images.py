import os.path as osp

import matplotlib.pyplot as plt
import numpy as np

from machinelearning.tools import get_files, is_dicom, process_image, drop_non_dicom

data_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'data')


def identify_images(folder, np_name=None, drop_non_dicoms=True):
    if drop_non_dicoms:
        drop_non_dicom(folder)
    print("Obtaining file paths...")
    filepaths = get_files(folder, is_dicom)
    print(len(filepaths), "found")
    feature_array = np.zeros((len(filepaths), 10000), dtype=np.float32)
    labels = np.zeros(len(filepaths))
    for idx, path in enumerate(filepaths):
        # try to automatically identify
        good_names = ('pf', 'vmat', '90', '270', '180', 'ra', 't1', 'picket')
        bad_names = ('open', 't2', 't3', 'speed', 'dose', 'rate', 'drgs', 'mlcs', 'jaw', 'coll', 'strip')
        auto_found = False
        basepath = osp.basename(path).lower()
        for name in good_names:
            if name in basepath:
                label_str = 1
                auto_found = True
                break
        for name in bad_names:
            if name in basepath:
                label_str = 0
                auto_found = True
                break
        if not auto_found:
            img = process_image(path).reshape(100, 100)
            plt.imshow(img, cmap=plt.cm.viridis)
            plt.axis('off')
            plt.title(osp.basename(path))
            plt.show()
            input_invalid = True
            break_out = False
            while input_invalid:
                label_str = input("Input the classification: 1 if a PF field, 0 otherwise. Enter 'done' to finish: ")
                if label_str == 'done':
                    break_out = True
                    input_invalid = False
                try:
                    int(label_str)
                except:
                    pass
                else:
                    if int(label_str) in [0, 1]:
                        input_invalid = False
            if break_out:
                break
        labels[idx] = int(label_str)
        feature_array[idx, :] = process_image(path)
    # trim feature array early if stopped early
    feature_array = feature_array[:idx, :]
    labels = labels[:idx]
    # scaled_features = preprocessing.minmax_scale(feature_array, axis=1)
    scaled_features = feature_array
    if np_name is None:
        np_name = osp.basename(folder)
    np.save(osp.join(data_dir, 'images_' + np_name + '_pf'), scaled_features)
    np.save(osp.join(data_dir, 'labels_' + np_name + '_pf'), labels)


if __name__ == '__main__':
    base_dir = r'C:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files\Picket Fences'
    img_dir = osp.join(base_dir, 'Katy iX')
    identify_images(img_dir, drop_non_dicoms=False)
