import concurrent.futures
import os
import gzip
import pickle
import os.path as osp
import random
import time

import pydicom
import numpy as np
from pylinac import image
from scipy.misc import imresize
from sklearn import svm, metrics, preprocessing, model_selection


def is_dicom(path):
    """Whether the file is a readable DICOM file via pydicom."""
    try:
        ds = pydicom.dcmread(path, force=True)
        ds.pixel_array
        return True
    except:
        return False


def drop_non_dicom(folder, use_pool=True):
    """Remove all files within a folder that are not DICOM images. Space-saving utility function."""
    print("Dropping non-DICOM files...")
    if not use_pool:
        for pdir, _, files in os.walk(folder):
            for file in files:
                file = osp.join(pdir, file)
                if not is_dicom(file):
                    os.remove(file)
                    print("Deleting", file)
    else:
        futures = {}
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, _, files in os.walk(folder):
                for file in files:
                    filepath = osp.join(pdir, file)
                    future = exec.submit(is_dicom, filepath)
                    futures[future] = filepath
        print("Queued {} file identifications".format(len(futures)))
        # filepaths = []
        for idx, future in enumerate(concurrent.futures.as_completed(futures)):
            if not future.result():
                os.remove(futures[future])
                # filepaths.append(futures[future])
        print("Done identifying files in {} in {:.2f}s".format(osp.basename(folder), time.time() - start))


def get_files(folder, func, use_pool=False, randomize=False, recursive=True):
    """Get a list of files that are valid images from the folder."""
    if not osp.isdir(folder):
        raise NotADirectoryError("{} is not a directory".format(folder))
    print("Grabbing file names...")
    # get filenames
    all_files = []
    if recursive:
        for pdir, _, files in os.walk(folder):
            for file in files:
                filepath = osp.join(pdir, file)
                all_files.append(filepath)
    else:
        files = os.listdir(folder)
        for file in files:
            filepath = osp.join(folder, file)
            if osp.isfile(filepath):
                all_files.append(filepath)

    if not use_pool:
        filepaths = []
        for file in all_files:
            if func(file):
                filepaths.append(file)
    else:
        futures = {}
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor() as exec:
            for pdir, _, files in os.walk(folder):
                for file in files:
                    filepath = osp.join(pdir, file)
                    future = exec.submit(func, filepath)
                    futures[future] = filepath
        print("Queued {} file identifications".format(len(futures)))
        filepaths = []
        for idx, future in enumerate(concurrent.futures.as_completed(futures)):
            if future.result():
                filepaths.append(futures[future])
        print("Done identifying files in {} in {:.2f}s".format(osp.basename(folder), time.time() - start))
    if randomize:
        random.shuffle(filepaths)
    return filepaths


def load_images(path):
    """Load the built images for training."""
    imgs = get_files(path, lambda x: 'images' in x, recursive=False)
    img_arr = np.vstack([np.load(f) for f in imgs])
    labels = get_files(path, lambda x: 'labels' in x, recursive=False)
    labels_arr = np.concatenate([np.load(f) for f in labels])
    return img_arr, labels_arr


def process_image(path):
    """Load and resize the images and return as flattened numpy array"""
    img = image.load(path, dtype=np.float32)
    resized_img = imresize(img.array, size=(100, 100), mode='F').flatten()
    rescaled_img = preprocessing.minmax_scale(resized_img)
    return rescaled_img


def train(path, train_size, parameters, clf_name):
    """Train an SVM classifier on a set of labeled images.

    Parameters
    ----------
    path : str
        Path to the folder containing the images and labels as numpy array files.
    train_size : float
        Training size proportion of input images. Must be between 0 and 1.
    parameters : dict
        Set of parameters to pass to the SMV grid search algorithm.
    clf_name : str
        Prefix name of classifier; e.g. 'vmat', 'cbct'.
    """
    data, labels = load_images(path)

    data_train, data_test, y_train, y_test = model_selection.train_test_split(data, labels, train_size=train_size)
    start = time.time()
    classifier = model_selection.GridSearchCV(svm.SVC(verbose=True), parameters)
    classifier.fit(data_train, y_train)
    print()
    print("Training took: {:.2f}s".format(time.time() - start))

    for params, mean_score, scores in classifier.grid_scores_:
        print("%0.3f (+/-%0.03f) for %r"
              % (mean_score, scores.std() * 2, params))
    print()
    print(classifier.best_estimator_)
    print("Best parameters found:")
    print(classifier.best_params_)
    print("With a training score of:")
    print(classifier.best_score_)
    print()
    print("Test data classification report:")
    print(metrics.classification_report(y_test, classifier.predict(data_test)))
    with gzip.open(clf_name + '_classifier.pkl.gz', mode='wb') as m:
        pickle.dump(classifier, m)


def strip(folder, classifier_prefix, correct_prediction, correct_names=None, incorrect_names=None, drop_non_dicoms=False):
    """Strip a folder of non-DICOM files and of image files that are not of the predicted classification."""
    if drop_non_dicoms:
        drop_non_dicom(folder)
    filepaths = get_files(folder, lambda x: True, randomize=False)
    # load classifier
    folder = osp.join(osp.dirname(__file__), classifier_prefix)
    with gzip.open(osp.join(folder, classifier_prefix + '_classifier.pkl.gz'), mode='rb') as m:
        clf = pickle.load(m)
    files2delete = []
    for file in filepaths:
        incorrect = True
        if incorrect_names is not None:
            for name in incorrect_names:
                if name in osp.basename(file).lower():
                    files2delete.append(file)
                    incorrect = False
                    break
        if incorrect:
            classify = True
            if correct_names is not None:
                for name in correct_names:
                    if name in osp.basename(file).lower():
                        classify = False
                        break
            if classify:
                img = process_image(file)
                prediction = clf.predict(img.reshape(1, -1))
                print("Prediction {} for file: {}".format(prediction, file))
                time.sleep(0.3)
                if prediction not in correct_prediction:
                    files2delete.append(file)
    for file in files2delete:
        os.remove(file)
    print("Done stripping")


if __name__ == '__main__':
    pass
    path = r'C:\Users\James\Dropbox\Programming\Python\Projects\pylinac test files\Picket Fences'
    drop_non_dicom(path)
