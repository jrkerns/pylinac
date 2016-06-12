import concurrent.futures
import os
import gzip
import pickle
import os.path as osp
import time

import dicom
import numpy as np
from pylinac import image
from scipy.misc import imresize
from sklearn import svm, metrics, cross_validation, grid_search, preprocessing


def is_dicom(path):
    """Whether the file is a readable DICOM file via pydicom."""
    try:
        ds = dicom.read_file(path, force=True)
        ds.pixel_array
        return True
    except:
        return False


def drop_non_dicom(folder):
    """Remove all files within a folder that are not DICOM images. Space-saving utility function."""
    print("Dropping non-DICOM files...")
    for pdir, _, files in os.walk(folder):
        for file in files:
            file = osp.join(pdir, file)
            if not is_dicom(file):
                os.remove(file)


def get_files(folder, func, use_pool=False):
    """Get a list of files that are valid images from the folder."""
    if not osp.isdir(folder):
        raise NotADirectoryError("{} is not a directory".format(folder))
    print("Grabbing file names...")
    if not use_pool:
        paths = []
        for pdir, _, files in os.walk(folder):
            for file in files:
                filepath = osp.join(pdir, file)
                if func(filepath):
                    paths.append(filepath)
        return paths
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
        return filepaths


def load_images(path):
    """Load the built images for training."""
    imgs = get_files(path, lambda x: 'images' in x)
    img_arr = np.vstack([np.load(f) for f in imgs])
    labels = get_files(path, lambda x: 'labels' in x)
    labels_arr = np.concatenate([np.load(f) for f in labels])
    return img_arr, labels_arr


def process_image(path):
    """Load and resize the images and return as flattened numpy array"""
    img = image.load(path, dtype=np.float32)
    resized_img = imresize(img.array, size=(100, 100), mode='F').flatten()
    rescaled_img = preprocessing.minmax_scale(resized_img, axis=1)
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

    data_train, data_test, y_train, y_test = cross_validation.train_test_split(data, labels, train_size=train_size)
    start = time.time()
    classifier = grid_search.GridSearchCV(svm.SVC(verbose=True), parameters)
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
    print("Classification report:")
    print(metrics.classification_report(y_test, classifier.predict(data_test)))
    with gzip.open(clf_name + '_classifier.pkl.gz', mode='wb') as m:
        pickle.dump(classifier, m)


def strip(folder, classifier_prefix, correct_prediction):
    """Strip a folder of non-DICOM files and of image files that are not of the predicted classification."""
    drop_non_dicom(folder)
    filepaths = get_files(folder, lambda x: True)
    # load classifier
    folder = osp.join(osp.dirname(__file__), classifier_prefix)
    with gzip.open(osp.join(folder, classifier_prefix + '_classifier.pkl.gz'), mode='rb') as m:
        clf = pickle.load(m)
    for file in filepaths:
        img = process_image(file)
        img = preprocessing.minmax_scale(img.reshape(1, -1), axis=1)
        prediction = clf.predict(img)
        print("Prediction {} for file: {}".format(prediction, osp.basename(file)))
        time.sleep(0.3)
        if prediction == correct_prediction:
            os.remove(file)
    print("Done stripping")


if __name__ == '__main__':
    pass
