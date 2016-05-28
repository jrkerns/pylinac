import gzip
import pickle
import time

from sklearn import svm, metrics, cross_validation, grid_search

from machinelearning.cbct.build_cbct_images import load_images

data, labels = load_images()

data_train, data_test, y_train, y_test = cross_validation.train_test_split(data, labels, train_size=0.85)

parameters = {
    'kernel': ['linear'],
    'C': [1, 0.1, 5],
}
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
with gzip.open('cbct_classifier.pkl.gz', mode='wb') as m:
    pickle.dump(classifier, m)
