from sklearn import svm, metrics, cross_validation, grid_search
import pickle
import gzip

from pylinac_ml.pylinac_images import load_images

# build_images()
data, labels = load_images()

data_train, data_test, y_train, y_test = cross_validation.train_test_split(data, labels, train_size=0.8)

parameters = {
    'kernel': ['rbf'],
    'C': [10],
    'gamma': [0.001]
}
classifier = grid_search.GridSearchCV(svm.SVC(verbose=True), parameters)
# parameters = {
#     'hidden_layer_sizes': [(2500,)],
#     'activation': ['relu'],
#     'alpha': [0.1, 0.01, 10],
#     'algorithm': ['adam'],
#     'tol': [0.01],
#     'learning_rate': ['invscaling']
# }
# classifier = model_selection.GridSearchCV(neural_network.MLPClassifier(verbose=True), parameters)

classifier.fit(data_train, y_train)

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
# print(metrics.f1_score(y_train, classifier.predict(iris_train), average='binary'))
# print("And test score of")
# print(metrics.f1_score(y_test, classifier.predict(data_test), average='binary'))
with gzip.open('pylinac_model.pkl.gz', mode='wb') as m:
    pickle.dump(classifier, m)

# with gzip.open('pylinac_model.pkl.gz', mode='rb') as m:
#     unp_clf = pickle.load(m)

# print("Classification report after pickling/unpickling:")
# print(metrics.classification_report(y_test, unp_clf.predict(data_test)))