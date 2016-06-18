import os.path as osp

from machinelearning.tools import train

path = osp.join(osp.dirname(__file__), 'data')

parameters = {
    'kernel': ['linear'],
    'C': [1, 5],
    # 'gamma': [0.001, 0.01, 0.1]
}

train(path, train_size=0.85, parameters=parameters, clf_name='singleimage')
