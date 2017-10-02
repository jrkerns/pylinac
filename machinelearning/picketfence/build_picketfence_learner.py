import os.path as osp

from machinelearning.tools import train

path = osp.join(osp.dirname(__file__), 'data')

parameters = {
    'kernel': ['linear'],
    'C': [1, 0.5, 5],
}

train(path, train_size=0.98, parameters=parameters, clf_name='picketfence')
