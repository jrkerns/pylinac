import os.path as osp

from machinelearning.tools import train

path = osp.join(osp.dirname(__file__), 'data', 'CatPhan 600')

parameters = {
    'kernel': ['linear'],
    'C': [1, 0.1, 5, 10, 50],
}

train(path, train_size=0.95, parameters=parameters, clf_name='catphan600')
