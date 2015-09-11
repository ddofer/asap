'''
A script to train a classifier from NueroPred's dataset. Will log performance to stdout.
Arguments:
- advanced (boolean): Whether to use the advanced features (extracted with the extra tracks), or just the simple ones.
- predictor_dump_path (file path, optional): The path to dump the trained PeptidePredictor into. If not provided, will not dump it at all.
'''

import sys
import pickle

import pandas as pd

from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC

from . import util
from . import project_paths
from .common import window_extraction_params
from ..asap import train_window_classifier, PeptidePredictor

### Parse arguments ###

advanced = util.parse_bool(sys.argv[1])

if len(sys.argv) > 2:
    predictor_dump_path = sys.argv[2]
else:
    predictor_dump_path = None

### Configuration ###

# We use NeuroPred's dataset for training/validation of our predictors.
project_paths.dataset_name = 'neuropred'

classifiers = [
    LogisticRegressionCV(Cs = 20, n_jobs = -1, class_weight = 'auto'),
    RandomForestClassifier(max_features = 92, n_estimators = 150, n_jobs = -1, class_weight = 'auto'),
    SVC(kernel = 'rbf', C = 3.798, probability = True, cache_size = 1600, class_weight = 'auto'),
]
classifiers = [classifiers[1]] # XXX In the meantime let's just use Random Forest, for the sake of speed

### Train the classifier and dump the predictor ###

windows_file = None
predictor_dump_file = None

def open_files():
    windows_file = open(project_paths.get_window_features_file_path(advanced), 'rb')
    predictor_dump_file = util.open_file(predictor_dump_path, 'wb')
    
def close_files():
    util.close_files([windows_file, predictor_dump_file])
    
def dump_predictor(predictor):
    if predictor_dump_file is not None:
        pickle.dump(predictor, predictor_dump_file)
    
def train_classifier():
    windows_data_frame = pd.read_csv(windows_file)
    window_classifier, classifier_performance = train_window_classifier(windows_data_frame, classifiers = classifiers, \
            drop_only_almost_positives = True)
    peptide_predictor = PeptidePredictor(window_classifier, window_extraction_params = window_extraction_params)
    dump_predictor(peptide_predictor)
    
if __name__ == '__main__':
    try:
        open_files()
        train_classifier()
    finally:
        close_files()
