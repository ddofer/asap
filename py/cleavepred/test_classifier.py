'''
A script to test a classifier that was trained on NeuroPred's dataset against UniProt's dataset
Arguments:
- advanced (boolean): Whether to test the advanced or simple classifier.
'''

import sys
import pickle
import logging

import pandas as pd

from cleavepred import util
from cleavepred import project_paths

### Parse arguments ###

advanced = util.parse_bool(sys.argv[1])

### Configuration ###

# We use UniProt's dataset for testing our predictors.
project_paths.dataset_name = 'uniprot'

logger = logging.getLogger('TEST')

### Test the classifier ###

predictor_dump_file = None
windows_file = None

def open_files():
    global predictor_dump_file, windows_file
    predictor_dump_file = open(project_paths.get_peptide_predictor_dump_file_path(advanced), 'rb')
    windows_file = open(project_paths.get_window_features_file_path(advanced), 'rb')
    
def close_files():
    util.close_files([predictor_dump_file, windows_file])

def test_classifier():
    peptide_predictor = pickle.load(predictor_dump_file)
    windows_data_frame = pd.read_csv(windows_file)
    score, roc, sensitivity, precision, specificity, cm = peptide_predictor.window_classifier.test_performance(windows_data_frame, \
            drop_only_almost_positives = True)
    logger.info('score = %f, roc = %f, sensitivity = %f, precision = %f, specificity = %f' % (score, roc, sensitivity, precision, specificity))
    logger.info('Confusion matrix:' + '\n' + str(cm))
    
if __name__ == '__main__':
    try:
        open_files()
        test_classifier()
    finally:
        close_files()
