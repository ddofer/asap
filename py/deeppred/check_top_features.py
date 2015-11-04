'''
Checks the top features predicted for a given dataset.
Arguments:
- dataset_name (string): The name of the dataset to use.
- advanced (boolean): Whether to use all or only some features, corresponding to simple and advanced classifiers respectively.
'''

import sys
import logging

import pandas as pd

from asap import get_top_features

from cleavepred import util
from cleavepred import project_paths

logger = logging.getLogger('FEATURES')

### Parse arguments ###

project_paths.dataset_name = sys.argv[1].lower()
advanced = util.parse_bool(sys.argv[2])

### Get top features ###

windows_file = None

def open_files():
    global windows_file
    windows_file = open(project_paths.get_window_features_file_path(advanced), 'rb')
    
def close_files():
    util.close_files([windows_file])
    
def get_advanced_label():
    if advanced:
        return 'advanced'
    else:
        return 'simple'
    
def check_top_features():
    windows_data_frame = pd.read_csv(windows_file)
    logger.info('Checking top features over %s dataset with %s features...' % (project_paths.dataset_name, get_advanced_label()))
    top_features = get_top_features(windows_data_frame, drop_only_almost_positives = True)
    logger.info('Top features: ' + ', '.join(top_features))
    
if __name__ == '__main__':
    try:
        open_files()
        check_top_features()
    finally:
        close_files()
