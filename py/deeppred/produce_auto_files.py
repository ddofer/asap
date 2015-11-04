'''
A master script to produce all the required auto-generated files:
1. Uniprot's annotated seqs .lf file
2. CSVs of the windows with features
3. Pickle dump files of trained predictors.
Just execute this script as it is, with no arguments. Make sure to be in the py/ directory when running it.
'''

import sys
import logging

# For logger initialization
import asap

LOGGING_PREFIX = '********** '

logger = logging.getLogger('EXEC')

# Uniprot's .lf file
logger.info(LOGGING_PREFIX + 'Running extract_uniprot_annotated_seqs_from_xml.py')
sys.argv = ['']
execfile('cleavepred/extract_uniprot_annotated_seqs_from_xml.py')

# Create CSVs
for dataset in ['neuropred', 'uniprot']:
    for advanced in ['false', 'true']:
        logger.info(LOGGING_PREFIX + 'Running extract_windows.py with dataset="%s" and advanced="%s"' % (dataset, advanced))
        sys.argv = ['', dataset, advanced]
        execfile('cleavepred/extract_windows.py')

# # Create dump files
# for advanced in ['false', 'true']:
#     logger.info(LOGGING_PREFIX + 'Running train_classifier.py with advanced="%s"' % advanced)
#     sys.argv = ['', advanced, 'auto']
#     execfile('cleavepred/train_classifier.py')

logger.info(LOGGING_PREFIX + 'Finished executing scripts. All auto-generated files should now be updated.')
