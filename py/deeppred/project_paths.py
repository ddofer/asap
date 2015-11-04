import os

from .common import AVAILABLE_TRACKS

# A global variable to update whenever looking to work on another dataset
dataset_name = None

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# DATA_DIR = os.path.join(BASE_DIR, 'data/cleavage')
DATA_DIR = os.path.join(BASE_DIR, 'data/deep/cleavage')

def get_dataset_dir():
    return os.path.join(DATA_DIR, '%s_dataset' % dataset_name)

def get_peptide_predictor_dump_file_path(advanced):
    if advanced:
        return os.path.join(DATA_DIR, 'advanced_peptide_predictor.pkl')
    else:
        return os.path.join(DATA_DIR, 'simple_peptide_predictor.pkl')

def get_window_features_file_path(advanced):
    if advanced:
        return os.path.join(get_dataset_dir(), 'window_advanced_features.csv')
    else:
        return os.path.join(get_dataset_dir(), 'window_simple_features.csv')

def get_raw_data_xml_file_path():
    # Relevant only for when dataset_name = 'uniprot'
    return os.path.join(get_dataset_dir(), 'raw_data.xml')

def get_annotated_seqs_file_path():
    return os.path.join(get_dataset_dir(), 'annotated_seqs.lf')

def get_filtered_seqs_file_path():
    return os.path.join(get_dataset_dir(), 'filtered_seqs.fasta')

def get_track_file_paths():
    return {track: os.path.join(get_dataset_dir(), 'extra_tracks/seqs.%s' % track) for track  in AVAILABLE_TRACKS}