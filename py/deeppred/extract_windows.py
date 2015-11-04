'''
A script to extract the window features for a given dataset.
Arguments:
- dataset_name (string): The name of the dataset to extract the window features for.
- advanced (boolean): Whether to use the extra tracks when extracting the windows, or extracting only the simple sequence-based features.
- output_file_path (file path, optional): The path to write the output CSV to. If not provided, will update the project's relevant file.
'''

import sys

import asap

from cleavepred import util
from cleavepred import project_paths
from cleavepred.common import window_extraction_params

### Parse arguments ###

project_paths.dataset_name = sys.argv[1]
advanced = util.parse_bool(sys.argv[2])

if len(sys.argv) > 3:
    output_file_path = sys.argv[3]
else:
    output_file_path = project_paths.get_window_features_file_path(advanced)

### Extract the windows ###

annotated_seqs_file = None
seqs_filtration_file = None
csv_output_file = None
extra_tracks_files = {}

def open_files():

    global annotated_seqs_file, seqs_filtration_file, csv_output_file, extra_tracks_files

    annotated_seqs_file = open(project_paths.get_annotated_seqs_file_path(), 'rb')
    seqs_filtration_file = open(project_paths.get_filtered_seqs_file_path(), 'rb')
    csv_output_file = open(output_file_path, 'wb')

    if advanced:
        for track_name, track_file_path in project_paths.get_track_file_paths().items():
            extra_tracks_files[track_name] = open(track_file_path, 'rb')

def close_files():
    util.close_files([annotated_seqs_file, seqs_filtration_file, csv_output_file])
    util.close_files(extra_tracks_files.values())

def extract_windows():
    asap.extract_windows_from_file(annotated_seqs_file, extract_annotations = True, seqs_filtration_file = seqs_filtration_file, \
            extra_tracks_files = extra_tracks_files, csv_output_file = csv_output_file, window_extraction_params = window_extraction_params)

if __name__ == '__main__':
    try:
        open_files()
        extract_windows()
    finally:
        close_files()
