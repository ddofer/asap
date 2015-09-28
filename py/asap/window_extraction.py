from StringIO import StringIO
import datetime
import csv
import logging

from . import util
from . import features
from . import parse

BASIC_HEADERS = [
    'peptide_id',
    'window_hot_index',
    'window_seq',
    'window_neighbourhood',
]

ANNOTATION_HEADERS = [
    'window_annotation_mask',
    'window_label',
    'window_only_almost_positive',
]

META_WINDOW_HEADERS = BASIC_HEADERS + ANNOTATION_HEADERS

LOGGER = logging.getLogger('EXTRACTION')

class WindowExtractionParams(object):

    '''
    Parameters that should be used when extracting windows and their features from full records.
    '''

    def __init__(self, window_prefix = 9, window_suffix = 9, neighbourhood_prefix = 5, neighbourhood_suffix = 5, \
            windows_filter = None, feature_keys = features.DEFAULT_FEATURE_KEYS):

        '''
        @param window_prefix, window_suffix (int, both default 7):
            The number of residues before and after the hot index in each window, where the hot index is the position determining the label of
            the window (i.e. the window's label is the value of the annotation mask in the hot index). It follows that the total window size
            is (window_prefix + window_suffix + 1).
        @param neighbourhood_prefix, neighbourhood_suffix (int, both default 5):
            The number of residues before and after the hot index in determining the neighbourhood of the window. The neighbourhood can be used
            during the training process in order to avoid duplicates of very similar windows.
        @param windows_filter (function, optional):
            A function to filter the extracted windows by. The function will receive Window objects and should return a bool stating whether to
            include them or not. If not provided, no filtration will take place in the windows level, and all windows will be extracted.
        @param feature_keys (list, optional, default features.DEFAULT_FEATURE_KEYS):
            A list of features to extract for each window. You can see the full list of optional keywords in features.FEARURE_KEY_OPTIONS, where
            documentation is also provided for each of the feature keys. If not provided, all features will be extracted by default (i.e. will
            use all of the features in features.FEARURE_KEY_OPTIONS). By default will use features.DEFAULT_FEATURE_KEYS, which is another list
            containing most of the features, but not all of them, as there are some features that it doesn't make much sense to use at the same
            time (e.g. both 'aa' and 'aa_reduced'). We anticipate that using the default features should give pretty good results in most
            scenarios, so fine-tuning the exact used features can be left for late stages of a project.
            Important note: Features that rely on extra tracks (ss, acc, disorder, pssm) will not be extracted if the tracks are not provided,
            even if those features are explicitly given in this list.
        '''

        self.window_prefix = window_prefix
        self.window_suffix = window_suffix
        self.neighbourhood_prefix = neighbourhood_prefix
        self.neighbourhood_suffix = neighbourhood_suffix
        self.windows_filter = windows_filter
        self.feature_keys = feature_keys

        self.window_size = window_prefix + window_suffix + 1
        self.window_hot_index = window_prefix
        self.neighbourhood_size = neighbourhood_prefix + neighbourhood_suffix + 1

def extract_windows_from_file(source, extract_annotations = False, seqs_filtration_file = None, \
        extra_tracks_files = {}, csv_output_file = None, window_extraction_params = WindowExtractionParams()):

    '''
        Parses a given file with peptide sequences and breaks it into windows with features, outputs a CSV with a row for each window
        and a column for each feature (along with a few other meta headers).
    @param source (file):
        A file handle to parse the peptide sequences from. Can be either a .fasta or .lf format, depending on the extract_annotations
        parameter.
    @param extract_annotations (boolean, default False):
        Whether to expect finding annotation masks inside the given file of sequences. If set to True, will expect getting a .lf file
        which also contains annotations. If set to False, will expect getting a .fasta file that contains only the sequences. Annotations
        are required only if one plans using the extracted windows to train a new classifier, rather than using an existing one.
    @param seqs_filtration_file (file, optional):
        A fasta format file handle to use for filtering records. If given, will use only sequences with an ID that is also present in the
        ids of this FASTA file. If not given, will not perform any filtration.
    @param extra_tracks_files (dict, empty by default):
        A dictionary for providing extra tracks to extract data from (beside the actual amino-acid sequence and annotations mask). The
        given dictionary should map from a track name to a file handle containing the track data for each of the records in a FASTA
        format. The currently supported extra tracks are: ss (secondary-structure), acc (accessibility), disorder and pssm (position-specific
        scoring matrix).
    @param csv_output_file (file, optional):
        A file handle with writing permissions to write the output CSV into. If not provided, will return a StringIO object from which the
        output CSV can be read.
    @param window_extraction_params (WindowExtractionParams, default params by default):
        Parameters to use for extracting the windows.
    @return:
        If csv_output_file is given, will return nothing. If csv_output_file is not given, will return a a StringIO object from which the
        output CSV can be read.
    '''

    relevant_ids = _get_relevant_ids(seqs_filtration_file)
    extra_tracks = _get_extra_tracks_from_files(extra_tracks_files)
    full_records = list(parse.parse_records_from_file(source, extract_annotations, relevant_ids, extra_tracks))
    parse.LOGGER.info('Final records: %d' % len(full_records))
    _pad_records(full_records, window_extraction_params)
    return _extract_windows(full_records, csv_output_file, window_extraction_params)

def extract_windows_from_seq(seq, annotation_mask = None, extra_tracks_data = {}, csv_output_file = None, \
        window_extraction_params = WindowExtractionParams()):

    '''
    Breaking a peptide sequence into windows with features, outputting a CSV with a row for each window and a column for each feature
    (along with a few other meta headers).
    @param seq (string):
        The peptide sequence to use, given in a 20 amino-acid alphabet.
    @param annotation_mask (string, optional):
        An annotation mask to use as a labeling for each position along the sequence. Expecting a binary sequence (of 0's and 1's) in
        the same length of the given amino-acid sequence. If not provided, the extracted windows won't have labels, meaning they cannot
        be used for training a new classifier (only fed to an already trained classifier).
    @param extra_tracks_data (dict, empty by default):
        A dictionary for providing extra tracks to extract data from (beside the actual amino-acid sequence and annotations mask). The
        given dictionary should map from track names to their sequence. Currently supported extra tracks are: ss (secondary-structure),
        acc (accessibility), disorder and pssm (position-specific scoring matrix).
    @param csv_output_file (file, optional):
        A file handle with writing permissions to write the output CSV into. If not provided, will return a StringIO object from which
        the output CSV can be read.
    @param window_extraction_params (WindowExtractionParams, default params by default):
        Parameters to use for extracting the windows.
    @return:
        If csv_output_file is given, will return nothing. If csv_output_file is not given, will return a a StringIO object from which the
        output CSV can be read.
    '''

    extra_tracks = _get_extra_tracks_from_raw_data(extra_tracks_data)
    full_record = parse.get_record_from_seq(seq, annotation_mask, extra_tracks)
    _pad_record(full_record, window_extraction_params)
    return _extract_windows([full_record], csv_output_file, window_extraction_params)

def _get_relevant_ids(seqs_filtration_file):
    if seqs_filtration_file is None:
        return None
    else:
        relevant_ids = parse.parse_track_from_file(seqs_filtration_file, 'seq').keys()
        parse.LOGGER.info('%d records are in the filtration FASTA file' % len(relevant_ids))
        return relevant_ids

def _get_extra_tracks_from_files(extra_tracks_files):

    extra_tracks = {}

    for track_name, track_source in extra_tracks_files.items():
        if track_name in _TRACK_NAME_TO_TYPE:
            track_type = _TRACK_NAME_TO_TYPE[track_name]
            extra_tracks[track_name] = parse.parse_track_from_file(track_source, track_type)
        else:
            raise Exception('Unknown track name: ' + str(track_name))

    return extra_tracks

def _get_extra_tracks_from_raw_data(extra_tracks_data):

    extra_tracks = {}

    for track_name, raw_track_seq in extra_tracks_data.items():
        if track_name in _TRACK_NAME_TO_TYPE:
            track_type = _TRACK_NAME_TO_TYPE[track_name]
            extra_tracks[track_name] = parse.parse_track_from_seq(raw_track_seq, track_type)
        else:
            raise Exception('Unknown track name: ' + str(track_name))

    return extra_tracks

def _pad_records(records, window_extraction_params):
    for record in records:
        _pad_record(record, window_extraction_params)

def _pad_record(record, window_extraction_params):
    record.pad(window_extraction_params.window_prefix - 1, window_extraction_params.window_suffix - 1)

def _extract_windows(full_records, csv_output_file, window_extraction_params):
    if csv_output_file is None:
        csv_buffer = StringIO()
        _extract_windows_to_csv(full_records, csv_buffer, window_extraction_params)
        csv_buffer.seek(0)
        return csv_buffer
    else:
        _extract_windows_to_csv(full_records, csv_output_file, window_extraction_params)

def _extract_windows_to_csv(full_records, output_file, window_extraction_params):

    LOGGER.info('Extracting windows with features in CSV format...')
    start = datetime.datetime.now()
    csv_writer = csv.writer(output_file)
    feature_headers = None
    include_annotations = None

    for record in full_records:
        for window in record.get_windows(window_extraction_params.window_size):
            feature_headers, include_annotations = _process_window_to_csv(window, csv_writer, window_extraction_params, \
                        feature_headers, include_annotations)

    time_diff = datetime.datetime.now() - start
    LOGGER.info('Done. Extraction took %d seconds.' % time_diff.total_seconds())

def _process_window_to_csv(window, csv_writer, window_extraction_params, feature_headers, include_annotations):

    if feature_headers is None:
        features = window.get_features(window_extraction_params.window_hot_index, window_extraction_params.feature_keys)
        feature_headers = list(sorted(features.keys()))
        include_annotations = window.has_annotation_mask()
        util.write_csv_line(csv_writer, _get_meta_headers(include_annotations) + feature_headers)

    if window_extraction_params.windows_filter is None or window_extraction_params.windows_filter(window):
        meta_values = _get_window_meta_values(window, window_extraction_params, include_annotations)
        features = window.get_features(window_extraction_params.window_hot_index, window_extraction_params.feature_keys)
        feature_values = [features[header] for header in feature_headers]
        util.write_csv_line(csv_writer, meta_values + feature_values)

    return feature_headers, include_annotations

def _get_meta_headers(include_annotations):
    if include_annotations:
        return BASIC_HEADERS + ANNOTATION_HEADERS
    else:
        return BASIC_HEADERS

def _get_window_meta_values(window, window_extraction_params, include_annotations):

    hot_index = window.original_index + window_extraction_params.window_hot_index
    neighbourhood = window.get_neighbourhood(window_extraction_params.window_hot_index, window_extraction_params.neighbourhood_prefix, \
            window_extraction_params.neighbourhood_suffix)
    meta_values = [window.full_record.id, hot_index, window.get_aa_seq(), neighbourhood]

    if include_annotations:
        label = window.get_label(window_extraction_params.window_hot_index)
        is_only_almost_positive = window.is_only_almost_positive(window_extraction_params.window_hot_index)
        meta_values += [window.get_annotation_mask(), label, is_only_almost_positive]

    return meta_values

_TRACK_NAME_TO_TYPE = {
    'ss': 'seq',
    'acc': 'seq',
    'disorder': 'disorder',
    'pssm': 'pssm',
}
