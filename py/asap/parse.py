import logging

from Bio import SeqIO

from . import util
from . import config
from . import data

LOGGER = logging.getLogger('PARSE')

def convert_lf_to_fasta(source, output_file):

    '''
    Converts a .lf file, which also contains annotations, to a .fasta file, which contains only the amino-acid sequences
    of the records.
    @param source (file handle):
        The source .lf file to read.
    @param output_file (file handle):
        A file handlw with writing permissions to write the output FASTA into.
    '''

    full_records = parse_records_from_file(source, extract_annotations = True)
    fasta_records = [full_record.to_fasta_record() for full_record in full_records]
    SeqIO.write(fasta_records, output_file, 'fasta')

def parse_records_from_file(source, extract_annotations, relevant_ids = None, extra_tracks = {}):

    '''
    Parses full data records from a file (either .lf format, which also includes annotation masks,
        or a simple fasta )
    @param source (file):
        A file handle to parse the records from.
    @param extract_annotations (bool):
        Whether to expect a .lf format which also contains annotations, or a simple .fasta format.
    @param relevant_ids (collection, optional):
        An optional list of ids. If None, will do nothing. If provided with a collection,
    will return only records with the given ids.
    @param extra_tracks (dict, empty by default):
        Extra tracks to give the records, given in the following format:
        {
            track_name : {
                record_id: (seq, padding_value),
                ...
            }
            ...
        }
    @return:
        A generator for the parsed records (each of type FullDataRecord).
    '''

    seqs = list(SeqIO.parse(source, 'fasta'))
    LOGGER.info('Parsing %d sequencess...' % len(seqs))

    for seq in seqs:
        if relevant_ids is None or _format_id(seq.id) in relevant_ids:
            yield _parse_fasta_record(seq, extra_tracks, extract_annotations)

def get_record_from_seq(seq, annotation_mask = None, extra_tracks = {}):

    '''
    Creates a full data record from sequences.
    @param seq (string):
        The amino-acid sequence of the record
    @param annotation_mask (string, optional):
        A binary mask (made of 0's and 1's) in the same size of the given sequence to use
        as an annotation mask. If not provided, the record will not have an annotation mask.
    @param extra_tracks (dict, empty by default):
    Extra tracks to give the record, given in the following format:
        {
            track_name: (seq, padding_value),
            ...
        }
    @return:
        A FullDataRecord created from the provided data.
    '''

    sequence_tracks = data.SequenceTracks()
    sequence_tracks.add_track(data.SequenceTrack('aa', seq))

    if annotation_mask is not None:
        sequence_tracks.add_track(data.SequenceTrack('annotation', annotation_mask))

    for track_name, track_data in extra_tracks.items():
        track_seq, track_padding_value = track_data
        sequence_tracks.add_track(data.SequenceTrack(track_name, track_seq, track_padding_value))

    return data.FullDataRecord('N/A', 'N/A', 'N/A', sequence_tracks)

def parse_track_from_file(source, type):

    '''
    Parses the track data of multiple records from a FASTA file .
    @param source (file):
        The file handle to parse (in FASTA format)
    @param type (string):
        The type of the track to parse (options: seq, disorder, pssm)
    @return:
        A dictionary of the following format:
        {
            record_id: (seq, padding_value),
            ...
        }
    '''

    track_file_parser, track_seq_parser, padding_value = _TRACK_TYPE_TO_PARSERS_AND_PADDING[type]
    track_data = {}

    for record_id, seq in track_file_parser(source):
        track_data[_format_id(record_id)] = (seq, padding_value)

    return track_data

def parse_track_from_seq(seq, type):

    '''
    Parses the track data of a single record from a raw sequence.
    @param seq (string):
        The raw sequence to parse
    @param type (string):
        The type of the track to parse (options: seq, disorder, pssm)
    @return:
        A tuple containing the parsed track sequence and its padding value.
    '''

    track_file_parser, track_seq_parser, padding_value = _TRACK_TYPE_TO_PARSERS_AND_PADDING[type]
    return track_seq_parser(seq), padding_value

def _parse_fasta_record(fasta_seq, extra_tracks, extract_annotations):

    record_id = _format_id(fasta_seq.id)

    if extract_annotations:
        raw_seq_and_mask = str(fasta_seq.seq)
        mask_start_index = util.find_first_index_of(raw_seq_and_mask, '01')
        aa_seq = _fix_aa_seq(raw_seq_and_mask[:mask_start_index])
        annotation_mask = raw_seq_and_mask[mask_start_index:]
    else:
        aa_seq = fasta_seq.seq
        annotation_mask = None

    sequence_tracks = data.SequenceTracks()
    sequence_tracks.add_track(data.SequenceTrack('aa', aa_seq))

    if annotation_mask is not None:
        sequence_tracks.add_track(data.SequenceTrack('annotation', annotation_mask))

    for track_name, extra_track_data in extra_tracks.items():
        if record_id in extra_track_data:
            track_seq, track_padding_value = extra_track_data[record_id]
            sequence_tracks.add_track(data.SequenceTrack(track_name, track_seq, track_padding_value))
        else:
            raise Exception('No record for %s in track %s' % (record_id, track_name))

    return data.FullDataRecord(record_id, fasta_seq.name, fasta_seq.description, sequence_tracks)

def _parse_seq_track_from_file(source):
    for seq in SeqIO.parse(source, 'fasta'):
        yield seq.id, str(seq.seq)

def _parse_seq_track_from_seq(seq):
    return seq

def _parse_disorder_track_from_file(source):
    for seq in SeqIO.parse(source, 'fasta'):
        yield seq.id, _parse_disorder_track_from_seq(str(seq.seq))

def _parse_disorder_track_from_seq(seq):
    disorder_start_index = util.find_first_index_of(seq, config.DISORDER_OPTIONS)
    return seq[disorder_start_index:]

def _parse_pssm_track_from_file(source):
    for raw_record in source.read().split('>')[1:]:
        lines = raw_record.splitlines()
        record_id = lines[0].split(' ')[0]
        pssm = _parse_pssm(lines[1:])
        yield record_id, pssm

def _parse_pssm_track_from_seq(seq):
    return _parse_pssm(seq.splitlines())

def _parse_pssm(lines):

    pssm = []

    for line in lines:
        freqs_vector = map(float, line.split(' ')[1:])
        freqs_dict = dict(zip(config.PSSM_AMINO_ACIDS, freqs_vector))
        freqs_dict['_'] = 0.0
        pssm += [freqs_dict]

    return pssm

def _fix_aa_seq(seq):

    fixed_seq = ''

    for aa in seq:
        if aa in config.AMINO_ACIDS:
            fixed_seq += aa
        else:
            fixed_seq += '_'

    return fixed_seq

def _format_id(id):
    return id.replace('|', '__').replace('-', '_')

_TRACK_TYPE_TO_PARSERS_AND_PADDING = {
    'seq': (_parse_seq_track_from_file, _parse_seq_track_from_seq, '_'),
    'disorder': (_parse_disorder_track_from_file, _parse_disorder_track_from_seq, '_'),
    'pssm': (_parse_pssm_track_from_file, _parse_pssm_track_from_seq, [dict([(aa, 0.0) for aa in config.PSSM_AMINO_ACIDS] + [('_', 1.0)])]),
}
