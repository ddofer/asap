from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from . import features

class SequenceTrack(object):

    def __init__(self, name, seq, padding_value = '_'):
        '''
        @param name:
            The name of the track (string).
        @param seq:
            The actual sequence (string).
        '''
        self.name = name
        self.seq = seq
        self.padding_value = padding_value

    def length(self):
        return len(self.seq)

    def get_subsequence(self, start, length):
        return SequenceTrack(self.name, self.seq[start:(start + length)])

    def pad(self, prefix_length, suffix_length):
        prefix = prefix_length * self.padding_value
        suffix = suffix_length * self.padding_value
        self.seq = prefix + self.seq + suffix

    def __repr__(self):
        return '%s: %s' % (self.name, self.seq)

class SequenceTracks(object):

    def __init__(self):
        self.tracks = {}

    def add_track(self, sequence_track):
        self.tracks[sequence_track.name] = sequence_track

    def get_track(self, name):
        return self.tracks[name]

    def get_subsequence(self, start, length):

        subsequence = SequenceTracks()

        for track in self.tracks.values():
            subsequence.add_track(track.get_subsequence(start, length))

        return subsequence

    def pad(self, prefix_length, suffix_length):
        for track in self.tracks.values():
            track.pad(prefix_length, suffix_length)

    def length(self):

        used_track_name = None
        length = None

        for track in self.tracks.values():
            if length is None:
                used_track_name = track.name
                length = track.length()
            elif length != track.length():
                raise Exception('Track lengths don\'t match (%s: %d, %s: %d)' % (used_track_name, length, track.name, track.length()))

        return length

class DataRecord(object):

    def __init__(self, sequence_tracks):
        '''
        @param sequence_tracks (SequenceTracks):
            All the sequence tracks used for this data record.
        '''
        self.sequence_tracks = sequence_tracks
        self.padding_prefix_length = 0
        self.padding_suffix_length = 0
        
    def length(self):
        return self.sequence_tracks.length()

    def pad(self, prefix_length, suffix_length):
        self.padding_prefix_length += prefix_length
        self.padding_suffix_length += suffix_length
        self.sequence_tracks.pad(prefix_length, suffix_length)

    def get_track_seq(self, name):
        return self.sequence_tracks.get_track(name).seq

    def get_aa_seq(self):
        return self.get_track_seq('aa')

    def get_annotation_mask(self):
        return self.get_track_seq('annotation')

    def get_available_tracks(self):
        return self.sequence_tracks.tracks.keys()

    def has_annotation_mask(self):
        return 'annotation' in self.get_available_tracks()

class FullDataRecord(DataRecord):

    def __init__(self, id, name, description, sequence_tracks):
        '''
        @see DataRecord
        @param id (string):
            The record's ID from FASTA
        @param name (string):
            The record's name from FASTA
        @param description (string):
            The record's description from FASTA
        '''
        DataRecord.__init__(self, sequence_tracks)
        self.id = id
        self.name = name
        self.description = description

    def to_fasta_record(self):
        return SeqRecord(Seq(self.get_aa_seq(), IUPAC.protein),
                         id = self.id,
                         name = self.name,
                         description = self.description)

    def get_windows(self, window_size):
        for i in range(self.length() - window_size + 1):
            yield Window(self, i, i - self.padding_prefix_length, self.sequence_tracks.get_subsequence(i, window_size))

    def __repr__(self):
        return 'Record %s' % self.id

class Window(DataRecord):

    def __init__(self, full_record, offset, original_index, sequence_tracks):
        DataRecord.__init__(self, sequence_tracks)
        self.full_record = full_record
        self.offset = offset
        self.original_index = original_index

    def get_left_context_track_seq(self, name):
        return self.full_record.get_track_seq(name)[:self.offset]

    def get_right_context_track_seq(self, name):
        return self.full_record.get_track_seq(name)[(self.offset + self.length()):]

    def get_neighbourhood(self, hot_index, neighbourhood_prefix, neighbourhood_suffix):
        return self.get_aa_seq()[(hot_index - neighbourhood_prefix):(hot_index + neighbourhood_suffix + 1)]

    def get_label(self, hot_index):
        return self.get_annotation_mask()[hot_index] == '1'

    def is_only_almost_positive(self, hot_index):
        '''
        @return:
            whether the hot index is negative, but one of the flanking indices is positive. If so, this window shouldn't
            be considered during the learning process (we treat it neither positive nor negative), assuming that default
            configuration is used.
        '''
        mask = self.get_annotation_mask()
        return mask[hot_index] == '0' and (mask[hot_index - 1] == '1' or mask[hot_index + 1] == '1')

    def get_features(self, hot_index, feature_keys = features.DEFAULT_FEATURE_KEYS):
        return features.get_features(self, hot_index, feature_keys)

    def __repr__(self):
        return 'Window %d of %s' % (self.offset, self.full_record.id)
