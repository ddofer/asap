from StringIO import StringIO
import logging

from Bio import SeqIO

from django.template import RequestContext
from django.shortcuts import render_to_response

from cleavepred import simple_cleavage_predictor
from cleavepred.util import split_to_chunks

LOGGER = logging.getLogger('WEB')

def cleavage_prediction(request):
    seqs_data = _get_seqs_data(_get_raw_seqs(request))
    return render_to_response('cleavage-prediction.html', {'seqs_data': seqs_data}, context_instance = RequestContext(request))
    
def _get_raw_seqs(request):
    if 'seqs-file' in request.FILES:
        return request.FILES['seqs-file'].read()
    else:
        return request.POST.get('seqs-text')
    
def _get_seqs_data(raw_seqs):
    
    raw_seqs = _fix_string_newlines(unicode(raw_seqs)) # This patch is required for some reason, because Django somehow corrupts uploaded files
    LOGGER.info('Received a %d bytes long FASTA' % len(raw_seqs))
    seqs_fasta = StringIO(_fix_fasta_if_needed(raw_seqs))
    records = list(SeqIO.parse(seqs_fasta, 'fasta'))
    LOGGER.info('About to process %d records' % len(records))
    
    for record in records:
        
        id = record.id
        seq = str(record.seq)
        LOGGER.info('Processing record %s: %s' % (id, seq))
        
        cleavage_mask, cleavage_products = simple_cleavage_predictor.predict(seq, proba = True)
        labeled_aa_chunks = split_to_chunks([_LabeledAminoAcid(aa, label) for aa, label in zip(seq, cleavage_mask)], _RESIDUES_TO_DISPAY_PER_ROW)
        yield id, labeled_aa_chunks, cleavage_products
        
def _fix_fasta_if_needed(raw_fasta):
    if _has_fasta_headers(raw_fasta):
        return raw_fasta
    else:
        return _DEFAULT_FASTA_HEADER + '\n' + raw_fasta
    
def _has_fasta_headers(raw_fasta):
    
    for line in raw_fasta.splitlines():
        if line.startswith('>'):
            return True
            
    return False
    
def _fix_string_newlines(string):
    
    fixed_string = ''
    
    for i in xrange(len(string)):
        if string[i] == '\r' and i < len(string) - 1 and string[i + 1] != '\n':
            fixed_string += '\r\n'
        else:
            fixed_string += string[i]
            
    return fixed_string
    
class _LabeledAminoAcid(object):
    
    def __init__(self, aa, cleavage_probability):
        self.aa = aa
        self.cleavage_probability = cleavage_probability
        
    def is_cleavage(self):
        return self.cleavage_probability >= 0.5
        
    def background_color(self):
        if self.is_cleavage():
            return '#ff0000'
        else:
            return '#ffffff'
        
    def probability_color(self):
        if self.cleavage_probability <= 0.0:
            return '#aaaaaa'
        if self.cleavage_probability < 0.1:
            return '#666666'
        if self.cleavage_probability < 0.3:
            return '#aaaa00'
        if self.cleavage_probability < 0.5:
            return '#ffaa00'
        else:
            return '#ff0000'
            
    def __repr__(self):
        return '<%s>' % self.aa

_DEFAULT_FASTA_HEADER = '>input_seq'
_RESIDUES_TO_DISPAY_PER_ROW = 20
