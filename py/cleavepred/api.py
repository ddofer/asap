import pickle

from . import project_paths

class CleavagePredictor(object):

    '''
    A predictor trained to predict the cleavage of peptides.
    There should be only two instances of this class:
    1. simple_cleavage_predictor - Uses only the basic features derived from the amino-acid sequence of peptides.
    2. advanced_cleavage_predictor - Used also features derived from external tools (ss, acc, disorder and pssm).
    '''

    def __init__(self, advanced):
        self.advanced = advanced
        self._peptide_predictor = None

    def predict(self, seq, extra_tracks_data = {}, proba = False):
        '''
        Predicts cleavage for a given peptide.
        @param seq (string):
            The amino-acid sequence of the peptide to predict the annotations for, given in a 20 amino-acid alphabet.
        @param extra_tracks_data (dict, empty by default):
            A dictionary for providing extra tracks of the given peptide. If using the simple predictor (i.e. advanced = False), it can
            be left empty. If using the advanced predictor (i.e. advanced = True), must receive all tracks (i.e. ss, acc, disorder
            and pssm). The given dictionary should map from track names to their sequence.
        @param proba (default False):
            Whether to return mask of predicted probabilities (floats from between 0 to 1) or binary labels (0s or 1s).
        @return:
            A tuple composed of:
            1. cleavage_mask - If proba = False, it will be a binary string (0's and 1's) representing whether each residue is a cleavage
            site (1) or not (0). If proba = True, it will be a list of floats (between 0 to 1) representing the probability of each residue
            to be a cleavage site. Either way, the length of the returned string/list will correspond to the length of the provided peptide
            sequence.
            2. cleavage_products - A list of strings, each representing the amino-acid sequence of a predicted cleavage product.
        '''
        cleavage_mask = self.get_peptide_predictor().predict_annotations(seq, extra_tracks_data = extra_tracks_data, proba = proba)
        cleavage_products = _get_cleavage_products(seq, cleavage_mask)
        return cleavage_mask, cleavage_products

    def get_peptide_predictor(self):

        '''
        @return:
            The PeptidePredictor object associated with this cleavage predictor.
        '''

        if self._peptide_predictor is None:
            self._peptide_predictor = self._load_peptide_predictor()

        return self._peptide_predictor

    def _load_peptide_predictor(self):

        predictor_dump_file = open(project_paths.get_peptide_predictor_dump_file_path(self.advanced), 'rb')

        try:
            return pickle.load(predictor_dump_file)
        finally:
            predictor_dump_file.close()

simple_cleavage_predictor = CleavagePredictor(False)
advanced_cleavage_predictor = CleavagePredictor(True)

def _get_cleavage_products(seq, cleavage_mask):

    products = []
    current_product = ''

    for i in range(len(seq)):
    
        current_product += seq[i]
    
        # When we have continuous positive cleavage sites, we consider only the most C-terminus one.
        if _is_cleavage(cleavage_mask[i]) and (i >= len(seq) - 1 or not _is_cleavage(cleavage_mask[i + 1])):
            _add_if_not_empty(products, current_product)
            current_product = ''

    _add_if_not_empty(products, current_product)
    return products
    
def _is_cleavage(label):
    if isinstance(label, str):
        return label == '1'
    elif isinstance(label, int) or isinstance(label, float):
        return int(round(label)) == 1
    else:
        raise Exception('Unknown label type: ' + str(type(label)))

def _add_if_not_empty(array, string):
    if len(string) > 0:
        array += [string]
