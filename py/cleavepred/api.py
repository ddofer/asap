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

    def predict(self, seq, extra_tracks_data = {}):
        '''
        Predicts cleavage for a given peptide.
        @param seq (string):
         The amino-acid sequence of the peptide to predict the annotations for, given in a 20 amino-acid alphabet.
        @param extra_tracks_data (dict, empty by default):
            A dictionary for providing extra tracks of the given peptide. If using the simple
            predictor (i.e. advanced = False), it can be left empty. If using the advanced predictor (i.e. advanced = True), must receive all
            tracks (i.e. ss, acc, disorder and pssm). The given dictionary should map from track names to their sequence.
        @return:
            A tuple composed of:
            1. cleavage_mask - A binary string (0's and 1's) representing whether each residue is a cleavage site (1) or not (0). The length
            of the returned string corresponds to the length of the provided peptide sequence.
            2. cleavage_products - A list of strings, each representing the amino-acid sequence of a predicted cleavage product.
        '''
        cleavage_mask = self.get_peptide_predictor().predict_annotations(seq, extra_tracks_data = extra_tracks_data)
        cleavage_products = _get_cleavage_products(seq, cleavage_mask)
        return cleavage_mask, cleavage_products

    def get_peptide_predictor(self):

        '''
        @return: The PeptidePredictor object associated with this cleavage predictor.
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

    for aa, label in zip(seq, cleavage_mask):
        if label == '1':
            _add_if_not_empty(products, current_product)
            current_product = ''
        else:
            current_product += aa

    _add_if_not_empty(products, current_product)
    return products

def _add_if_not_empty(array, string):
    if len(string) > 0:
        array += [string]
