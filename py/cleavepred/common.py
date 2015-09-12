from asap import WindowExtractionParams
from asap.config import POSITIVE_AMINO_ACIDS

AVAILABLE_TRACKS = [
    'ss',
    'acc',
    'disorder',
    'pssm',
]

def windows_filter(window):
    # We consider only interesting windows with a positively charged amino-acid (i.e. K/R) in the hot index (only then it can be a
    # cleavage candidate).
    return window.get_aa_seq()[window_extraction_params.window_hot_index] in POSITIVE_AMINO_ACIDS
    
window_extraction_params = WindowExtractionParams(window_prefix = 11, window_suffix = 8, neighbourhood_prefix = 5, \
        neighbourhood_suffix = 5, windows_filter = windows_filter)
