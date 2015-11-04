from asap import FEATURE_KEY_OPTIONS, WindowExtractionParams
from asap.config import POSITIVE_AMINO_ACIDS

AVAILABLE_TRACKS = [
    'ss',
    'acc',
    'disorder',
    'pssm',
]

# Here we prefer using 'aa_reduced' over 'aa'. We give up on some other features.
USED_FEATURES = set(FEATURE_KEY_OPTIONS).difference(['aa', 'accum_charge_left', 'accum_charge_right', 'accum_pos_charge_left', 'accum_pos_charge_right'])

def windows_filter(window):
    '''
    We consider only windows with a positively charged amino-acid (i.e. K/R) in the hot index (only then it can be a
    cleavage candidate).
    '''
    return window.get_aa_seq()[window_extraction_params.window_hot_index] in POSITIVE_AMINO_ACIDS

window_extraction_params = WindowExtractionParams(window_prefix = 11, window_suffix = 8, neighbourhood_prefix = 5, \
        neighbourhood_suffix = 5, windows_filter = windows_filter, feature_keys = USED_FEATURES)
