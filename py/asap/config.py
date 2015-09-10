import logging

AMINO_ACIDS = '_ACDEFGHIKLMNPQRSTVWY'
REDUCED_AMINO_ACIDS = 'ACDGFIHKMNQPSRW'
POSITIVE_AMINO_ACIDS = 'KR'
NEGATIVE_AMINO_ACIDS = 'DE'

SS_OPTIONS = '_HCE' # 3 state 2D structure
ACC_OPTIONS = '_-e' # Binary
DISORDER_OPTIONS = '.-^'
PSSM_AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV' # The order is super-important here!

# Init logger
LOG_FORMAT = '%(asctime)s [%(name)s:%(levelname)s] %(message)s'
logging.basicConfig(format = LOG_FORMAT, level = 'INFO')
