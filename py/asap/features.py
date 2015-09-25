from . import config
from .features_deps.Disorder import getDisordered
from .features_deps.AAScales import *

import math
from collections import Counter, defaultdict

FEATURE_KEY_OPTIONS = [
    'position', # Features related to the position of the window within the whole protein
    'basic_properties', # PI, MW, GRAVY, aromaticity, aliphaticness, Net Charge.
    'scales', # Extracts quantitative features (e.g. Average at each position) using multiple amino acid propensity scales, for a given window size
    'cysteine_motifs', # Check for putative Cysteine spacer motifs. ('C[^C]{0,3}C' , 'C[^C]{0,3}C[^C]{15,40}C')
    'fid_disorder', # Get for putative disordered segment(s) according to FoldIndex method
    'kr_motifs', # Get putative canonical/"Known motif" cleavage sites, according to a regular expression.
    'charge', # The charge (+1, -1 or 0) of each residue in the window
    # 'accum_charge_left', # Accumulating charge from the left of the window
    # 'accum_charge_right', # Accumulating charge from the right of the window
    # 'accum_pos_charge_left', # Accumulating charge from the left of the window where only positive amino-acids (K and R) are considered
    # 'accum_pos_charge_right', # Accumulating charge from the right of the window where only positive amino-acids (K and R) are considered
    # 'aa', # the actual amino-acid at each position given by one-hot encoding
    'aa_reduced', # Same as 'aa', after using a reduced alphabet of 15 amino-acids
    'aa_context_count', # Counting the number of occurrences of each amino-acid in the part of the protein to the left/right of the window
    'aa_counts', # Counting the number of occurrences of each amino-acid in the window
    'ss', # The given secondary-structure prediction at each position in the window given by one-hot encoding
    'ss_context_count', # Like 'aa_context_count', only for ss instead of aa
    'ss_segment', # Get the total amount of each type of "letter/state" in each subsegment of the sequence. Default is dividing the sequence into 3 segments.
    'acc', # The given accessibility prediction at each position in the window given by one-hot encoding
    'acc_context_count', # Like 'aa_context_count', only for acc instead of aa
    'acc_segment', # Get the total amount of each type of "letter/state" in each subsegment of the sequence. Default is dividing the sequence into 3 segments.
    'disorder', # The given disorder prediction at each position in the window given by one-hot encoding
    'disorder_context_count', # Like 'aa_context_count', only for disorder instead of aa
    'disorder_segment', # Get the total amount of each type of "letter/state" in each subsegment of the sequence. Default is dividing the sequence into 3 segments.
    'pssm', # AA Frequency from the PSSM at each position in the window
    'pssm_entropy', # Entropy of the PSSM at each position
]

DEFAULT_FEATURE_KEYS = [
    'position',
    'basic_properties',
    'scales',
    'fid_disorder',
    'charge',
    'accum_charge_left',
    'accum_charge_right',
    'accum_pos_charge_left',
    'accum_pos_charge_right',
    'aa',
    'aa_context_count',
    'aa_counts',
    'ss',
    'ss_context_count',
    'ss_segment',
    'acc',
    'acc_context_count',
    'acc_segment',
    'disorder',
    'disorder_context_count',
    'disorder_segment',
    'pssm',
    'pssm_entropy',
]

def get_features(window, hot_index, feature_keys = DEFAULT_FEATURE_KEYS):

    if feature_keys is None:
        feature_keys = FEATURE_KEY_OPTIONS
    elif not set(feature_keys).issubset(FEATURE_KEY_OPTIONS):
        raise Exception('Unknown feature keys: ' + ', '.join(map(str, set(feature_keys).difference(FEATURE_KEY_OPTIONS))))

    features = {}

    if 'position' in feature_keys:
        features.update(get_position_features(window))
    if 'basic_properties' in feature_keys:
        features.update(get_basic_properties_features(window.get_aa_seq().replace('_', '').replace('X', '')))
    if 'scales' in feature_keys:
        features.update(get_scales_features(window.get_aa_seq()))
    if 'cysteine_motifs' in feature_keys:
        features.update(getCysteineMotifs(window.get_aa_seq()))
    if 'fid_disorder' in feature_keys:
        features.update(getFIDisorder(window.get_aa_seq()))
    if 'kr_motifs' in feature_keys:
        features.update(GetKRMotifCounts(window.get_aa_seq(), hot_index))

    if 'charge' in feature_keys:
        features.update(get_charge_features(window.get_aa_seq()))
    if 'accum_charge_left' in feature_keys:
        features.update(get_accumulating_charge_features(window.get_aa_seq(), 'left'))
    if 'accum_charge_right' in feature_keys:
        features.update(get_accumulating_charge_features(window.get_aa_seq()[::-1], 'right'))
    if 'accum_pos_charge_left' in feature_keys:
        features.update(get_accumulating_positive_charge_features(window.get_aa_seq(), 'left'))
    if 'accum_pos_charge_right' in feature_keys:
        features.update(get_accumulating_positive_charge_features(window.get_aa_seq()[::-1], 'right'))

    if 'aa' in feature_keys:
        features.update(get_aa_features(window.get_aa_seq()))
    if 'aa_reduced' in feature_keys:
        features.update(get_reduced_aa_features(window.get_aa_seq()))
    if 'aa_context_count' in feature_keys:
        features.update(get_context_count_features('aa', config.AMINO_ACIDS, window, 'aa'))
    if 'aa_counts' in feature_keys:
        features.update(get_aa_counts_features(window.get_aa_seq()))

    if 'ss' in window.get_available_tracks():
        if 'ss' in feature_keys:
            features.update(get_ss_features(window.get_track_seq('ss')))
        if 'ss_context_count' in feature_keys:
            features.update(get_context_count_features('ss', config.SS_OPTIONS, window, 'ss'))
        if 'ss_segment' in feature_keys:
            features.update(get_segment_features(window.get_track_seq('ss'), hot_index, 'ss'))

    if 'acc' in window.get_available_tracks():
        if 'acc' in feature_keys:
            features.update(get_acc_features(window.get_track_seq('acc')))
        if 'acc_context_count' in feature_keys:
            features.update(get_context_count_features('acc', config.ACC_OPTIONS, window, 'acc'))
        if 'acc_segment' in feature_keys:
            features.update(get_segment_features(window.get_track_seq('acc'), hot_index, 'acc'))

    if 'disorder' in window.get_available_tracks():
        if 'disorder' in feature_keys:
            features.update(get_disorder_features(window.get_track_seq('disorder')))
        if 'disorder_context_count' in feature_keys:
            features.update(get_context_count_features('disorder', config.DISORDER_OPTIONS, window, 'disorder'))
        if 'disorder_segment' in feature_keys:
            features.update(get_segment_features(window.get_track_seq('disorder'), hot_index, 'disorder'))

    if 'pssm' in window.get_available_tracks():
        if 'pssm' in feature_keys:
            features.update(get_pssm_features(window.get_track_seq('pssm')))
        if 'pssm_entropy' in feature_keys:
            features.update(Entropy_pssm(window.get_track_seq('pssm'), hot_index))

    return features

def get_position_features(window):

    total_size = window.full_record.length
    left_size = window.offset
    right_size = total_size - left_size - window.length

    return {
        'total_size': total_size,
        'left_size': left_size,
        'right_size': right_size,
        'left_size_percentage': float(left_size) / float(total_size),
        'right_size_percentage': float(right_size) / float(total_size),
    }

def get_context_count_features(feature_name, options, window, track_name):
    features = {}
    features.update(get_count_features(feature_name + '_left_context', options, window.get_left_context_track_seq(track_name)))
    features.update(get_count_features(feature_name + '_right_context', options, window.get_right_context_track_seq(track_name)))
    return features

def get_count_features(feature_name, options, seq):

    features = {}

    for option in options:
        features[feature_name + ':' + option] = seq.count(option)

    return features

def get_reduced_aa_seq(seq):
    '''
    Get a reduced alphabet version of a sequence.
    (Could be used also for SS8 to SS5 !).
    This uses the OFER13KR alphabet by default, due to differences in
    python's Translate between Py3 and 2.7! *(Otherwise, AAlphabets.py could be used)
    '''

    from string import maketrans

    BaseAlph = 'ACEDGFIHKMLNQPSRTWVY'
    OutAlph =  'ACDDGFIHKMINQPSRSWIF'
    tran = maketrans(BaseAlph,OutAlph)

    return seq.translate(tran)

def get_aa_features(seq):
    return get_multi_option_features('aa', config.AMINO_ACIDS, seq)

def get_reduced_aa_features(seq):
    return get_multi_option_features('r_aa', config.REDUCED_AMINO_ACIDS, get_reduced_aa_seq(seq))

def get_ss_features(ss):
    return get_multi_option_features('ss', config.SS_OPTIONS, ss)

def get_acc_features(acc):
    return get_multi_option_features('acc', config.ACC_OPTIONS, acc)

def get_disorder_features(disorder):
    return get_multi_option_features('disorder', config.DISORDER_OPTIONS, disorder)

def get_aa_counts_features(seq):
    features = {}

    for aa in config.AMINO_ACIDS:
        feature_name = 'aa_count:' + aa
        features[feature_name] = seq.count(aa)

    return features


def Dict_Keys_prefix(multilevelDict,PrefixStr):
    '''
    given a dict, Returns a new dict with prefix_str added before each key.
    i.e:

    >>> aa_dict = {'fruit':'orange','K':0.23}
    >>>  print(cds.transform(aa_dict,"entropy"))
    '''
    return {PrefixStr+":"+str(key): (Dict_Keys_prefix(value,PrefixStr) if isinstance(value, dict) else value) for key, value in multilevelDict.items()}

# options_dict = elegant trick :) Dan
def get_segment_features(seq, hot_index, seqtype):

    '''
    Get the total amount of each type of "letter/state" ,
    in each subsegment of the sequence.
    e.g. if seq is the disordered sequence representation:
    Returns the total amount of disorder, non-disordered & predicted-binding in each seg.

    By Default, 3 segments are defined:
    sequence start - 5 positions prior to cleavage;
    cleavege-5 until the cleavage site. (i.e strongest recognition/cleavage site)
    After the cleavage site - until the end of the seq (i.e the putative peptides region)
    '''

    options_dict = {
        'ss': config.SS_OPTIONS,
        'acc': config.ACC_OPTIONS,
        'disorder': config.DISORDER_OPTIONS
    }

    options = options_dict[str(seqtype)] #Get appropiate options

    seg_features = defaultdict(int)
    segments = [seq[0:hot_index-5],seq[hot_index-5:hot_index+1],seq[hot_index+1:]]

    for i,seg in enumerate(segments,start=1):
        for letter in options:
            feature_name = str(seqtype)+'-seg'+str(i)+':'+letter
            seg_features[feature_name] = seg.count(letter)
    return seg_features

def get_charge_features(seq):

    features = {}

    for i, seq_aa in enumerate(seq):
        feature_name = 'charge@' + str(i)
        features[feature_name] = get_aa_charge(seq_aa)

    return features

def get_accumulating_charge_features(seq, side):

    features = {}
    charge = 0

    for i, seq_aa in enumerate(seq):
        charge += get_aa_charge(seq_aa)
        feature_name = 'acc_charge_' + side + '@' + str(i)
        features[feature_name] = charge

    return features

def get_accumulating_positive_charge_features(seq, side):
    features = {}
    charge = 0

    for i, seq_aa in enumerate(seq):

        if seq_aa in config.POSITIVE_AMINO_ACIDS:
            charge += 1

        feature_name = 'acc_pos_charge_' + side + '@' + str(i)
        features[feature_name] = charge

    return features

def get_pssm_features(pssm):
    '''
    PSSM profile AA Order: {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', '_'}
    '''

    features = {}

    for i, profile in enumerate(pssm):
        for aa, value in profile.items():
            feature_name = 'pssm@' + str(i) + ':' + aa
            features[feature_name] = value

    return features

def get_multi_option_features(feature_name, options, values):

    features = {}

    for i, value in enumerate(values):
        for option in options:
            full_feature_name = feature_name + '@' + str(i) + ':' + option
            features[full_feature_name] = (value == option)

    return features

def get_aa_charge(aa):
    if aa in config.POSITIVE_AMINO_ACIDS:
        return 1
    elif aa in config.NEGATIVE_AMINO_ACIDS:
        return -1
    else:
        return 0

# -------------------------------------------------- #
#   Following includes reused code from ProFET       #
# -------------------------------------------------- #

from collections import Counter, defaultdict
import re

from Bio.SeqUtils import ProtParam as pp

'Used to calc protein net charge'
'TODO: Modify PI, Netcharge calcs, to use mofidied deep copuies of these  when no N or C termini'
pKa     = {'D':3.9, 'E':4.3, 'H':6.1, 'C':8.3, 'Y':10.1, 'K':10.5, 'R':12, 'N-term':8, 'C-term':3.1}
charges = {'D':-1,  'E':-1,  'H':+1,  'C':-1,  'Y':-1,   'K':1,    'R':1,  'N-term':1, 'C-term':-1}


def get_scales_features(seq):
    '''
    TODO: Add Average, max, min feature.
    TODO - reimplement window scaling.
    '''
    SCALE_WINDOW_SIZE = 4

    features = {}

    for scale in PTMScales_Dict:

        aa_dict = PTMScales_Dict[scale]
        avg = PTMScales_Avg[scale]
        values = [aa_dict.get(aa, avg) for aa in seq]

        for i in range(len(seq) - SCALE_WINDOW_SIZE + 1):
            feature_name = 'scale_%s_window_%d' % (scale, i)
            feature_value = float(np.average(values[i:(i + SCALE_WINDOW_SIZE)]))
            features[feature_name] = feature_value

    return features

def GetFlex(Bio_ProtParam):
    '''
    Get parameters for B-values / flexibility. Built for a 9 window window;
    alters values from a BioPy ProtParam module.
    '''
    #PP = self.Bio_ProtParam
    flex = Bio_ProtParam.flexibility()
    Flex_mean =np.mean(flex)
    Flex_max=np.amax(flex)
    Flex_min=np.amin(flex)
    res = {'Flex_mean':Flex_mean,'Flex_max':Flex_max,'Flex_min':Flex_min}
    return res


def GetKRMotifLocations(seq):
    '''
    Get potential matches + locations to the Known motif cleave RegEx,
     along the sequence itself.

     known motif model:
     Xx-Xx-K-K# , Xx-Xx-K-R# , Xx-Xx-R-R# ,  R-Xx-Xx-K# , R-Xx-Xx-R#
    lys: K
    aRg: R

     Useful to ID/filter consecutive "overlapping" sites
      (i.e K*R {false},  KR* {true} , where * is cleavage)

        KM_RE_1 = 'R.{2}[RK]'2
        KM_RE_2 = '[^R].K[RK]'
        KM_RE_3 = '[^R].R{2}'

    Getting overlapping matches:
    ~ finditer() won't find overlapping strings ! (Updated RegEx package can)
    http://stackoverflow.com/questions/5616822/python-regex-find-all-overlapping-matches?rq=1
    https://mail.python.org/pipermail/tutor/2005-September/041126.html

    http://code.activestate.com/recipes/499314-find-all-indices-of-a-substring-in-a-given-string/
    '''

    #http://stackoverflow.com/questions/5616822/python-regex-find-all-overlapping-matches?rq=1
    # (?=...) is a lookahead assertion. (Doesn't consume string)
    # matches = re.finditer(r'(?=(\d{10}))',s)
    matches = re.finditer(r'(?=(\d{10}))',s)
    results = [int(match.group(1)) for match in matches]

    Motif_1 = [m.start(0) for m in re.finditer(r'(?=(K[RK]))',seq)]
    Motif_2 = [m.start(0) for m in re.finditer(r'(?=(R.{2}[RK]))',seq)]
    Motif_3 = [m.start(0) for m in re.finditer(r'(?=(R{2}))',seq)]


    Motif_count1 = len(re.findall(KM_RE_1, km_flank)) #Known motif model. Arg-Xxx-Xxx-[Arg|Lys]
    Motif_count2 = len(re.findall('K[RK]', seq)) #Known motif model. Not R at P1, to avoid overlap with other motif!
    Motif_count3 = len(re.findall(KM_RE_3, km_flank)) #Known motif model - Xxx-Xxx-Arg-Arg
    return ({"Motif_count1":Motif_count1,"Motif_count2":Motif_count2,"Motif_count3":Motif_count3})

def GetKRMotifCounts(seq, hot_index) :
    '''
    counts # of suspected cleavage sites according to known motif model:
    Xxx-Xxx-Lys-Lys# , Xxx-Xxx-Lys-Arg# , Xxx-Xxx-Arg-Arg# ,  Arg-Xxx-Xxx-Lys# , # # Arg-Xxx-Xxx-Arg#
    lysine: K. arginine: R.   #: Cleavage site.  "Not Proline" = [^P].
    TODO: CHECK that positions and motif are right! (e.g [hot_index-3] ; ^P ? ...)
    '''
    KnownMotifFound = 0

    #Look at the  4 positions prior to and including cleavage site (i.e length 4):
    km_flank = str(seq[(hot_index-3):(hot_index+1)])

    before_cleavewin_seq =  str(seq[0:(hot_index-3)])
    #Look at a wide potential window - non canonical
    potential_flank = str(seq[(hot_index-5):(hot_index+1)])

    "Check for presence of canonical known motif right before potential cleavage site"

    #Cleavage sites Based on to Known motif model. (+Strong likelihood of NO Proline adjacent)
    #Pay attention when building RegEx to avoid "double counting"!
    #Replaced "seq" with "km_flank" - to count only at the potential cleavage site "legal" location.
    'Replace strings with KM_RE_1/2/3 - TODO - Dan'
    Motif_count1 = len(re.findall('R.{2}[RK]', km_flank)) #Known motif model. Arg-Xxx-Xxx-[Arg|Lys]
    Motif_count2 = len(re.findall('[^R].K[RK]', km_flank)) #Known motif model. Not R at P1, to avoid overlap with other motif!
    Motif_count3 = len(re.findall('[^R].R{2}', km_flank)) #Known motif model - Xxx-Xxx-Arg-Arg

    KR_motif_counts = Motif_count1+Motif_count2+Motif_count3

    if KR_motif_counts>0:
        KnownMotifFound = 1

    #Get counts of possible cleavage sites (KM based) in regions prior to desired cleavage site - might help with noise?
    count1 = len(re.findall('R..[RK]', before_cleavewin_seq))
    #Arg-Xxx-Xxx-Arg|Lys
    count2 = (len(re.findall('.[^P][RK][RK]', before_cleavewin_seq))) #The ^P is not quite necessarily legal always.

    Potential_KR_counts = count1+count2

    #Count amount of positive charged AA in 6 positions prior to cleavage site.
    basic_counts = (potential_flank.count('K')+potential_flank.count('R'))
    # #See if there's at least 2 charged AA before the cleavage site. Good for lax filtering.
    # minimalBasicCount = 0
    # if basic_counts>1:
    #     minimalBasicCount = 1

    '''
    'RR, KK, KR or RxxR where x is any amino acid except Lys or Arg) - Immeditely before the cleavage sites'
    lst = ["RR", "KK", "KR"]
    # #Alt:
    "Is This version is NOT correct ? -  looks at too wide a window (should be just the 4 before cleavage, not 5)"
    km_window_1 = str(seq[(hot_index-4):(hot_index+1)])

    if any(s in km_window_1 for s in lst):
        KnownMotifFound = 1
    if (len(re.findall('R[^RK][^RK]R', km_window_2)) >0):
        KnownMotifFound=1
    '''

    return ({"Potential_KR_counts":Potential_KR_counts,
             "KnownMotifFound":KnownMotifFound,
             "K+R_found_FlankingCleavage":basic_counts,
             "minimalBasicCount":basic_counts
             }
             )

def get_basic_properties_features(seq):
    '''
    Get basic physical properties as in BioPython/ExPasy ProtParam
    module.
    Returns: PI, MW, GRAVY, aromaticity,aliphaticness,Net Charge.

    Note: These methods all assume a standard AA Alphabet.
    Warning! Returned PI is INNACCURATE For a parsed (Tail(s) removed) subseq.
    (BioPy-ProtParam.isoelectric_point assumes N,C terminii!)
    '''
    Bio_ProtParam = pp.ProteinAnalysis(seq) #BioPython SequenceAnalysis object from str

    PI = Bio_ProtParam.isoelectric_point()
    MW=Bio_ProtParam.molecular_weight()
    GRAVY=Bio_ProtParam.gravy()
    aromaticity=Bio_ProtParam.aromaticity()
    aliphaticness = GetAliphaticness(seq)
    NetCharges = get_netCharge(seq)

    prot_pp = {'PI':PI,
    'Molecular_Weight':round(MW,4), 'GRAVY':round(GRAVY,4),
    'Aromaticity':round(aromaticity,4)}

    # #Added now - Dan.
    # flex = GetFlex(Bio_ProtParam) #Returns 3 keys/values
    # prot_pp.update(flex) # Problem with mpty window

    prot_pp.update(aliphaticness)
    prot_pp.update(NetCharges)

    return prot_pp

def cysteineMotif(seq,segDivide=1):
    cysPattern = r'C[^C]{0,3}C'
    length = float(len(seq))  ##If used as part of regular package,maybe use length=self.length instead; seq=self.seq
    # window_size = length / 5  # window size 20% of the protein length #Orig
    window_size = int(length /segDivide)
    # scores = [0 for _ in range(5)]
    scores = [0 for _ in range(segDivide)]
    prog = re.compile(cysPattern)
    pos = 0
    for i in range(segDivide): #Changed to segDivide
        scores[i] = len(prog.findall(seq[pos:pos + window_size]))
        pos += window_size
    scores[-1] = len(prog.findall(seq[pos:])) #-1? Why not do scores.append? - D
    res = {}
    key = "cysteine_window_"
    for i, score in enumerate(scores):
        #Changed if score > 2 to: if score > 1
        res[key + str(i)] = 1 if score > 1 else 0   #Why not >= ? 3 means at least 6 cysteins! That's a lot
    return res

def cysteineSpaceMotif(seq):
    cysPattern = r'C[^C]{0,3}C[^C]{15,40}C'
    match = re.search(cysPattern, seq)
    score = 1 if match else 0 # if match isn't none, the score is 1 ("space motif" was found)
    return {"cysteine_space": score}

def getCysteineMotifs(seq): #(Why have called submethods as static (and not calling to "self")? Less consistant; .)
    '''
    Get Cysteine spacer motif, and counts of frequent CxxC motifs over protein sequence
    '''
    res = {}
    res.update(cysteineMotif(seq))
    res.update(cysteineSpaceMotif(seq))
    return res

def getFIDisorder(seq,segments=3):
    '''
    Divide protein sequence into segments, and returns for each seg.
    predicted disorder (Y/N) according to FoldIndex method (Uversky et al).
    Method implemented with help of Tal Arian.
    '''
    return getDisordered(seq, segments=segments)

def calculateProteinCharge(sequence, pH=7.2): # Has_N_Terminal=True,Has_C_Terminal=True):
    '''
    Get Net-Protein charge at given PH. Can alse be Used to predict PI of a protein.
    Be careful if used on a "partial" seq - +-N/C termini..
    If using only part of a sequence (Just the N-tail etc'), then modify call appropaitely.
    http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt
    We could also try using:
    http://pythonhosted.org/pyteomics/_modules/pyteomics/electrochem.html#charge
    '''
    AA_Counts=Counter(sequence)
    AA_Counter = AA_Counts

    # TODO: this can actually be important. fix it
    #Has_N_Terminal = self.HAS_N
    #Has_C_Terminal = self.HAS_N
    Has_N_Terminal = False
    Has_C_Terminal = False
    N_charge = 0
    C_charge = 0

    # @jit  #autojitCauses slow down for some reason
    def calculateAminoAcidCharge(amino_acid, pH):
        ratio = 1 / (1 + 10**(pH - pKa[amino_acid]))

        if charges[amino_acid] == 1:
            return ratio
        else:
            return ratio - 1

    if Has_N_Terminal:
        N_charge = calculateAminoAcidCharge('N-term', pH)
    if Has_C_Terminal:
        C_charge = calculateAminoAcidCharge('C-term', pH)
    protein_charge = N_charge+C_charge

    for amino_acid in pKa.keys():
        # protein_charge += sequence.count(amino_acid) * calculateAminoAcidCharge(amino_acid, pH)
        protein_charge += AA_Counter[amino_acid] * calculateAminoAcidCharge(amino_acid, pH)
    return protein_charge

def get_netCharge(seq,PH_ranges = [5.5,6.8,7.3,8.1]):
        st = "PH_Charge_"
        res = defaultdict(float)
        # res = Counter()
        'add a map(round(x,3),x)'
        for PH in PH_ranges:
            res [st+str(PH)] = calculateProteinCharge(seq,PH)
        return res

def GetAliphaticness(seq) :
    a = 2.9
    b = 3.9
    #seq = self.seq
    length = float(len(seq))

    #AA_Counts = self.AA_Counts
    AA_Counts=Counter(seq)
    AA_Counter = AA_Counts

    alanine_per = (AA_Counts['A'] / length )
    valine_per = (AA_Counts['V'] / length )
    isoleucine_per = (AA_Counts['I'] / length )
    leucine_per = (AA_Counts['L'] / length )
    # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    aliphatic_index = (100 * (alanine_per + a * valine_per + b * (isoleucine_per + leucine_per )))
    return {'Aliphaticness':aliphatic_index}

'TODO: Maximally independent ICA from current scales. +  Kidera factors..'
def Get_ParamScales(Bio_PP,window=5,edge=0.9,PickScales = PTMScales_Dict): ##,SA_window = 9, TMD_window = 19):
    """
    Transform the overall sequence to various AA scales based representations, and returns them.
    (LATER, we extract features according to locations/masks of subsequences)
    """
    '''
    Gets numerical represention of sequence,
    for each amino acid propensity scale (and window size). This is then
    used (later) to Get values for sequences using different amino propensities,
    via def Get_paramScales.

    Default is "built in" scales, but can be expanded easily.
    Default window size from literature is !17-19 for detecting TMDs' using hydrophobicity.

    Returns a list of (string:list) tuples.
    Uses:  Bio.SeqUtils.ProtParam.ProteinAnalysis(protein_scale(self, param_dict, window))
    http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-module.html - builtin scales
    '''
    # Bio_PP = self.Bio_ProtParam
    PP_scales = []
    if PickScales is None: #Default is to use all preloaded scales from AAScales.py
       aaScales = Scales_Dict
       'Scales_Dict is a dict of dicts from AAScales.py'
    else:
#Use user defined list of scales
        aaScales=PickScales
    for scaleName,v in aaScales.items():
        # PP_scales [str(scaleName)] = Bio_PP.protein_scale(v,window,edge)
        s = str(scaleName)
        v = Bio_PP.protein_scale(aaScales[str(scaleName)],window=window,edge=edge)
        t = (s,v)
        PP_scales.append (t)
    return PP_scales

"TODO: Mod this to work with a given WINDOW, and it's parent sequence as a scale!!"
def Get_ParamScales_Features(seq,window=4,PickScales = MinScales_Dict,segs=None):
    '''
    Similar to "Get_ParamScales_Features", but gets a small(er) set of features (per scale),
    for a smaller set (by efault) of AA scales, extracted from multiple segments of the sequence.

    Extract features for given scale/propensity represention of a protein sequence.
    By default, gets less features (min, max, average), with a minSet of scales. (Includes Atchley).
    An alternative,implementation is to divide into  segments - e.g thirds
    '''
   #  Bio_PP = pp.ProteinAnalysis(seq)
   # #Get list of scales to use.
   #  PP_scales=Get_ParamScales(Bio_PP,window,edge,PickScales=MinScales_Dict)
    # res = defaultdict(float)
    res = {}
    window_prefix = ('_'+str(window)+'_')
    'Check that calced seg length is right!!?'
    if segs is None:
        segs = len(seq)-1 #Check right length. Depends on window size!! Dan

    for scale in PP_scales:
            name=scale[0]
            arr = np.asarray(scale[1])
            seg_size = int(len(arr) / segs)  # window size of each segment of the protein
            pos = 0

            res[str('AA-Scale-'+window_prefix+str(name))+':AVerage'] = np.mean(Arr)
            res[str('AA-Scale-'+window_prefix+str(name))+':MAX'] = np.amax(Arr)
            res[str('AA-Scale-'+window_prefix+str(name))+':MIN'] = np.amin(Arr)
            for i in range(segs-1):
                'Get aa scale feature per individual position'
                res[str('AA-Scale-'+str(i)+"-"+window_prefix+str(name))] = float(seq[i])
                # #subArr = arr[pos:pos + seg_size]
                # res[FeatPrefix+'AV'] = np.mean(subArr)
                # # res[FeatPrefix+'MAX'] = np.amax(subArr)
                # # res[FeatPrefix+'MIN'] = np.amin(subArr)
    return res


##################################################################
# Pick a way of getting entropy -
def log2(number):
    """Shorthand for base-2 logarithms"""
    return math.log(number, 2) if number > 0 else math.log(10E-50, 2)

VERT_AA_FREQ = {
    'A':0.074,'R':0.042,'N':0.044,'D':0.059,
    'C': 0.033 ,'Q': 0.037 ,'E': 0.058 ,
'G': 0.074 ,'H': 0.029 ,'I': 0.038 ,'L': 0.076 ,
'K': 0.072 ,'M': 0.018 ,'F': 0.04 ,'P': 0.05 ,
'S': 0.081 ,'T': 0.062 ,'W': 0.013 ,'Y': 0.033 ,
'V': 0.068,'_':1
}

def Entropy_pssm(pssm, hot_index, get_entropy_segs = True,background_freq_dict = VERT_AA_FREQ):
    '''
    PSSM profile AA Order: {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', '_'}
    TODO: Background frequency correction.
    Background frequencies taken from default for Vertebrates:
    http://www.tiem.utk.edu/~gross/bioed/webmodules/aminoacid.html
    It might be useful to use other background frequencies - e.g. from the training set proteins, or default BloSum.
    Additional sources of background frequencies: http://www.bioinformatics.org/pipermail/ssml-general/2005-July/000203.html

    Calculate the Relative entropy (KL-divergence, with relative frequency accounted for)
     at each position, from a PSSM profile.

    Entropy is the expected value of the measure of information content in system.
    http://rosettacode.org/wiki/Entropy#Python
    https://github.com/talmo/MotifSearch/blob/master/genomics.py

    '''

    features = {}
    # background_freq_uniform = 0.047 #http://bioinf.uab.es/transcout/entropy.html
    entropy_seq = [] #store entropy at each position along the sequence

    for i, profile in enumerate(pssm):

        feature_name = 'EntropyPssm@' + str(i)
        #ORIG
        # features[feature_name] = (-sum( (count * log2(count/background_freq_uniform)) for count in profile.values())) #Pre Change to new background frequencies
        pssm_keys = profile.keys()
        features[feature_name] = (-sum((profile[aa_letter] * log2(profile[aa_letter]/background_freq_dict[aa_letter]) for aa_letter in pssm_keys)))

        entropy_seq.append(features[feature_name])

    if get_entropy_segs:
        entropy_seg_1 = math.fsum(entropy_seq[0:hot_index-4])
        entropy_seg_2 = math.fsum(entropy_seq[hot_index-4:hot_index+1])
        entropy_seg_3 = math.fsum(entropy_seq[hot_index+1:])

        entropy_seg_peptide = math.fsum(entropy_seq[hot_index-2:-2])

        features['entropy_seg_1']=entropy_seg_1
        features['entropy_seg_2']=entropy_seg_2
        features['entropy_seg_3']=entropy_seg_3
        features['entropy_seg_peptide']=entropy_seg_peptide
        #The cleavage site may be before or after the peptide. We attempt to find a conserved stretch:
        features['Max_Local_entropy_segment']=max(entropy_seg_1,entropy_seg_3)

        # for aa, value in profile.items():
    # s=self.seq
    # p=self.AA_Counts
    # lns = self.length
    # def count_e():
    #     # return(-sum( count/lns * log(count/lns, 2) for count in p.values()))
    #     return(-sum( count * log(count, 2) for count in p.values()))
    # entropy=count_e()
    # return {'Sequence Entropy: ':count_e()}
    return features


######################

    # def compute_entropy(motif):
    #     '''
    #     http://stackoverflow.com/questions/23480002/shannon-entropy-of-data-in-this-format-dna-motif?rq=1
    #     '''
    #     arr = np.array(motif)
    #     H = (arr[arr > 0] * np.log2(arr[arr > 0])).sum(axis=1)
    #     print 'entropy:', -H.mean(), 'bits'

    # motif = []
    # for line in sys.stdin:
    #     line = line.strip().lower()
    #     if line.startswith('de'):
    #         print line
    #     elif line == 'xx':
    #         if motif:
    #             compute_entropy(motif)
    #         motif = []
    #     else:
    #         motif.append(map(float, line.split()[1:-1]))

    #     def GetEntropy (self,normalizeTotEntropy=False, getLettersEntropy=True):
    #     '''
    #     http://bugra.github.io/work/notes/2014-05-16/entropy-perplexity-image-text/
    #     VS Entropy_Kap_AA , Entropy_seq..

    #     normalizeTotEntropy : If True, then divide total sequence entropy by log2(length).
    #                             (Check if correct!!)
    #     getLettersEntropy : If True, also return  entropy per letter; If false, then
    #     return only the "total entropy".
    #     '''
    #     # word_count_information = []
    #     AA_information = {}

    #     seq = self.seq
    #     length = float(len(seq))
    #     wordset = set(self.alph)
    #     # freq = self.AA_Counts
    #     "Don't count entropy for letters which don't appear:"
    #     freq = {k:v for (k, v) in self.AA_Counts.items() if v != 0}
    #     entropy = 0
    #     # for word in wordset:
    #     for word in freq.keys():
    #         probability = freq[word] / (1.0 * length)  #Could be replaced with copied use of AA_Freq
    #         self_information = np.log2(1.0/probability)
    #         entropy += (probability * self_information)
    #         # word_count_information.append([word, freq[word], self_information])
    #         if getLettersEntropy==True:
    #             AA_information[str(word)+' Entropy']=self_information

    #     if normalizeTotEntropy is True:
    #         from math import log2
    #         l = log2(length)
    #         AA_information['Total Entropy - Normalized By Length']=(float(entropy)/l)
    #     # else:
    #     AA_information['Total Entropy']=entropy # Equivalent to old Entropy_Seq.

    #     return self.alphabet_prefix(AA_information)


KM_RE_1 = 'R.{2}[RK]'
KM_RE_2 = '[^R].K[RK]'
KM_RE_3 = '[^R].R{2}'
