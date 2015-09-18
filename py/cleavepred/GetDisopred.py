'''
Script to run disopred3 on a multifasta file, then parse and collate the output.
Steps:
1. Extract each fasta from a copy of original multifasta file. (Stored in a seperate directory)
2. (Opt?) Remove original multifasta file.
3. Run disopred on each fasta.
4. i. gather the data/output for all the fastas. (output file name is the fasta's name).
4. ii. Clean the format (for each file), and save to file: "output_feat.diso", (in a format like lf/ss/acc)
4.iii. Save this output to the external features folder.
5. Import from the standard pipeline (config) / not here.
'''

import os
import subprocess
from subprocess import call
import sys
import csv
import glob
import pandas as pd
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


#Location of the directory containing "run_disopred.pl"
DISOPRED_LOCATION = r'/cs/stud/danofer/Desktop/danofer/Software/DISOPRED'
DISO_PROG = 'run_disopred.pl'

# FASTA_LOCATION = '/cs/prt3/danofer/CleavePred/Dataset/Neuropred/V3/'
# FASTA_TARGET = 'NeuroPred_nNames_70ID.fasta'

#NEW:
# FASTA_LOCATION = '/cs/prt3/danofer/CleavePred/data/uniprot/UniProtTestSeqs/D3/'
FASTA_LOCATION = '/a/fr-05/vol/protein/danofer/imac/Desktop/DFTP/'

# FASTA_TARGET = 'D3_TEST_50_FILT.fasta'
FASTA_TARGET = 'D3_TEST_50_FILT_mod.fasta'

SPLIT_FASTAS_DIR = 'splitfasta1'
split_fastas_dir = os.path.join(FASTA_LOCATION,SPLIT_FASTAS_DIR)


file_in = os.path.join(FASTA_LOCATION,FASTA_TARGET)

ALT_DIR_OUTPUT = os.path.join(FASTA_LOCATION,'altFiltered')



def parse_DISOPRED():
    '''
    Parse all pbdat files in a dir, (Output of DISOPRED)
    save their output in a format like hssp/ss/acc
    '''
    # os.chdir(os.path.dirname(split_fastas_dir))
    os.chdir(split_fastas_dir)

    files = glob.glob('*.pbdat')
    print('amount of .pbdat files:',len(files))
    # output_file = open('/cs/prt3/danofer/CleavePred/Dataset/Neuropred/V2_ExternalFeat_NP/output_feat.DISO', 'w')
    # output_file = open(FASTA_LOCATION+'output_feat.DISO', 'w')
    output_file = open(os.path.join(FASTA_LOCATION,'output_feat.DISO'), 'w')

    # print('Joined results will be saved to ',(FASTA_LOCATION+'output_feat.DISO'))
    print('Joined results will be saved to ',os.path.join(FASTA_LOCATION+'output_feat.DISO'))

    for f in files:
        accession = str('>'+os.path.splitext(os.path.basename(f))[0])
        f=open(f)
        lines = f.readlines()
        seq = []
        diso_state=[]
        for line in lines:
            parts = line.strip(' \n \t').split(' ')
            if parts[0] != '#':
                seq += parts[1]
                diso_state += parts[2]
        seq = ''.join(seq)
        diso_state = ''.join(diso_state)
        output_file.write(accession+'\n')
        output_file.write((seq)+'\n')
        output_file.write((diso_state)+'\n')
        # print(accession)
    print("Saved to output_feat.DISO")
    output_file.close()


def split_fasta(filter = False):
    '''
    https://py4bio.wordpress.com/2009/07/22/split_fasta_file/
    This script takes a fasta file and split it in one file per fasta entry.
    It saves the outputs fastas in a new directory
    '''
    os.chdir(os.path.dirname(FASTA_LOCATION))
    print("Current working Directory:",os.getcwd())
    filter_fastas = []
    file_in = os.path.join(FASTA_LOCATION,FASTA_TARGET)
    split_output_dir = SPLIT_FASTAS_DIR

    if filter == True:
        filter_fastas = filter_fasta_queries(fastas_dir=split_fastas_dir)
        split_output_dir = ALT_DIR_OUTPUT

    if not os.path.exists(split_output_dir):
        os.makedirs(split_output_dir)

    os.chdir(split_output_dir)

    i = 0
    read_counts = 0
    for record in SeqIO.parse(open(file_in), "fasta"):
        read_counts += 1
        if not (record.id in filter_fastas):
            # f_out = os.path.join(split_output_dir,record.id+'.fasta')
            f_out = (record.id+'.fasta')
            # f_out =(split_output_dir+record.id+'.fasta')
            print('save to:',f_out)
            # SeqIO.write([record],open(f_out,'w'),"fasta")
            with open(f_out, "w") as handle:
                SeqIO.write([record], handle, "fasta")
            i += 1

    print(read_counts," = Fastas in the original multifasta-file")
    print(i," = # Splitted Fasta files made")


def call_DISOPRED(split_fastas_list):
    os.chdir(os.path.dirname(DISOPRED_LOCATION))
    print(os.getcwd())
    print("In DisoPred folder")
    for i, fasta in enumerate(split_fastas_list):
        print(i)
        print(fasta)
        subprocess.call([DISOPRED_LOCATION+'/'+DISO_PROG,fasta])
        print()

def filter_fasta_queries(fastas_dir=split_fastas_dir):
    '''
    If Disopred job was interrupted in midway -
    This lets us continue (in a new dir) for
    only those sequences that do not have
    Disopred predictions = *.pbdat
    '''
    os.chdir(fastas_dir)
    files = glob.glob('*.pbdat')
    print('# .pbdat files = fastas that were processed succesfully, previously:',len(files))
    ids = [os.path.splitext(os.path.basename(f))[0] for f in files]
    print('len(ids)',len(ids))
    return ids

def find_missing():
    '''
    Disopred seems to "miss" some fastas.
    This helps us find them, assuming the
    disopred output is in the same dirr as
    the (split / "filtered") our fasta candidates
    '''
    # os.chdir('/cs/prt3/danofer/CleavePred/Dataset/Uniprot/altFiltered')
    os.chdir(split_fastas_dir)

    fastas = [f for f in os.listdir('.') if f.endswith('.fasta')]
    ids = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    print('# Fastas present:',str(len(ids)))

    preds = fastas = [f for f in os.listdir('.') if f.endswith('.pbdat')]
    dis_ids = [os.path.splitext(os.path.basename(f))[0] for f in preds]
    print('# Predictions present:',str(len(dis_ids)))

    missing = [a for a in ids if a not in dis_ids]
    print('Missing IDs:')
    print(missing)
    return missing


if __name__ == '__main__':

    # fe = filter_fasta_queries()
    # print('len filter_fasta_queries()',len(fe))
    # split_fasta(filter = True)

    SPLIT_F = False
    CALL_DISO = False
    USE_FILTERED = False

    RESUME_DISO_PARTIAL = True

    if SPLIT_F == True:
        split_fasta()

    # split_fastas_list = glob.glob(split_fastas_dir+'/*.fasta')
    split_fastas_list = glob.glob(os.path.join(split_fastas_dir+'/*.fasta'))
    print('\n split_fastas_list in dir: ',len(split_fastas_list))

    if USE_FILTERED == True:
        # filt_split_fastas_list = glob.glob(ALT_DIR_OUTPUT+'/*.fasta') #ORIGINAL
        filt_split_fastas_list = glob.glob(split_fastas_dir+'/*.fasta')  #CHANGED #D
        print('\n filt_split_fastas_list ',len(filt_split_fastas_list))
        print(filt_split_fastas_list[0])
        print(filt_split_fastas_list[1])
        call_DISOPRED(filt_split_fastas_list)


    if CALL_DISO == True:
        call_DISOPRED(split_fastas_list)
        print("\n DISOPRED DONE! \n")

    parse_DISOPRED()

    missing_ID = find_missing()
    if RESUME_DISO_PARTIAL:
        # fe = filter_fasta_queries()
        print("\n Calling disopred on missing IDs")
        print("Missing: \n",missing_ID)
        call_DISOPRED(split_fastas_list)
        print("\n DISOPRED DONE! \n")


