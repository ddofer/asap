'''
Extract a .lf file, containing sequences and cleavage annotation masks, out of a UniProt's XML file.
Arguments:
- output_file_path (file path, optional): The path to write the output .lf file to. If not provided, will update the project's relevant file.
'''

import sys
import re
import xml.etree.ElementTree as et
from StringIO import StringIO

from cleavepred import util
from cleavepred import project_paths

project_paths.dataset_name = 'uniprot'

if len(sys.argv) > 1:
    output_file_path = sys.argv[1]
else:
    output_file_path = project_paths.get_annotated_seqs_file_path()

def get_unique(element, xpath):

    subelements = element.findall(xpath)

    if len(subelements) == 0:
        return None
    if len(subelements) == 1:
        return subelements[0]
    else:
        raise Exception('%d subelements: %s' % (len(subelements), xpath))

def parse_uniprot_xml(raw_xml_path):
    raw = util.read_file(raw_xml_path)
    fixed_raw = re.sub(r'xmlns="[^"]*"', '', raw)
    return et.fromstring(fixed_raw)

def get_proteins_with_cleavage_sites(raw_xml_path):

    root = parse_uniprot_xml(raw_xml_path)

    for entry in root.findall('./entry'):

        accession = entry.findall('./accession')[0].text
        raw_seq = get_unique(entry, './sequence').text
        seq = re.sub(r'\s', '', raw_seq)

        signal_peptide_end = 0
        cleavage_sites = set()
        skip_protein = False

        for feature in entry.findall('./feature'):

            type = feature.get('type').lower()

            if type in ['peptide', 'chain', 'propeptide', 'signal peptide']:

                try:
                    begin = int(get_unique(feature, './location/begin').get('position'))
                except:
                    begin = None

                try:
                    end = int(get_unique(feature, './location/end').get('position'))
                except:
                    end = None

                if type == 'signal peptide':
                    if end is None:
                        print ('%s: no end to signal peptide. We will ignore this protein.' % accession)
                        skip_protein = True
                        break
                    else:
                        signal_peptide_end = max(signal_peptide_end, end)
                else:

                    if begin is not None:
                        cleavage_sites.add(begin - 1)
                        cleavage_sites.add(begin - 2)

                    if end is not None:
                        if type == 'propeptide':
                            cleavage_sites.add(end - 1)
                        else:
                            cleavage_sites.add(end)

        if skip_protein:
            continue
            
        cleavage_sites = set([i for i in cleavage_sites if i >= signal_peptide_end + 3 and i < len(seq) - 3 and seq[i] in 'KR'])
        cleavage_sites_to_remove = set([i - 1 for i in cleavage_sites]) # If 11, we take only the second
        cleavage_sites = cleavage_sites.difference(cleavage_sites_to_remove)

        if cleavage_sites: # we don't want samples with no cleavages at all - it's probably a mistake
            yield accession, seq, cleavage_sites, signal_peptide_end

def cleavage_sites_to_mask(seq_length, cleavage_sites):

    mask = ['0'] * seq_length

    for cleavage_site in cleavage_sites:
        mask[cleavage_site] = '1'

    return ''.join(mask)
    
def remove_xs(seq, mask):
    
    revised_seq = ''
    revised_mask = ''
    
    for aa, label in zip(seq, mask):
        if aa.lower() != 'x':
            revised_seq += aa
            revised_mask += label
            
    return revised_seq, revised_mask

def space_seq(seq, chunk_length = 10):
    return ' '.join(util.split_to_chunks(seq, chunk_length))

def write_fasta_like_record(file, accession, seq, mask):
    file.write('>' + accession + '\n')
    file.write(space_seq(seq) + '\n')
    file.write(space_seq(mask) + '\n')
    file.write('\n')

if __name__ == '__main__':

    output_file = open(output_file_path, 'wb')

    try:
        for accession, seq, cleavage_sites, signal_peptide_end in get_proteins_with_cleavage_sites(project_paths.get_raw_data_xml_file_path()):
            mask_to_write = cleavage_sites_to_mask(len(seq), cleavage_sites)[signal_peptide_end:]
            seq_to_write = seq[signal_peptide_end:]
            seq_to_write, mask_to_write = remove_xs(seq_to_write, mask_to_write)
            write_fasta_like_record(output_file, accession, seq_to_write, mask_to_write)
    finally:
        output_file.close()

    print 'Done.'
