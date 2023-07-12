# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import os
import time
import argparse
from subprocess import Popen, PIPE, CalledProcessError
from parse_reference_files import readFasta
import re

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input", dest = "input",\
                    help="Directory containing downloaded protein.faa.gz files from RefSeq")
parser.add_argument("-ouput", "--ouput", dest = "ouput",\
                    default = 'data/reference/Actinopterygii/Actinopterygii_refseq',\
                    help="")
parser.add_argument("-t", "--taxon_list", dest = "taxon_list",\
                    default = 'data/reference/Actinopterygii/ancestor7898_proteomes.tsv',\
                    help="list of taxon ids with ancestor:7898 from UniProt")
parser.add_argument("-makeblastdb", "--makeblastdb", dest = "makeblastdb", \
                    default = 'src/ncbi-blast-2.14.0+/bin/makeblastdb', \
                    help="makeblastdb path")
args = parser.parse_args()
 

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

uniprot_taxlist_fn = os.path.join(base_dir, os.path.normpath(args.taxon_list))
if not os.path.isfile(uniprot_taxlist_fn):
    raise ValueError('DNE: ', uniprot_taxlist_fn)

blast_path = os.path.join(base_dir, os.path.normpath(args.makeblastdb))
if not os.path.isfile(blast_path):
    err_msg = 'blast path DNE'
    raise ValueError(err_msg)

input_dir = os.path.normpath(args.input)
if not os.path.isdir(input_dir):
    raise ValueError('Check input directory path...')

ouput_path = os.path.join(base_dir, os.path.normpath(args.ouput))
fasta_directory = os.path.join(ouput_path, 'fasta')
blast_db_dir = os.path.join(ouput_path, 'blast_db')
for directory in [fasta_directory, blast_db_dir]:
    if not os.path.isdir(directory):
        os.makedirs(directory)


sciname_to_taxon = {}
scientific_names_list = set()
with open(uniprot_taxlist_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            taxon = sp[colidx['Taxon Id']]
            Scientific_name = sp[colidx['Scientific name']]
            scientific_names_list.add(Scientific_name)
            sciname_to_taxon[Scientific_name] = taxon
print('taxons: ', len(scientific_names_list))


def readFasta_taxon(fastafn):
    fasta_proteins = {}
    sequence = []
    accession = ''
    header = ''
    skip_entry = False
    organism = ''
    with open(fastafn,'r') as infile:
        for line in infile:
            if line[0]=='>':
                skip_entry = False
                #>NP_001012512.1 peptide chain release factor 1, mitochondrial [Danio rerio]
                if sequence:
                    if organism not in fasta_proteins:
                        fasta_proteins[organism] = {}
                    fasta_proteins[organism][accession] = {'seq':sequence,'header':header}
                sequence = []
                accession = line.strip()[1:].split()[0]
                header = line.strip()
                organism = re.search(r"\[.*?\]", header)[0]
                organism = organism.replace('[','').replace(']','')
                if organism not in scientific_names_list:
                    skip_entry = True
                    organism = ''
                    accession = ''
                    header = ''
                    continue
            else:
                if skip_entry==False:
                    sequence.append(line.strip())
        if sequence:
            if organism not in fasta_proteins:
                fasta_proteins[organism] = {}
            fasta_proteins[organism][accession] = {'seq':sequence,'header':header}
    return fasta_proteins

for fn in os.listdir(input_dir):
    print(fn)
    fasta_entries = readFasta_taxon(os.path.join(input_dir, fn))
    for organism in fasta_entries:
        taxon = sciname_to_taxon[organism]
        outfn = os.path.join(fasta_directory, taxon+'_refseq.fasta')
        if os.path.exists(outfn):
            outfile = open(outfn,'a')
        if not os.path.exists(outfn):
            outfile = open(outfn,'w')
        for accession in fasta_entries[organism]:
            header = fasta_entries[organism][accession]['header']
            seq = fasta_entries[organism][accession]['seq']
            outfile.write(header+'\n')
            outfile.write('\n'.join(seq)+'\n')
        outfile.close()

                                        
print('prepare files for blast db creation...')
combined_fn = os.path.join(blast_db_dir, 'Actinopterygii.fsa')
map_fn = os.path.join(blast_db_dir, 'Actinopterygii_map.txt')
with open(combined_fn,'w') as outfile, open(map_fn,'w') as outfile2:
    for fn in os.listdir(fasta_directory):
        taxon = str(fn).replace('_refseq.fasta','')
        fasta_fn = readFasta(os.path.join(fasta_directory, fn))
        for accession in fasta_fn:
            protseq = fasta_fn[accession]['seq']
            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            outfile.write('>'+accession+'\n'+'\n'.join(sequence_split)+'\n')
            outfile2.write(accession+' '+taxon+'\n')

#create blast db
def has_bad_exception(content):
    lines = content.splitlines()
    for i, line in enumerate(lines):
        lower_case_line = line.decode('utf8').lower()
        if "[error]" in lower_case_line:
            print(lower_case_line)
            return True
        elif "exception" in lower_case_line:
            print(lower_case_line)
            return True
    return False

def cmdrun(cmd_listform):
    start_time = time.time()
    process = Popen(cmd_listform, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    # check the output to determine if a silent error occurred
    if process.returncode == 0:
        if (has_bad_exception(stdout) or has_bad_exception(stderr)):
            raise CalledProcessError(process.returncode, cmd_listform, output=stdout)
        else:
            msg = "complete: "+str(time.time()-start_time)+"s"
            return msg
    else:
        raise CalledProcessError(process.returncode, cmd_listform, output=stdout)


cmd = [blast_path,
       '-in',
       combined_fn,
       '-parse_seqids',
       '-blastdb_version',
       '5',
       '-taxid_map',
       map_fn,
       '-title',
       'Actinopterygii_refseq',
       '-dbtype',
       'prot']

result = cmdrun(cmd)
print(result)
