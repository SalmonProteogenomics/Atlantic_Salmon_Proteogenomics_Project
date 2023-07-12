# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import argparse
import os
from parse_reference_files import readFasta

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23',\
                    help="")
args = parser.parse_args()

print(args.dataset)
fdr = 0.01

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)
##reference files
reference_dir = os.path.join(base_dir, 'data', 'reference')
#RefSeq + Ensembl + Augustus Search directory
RefSeq_Ensembl_Augustus_Search_dir = os.path.join(data_dir,\
                                        '5_RefSeq_Ensembl_Augustus_Search')
#msms meta info  
spectra_meta_fn = os.path.join(data_dir, 'msms_meta_info',\
                               'msms_meta.txt')

spectrumfiles = {}
specfile_order = []
datasets = set()
with open(spectra_meta_fn,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            specfile = sp[colidx['#Sample']].replace(' ','-')
            specfile_order.append(specfile)
            ds = sp[colidx['Dataset']]
            heatmap_label = sp[colidx['heatmap_label']]
            spectrumfiles[specfile] = {'heatmap_label':heatmap_label,
                                       'Dataset':ds}
            datasets.add(ds)
datasets = sorted(datasets)

#RefSeq
print('Reading RefSeq fasta files...')
refseq_dir = os.path.join(reference_dir, 'RefSeq')
refseq_prot_fasta = os.path.join(refseq_dir,\
                                 'GCF_000233375.1_ICSASG_v2_protein.faa')
refseq_proteins = readFasta(refseq_prot_fasta)

#Ensembl
print('Reading Ensembl fasta files...')
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
ensembl_prot_fasta = os.path.join(ensembl_dir,\
                                  'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta) 


#GPM contaminants database
print('Reading GPM contaminants fasta file...')
gpm_fn = os.path.join(reference_dir, 'contaminant_db','gpm_cRAP_012915.fasta')
contaminant_proteins = readFasta(gpm_fn)

#RefSeq Ensembl Augustus DB
print('Reading RefSeq+Ensembl+Augustus fasta file...')
refseq_ensembl_augustus_db = os.path.join(RefSeq_Ensembl_Augustus_Search_dir,\
                                          'fasta_db',\
                                          'RefSeq+Ensembl+Augustus.fasta')
refseq_ensembl_augustus_proteins = readFasta(refseq_ensembl_augustus_db)


#MSGF peptidematch CMD output directories
FDR_dir = os.path.join(RefSeq_Ensembl_Augustus_Search_dir,\
                       'workflow_output', 'MSGF_fdrdir')
Proteinlevel_psms = os.path.join(FDR_dir, 'Proteinlevel_'+str(fdr)+'_FDR_PSMs.fdr')
Proteinlevel_proteins_txt = os.path.join(FDR_dir, 'Proteinlevel_'+str(fdr)+'_FDR_proteins.txt')
exactmatch_knowndir = os.path.join(RefSeq_Ensembl_Augustus_Search_dir,\
                                   'workflow_output', 'ExactMatch_output')

matrix_dir = os.path.join(RefSeq_Ensembl_Augustus_Search_dir,\
                          'knowndb_matrix_proteinlevel_'+str(fdr)+'_FDR')
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

matrix_fn = os.path.join(matrix_dir, 'knowndb_protein_matrix.txt')
matrix_protein_fn = os.path.join(matrix_dir, 'knowndb_protein_matrix_proteins.txt')
#Normalize_Count_Output_File : Normalize PSM count - Protein Spectra per 100 aminoacid per 1000 spectra (PSK)
matrix_fn_psk = os.path.join(matrix_dir, 'knowndb_protein_matrix_psk.txt')
matrix_fn_psk_top_1000 = os.path.join(matrix_dir, 'knowndb_protein_matrix_psk_1000.txt')            
    
 
###########################################################################
#MSGF+ identifications
print('Read Proteinlevel_0.01_FDR_proteins.txt...')
Protein_levelFDR_proteins = {}
with open(Proteinlevel_proteins_txt,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            col = {x:xi for xi,x in enumerate(sp)}
        else:
            proteins = [x.split('/')[0] for x in sp[col['Proteins']].split(';')]
            Peptides = sp[col['Peptides']].split(';')
            Peptides = set([x.replace('I','L') for x in Peptides])            
            for protein in proteins:
                if protein not in Protein_levelFDR_proteins:
                    Protein_levelFDR_proteins[protein] = set()
                Protein_levelFDR_proteins[protein].update(Peptides)

#all proteins passing 1% protein-level FDR 
all_proteins = set()
Known_peptide_proteins = {}
for prot in Protein_levelFDR_proteins:
    peptides = Protein_levelFDR_proteins[prot]
    for p in peptides:
        if p not in Known_peptide_proteins:
            Known_peptide_proteins[p] = set()
        Known_peptide_proteins[p].add(prot)
    all_proteins.add(prot)

print('Read Proteinlevel_0.01_FDR_PSMs.fdr...')
protein_level_fdr_precursors = {x:set() for x in datasets}
protein_level_fdr_peptides = {x:set() for x in datasets}
total_psms = 0
Known_peptides = {}
protein_level_psms = {x:set() for x in datasets}
with open(Proteinlevel_psms,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            col = {x:xi for xi,x in enumerate(sp)}
        else:
            SpecFile = os.path.splitext(sp[col['#SpecFile']])[0]
            if SpecFile not in spectrumfiles:
                raise ValueError(SpecFile)
            ScanNum = sp[col['ScanNum']]
            total_psms+=1
            Charge = sp[col['Charge']]
            SpecEValue = sp[col['SpecEValue']]
            peptide = sp[col['Peptide']][2:-2]
            nomod_pep = ''.join([aa for aa in peptide if aa.isalpha()]).replace('I','L')
            Protein = sp[col['Protein']].split(';')
            precursor = sp[col['precursor_pepcharge']].replace('I','L')
                        
            if nomod_pep not in Known_peptides:
                Known_peptides[nomod_pep] = {}
            if SpecFile not in Known_peptides[nomod_pep]:
                Known_peptides[nomod_pep][SpecFile] = {}
            if ScanNum not in Known_peptides[nomod_pep][SpecFile]:
                Known_peptides[nomod_pep][SpecFile][ScanNum] = {'Charge':Charge, 'SpecEValue':SpecEValue, 'peptide':peptide}

            exp_type = spectrumfiles[SpecFile]['Dataset']
            protein_level_psms[exp_type].add((SpecFile,ScanNum,peptide))
            
            contaminant_p = True
            for prot in Protein:
                if prot not in contaminant_proteins:
                    contaminant_p = False
                    break
            if contaminant_p==False:
                protein_level_fdr_precursors[exp_type].add(precursor)
                protein_level_fdr_peptides[exp_type].add(nomod_pep)


###########################################################################
print('Write knowndb_protein_matrix_proteins.txt...')
prot_lengths = {}
all_prot_sequences = set()
prot_ref = {}
with open(matrix_protein_fn,'w') as outfile:
    outfile.write('\t'.join(['#Protein','Description','Protein_sequence'])+'\n')
    for p in all_proteins:
        if p in ensembl_proteins:
            prot_ref[p] = 'Ensembl'
            all_prot_sequences.add(ensembl_proteins[p]['seq'])
            plen = len(ensembl_proteins[p]['seq'])*0.01
            prot_lengths[p] = plen
            outfile.write('\t'.join([p, ensembl_proteins[p]['def'],\
                                     ensembl_proteins[p]['seq']])+'\n')
        elif p in refseq_proteins:
            prot_ref[p] = 'RefSeq'
            all_prot_sequences.add(refseq_proteins[p]['seq'])
            plen = len(refseq_proteins[p]['seq'])*0.01
            prot_lengths[p] = plen
            outfile.write('\t'.join([p, refseq_proteins[p]['def'],\
                                     refseq_proteins[p]['seq']])+'\n')
        elif p in contaminant_proteins:
            prot_ref[p] = 'Contaminant'
            all_prot_sequences.add(contaminant_proteins[p]['seq'])
            plen = len(contaminant_proteins[p]['seq'])*0.01
            prot_lengths[p] = plen
            outfile.write('\t'.join([p, contaminant_proteins[p]['def'],\
                                     contaminant_proteins[p]['seq']])+'\n')
        else:
            prot_ref[p] = 'Augustus'
            all_prot_sequences.add(refseq_ensembl_augustus_proteins[p]['seq'])
            plen = len(refseq_ensembl_augustus_proteins[p]['seq'])*0.01
            prot_lengths[p] = plen
            outfile.write('\t'.join([p, refseq_ensembl_augustus_proteins[p]['def'],\
                                     refseq_ensembl_augustus_proteins[p]['seq']])+'\n')

print('Write knowndb_protein_matrix.txt...')
prot_sample = {x:{z:0 for z in specfile_order} for x in all_proteins}
sample = {x:0 for x in specfile_order}
with open(matrix_fn,'w') as outfile:
    outfile.write('#Protein'+'\t'+'\t'.join([spectrumfiles[x]['heatmap_label'] for x in specfile_order])+'\n')
    for prot in all_proteins:
        prot_len = prot_lengths[prot]
        row_values = {x:0 for x in specfile_order}
        for pep in Protein_levelFDR_proteins[prot]:
            num_prot = float(len(Known_peptide_proteins[pep]))
            for spfn in specfile_order:
                if spfn in Known_peptides[pep]:
                    psms = len(Known_peptides[pep][spfn])/num_prot
                    row_values[spfn] += psms
        for spfn in row_values:
            prot_sample[prot][spfn] = row_values[spfn]/prot_len
            sample[spfn]+=prot_sample[prot][spfn]
        outfile.write(prot+'\t'+'\t'.join([str(row_values[x]) for x in specfile_order])+'\n')

print('Computing psk values, writing knowndb_protein_matrix_psk.txt...')
matrix = []
for prot in all_proteins:
    row_value = []
    for spfn in specfile_order:
        psk_value = (prot_sample[prot][spfn]/sample[spfn])*1000
        row_value.append(psk_value)
    matrix.append((sum(row_value), prot+'\t'+'\t'.join([str(x) for x in row_value])+'\n'))    
matrix = sorted(matrix, reverse=True)

with open(matrix_fn_psk,'w') as outfile_psk, open(matrix_fn_psk_top_1000,'w') as top_outfile_psk:
    outfile_psk.write('#Protein'+'\t'+'\t'.join([spectrumfiles[x]['heatmap_label'] for x in specfile_order])+'\n')
    top_outfile_psk.write('#Protein'+'\t'+'\t'.join([spectrumfiles[x]['heatmap_label'] for x in specfile_order])+'\n')
    c=0
    for entry in matrix:
        outfile_psk.write(entry[1])
        c+=1
        if c<=1000:
            top_outfile_psk.write(entry[1])

# print('total proteins: ', len(all_proteins))
# print(len([x for x in prot_ref if prot_ref[x]=='RefSeq']), ' RefSeq')
# print(len([x for x in prot_ref if prot_ref[x]=='Ensembl']), ' Ensembl')
# print(len([x for x in prot_ref if prot_ref[x]=='Augustus']), ' Augustus')
# print(len([x for x in prot_ref if prot_ref[x]=='Contaminant']), ' Contaminant')
# print('total protein sequences: ', len(all_prot_sequences))

# total_precursors = set()
# total_peptides = set()
# print('tissue\tprecursors\tpeptides')
# for t in protein_level_fdr_precursors:
#     print(t, len(protein_level_fdr_precursors[t]), len(protein_level_fdr_peptides[t]))
#     total_precursors.update(protein_level_fdr_precursors[t])
#     total_peptides.update(protein_level_fdr_peptides[t])

# print('total precursors:', len(total_precursors))
# print('total peptides:', len(total_peptides))