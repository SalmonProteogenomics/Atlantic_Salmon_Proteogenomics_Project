# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import os
from operator import itemgetter
import re
from Bio import Entrez

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

reference_dir = os.path.join(base_dir , 'data', 'reference')

Ensembl_dir = os.path.join(reference_dir, 'Ensembl')
RefSeq_dir = os.path.join(reference_dir, 'RefSeq')

Ensembl_gff = os.path.join(Ensembl_dir, 'release-99', 'Salmo_salar.ICSASG_v2.99.chr.gff3')
RefSeq_gff = os.path.join(RefSeq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')

def ncbi_chrname_mapping():
    ncbi_chrname_map = {'NC_027300.1': 'chrssa01', 'NC_027301.1': 'chrssa02', 'NC_027302.1': 'chrssa03',
     'NC_027303.1': 'chrssa04', 'NC_027304.1': 'chrssa05', 'NC_027305.1': 'chrssa06', 'NC_027306.1': 'chrssa07',
     'NC_027307.1': 'chrssa08', 'NC_027308.1': 'chrssa09', 'NC_027309.1': 'chrssa10', 'NC_027310.1': 'chrssa11',
     'NC_027311.1': 'chrssa12', 'NC_027312.1': 'chrssa13', 'NC_027313.1': 'chrssa14', 'NC_027314.1': 'chrssa15',
     'NC_027315.1': 'chrssa16', 'NC_027316.1': 'chrssa17', 'NC_027317.1': 'chrssa18', 'NC_027318.1': 'chrssa19',
     'NC_027319.1': 'chrssa20', 'NC_027320.1': 'chrssa21', 'NC_027321.1': 'chrssa22', 'NC_027322.1': 'chrssa23',
     'NC_027323.1': 'chrssa24', 'NC_027324.1': 'chrssa25', 'NC_027325.1': 'chrssa26', 'NC_027326.1': 'chrssa27',
     'NC_027327.1': 'chrssa28', 'NC_027328.1': 'chrssa29'}
    return ncbi_chrname_map

ncbi_chrname_map = ncbi_chrname_mapping()
chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}

def readFasta(fastafn):
    fasta_proteins = {}
    sequence = []
    accession = ''
    definition = ''
    with open(fastafn,'r') as infile:
        for line in infile:
            if line[0]=='>':
                if sequence:
                    fasta_proteins[accession] = {'seq':''.join(sequence),'def':definition}
                accession = line.strip()[1:].split()[0]
                definition = line.strip()[1:].replace(accession,'').strip()
                sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            fasta_proteins[accession] = {'seq':''.join(sequence),'def':definition}
    return fasta_proteins

def salmosalar_genomic_fna(fna_dir, chromosome):
    with open(os.path.join(fna_dir, chromosome+'_Salmo_salar.fa'),'r') as infile:
        dna_seq = []
        for line in infile:
            if line[0]!='>':
                dna_seq.append(line.strip().upper())
    return ''.join(dna_seq)

def interval_overlap(b, int_list):
    overlap_list = set()
    for a in int_list:
        #true if overlapping intervals
        if a[0] <= b[1] and b[0] <= a[1]:
            s_e = str(a[0])+'-'+str(a[1])
            overlap_list.add(s_e)
    return list(overlap_list)

def combine_indices(start_end_list):
    start_end_list = sorted(start_end_list,key=itemgetter(0))
    if len(start_end_list)==1:
        return start_end_list
    else:
        combined_list = []
        t = list(start_end_list[0])
        for i, entry in enumerate(start_end_list):
            s = int(entry[0])
            e = int(entry[1])
            if s<=t[1] and e>=t[1]:
                t[1] = e
            elif s<=t[1] and e<t[1]:
                continue
            else:
                combined_list.append(t)
                t = [s,e]
        combined_list.append(t)
        return combined_list

def combine_int_into_interval(list_int):
    list_int = sorted(list_int)
    combine_list = []
    s = list_int[0]
    for xi,x in enumerate(list_int):
        if xi==len(list_int)-1:
            combine_list.append((s,x))
            continue
        if list_int[xi+1]-x>1:
            combine_list.append((s,x))
            s = list_int[xi+1]    
    return combine_list

def translate(seq): 
    seq = seq.upper()
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein

def rev_comp(dna_str):
    cmap = {'G':'C','C':'G','A':'T','T':'A','N':'N'}
    return ''.join([cmap[nt] for nt in dna_str.upper()])[::-1]

def comp(dna_str):
    cmap = {'G':'C','C':'G','A':'T','T':'A','N':'N'}
    return ''.join([cmap[nt] for nt in dna_str.upper()])

def find_eventpep(protein_sequence, event_peptides):
    found_peps = {}
    for event_pep in event_peptides: 
        index = [m.start() for m in re.finditer(event_pep.replace(':',''), protein_sequence)]
        if index:
            found_peps[event_pep] = index
    return found_peps

gff3_columns = ['seqid','source','type','start','end','score','strand','phase','attributes']

def read_ensembl_gff(fn):
    ensembl_gff_entries = {}
    ensembl_gff_start_end = {}
    ensembl_mRNA_protein = {}
    ensembl_mRNA_gene = {}
    ensembl_gene_coor = {}
    with open(fn,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                gff_values = {gff3_columns[xi]:x for xi,x in enumerate(sp)} 
                chrname = 'chr'+gff_values['seqid']
                strand = gff_values['strand']
                
                if chrname not in chrname_ncbi_map:
                    continue
                
                if chrname not in ensembl_gff_entries:
                    ensembl_gff_entries[chrname] = {}
                if strand not in ensembl_gff_entries[chrname]:
                    ensembl_gff_entries[chrname][strand] = {}
                
                feature = gff_values['type']
                
                attributes = {x.split('=')[0]:x.split('=')[1] for x in gff_values['attributes'].split(';')}
                ID = ''
                if 'ID' in attributes:
                    ID = attributes['ID']
        
                if feature in ['gene','mRNA','CDS','five_prime_UTR','three_prime_UTR','V_gene_segment','D_gene_segment','J_gene_segment','pseudogene']:
                    start = gff_values['start']
                    end = gff_values['end']
                    key = start+'-'+end
                    
                    if feature in ['gene', 'pseudogene']:
                        if chrname not in ensembl_gff_start_end:
                            ensembl_gff_start_end[chrname] = {}
                        if strand not in ensembl_gff_start_end[chrname]:
                            ensembl_gff_start_end[chrname][strand] = {}
                        
                        if key in ensembl_gff_start_end[chrname][strand]:
                            print(ensembl_gff_start_end[chrname][strand][key])
                            
                        if key not in ensembl_gff_start_end[chrname][strand]:
                            ensembl_gff_start_end[chrname][strand][key] = {}
                        
                        ensembl_gff_start_end[chrname][strand][key] = {'start':start, 'end':end, 'attributes':attributes}
                        if not ID:
                            raise ValueError('')
                        
                        gene_name = ''
                        if 'Name' in attributes:
                            gene_name = attributes['Name']
                        gene_id = ''
                        if 'gene_id' in attributes:
                            gene_id = attributes['gene_id']
                        
                        ensembl_gene_coor[ID] = {'gene_name':gene_name, 'gene_id':gene_id, 'start':start, 'end':end, 'strand': strand, 'chrom':chrname, 'attributes':attributes}
                                
                    elif feature in ['mRNA','V_gene_segment','D_gene_segment','J_gene_segment']:
                        parent_gene = attributes['Parent']
                        ensembl_mRNA_gene[ID] = parent_gene
                        if parent_gene not in ensembl_gff_entries[chrname][strand]:
                            ensembl_gff_entries[chrname][strand][parent_gene] = {}#transcript ID as key
                        if ID not in ensembl_gff_entries[chrname][strand][parent_gene]:
                            ensembl_gff_entries[chrname][strand][parent_gene][ID] = {'strand':strand, 'start':start, 'end':end, 'attributes':attributes, 'CDS':[], '5UTR':[], '3UTR':[]}
                            
                    elif feature == 'CDS':
                        parent_type = attributes['Parent'].split(':')[0]
                        parent_transcript = attributes['Parent']
                        parent_gene = ''
                        if parent_type == 'transcript':  
                            parent_gene = ensembl_mRNA_gene[parent_transcript]
                        else:
                            print('parent_type not transcript')
                            print(attributes)
                        protein_id = attributes['protein_id']
    
                        pt = attributes['Parent'].split(':')[1]
                        transcript_name = ''
                        if pt not in ensembl_mRNA_protein:
                            ensembl_mRNA_protein[pt] = set()
                        if 'Name' in ensembl_gff_entries[chrname][strand][parent_gene][parent_transcript]['attributes']:
                            transcript_name = ensembl_gff_entries[chrname][strand][parent_gene][parent_transcript]['attributes']['Name']                        
                            if transcript_name not in ensembl_mRNA_protein:
                                ensembl_mRNA_protein[transcript_name] = set()
          
                        if protein_id:
                            ensembl_mRNA_protein[pt].add(protein_id)
                            if transcript_name:
                                ensembl_mRNA_protein[transcript_name].add(protein_id)
                        
                        ensembl_gff_entries[chrname][strand][parent_gene][parent_transcript]['CDS'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes, 'protein_id':protein_id})
                        
                    elif feature == 'five_prime_UTR':
                        if not ID:
                            ID = '5utr'
                        parent_type = attributes['Parent'].split(':')[0]
                        parent_transcript = attributes['Parent']
                        parent_gene = ''
                        if parent_type == 'transcript':
                            parent_gene = ensembl_mRNA_gene[parent_transcript]
                            ensembl_gff_entries[chrname][strand][parent_gene][parent_transcript]['5UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes, 'ID':ID})
                        else:
                            print('parent_type not transcript')
                            print(attributes)
                         
                    elif feature == 'three_prime_UTR':
                        if not ID:
                            ID = '3utr'
                        parent_type = attributes['Parent'].split(':')[0]
                        parent_transcript = attributes['Parent']
                        parent_gene = ''
                        if parent_type == 'transcript':
                            parent_gene = ensembl_mRNA_gene[parent_transcript]
                            ensembl_gff_entries[chrname][strand][parent_gene][parent_transcript]['3UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes, 'ID':ID})
                        else:
                            print('parent_type not transcript')
                            print(attributes)

    return {'ensembl_gff_entries':ensembl_gff_entries, 'ensembl_gff_start_end': ensembl_gff_start_end,
            'ensembl_mRNA_gene':ensembl_mRNA_gene, 'ensembl_mRNA_protein':ensembl_mRNA_protein,
            'ensembl_gene_coor':ensembl_gene_coor}


def read_refseq_gff(fn):
    refseq_gff_entries = {}
    refseq_gff_start_end = {}
    refseq_mRNA_gene = {}
    refseq_mRNA_protein = {}
    refseq_ID_type = {}
    refseq_gene_coor = {}
    refseq_protein_cds = {}
    with open(fn,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                gff_values = {gff3_columns[xi]:x for xi,x in enumerate(sp)} 
                chrname = gff_values['seqid']
                strand = gff_values['strand']
                
                if chrname not in ncbi_chrname_map:
                    continue
                
                chrname = ncbi_chrname_map[chrname]
                
                if chrname not in refseq_gff_entries:
                    refseq_gff_entries[chrname] = {}
                if strand not in refseq_gff_entries[chrname]:
                    refseq_gff_entries[chrname][strand] = {}
                
                feature = gff_values['type']            
                attributes = {x.split('=')[0]:x.split('=')[1] for x in gff_values['attributes'].split(';')}
                            
                if feature in ['gene','mRNA', 'transcript', 'CDS','five_prime_UTR','three_prime_UTR','V_gene_segment','C_gene_segment','pseudogene']:                
                    start = gff_values['start']
                    end = gff_values['end']
                    key = start+'-'+end
                    ID = ''                    
                    if feature in ['gene', 'pseudogene']:
                        ID = attributes['ID']
                        refseq_ID_type[ID] = feature
                            
                        if chrname not in refseq_gff_start_end:
                            refseq_gff_start_end[chrname] = {}
                        if strand not in refseq_gff_start_end[chrname]:
                            refseq_gff_start_end[chrname][strand] = {}
                        
                        if key in refseq_gff_start_end[chrname][strand]:
                            print('key in refseq_gff_start_end[chrname][strand]:', refseq_gff_start_end[chrname][strand][key])
                            
                        if key not in refseq_gff_start_end[chrname][strand]:
                            refseq_gff_start_end[chrname][strand][key] = {}
                        
                        refseq_gff_start_end[chrname][strand][key] = {'start':start, 'end':end, 'attributes':attributes}
                        
                        if not ID:
                            raise ValueError('')

                        gene_name = ''
                        if 'Name' in attributes:
                            gene_name = attributes['Name']
                        gene_id = ''
                        if 'gene_id' in attributes:
                            gene_id = attributes['gene_id']
                            
                        if ID not in refseq_gff_entries[chrname][strand]:
                            refseq_gff_entries[chrname][strand][ID] = {}
        
                        refseq_gene_coor[ID] = {'gene_name':gene_name, 'gene_id':gene_id, 'start':start, 'end':end, 'strand': strand, 'chrom':chrname, 'attributes':attributes}
                                
                    elif feature in ['mRNA','transcript']:
                        ID = attributes['ID']
                        refseq_ID_type[ID] = feature
                        
                        parent_type = attributes['Parent'].split('-')[0]
                        parent_gene = attributes['Parent']
                            
                        refseq_mRNA_gene[ID] = parent_gene
                        if ID not in refseq_gff_entries[chrname][strand][parent_gene]:
                            refseq_gff_entries[chrname][strand][parent_gene][ID] = {'strand':strand, 'start':start, 'end':end, 'attributes':attributes, 'CDS':[], '5UTR':[], '3UTR':[]}
                            
                    elif feature == 'CDS':
                        ID = attributes['ID']                            
                        refseq_ID_type[ID] = feature
    
                        protein_id = ''
                        if 'protein_id' in attributes:
                            protein_id = attributes['protein_id']
                        product = attributes['product']
                        
                        parent_type = attributes['Parent'].split('-')[0]
                        parent_transcript = attributes['Parent']
                        
                        pt = attributes['Parent']
                        if pt not in refseq_mRNA_protein:
                            refseq_mRNA_protein[pt] = set()
    
                        if protein_id:
                            refseq_mRNA_protein[pt].add(protein_id)
                        
                        parent_gene = ''
                        if parent_type == 'rna':  
                            parent_gene = refseq_mRNA_gene[parent_transcript]
                            refseq_gff_entries[chrname][strand][parent_gene][parent_transcript]['CDS'].append({'strand': strand, 'start':start, 'end':end, \
                                                                                                              'attributes':attributes, 'protein_id':protein_id, 'product':product})                                           
                            transcript_name = ''
                            if 'Name' in refseq_gff_entries[chrname][strand][parent_gene][parent_transcript]['attributes']:
                                transcript_name = refseq_gff_entries[chrname][strand][parent_gene][parent_transcript]['attributes']['Name']
                            if protein_id:
                                if transcript_name not in refseq_mRNA_protein:
                                    refseq_mRNA_protein[transcript_name] = set()
                                refseq_mRNA_protein[transcript_name].add(protein_id)
                                if protein_id not in refseq_protein_cds:
                                    refseq_protein_cds[protein_id] = {'CDS':[]}
                                refseq_protein_cds[protein_id]['CDS'].append({'chromosome':chrname, 'strand':strand, 'start':start, 'end':end, 'attributes':attributes,\
                                                                              'parent_gene':parent_gene, 'parent_transcript':parent_transcript})                        
                        elif parent_type == 'gene':
                            parent_gene = attributes['Parent']
                            if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                            if 'CDS' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['CDS'] = []
                            refseq_gff_entries[chrname][strand][parent_gene]['other']['CDS'].append({'strand': strand, 'start':start, 'end':end, \
                                                                                                      'attributes':attributes, 'protein_id':protein_id, 'product':product})
                        elif parent_type == 'id':
                            if refseq_ID_type[attributes['Parent'].replace('id-','gene-')]=='gene':
                                parent_gene = attributes['Parent'].replace('id-','gene-')
                                if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                                if 'CDS' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other']['CDS'] = []
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['CDS'].append({'strand': strand, 'start':start, 'end':end, \
                                                                                                          'attributes':attributes, 'protein_id':protein_id, 'product':product})
                        else:        
                            print('parent_type not transcript', parent_type)
    
                    elif feature == 'five_prime_UTR':
                        ID = attributes['ID'] 
                        refseq_ID_type[ID] = feature                        
                        parent_type = attributes['Parent'].split('-')[0]
                        parent_transcript = attributes['Parent']                        
                        parent_gene = ''
                        if parent_type == 'rna':
                            parent_gene = refseq_mRNA_gene[parent_transcript]
                            refseq_gff_entries[chrname][strand][parent_gene][parent_transcript]['5UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes})                            
                            if parent_transcript in refseq_mRNA_protein:
                                protid = list(refseq_mRNA_protein[parent_transcript])[0]
                                if protid in refseq_protein_cds:
                                    if '5UTR' not in refseq_protein_cds:
                                        refseq_protein_cds[protid]['5UTR'] = []
                                    refseq_protein_cds[protid]['5UTR'].append({'chromosome':chrname, 'strand':strand, 'start':start, 'end':end, 'attributes':attributes,\
                                                                              'parent_gene':parent_gene, 'parent_transcript':parent_transcript})
                        elif parent_type == 'gene':
                            parent_gene = attributes['Parent']
                            if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                            if '5UTR' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['5UTR'] = []
                            refseq_gff_entries[chrname][strand][parent_gene]['other']['5UTR'].append({'strand': strand,'start':start, 'end':end, 'attributes':attributes})
                        elif parent_type == 'id':
                            if refseq_ID_type[attributes['Parent'].replace('id-','gene-')]=='gene':
                                parent_gene = attributes['Parent'].replace('id-','gene-')
                                if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                                if '5UTR' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other']['5UTR'] = []
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['5UTR'].append({'strand': strand,'start':start, 'end':end, 'attributes':attributes})
                        else:        
                            print('parent_type not transcript', parent_type)
    
                    elif feature == 'three_prime_UTR':
                        ID = attributes['ID']
                        refseq_ID_type[ID] = feature                        
                        parent_type = attributes['Parent'].split('-')[0]
                        parent_transcript = attributes['Parent']                        
                        parent_gene = ''
                        if parent_type == 'rna':
                            parent_gene = refseq_mRNA_gene[parent_transcript]
                            refseq_gff_entries[chrname][strand][parent_gene][parent_transcript]['3UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes})
                            if parent_transcript in refseq_mRNA_protein:
                                protid = list(refseq_mRNA_protein[parent_transcript])[0]
                                if protid in refseq_protein_cds:
                                    if '3UTR' not in refseq_protein_cds:
                                        refseq_protein_cds[protid]['3UTR'] = []
                                    refseq_protein_cds[protid]['3UTR'].append({'chromosome':chrname, 'strand':strand, 'start':start, 'end':end, 'attributes':attributes, \
                                                                              'parent_gene':parent_gene, 'parent_transcript':parent_transcript})
                        elif parent_type == 'gene':
                            parent_gene = attributes['Parent']
                            if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                            if '3UTR' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['3UTR'] = []
                            refseq_gff_entries[chrname][strand][parent_gene]['other']['3UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes})
                        elif parent_type == 'id':
                            if refseq_ID_type[attributes['Parent'].replace('id-','gene-')]=='gene':
                                parent_gene = attributes['Parent'].replace('id-','gene-')
                                if 'other' not in refseq_gff_entries[chrname][strand][parent_gene]:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other'] = {}
                                if '3UTR' not in refseq_gff_entries[chrname][strand][parent_gene]['other']:
                                    refseq_gff_entries[chrname][strand][parent_gene]['other']['3UTR'] = []
                                refseq_gff_entries[chrname][strand][parent_gene]['other']['3UTR'].append({'strand': strand, 'start':start, 'end':end, 'attributes':attributes})
                        else:        
                            print('parent_type not transcript', parent_type    )

    return {'refseq_gff_entries':refseq_gff_entries, 'refseq_gff_start_end':refseq_gff_start_end,
            'refseq_mRNA_gene':refseq_mRNA_gene, 'refseq_mRNA_protein':refseq_mRNA_protein,
            'refseq_gene_coor':refseq_gene_coor, 'refseq_protein_cds':refseq_protein_cds}


def reference_gene_overlap(reference_se, reference_gff, event_min_loc,\
                           event_max_loc, events_splice_loc, events_other_loc,\
                              ensembl_mRNA_protein, ensembl_proteins):
    #overlapping gene
    overlapping_gene = interval_overlap((event_min_loc, event_max_loc), sorted([sorted((int(reference_se[key]['start']), int(reference_se[key]['end']))) for key in reference_se]))
    #splice events at CDS region?
    ref_gene_supported_by_event = {}
    for rg in overlapping_gene:
        ref_gene_id = reference_se[rg]['attributes']['ID']
        
        ref_gene_name = ''
        if 'Name' in reference_se[rg]['attributes']:
            ref_gene_name = reference_se[rg]['attributes']['Name']
        if not ref_gene_name:
            ref_gene_name = reference_se[rg]['attributes']['ID']

        for rt in reference_gff[ref_gene_id]:
            if rt=='other':
                continue            
            CDS = reference_gff[ref_gene_id][rt]['CDS']
            transcript_name = ''
            if 'Name' in reference_gff[ref_gene_id][rt]['attributes']:
                transcript_name = reference_gff[ref_gene_id][rt]['attributes']['Name']
            if not transcript_name:
                transcript_name = rt
            cds_loc = sorted([sorted((int(cds['start']), int(cds['end']))) for cds in CDS])
            #
            transcriptID = transcript_name
            if 'transcript:' in transcript_name:
                transcriptID = transcript_name.split('transcript:')[1]
            if transcriptID in ensembl_mRNA_protein:
                protein_id = ensembl_mRNA_protein[transcriptID]
                if protein_id:
                    for p in protein_id:
                        protid_v = p+'.1'
                        prot_seq = ensembl_proteins[protid_v]['seq']
                        for pep in events_other_loc:
                            if pep.replace(':','') in prot_seq:
                                if ref_gene_name not in ref_gene_supported_by_event:
                                    ref_gene_supported_by_event[ref_gene_name] = {}
                                if transcript_name not in ref_gene_supported_by_event[ref_gene_name]:
                                    ref_gene_supported_by_event[ref_gene_name][transcript_name] = set()
                                ref_gene_supported_by_event[ref_gene_name][transcript_name].add(pep)
            #
            CDS_borders = []            
            if len(cds_loc)==1:
                continue
            if len(cds_loc)==2:
                CDS_borders = [(cds_loc[0][1], cds_loc[1][0])]
            if len(cds_loc)>2:
                CDS_borders = []
                for i,il in enumerate(cds_loc):
                    if i==len(cds_loc)-1:
                        continue
                    CDS_borders.append((cds_loc[i][1], cds_loc[i+1][0]))
            CDS_borders = set(CDS_borders)
            for pep in events_splice_loc:
                for pep_locations in events_splice_loc[pep]:
                    overlapping = [el for el in pep_locations if el in CDS_borders]
                    if len(overlapping)==len(pep_locations):
                        if ref_gene_name not in ref_gene_supported_by_event:
                            ref_gene_supported_by_event[ref_gene_name] = {}
                        if transcript_name not in ref_gene_supported_by_event[ref_gene_name]:
                            ref_gene_supported_by_event[ref_gene_name][transcript_name] = set()
                        ref_gene_supported_by_event[ref_gene_name][transcript_name].add(pep)

    return ref_gene_supported_by_event

def reference_augustus_overlap(reference_gff, augustus_cds_SE, prot_id):
    ref_gene_list = {}
    parent_transcripts = set()
    parent_gene = set()
    
    #CDS
    CDS_entries = {}
    CDS_se_list = []
    if 'CDS' in reference_gff:
        for cds in reference_gff['CDS']:
            cds_s_e = sorted([int(cds['start']), int(cds['end'])])
            parent_transcripts.add(cds['parent_transcript'])
            parent_gene.add(cds['parent_gene'])
            if cds_s_e not in CDS_se_list:
                CDS_se_list.append(cds_s_e)
            cds_s_e = str(cds_s_e[0])+'-'+str(cds_s_e[1])
            CDS_entries[cds_s_e] = cds

    #5UTR
    fiveUTR_entries = {}
    fiveUTR_se_list = []
    if '5UTR' in reference_gff:
        for utr in reference_gff['5UTR']:
            utr_s_e = sorted([int(utr['start']), int(utr['end'])])
            parent_transcripts.add(utr['parent_transcript'])
            parent_gene.add(utr['parent_gene'])
            if utr_s_e not in fiveUTR_se_list:
                fiveUTR_se_list.append(utr_s_e)
            utr_s_e = str(utr_s_e[0])+'-'+str(utr_s_e[1])
            fiveUTR_entries[utr_s_e] = utr

    #3UTR
    threeUTR_entries = {}
    threeUTR_se_list = []
    if '3UTR' in reference_gff:
        for utr in reference_gff['3UTR']:
            utr_s_e = sorted([int(utr['start']), int(utr['end'])])
            parent_transcripts.add(utr['parent_transcript'])
            parent_gene.add(utr['parent_gene'])
            if utr_s_e not in threeUTR_se_list:
                threeUTR_se_list.append(utr_s_e)
            utr_s_e = str(utr_s_e[0])+'-'+str(utr_s_e[1])
            threeUTR_entries[utr_s_e] = utr
    
    if len(parent_transcripts)!=1:
        raise ValueError('')

    rt = list(parent_transcripts)[0]
    rg = list(parent_gene)[0]
    # reference CDS,UTR and augustus CDS overlap
    for aug_cds_i, aug_cds in enumerate(augustus_cds_SE):
        overlapping_CDS = [CDS_entries[x] for x in interval_overlap((aug_cds[0], aug_cds[1]), sorted(CDS_se_list))]
        overlapping_5UTR = [fiveUTR_entries[x] for x in interval_overlap((aug_cds[0], aug_cds[1]), sorted(fiveUTR_se_list))]
        overlapping_3UTR = [threeUTR_entries[x] for x in interval_overlap((aug_cds[0], aug_cds[1]), sorted(threeUTR_se_list))]

        if not (overlapping_CDS or overlapping_5UTR or overlapping_3UTR):
            continue
        
        ref_gene_list['parent_transcript'] = rt
        ref_gene_list['parent_gene'] = rg
        
        for CDS in overlapping_CDS:
            ref_Id = CDS['attributes']['ID']
            ref_s_e = sorted([int(CDS['start']), int(CDS['end'])])
            if 'CDS' not in ref_gene_list:
                ref_gene_list['CDS'] = []            
            ref_gene_list['CDS'].append({'protein_id':prot_id, 'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})

        for utr in overlapping_5UTR:
            ref_Id = '5utr'
            if 'ID' in utr['attributes']:
                ref_Id = utr['attributes']['ID']                   
            ref_s_e = sorted([int(utr['start']), int(utr['end'])])
            if '5UTR' not in ref_gene_list:
                ref_gene_list['5UTR'] = []
            ref_gene_list['5UTR'].append({'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})
            
        for utr in overlapping_3UTR:
            ref_Id = '3utr'
            if 'ID' in utr['attributes']:
                ref_Id = utr['attributes']['ID']
            ref_s_e = sorted([int(utr['start']), int(utr['end'])])
            if '3UTR' not in ref_gene_list:
                ref_gene_list['3UTR'] = []
            ref_gene_list['3UTR'].append({'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})

    return ref_gene_list

def read_ref_prot_genomic(fn):
    known_prot_genomic_coor = {}
    with open(fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                #['#Chromosome','strand','protein', 'genomic_coordinates']
                colidx = {x:xi for xi,x in enumerate(sp)}
            else:
                chrom = sp[colidx['#Chromosome']]
                strand = sp[colidx['strand']]
                protein = sp[colidx['protein']].split('|')[0]
                genomic_coordinates = sp[colidx['genomic_coordinates']].split(';')
                if protein not in known_prot_genomic_coor:
                    known_prot_genomic_coor[protein] = {}
                known_prot_genomic_coor[protein] = {'chrom':chrname_ncbi_map[chrom], \
                                                           'coor':{xi:[int(z) for z in x.split(':')[1].split(',')] for xi,x in enumerate(genomic_coordinates)},
                                                           'strand':strand}
    return known_prot_genomic_coor
    
def read_ref_prot_cds(fn):
    known_prot_genomic_coor_incomplete = {}
    with open(fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                #['#Chromosome','strand','protein', 'CDS_coordinates']
                colidx = {x:xi for xi,x in enumerate(sp)}
            else:
                chrom = sp[colidx['#Chromosome']]
                strand = sp[colidx['strand']]
                protein = sp[colidx['protein']].split('|')[0]
                CDS_coordinates = sp[colidx['CDS_coordinates']].split(';')
                if protein not in known_prot_genomic_coor_incomplete:
                    known_prot_genomic_coor_incomplete[protein] = {}
                #1-index cds coor
                known_prot_genomic_coor_incomplete[protein] = {'chrom':chrname_ncbi_map[chrom], \
                                                           'coor':[(int(x.split(':')[1].split('-')[0]),int(x.split(':')[1].split('-')[1])) for x in CDS_coordinates],
                                                           'strand':strand}    
    return known_prot_genomic_coor_incomplete
    
def read_ref_prot_gff(fn):
    Salmo_salar_ref_prot_gff = {}
    with open(fn,'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            if line.strip().split()[0]=='LOCUS':
                prot_def = ''
                prot_acc = ''
                CDS_start = False
                gene_source = False 
            if line.strip().split()[0] == 'DEFINITION':
                prot_def = line.strip().split('DEFINITION')[1].strip() 
            if line.strip().split()[0] == 'VERSION':
                prot_acc = line.strip().split('VERSION')[1].strip()
                Salmo_salar_ref_prot_gff[prot_acc] = {'def':prot_def,'dbsource':'',
                                                        'chrom':'','gene':'',
                                                        'coded_by':'','db_xref':''}   
            if line.strip().split()[0] == 'DBSOURCE':
                Salmo_salar_ref_prot_gff[prot_acc]['dbsource'] = line.strip().split('DBSOURCE')[1].strip()         
            if 'source          1..' in line:
                gene_source = True   
            if 'CDS             1..' in line:
                CDS_start = True
            if gene_source==True:
                if '/chromosome="' in line:
                    Salmo_salar_ref_prot_gff[prot_acc]['chrom'] = 'chr'+line.strip().split('/chromosome="')[1].split('"')[0]  
            if CDS_start==True:
                if '/gene="' in line:
                    Salmo_salar_ref_prot_gff[prot_acc]['gene'] = line.strip().split('/gene="')[1].split('"')[0]
                if '/coded_by="' in line:
                    Salmo_salar_ref_prot_gff[prot_acc]['coded_by'] = line.strip().split('/coded_by="')[1].split('"')[0]         
                if '/db_xref="' in line: 
                    Salmo_salar_ref_prot_gff[prot_acc]['db_xref'] = line.strip().split('/db_xref="')[1].split('"')[0]
    return Salmo_salar_ref_prot_gff

def ncbi_entrez_genpept(accession_str, GenPept_dir, entrez_email):
    Entrez.email = entrez_email
    genpept_fn = os.path.join(GenPept_dir, accession_str+'.txt')
    if '|' in accession_str:
        genpept_fn = os.path.join(GenPept_dir, accession_str.replace('|','-')+'.txt')

    entrez_lines = ''
    
    if not os.path.isfile(genpept_fn):
        try:
            handle = Entrez.efetch(db="protein", id=accession_str, rettype="gp", retmode= "text")
            entrez_lines = handle.readlines()
            handle.close()
        except:
            print('Could not retrieve: ', accession_str)
            return 'Could not retrieve'

        with open(genpept_fn,'w') as outf:
            outf.write(''.join(entrez_lines))
    
    TAGS = ['DEFINITION','ACCESSION','VERSION','DBSOURCE','SOURCE','ORGANISM','FEATURES', 'CDS', 'ORIGIN']
    ABC_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integers = set([x for x in '0123456789'])
    ABC_upper = set([x for x in ABC_upper])

    '''
    The location of each feature is provided as well, and can be a single base, a contiguous span of bases, a joining of sequence spans, and other representations. 
    If a feature is located on the complementary strand, the word "complement" will appear before the base span. 
    If the "<" symbol precedes a base span, the sequence is partial on the 5' end (e.g., CDS  <1..206).  
    If the ">" symbol follows a base span, the sequence is partial on the 3' end (e.g., CDS   435..915>).
    Base span of the biological feature indicated to the left, in this case, a CDS feature. 
    (The CDS feature is described above, and its base span includes the start and stop codons.) 
    Features can be complete, partial on the 5' end, partial on the 3' end, and/or on the complementary strand. 
    '''
    with open(genpept_fn,'r') as entrez_lines:
        tag_lines = {x:[] for x in TAGS}
        tag_begin = False
        for line in entrez_lines:
            if line.strip():
                sp = set([x for x in line.split()[0]])-ABC_upper
                if not sp:
                    tag = line.split()[0]
                    tag_begin = False
                    if tag in TAGS:                    
                        tag_lines[tag].append(line.strip())
                        tag_begin = True
                else:
                    if tag_begin == True:
                        tag_lines[tag].append(line.strip())
    
    genpept_features = tag_lines['FEATURES']
    '''
    FEATURES             Location/Qualifiers
     source          1..262
                     /organism="Salmo salar"
                     /db_xref="taxon:8030"
                     /chromosome="ssa23"
    '''

    feature_index = {}
    chromosome = ''
    feature_count = 0
    source_fc = ''
    for i in range(1,len(genpept_features)):     
        current_line = genpept_features[i]
        if current_line.strip()[0]!='/':
            if len(current_line.strip().split())==2 and not set([x for x in current_line.strip().split()[1].replace('.','')])-integers:   
                loc_tag = current_line.strip().split()[0]
                feature_count+=1
                feature_index[feature_count] = {'feature':loc_tag, 'index':i}
                if loc_tag == 'source':
                    source_fc = feature_count
        else:
            if '/chromosome=' in current_line:
                chromosome = current_line.strip().split('/chromosome=')[1].replace('"','')
    
    if not feature_index:
        print('Unexpected source line: ', accession_str)
        return 'Unexpected source line'
    
    if feature_count==1:
        source_lines = ' '.join(genpept_features[feature_index[source_fc]['index']:len(genpept_features)][1:]).split('"')
    else:
        source_lines = ' '.join(genpept_features[feature_index[source_fc]['index']:feature_index[source_fc+1]['index']][1:]).split('"')
    
    source_lines = ''.join([x.replace('/','|') if '=' not in x else x for x in source_lines if x.strip()]).split('/')
    source_lines = {x.strip().split('=')[0]:x.strip().split('=')[1].strip() for x in source_lines if x.strip() and '=' in x}
    
    DEFINITION = ' '.join(tag_lines['DEFINITION']).lstrip('DEFINITION').strip()
    #ACCESSION = ' '.join(tag_lines['ACCESSION']).lstrip('ACCESSION').strip().split(' ')
    VERSION = ' '.join(tag_lines['VERSION']).lstrip('VERSION').strip().split(' ')
    ORGANISM = tag_lines['ORGANISM']
    ORGANISM = ORGANISM[0].replace('ORGANISM','').strip()    
    
    DBSOURCE = ''
    if tag_lines['DBSOURCE']:
        DBSOURCE = ' '.join(tag_lines['DBSOURCE']).lstrip('DBSOURCE').strip()

    gene = 'NA'
    gene_id = 'NA'
    coded_by = 'NA'
    if tag_lines['CDS']:
        CDS = ' '.join(tag_lines['CDS']).split('"')
        CDS = ''.join([x.replace('/','|') if '=' not in x else x for x in CDS if x.strip()]).split('/')
        CDS[0] = CDS[0].replace('CDS','CDS=').strip()
        CDS = [x.replace('db_xref=','').replace(':','=') if 'db_xref=' in x else x for x in CDS]
        CDS = {x.strip().split('=')[0]:''.join(x.strip().split('=')[1:]).replace(x.strip().split('=')[0],'').strip() for x in CDS}
        if 'gene' in CDS:
            gene = CDS['gene']
        if 'GeneID' in CDS:
            gene_id = CDS['GeneID']
        if 'coded_by' in CDS:
            if 'join' in CDS['coded_by']:  
                coded_by = []
                CDS_parts = {}
                for x in CDS['coded_by'].split('join(')[1].split(')')[0].replace(' ','').split(','):
                    mrna, mrna_loc = x.split(':')
                    if mrna not in CDS_parts:
                        CDS_parts[mrna] = []
                    CDS_parts[mrna].append(mrna_loc)
                for x in CDS_parts:
                    coded_by.append(x+':'+','.join(CDS_parts[mrna]))
                coded_by = '; '.join(coded_by)
            else:
                coded_by = CDS['coded_by']
            
    SEQUENCE = ''
    if tag_lines['ORIGIN']:
        SEQUENCE = ''.join([''.join([aa.upper() for aa in x if aa!=' ' and aa not in integers]) for x in tag_lines['ORIGIN'][1:-1]])

    '''
    XP_054594841.1
    {'mrna': 'XM_054738866.1:352..2493', 'gene': 'LOC129162763',
     'gene_id': '129162763', 
     'DEFINITION': 'uncharacterized protein LOC129162763 [Nothobranchius furzeri].', 
     'VERSION': ['XP_054594841.1'], 'ORGANISM': 'Nothobranchius furzeri',
     'SEQUENCE': 'MYQNLNQKQACVFYAVRDWCIKRVCGLNPDPFFFFLEGGAGTGKSMVVRCIHSEASKILSRLPAEDVDLSNPTVLLTSFTGTAAFNIGGTTLHSLLKLPRSLKPPIQGLGNQLDEVRCELLNAEILVIDEISMVSKPLFAYVDARLKQIKGNQRPFGGMSVLAVGDFYQLPPVRQSKPLCVYDPSDIDLWQPYFQMASLTEIMRQKDDVAFAEMLNRIRVKEKTDELSPADRDMLSRTITEPELCPSDVLHIFATNKDVEAHNSATLERLHDNIITIDADDFQKDPRTGRMERKDTPLKGGRGELSDCLKVAEGARVMLTRNINVQDGLVNGAFGKLVRVIASEIDPQHIFKLALRMDNQSSVRSNRHGASGSDDLVYIERAEDSLKQRGGVRRQFPVRLAFSCTTHKTQGLTTHAAVVSLKNIFEPGMAYVALSRVTSLGGLYLLGLDERKIYANPDVTAALQSMRQASVDQMMPLLLVREAVSRPDTLTLIHHNTEGLPAHINDITSHHELSLADVLCLTETHLQGSFVAESLVLDGYTMFKRSRHQSYTKFPQMACKSGGGVAVYVKNHIHVREKRYVHNVTDLEFLVLKVETPFPALIAVIYRPPDYIMRPFMENLVSLLDSLEVMDCHPVIVCGDFNENQFSGGRKQIVEHFQSRGYAQMITSATTDKNTLLDLVFISQPQRSLHCGVLRTYYSYHDPVFCVLSSSQP',
     'DBSOURCE': 'REFSEQ: accession XM_054738866.1'}
    '''

    return {'chromosome':chromosome, 'mrna':coded_by,'gene':gene, 'gene_id':gene_id,\
            'DEFINITION':DEFINITION, 'VERSION':VERSION, 'ORGANISM':ORGANISM,\
            'SEQUENCE':SEQUENCE, 'DBSOURCE':DBSOURCE}