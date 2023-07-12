# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import os
from parse_reference_files import salmosalar_genomic_fna, combine_indices

def generate_blastx_hints(blastx_dir, reference_dir):
    blastx_input_dir = os.path.join(blastx_dir, 'input')
    blastx_output_dir = os.path.join(blastx_dir, 'output') 
    
    # refseq_dir = os.path.join(reference_dir, 'RefSeq')
    # chromosome_genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')    
    
    hsp_cutoff = 1e-10
 
    blastx_hints = {}
    EG_blastx_query_range = {}
    for blastx_outfn in os.listdir(blastx_output_dir):
        chromosome = '' # blastx_outfn.split('_')[0]
        EG = '' #blastx_outfn.split('.blastx.txt')[0].split(chromosome+'_')[1]    
        all_hits = {}
        hit = ''
        with open(os.path.join(blastx_output_dir, blastx_outfn),'r') as infile:
            for line in infile:
                if '<BlastOutput_query-def>' in line:
                    #<BlastOutput_query-def>chrssa06_EG_1</BlastOutput_query-def>
                    tag = line.strip().split('<')[1].split('>')[0]
                    blastx_input_fn = line.strip().split('<'+tag+'>')[1].split('</'+tag+'>')[0]
                    chromosome = blastx_input_fn.split('_')[0]
                    EG = blastx_input_fn
                    #>chrssa15_EG_5	chrssa15[102255703-102350363]	FVECFIAEQNMV:SVAIGCATR/102256053-102256085;102349984-102350013|IAEQNMV:SVAIGCATR/102256068-102256085;102349984-102350013
                    input_header = ''
                    blastin_seq = []
                    with open(os.path.join(blastx_input_dir, blastx_input_fn+'.blastx.txt'),'r') as blastinput:
                        for query in blastinput:
                            if query[0]=='>':
                                input_header = query[1:].split('\t')
                                #dna_seq[EG_group_min-window_bp:EG_group_max+window_bp] 0-index
                                chromosome_location = input_header[1].split('[')[1].split(']')[0].split('-')
                            else:
                                blastin_seq.append(query.strip())
                    blastin_seq = ''.join(blastin_seq)
                if '<Hit>' in line:
                    hit = {'Hsp':[]}
                    continue            
                if hit:
                    if '<Hit_hsps>' in line or '</Hit_hsps>' in line:
                        continue                
                    if '<Hsp>' in line:
                        hsp = {}
                        continue    
                    if '</Hsp>' in line:
                        hit['Hsp'].append(hsp)
                        hsp = {}
                        continue                
                    if '</Hit>' in line:
                        all_hits[hit['Hit_num']] = hit
                        hit = ''
                        continue    
                    tag = line.strip().split('<')[1].split('>')[0]
                    text = line.strip().split('<'+tag+'>')[1].split('</'+tag+'>')[0]                
                    if 'Hsp_' in tag:
                        hsp[tag] = text
                    else:
                        hit[tag] = text                    
    
        #chr_dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, chromosome)
        
        low_evalue_query_coor = {'+':set(),'-':set()}
        min_query_coor = ''
        max_query_coor = ''
        no_hits = True
        for hit in all_hits.values():
            Hit_num = hit['Hit_num']
            Hit_hsps = hit['Hsp']
            for hsp in Hit_hsps:
                Hsp_num = hsp['Hsp_num']
                #query 
                Hsp_query_from = int(hsp['Hsp_query-from'])
                Hsp_query_to = int(hsp['Hsp_query-to'])
                if Hsp_query_to<Hsp_query_from:
                    raise ValueError('Hsp_query_to<Hsp_query_from')
                blast_input_query = blastin_seq[int(hsp['Hsp_query-from'])-1:int(hsp['Hsp_query-to'])]
                Hsp_query_length = Hsp_query_to-Hsp_query_from
                if Hsp_query_length/abs(Hsp_query_length)==-1:
                    raise ValueError(Hsp_query_length)
                Hsp_query_frame = int(hsp['Hsp_query-frame'])
                strand = '+'
                if Hsp_query_frame/abs(Hsp_query_frame)==-1:
                    strand = '-'
    
                Hsp_query_from = int(chromosome_location[0])+(Hsp_query_from-1)+1 #1 index
                Hsp_query_to = Hsp_query_from+Hsp_query_length-1 #1 index                
            
    #            chr_subsequence = chr_dna_seq[Hsp_query_from-1:Hsp_query_to+1] #last coordinate is inclusive
    #            Hsp_qseq = hsp['Hsp_qseq'].replace('-','')
    #            if strand == '+':
    #                expected_seq = str(Seq(chr_subsequence).translate().rstrip('*'))
    #            else:
    #                expected_seq = str(Seq(rev_comp(chr_subsequence)).translate().rstrip('*'))
    #            
    #            if Hsp_qseq!= expected_seq:
    #                raise ValueError('Hsp_qseq!= expected_seq')
    
                #evalue
                Hsp_evalue = float(hsp['Hsp_evalue'])
                #used 0.001 for cutoff when running blastx
                if Hsp_evalue<=hsp_cutoff:
                    low_evalue_query_coor[strand].add((Hsp_query_from, Hsp_query_to))
                    no_hits = False
                    if not min_query_coor:
                        min_query_coor = Hsp_query_from
                        max_query_coor = Hsp_query_to
                    else:
                        if Hsp_query_from<min_query_coor:
                            min_query_coor = Hsp_query_from
                        if Hsp_query_to>max_query_coor:
                            max_query_coor = Hsp_query_to                    
                        
        if no_hits==True:
            continue
        
        if not min_query_coor:
            print(EG_blastx_query_range)
            raise ValueError(blastx_outfn)
        
        if low_evalue_query_coor['+']:
            low_evalue_query_coor['+'] = sorted(combine_indices(low_evalue_query_coor['+']))
        if low_evalue_query_coor['-']:
            low_evalue_query_coor['-'] = sorted(combine_indices(low_evalue_query_coor['-']))
    
        if chromosome not in EG_blastx_query_range:
            EG_blastx_query_range[chromosome] = {}
            
        EG_blastx_query_range[chromosome][EG] = (min_query_coor, max_query_coor)
    
        if chromosome not in blastx_hints:
            blastx_hints[chromosome] = {}
        blastx_hints[chromosome][EG] = []
    
        source = 'blastx'
        feature = 'CDSpart'
        score = '0'
        frame = '.' # . if unknown or irreelevant
        attribute = 'pri=4;src=P'       
        for strand in low_evalue_query_coor:
            if not low_evalue_query_coor[strand]:
                continue
            for coor in low_evalue_query_coor[strand]:
                start = coor[0]
                stop = coor[1]
                blastx_hints[chromosome][EG].append([chromosome, source+'@'+EG, feature, str(start), str(stop), score, strand, frame, attribute])

    return {'blastx_hints':blastx_hints, 'EG_blastx_query_range':EG_blastx_query_range}