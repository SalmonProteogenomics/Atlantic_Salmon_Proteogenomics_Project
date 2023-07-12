# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import argparse
import os
import time
from operator import itemgetter
import numpy as np
import pickle
from parse_reference_files import read_refseq_gff, interval_overlap,\
    find_eventpep, combine_indices, readFasta, ncbi_entrez_genpept,\
    combine_int_into_interval,\
    read_ref_prot_genomic, read_ref_prot_cds, read_ref_prot_gff,\
    ncbi_chrname_mapping
 
base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23',\
                    help="")
parser.add_argument("-email", "--email", dest = "email",\
                    default = '',\
                    help="email to use for entrez")    
args = parser.parse_args()

print(args.dataset)

script_start_time = time.time()

max_hits_per_prot = 500

ncbi_chrname_map = ncbi_chrname_mapping()
chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}

entrez_email = str(args.email)
print(entrez_email)
if not entrez_email:
    raise ValueError('no entrez email provided')

##reference files
reference_dir = os.path.join(base_dir , 'data', 'reference')
#GenPept directory
GenPept_dir = os.path.join(reference_dir, 'GenPept')
if not os.path.isdir(GenPept_dir):
    os.mkdir(GenPept_dir)
            
#no identity cutoff, but add areas with midline + / ' ' continuous>8 to non align region
alignment_gap_nonalgn = 10
hsp_cutoff = 1e-10

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

#RefSeq Ensembl db search output
RefSeq_Ensembl_dir = os.path.join(data_dir,'2_RefSeq_Ensembl_Search')
workflow_output_dir = os.path.join(RefSeq_Ensembl_dir, 'workflow_output')
genomic_coor_output_dir = os.path.join(RefSeq_Ensembl_dir, 'genomic_coordinates')

#Augustus
augustus_dir = os.path.join(data_dir, '3_Augustus')
augustus_parsed_dir = os.path.join(augustus_dir, 'parsed_output')
augustus_ref_overlap = os.path.join(augustus_parsed_dir,'Augustus_reference_overlap.pickle')
#Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]
augustus_ref_overlap_entries = pickle.load(open(augustus_ref_overlap,'rb'))

'#Augustus_Predicted_Prot	Chromosome'
'ProtSeq1|chrssa21_EG_1	chrssa21'
'Augustus_pred_prot_summary.txt'

#blastp
blastp_dir = os.path.join(data_dir, '4_blastp')
blastp_out_dir = os.path.join(blastp_dir, 'output')
blastp_parsed_dir = os.path.join(blastp_dir, 'parsed_output')


all_blastp_hits_pickle = os.path.join(blastp_parsed_dir, 'all_blastp_hits.pickle')
all_accessions_pickle = os.path.join(blastp_parsed_dir, 'all_accessions.pickle')

augustus_parsed_ss = os.path.join(blastp_parsed_dir, 'Salmo_salar')
augustus_parsed_other = os.path.join(blastp_parsed_dir, 'Other')
for dirn in [augustus_parsed_ss, augustus_parsed_other]:
    if not os.path.isdir(dirn):
        os.makedirs(dirn)

augustus_predprot_info = {}
with open(os.path.join(augustus_parsed_dir, 'Augustus_pred_prot_summary.txt'),'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            #Augustus_Predicted_Prot	Chromosome	jbrowse_coor	Augustus_gene_id	Augustus_transcript_id	Augustus_transcript_cds
            prot = sp[colidx['#Augustus_Predicted_Prot']].split('|')[0] #ProtSeq1|chrssa21_EG_1
            Chromosome = sp[colidx['Chromosome']]
            jbrowse_coor = sp[colidx['jbrowse_coor']]
            Augustus_transcript_cds = sp[colidx['Augustus_transcript_cds']]
            augustus_predprot_info[prot] = Chromosome
    
        
#RefSeq
refseq_dir = os.path.join(reference_dir, 'RefSeq')
genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')
Protein_genomic_coordiates = os.path.join(genomic_coor_output_dir,\
                                'RefSeq_Ensembl_Protein_genomic_coordinates.txt')
Protein_CDS_coordiates = os.path.join(genomic_coor_output_dir,\
                                'RefSeq_Ensembl_Protein_CDS_coordinates.txt')

refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
refseq_gff_results = read_refseq_gff(refseq_gff)
refseq_gene_coor = refseq_gff_results['refseq_gene_coor']

refseq_prot_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.gpff')
Salmo_salar_ref_prot_gff = read_ref_prot_gff(refseq_prot_gff)

refseq_prot_fasta = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.faa')
refseq_proteins = readFasta(refseq_prot_fasta)

#Ensembl
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')
ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta)

refseq_ensembl_proteins = list(refseq_proteins.keys())
refseq_ensembl_proteins.extend(ensembl_proteins.keys())
 
#####################################################################

print('Reading Known protein genomic coordinates...')
known_prot_genomic_coor_fn = os.path.join(genomic_coor_output_dir, 'RefSeq_Ensembl_Protein_genomic_coordiates.pickle')
if not os.path.isfile(known_prot_genomic_coor_fn):
    known_prot_genomic_coor = read_ref_prot_genomic(Protein_genomic_coordiates)
    with open(known_prot_genomic_coor_fn, 'wb') as handle:
        pickle.dump(known_prot_genomic_coor , handle, protocol=pickle.HIGHEST_PROTOCOL)
known_prot_genomic_coor = pickle.load(open(known_prot_genomic_coor_fn,'rb'))        
print(len(known_prot_genomic_coor),  'reference proteins with genomic coordinates...')

print( 'Reading Known proteins with incomplete CDS info...')
known_prot_genomic_coor_incomplete_fn = os.path.join(genomic_coor_output_dir, 'RefSeq_Ensembl_Protein_CDS_coordiates.pickle')
if not os.path.isfile(known_prot_genomic_coor_incomplete_fn):
    known_prot_genomic_coor_incomplete = read_ref_prot_cds(Protein_CDS_coordiates)
    with open(known_prot_genomic_coor_incomplete_fn, 'wb') as handle:
        pickle.dump(known_prot_genomic_coor_incomplete , handle, protocol=pickle.HIGHEST_PROTOCOL)
known_prot_genomic_coor_incomplete = pickle.load(open(known_prot_genomic_coor_incomplete_fn,'rb'))
print(len(known_prot_genomic_coor_incomplete), 'reference proteins with INCOMPLETE genomic coordinates...')

all_accessions = {}    
all_blastp_hits = {}

def read_blastp_xml(outfn):
    all_hits = {}
    hit = ''
    with open(outfn,'r') as infile:
        for line in infile:
            if '<BlastOutput_query-def>' in line:
                tag = line.strip().split('<')[1].split('>')[0]
                if os.path.basename(outfn).split('.')[0] not in line.strip().split('<'+tag+'>')[1].split('</'+tag+'>')[0]:
                    raise ValueError(outfn, line)
                
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
    return all_hits

def parse_blastp(blastp_fn, predprotid, augustus_pp_sequence):    
    #read xml
    blastp_hits = read_blastp_xml(os.path.join(blastp_out_dir, blastp_fn))
    if not blastp_hits:
        return None     

    #determine which hits are salmon salar hits, and keep top 5 lowest evalue non salmon salar hits
    Salmon_hits = {}
    Non_Salmon_hits = {}
    Salmon = False
    count_hits = 0
    for hit in blastp_hits.values():
        Salmon = False
        skip_pred_prot = False
        
        alignments = [(float(hsp['Hsp_evalue']), 100*float(hsp['Hsp_identity'])/float(hsp['Hsp_align-len']))\
                      for hsp in hit['Hsp'] if float(hsp['Hsp_evalue'])<=hsp_cutoff]
        if not alignments:
            continue            
        
        Hit_num = hit['Hit_num']
        ########## Accessions associated with hit ###########
        #<Hit_id>ref|XP_036820350.1|</Hit_id>
        Hit_id = hit['Hit_id']
        accession = ''
        if not Hit_id.split('|')[2]:
            accession = Hit_id.split('|')[1]
        else:
            accession = Hit_id[0:]
        #<Hit_accession>XP_036820350</Hit_accession>
        #Hit_accession = hit['Hit_accession']
            
        if accession not in all_accessions: 
            entrez_result = ncbi_entrez_genpept(accession, GenPept_dir, entrez_email)
            if entrez_result == 'Could not retrieve':
                #try again
                entrez_result = ncbi_entrez_genpept(accession, GenPept_dir, entrez_email)
                if entrez_result == 'Could not retrieve':
                    print('Could not retrieve ', accession, entrez_result)
                    continue
            if entrez_result == 'Unexpected source line':
                continue
            if entrez_result == 'skip':
                continue
            all_accessions[accession] = entrez_result

        ss_accessions = {}
        other_hits = {}            
        organism = all_accessions[accession]['ORGANISM']
        chromosome = all_accessions[accession]['chromosome']
        if organism == 'Salmo salar':
            Salmon = True
            #some accessions don't have genomic_coor due to incomplete annotations
            if accession in refseq_ensembl_proteins:
                ss_accessions[accession] = 'Ref'
            else:
                ss_accessions[accession] = 'non-reference'
            
            if 'chr'+chromosome in chrname_ncbi_map:
                chromosome = chrname_ncbi_map['chr'+chromosome]
                augustus_chromosome = augustus_predprot_info[predprotid]
                #skip if predicted augustus protein is the same as 
                if augustus_pp_sequence == all_accessions[accession]['SEQUENCE']:
                    if augustus_chromosome==chromosome:
                        print(accession, chromosome)
                        skip_pred_prot = True
        else:
            other_hits[accession] = organism
        
        if skip_pred_prot:
            continue
        
        if Salmon == True:
            Salmon_hits[Hit_num] = {'ss':ss_accessions,\
                                    'other':[x+'['+other_hits[x]+']' for x in other_hits]}
        else:
            if not other_hits:
                continue
            lowest_evalue = min([x[0] for x in alignments])
            aligned_identities = [x[1] for x in alignments]
            Non_Salmon_hits[Hit_num] = {'other':other_hits,\
                                        'evalue':lowest_evalue,\
                                        'aligned_identities':aligned_identities}

        count_hits+=1
        if count_hits>max_hits_per_prot:
            break

    if Non_Salmon_hits:
        lowest_evalue_non_ss = [(Non_Salmon_hits[x]['evalue'],\
                                 np.mean(Non_Salmon_hits[x]['aligned_identities']),\
                                 len(Non_Salmon_hits[x]['aligned_identities']),\
                                 x)\
                                 for x in Non_Salmon_hits]
        lowest_evalue_non_ss.sort(key=itemgetter(1),reverse=True)
        lowest_evalue_non_ss.sort(key=itemgetter(2),reverse=True)
        lowest_evalue_non_ss = [x[3] for x in lowest_evalue_non_ss[:max_hits_per_prot]]
        
        non_ss_hits = list(Non_Salmon_hits.keys())
        for Hit_num in non_ss_hits:
            if Hit_num not in lowest_evalue_non_ss:
                del Non_Salmon_hits[Hit_num]
    
    all_hitnums = set(Salmon_hits.keys())
    all_hitnums.update(Non_Salmon_hits.keys())
    blastp_hits = {x:blastp_hits[x] for x in blastp_hits if x in all_hitnums}

    return {'Salmon_hits':Salmon_hits,
            'Non_Salmon_hits':Non_Salmon_hits,
            'blastp_hits':blastp_hits} 

print('parsing blastp output...')
if not os.path.isfile(all_blastp_hits_pickle):
    for outfn in os.listdir(blastp_out_dir):
        if '.blastp.log' in outfn:
            continue
        predprotid = outfn.split('.blastp.txt')[0]
        augustus_pp_sequence = augustus_ref_overlap_entries[predprotid]['aug_prot_seq']
        print(predprotid)
        results = parse_blastp(outfn, predprotid, augustus_pp_sequence)
        
        if results:
            all_blastp_hits[predprotid] = results
        else:
            print('no blastp results: ',predprotid)

    with open(all_accessions_pickle, 'wb') as handle:
        pickle.dump(all_accessions, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open(all_blastp_hits_pickle, 'wb') as handle:
        pickle.dump(all_blastp_hits, handle, protocol=pickle.HIGHEST_PROTOCOL)
 
all_blastp_hits = pickle.load(open(all_blastp_hits_pickle,'rb'))
all_accessions = pickle.load(open(all_accessions_pickle,'rb'))


logfn = open(os.path.join(blastp_dir,'Parse_blastp_XML_output_log.txt'),'w')

for outfn in os.listdir(blastp_out_dir):
    if '.blastp.log' in outfn:
        continue
    
    predprotid = outfn.split('.blastp.txt')[0]
    
    if predprotid not in all_blastp_hits:
        continue
    
    ss_hits_fn = os.path.join(augustus_parsed_ss, predprotid+'_Salmon_hits_info.pickle')
    ss_ns_hits_fn = os.path.join(augustus_parsed_ss, predprotid+'_Non_Salmon_hits_info.pickle')
    ns_hits_fn = os.path.join(augustus_parsed_other, predprotid+'_Non_Salmon_hits_info.pickle')

    augustus_pp_sequence = augustus_ref_overlap_entries[predprotid]['aug_prot_seq']   
    augustus_pp_length = len(augustus_pp_sequence)
    
    prediction = list(augustus_ref_overlap_entries[predprotid]['predictions'].keys())[0]
    logfn.write('##'+prediction+'\n')
    
    EG, aug_geneid, aug_transcriptid = prediction.split('.')
    logfn.write('#'+'\t'.join([predprotid, EG, aug_geneid, aug_transcriptid])+'\n')
    
    chromosome = EG.split('_')[0]
    
    Entry = augustus_ref_overlap_entries[predprotid]['predictions'][prediction]
    aug_gene_s_e = Entry['aug_gene_s_e'] #1-index
    aug_cds_s_e = Entry['aug_cds_s_e'] #1-index
    aug_strand = Entry['aug_strand']
    
    #{'AAQKPQPVFPSPEQMLFDHEFTK|chrssa08|+|16570207-16570275', 'ATGVPADILTETINTVSEVIR|chrssa08|+|16564723-16564755,16564896-16564925'}
    known_peps_genomic_coordinates = Entry['known_peps_genomic_coordinates']
    #TVWDWELMNDIK[316]
    Known_Peptides_Prot_Coor = Entry['known_peps_protein_coordinates']

    #{'LMDLLADSR|chrssa08|+|5424|16563576-16563602'}
    event_peps_genomic_coordinates = Entry['event_peps_genomic_coordinates']
    predicted_prot_event_peptides = set([x.split('|')[0] for x in event_peps_genomic_coordinates])
        
    Aug_Known_overlap = {}
    overlap_ref_genes = {}
    #CDS 5UTR 3UTR genomic coordinate overlap between reference annotation and augustus predicted transcript
    refseq_overlap = Entry['refseq_overlap']
    for x in refseq_overlap:
        rg = x.replace('gene-','')
        overlap_ref_genes[rg] = {'start':int(refseq_gene_coor[x]['start']),'end':int(refseq_gene_coor[x]['end'])}
        Aug_Known_overlap[rg] = {}
        for z in refseq_overlap[x]:
            rt = z.replace('rna-','')
            Aug_Known_overlap[rg][rt] = refseq_overlap[x][z]
            for t in ['CDS','5UTR','3UTR']:
                if t not in Aug_Known_overlap[rg][rt]:
                    Aug_Known_overlap[rg][rt][t] = []
    
    related_genes = Entry['related_genes_RefSeq']
    overlap_ref_protein_genes = {}
    overlap_ref_mrna_prot = {}
    for x in related_genes:
        rg,rt,rp = x.split('|')
        rg = rg.replace('gene-','')
        rt = rt.replace('rna-','')
        if rp in ['5UTR','3UTR']:
            continue
        if rp in overlap_ref_protein_genes:
            raise ValueError('')
        overlap_ref_protein_genes[rp] = {'gene':rg,'mrna':rt}
        overlap_ref_mrna_prot[rt] = rp

    #keys are protein aa index (integer):{'aa': 'M', 'aa_coor': [16559841, 16559842, 16559843]}
    aug_prot_aa_genomic_coor = Entry['aug_prot_aa_genomic_coor'] #1 index

    blastp_results = all_blastp_hits[predprotid]
    Non_Salmon_hits = blastp_results['Non_Salmon_hits']
    Salmon_hits = blastp_results['Salmon_hits']
    blastp_hits = blastp_results['blastp_hits']

   
    ###################  Non_Salmon_hits ###################
    #           #report aligned regions
    #           #event peptide in predicted protein's aligned regions?
    #           #alignment qseq,midline,hseq
    Non_Salmon_hits_info = {}
    if Non_Salmon_hits:
        #top5 non-salmo salar hits
        for Hit_num in Non_Salmon_hits:
            main_accession = list(Non_Salmon_hits[Hit_num]['other'].keys())[0]
            organism = Non_Salmon_hits[Hit_num]['other'][main_accession]
            hit = blastp_hits[Hit_num]
            ########## HSPs ###########
            alignment_text = []
            query_hsp_start_end_prot_coor = []
            query_aligned_full_genomic_coor = [] # aligned + 100p identity
            query_aligned_other_genomic_coor = [] # aligned + <100p identity 
            evalue = []
            for hsp in hit['Hsp']:
                if float(hsp['Hsp_evalue'])<=hsp_cutoff:
                    evalue.append(float(hsp['Hsp_evalue']))      
                    hsp_identityp = round(100*float(hsp['Hsp_identity'])/float(hsp['Hsp_align-len']),1)
                    #keep alignment if pass evalue cutoff
                    query_pre_txt = 'Query'+' '+hsp['Hsp_query-from']+' '
                    hit_pre_txt = 'Sbjct'+' '+hsp['Hsp_hit-from']+' '
                    midline_pre_txt = ''
                    if len(query_pre_txt)>len(hit_pre_txt):
                        hit_pre_txt += ''.join(['\s' for x in range(len(query_pre_txt)-len(hit_pre_txt))])
                        midline_pre_txt += ''.join(['\s' for x in range(len(query_pre_txt))])
                    elif len(query_pre_txt)<len(hit_pre_txt):
                        query_pre_txt += ''.join(['\s' for x in range(len(hit_pre_txt)-len(query_pre_txt))])
                        midline_pre_txt += ''.join(['\s' for x in range(len(hit_pre_txt))])
                    else:
                        midline_pre_txt = ''.join(['\s' for x in range(len(hit_pre_txt))])          
                    query_txt = query_pre_txt+hsp['Hsp_qseq']+' '+hsp['Hsp_query-to']
                    hit_txt = hit_pre_txt+hsp['Hsp_hseq']+' '+hsp['Hsp_hit-to']
                    midline_txt = midline_pre_txt+hsp['Hsp_midline']                            
                    alignment_text.append(['#identity:'+str(hsp_identityp), query_txt,midline_txt,hit_txt])

                    #find large gaps                         
                    non_match_pos = []
                    for ci,c in enumerate(hsp['Hsp_midline']):
                        if not c.isalpha():
                            non_match_pos.append(ci)
                    if non_match_pos:
                        non_match_pos = combine_int_into_interval(non_match_pos) #0 index
                        
                    add_nonalign = []
                    for c in non_match_pos:
                        #if gap region is >=alignment_gap_nonalgn
                        if len(hsp['Hsp_midline'][c[0]:c[1]+1])>=alignment_gap_nonalgn:
                            #print hsp['Hsp_midline'][c[0]:c[1]+1]
                            add_nonalign.append(c)
                    match_pos = []
                    for i in range(int(hsp['Hsp_query-from']),int(hsp['Hsp_query-to'])+1):
                        if not [c for c in add_nonalign if c[0]<=i-int(hsp['Hsp_query-from'])<=c[1]]:
                            match_pos.append(i)
                    if match_pos:
                        match_pos = combine_int_into_interval(match_pos) #1index protein coor
                    
                    for mr in match_pos:
                        midline = hsp['Hsp_midline'][mr[0]-1:mr[1]]
                        #if match region < alignment_gap_nonalgn, not a match region?
                        if len(midline)<alignment_gap_nonalgn:
                            continue
                        #identity: length after removing /s and + in Hsp_midline   
                        hsp_identity = len(midline.replace(' ','').replace('-',''))
                        hsp_identityp = round(100*float(hsp_identity)/float(len(midline)),1) 
                        
                        query_hsp_start_end_prot_coor.append(mr)
                        
                        query_genomic_coordinates = set()
                        for i in range(mr[0]-1,mr[1]): #change to 0 index
                            query_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor']) #1 index
                        query_genomic_coordinates = (min(query_genomic_coordinates),max(query_genomic_coordinates))
                        
                        if hsp_identityp==100.0:
                            query_aligned_full_genomic_coor.append(query_genomic_coordinates)
                        else:
                            query_aligned_other_genomic_coor.append(query_genomic_coordinates)
            
            if not query_hsp_start_end_prot_coor:
                continue
            
            query_aligned_regions = combine_indices(query_hsp_start_end_prot_coor) # 1-index
            #event peptide in aligned region of pred prot?
            query_aligned_region_eventpeps = set()
            for s,e in query_aligned_regions:
                found = find_eventpep(augustus_pp_sequence[s-1:e], predicted_prot_event_peptides)
                #found[ep] = index 
                if found:
                    for event_peptide in found:
                        prot_coor = found[event_peptide]
                        prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                        for pc in prot_coor:
                            query_aligned_region_eventpeps.add(((s,e), pc, event_peptide)) # 1index
                
            #when aligned to main_accession
            if query_aligned_full_genomic_coor:
                query_aligned_full_genomic_coor = combine_indices(query_aligned_full_genomic_coor)
            if query_aligned_other_genomic_coor:
                query_aligned_other_genomic_coor = combine_indices(query_aligned_other_genomic_coor)
            #1 index
            Non_Salmon_hits_info[main_accession] = {'query_aligned_region_eventpeps':query_aligned_region_eventpeps, #1index
                                                    'query_aligned_full_genomic_coor':query_aligned_full_genomic_coor, #1index
                                                    'query_aligned_other_genomic_coor':query_aligned_other_genomic_coor, #1index
                                                    'alignment_text':alignment_text,
                                                    'organism':organism,
                                                    'evalue':min(evalue)}
            
    ################### Salmon_hits ###################
    #           #min(Hsp_query-from) < or > min(Hsp_hit-from)
    #           #max(Hsp_query-to) < or > max(Hsp_hit-to)
    #           #
    #           #known peptides in pred prot
    #           #event peptides in pred prot
    #           #100p aligned
    #           #   #augustus CDS == reference CDS difference?
    #           #non-aligned
    #           #   #event peptide in non-aligned region?
    #           #   #reference UTR in non-aligned region?
    #           #   #augustus CDS != reference CDS difference? augustus CDS shorter? augustus CDS longer?
    #           #<100p aligned
    #           #   #event peptide in non-aligned region?
    #           #   #augustus CDS and reference CDS difference?
    #           #   #
    
    Salmon_hits_info = {}        
    if Salmon_hits:
        for Hit_num in Salmon_hits:
            ###########################################################################
            ss_accessions = Salmon_hits[Hit_num]['ss']
            ref_accessions = [a for a in ss_accessions if ss_accessions[a]=='Ref']  
            ref_accessions_unplaced = [a for a in ref_accessions if a not in Salmo_salar_ref_prot_gff]
            ref_accessions_placed = [a for a in ref_accessions if a in Salmo_salar_ref_prot_gff]
            ref_accessions_same_chrom = [a for a in ref_accessions_placed if Salmo_salar_ref_prot_gff[a]['chrom']==chromosome]
            ref_accessions_diff_chrom = [a for a in ref_accessions_placed if Salmo_salar_ref_prot_gff[a]['chrom']!=chromosome]
            
            non_ref_accessions = [a for a in ss_accessions if ss_accessions[a]=='non-reference']
            non_ref_accessions_unplaced = [a for a in non_ref_accessions if a not in Salmo_salar_ref_prot_gff]
            non_ref_accessions_placed = [a for a in non_ref_accessions if a in Salmo_salar_ref_prot_gff]
            non_ref_accessions_same_chrom = [a for a in non_ref_accessions_placed if Salmo_salar_ref_prot_gff[a]['chrom']==chromosome]
            non_ref_accessions_diff_chrom = [a for a in non_ref_accessions_placed if Salmo_salar_ref_prot_gff[a]['chrom']!=chromosome]
            
            ref_ss_accessions = {}
            for a in ref_accessions_same_chrom:
                if a in known_prot_genomic_coor:
                    if 'g_coor' not in ref_ss_accessions:
                        ref_ss_accessions['g_coor'] = set()
                    ref_ss_accessions['g_coor'].add(a)
                elif a in known_prot_genomic_coor_incomplete:
                    if 'cds_coor' not in ref_ss_accessions:
                        ref_ss_accessions['cds_coor'] = set()
                    ref_ss_accessions['cds_coor'].add(a)  
            
            #accessions whose reference gene has overlap with augustus predicted gene
            main_accessions = {a:overlap_ref_protein_genes[a] for a in ss_accessions if a in overlap_ref_protein_genes}
            main_accessions_genes = set([overlap_ref_protein_genes[a]['gene'] for a in ss_accessions if a in overlap_ref_protein_genes])
            
#            if main_accessions:
#                print chromosome, prot_seq_i
#                print main_accessions

            
            ###########################################################################
            hit = blastp_hits[Hit_num]
            ########## HSPs ###########
            alignment_text = []
            query_hsp_start_end_prot_coor = []
            hit_hsp_start_end_prot_coor = []

            query_aligned_full_genomic_coor = [] # aligned + 100p identity
            query_aligned_full_prot_coor = []
            hit_aligned_full_prot_coor = []
            query_aligned_other_genomic_coor = [] # aligned + <100p identity
            query_aligned_other_prot_coor = []
            hit_aligned_other_prot_coor = []
            evalue = []                                
            for hsp in hit['Hsp']:
                if float(hsp['Hsp_evalue'])<=hsp_cutoff:
                    evalue.append(float(hsp['Hsp_evalue']))      
                    hsp_identityp = round(100*float(hsp['Hsp_identity'])/float(hsp['Hsp_align-len']),1)
                    #keep alignment if pass evalue cutoff
                    query_pre_txt = 'Query'+' '+hsp['Hsp_query-from']+' '
                    hit_pre_txt = 'Sbjct'+' '+hsp['Hsp_hit-from']+' '
                    midline_pre_txt = ''
                    if len(query_pre_txt)>len(hit_pre_txt):
                        hit_pre_txt += ''.join(['\s' for x in range(len(query_pre_txt)-len(hit_pre_txt))])
                        midline_pre_txt += ''.join(['\s' for x in range(len(query_pre_txt))])
                    elif len(query_pre_txt)<len(hit_pre_txt):
                        query_pre_txt += ''.join(['\s' for x in range(len(hit_pre_txt)-len(query_pre_txt))])
                        midline_pre_txt += ''.join(['\s' for x in range(len(hit_pre_txt))])
                    else:
                        midline_pre_txt = ''.join(['\s' for x in range(len(hit_pre_txt))])          
                    query_txt = query_pre_txt+hsp['Hsp_qseq']+' '+hsp['Hsp_query-to']
                    hit_txt = hit_pre_txt+hsp['Hsp_hseq']+' '+hsp['Hsp_hit-to']
                    midline_txt = midline_pre_txt+hsp['Hsp_midline']                            
                    alignment_text.append(['#identity:'+str(hsp_identityp), query_txt,midline_txt,hit_txt])

                    #query
                    hsp['Hsp_query-from'] = int(hsp['Hsp_query-from']) #1 index
                    hsp['Hsp_query-to'] = int(hsp['Hsp_query-to'])
                    
                    #hit
                    hsp['Hsp_hit-from'] = int(hsp['Hsp_hit-from'])
                    hsp['Hsp_hit-to'] = int(hsp['Hsp_hit-to'])

                    #find large gaps                         
                    non_match_pos = []
                    for ci,c in enumerate(hsp['Hsp_midline']):
                        if not c.isalpha():
                            non_match_pos.append(ci)
                    if non_match_pos:
                        non_match_pos = combine_int_into_interval(non_match_pos) #0index
                        
                    add_nonalign = []
                    for c in non_match_pos:
                        if len(hsp['Hsp_midline'][c[0]:c[1]+1])>=alignment_gap_nonalgn:
                            #print hsp['Hsp_midline'][c[0]:c[1]+1]
                            add_nonalign.append(c)
                    match_pos = []
                    for i in range(int(hsp['Hsp_query-from']),int(hsp['Hsp_query-to'])+1):
                        if not [c for c in add_nonalign if c[0]<=i-int(hsp['Hsp_query-from'])<=c[1]]:
                            match_pos.append(i)
                    if match_pos:
                        match_pos = combine_int_into_interval(match_pos) #1index protein coor
                    
                    for mr in match_pos:
                        midline = hsp['Hsp_midline'][mr[0]-1:mr[1]]
                        if len(midline)<alignment_gap_nonalgn:
                            continue
                        #identity: length after removing /s and + in Hsp_midline   
                        hsp_identity = len(midline.replace(' ','').replace('-',''))
                        hsp_identityp = round(100*float(hsp_identity)/float(len(midline)),1) 
                        
                        query_hsp_start_end_prot_coor.append(mr)
                        h_mr_s = hsp['Hsp_hit-from']+(mr[0]-hsp['Hsp_query-from'])
                        h_mr_e = h_mr_s+(mr[1]-mr[0])
                        hit_hsp_start_end_prot_coor.append((h_mr_s, h_mr_e))
                        
                        query_genomic_coordinates = set()
                        for i in range(mr[0]-1,mr[1]):
                            query_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor'])
                        query_genomic_coordinates = (min(query_genomic_coordinates),max(query_genomic_coordinates))
                        
                        if hsp_identityp==100.0:
                            query_aligned_full_genomic_coor.append(query_genomic_coordinates)
                            query_aligned_full_prot_coor.append(mr) # 1 index
                            hit_aligned_full_prot_coor.append((h_mr_s, h_mr_e))
                        else:
                            query_aligned_other_genomic_coor.append(query_genomic_coordinates)
                            query_aligned_other_prot_coor.append(mr)
                            hit_aligned_other_prot_coor.append((h_mr_s, h_mr_e))
                    
                    #don't have genomic info for salmo salar accessions? treat as non salmo salar
                    if ref_ss_accessions:
                        query_genomic_coordinates = set()
                        for i in range(hsp['Hsp_query-from']-1,hsp['Hsp_query-to']):
                            query_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor'])
                        query_genomic_coordinates = (min(query_genomic_coordinates),max(query_genomic_coordinates))

                        if aug_strand == '-':
                            query_genomic_coordinates = query_genomic_coordinates[::-1] 
                            
            if not query_hsp_start_end_prot_coor:
                continue
            #1-index
            query_aligned_regions = combine_indices(query_hsp_start_end_prot_coor)
            if query_aligned_full_genomic_coor:
                query_aligned_full_genomic_coor = combine_indices(query_aligned_full_genomic_coor)
            if query_aligned_other_genomic_coor:
                query_aligned_other_genomic_coor = combine_indices(query_aligned_other_genomic_coor) 
            
            #?#
            if not ref_ss_accessions:
                #event peptide in aligned region of pred prot?
                query_aligned_region_eventpeps = set()
                for s,e in query_aligned_regions:
                    found = find_eventpep(augustus_pp_sequence[s-1:e], predicted_prot_event_peptides)
                    if found:
                        for event_peptide in found:
                            prot_coor = found[event_peptide]
                            prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                            for pc in prot_coor:
                                query_aligned_region_eventpeps.add(((s,e), pc, event_peptide))
                    
                #when aligned to main_accession
                Salmon_hits_info[Hit_num] = {'main_accessions':main_accessions,
                                             'ref_accessions_unplaced':ref_accessions_unplaced,
                                             'non_ref_accessions_unplaced':non_ref_accessions_unplaced,
                                             'ref_accessions_same_chrom':ref_accessions_same_chrom,
                                             'ref_accessions_diff_chrom':ref_accessions_diff_chrom,
                                             'non_ref_accessions_same_chrom':non_ref_accessions_same_chrom,
                                             'non_ref_accessions_diff_chrom':non_ref_accessions_diff_chrom,
                                             'query_aligned_region_eventpeps':query_aligned_region_eventpeps, #1 index
                                             'query_aligned_full_genomic_coor':query_aligned_full_genomic_coor, #1 index
                                             'query_aligned_other_genomic_coor':query_aligned_other_genomic_coor, #1 index
                                             'alignment_text':alignment_text,
                                             'identical_accessions':Salmon_hits[Hit_num],
                                             'evalue':min(evalue)}
                continue

            Salmon_hits_info[Hit_num] = {'ref_ss_accessions':ref_ss_accessions,
                                        'main_accessions':main_accessions,
                                        'ref_accessions_unplaced':ref_accessions_unplaced,
                                        'non_ref_accessions_unplaced':non_ref_accessions_unplaced,
                                        'ref_accessions_same_chrom':ref_accessions_same_chrom,
                                        'ref_accessions_diff_chrom':ref_accessions_diff_chrom,
                                        'non_ref_accessions_same_chrom':non_ref_accessions_same_chrom,
                                        'non_ref_accessions_diff_chrom':non_ref_accessions_diff_chrom,
                                        'non-aligned':{'query_events':[], 'overlap_ref_geneids':[]},
                                        '100p-aligned':{'query_events':[],'hit_events':[], 'overlap_ref_geneids':[]},
                                        '<100p-aligned':{'query_events':[],'hit_events':[], 'overlap_ref_geneids':[]},
                                        'alignment_text':alignment_text,
                                        'query_aligned_full_genomic_coor':query_aligned_full_genomic_coor, #1 index
                                        'query_aligned_other_genomic_coor':query_aligned_other_genomic_coor, #1 index
                                        'identical_accessions':Salmon_hits[Hit_num],
                                        'evalue':min(evalue)}

            ################ NON-ALIGNED regions ################
            # a. event peptide in predicted protein? 
            # b1. ref gene in region? UTR regions?
            # b2. if ref gene, overlap with augustus CDS? augustus CDS shorter? augustus CDS longer?
            #Determine query nonaligned regions given aligned regions   
            query_nonaligned_p = []                
            for i in range(1,augustus_pp_length+1):
                if not [c for c in query_aligned_regions if c[0]<=i<=c[1]]:
                    query_nonaligned_p.append(i)
            if query_nonaligned_p:
                query_nonaligned_p = combine_int_into_interval(query_nonaligned_p) #1-index
            
            #make 1 index
            query_nonaligned_regions = {}
            for s,e in query_nonaligned_p:
                if s==e:
                    continue
                query_nonaligned_regions[s] = {e:{}}
            
            nonaligned_overlapping_genes = set()
            #protein coordinates, 1 index
            for s in query_nonaligned_regions:
                for e in query_nonaligned_regions[s]:
                    #a. events in nonaligned region of predicted protein seq
                    found = find_eventpep(augustus_pp_sequence[s-1:e], predicted_prot_event_peptides)
                    if found:
                        query_nonaligned_regions[s][e]['events'] = set()
                        for event_peptide in found:
                            prot_coor = found[event_peptide]
                            prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                            for pc in prot_coor:
                                query_nonaligned_regions[s][e]['events'].add((pc, event_peptide))
                
                    #b1. check area for reference gene
                    # hit proteins should belong to reference gene
                    non_aligned_genomic_coordinates = set()
                    #0 index
                    for i in range(s-1,e):
                        non_aligned_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor'])
                    genomic_s, genomic_e = (min(non_aligned_genomic_coordinates), max(non_aligned_genomic_coordinates)) #1 index
                    
                    nonaligned_overlapping_genes = set()
                    #entry is reference gene that overlaps (genomic coordinates) with predicted augustus gene but non-aligned??
                    for ref_gene in overlap_ref_genes:
                        if ref_gene not in main_accessions_genes:
                            continue
                        #overlap gene start end (genomic)
                        b = (overlap_ref_genes[ref_gene]['start'], overlap_ref_genes[ref_gene]['end'])
                        if genomic_s <= b[1] and b[0] <= genomic_e:
                            nonaligned_overlapping_genes.add(ref_gene)
                    
                    #b2. if ref gene, overlap with augustus CDS or UTR?
                    if nonaligned_overlapping_genes:
                        for refseq_geneid in Aug_Known_overlap:
                            if refseq_geneid in nonaligned_overlapping_genes:
                                ref_gene_entry = Aug_Known_overlap[refseq_geneid]
                                #['3UTR', '5UTR', 'CDS']
                                #{'augustus': ('cds1', '2393082', '2393210'),
                                # 'protein_id': 'ENSSSAP00000054587',
                                # 'ref': ('CDS:ENSSSAP00000054587', 2393082, 2393233)}
                                for ref_gene_transcript in ref_gene_entry:
                                    ref_gene_CDS = ref_gene_entry[ref_gene_transcript]['CDS']
                                    ref_gene_5UTR = ref_gene_entry[ref_gene_transcript]['5UTR']
                                    ref_gene_3UTR = ref_gene_entry[ref_gene_transcript]['3UTR']
                                    if ref_gene_CDS or ref_gene_5UTR or ref_gene_3UTR:
                                        logfn.write(refseq_geneid+'\t'+ref_gene_transcript+'\n')
                                        if 'overlap_ref_geneids' not in query_nonaligned_regions[s][e]:
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'] = {}
                                        if refseq_geneid not in query_nonaligned_regions[s][e]['overlap_ref_geneids']:
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid] = {}
                                        if ref_gene_transcript not in query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid]:
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript] = {'CDS':[],\
                                                  '5UTR':[],'3UTR':[]}
                                        if ref_gene_CDS:
                                            logfn.write('CDS: '+str(ref_gene_CDS)+'\n')
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['CDS'].extend(ref_gene_CDS)
                                        elif ref_gene_5UTR:
                                            logfn.write('5UTR: '+str(ref_gene_5UTR)+'\n')
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['5UTR'].extend(ref_gene_5UTR)
                                        elif ref_gene_3UTR:
                                            logfn.write('3UTR: '+str(ref_gene_3UTR)+'\n')
                                            query_nonaligned_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['3UTR'].extend(ref_gene_3UTR)
            #1 index
            for s in query_nonaligned_regions:
                for e in query_nonaligned_regions[s]:
                    if 'events' in query_nonaligned_regions[s][e]:
                        Salmon_hits_info[Hit_num]['non-aligned']['query_events'].append((s,e,query_nonaligned_regions[s][e]['events']))
                    if 'overlap_ref_geneids' in query_nonaligned_regions[s][e]:
                        Salmon_hits_info[Hit_num]['non-aligned']['overlap_ref_geneids'].append((s,e,query_nonaligned_regions[s][e]['overlap_ref_geneids']))

            ################ 100p ALIGNED regions ################
            #a. event peptide in predicted protein but not ref protein? event peptide in aligned hit sequence? 
            #b. augustus pred CDS and ref CDS overlap?
            # 1 index
            if query_aligned_full_prot_coor:
                query_aligned_full_prot_coor = combine_indices(query_aligned_full_prot_coor)
                hit_aligned_full_prot_coor = combine_indices(hit_aligned_full_prot_coor)

                #a.
                hits_aligned_full_regions_ep = {}
                for t in ref_ss_accessions:
                    for accession in ref_ss_accessions[t]: 
                        found = find_eventpep(all_accessions[accession]['SEQUENCE'], predicted_prot_event_peptides)
                        if found:
                            for event_peptide in found:
                                prot_coor = found[event_peptide]
                                prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                                for pc in prot_coor:
                                    ovp = interval_overlap(pc, hit_aligned_full_prot_coor)
                                    if ovp:
                                        if accession not in hits_aligned_full_regions_ep:
                                            hits_aligned_full_regions_ep[accession]= set()
                                        hits_aligned_full_regions_ep[accession].add((pc, event_peptide))
                
                if hits_aligned_full_regions_ep:
                    Salmon_hits_info[Hit_num]['100p-aligned']['hit_events'] = hits_aligned_full_regions_ep

                query_aligned_full_regions = {}
                #1 index
                for s,e in query_aligned_full_prot_coor:
                    found = find_eventpep(augustus_pp_sequence[s-1:e], predicted_prot_event_peptides)
                    if found:
                        if s not in query_aligned_full_regions:
                            query_aligned_full_regions[s] = {e:{}}
                        query_aligned_full_regions[s][e]['events'] = set()
                        for event_peptide in found:
                            prot_coor = found[event_peptide]
                            prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                            for pc in prot_coor:
                                query_aligned_full_regions[s][e]['events'].add((pc, event_peptide))

                    #b1. check area for reference gene     
                    aligned_full_genomic_coordinates = set()
                    for i in range(s-1,e):
                        aligned_full_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor'])
                    genomic_s, genomic_e = (min(aligned_full_genomic_coordinates), max(aligned_full_genomic_coordinates))
                    
                    aligned_full_overlapping_genes = set()
                    #entry is reference gene that overlaps with predicted augustus gene
                    for ref_gene in overlap_ref_genes:
                        if ref_gene not in main_accessions_genes:
                            continue
                        #overlap gene start end (genomic)
                        b = (overlap_ref_genes[ref_gene]['start'], overlap_ref_genes[ref_gene]['end'])
                        if genomic_s <= b[1] and b[0] <= genomic_e:
                            aligned_full_overlapping_genes.add(ref_gene)
                    
                    #b2. if ref gene, overlap with augustus CDS or UTR?
                    if aligned_full_overlapping_genes:
                        for refseq_geneid in Aug_Known_overlap:
                            #
                            if refseq_geneid in nonaligned_overlapping_genes:
                                ref_gene_entry = Aug_Known_overlap[refseq_geneid]
                                #['3UTR', '5UTR', 'CDS']
                                #{'augustus': ('cds1', '2393082', '2393210'),
                                # 'protein_id': 'ENSSSAP00000054587',
                                # 'ref': ('CDS:ENSSSAP00000054587', 2393082, 2393233)}
                                for ref_gene_transcript in ref_gene_entry:
                                    ref_gene_CDS = ref_gene_entry[ref_gene_transcript]['CDS']
                                    ref_gene_5UTR = ref_gene_entry[ref_gene_transcript]['5UTR']
                                    ref_gene_3UTR = ref_gene_entry[ref_gene_transcript]['3UTR']
                                    #if not empty, list of overlapping features
                                        #[{'augustus': ('cds0', 44633435, 44633437),
                                        #  'ref': ('utr123579', 44632626, 44633448)},
                                        # {'augustus': ('cds1', 44633435, 44633448),
                                        #  'ref': ('utr123579', 44632626, 44633448)}]
                                    if ref_gene_CDS or ref_gene_5UTR or ref_gene_3UTR:
                                        logfn.write(refseq_geneid+'\t'+ref_gene_transcript+'\n')
                                        if s not in query_aligned_full_regions:
                                            query_aligned_full_regions[s] = {e:{}}    
                                        if 'overlap_ref_geneids' not in query_aligned_full_regions[s][e]:
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'] = {}
                                        if refseq_geneid not in query_aligned_full_regions[s][e]['overlap_ref_geneids']:
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid] = {}
                                        if ref_gene_transcript not in query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid]:
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript] = {'CDS':[],\
                                                  '5UTR':[],'3UTR':[]}
                                        if ref_gene_CDS:
                                            logfn.write('CDS: '+str(ref_gene_CDS)+'\n')
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['CDS'].extend(ref_gene_CDS)
                                        elif ref_gene_5UTR:
                                            logfn.write('5UTR: '+str(ref_gene_5UTR)+'\n')
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['5UTR'].extend(ref_gene_5UTR)
                                        elif ref_gene_3UTR:
                                            logfn.write('3UTR: '+str(ref_gene_3UTR)+'\n')
                                            query_aligned_full_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['3UTR'].extend(ref_gene_3UTR)           
                #1 index
                for s in query_aligned_full_regions:
                    for e in query_aligned_full_regions[s]:
                        if 'events' in query_aligned_full_regions[s][e]:
                            Salmon_hits_info[Hit_num]['100p-aligned']['query_events'].append((s,e,query_aligned_full_regions[s][e]['events']))
                        if 'overlap_ref_geneids' in query_aligned_full_regions[s][e]:
                            Salmon_hits_info[Hit_num]['100p-aligned']['overlap_ref_geneids'].append((s,e,query_aligned_full_regions[s][e]['overlap_ref_geneids']))
                
            ################ <100p ALIGNED regions ################
            #a. event peptide in <100p region of pred prot?
            #b. augustus CDS and reference CDS difference?
            # 1 index
            if query_aligned_other_prot_coor:
                query_aligned_other_prot_coor = combine_indices(query_aligned_other_prot_coor)
                hit_aligned_other_prot_coor = combine_indices(hit_aligned_other_prot_coor)

                #a.
                hits_aligned_other_regions_ep = {}
                for t in ref_ss_accessions:
                    for accession in ref_ss_accessions[t]: 
                        found = find_eventpep(all_accessions[accession]['SEQUENCE'], predicted_prot_event_peptides)
                        if found:
                            for event_peptide in found:
                                prot_coor = found[event_peptide]
                                prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                                for pc in prot_coor:
                                    ovp = interval_overlap(pc, hit_aligned_other_prot_coor)
                                    if ovp:
                                        if accession not in hits_aligned_other_regions_ep:
                                            hits_aligned_other_regions_ep[accession]= set()
                                        hits_aligned_other_regions_ep[accession].add((pc, event_peptide))
                
                if hits_aligned_other_regions_ep:
                    Salmon_hits_info[Hit_num]['<100p-aligned']['hit_events'] = hits_aligned_other_regions_ep

                #b 1 index
                query_aligned_other_regions = {}
                for s,e in query_aligned_other_prot_coor:
                    found = find_eventpep(augustus_pp_sequence[s-1:e], predicted_prot_event_peptides)
                    if found:
                        if s not in query_aligned_other_regions:
                            query_aligned_other_regions[s] = {e:{}}
                        query_aligned_other_regions[s][e]['events'] = set()
                        for event_peptide in found:
                            prot_coor = found[event_peptide]
                            prot_coor = [(i+1, i+1+len(event_peptide.replace(':',''))-1) for i in prot_coor]
                            for pc in prot_coor:
                                query_aligned_other_regions[s][e]['events'].add((pc, event_peptide))
                    
                    #b1. check area for reference gene     
                    aligned_other_genomic_coordinates = set()
                    for i in range(s-1,e):
                        aligned_other_genomic_coordinates.update(aug_prot_aa_genomic_coor[i]['aa_coor'])
                    genomic_s, genomic_e = (min(aligned_other_genomic_coordinates), max(aligned_other_genomic_coordinates))
                    
                    aligned_other_overlapping_genes = set()
                    #entry is reference gene that overlaps with predicted augustus gene
                    for ref_gene in overlap_ref_genes:
                        if ref_gene not in main_accessions_genes:
                            continue
                        #overlap gene start end (genomic)
                        b = (overlap_ref_genes[ref_gene]['start'], overlap_ref_genes[ref_gene]['end'])
                        if genomic_s <= b[1] and b[0] <= genomic_e:
                            aligned_other_overlapping_genes.add(ref_gene)
                    
                    #b2. if ref gene, overlap with augustus CDS or UTR?
                    if aligned_other_overlapping_genes:
                        for refseq_geneid in Aug_Known_overlap:
                            if refseq_geneid in nonaligned_overlapping_genes:
                                ref_gene_entry = Aug_Known_overlap[refseq_geneid]
                                #['3UTR', '5UTR', 'CDS']
                                #{'augustus': ('cds1', '2393082', '2393210'),
                                # 'protein_id': 'ENSSSAP00000054587',
                                # 'ref': ('CDS:ENSSSAP00000054587', 2393082, 2393233)}
                                for ref_gene_transcript in ref_gene_entry:
                                    ref_gene_CDS = ref_gene_entry[ref_gene_transcript]['CDS']
                                    ref_gene_5UTR = ref_gene_entry[ref_gene_transcript]['5UTR']
                                    ref_gene_3UTR = ref_gene_entry[ref_gene_transcript]['3UTR']
                                    #if not empty, list of overlapping features
                                        #[{'augustus': ('cds0', 44633435, 44633437),
                                        #  'ref': ('utr123579', 44632626, 44633448)},
                                        # {'augustus': ('cds1', 44633435, 44633448),
                                        #  'ref': ('utr123579', 44632626, 44633448)}]
                                    if ref_gene_CDS or ref_gene_5UTR or ref_gene_3UTR:
                                        logfn.write(refseq_geneid+'\t'+ref_gene_transcript+'\n')
                                        if s not in query_aligned_other_regions:
                                            query_aligned_other_regions[s] = {e:{}}
                                        if 'overlap_ref_geneids' not in query_aligned_other_regions[s][e]:
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'] = {}
                                        if refseq_geneid not in query_aligned_other_regions[s][e]['overlap_ref_geneids']:
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid] = {}
                                        if ref_gene_transcript not in query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid]:
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript] = {'CDS':[],\
                                                  '5UTR':[],'3UTR':[]}
                                        if ref_gene_CDS:
                                            logfn.write('CDS: '+str(ref_gene_CDS)+'\n')
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['CDS'].extend(ref_gene_CDS)
                                        elif ref_gene_5UTR:
                                            logfn.write('5UTR: '+str(ref_gene_5UTR)+'\n')
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['5UTR'].extend(ref_gene_5UTR)
                                        elif ref_gene_3UTR:
                                            logfn.write('3UTR: '+str(ref_gene_3UTR)+'\n')
                                            query_aligned_other_regions[s][e]['overlap_ref_geneids'][refseq_geneid][ref_gene_transcript]['3UTR'].extend(ref_gene_3UTR)           

                for s in query_aligned_other_regions:
                    for e in query_aligned_other_regions[s]:
                        if 'events' in query_aligned_other_regions[s][e]:
                            Salmon_hits_info[Hit_num]['<100p-aligned']['query_events'].append((s,e,query_aligned_other_regions[s][e]['events']))
                        if 'overlap_ref_geneids' in query_aligned_other_regions[s][e]:
                            Salmon_hits_info[Hit_num]['<100p-aligned']['overlap_ref_geneids'].append((s,e,query_aligned_other_regions[s][e]['overlap_ref_geneids']))
    
    if Salmon_hits_info:
        with open(os.path.join(augustus_parsed_ss, predprotid+'_Salmon_hits_info.pickle'), 'wb') as handle:
            pickle.dump(Salmon_hits_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if Non_Salmon_hits_info:
            with open(os.path.join(augustus_parsed_ss, predprotid+'_Non_Salmon_hits_info.pickle'), 'wb') as handle:
                pickle.dump(Non_Salmon_hits_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:                
        if Non_Salmon_hits_info:
            if Aug_Known_overlap:
                logfn.write(str(Aug_Known_overlap)+'\n')  
            with open(os.path.join(augustus_parsed_other, predprotid+'_Non_Salmon_hits_info.pickle'), 'wb') as handle:
                pickle.dump(Non_Salmon_hits_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
   
logfn.close()

print('\n\nComplete ('+str(time.time()-script_start_time)+' seconds)')