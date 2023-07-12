# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import os
from parse_reference_files import salmosalar_genomic_fna,\
    combine_indices, readFasta, rev_comp, comp,\
    combine_int_into_interval, ncbi_chrname_mapping
from Bio.Seq import Seq

def generate_event_peptide_hints(combined_Enosi_dir, EG_blastx_query_range):
    
    event_fn = os.path.join(combined_Enosi_dir, 'event_parsed.txt')
    eventgroup_summary = os.path.join(combined_Enosi_dir, 'event_group_summary.txt')
    
    source = 'enosi'
    feature = 'CDSpart'
    score = '0'
    frame = '.' # . if unknown or irreelevant
    attribute = 'pri=6;src=P'    
    enosi_hints = {}
    with open(event_fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                col = {x:xi for xi,x in enumerate(sp)}
            else:
                Event_group = sp[col['#EventGroup']]
                chromosome = sp[col['Chromosome']]
                if chromosome not in enosi_hints:
                    enosi_hints[chromosome] = {}
                if Event_group not in enosi_hints[chromosome]:
                    enosi_hints[chromosome][Event_group] = []
                strand = sp[col['Strand']]
                Peptides = sp[col['Peptides']].split('|')
                for entry in Peptides:              
                    pep, loc = entry.split('/') #1 index
                    for coordinates in loc.split(';'):
                        coordinates = coordinates.split('-')
                        start = int(coordinates[0])
                        stop = int(coordinates[1])
                        enosi_hints[chromosome][Event_group].append([chromosome,\
                                source+'@'+Event_group+'/'+pep, feature,\
                                str(start), str(stop), score, strand,\
                                frame, attribute])

    with open(eventgroup_summary,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                col = {x:xi for xi,x in enumerate(sp)}
            else:
                Event_group = sp[col['#EventGroup']]
                EG_chr = sp[col['Chromosome']]
                min_loc = int(sp[col['Location(min)']]) #1 index
                max_loc = int(sp[col['Location(max)']])
                if EG_chr not in EG_blastx_query_range:
                    EG_blastx_query_range[EG_chr] = {}
                if Event_group not in EG_blastx_query_range[EG_chr]:
                    EG_blastx_query_range[EG_chr][Event_group] = (min_loc, max_loc)
                else:
                    if min_loc< EG_blastx_query_range[EG_chr][Event_group][0]:
                        EG_blastx_query_range[EG_chr][Event_group] = (min_loc, EG_blastx_query_range[EG_chr][Event_group][1])
                    if max_loc> EG_blastx_query_range[EG_chr][Event_group][1]:
                        EG_blastx_query_range[EG_chr][Event_group] = (EG_blastx_query_range[EG_chr][Event_group][0], max_loc)
                       

    return {'enosi_hints':enosi_hints,'EG_blastx_query_range':EG_blastx_query_range}




def generate_known_peptide_hints(workflow_output_dir, \
                                 genomic_coor_output_dir,\
                                 reference_dir):
    #RefSeq
    refseq_dir = os.path.join(reference_dir, 'RefSeq')
    chromosome_genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')
    refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
    refseq_prot_fasta = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.faa')
    refseq_proteins = readFasta(refseq_prot_fasta)   
    #Ensembl
    ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
    ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')
    ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
    ensembl_proteins = readFasta(ensembl_prot_fasta)
    #contaminants
    gpm_fn = os.path.join(reference_dir, 'contaminant_db', 'gpm_cRAP_012915.fasta')
    contaminant_prot = set(readFasta(gpm_fn).keys())

    Protein_genomic_coordiates = os.path.join(genomic_coor_output_dir,\
                                    'RefSeq_Ensembl_Protein_genomic_coordinates.txt')
    Protein_CDS_coordiates = os.path.join(genomic_coor_output_dir,\
                                    'RefSeq_Ensembl_Protein_CDS_coordinates.txt')
    Peptide_genomic_coordinates = os.path.join(genomic_coor_output_dir,\
                                    'RefSeq_Ensembl_Peptide_genomic_coordinates.txt')
    
    prot_genomic_coor_fn = open(Protein_genomic_coordiates,'w')
    prot_genomic_coor_fn.write('\t'.join(['#Chromosome', 'strand',\
                                          'protein', 'genomic_coordinates'])+'\n')
    
    prot_genomic_coor_incomplete_fn = open(Protein_CDS_coordiates,'w')
    prot_genomic_coor_incomplete_fn.write('\t'.join(['#Chromosome', 'strand',\
                                          'protein', 'CDS_coordinates'])+'\n')
    
    peptide_genomic_coor_fn = open(Peptide_genomic_coordinates,'w')
    peptide_genomic_coor_fn.write('\t'.join(['#Chromosome','strand','genomic_coordinate_intervals',\
                                                'peptide','protein','protein_coordinates'])+'\n')
    
    ncbi_chrname_map = ncbi_chrname_mapping()
    chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}
        
    def KnownPeptideProtein_hints(ref_protein_CDS):
        kp_hints = {}
        for chrom in ref_protein_CDS:
            chrom_dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, chrom)
            kp_hints[chrom] = []
            for prot in ref_protein_CDS[chrom]:
                entry = ref_protein_CDS[chrom][prot]
                
                protein_seq = ''
                if prot in refseq_proteins:
                    protein_seq = refseq_proteins[prot]['seq']
                else:
                    prot = prot+'.1'
                    protein_seq = ensembl_proteins[prot]['seq']
                if not protein_seq:
                    raise ValueError('no protein_seq')
                
                prot_aacoor = {}
                #may have multiple parent rna but same protein
                cds_ids = set([x[3]['ID'] for x in entry])
                for cds_id in cds_ids:
                    entries = [x for x in entry if x[3]['ID']==cds_id]
                    
                    strand = [x[2] for x in entries][0]
                
                    cds_coor = [(x[0],x[1]) for x in entries] #1-index
                    cds_coor = combine_indices(cds_coor)
               
                    cds_seq = []
                    CDS_coor = []
                    if strand == '-':
                        cds_coor = cds_coor[::-1]
                    for b in cds_coor:
                        subseq = chrom_dna_seq[b[0]-1:b[1]]
                        if strand == '-':
                            cds_seq.append(rev_comp(subseq))
                            CDS_coor.extend(range(b[0]-1,b[1])[::-1])
                        else:
                            cds_seq.append(subseq)
                            CDS_coor.extend(range(b[0]-1,b[1]))
                    cds_join = ''.join(cds_seq).upper()
                    translated_cds = str(Seq(cds_join).translate().rstrip('*'))
                    
                    if translated_cds!=protein_seq:
                        prot_genomic_coor_incomplete_fn.write('\t'.join([chrom, strand, prot+'|'+cds_id, ';'.join([cds_seq[vi]+':'+str(v[0])+'-'+str(v[1]) for vi,v in enumerate(cds_coor)])])+'\n')
                        continue
                    
                    prot_aacoor[cds_id] = []
                    #write genomic coordinates of refseq proteins to file for use later when parsing blastp xml output
                    genomic_coordinates = []
                    for i in range(0,len(CDS_coor),3):
                        prot_aacoor[cds_id].append(CDS_coor[i:i+3])
                        genomic_coordinates.append(cds_join[i:i+3]+':'+','.join([str(v) for v in CDS_coor[i:i+3]]))
                    
                    prot_genomic_coor_fn.write('\t'.join([chrom, strand, prot+'|'+cds_id, ';'.join(genomic_coordinates)])+'\n')
        
                for peptide in Known[prot]:
                    for prot_coor in Known[prot][peptide]:
                        prot_s, prot_e = prot_coor #1 index
                        for cds_id in prot_aacoor:
                            #pep_seq = ref_proteins[prot_acc]['seq'][prot_coor[0]:prot_coor[1]+1]
                            genomic_coor = [] 
                            pep_cds = []
                            #0 index
                            for i in range(prot_s-1, prot_e):
                                ss = prot_aacoor[cds_id][i]
                                genomic_coor.extend(ss)
                                
                                subseq = ''.join([chrom_dna_seq[n] for n in ss])
                                
                                if strand == '-':
                                    subseq = comp(subseq)
                                    
                                pep_cds.append(subseq)
                                
                            pep_cds = ''.join(pep_cds)
                            translated_pep = str(Seq(pep_cds).translate())
                            if translated_pep.replace('L','I')!=peptide.replace('L','I'):
                                print( translated_pep, peptide)
                                raise ValueError('translated_pep!=peptide')
                            
                            #genomic_coor are 0index based
                            genomic_coor_intervals = combine_int_into_interval(genomic_coor)
                            # change back to 1 index
                            for gen_int in genomic_coor_intervals:
                                hint_line = [chrom, 'MSGF+@'+peptide+'['+prot+'|'+cds_id+':'+str(prot_s)+'-'+str(prot_e)+']', 'CDSpart',\
                                             gen_int[0]+1, gen_int[1]+1, \
                                             '0', strand, '.', 'pri=7;src=P']
                                kp_hints[chrom].append(hint_line)
                            
                            genomic_coor_intervals = ';'.join([str(x[0]+1)+'-'+str(x[1]+1) for x in genomic_coor_intervals])                
                            peptide_genomic_coor_fn.write('\t'.join([chrom, strand, genomic_coor_intervals, peptide,\
                                                                        prot+'|'+cds_id, str(prot_s)+'-'+str(prot_e)])+'\n')
        return kp_hints

    #exact match peptide
    Known = {}
    protein_peptide = {}
    identified_known_peptides = set()
    
    #exact match peptide
    knowndb_exact_match_fn = os.path.join(workflow_output_dir, 'ExactMatch_output',\
                   'knowndb_00000_0.upload.revCat.txt')
    known_db_target_precursor_fn = os.path.join(workflow_output_dir, 'MSGF_fdrdir',\
                   'MSGF_tsvdir_PrecursorLevelFDR_targets.fdr')
    
    known_db_target_precursors = set()
    with open(known_db_target_precursor_fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                colidx = {x:xi for xi,x in enumerate(sp)}
            else:
                known_db_target_precursors.add(''.join([aa for aa in sp[colidx['Peptide']][2:-2] if aa.isalpha()]).replace('L','I'))
    print(len(known_db_target_precursors), ' peptides passing 1% precursor level fdr')
     
    with open(knowndb_exact_match_fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0:2]=='##':
                col = {x.replace('#',''):xi for xi,x in enumerate(sp)}
            elif line[0]!='#':            
                protein = sp[col['Subject']].split('/')[0].strip()
                if protein=='No match':
                    continue
                if protein in contaminant_prot:
                    continue
                protseq = ''
                if protein in refseq_proteins:
                    protseq = refseq_proteins[protein]['seq']
                if protein in ensembl_proteins:
                    protseq = ensembl_proteins[protein]['seq']
                if not protseq:
                    if 'XXX_' not in protein:
                        print(protein)
                    continue
                match_start = int(sp[col['MatchStart']]	)
                match_end = int(sp[col['MatchEnd']])
                pepseq = protseq[match_start-1:match_end]            
                peptide = sp[col['Query']]            
                if pepseq.replace('L','I')!=peptide.replace('L','I'):
                    raise ValueError('')
                if peptide.replace('L','I') not in known_db_target_precursors:
                    continue
                if protein not in protein_peptide:
                    protein_peptide[protein] = set()
                    Known[protein] = {}
                protein_peptide[protein].add(peptide)
                if peptide not in Known[protein]:
                    Known[protein][peptide] = set()
                Known[protein][peptide].add((match_start, match_end)) # 1 index
                identified_known_peptides.add(peptide.replace('L','I'))
    
    identified_known_proteins = sorted(protein_peptide.keys())
    print('identified_known_proteins: ', len(identified_known_proteins))
    print('identified_known_peptides: ', len(identified_known_peptides))

    with open(os.path.join(genomic_coor_output_dir, 'RefSeq_Ensembl_peptides.txt'),'w') as outfile:
        outfile.write('\n'.join(identified_known_peptides))
    with open(os.path.join(genomic_coor_output_dir, 'RefSeq_Ensembl_proteins.txt'),'w') as outfile:
        outfile.write('\n'.join(identified_known_proteins))
    
    gff3_columns = ['seqid','source','type','start','end','score','strand','phase','attributes']     
    ensembl_protein_CDS = {}
    with open(ensembl_gff,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                gff_values = {gff3_columns[xi]:x for xi,x in enumerate(sp)} 
                chrname = 'chr'+gff_values['seqid']            
                if chrname not in chrname_ncbi_map:
                    continue
                feature = gff_values['type']
                if feature == 'CDS':
                    attributes = {x.split('=')[0]:x.split('=')[1] for x in gff_values['attributes'].split(';')}
                    if 'protein_id' in attributes:
                        strand = gff_values['strand']
                        start = int(gff_values['start'])
                        end = int(gff_values['end'])
                        protein_id = attributes['protein_id']
                        if protein_id+'.1' not in identified_known_proteins:
                            continue
                        if chrname not in ensembl_protein_CDS:
                            ensembl_protein_CDS[chrname] = {} 
                        if protein_id not in ensembl_protein_CDS[chrname]:
                            ensembl_protein_CDS[chrname][protein_id] = []
                        ensembl_protein_CDS[chrname][protein_id].append((start, end, strand, attributes))
    
    refseq_protein_CDS = {}
    with open(refseq_gff,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                gff_values = {gff3_columns[xi]:x for xi,x in enumerate(sp)} 
                chrname = gff_values['seqid']
                if chrname not in ncbi_chrname_map:
                    continue  
                chrname = ncbi_chrname_map[chrname]
                feature = gff_values['type']
                if feature == 'CDS':
                    attributes = {x.split('=')[0]:x.split('=')[1] for x in gff_values['attributes'].split(';')}
                    if 'protein_id' in attributes:
                        strand = gff_values['strand']
                        start = int(gff_values['start'])
                        end = int(gff_values['end'])
                        protein_id = attributes['protein_id']
                        if protein_id not in identified_known_proteins:
                            continue
                        if chrname not in refseq_protein_CDS:
                            refseq_protein_CDS[chrname] = {} 
                        if protein_id not in refseq_protein_CDS[chrname]:
                            refseq_protein_CDS[chrname][protein_id] = []
                        refseq_protein_CDS[chrname][protein_id].append((start, end, strand, attributes))

    knownpep_hints = KnownPeptideProtein_hints(refseq_protein_CDS)
    ensembl_kp_hints =  KnownPeptideProtein_hints(ensembl_protein_CDS)
    for chromosome in ensembl_kp_hints:
        if chromosome not in knownpep_hints:
            knownpep_hints[chromosome] = ensembl_kp_hints[chromosome]
        else:
            knownpep_hints[chromosome].extend(ensembl_kp_hints[chromosome])
    
    prot_genomic_coor_fn.close()
    prot_genomic_coor_incomplete_fn.close()
    peptide_genomic_coor_fn.close()


    return {'knownpep_hints':knownpep_hints,
            'ensembl_kp_hints':ensembl_kp_hints}
