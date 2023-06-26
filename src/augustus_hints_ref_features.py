# -*- coding: utf-8 -*-
import os 
from parse_reference_files import ncbi_chrname_mapping

def generate_ref_features_hints(reference_dir):
    #RefSeq
    refseq_dir = os.path.join(reference_dir, 'RefSeq')
    refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
    #Ensembl
    ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
    ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')

    augustus_feature_term = set(['CDS','exon','UTRpart'])
    gff3_columns = ['seqid','source','type','start','end','score','strand','phase','attributes']

    ncbi_chrname_map = ncbi_chrname_mapping()
    chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}
     
    ##################### All GFF hints #####################
    refseq_gff_entries = {}
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
                if feature in ['five_prime_UTR','three_prime_UTR']:
                    feature = 'UTRpart'                         
                if feature not in augustus_feature_term:
                    continue    
                strand = gff_values['strand']
                start = int(gff_values['start'])
                end = int(gff_values['end'])    
                source = 'RefSeqGFF@'+gff_values['source']
                score = '0'
                frame = gff_values['phase'] # . if unknown or irreelevant
                attribute = 'pri=5;src=M'
                if chrname not in refseq_gff_entries:
                    refseq_gff_entries[chrname] = []
                refseq_gff_entries[chrname].append((chrname, source, feature, start, end, score, strand, frame, attribute))
    
    ensembl_gff_entries = {}
    with open(ensembl_gff,'r') as infile:
        for line in infile:
            if line[0]!='#':
                sp = line.strip().split('\t')
                gff_values = {gff3_columns[xi]:x for xi,x in enumerate(sp)}                  
                chrname = 'chr'+gff_values['seqid']
                if chrname not in chrname_ncbi_map:
                    continue                
                feature = gff_values['type']    
                if feature in ['five_prime_UTR','three_prime_UTR']:
                    feature = 'UTRpart'                         
                if feature not in augustus_feature_term:
                    continue    
                strand = gff_values['strand']
                start = int(gff_values['start'])
                end = int(gff_values['end'])    
                source = 'EnsemblGFF@'+gff_values['source']
                score = '0'
                frame = gff_values['phase'] # . if unknown or irreelevant
                attribute = 'pri=5;src=M'
                if chrname not in ensembl_gff_entries:
                    ensembl_gff_entries[chrname] = []
                ensembl_gff_entries[chrname].append((chrname, source, feature, start, end, score, strand, frame, attribute))    
    
    return {'refseq_gff_entries':refseq_gff_entries,
            'ensembl_gff_entries':ensembl_gff_entries}