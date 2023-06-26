# -*- coding: utf-8 -*-
import argparse
import os
import pickle
from parse_reference_files import readFasta, ncbi_chrname_mapping,\
    read_ensembl_gff, read_refseq_gff, read_ref_prot_gff,\
    ncbi_entrez_genpept, reference_augustus_overlap

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23',\
                    help="")
args = parser.parse_args()

print(args.dataset)


ncbi_chrname_map = ncbi_chrname_mapping()
chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

##reference files
reference_dir = os.path.join(base_dir, 'data', 'reference')

#GenPept directory
GenPept_dir = os.path.join(reference_dir, 'GenPept')

#RefSeq
print('Reading RefSeq reference files...')
refseq_dir = os.path.join(reference_dir, 'RefSeq')
 
refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
refseq_gff_results = read_refseq_gff(refseq_gff)
refseq_gff_entries = refseq_gff_results['refseq_gff_entries']
refseq_gff_start_end = refseq_gff_results['refseq_gff_start_end']
refseq_gene_coor = refseq_gff_results['refseq_gene_coor']
refseq_protein_cds = refseq_gff_results['refseq_protein_cds']

refseq_prot_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.gpff')
Salmo_salar_ref_prot_gff = read_ref_prot_gff(refseq_prot_gff)

refseq_prot_fasta = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.faa')
refseq_proteins = readFasta(refseq_prot_fasta)


#Ensembl
print('Reading Ensembl reference files...')
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')

ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')
ensembl_gff_results = read_ensembl_gff(ensembl_gff)
ensembl_gff_entries = ensembl_gff_results['ensembl_gff_entries']
ensembl_gff_start_end = ensembl_gff_results['ensembl_gff_start_end']
ensembl_gene_coor = ensembl_gff_results['ensembl_gene_coor']

ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta) 

#RefSeq and Ensembl proteins
ref_proteins = {}
for x in refseq_proteins:
    ref_proteins[x] = refseq_proteins[x]
for x in ensembl_proteins:
    ref_proteins[x] = ensembl_proteins[x]

#Augustus
print('Reading Augustus parsed files...')
augustus_dir = os.path.join(data_dir, '3_Augustus')
augustus_parsed_dir = os.path.join(augustus_dir, 'parsed_output')
augustus_pred_prot_summary_fn = os.path.join(augustus_parsed_dir,\
                                             'Augustus_pred_prot_summary.txt')
augustus_ref_overlap = os.path.join(augustus_parsed_dir,\
                                    'Augustus_reference_overlap.pickle') 
augustus_ref_overlap_entries = pickle.load(open(augustus_ref_overlap,'rb'))

#blastp
blastp_dir = os.path.join(data_dir, '4_blastp')
blastp_parsed_dir = os.path.join(blastp_dir, 'parsed_output') 
blastp_parsed_ss = os.path.join(blastp_parsed_dir, 'Salmo_salar')
blastp_parsed_other = os.path.join(blastp_parsed_dir, 'Other')
all_accessions_pickle = os.path.join(blastp_parsed_dir, 'all_accessions.pickle') 
all_accessions = pickle.load(open(all_accessions_pickle,'rb'))


RefSeq_Ensembl_Augustus_Search_dir = os.path.join(data_dir,\
                                        '5_RefSeq_Ensembl_Augustus_Search')
evidence_table_dir = os.path.join(RefSeq_Ensembl_Augustus_Search_dir,\
                                         'evidence_table')  
Non_salmo_salar_alignment_dir = os.path.join(evidence_table_dir,\
                                         'Non_salmo_salar_alignments')
Salmo_salar_alignment_dir = os.path.join(evidence_table_dir,\
                                         'Salmo_salar_alignments')
for dirnm in [Non_salmo_salar_alignment_dir, Salmo_salar_alignment_dir]:
    if not os.path.isdir(dirnm):
        os.makedirs(dirnm)



#############################################################################
logfn = open(os.path.join(evidence_table_dir, 'Write_Tables.log'),'w')

print('Read augustus_pred_prot_summary.txt...') 
augustus_pred_prot_summary = {}
Table1_header = ''
with open(augustus_pred_prot_summary_fn,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            col = {x:xi for xi,x in enumerate(sp)}
            Table1_header = line.strip().split('\t')
        else:
            proti = sp[col['#Augustus_Predicted_Prot']].split('|')[0]
            augustus_pred_prot_summary[proti] = {x:sp[col[x]] for x in col}
            #'LMDLLADSR|chrssa08|+|5424|16563576-16563602'
            augustus_pred_prot_summary[proti]['Related_Genes_CDS'] = {}
            augustus_pred_prot_summary[proti]['Reference_Overlap'] = {}
            augustus_pred_prot_summary[proti]['Aligned_region_Events'] = {}
            augustus_pred_prot_summary[proti]['Aligned_region_Ref_CDS'] ={}
            augustus_pred_prot_summary[proti]['NonAligned_region_Events'] = {}
            augustus_pred_prot_summary[proti]['NonAligned_region_Ref_CDS'] = {}
            augustus_pred_prot_summary[proti]['SemiAligned_region_Events'] = {}
            augustus_pred_prot_summary[proti]['SemiAligned_region_Ref_CDS'] = {}
            augustus_pred_prot_summary[proti]['Top5_Blastp_Hits'] = {}           
            augustus_pred_prot_summary[proti]['Hit_Event_Peptides'] = {}       
            augustus_pred_prot_summary[proti]['gene.transcript.id'] = sp[col['Augustus_gene_id']].split('|')[0]+'.'+sp[col['Augustus_transcript_id']].split('|')[0]

############################ Non_Salmon_hits ###########################################
print('Add non-Salmon hits...')
def Non_Salmo_salar_hits(input_dir, alignment_dir):
    for fn in os.listdir(input_dir):
        proti = fn.split('_')[0]
        if not os.path.isfile(os.path.join(input_dir, proti+'_Non_Salmon_hits_info.pickle')):
            continue
        pred_prot_seq = augustus_pred_prot_summary[proti]['Predicted_Prot_Seq']
        Non_salmo_salar_alignment_fn = open(os.path.join(alignment_dir, proti+'_Non_Salmon_alignments.txt'),'w')
        Non_salmo_salar_alignment_fn.write('#'+proti+'|'+str(len(pred_prot_seq))+'aa'+'\n')
        Non_Salmon_hits = pickle.load(open(os.path.join(input_dir, proti+'_Non_Salmon_hits_info.pickle'),'rb'))
        for accession in Non_Salmon_hits:            
            organism = Non_Salmon_hits[accession]['organism']
            evalue = Non_Salmon_hits[accession]['evalue']
            alignment_text = Non_Salmon_hits[accession]['alignment_text']

            if accession not in all_accessions: 
                entrez_result = ncbi_entrez_genpept(accession, GenPept_dir)
                all_accessions[accession] = entrez_result
                
            accession_protseq = all_accessions[accession]['SEQUENCE']
            #['DEFINITION', 'SEQUENCE', 'gene_id', 'VERSION', 'DBSOURCE', 'mrna', 'gene', 'ORGANISM']
            accession_def = all_accessions[accession]['DEFINITION']
            
            query_aligned_region_eventpeps = {}
            #query_aligned_region_eventpeps.add(((s,e), pc, event_peptide, EG))
            for ep in Non_Salmon_hits[accession]['query_aligned_region_eventpeps']:            
                query_align_s = ep[0][0]
                ep_s = (query_align_s-1)+ep[1][0] #1-index
                event_pep = ep[2]
                if pred_prot_seq[ep_s-1:ep_s-1+len(event_pep.replace(':',''))]!=event_pep.replace(':',''):
                    raise ValueError('pep seq error')
                if event_pep not in query_aligned_region_eventpeps:
                    query_aligned_region_eventpeps[event_pep] = set()
                query_aligned_region_eventpeps[event_pep].add((ep_s-1, ep_s-1+len(event_pep.replace(':',''))))
            
            #start query, start hit
            #end query, end hit
            augustus_pred_prot_summary[proti]['Top5_Blastp_Hits'][accession] = {'def':accession_def.split('[')[0],
                                      'organism':organism, 'evalue':evalue, 'qh_start':[], 'qh_end':[], 'accession_len':str(len(accession_protseq))+'aa'}

            if query_aligned_region_eventpeps:
                augustus_pred_prot_summary[proti]['Aligned_region_Events'][accession] = query_aligned_region_eventpeps
             
            #alignment
            Non_salmo_salar_alignment_fn.write('##Begin##\n')
            Non_salmo_salar_alignment_fn.write('##'+'|'.join([accession+'['+organism+']', str(len(accession_protseq))+'aa', 'e-value:'+str(evalue)])+'\n')
            Non_salmo_salar_alignment_fn.write('##'+accession_protseq+'\n')
            min_query_start = ''
            max_query_end = ''
            min_hit_start = ''
            max_hit_end = ''
            for hsp in alignment_text:
                for line in hsp:
                    Non_salmo_salar_alignment_fn.write(line.replace('\s',' ')+'\n')
                    if line[0]!='#':
                        if 'Query' in line:
                            q_s = int(line.split('Query')[1].strip().split()[0])
                            q_e = int(line.split('Query')[1].strip().split()[-1])
                            if not min_query_start:
                                min_query_start = q_s
                                max_query_end = q_e
                            else:
                                if q_s<min_query_start:
                                    min_query_start = q_s
                                if q_e>max_query_end:
                                    max_query_end = q_e                        
                        elif 'Sbjct' in line:
                            h_s = int(line.split('Sbjct')[1].strip().split()[0])
                            h_e = int(line.split('Sbjct')[1].strip().split()[-1])
                            if not min_hit_start:
                                min_hit_start = h_s
                                max_hit_end = h_e
                            else:
                                if h_s<min_hit_start:
                                    min_hit_start = h_s
                                if h_e>max_hit_end:
                                    max_hit_end = h_e                     
            Non_salmo_salar_alignment_fn.write('##End##'+'\n\n')
            
            augustus_pred_prot_summary[proti]['Top5_Blastp_Hits'][accession]['qh_start'] = 'Q:'+str(min_query_start)+'|H:'+str(min_hit_start)
            augustus_pred_prot_summary[proti]['Top5_Blastp_Hits'][accession]['qh_end'] = 'Q:'+str(max_query_end)+'|H:'+str(max_hit_end)
            
        Non_salmo_salar_alignment_fn.close() 
    
    return None
    
Non_Salmo_salar_hits(blastp_parsed_other, Non_salmo_salar_alignment_dir)
Non_Salmo_salar_hits(blastp_parsed_ss, Salmo_salar_alignment_dir)


############################ Salmon_hits ###############################################
print('Add Salmon hits...')
#min(Hsp_query-from) < or > min(Hsp_hit-from)
#max(Hsp_query-to) < or > max(Hsp_hit-to)
for fn in os.listdir(blastp_parsed_ss):
    proti = fn.split('_')[0]
    logfn.write('#'+proti+'\n')
    entry = augustus_pred_prot_summary[proti]

    chromosome = entry['Chromosome']
    ncbi_chrom = chrname_ncbi_map[chromosome]    
    Augustus_predictionID = entry['#Augustus_Predicted_Prot']        
    pred_prot_seq = entry['Predicted_Prot_Seq']
    
    Salmo_salar_alignment_fn = open(os.path.join(Salmo_salar_alignment_dir, proti+'_Salmon_alignments.txt'),'w')
    Salmo_salar_alignment_fn.write('#'+proti+'|'+str(len(pred_prot_seq))+'aa'+'\n')
    Salmon_hits = pickle.load(open(os.path.join(blastp_parsed_ss, proti+'_Salmon_hits_info.pickle'),'rb'))
    
    blastp_related_proteins = set()
    augustus_pred_prot_summary[proti]['Salmo_salar_Blastp_Hits'] = {}
    for Hit_num in Salmon_hits:
        add_accessions = {}
        logfn.write('Hit_num:'+str(Hit_num)+'\n')
        identical_accessions = Salmon_hits[Hit_num]['identical_accessions']
        for accession in identical_accessions['ss']:
            if accession not in all_accessions: 
                entrez_result = ncbi_entrez_genpept(accession, GenPept_dir)
                all_accessions[accession] = entrez_result
        Hit_prot_sequence = [ref_proteins[a]['seq'] if a in ref_proteins else all_accessions[a]['SEQUENCE'] for a in identical_accessions['ss']][0]
        evalue = Salmon_hits[Hit_num]['evalue']

        #main_accessions = {a:overlap_ref_gene_proteins[a] for a in ref_accessions if a in overlap_ref_gene_proteins}
        #{'gene':ref_g,'mrna':ref_t}
        main_accessions = Salmon_hits[Hit_num]['main_accessions']
        logfn.write('main_accessions '+str(main_accessions)+'\n')
        main_accession = ''
        main_accession_protseq = ''
        main_accession_def = ''
        main_accession_gene = '' #gene may be NA 
        main_accession_mrna = ''
        main_accession_chrom = ''        
        if not main_accessions:
            ref_accessions_same_chrom = Salmon_hits[Hit_num]['ref_accessions_same_chrom'] 
            non_ref_accessions_same_chrom = Salmon_hits[Hit_num]['non_ref_accessions_same_chrom'] 
            ref_accessions_diff_chrom = Salmon_hits[Hit_num]['ref_accessions_diff_chrom']  
            ref_accessions_unplaced = Salmon_hits[Hit_num]['ref_accessions_unplaced']     
            non_ref_accessions_diff_chrom = Salmon_hits[Hit_num]['non_ref_accessions_diff_chrom']
            non_ref_accessions_unplaced = Salmon_hits[Hit_num]['non_ref_accessions_unplaced'] 
            
            while not main_accession:
                if ref_accessions_same_chrom:
                    main_accession = ref_accessions_same_chrom[0]
                    main_accession_protseq = ref_proteins[main_accession]['seq']
                    main_accession_def = ref_proteins[main_accession]['def']
                    main_accession_gene = Salmo_salar_ref_prot_gff[main_accession]['gene']
                    main_accession_mrna = Salmo_salar_ref_prot_gff[main_accession]['coded_by'].split(':')[0]
                    main_accession_chrom = chromosome
                    continue
                if non_ref_accessions_same_chrom:
                    main_accession = non_ref_accessions_same_chrom[0]
                    main_accession_protseq = all_accessions[main_accession]['SEQUENCE']
                    main_accession_def = all_accessions[main_accession]['DEFINITION']
                    main_accession_gene = all_accessions[main_accession]['gene']
                    main_accession_mrna = all_accessions[main_accession]['mrna'].split(':')[0]
                    main_accession_chrom = chromosome
                    continue
                if ref_accessions_diff_chrom:
                    main_accession = ref_accessions_diff_chrom[0]
                    main_accession_protseq = ref_proteins[main_accession]['seq']
                    main_accession_def = ref_proteins[main_accession]['def']
                    main_accession_gene = Salmo_salar_ref_prot_gff[main_accession]['gene']
                    main_accession_mrna = Salmo_salar_ref_prot_gff[main_accession]['coded_by'].split(':')[0]
                    main_accession_chrom = Salmo_salar_ref_prot_gff[main_accession]['chrom']
                    continue
                if ref_accessions_unplaced:
                    main_accession = ref_accessions_unplaced[0]
                    main_accession_protseq = all_accessions[main_accession]['SEQUENCE']
                    main_accession_def = all_accessions[main_accession]['DEFINITION']
                    main_accession_gene = all_accessions[main_accession]['gene']
                    main_accession_mrna = all_accessions[main_accession]['mrna'].split(':')[0]
                    main_accession_chrom = 'Unplaced'
                    continue                
                if non_ref_accessions_diff_chrom:
                    main_accession = non_ref_accessions_diff_chrom[0]
                    main_accession_protseq = all_accessions[main_accession]['SEQUENCE']
                    main_accession_def = all_accessions[main_accession]['DEFINITION']
                    main_accession_gene = all_accessions[main_accession]['gene']
                    main_accession_mrna = all_accessions[main_accession]['mrna'].split(':')[0]
                    main_accession_chrom = Salmo_salar_ref_prot_gff[main_accession]['chrom']
                    continue
                if non_ref_accessions_unplaced:
                    main_accession = non_ref_accessions_unplaced[0]
                    main_accession_protseq = all_accessions[main_accession]['SEQUENCE']
                    main_accession_def = all_accessions[main_accession]['DEFINITION']
                    main_accession_gene = all_accessions[main_accession]['gene']
                    main_accession_mrna = all_accessions[main_accession]['mrna'].split(':')[0]
                    main_accession_chrom = 'Unplaced'
                    continue
                raise ValueError('main_accession?')
            
            if main_accession_gene=='NA':
                main_accession_gene = 'Unknown_'+main_accession_mrna
            
            add_accessions[main_accession] = {'gene':main_accession_gene, 'mrna':main_accession_mrna, 'chrom':main_accession_chrom, 'protseq':main_accession_protseq}
            blastp_related_proteins.add(main_accession)

        else:                
            for acc in main_accessions:
                main_accession = acc
                main_accession_def = ref_proteins[acc]['def']
                main_accession_gene = main_accessions[acc]['gene']
                main_accession_mrna = main_accessions[acc]['mrna']
                main_accession_protseq = ref_proteins[acc]['seq']
                main_accession_chrom = Salmo_salar_ref_prot_gff[acc]['chrom']
                
                if main_accession_gene=='NA':
                    print(main_accession_gene)
                    main_accession_gene = 'Unknown_'+main_accession_mrna
                
                add_accessions[main_accession] = {'gene':main_accession_gene, 'mrna':main_accession_mrna, 'chrom':main_accession_chrom, 'protseq':main_accession_protseq}
                blastp_related_proteins.add(main_accession)
        
        ####################################################################################
        #alignment            
        alignment_text = Salmon_hits[Hit_num]['alignment_text']
        ss_accessions = ','.join([a+'['+identical_accessions['ss'][a]+']' for a in identical_accessions['ss']])
        other_accessions = ', '.join(identical_accessions['other'])
        Salmo_salar_alignment_fn.write('##Begin##\n')
        Salmo_salar_alignment_fn.write('##Hit_sequence:'+main_accession_protseq+'\n')
        Salmo_salar_alignment_fn.write('##'+'|'.join([str(len(main_accession_protseq))+'aa', 'e-value:'+str(evalue)])+'\n')
        if ss_accessions:
            Salmo_salar_alignment_fn.write('##Salmo_salar:'+ss_accessions+'\n')
        if other_accessions:
            Salmo_salar_alignment_fn.write('##other:'+other_accessions+'\n')

        min_query_start = ''
        max_query_end = ''
        min_hit_start = ''
        max_hit_end = ''
        for hsp in alignment_text:
            for line in hsp:
                Salmo_salar_alignment_fn.write(line.replace('\s',' ')+'\n')
                if line[0]!='#':
                    if 'Query' in line:
                        q_s = int(line.split('Query')[1].strip().split()[0])
                        q_e = int(line.split('Query')[1].strip().split()[-1])
                        if not min_query_start:
                            min_query_start = q_s
                            max_query_end = q_e
                        else:
                            if q_s<min_query_start:
                                min_query_start = q_s
                            if q_e>max_query_end:
                                max_query_end = q_e                        
                    elif 'Sbjct' in line:
                        h_s = int(line.split('Sbjct')[1].strip().split()[0])
                        h_e = int(line.split('Sbjct')[1].strip().split()[-1])
                        if not min_hit_start:
                            min_hit_start = h_s
                            max_hit_end = h_e
                        else:
                            if h_s<min_hit_start:
                                min_hit_start = h_s
                            if h_e>max_hit_end:
                                max_hit_end = h_e                     
        Salmo_salar_alignment_fn.write('##End##'+'\n\n')

        ############## salmo salar hit with incomplete genomic coor info/non-ref prot only ###########
        if 'ref_ss_accessions' not in Salmon_hits[Hit_num]:
            #query_aligned_region_eventpeps.add(((s,e), pc, event_peptide, EG))
            query_aligned_region_eventpeps = {}
            for ep in Salmon_hits[Hit_num]['query_aligned_region_eventpeps']:       
                query_align_s = ep[0][0]
                query_align_e = ep[0][1]
                ep_s = (query_align_s-1)+ep[1][0] #1-index
                event_pep = ep[2]
                if pred_prot_seq[ep_s-1:ep_s-1+len(event_pep.replace(':',''))]!=event_pep.replace(':',''):
                    raise ValueError('pep seq error')
                if event_pep not in query_aligned_region_eventpeps:
                    query_aligned_region_eventpeps[event_pep] = set()
                query_aligned_region_eventpeps[event_pep].add((ep_s-1, ep_s-1+len(event_pep.replace(':',''))))
            
            if query_aligned_region_eventpeps:
                for acc in add_accessions:
                    augustus_pred_prot_summary[proti]['Aligned_region_Events'][acc] = query_aligned_region_eventpeps.keys()

        ########################### salmo salar with genomic coor info ###############################
        else:            
            #events
            #############
            p100_aligned_query_events = Salmon_hits[Hit_num]['100p-aligned']['query_events']
            if p100_aligned_query_events:
                for acc in add_accessions:
                    pre_e = set()
                    for x in p100_aligned_query_events:
                        for z in x[2]:
                            pre_e.add(z[1])
                    augustus_pred_prot_summary[proti]['Aligned_region_Events'][acc] = pre_e
                    
            non_aligned_query_events = Salmon_hits[Hit_num]['non-aligned']['query_events']
            if non_aligned_query_events:
                for acc in add_accessions:
                    pre_e = set()
                    for x in non_aligned_query_events:
                        for z in x[2]:
                            pre_e.add(z[1])
                    augustus_pred_prot_summary[proti]['NonAligned_region_Events'][acc] = pre_e
            
            p_1_99_aligned_query_events = Salmon_hits[Hit_num]['<100p-aligned']['query_events']
            if p_1_99_aligned_query_events:
                for acc in add_accessions:
                    pre_e = set()
                    for x in p_1_99_aligned_query_events:
                        for z in x[2]:
                            pre_e.add(z[1])
                    augustus_pred_prot_summary[proti]['SemiAligned_region_Events'][acc] = pre_e
                      
            #ref cds
            #############
            non_aligned_overlap_ref_geneids = Salmon_hits[Hit_num]['non-aligned']['overlap_ref_geneids']                
            if non_aligned_overlap_ref_geneids:
                for acc in add_accessions:
                    acc_gene = add_accessions[acc]['gene']
                    acc_mrna = add_accessions[acc]['mrna']
                    for en in non_aligned_overlap_ref_geneids:
                        ref_aug = en[2]
                        if acc_gene in ref_aug:
                            augustus_pred_prot_summary[proti]['NonAligned_region_Ref_CDS'][acc] = set()
                            for feature in ref_aug[acc_gene][acc_mrna]:
                                augustus_pred_prot_summary[proti]['NonAligned_region_Ref_CDS'][acc].update([z['ref'][0] for z in ref_aug[acc_gene][acc_mrna][feature]])
                
            p100_aligned_overlap_ref_geneids = Salmon_hits[Hit_num]['100p-aligned']['overlap_ref_geneids']
            if p100_aligned_overlap_ref_geneids:
                for acc in add_accessions:
                    acc_gene = add_accessions[acc]['gene']
                    acc_mrna = add_accessions[acc]['mrna']
                    for en in p100_aligned_overlap_ref_geneids:
                        ref_aug = en[2]
                        if acc_gene in ref_aug:
                            augustus_pred_prot_summary[proti]['Aligned_region_Ref_CDS'][acc] = set()
                            for feature in ref_aug[acc_gene][acc_mrna]:
                                augustus_pred_prot_summary[proti]['Aligned_region_Ref_CDS'][acc].update([z['ref'][0] for z in ref_aug[acc_gene][acc_mrna][feature]])

            p_1_99_aligned_overlap_ref_geneids = Salmon_hits[Hit_num]['<100p-aligned']['overlap_ref_geneids']
            if p_1_99_aligned_overlap_ref_geneids:
                for acc in add_accessions:
                    acc_gene = add_accessions[acc]['gene']
                    acc_mrna = add_accessions[acc]['mrna']
                    for en in p_1_99_aligned_overlap_ref_geneids:
                        ref_aug = en[2]
                        if acc_gene in ref_aug:
                            augustus_pred_prot_summary[proti]['SemiAligned_region_Ref_CDS'][acc] = set()
                            for feature in ref_aug[acc_gene][acc_mrna]:
                                augustus_pred_prot_summary[proti]['SemiAligned_region_Ref_CDS'][acc].update([z['ref'][0] for z in ref_aug[acc_gene][acc_mrna][feature]])
                                
        augustus_pred_prot_summary[proti]['Salmo_salar_Blastp_Hits'][Hit_num] = add_accessions
        #logfn.write(str(add_accessions)+'\n')

    prediction = list(augustus_ref_overlap_entries[proti]['predictions'].keys())[0]
    prediction = augustus_ref_overlap_entries[proti]['predictions'][prediction]
    aug_gene_s_e = prediction['aug_gene_s_e']
    aug_cds_s_e = prediction['aug_cds_s_e']
    aug_strand = prediction['aug_strand']
    related_proteins = set([x.split('|')[2] for x in prediction['related_genes_RefSeq'] if x.split('|')[2] not in ['5UTR','3UTR']])
    related_proteins.update([x.split('|')[2] for x in prediction['related_genes_Ensembl'] if x.split('|')[2] not in ['5UTR','3UTR']])

    blastp_protein_refoverlap = {}
    blastp_protein_norefcds = set()
    blastp_protein_sameseq = set()
    #blastp proteins. may not have checked for overlap CDS, 5UTR, 3UTR
    for p in blastp_related_proteins:
        if p not in related_proteins:
            if p in refseq_protein_cds:
                if pred_prot_seq==refseq_proteins[p]:
                    blastp_protein_sameseq.add(p)
                else:
                    refseq_overlap = reference_augustus_overlap(refseq_protein_cds[p], aug_cds_s_e, p)
                    if refseq_overlap:
                        blastp_protein_refoverlap[p] = refseq_overlap     
            else:
                blastp_protein_norefcds.add(p)
    
    augustus_pred_prot_summary[proti]['blastp_protein_refoverlap'] = blastp_protein_refoverlap
    augustus_pred_prot_summary[proti]['blastp_protein_norefcds'] = blastp_protein_norefcds
    augustus_pred_prot_summary[proti]['blastp_protein_sameseq'] = blastp_protein_sameseq

logfn.close()

######################## WRITE to table ###########################
############################# Table 1 #############################
print ('writing Table 1...')
with open(os.path.join(evidence_table_dir,'Table1.txt'),'w') as Table1:
    Table1.write('\t'.join(Table1_header)+'\t'+'Top5_Blastp_Hits (non-Salmo salar)'+'\n')
    for proti in augustus_pred_prot_summary:
        entry = augustus_pred_prot_summary[proti]
        
        Top5_Blastp_Hits = {}
        for x in entry['Top5_Blastp_Hits']:
            organism = entry['Top5_Blastp_Hits'][x]['organism']
            if organism not in Top5_Blastp_Hits:
                Top5_Blastp_Hits[organism] = set()
            Top5_Blastp_Hits[organism].add(x)
        Top5_Blastp_Hits = ';'.join([x+':'+','.join(Top5_Blastp_Hits[x]) for x in Top5_Blastp_Hits])
        if not Top5_Blastp_Hits:
            Top5_Blastp_Hits = 'NA'
            
        Table1.write('\t'.join([entry[x] for x in Table1_header])+\
                     '\t'+Top5_Blastp_Hits+'\n')

#####################################################################
############################# Table 3 ###############################
print ('writing Table 3...')
with open(os.path.join(evidence_table_dir,'Table3.txt'),'w') as Table3:
    Table3.write('\t'.join(['#Chromosome', 'Augustus_Predicted_Prot',\
                            'Organism',\
                            'Hit_Accession',\
                            'Hit_Def',\
                            'Query_Hit_Alignment_Start',\
                            'Query_Hit_Alignment_End',\
                            'Hit_evalue',\
                            'Aligned_region_Events'])+'\n')
    
    for proti in augustus_pred_prot_summary:
        entry = augustus_pred_prot_summary[proti]  
        chromosome = entry['Chromosome']
        ncbi_chrom = chrname_ncbi_map[chromosome]
        Augustus_predictionID = entry['#Augustus_Predicted_Prot']
        Predicted_Prot_Seq = entry['Predicted_Prot_Seq']
        
        Aligned_region_Events = entry['Aligned_region_Events']
        for accession in entry['Top5_Blastp_Hits']:
            organism = entry['Top5_Blastp_Hits'][accession]['organism']
            acc_def = entry['Top5_Blastp_Hits'][accession]['def']
            evalue = entry['Top5_Blastp_Hits'][accession]['evalue']
            qh_start = entry['Top5_Blastp_Hits'][accession]['qh_start']
            qh_end = entry['Top5_Blastp_Hits'][accession]['qh_end']
            hit_length = entry['Top5_Blastp_Hits'][accession]['accession_len']
            
            Ar_Events = []
            if accession in Aligned_region_Events:
                Ar_Events = Aligned_region_Events[accession]
            
            Table3.write('\t'.join([chromosome, Augustus_predictionID+'('+str(len(Predicted_Prot_Seq))+'aa)',\
                                organism, accession+'('+hit_length+')', acc_def,\
                                qh_start, qh_end, str(evalue),\
                                ';'.join(Ar_Events)])+'\n') 

#####################################################################
############################### Table 2 #############################
print( 'writing Table 2...')
with open(os.path.join(evidence_table_dir,'Table2.txt'),'w') as Table2:
    Table2.write('\t'.join(['#Chromosome', 'Augustus_Predicted_Prot',\
                            'Reference',\
                            'Related_Gene',\
                            "Related_Transcript",\
                            "Related_Protein",\
                            "Related_Transcript_5'UTR",\
                            'Related_Transcript_CDS',\
                            "Related_Transcript_3'UTR",\
                            'Overlapping_5UTR',\
                            'Overlapping_CDS',\
                            'Overlapping_3UTR',\
                            'Aligned_region_Events',\
                            'Aligned_region_Ref_CDS',\
                            'NonAligned_region_Events',\
                            'NonAligned_region_Ref_CDS',\
                            'SemiAligned_region_Events',\
                            'SemiAligned_region_Ref_CDS'])+'\n')
    
    for proti in augustus_pred_prot_summary:
        entry = augustus_pred_prot_summary[proti]  
        #if 'Salmo_salar_Blastp_Hits' in entry:
        chromosome = entry['Chromosome']
        Augustus_predictionID = entry['#Augustus_Predicted_Prot']            
        aug_prot_seq = entry['Predicted_Prot_Seq']
        Aligned_region_Events = entry['Aligned_region_Events']
        Aligned_region_Ref_CDS = entry['Aligned_region_Ref_CDS']
        NonAligned_region_Events = entry['NonAligned_region_Events']
        NonAligned_region_Ref_CDS = entry['NonAligned_region_Ref_CDS']
        SemiAligned_region_Events = entry['SemiAligned_region_Events']
        SemiAligned_region_Ref_CDS = entry['SemiAligned_region_Ref_CDS']

        prediction = list(augustus_ref_overlap_entries[proti]['predictions'].keys())[0]
        prediction = augustus_ref_overlap_entries[proti]['predictions'][prediction]
        aug_gene_s_e = prediction['aug_gene_s_e']
        aug_cds_s_e = prediction['aug_cds_s_e']
        aug_strand = prediction['aug_strand']

        #augustus CDS overlap with reference 5UTR 3UTR CDS
        Reference_CDS_Overlap = {}
        overlap_ref_genes = {}
        refseq_overlap = prediction['refseq_overlap']
        for x in refseq_overlap:
            rg = x.replace('gene-','')
            overlap_ref_genes[rg] = {'Ref':'RefSeq', 'chromosome':chromosome, 'strand':aug_strand, \
                                     'transcripts':{}, 'start':int(refseq_gene_coor[x]['start']),'end':int(refseq_gene_coor[x]['end'])}
            Reference_CDS_Overlap[rg] = {}
            for z in refseq_overlap[x]:
                rt = z.replace('rna-','')
                if rt not in overlap_ref_genes[rg]['transcripts']:
                    Ref_seq_entry = refseq_gff_entries[chromosome][aug_strand][x][z]
                    mrna_start = Ref_seq_entry['start']
                    mrna_end = Ref_seq_entry['end']
                    mrna_strand = Ref_seq_entry['strand']
                    mrna_cds = []
                    fiveprimeUTR = []
                    threeprimeUTR = []
                    if 'CDS' in Ref_seq_entry:
                        mrna_cds = sorted([(cds['start'], cds['end'], cds['strand'], cds['attributes']['ID'], cds['protein_id']) for cds in Ref_seq_entry['CDS']])
                    if '5UTR' in Ref_seq_entry:
                        fiveprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['attributes']['ID']) for utr in Ref_seq_entry['5UTR']])
                    if '3UTR' in Ref_seq_entry:
                        threeprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['attributes']['ID']) for utr in Ref_seq_entry['3UTR']])                            
                    overlap_ref_genes[rg]['transcripts'][rt] = {'strand':mrna_strand, 'start':mrna_start, 'end':mrna_end, 'CDS':mrna_cds, '5UTR':fiveprimeUTR, '3UTR':threeprimeUTR}
                Reference_CDS_Overlap[rg][rt] = refseq_overlap[x][z]
                for t in ['CDS','5UTR','3UTR']:
                    if t not in Reference_CDS_Overlap[rg][rt]:
                        Reference_CDS_Overlap[rg][rt][t] = []                 

        ensembl_overlap = prediction['ensembl_overlap']
        for x in ensembl_overlap:
            rg = x.replace('gene:','')
            overlap_ref_genes[rg] = {'Ref':'Ensembl','chromosome':chromosome, 'strand':aug_strand, 'transcripts':{},\
                                     'start':int(ensembl_gene_coor[x]['start']),'end':int(ensembl_gene_coor[x]['end'])}
            Reference_CDS_Overlap[rg] = {} 
            for z in ensembl_overlap[x]:
                rt = z.replace('transcript:','')
                if rt not in overlap_ref_genes[rg]['transcripts']:
                    Ensembl_seq_entry = ensembl_gff_entries[chromosome][aug_strand][x][z]
                    mrna_start = Ensembl_seq_entry['start']
                    mrna_end = Ensembl_seq_entry['end']
                    mrna_strand = Ensembl_seq_entry['strand']
                    mrna_cds = []
                    fiveprimeUTR = []
                    threeprimeUTR = []
                    if 'CDS' in Ensembl_seq_entry:
                        mrna_cds = sorted([(cds['start'], cds['end'], cds['strand'], cds['attributes']['ID'], cds['protein_id']) for cds in Ensembl_seq_entry['CDS']])
                    if '5UTR' in Ensembl_seq_entry:
                        fiveprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['ID']) for utr in Ensembl_seq_entry['5UTR']])
                    if '3UTR' in Ensembl_seq_entry:
                        threeprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['ID']) for utr in Ensembl_seq_entry['3UTR']])                            
                    overlap_ref_genes[rg]['transcripts'][rt] = {'strand':mrna_strand, 'start':mrna_start, 'end':mrna_end, 'CDS':mrna_cds, '5UTR':fiveprimeUTR, '3UTR':threeprimeUTR}
                Reference_CDS_Overlap[rg][rt] = ensembl_overlap[x][z]
                for t in ['CDS','5UTR','3UTR']:
                    if t not in Reference_CDS_Overlap[rg][rt]:
                        Reference_CDS_Overlap[rg][rt][t] = []
        
        #blastp proteins. may not have checked for overlap CDS, 5UTR, 3UTR                    
        #{'chromosome':chrname, 'strand':strand, 'start':start, 'end':end, 'attributes':attributes,
        # 'parent_gene':parent_gene, 'parent_transcript':parent_transcript}
        if 'blastp_protein_refoverlap' in entry:
            blastp_protein_refoverlap = entry['blastp_protein_refoverlap']
            for p in blastp_protein_refoverlap:
                p_entry = blastp_protein_refoverlap[p]
                x = p_entry['parent_gene']
                rg = x.replace('gene-','')
                rg_chrom = refseq_gene_coor[x]['chrom']
                rg_strand = refseq_gene_coor[x]['strand']
                #refseq_gene_coor[gene_id] = {'start':start, 'end':end, 'strand': strand, 'chrom':chrname, 'attributes':attributes}
                overlap_ref_genes[rg] = {'Ref':'Blastp','chromosome':rg_chrom, 'strand':rg_strand, 'transcripts':{}, 'start':int(refseq_gene_coor[x]['start']),'end':int(refseq_gene_coor[x]['end'])}
                z = p_entry['parent_transcript']
                rt = z.replace('rna-','')
                if rt not in overlap_ref_genes[rg]['transcripts']:
                    Ref_seq_entry = refseq_gff_entries[rg_chrom][rg_strand][x][z]
                    mrna_start = Ref_seq_entry['start']
                    mrna_end = Ref_seq_entry['end']
                    mrna_strand = Ref_seq_entry['strand']
                    mrna_cds = []
                    fiveprimeUTR = []
                    threeprimeUTR = []
                    if 'CDS' in Ref_seq_entry:
                        mrna_cds = sorted([(cds['start'], cds['end'], cds['strand'], cds['attributes']['ID'], cds['protein_id']) for cds in Ref_seq_entry['CDS']])
                    if '5UTR' in Ref_seq_entry:
                        fiveprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['attributes']['ID']) for utr in Ref_seq_entry['5UTR']])
                    if '3UTR' in Ref_seq_entry:
                        threeprimeUTR = sorted([(utr['start'], utr['end'], utr['strand'], utr['attributes']['ID']) for utr in Ref_seq_entry['3UTR']])                            
                    overlap_ref_genes[rg]['transcripts'][rt] = {'strand':mrna_strand, 'start':mrna_start, 'end':mrna_end, 'CDS':mrna_cds, '5UTR':fiveprimeUTR, '3UTR':threeprimeUTR}
                if rg not in Reference_CDS_Overlap:
                    Reference_CDS_Overlap[rg] = {}
                if rt not in Reference_CDS_Overlap[rg]:
                    Reference_CDS_Overlap[rg][rt] = {}
                for t in ['CDS','5UTR','3UTR']:
                    if t in p_entry:
                        Reference_CDS_Overlap[rg][rt][t] = p_entry[t]
                    else:
                        Reference_CDS_Overlap[rg][rt][t] = []
        
        for rg in overlap_ref_genes:
            gene_chr = overlap_ref_genes[rg]['chromosome']
            gene_strand = overlap_ref_genes[rg]['strand']
            gene_start = overlap_ref_genes[rg]['start']
            gene_end = overlap_ref_genes[rg]['end']
            Ref = overlap_ref_genes[rg]['Ref']
            gene_ncbi_chr = chrname_ncbi_map[gene_chr]
            gene = rg+'@'+gene_chr+'/'+gene_ncbi_chr+'('+gene_strand+')'+':'+str(gene_start)+'-'+str(gene_end)
            for rt in overlap_ref_genes[rg]['transcripts']:
                ref_entry = overlap_ref_genes[rg]['transcripts'][rt]
                #CDS 5UTR 3UTR reference coordinates
                CDS = ''
                fiveUTR = ''
                threeUTR = ''
                protein = ''
                mrna_start = ref_entry['start']
                mrna_end = ref_entry['end']
                mrna_strand = ref_entry['strand']
                transcript = rt+'@('+mrna_strand+'):'+str(mrna_start)+'-'+str(mrna_end)
                if ref_entry['CDS']:
                    CDS = ref_entry['CDS']
                    cds_strand_id = set([cds[3]+'('+cds[2]+')' for cds in CDS])
                    if len(cds_strand_id)>1:
                        raise ValueError('')
                    protein = list(set([cds[4] for cds in CDS]))[0]
                    cds_strand_id = list(cds_strand_id)[0]
                    CDS = cds_strand_id+':'+''.join(['@'+str(cds[0])+'-'+str(cds[1]) for cds in CDS])
                if ref_entry['5UTR']:
                    fiveUTR = ''.join(['@'+utr[3]+'('+utr[2]+'):'+str(utr[0])+'-'+str(utr[1]) for utr in ref_entry['5UTR']])
                if ref_entry['3UTR']:
                    threeUTR = ''.join(['@'+utr[3]+'('+utr[2]+'):'+str(utr[0])+'-'+str(utr[1]) for utr in ref_entry['3UTR']])                   
                
                #Augustus overlap with reference CDS 5UTR 3UTR
                overlap_CDS = []
                overlap_5UTR = []
                overlap_3UTR = []
                if Reference_CDS_Overlap[rg][rt]['CDS']:
                    for ann in Reference_CDS_Overlap[rg][rt]['CDS']:
                        aug_s = int(ann['augustus'][1])
                        aug_e = int(ann['augustus'][2])
                        ref_s = int(ann['ref'][1])
                        ref_e = int(ann['ref'][2])
                        overlap_int = (max(aug_s, ref_s), min(aug_e, ref_e))
                        overlap_CDS.append((overlap_int[0], ann['ref'][0]+'@'+str(overlap_int[0])+'-'+str(overlap_int[1])+'('+str(overlap_int[1]-overlap_int[0])+' bp)'))
                    overlap_CDS = [cds[1] for cds in sorted(overlap_CDS)]
                if Reference_CDS_Overlap[rg][rt]['5UTR']:
                    for ann in Reference_CDS_Overlap[rg][rt]['5UTR']:
                        aug_s = int(ann['augustus'][1])
                        aug_e = int(ann['augustus'][2])
                        ref_s = int(ann['ref'][1])
                        ref_e = int(ann['ref'][2])
                        overlap_int = (max(aug_s, ref_s), min(aug_e, ref_e))
                        overlap_5UTR.append((overlap_int[0], ann['ref'][0]+'@'+str(overlap_int[0])+'-'+str(overlap_int[1])+'('+str(overlap_int[1]-overlap_int[0])+' bp)'))
                    overlap_5UTR = [utr[1] for utr in sorted(overlap_5UTR)]
                if Reference_CDS_Overlap[rg][rt]['3UTR']:
                    for ann in Reference_CDS_Overlap[rg][rt]['3UTR']:
                        aug_s = int(ann['augustus'][1])
                        aug_e = int(ann['augustus'][2])
                        ref_s = int(ann['ref'][1])
                        ref_e = int(ann['ref'][2])
                        overlap_int = (max(aug_s, ref_s), min(aug_e, ref_e))
                        overlap_3UTR.append((overlap_int[0], ann['ref'][0]+'@'+str(overlap_int[0])+'-'+str(overlap_int[1])+'('+str(overlap_int[1]-overlap_int[0])+' bp)'))
                    overlap_3UTR = [utr[1] for utr in sorted(overlap_3UTR)]
                
                Ar_Events = set()
                if protein in Aligned_region_Events:
                    Ar_Events.update(Aligned_region_Events[protein]) 
                Nr_Events = set()
                if protein in NonAligned_region_Events:
                    Nr_Events.update(NonAligned_region_Events[protein])
                SAr_Event = set()                            
                if protein in SemiAligned_region_Events:
                    SAr_Event.update(SemiAligned_region_Events[protein])
            
                Ar_RefCDS = set()                            
                if protein in Aligned_region_Ref_CDS:
                    Ar_RefCDS.update(Aligned_region_Ref_CDS[protein])
                Nr_RefCDS = set()
                if protein in NonAligned_region_Ref_CDS:
                    Nr_RefCDS.update(NonAligned_region_Ref_CDS[protein])
                SAr_RefCDS = set()                                                                                               
                if protein in SemiAligned_region_Ref_CDS:
                    SAr_RefCDS.update(SemiAligned_region_Ref_CDS[protein])

                Table2.write('\t'.join([chromosome, Augustus_predictionID, Ref,\
                                        gene, transcript, protein,\
                                        fiveUTR, CDS, threeUTR,\
                                        ';'.join(overlap_5UTR),\
                                        ';'.join(overlap_CDS),\
                                        ';'.join(overlap_3UTR),\
                                        ';'.join(Ar_Events),\
                                        ';'.join(Ar_RefCDS),\
                                        ';'.join(Nr_Events),\
                                        ';'.join(Nr_RefCDS),\
                                        ';'.join(SAr_Event),\
                                        ';'.join(SAr_RefCDS)])+'\n')

        if 'blastp_protein_sameseq' in entry:
            blastp_protein_sameseq = entry['blastp_protein_sameseq']
            if blastp_protein_sameseq:
                print (proti, blastp_protein_sameseq)
        
        if 'blastp_protein_norefcds' in entry:
            blastp_protein_norefcds = entry['blastp_protein_norefcds']            
            for protein in blastp_protein_norefcds:
                Ar_Events = set()
                if protein in Aligned_region_Events:
                    Ar_Events.update(Aligned_region_Events[protein]) 
                Nr_Events = set()
                if protein in NonAligned_region_Events:
                    Nr_Events.update(NonAligned_region_Events[protein])
                SAr_Event = set()                            
                if protein in SemiAligned_region_Events:
                    SAr_Event.update(SemiAligned_region_Events[protein])
    
                Ar_RefCDS = set()                            
                if protein in Aligned_region_Ref_CDS:
                    Ar_RefCDS.update(Aligned_region_Ref_CDS[protein])
                Nr_RefCDS = set()
                if protein in NonAligned_region_Ref_CDS:
                    Nr_RefCDS.update(NonAligned_region_Ref_CDS[protein])
                SAr_RefCDS = set()                                                                                               
                if protein in SemiAligned_region_Ref_CDS:
                    SAr_RefCDS.update(SemiAligned_region_Ref_CDS[protein])
                
                Table2.write('\t'.join(['Unplaced?', Augustus_predictionID, 'Blastp',\
                                        'Unplaced?', 'Unplaced?', protein,\
                                        '', '', '',\
                                        '',\
                                        '',\
                                        '',\
                                        ';'.join(Ar_Events),\
                                        ';'.join(Ar_RefCDS),\
                                        ';'.join(Nr_Events),\
                                        ';'.join(Nr_RefCDS),\
                                        ';'.join(SAr_Event),\
                                        ';'.join(SAr_RefCDS)])+'\n')
