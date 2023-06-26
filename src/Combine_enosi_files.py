# -*- coding: utf-8 -*-
import time
import os
import re
import shutil
import pickle
import argparse
from parse_reference_files import read_refseq_gff, read_ensembl_gff,\
    readFasta, salmosalar_genomic_fna, reference_gene_overlap,\
        translate, rev_comp, find_eventpep

ncbi_range = 60000
base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset", default = '2019_tissues_n23', \
                    help="")
args = parser.parse_args()

print(args.dataset)

start_time = time.time()

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

dataset_dir = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search', 'workflow_output')
output_dir = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search', 'combined_Enosi_Output')

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

#output files     
RefSeq_with_supported_ensembl_fasta = os.path.join(output_dir,\
                                        'GCF_000233375.1_ICSASG_v2_protein+Ensembl.fasta')
event_fn = os.path.join(output_dir, 'event_parsed.txt')
eventgroup_summary = os.path.join(output_dir, 'event_group_summary.txt')
agg_events = os.path.join(output_dir, 'all_events.pickle')

##reference files
reference_dir = os.path.join(base_dir, 'data', 'reference')
print('Reading RefSeq reference files...')
#RefSeq
refseq_dir = os.path.join(reference_dir, 'RefSeq')
genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')
refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
refseq_prot_fasta = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.faa')
refseq_proteins = readFasta(refseq_prot_fasta)

refseq_gff_results = read_refseq_gff(refseq_gff)
refseq_gff_entries = refseq_gff_results['refseq_gff_entries']
refseq_gff_start_end = refseq_gff_results['refseq_gff_start_end']
refseq_mRNA_protein = refseq_gff_results['refseq_mRNA_protein']

print('Reading Ensembl reference files...')
#Ensembl
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')
ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta)

ensembl_gff_results = read_ensembl_gff(ensembl_gff)
ensembl_gff_entries = ensembl_gff_results['ensembl_gff_entries']
ensembl_gff_start_end = ensembl_gff_results['ensembl_gff_start_end']
ensembl_mRNA_protein = ensembl_gff_results['ensembl_mRNA_protein']


print('Parsing proteogenomic events files...')
#read event.txt files from workflow output
#output combined events file
combined_events = {}
header = ''
colidx = ''
for subdir in os.listdir(dataset_dir):
    #copy and rename to event.txt
    Output_directory = os.path.join(dataset_dir, subdir)
    event_txt_fn = os.path.join(Output_directory, 'event.txt')
    if not os.path.isfile(event_txt_fn):
        downloaded_event_txt_fn = [x for x in os.listdir(Output_directory)\
                                   if re.search("^CUSTOM_COMPREHENSIVE_DB_SEARCH-.*all_events-main.tsv$", x)]
        if not downloaded_event_txt_fn:
            raise ValueError('Place CUSTOM_COMPREHENSIVE_DB_SEARCH-xxxxxxxx-all_events-main.tsv in Enosi_Output directory!')
        print('Copying '+downloaded_event_txt_fn[0]+' as event.txt...')
        downloaded_event_txt_fn = os.path.join(Output_directory, downloaded_event_txt_fn[0])
        shutil.copy2(downloaded_event_txt_fn, event_txt_fn)
    #parse event.txt
    with open(event_txt_fn,'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                colidx = {x:xi for xi,x in enumerate(sp)}
                header = line.strip()
            else:
                event_num = subdir+'_'+sp[colidx['#Num']]
                sp[colidx['#Num']] = event_num
                chromosome = sp[colidx['Chromosome']]
                if chromosome not in combined_events:
                    combined_events[chromosome] = []
                combined_events[chromosome].append(sp)

print('Checking event overlap with Ensembl annotations...')
#check if events overlap with Ensembl annotations
#logfn = open(os.path.join(output_dir, 'all_events.log'),'w')
related_genes_not_in_refseq = set()
enosi_events = {}
for chromosome in combined_events:         
    if chromosome not in enosi_events:
        enosi_events[chromosome] = {}
    dna_seq = salmosalar_genomic_fna(genomic_fna_dir, chromosome)
    for entry in combined_events[chromosome]:
        event_num = entry[colidx['#Num']]
        event_type = entry[colidx['Event']]
        location = entry[colidx['Location']]
        strand = entry[colidx['Strand']]
        peptide = entry[colidx['Peptide']]
        spec_count = entry[colidx['spec_count']]
        #issue with RelatedGene column created by Enosi script
        related_genes = [x.replace('NA','') for x in entry[colidx['RelatedGene']].rstrip(';').split(';') if x.replace('NA','').strip()]
        
        groupinfo = {peptide:[[sorted((int(z.split('-')[0]),int(z.split('-')[1])))\
                               for z in location.split(';')]]}
            
        for x in entry[colidx['GroupInfo']].rstrip('|').split('|'):
            if not x:
                continue
            pep,loc = x.split('/')
            if pep not in groupinfo:
                groupinfo[pep] = []
            groupinfo[pep].append([sorted((int(z.split('-')[0]),int(z.split('-')[1])))\
                                   for z in loc.split(';')])
            
        relatedgene = []
        for rg in related_genes:
            #issue with RelatedGene column created by Enosi script
            if rg not in refseq_mRNA_protein:
                #logfn.write(event_num+'\tgene not in reference: '+rg+'\n')
                related_genes_not_in_refseq.add(rg)
            #look for peptide in reference protein sequence. If present, remove.
            else:
                relatedgene.append(rg)
                prot_name = refseq_mRNA_protein[rg]
                prot_name = list(prot_name)[0]
                found = find_eventpep(refseq_proteins[prot_name]['seq'], groupinfo.keys())
                if found:
                    for p in found:
                        #logfn.write(event_num+'\t'+p+'\t'+'found event peptide in: '+prot_name+'\n')
                        del groupinfo[p]
        
        EventNum_peps = {}
        min_max_location = []
        #using locations of peptides, check annotations
        for pep in groupinfo:
            #disregard peptide if length of splice part <=3 amino acids
            #SALS:ATGGMTATGQ:AK
            pep_parts = set(['remove' if len(x)<=3 else '' for x in pep.split(':')])
            if 'remove' in pep_parts:
                #logfn.write(event_num+'\t'+pep+'\n')
                continue
            
            for peptide_location in groupinfo[pep]:
                #disregard peptide if length of splice part <=3 amino acids
                loc_parts = set(['remove' if abs(x[1]-x[0])<=9 else '' for x in peptide_location])
                if 'remove' in loc_parts:
                    #logfn.write(event_num+'\t'+pep+'\t'+str(peptide_location)+'\n')
                    continue
    
                #check peptide coordinates
                dseq = ''.join([dna_seq[l[0]:l[1]] for l in peptide_location])
                if strand == '+':
                    expected_pepseq = translate(dseq)
                else:
                    expected_pepseq = translate(rev_comp(dseq))
                if expected_pepseq!=pep.replace(':',''):
                    print( pep, expected_pepseq)
                    err_msg = 'pep seq doesnt match'
                    raise ValueError(err_msg)
                
                #use 1 index
                index1_peptide_location = [(x[0]+1,x[1]) for x in peptide_location]
                loc_start = min([min(x) for x in index1_peptide_location])
                loc_end = max([max(x) for x in index1_peptide_location])
    
                if not min_max_location:
                    min_max_location = [loc_start, loc_end]
                else:
                    if loc_start<min_max_location[0]:
                        min_max_location[0] = loc_start
                    if loc_end>min_max_location[1]:
                        min_max_location[1] = loc_end
                
                if pep not in EventNum_peps:
                    EventNum_peps[pep] = []
                EventNum_peps[pep].append(index1_peptide_location)
        
        if not EventNum_peps:
            #logfn.write('EventNum_peps empty: '+event_num+'\n')
            continue
        
        splice_coor = {}
        other_coor = {}
        for pep in EventNum_peps:
            for peploc in EventNum_peps[pep]:
                if len(peploc)==1:
                    if pep not in other_coor:
                        other_coor[pep] = []
                    other_coor[pep].append(peploc)                         
                elif len(peploc)>=2:
                    if pep not in splice_coor:
                        splice_coor[pep] = []
                    if len(peploc)==2:
                        splice_coor[pep].append([(peploc[0][1], peploc[1][0])])
                    if len(peploc)>2:
                        locs = []
                        for i,il in enumerate(peploc):
                            if i==len(peploc)-1:
                                continue
                            locs.append((peploc[i][1], peploc[i+1][0]))
                        splice_coor[pep].append(locs)
        
        ensembl_overlap = []
        refseq_overlap = []
        if splice_coor:                    
            ensembl_overlap = reference_gene_overlap(ensembl_gff_start_end[chromosome][strand],\
                                                      ensembl_gff_entries[chromosome][strand],\
                                                          min_max_location[0], min_max_location[1],\
                                                              splice_coor, other_coor,\
                                                                  ensembl_mRNA_protein, ensembl_proteins)
            refseq_overlap = reference_gene_overlap(refseq_gff_start_end[chromosome][strand],\
                                                    refseq_gff_entries[chromosome][strand],\
                                                        min_max_location[0], min_max_location[1],\
                                                            splice_coor, other_coor,\
                                                                ensembl_mRNA_protein, ensembl_proteins)
            
        enosi_events[chromosome][event_num]= {'PSMs':spec_count, 'min_max_location':min_max_location,
                                                'strand':strand, 'event_type':event_type,
                                                'peptides':EventNum_peps,'RelatedGene':relatedgene,
                                                'Ensembl_overlap':ensembl_overlap,\
                                                'Refseq_overlap':refseq_overlap,\
                                                'splice_coor':splice_coor}
#logfn.close()
with open(agg_events, 'wb') as handle:
    pickle.dump(enosi_events , handle, protocol=pickle.HIGHEST_PROTOCOL)

print('Writing event group files...')
all_cmds = []
number_of_commands = 0
EG_group_info = {}
with open(event_fn,'w') as outfile, open(eventgroup_summary,'w') as EGout:
    outfile.write('\t'.join(['#EventGroup','EventNum', 'Event', 'Chromosome',\
                             'Strand','Location(min)','Location(max)','Peptides',\
                                 'RelatedGene','Ensembl_overlap'])+'\n')
    EGout.write('\t'.join(['#EventGroup','EventNum', 'Chromosome','Strand',\
                           'Location(min)','Location(max)','Peptides'])+'\n')
    for chromosome in enosi_events:
        dna_seq = salmosalar_genomic_fna(genomic_fna_dir, chromosome)
        
        related_events = {}
        for event_num in enosi_events[chromosome]:
            event_loc_min = enosi_events[chromosome][event_num]['min_max_location'][0]
            if event_loc_min<ncbi_range:
                event_loc_min = 0
            else:
                event_loc_min = event_loc_min-ncbi_range
            event_loc_max = enosi_events[chromosome][event_num]['min_max_location'][1]+ncbi_range 
            related_events[event_num] = set([event_num])
            for other_event in enosi_events[chromosome]:
                if other_event!=event_num:
                    loc_min = enosi_events[chromosome][other_event]['min_max_location'][0]
                    loc_max = enosi_events[chromosome][other_event]['min_max_location'][1]
                    if event_loc_min<=loc_min<=event_loc_max and event_loc_min<=loc_max<=event_loc_max:
                        related_events[event_num].add(other_event)
        
        #combine subsets
        sorted_rg = set([';'.join(sorted(related_events[en])) for en in related_events])
        sorted_rg = sorted([sorted(x.split(';')) for x in sorted_rg])
        groups = sorted_rg[0:]
        for g in sorted_rg:
            for x in sorted_rg:
                if g!=x:
                    if len(set(g)-set(x))==0:
                        if g in groups:
                            groups.remove(g)
        
        EG_num = 0
        for group_events in groups:
            EG_num+=1
            EG = chromosome+'_EG_'+str(EG_num)
            EG_eventnums = ';'.join(group_events)
            EG_group_min = ''
            EG_group_max = ''
            EG_group_strand = set()
            EG_group_peptides = set()
            for event_num in group_events:
                event_type = enosi_events[chromosome][event_num]['event_type']
                Ensembl_overlap = []                
                for rg in enosi_events[chromosome][event_num]['Ensembl_overlap']:
                    for rt in enosi_events[chromosome][event_num]['Ensembl_overlap'][rg]:
                        overlapping_events = enosi_events[chromosome][event_num]['Ensembl_overlap'][rg][rt]
                        transcriptID = rt
                        if 'transcript:' in rt:
                            transcriptID = rt.split('transcript:')[1]
                        protein_id = ensembl_mRNA_protein[transcriptID]
                        if protein_id:
                            for p in protein_id:
                                protid_v = p+'.1'
                                prot_seq = ensembl_proteins[protid_v]['seq']
                                #check if peptide in prot_seq
                                found = set()
                                for ep in overlapping_events:
                                    index = [str(m.start()) for m in re.finditer(ep.replace(':',''), prot_seq)]
                                    if index:
                                        found.add(ep)
                                #as long as an event peptide is found, add protein to fasta
                                if found:
                                    if EG not in EG_group_info:
                                        EG_group_info[EG] = {}
                                    if protid_v not in EG_group_info[EG]:
                                        EG_group_info[EG][protid_v] = set()
                                    EG_group_info[EG][protid_v].update(found)
                        else:
                            print( rt)
                        
                        Ensembl_overlap.append(rg+'|'+rt+'|'+','.join(sorted(protein_id))+'|'+','.join(sorted(overlapping_events)))                              

                strand = enosi_events[chromosome][event_num]['strand']
                loc_min = str(enosi_events[chromosome][event_num]['min_max_location'][0])
                loc_max = str(enosi_events[chromosome][event_num]['min_max_location'][1])
                if not EG_group_min:
                    EG_group_min = int(loc_min)
                    EG_group_max = int(loc_max)
                else:
                    if int(loc_min)<EG_group_min:
                        EG_group_min = int(loc_min)
                    if int(loc_max)>EG_group_max:
                        EG_group_max = int(loc_max)
                if not EG_group_strand:
                    EG_group_strand = set([strand])
                else:
                    if strand!=EG_group_strand:
                        EG_group_strand.add(strand)
                RelatedGene = ';'.join([z for z in enosi_events[chromosome][event_num]['RelatedGene'] if z and z!='LOC'])
                peptides = enosi_events[chromosome][event_num]['peptides']
                pep_loc = []
                for pep in peptides:
                    for index1_peptide_location in peptides[pep]:
                        pep_loc.append(pep+'/'+';'.join([str(x[0])+'-'+str(x[1]) for x in index1_peptide_location]))
                EG_group_peptides.update(pep_loc)
                              
                if not Ensembl_overlap:
                    Ensembl_overlap = ['NA']
                if not pep_loc:
                    pep_loc = ['NA']
                
                outfile.write('\t'.join([EG, event_num, event_type, chromosome,\
                                         strand, loc_min, loc_max,'|'.join(pep_loc),\
                                             RelatedGene, ';'.join(Ensembl_overlap)])+'\n')
            
            EGout.write('\t'.join([EG, EG_eventnums, chromosome, ','.join(EG_group_strand),\
                                   str(EG_group_min), str(EG_group_max), '|'.join(EG_group_peptides)])+'\n')


print('Creating RefSeq+Ensembl database...')
# write databse with event supported Ensembl sequences
all_refseq_sequences = set([refseq_proteins[x]['seq'] for x in refseq_proteins])
ensembl_added = []
with open(RefSeq_with_supported_ensembl_fasta,'w') as outfile:
    with open(refseq_prot_fasta,'r') as infile:
        for line in infile:
            outfile.write(line)
    for EG in EG_group_info:
        chrom = EG.split('_')[0]
        #print EG, EG_group_info[EG]
        sequences = set([ensembl_proteins[v]['seq'] for v in EG_group_info[EG]])
        for prot_seq in sequences:
            prot_header = set([v+' /'+','.join(EG_group_info[EG][v]) for v in EG_group_info[EG] if ensembl_proteins[v]['seq']==prot_seq])
            prot_header = '>'+';'.join(prot_header)
            if prot_seq in all_refseq_sequences:
                #print 'skipping ',prot_header
                continue
            sequence_split = [prot_seq[i:i+80] for i in range(0, len(prot_seq), 80)]
            if ''.join(sequence_split)!=prot_seq:
                raise ValueError('')
            outfile.write(prot_header+'\n'+'\n'.join(sequence_split)+'\n')
            ensembl_added.append(prot_seq)
print (len(ensembl_added), len(set(ensembl_added)))

print('Combine Enosi Files Complete ('+str(time.time()-start_time)+' seconds)')