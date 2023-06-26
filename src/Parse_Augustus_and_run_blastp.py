# -*- coding: utf-8 -*-
import argparse
from subprocess import Popen, PIPE, CalledProcessError
import os
import time
import re
import pickle
from parse_reference_files import ncbi_chrname_mapping,\
      readFasta, combine_indices,\
      interval_overlap, read_refseq_gff, read_ensembl_gff,\
      salmosalar_genomic_fna, rev_comp

 
base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23',\
                    help="")
parser.add_argument("-blastp", "--blastp", dest = "blastp", \
                    default = 'src/ncbi-blast-2.14.0+/bin/blastp', \
                    help="blastp path")
parser.add_argument("-run", "--run", dest = "run", \
                    default = 1, \
                    help="1: run blastp commands. 0: do not run blastp commands.")
args = parser.parse_args()

print(args.dataset)

script_start_time = time.time()

runblast = int(args.run)
blastp_path = os.path.join(base_dir, os.path.normpath(args.blastp))
if not os.path.isfile(blastp_path):
    err_msg = 'blastp path DNE'
    raise ValueError(err_msg)

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)    

ncbi_chrname_map = ncbi_chrname_mapping()
chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}


#reference files
reference_dir = os.path.join(base_dir , 'data', 'reference')
#RefSeq
refseq_dir = os.path.join(reference_dir, 'RefSeq')
chromosome_genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna') 
 
refseq_gff = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic_with_utrs.gff')
refseq_gff_results = read_refseq_gff(refseq_gff)
refseq_gff_entries = refseq_gff_results['refseq_gff_entries']
refseq_gff_start_end = refseq_gff_results['refseq_gff_start_end']

refseq_prot_fasta = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_protein.faa')
refseq_proteins = readFasta(refseq_prot_fasta)
 
#Ensembl
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')

ensembl_gff = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.99.chr.gff3')
ensembl_gff_results = read_ensembl_gff(ensembl_gff)
ensembl_gff_entries = ensembl_gff_results['ensembl_gff_entries']
ensembl_gff_start_end = ensembl_gff_results['ensembl_gff_start_end']

ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta) 

#Augustus
augustus_dir = os.path.join(data_dir, '3_Augustus')
hints_dir = os.path.join(augustus_dir, 'hints_dir')
augustus_outdir = os.path.join(augustus_dir, 'out_dir')   
augustus_genome_fna = os.path.join(augustus_dir,'genome_fna')

#output
blastp_dir = os.path.join(data_dir, '4_blastp')
blastp_in_dir = os.path.join(blastp_dir, 'input')
blastp_out_dir = os.path.join(blastp_dir, 'output')
blastp_cmds_dir = os.path.join(blastp_dir, 'cmds')
augustus_parsed_dir = os.path.join(augustus_dir, 'parsed_output')
Augustus_predicted_proteins = os.path.join(augustus_parsed_dir,\
                                           'Augustus_predicted_proteins.pickle')
for directory in [blastp_in_dir, blastp_out_dir,\
                  blastp_cmds_dir, augustus_parsed_dir]:
    if not os.path.isdir(directory):
        os.makedirs(directory)


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

def cmdrun(cmd_listform, outfn):
    start_time = time.time()
    process = Popen(cmd_listform, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    # check the output to determine if a silent error occurred
    if process.returncode == 0:
        if (has_bad_exception(stdout) or has_bad_exception(stderr)):
            raise CalledProcessError(process.returncode, cmd_listform, output=stdout)
        else:
            msg = "complete: "+str(time.time()-start_time)+"s"
            outfile = open(outfn, 'w')
            outfile.write(msg+'\n')
            outfile.write(' '.join(cmd_listform))
            outfile.close()
            return msg
    else:
        raise CalledProcessError(process.returncode, cmd_listform, output=stdout)

####################### Read Augustus output #######################
if not os.path.isfile(Augustus_predicted_proteins):
    print( 'Reading Augustus output...')
    nucleotides = set(['a','c','t','g'])
    amino_acids = set(['A', 'C', 'E', 'D', 'G',
                        'F', 'I', 'H', 'K', 'M',
                         'L', 'N', 'Q', 'P', 'S',
                          'R', 'T', 'W', 'V', 'Y','X'])
    aa_codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 
    'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W'} 
        
    def read_augustus_coor():
        fna_coordinates = {}
        for fn in os.listdir(augustus_genome_fna):
            from_loc = ''
            to_loc = ''
            seq = []
            with open(os.path.join(augustus_genome_fna, fn),'r') as infile:
                for line in infile:
                    if line[0]=='>':
                        #chrssa02_EG_38_Salmo_salar.fa
                        #>chrssa02 [42964064-43086419]         
                        coor = line.strip().split('[')[1].split(']')[0].split('-')
                        from_loc = int(coor[0])
                        to_loc = int(coor[1])
                    else:
                        seq.append(line.strip())
            fna_coordinates[fn.replace('_Salmo_salar.fa','')] = (from_loc, to_loc)
        return fna_coordinates
    
    augustus_log = open(os.path.join(augustus_parsed_dir, 'Parse_Augustus_output.log'),'w')   
    
    Augustus_genome_coordinates = read_augustus_coor()
    gff_tags = ['seqname',  'source',  'feature',  'start',  'end',\
                'score',  'strand',  'frame',  'transcript and gene name']
    gff_columns = {x:xi for xi,x in enumerate(gff_tags)}
    tag_list = ['coding sequence', 'protein sequence', 'Evidence for and against this transcript', \
                "% of transcript supported by hints (any source)",\
               "CDS exons", "CDS introns", "5'UTR exons and introns", "3'UTR exons and introns",\
               "hint groups fully obeyed", "incompatible hint groups"] 
    
    predicted_prots_event_hints = {}
    total_EG_peptides = set()
    total_known_peptides = set()
    augustus_log.write('###Reading Augustus Output###\n')
    for outfn in os.listdir(augustus_outdir):
        augustus_log.write('#'+outfn+'\n')
        EG = outfn.split('.augustus.hints.gff')[0]
        aug_genome_fna_coor = Augustus_genome_coordinates[EG]
        chromosome = outfn.split('_')[0]
        chrom_dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, chromosome)
        
        #read augustus hints
        aug_hints = {'enosi@':set(),'MSGF+@':{}}
        #augustus output doesn't tell you which hint exactly, just the source (P or M)
        with open(os.path.join(hints_dir, EG+'.hints.gff'),'r') as infile:
            for line in infile:
                sp = line.strip().split('\t')
                if 'enosi@' in sp[1]:
                    #enosi@chrssa02_EG_1/HQSEIRPG:EAEPQIMDYQTQQYK
                    EG_group,peptide = sp[1].split('@')[1].split('/')
                    aug_hints['enosi@'].add(peptide)
                elif 'MSGF+' in sp[1]:
                    #MSGF+@TNEETDDDGWTTVAR[XP_014055095.1|cds-XP_014055095.1:1243-1257]
                    known_pep = sp[1].split('MSGF+@')[1].split('[')[0]
                    if known_pep not in aug_hints['MSGF+@']:
                        aug_hints['MSGF+@'][known_pep] = set()
                    aug_hints['MSGF+@'][known_pep].add(sp[1].split('MSGF+@')[1])
    
        parse_start = False
        gff_list = {}
        #each predicted gene from augustus may have multiple transcripts and/or multiple protein sequences.
        with open(os.path.join(augustus_outdir, outfn),'r') as infile:
            for line in infile:
                if '# Predicted genes for sequence number' in line:
                    parse_start = True
                    continue
                if parse_start==True:
                    #each gene may have multiple transcripts
                    if '# start gene ' in line:
                        seq_info = []
                        indexes = {x:[] for x in tag_list}
                    elif '# end gene ' in line:                  
                        gene_id = ''
                        transcript_id = ''
                        strand = ''                    
                        for xi,x in enumerate(seq_info):
                            if x[0]=='#':
                                tag = x.lstrip('# ').split(':')[0].strip().split('=')[0].strip() 
                                if tag in indexes:
                                    indexes[tag].append(xi)
                                else:
                                    if set([ch for ch in tag])-nucleotides-amino_acids-set(['[',']','n']):
                                        print(tag)
                            else:
                                if x.strip().split()[0]==chromosome:
                                    sp = x.strip().split('\t')
                                    gff_line = {v:sp[gff_columns[v]] for v in gff_columns}
                                    strand = gff_line['strand']
                                    feature = gff_line['feature']
                                    attributes = gff_line['transcript and gene name']
                                    
                                    #readjust coordinates
                                    b = sorted([int(gff_line['start'])+aug_genome_fna_coor[0]-1, int(gff_line['end'])+aug_genome_fna_coor[0]-1]) #1-index
                                    gff_line['start'] = str(b[0])
                                    gff_line['end'] = str(b[1]) #exact end
                                    gff_line = '\t'.join([gff_line[v] for v in gff_tags])
                                    
                                    if feature == 'gene':
                                        gene_id = attributes.strip()
                                        if gene_id not in gff_list:
                                            gff_list[gene_id] = {'i':xi,
                                                                 'main':gff_line,
                                                                 'transcripts':{}}
                                    elif feature == 'transcript':
                                        gene_id = attributes.split('.')[0]
                                        transcript_id = attributes.split('.')[1]
                                        if transcript_id not in gff_list[gene_id]['transcripts']:
                                            gff_list[gene_id]['transcripts'][transcript_id] = {'i':xi, 'strand':strand,
                                                                                    'main':gff_line,'CDS':[],
                                                                                    'start_codon':[],'stop_codon':[],
                                                                                    'cds_genomic_coor':[]}
                                    elif feature in ['CDS', 'start_codon', 'stop_codon']:
                                        #transcript_id "g7.t1"; gene_id "g7";
                                        gene_id = attributes.split('transcript_id "')[1].split('"')[0].split('.')[0]
                                        transcript_id = attributes.split('transcript_id "')[1].split('"')[0].split('.')[1]
                                        gff_list[gene_id]['transcripts'][transcript_id]['cds_genomic_coor'].append(b)
                                        if feature == 'CDS':
                                            gff_list[gene_id]['transcripts'][transcript_id][feature].append(gff_line)
                                        else:
                                            gff_list[gene_id]['transcripts'][transcript_id][feature] = {'main':gff_line,'i':xi}
                        
                        for transcript_id in gff_list[gene_id]['transcripts']:                        
                            strand = gff_list[gene_id]['transcripts'][transcript_id]['strand']
                            i = int(transcript_id.replace('t',''))-1
                            final_end = len(seq_info)
                            next_t = 't'+str(int(transcript_id.replace('t',''))+1)
                            if next_t in gff_list[gene_id]['transcripts']:
                                final_end = gff_list[gene_id]['transcripts'][next_t]['i']
                            
                            coding_sequence = seq_info[indexes['coding sequence'][i]:indexes['protein sequence'][i]]
                            coding_sequence = ''.join([x.lstrip('# ') for x in coding_sequence]).split('[')[1].split(']')[0].upper()
                            protein_sequence = seq_info[indexes['protein sequence'][i]:indexes['Evidence for and against this transcript'][i]]
                            protein_sequence = ''.join([x.lstrip('# ') for x in protein_sequence]).split('[')[1].split(']')[0].upper()                 
                            
                            if not protein_sequence:
                                augustus_log.write('\t'+gene_id+'.'+transcript_id+\
                                                   '\t'+'no protein sequence!'+'\n')
                                continue
                                                       
                            CDS_exons = [x.lstrip('# ') for x\
                                         in seq_info[indexes['CDS exons'][i]:indexes['CDS introns'][i]]]
                            CDS_exons = {x.split(':')[0].strip():x.split(':')[1].strip() for x in CDS_exons}
                            
                            five_prime_UTR = [x.lstrip('# ') for x\
                                              in seq_info[indexes["5'UTR exons and introns"][i]:indexes["3'UTR exons and introns"][i]]]
                            five_prime_UTR = {x.split(':')[0].strip():x.split(':')[1].strip()\
                                              for x in five_prime_UTR}
                            
                            three_prime_UTR = [x.lstrip('# ') for x\
                                               in seq_info[indexes["3'UTR exons and introns"][i]:indexes["hint groups fully obeyed"][i]]]
                            three_prime_UTR = {x.split(':')[0].strip():x.split(':')[1].strip()\
                                               for x in three_prime_UTR}
                            
                            hint_groups_fully_obeyed = [x.lstrip('# ') for x\
                                                        in seq_info[indexes['hint groups fully obeyed'][i]:indexes['incompatible hint groups'][i]]]
                            hint_groups_fully_obeyed = {x.split(':')[0].strip():x.split(':')[1].strip()\
                                                        for x in hint_groups_fully_obeyed}
                            
                            hint_groups_obeyed = int(hint_groups_fully_obeyed['hint groups fully obeyed'])
                            if hint_groups_obeyed==0:
                                augustus_log.write('\t'+gene_id+'.'+transcript_id+'\t'+'hint groups fully obeyed=0'+'\n')
                                continue                            
                            
                            ##skip prediction if no 'P' hints
                            P_hint = False   
                            if 'P' in CDS_exons:
                                P_hint = True
                            if 'P' in five_prime_UTR:
                                P_hint = True
                            if 'P' in three_prime_UTR:
                                P_hint = True
                            if 'P' in hint_groups_fully_obeyed:
                                P_hint = True
                                
                            if P_hint ==False:
                                augustus_log.write('\t'+gene_id+'.'+transcript_id+\
                                                   '\t'+'P_hint ==False'+'\n')
                                continue
                            
                            ##skip prediction if no enosi hints
                            found_peps = set()
                            for pep in aug_hints['enosi@']:
                                #find entire peptide in protein sequence
                                index = [str(m.start()) for m\
                                         in re.finditer(pep.replace(':',''), protein_sequence)]
                                if index:
                                    total_EG_peptides.add(pep)
                                    found_peps.add(pep+'['+','.join(index)+']')
    
                            if not found_peps:
                                augustus_log.write('\t'+gene_id+'.'+\
                                                transcript_id+'\t'+\
                                                'no event peptide in pred prot seq'+'\n')
                                continue
                                
                            ##check readjusted coordinates
                            cds_genomic_coor = gff_list[gene_id]['transcripts'][transcript_id]['cds_genomic_coor']
                            cds_genomic_coor = combine_indices(cds_genomic_coor)
                            if strand == '-':
                                cds_genomic_coor = cds_genomic_coor[::-1]
                            
                            Aug_CDS_coor = []    
                            cds_join = []
                            for b in cds_genomic_coor:
                                subseq = chrom_dna_seq[b[0]-1:b[1]]
                                if strand == '-':
                                    cds_join.append(rev_comp(subseq))
                                    Aug_CDS_coor.extend(range(b[0],b[1]+1)[::-1]) # 1-index
                                else:
                                    cds_join.append(subseq)
                                    Aug_CDS_coor.extend(range(b[0],b[1]+1))
    
                            if coding_sequence!=''.join(cds_join).upper():
                                raise ValueError('coding_sequence!=cds_join')
    
                            #1-index coor    
                            aug_prot_aa_genomic_coor = {}
                            start_cdi = 0
                            for ai,aa in enumerate(protein_sequence):
                                expected_aa = ''
                                for cds_i, cds in enumerate(coding_sequence):
                                    if cds_i<start_cdi:
                                        continue
                                    subseq = coding_sequence[cds_i:cds_i+3]
                                    expected_aa = aa_codon_table[subseq]
                                    if expected_aa==aa:
                                        break
                                if not expected_aa:
                                    raise ValueError('expected_aa empty')
                                aug_prot_aa_genomic_coor[ai] = {'aa':aa,
                                                                'aa_coor':Aug_CDS_coor[cds_i:cds_i+3]}
                                start_cdi = cds_i+3
                                    
                            if len(aug_prot_aa_genomic_coor)!=len(protein_sequence):
                                raise ValueError('')
                            
                            #known peptide in augustus prot seq
                            found_known = set()
                            for known_pep in aug_hints['MSGF+@']:
                                index = [str(m.start()) for m in re.finditer(known_pep.replace('L','I'), protein_sequence.replace('L','I'))]
                                if index:
                                    total_known_peptides.add(known_pep)
                                    found_known.add(known_pep+'['+','.join(index)+']')                         
                            
                            gff_lines = [gff_list[gene_id]['main']]
                            gff_lines.append(gff_list[gene_id]['transcripts'][transcript_id]['main'])
                            if strand == '-':
                                gff_lines.append(gff_list[gene_id]['transcripts'][transcript_id]['stop_codon']['main'])
                                gff_lines.extend(gff_list[gene_id]['transcripts'][transcript_id]['CDS'])
                                gff_lines.append(gff_list[gene_id]['transcripts'][transcript_id]['start_codon']['main'])
                            if strand == '+':
                                gff_lines.append(gff_list[gene_id]['transcripts'][transcript_id]['start_codon']['main'])
                                gff_lines.extend(gff_list[gene_id]['transcripts'][transcript_id]['CDS'])
                                gff_lines.append(gff_list[gene_id]['transcripts'][transcript_id]['stop_codon']['main'])
                            
                            evidence_lines = seq_info[indexes['Evidence for and against this transcript'][i]:final_end]
                            if '# end' in evidence_lines[-1]:
                                evidence_lines = evidence_lines[:-1]
    
                            if protein_sequence not in predicted_prots_event_hints:
                                predicted_prots_event_hints[protein_sequence] = {}
                            if chromosome not in predicted_prots_event_hints[protein_sequence]:
                                predicted_prots_event_hints[protein_sequence][chromosome] = {}                        
                            if EG not in predicted_prots_event_hints[protein_sequence][chromosome]:
                                predicted_prots_event_hints[protein_sequence][chromosome][EG] = []
                            predicted_prots_event_hints[protein_sequence][chromosome][EG].append(
                                                       {'event_peptides':found_peps,'known_peptides':found_known,\
                                                       'gene_id':gene_id, 'transcript_id':transcript_id,\
                                                       'evidence_lines':evidence_lines,\
                                                       'coding_sequence':coding_sequence,\
                                                       'gff_lines':gff_lines, 'aug_cds_s_e':cds_genomic_coor,\
                                                       'aug_prot_aa_genomic_coor':aug_prot_aa_genomic_coor})                        
                        seq_info = []
                    
                    else:
                        #lines in between # start gene and # end gene 
                        seq_info.append(line.strip())
    
    print ('Total Augustus predicted proteins: ', len(predicted_prots_event_hints))
    augustus_log.write('#Total Augustus predicted proteins: '+str(len(predicted_prots_event_hints))+'\n')
    augustus_log.write('###FINISHED Reading Augustus Output###\n')
    augustus_log.close()
    
    with open(Augustus_predicted_proteins, 'wb') as handle:
        pickle.dump(predicted_prots_event_hints , handle, protocol=pickle.HIGHEST_PROTOCOL)
    
else:
    predicted_prots_event_hints = pickle.load(open(Augustus_predicted_proteins,'rb'))
    
####################### blastp #######################
print ('blastp...')
pred_num = 0
protein_sequences = sorted(predicted_prots_event_hints.keys())
for protein_sequence in protein_sequences:
    pred_num+=1   
    blastp_infile = os.path.join(blastp_in_dir, 'ProtSeq'+str(pred_num)+'.blastp.txt')
                          
    with open(blastp_infile,'w') as outfile:
        split_60 = [protein_sequence[j:j+60] for j in range(0,len(protein_sequence),60)]
        outfile.write('>'+'ProtSeq'+str(pred_num)+'\n'+'\n'.join(split_60)+'\n')

for outfn in os.listdir(blastp_in_dir):
    blastp_infile = os.path.join(blastp_in_dir, outfn)
    blastp_output = os.path.join(blastp_out_dir, outfn)
    blastp_logfile = os.path.join(blastp_out_dir, outfn.replace('.blastp.txt','.blastp.log'))
    cmd = [blastp_path,
           "-query",
           blastp_infile,
           "-db",
           "Actinopterygii.fsa",
           "-evalue",
           str(1e-10),
           "-word_size",
           str(5),
           "-outfmt",
           "5", #BLAST 5 XML
           "-out",
           blastp_output]
    
    with open(os.path.join(blastp_cmds_dir, outfn.split('.')[0]+'_cmds.sh'),'w') as outfile:
        outfile.write(' '.join(cmd))

    if runblast==1:
        print(outfn)
        result = cmdrun(cmd, blastp_logfile)
        print(result)

#########################

##RefSeq Splice db search output
RefSeq_Splice_dir = os.path.join(data_dir,'1_RefSeq_SpliceDB_Search')
combined_enosi_output = os.path.join(RefSeq_Splice_dir, 'combined_Enosi_Output')
event_fn = os.path.join(combined_enosi_output, 'event_parsed.txt')
agg_events = os.path.join(combined_enosi_output, 'all_events.pickle')

#RefSeq Ensembl db search output
RefSeq_Ensembl_dir = os.path.join(data_dir,'2_RefSeq_Ensembl_Search')
genomic_coor_output_dir = os.path.join(RefSeq_Ensembl_dir, 'genomic_coordinates')
Peptide_genomic_coordinates = os.path.join(genomic_coor_output_dir,\
                                'RefSeq_Ensembl_Peptide_genomic_coordinates.txt')

#Augustus
augustus_dir = os.path.join(data_dir, '3_Augustus')
augustus_parsed_dir = os.path.join(augustus_dir, 'parsed_output')
Augustus_predicted_proteins = os.path.join(augustus_parsed_dir,\
                                           'Augustus_predicted_proteins.pickle')

augustus_ref_overlap = os.path.join(augustus_parsed_dir,\
                                    'Augustus_reference_overlap.pickle') 
    
#blastp
blastp_dir = os.path.join(data_dir, '4_blastp')
blastp_in_dir = os.path.join(blastp_dir, 'input')

#jbrowse gff
jbrowse_gff_dir = os.path.join(augustus_dir, 'jbrowse_gff')
if not os.path.isdir(jbrowse_gff_dir):
    os.mkdir(jbrowse_gff_dir)


def augustus_reference_gene_overlap(reference_se, reference_gff, reference_name,\
                                    augustus_prot_seq, augustus_gene_SE, augustus_cds_SE):
    #overlapping gene
    overlapping_gene = interval_overlap((augustus_gene_SE[0], augustus_gene_SE[1]),\
                                        sorted([sorted((int(reference_se[key]['start']), int(reference_se[key]['end']))) for key in reference_se]))
    
    ref_gene_list = {}
    for rg_s_e in overlapping_gene:
        if 'ID' not in reference_se[rg_s_e]['attributes']:
            print(reference_name, reference_se[rg_s_e]['attributes'])
        ref_gene_id = reference_se[rg_s_e]['attributes']['ID']
        
        #for transcript belonging to gene overlapping augustus gene
        for rt in reference_gff[ref_gene_id]:
            if rt=='other':
                continue       
            
            #CDS
            CDS_entries = {}
            CDS_se_list = []
            for cds in reference_gff[ref_gene_id][rt]['CDS']:
                cds_s_e = sorted([int(cds['start']), int(cds['end'])])
                if cds_s_e not in CDS_se_list:
                    CDS_se_list.append(cds_s_e)
                cds_s_e = str(cds_s_e[0])+'-'+str(cds_s_e[1])
                CDS_entries[cds_s_e] = cds

            #5UTR
            fiveUTR_entries = {}
            fiveUTR_se_list = []
            for utr in reference_gff[ref_gene_id][rt]['5UTR']:
                utr_s_e = sorted([int(utr['start']), int(utr['end'])])
                if utr_s_e not in fiveUTR_se_list:
                    fiveUTR_se_list.append(utr_s_e)
                utr_s_e = str(utr_s_e[0])+'-'+str(utr_s_e[1])
                fiveUTR_entries[utr_s_e] = utr
            
            #3UTR
            threeUTR_entries = {}
            threeUTR_se_list = []
            for utr in reference_gff[ref_gene_id][rt]['3UTR']:
                utr_s_e = sorted([int(utr['start']), int(utr['end'])])
                if utr_s_e not in threeUTR_se_list:
                    threeUTR_se_list.append(utr_s_e)
                utr_s_e = str(utr_s_e[0])+'-'+str(utr_s_e[1])
                threeUTR_entries[utr_s_e] = utr

            # reference CDS,UTR and augustus CDS overlap
            for aug_cds_i, aug_cds in enumerate(augustus_cds_SE):
                overlapping_CDS = [CDS_entries[x] for x\
                                   in interval_overlap((aug_cds[0], aug_cds[1]), sorted(CDS_se_list))]
                overlapping_5UTR = [fiveUTR_entries[x] for x\
                                    in interval_overlap((aug_cds[0], aug_cds[1]), sorted(fiveUTR_se_list))]
                overlapping_3UTR = [threeUTR_entries[x] for x\
                                    in interval_overlap((aug_cds[0], aug_cds[1]), sorted(threeUTR_se_list))]

                if not (overlapping_CDS or overlapping_5UTR or overlapping_3UTR):
                    continue
                    
                if ref_gene_id not in ref_gene_list:
                    ref_gene_list[ref_gene_id] = {}
                if rt not in ref_gene_list[ref_gene_id]:
                    ref_gene_list[ref_gene_id][rt] = {}
            
                for CDS in overlapping_CDS:
                    protein_id = CDS['protein_id']
                    ref_Id = CDS['attributes']['ID']
                    ref_s_e = sorted([int(CDS['start']), int(CDS['end'])])
                    if 'CDS' not in ref_gene_list[ref_gene_id][rt]:
                        ref_gene_list[ref_gene_id][rt]['CDS'] = []
                    protseq= ''
                    if reference_name == 'ensembl':
                        protseq = ensembl_proteins[protein_id+'.1']
                    else:
                        protseq = refseq_proteins[protein_id]
                    
                    if protseq==augustus_prot_seq:
                        continue
                    
                    ref_gene_list[ref_gene_id][rt]['CDS'].append({'protein_id':protein_id, 'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})

                for utr in overlapping_5UTR:
                    ref_Id = '5utr'
                    if 'ID' in utr['attributes']:
                        ref_Id = utr['attributes']['ID']                   
                    ref_s_e = sorted([int(utr['start']), int(utr['end'])])
                    if '5UTR' not in ref_gene_list[ref_gene_id][rt]:
                        ref_gene_list[ref_gene_id][rt]['5UTR'] = []
                    ref_gene_list[ref_gene_id][rt]['5UTR'].append({'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})
                    
                for utr in overlapping_3UTR:
                    ref_Id = '3utr'
                    if 'ID' in utr['attributes']:
                        ref_Id = utr['attributes']['ID']
                    ref_s_e = sorted([int(utr['start']), int(utr['end'])])
                    if '3UTR' not in ref_gene_list[ref_gene_id][rt]:
                        ref_gene_list[ref_gene_id][rt]['3UTR'] = []
                    ref_gene_list[ref_gene_id][rt]['3UTR'].append({'augustus':('cds'+str(aug_cds_i), str(aug_cds[0]), str(aug_cds[1])),'ref':(ref_Id, ref_s_e[0], ref_s_e[1])})

    return ref_gene_list
 

####################### Enosi events ############################
print( 'reading enosi events...')
enosi_events = pickle.load(open(agg_events,'rb'))

Event_Peptide_Coor = {}
with open(event_fn,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            col = {x:xi for xi,x in enumerate(sp)}
        else:
            Event_group = sp[col['#EventGroup']]
            EventNum = sp[col['EventNum']]
            chromosome = sp[col['Chromosome']]
            Ensembl_overlap = {}
            event_type = enosi_events[chromosome][EventNum]['event_type']
            for rg in enosi_events[chromosome][EventNum]['Ensembl_overlap']:
                for rt in enosi_events[chromosome][EventNum]['Ensembl_overlap'][rg]:
                    for p in enosi_events[chromosome][EventNum]['Ensembl_overlap'][rg][rt]:
                        if p not in Ensembl_overlap:
                            Ensembl_overlap[p] = {}
                        if rg not in Ensembl_overlap[p]:
                            Ensembl_overlap[p][rg] = set()
                        Ensembl_overlap[p][rg].add(rt.replace('transcript:',''))
            Ensembl_overlap_str = {}
            for p in Ensembl_overlap:
                Ensembl_overlap_str[p] = '/'.join([rg+'|transcript:'+','.join(Ensembl_overlap[p][rg]) for rg in Ensembl_overlap[p]])
            
            if Event_group not in Event_Peptide_Coor:
                Event_Peptide_Coor[Event_group] = {}
            strand = sp[col['Strand']]
            Peptides = sp[col['Peptides']].split('|')
            for entry in Peptides:              
                pep, loc = entry.split('/') #1 index
                loc_coor = []
                for coordinates in loc.split(';'):
                    coordinates = coordinates.split('-')
                    start = int(coordinates[0])
                    stop = int(coordinates[1])
                    loc_coor.append((start, stop))
                if strand=='-':
                    loc_coor = loc_coor[::-1]
                if pep not in Event_Peptide_Coor[Event_group]:
                    Event_Peptide_Coor[Event_group][pep] = []
                Ensembl_overlap = []
                if pep in Ensembl_overlap_str:
                    Ensembl_overlap = Ensembl_overlap_str[pep]
                Event_Peptide_Coor[Event_group][pep].append({'chromosome':chromosome,'strand':strand, 
                                                          'EventNum':EventNum, 'genomic_coordinates':loc_coor,
                                                          'Ensembl_overlap':Ensembl_overlap, 'event_type':event_type})

print('reading known peptide genomic coordinates...')
MSGF_peptides_gff = open(os.path.join(jbrowse_gff_dir, 'MSGF+_peptides.gff3'),'w')
MSGF_peptides_gff.write('##gff-version 3'+'\n')
##Chromosome	strand	genomic_coordinate_intervals	peptide	protein	protein_coordinates
#chrssa06	+	67189552-67189581	LSSTSGPGVK	XP_014061050.1|cds-XP_014061050.1	14-23
known_peptide_gcoor = {}
with open(Peptide_genomic_coordinates,'r') as infile:
    for line in infile:
        sp = line.strip().split('\t')
        if line[0]=='#':
            col = {x:xi for xi,x in enumerate(sp)}
        else:
            Chromosome = sp[col['#Chromosome']]
            nc_chrom = chrname_ncbi_map[Chromosome]
            strand = sp[col['strand']]
            genomic_coordinate_intervals = sp[col['genomic_coordinate_intervals']].replace(';',',')
            peptide = sp[col['peptide']]
            protein = sp[col['protein']].split('|')[0]
            protein_coordinates = sp[col['protein_coordinates']]
            if Chromosome not in known_peptide_gcoor:
                known_peptide_gcoor[Chromosome] = {}
            if peptide not in known_peptide_gcoor[Chromosome]:
                known_peptide_gcoor[Chromosome][peptide] = {}
            if protein not in known_peptide_gcoor[Chromosome][peptide]:
                known_peptide_gcoor[Chromosome][peptide][protein] = []
            location = {'strand':strand,'genomic_coordinate_intervals':genomic_coordinate_intervals}
            known_peptide_gcoor[Chromosome][peptide][protein].append(location)
            
            for loc in sp[col['genomic_coordinate_intervals']].split(';'):
                coor = loc.split('-')
                attributes = ['ID='+peptide+'['+protein+']','Name='+peptide+'['+protein+']']
                MSGF_peptides_gff.write('\t'.join([nc_chrom, 'MSGF+', 'CDSpart', str(coor[0]), str(coor[1]), '.', strand, '.', ';'.join(attributes)])+'\n')
MSGF_peptides_gff.close()

######################################################
print ('finding related genes and proteins...')
predicted_prots_event_hints = pickle.load(open(Augustus_predicted_proteins,'rb'))

Enosi_events_gff = open(os.path.join(jbrowse_gff_dir, 'Enosi_events.gff3'),'w')
Enosi_events_gff.write('##gff-version 3'+'\n')
Augustus_Genes_gff = open(os.path.join(jbrowse_gff_dir, 'Augustus_Genes.gff3'),'w')
Augustus_Genes_gff.write('##gff-version 3'+'\n')
Augustus_transcript_cds_gff = open(os.path.join(jbrowse_gff_dir, 'Augustus_CDS.gff3'),'w')
Augustus_transcript_cds_gff.write('##gff-version 3'+'\n')

augustus_pred_prot_summary = open(os.path.join(augustus_parsed_dir,\
                                               'Augustus_pred_prot_summary.txt'),'w')
augustus_pred_prot_summary.write('\t'.join(['#Augustus_Predicted_Prot',\
                                            'Chromosome',\
                                            'jbrowse_coor',\
                                            'Augustus_gene_id', 'Augustus_transcript_id',\
                                            'Augustus_transcript_cds',\
                                            'Related_RefSeq_Gene|Transcript|Protein',\
                                            'Related_Ensembl_Gene|Transcript|Protein',\
                                            'Event_Peptides', 'Event_Peptides_Prot_Coor',\
                                            'Known_Peptides', 'Known_Peptides_Prot_Coor',\
                                            'Event_Supported_Ensembl_Transcripts',\
                                            'Predicted_Prot_Seq'])+'\n')

EG_group_count = {}
Entries = {}
for fn in os.listdir(blastp_in_dir):
    predprotid = fn.split('.blastp.txt')[0]
    
    aug_prot_seq = []
    with open(os.path.join(blastp_in_dir, fn),'r') as infile:
        for line in infile:
            if line[0]!='>':
                aug_prot_seq.append(line.strip())
    aug_prot_seq = ''.join(aug_prot_seq)

    Entries[predprotid] = {'aug_prot_seq':aug_prot_seq, 'predictions':{}}

    augustus_entries = predicted_prots_event_hints[aug_prot_seq]
    
    c=0
    for chromosome in augustus_entries:
        nc_chrom = chrname_ncbi_map[chromosome]
        for EG in augustus_entries[chromosome]:
            if EG not in EG_group_count:
                EG_group_count[EG] = set()
            EG_event_peps = Event_Peptide_Coor[EG]
            for prediction in augustus_entries[chromosome][EG]:
                c+=1                
                jbrowse_coor = ''
                gene_id = prediction['gene_id']
                transcript_id = prediction['transcript_id']
                
                key = predprotid+'|'+EG                
                EG_group_count[EG].add(predprotid+'.'+gene_id+'.'+transcript_id)                
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id] = {}
                
                #writing gff files to import into jbrowse 
                Augustus_Genes_gff.write('#Augustus Prediction:'+key+'|'+gene_id+'.'+transcript_id+'\n')
                Augustus_transcript_cds_gff.write('#Augustus Prediction:'+key+'|'+gene_id+'.'+transcript_id+'\n')
                                    
                overlap_refseq_genes = {}
                overlap_ensembl_genes = {}
                aug_gene_s_e = []
                aug_transcript_s_e = []
                aug_cds_s_e = []
                aug_strand= ''
                for line in prediction['gff_lines']:
                    sp = line.split('\t')
                    gff_line = {x:sp[xi] for xi,x in enumerate(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript and gene name'])}
                    if gff_line['feature'] == 'gene':
                        aug_gene_s_e = sorted([int(gff_line['start']), int(gff_line['end'])])
                        jbrowse_coor = nc_chrom+':'+str(aug_gene_s_e[0])+'...'+str(aug_gene_s_e[1])
                        attributes = ';'.join(['ID=gene_'+key+'|'+gene_id, 'Name=Gene_'+key+'|'+gene_id, 'gbkey=Gene', 'gene='+key+'|'+gene_id, 'gene_biotype=protein_coding'])
                        Augustus_Genes_gff.write('\t'.join([nc_chrom, gff_line['source'], gff_line['feature'], gff_line['start'], gff_line['end'], gff_line['score'], gff_line['strand'], gff_line['frame'], attributes])+'\n')                        
                        
                    elif gff_line['feature'] == 'transcript':
                        aug_strand = gff_line['strand']
                        aug_transcript_s_e = sorted([int(gff_line['start']), int(gff_line['end'])])
                        attributes = ';'.join(['ID=rna_'+key+'|'+gene_id+'.'+transcript_id, 'Name=Transcript_'+key+'|'+gene_id+'.'+transcript_id, 'gbkey=mRNA', 'gene='+key+'|'+gene_id])
                        Augustus_transcript_cds_gff.write('\t'.join([nc_chrom, gff_line['source'], gff_line['feature'], gff_line['start'], gff_line['end'], gff_line['score'], gff_line['strand'], gff_line['frame'], attributes])+'\n')
                    
                    elif gff_line['feature'] in ['CDS','start_codon','stop_codon']:
                        aug_cds_s_e.append(sorted([int(gff_line['start']), int(gff_line['end'])]))
                        attributes = ';'.join(['ID='+gff_line['feature']+'_'+key+'|'+gene_id+'.'+transcript_id, 'Parent=rna_'+key+'|'+gene_id+'.'+transcript_id, 'Name='+key, 'gbkey=CDS', 'gene='+key+'|'+gene_id])
                        Augustus_transcript_cds_gff.write('\t'.join([nc_chrom, gff_line['source'], 'CDS', gff_line['start'], gff_line['end'], gff_line['score'], gff_line['strand'], gff_line['frame'], attributes])+'\n')
                
                aug_cds_s_e = combine_indices(aug_cds_s_e)             
                if aug_strand == '-':
                    aug_cds_s_e = aug_cds_s_e[::-1]
          
                if aug_cds_s_e!= prediction['aug_cds_s_e']:
                    raise ValueError('')

                #write augustus cds gcoor to file
                Augustus_Transcript_CDS = [str(x[0])+'-'+str(x[1]) for x in aug_cds_s_e]
               
                ensembl_overlap = augustus_reference_gene_overlap(ensembl_gff_start_end[chromosome][aug_strand], ensembl_gff_entries[chromosome][aug_strand], 'ensembl', aug_prot_seq, aug_gene_s_e, aug_cds_s_e)
                refseq_overlap = augustus_reference_gene_overlap(refseq_gff_start_end[chromosome][aug_strand], refseq_gff_entries[chromosome][aug_strand], 'refseq', aug_prot_seq, aug_gene_s_e, aug_cds_s_e)
                
                related_genes_Ensembl = set()
                related_plist = set()
                for rg in ensembl_overlap:
                    #gene:ENSSSAG00000051046
                    gene_name = rg.replace('gene:','')
                    for rt in ensembl_overlap[rg]:
                        #transcript:ENSSSAT00000082184                        
                        transcript_name = rt.replace('transcript:','')  
                        if 'CDS' in ensembl_overlap[rg][rt]:
                            protein_ids = set()
                            for cds in ensembl_overlap[rg][rt]['CDS']:
                                protein_ids.add(cds['protein_id']+'.1')
                                augustus = cds['augustus']
                                reference = cds['ref']
                            if len(protein_ids)>1:
                                print( protein_ids)
                            rp = list(protein_ids)[0]
                            related_genes_Ensembl.add(gene_name+'|'+transcript_name+'|'+rp)
                            related_plist.add(rp)
                        
                        if '3UTR' in ensembl_overlap[rg][rt]:
                            for utr in ensembl_overlap[rg][rt]['3UTR']:
                                augustus = utr['augustus']
                                reference = utr['ref']
                                related_genes_Ensembl.add(gene_name+'|'+transcript_name+'|'+'3UTR')

                        if '5UTR' in ensembl_overlap[rg][rt]:
                            for utr in ensembl_overlap[rg][rt]['5UTR']:
                                augustus = utr['augustus']
                                reference = utr['ref']
                                related_genes_Ensembl.add(gene_name+'|'+transcript_name+'|'+'5UTR')

                
                related_genes_RefSeq = set()
                for rg in refseq_overlap:
                    gene_name = rg.replace('gene-','')
                    for rt in refseq_overlap[rg]:
                        transcript_name = rt.replace('rna-','')
                        if 'CDS' in refseq_overlap[rg][rt]:
                            protein_ids = set()
                            for cds in refseq_overlap[rg][rt]['CDS']:
                                protein_ids.add(cds['protein_id'])
                                augustus = cds['augustus']
                                reference = cds['ref']
                            if len(protein_ids)>1:
                                print (protein_ids)
                            rp = list(protein_ids)[0]
                            related_genes_RefSeq.add(gene_name+'|'+transcript_name+'|'+rp)
                            related_plist.add(rp)
                        
                        if '3UTR' in refseq_overlap[rg][rt]:
                            for utr in refseq_overlap[rg][rt]['3UTR']:
                                augustus = utr['augustus']
                                reference = utr['ref']
                                related_genes_RefSeq.add(gene_name+'|'+transcript_name+'|'+'3UTR')

                        if '5UTR' in refseq_overlap[rg][rt]:
                            for utr in refseq_overlap[rg][rt]['5UTR']:
                                augustus = utr['augustus']
                                reference = utr['ref']
                                related_genes_RefSeq.add(gene_name+'|'+transcript_name+'|'+'5UTR')

                #event peps
                event_peps_genomic_coordinates = set()
                event_peps_protein_coordinates = prediction['event_peptides']
                event_supported_ensembl_transcripts = set()
                event_peptides = [x.split('[')[0] for x in event_peps_protein_coordinates]
                for pep in event_peptides:                    
                    for peploc in EG_event_peps[pep]:
                        ep_eventnum = peploc['EventNum']
                        ep_supported_ensembl = peploc['Ensembl_overlap']
                        ep_genomic_coor = sorted(peploc['genomic_coordinates'])
                        min_s = min([x[0] for x in ep_genomic_coor])
                        max_e = max([x[1] for x in ep_genomic_coor])
                        ep_genomic_coor = ','.join([str(x[0])+'-'+str(x[1]) for x in ep_genomic_coor])
                        if not (aug_gene_s_e[0]<min_s<aug_gene_s_e[1] and aug_gene_s_e[0]<max_e<aug_gene_s_e[1]):
                            ep_genomic_coor = 'UNDEFINED'
                            #raise ValueError('event peptide genomic coor not in augustus gene range')     
                        event_peps_genomic_coordinates.add(pep+'|'+peploc['chromosome']+'|'+peploc['strand']+'|'+ep_eventnum+'|'+ep_genomic_coor)
                        if ep_genomic_coor != 'UNDEFINED':
                            for pep_part_i,pep_part in enumerate(pep.split(':')):
                                coor = peploc['genomic_coordinates'][pep_part_i]
                                attributes = ['ID='+EG+'/'+ep_eventnum+'/'+pep_part, 'Name='+pep_part+' '+'['+'event_num:'+ep_eventnum+'/'+peploc['event_type']+'/'+pep+']']
                                if ep_supported_ensembl:
                                    attributes.append('Note=event supports Ensembl annotation')
                                #NC_027306.1	Enosi	CDSpart	1396047	1396074	.	+	.	ID=EG_3/EIAPAGTGMA;Name=EIAPAGTGMA [novel gene/EIAPAGTGMA:PGLMGIVGGR]
                                Enosi_events_gff.write('\t'.join([nc_chrom, 'EnoSi', 'CDSpart', str(coor[0]), str(coor[1]), '.', peploc['strand'], '.', ';'.join(attributes)])+'\n')
                        if ep_supported_ensembl:
                            event_supported_ensembl_transcripts.add(pep+'|'+ep_supported_ensembl)
                
                if not event_peps_genomic_coordinates:
                    raise ValueError('')
                
                #known peps
                known_peps_genomic_coordinates = set()
                known_peps_protein_coordinates = prediction['known_peptides']
                known_peptides = [x.split('[')[0] for x in known_peps_protein_coordinates]
                for kp in known_peptides:
                    for prot in known_peptide_gcoor[chromosome][kp]:
                        if prot not in related_plist:
                            continue
                        for loc in known_peptide_gcoor[chromosome][kp][prot]:
                            known_peps_genomic_coordinates.add(kp+'|'+chromosome+'|'+loc['strand']+'|'+loc['genomic_coordinate_intervals'])

                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['aug_strand'] = aug_strand
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['coding_sequence'] = prediction['coding_sequence']
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['aug_prot_aa_genomic_coor'] = prediction['aug_prot_aa_genomic_coor']
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['aug_cds_s_e'] = aug_cds_s_e
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['aug_gene_s_e'] = aug_gene_s_e                
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['ensembl_overlap'] = ensembl_overlap
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['refseq_overlap'] = refseq_overlap
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['related_genes_RefSeq'] = related_genes_RefSeq
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['related_genes_Ensembl'] = related_genes_Ensembl
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['event_peps_genomic_coordinates'] = event_peps_genomic_coordinates
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['event_peps_protein_coordinates'] = event_peps_protein_coordinates    
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['event_supported_ensembl_transcripts'] = event_supported_ensembl_transcripts   
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['known_peps_genomic_coordinates'] = known_peps_genomic_coordinates
                Entries[predprotid]['predictions'][EG+'.'+gene_id+'.'+transcript_id]['known_peps_protein_coordinates'] = known_peps_protein_coordinates
                
                augustus_pred_prot_summary.write('\t'.join([key, chromosome, jbrowse_coor,\
                                                            gene_id+'|'+str(aug_gene_s_e[0])+'-'+str(aug_gene_s_e[1]),\
                                                            transcript_id+'|'+str(aug_transcript_s_e[0])+'-'+str(aug_transcript_s_e[1])+'|'+aug_strand,\
                                                            ';'.join(Augustus_Transcript_CDS),\
                                                            ';'.join(related_genes_RefSeq),\
                                                            ';'.join(related_genes_Ensembl),\
                                                            ';'.join(event_peps_genomic_coordinates),\
                                                            ';'.join(event_peps_protein_coordinates),
                                                            ';'.join(known_peps_genomic_coordinates),\
                                                            ';'.join(known_peps_protein_coordinates),\
                                                            ';'.join(event_supported_ensembl_transcripts),\
                                                            aug_prot_seq])+'\n')

with open(augustus_ref_overlap, 'wb') as handle:
    pickle.dump(Entries , handle, protocol=pickle.HIGHEST_PROTOCOL)
    
augustus_pred_prot_summary.close()
Augustus_transcript_cds_gff.close()
Augustus_Genes_gff.close()
Enosi_events_gff.close()

print('\n\nComplete ('+str(time.time()-script_start_time)+' seconds)')
