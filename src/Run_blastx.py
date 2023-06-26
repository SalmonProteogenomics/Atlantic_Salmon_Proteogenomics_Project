# -*- coding: utf-8 -*-
import os
import math
import argparse
from parse_reference_files import readFasta, salmosalar_genomic_fna
from subprocess import Popen, PIPE, CalledProcessError
import time

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23', help="")
parser.add_argument("-par", "--parallel",\
                    default = 2, help="number of parallel processes")
parser.add_argument("-blastx", "--blastx", dest = "blastx", \
                    default = 'src/ncbi-blast-2.14.0+/bin/blastx', \
                    help="blastx path")
args = parser.parse_args()


window_bp = 350
ncbi_range = 60000
outfmt = '5' # 5 BLAST XML, 11 = BLAST archive (ASN.1)

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

#event peptides from RefSeq --> SpliceDB search
combined_Enosi_dir = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search',\
                                                'combined_Enosi_Output')
event_fn = os.path.join(combined_Enosi_dir, 'event_parsed.txt')
eventgroup_summary = os.path.join(combined_Enosi_dir, 'event_group_summary.txt')

#reference files
reference_dir = os.path.join(base_dir , 'data', 'reference')
#Ensembl
ensembl_dir = os.path.join(reference_dir, 'Ensembl', 'release-99')
ensembl_prot_fasta = os.path.join(ensembl_dir, 'Salmo_salar.ICSASG_v2.pep.all.fa')
ensembl_proteins = readFasta(ensembl_prot_fasta)
#RefSeq
refseq_dir = os.path.join(reference_dir, 'RefSeq')
chromosome_genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')

#blastx output
blastx_dir = os.path.join(data_dir,'1_RefSeq_SpliceDB_Search','blastx')
 
blastx_cmds_dir = os.path.join(blastx_dir, 'cmds')
blastx_input_dir = os.path.join(blastx_dir, 'input')
blastx_output_dir = os.path.join(blastx_dir, 'output')
for directory in [blastx_cmds_dir, blastx_input_dir,\
                  blastx_output_dir]:
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

   
all_cmds = []
number_of_commands = 0
current_chromosome = ''
dna_seq = ''
with open(eventgroup_summary,'r') as EGout:
    for i,line in enumerate(EGout):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            EG = sp[colidx['#EventGroup']]
            EG_eventnums = sp[colidx['EventNum']]
            chromosome = sp[colidx['Chromosome']]
            EG_group_strand = sp[colidx['Strand']].split(',')
            EG_group_min = int(sp[colidx['Location(min)']])
            EG_group_max = int(sp[colidx['Location(max)']])
            EG_group_peptides = sp[colidx['Peptides']]
            
            if not dna_seq:
                current_chromosome = str(chromosome)
                dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, current_chromosome)
            else:
                if chromosome!=current_chromosome:
                    current_chromosome = str(chromosome)
                    dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, current_chromosome)                    
                
            #blastx input
            Query_seq = dna_seq[EG_group_min-window_bp:EG_group_max+window_bp]
            blastx_infile = os.path.join(blastx_input_dir, EG+'.blastx.txt')
            blastx_outfile = os.path.join(blastx_output_dir, EG+'.blastx.txt')
            logfilename = os.path.join(blastx_output_dir, EG+'.blastx.log') 
            with open(blastx_infile,'w') as blastx_infile_out:
                blastx_infile_out.write('>'+EG+'\t'+\
                             chromosome+'['+str(EG_group_min-window_bp)+'-'+str(EG_group_max+window_bp)+']'+\
                             '\t'+EG_group_peptides+'\n'+Query_seq+'\n')
            
            all_cmds.append(EG)
            number_of_commands+=1

all_cmds = sorted(all_cmds)

#write parallel cmd files
parallel_n = int(args.parallel)
print('number of commands', number_of_commands)

num_per_par = int(math.ceil(number_of_commands/float(parallel_n)))
task_num = 0
for i in range(0, number_of_commands, num_per_par):
    task_idx = range(number_of_commands)[i:i+num_per_par]
    print('\tTask',task_num,':', len(task_idx),' processes')
    blastx_sh = os.path.join(blastx_cmds_dir, 'blastx_'+str(task_num)+'.sh')
    task_num+=1
    with open(blastx_sh,'w') as outfile:
        outfile.write('\n'.join([all_cmds[z] for z in task_idx]))


blastx_path = os.path.join(base_dir, os.path.normpath(args.blastx))
if not os.path.isfile(blastx_path):
    err_msg = 'blastx path DNE'
    raise ValueError(err_msg)

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

#blastx output
blastx_dir = os.path.join(data_dir,'1_RefSeq_SpliceDB_Search','blastx')
 
blastx_cmds_dir = os.path.join(blastx_dir, 'cmds')
blastx_input_dir = os.path.join(blastx_dir, 'input')
blastx_output_dir = os.path.join(blastx_dir, 'output')
for directory in [blastx_cmds_dir, blastx_input_dir,\
                  blastx_output_dir]:
    if not os.path.isdir(directory):
        os.makedirs(directory)

outfmt = '5' # 5 BLAST XML, 11 = BLAST archive (ASN.1)
 
for cmds_fn in os.listdir(blastx_cmds_dir):
    with open(os.path.join(blastx_cmds_dir, cmds_fn),'r') as infile:
        for line in infile:
            EG = line.rstrip('\n')
            blastx_infile = os.path.join(blastx_input_dir, EG+'.blastx.txt')
            blastx_outfile = os.path.join(blastx_output_dir, EG+'.blastx.txt')
            logfilename = os.path.join(blastx_output_dir, EG+'.blastx.log')  
            cmd = [blastx_path,
                   "-query",
                   blastx_infile,
                   "-db",
                   "Actinopterygii.fsa",
                   "-evalue",
                   str(1e-10),
                   "-word_size",
                   str(5),
                   "-outfmt",
                   outfmt, # 5 BLAST XML
                   "-out",
                   blastx_outfile]
            result = cmdrun(cmd, logfilename)
            print(EG, result)