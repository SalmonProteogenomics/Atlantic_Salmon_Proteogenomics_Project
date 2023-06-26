# -*- coding: utf-8 -*-
import os
from pathlib import Path
import time
import argparse
from parse_reference_files import salmosalar_genomic_fna,\
      ncbi_chrname_mapping
from augustus_hints_blastx import generate_blastx_hints
from augustus_hints_peptides import generate_event_peptide_hints, generate_known_peptide_hints
from augustus_hints_ref_features import generate_ref_features_hints
from subprocess import Popen, PIPE, CalledProcessError
        
base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", dest = "dataset",\
                    default = '2019_tissues_n23',\
                    help="")
parser.add_argument("-augustus", "--augustus", dest = "augustus",\
                    default = 'src/augustus-3.3.3/bin/augustus',\
                    help="augustus path")
parser.add_argument("-run", "--run", dest = "run",\
                    default = 1, \
                    help="1: run augustus. 0: do not run augustus.")
args = parser.parse_args()

print(args.dataset)

script_start_time = time.time()

run_augustus = int(args.run)

augustus_path = os.path.join(base_dir, os.path.normpath(args.augustus))
if not os.path.isfile(augustus_path):
    print(augustus_path)
    err_msg = 'augustus path DNE'
    raise ValueError(err_msg)

augusuts_parent_dir = Path(augustus_path).parent.parent

extrinsicCfgFile = os.path.join(augusuts_parent_dir, 'config',\
                        'extrinsic', 'extrinsic.M.RM.E.W.P.cfg')

    
ncbi_range = 60000
ncbi_chrname_map = ncbi_chrname_mapping()
chrname_ncbi_map = {ncbi_chrname_map[x]:x for x in ncbi_chrname_map}

##data directory
data_dir = os.path.join(base_dir, 'data', args.dataset)

#RefSeq Ensembl db search output
RefSeq_Ensembl_dir = os.path.join(data_dir,'2_RefSeq_Ensembl_Search')
workflow_output_dir = os.path.join(RefSeq_Ensembl_dir, 'workflow_output')

#reference directory
reference_dir = os.path.join(base_dir , 'data', 'reference')

#event peptides from RefSeq --> SpliceDB search
combined_Enosi_dir = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search',\
                                  'combined_Enosi_Output')

#blastx output
blastx_dir = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search','blastx')

#output
genomic_coor_output_dir = os.path.join(RefSeq_Ensembl_dir, 'genomic_coordinates')
augustus_dir = os.path.join(data_dir, '3_Augustus')
hints_dir = os.path.join(augustus_dir, 'hints_dir')
augustus_cmds_dir = os.path.join(augustus_dir, 'cmds_dir')
augustus_outdir = os.path.join(augustus_dir, 'out_dir')   
augustus_genome_fna = os.path.join(augustus_dir,'genome_fna')
for directory in [hints_dir, genomic_coor_output_dir,\
                  augustus_cmds_dir, augustus_outdir,\
                  augustus_genome_fna]:
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
            print(stderr)
            raise CalledProcessError(process.returncode, cmd_listform, output=stdout)
        else:
            msg = "complete: "+str(time.time()-start_time)+"s"
            outfile = open(outfn, 'w')
            lines = stdout.splitlines()
            for line in lines:
                outfile.write(line.decode('utf8').strip()+'\n')
            outfile.close()
            return msg
    else:
        raise CalledProcessError(process.returncode, cmd_listform, output=stdout)


###############  Generate hints for Augustus ###############
print('parse blastx output...')
#parse blastx output and generate hints
blastx_results = generate_blastx_hints(blastx_dir, reference_dir)
blastx_hints = blastx_results['blastx_hints']
EG_blastx_query_range = blastx_results['EG_blastx_query_range']    
    
print('parse proteogenomic events...')
#parse proteogenomic events and generate hints
event_peptide_results = generate_event_peptide_hints(combined_Enosi_dir, EG_blastx_query_range)
enosi_hints = event_peptide_results['enosi_hints']
EG_blastx_query_range = event_peptide_results['EG_blastx_query_range']

#check if augustus commands already written to file
write_hints = True
all_event_groups = set()
for chromosome in enosi_hints:
    all_event_groups.update(enosi_hints[chromosome].keys())
print('all_event_groups: ', len(all_event_groups))

augustus_event_group_cmds = set([x.replace('.augustus.sh','') for x in os.listdir(augustus_cmds_dir)])
if all_event_groups==augustus_event_group_cmds:
    write_hints = False

if write_hints==True:    
    print('parse RefSeq+Ensembl db search output...')
    #parse RefSeq+Ensembl DB search output and generate hints
    refseq_ensembl_results = generate_known_peptide_hints(workflow_output_dir,\
                                                        genomic_coor_output_dir,\
                                                        reference_dir)
    ensembl_kp_hints = refseq_ensembl_results['ensembl_kp_hints']
    knownpep_hints = refseq_ensembl_results['knownpep_hints']

    print('parse reference features...')
    #parse reference features and geneate hints
    reference_features_results = generate_ref_features_hints(reference_dir)
    refseq_gff_entries = reference_features_results['refseq_gff_entries'] 
    ensembl_gff_entries = reference_features_results['ensembl_gff_entries']
   
    logfn = open(os.path.join(augustus_dir, 'Create_Augustus_input.log'),'w')  

def GenerateAugustusInput(chromosome):
    #RefSeq
    refseq_dir = os.path.join(reference_dir, 'RefSeq')
    chromosome_genomic_fna_dir = os.path.join(refseq_dir, 'GCF_000233375.1_ICSASG_v2_genomic.fna')
    chr_dna_seq = salmosalar_genomic_fna(chromosome_genomic_fna_dir, chromosome)    
    
    #use 1 index
    for EG in enosi_hints[chromosome]:
        augustus_sh_fn = os.path.join(augustus_cmds_dir, EG+'.augustus.sh')
        hints_fn = os.path.join(hints_dir, EG+'.hints.gff')
        augustus_output = os.path.join(augustus_outdir, EG+'.augustus.hints.gff')
        genome_fa = os.path.join(augustus_genome_fna, EG+'_Salmo_salar.fa')
        
        if write_hints==True:
            logfn.write('##########################\n'+chromosome+'\n')
            
            current_range = EG_blastx_query_range[chromosome][EG]
            from_loc = int(current_range[0])
            to_loc = int(current_range[1])+ncbi_range
            if from_loc>ncbi_range:
                from_loc = from_loc-ncbi_range
            else:
                from_loc = 0
            
            sub_seq = chr_dna_seq[from_loc-1:to_loc]
            #generate genomic fna subsequences
            #write as 1-index
            with open(genome_fa,'w') as outfile:
                outfile.write('>'+chromosome+' ['+str(from_loc)+'-'+str(to_loc)+']'+'\n') #1 index
                split_60 = [sub_seq[j:j+60] for j in range(0,len(sub_seq),60)]
                outfile.write('\n'.join(split_60)+'\n')
                    
            #change coordinates based on genomic fna start end
            all_hints = []
            #already included event nums within ncbi range as same group
            for entry in enosi_hints[chromosome][EG]:
                hint = list(entry[0:])
                hint_len = int(hint[4]) - int(hint[3]) +1
                hint_start = int(hint[3]) - from_loc +1 #1 index
                hint_stop = (hint_start + hint_len)-1 #1 index
                hint[3] = str(hint_start)
                hint[4] = str(hint_stop)
                all_hints.append('\t'.join(hint))
                    
            #blastx
            num_blastx_hints = 0
            if chromosome in blastx_hints:
                if EG in blastx_hints[chromosome]:
                    for entry in blastx_hints[chromosome][EG]:
                        hint = list(entry[0:])
                        hint_len = int(hint[4]) - int(hint[3]) +1
                        hint_start = int(hint[3]) - from_loc +1 #1 index
                        hint_stop = (hint_start + hint_len)-1 #1 index
                        hint[3] = str(hint_start)
                        hint[4] = str(hint_stop)
                        all_hints.append('\t'.join(hint))
                        num_blastx_hints+=1
                if num_blastx_hints!=0:
                    #print(num_blastx_hints, ' blastx hints')
                    logfn.write(EG+' '+str(num_blastx_hints)+'  blastx hints'+'\n')
    
            #known peptides
            num_knownpep_hints = 0
            for entry in knownpep_hints[chromosome]:
                hint = list(entry[0:])
                if from_loc<int(hint[3])<to_loc and from_loc<int(hint[4])<to_loc:
                    hint_len = int(hint[4]) - int(hint[3]) +1
                    hint_start = int(hint[3]) - from_loc +1 #1 index
                    hint_stop = (hint_start + hint_len)-1 #1 index
                    hint[3] = str(hint_start)
                    hint[4] = str(hint_stop)
                    all_hints.append('\t'.join(hint))
                    num_knownpep_hints+=1
            if num_knownpep_hints!=0:
                #print(num_knownpep_hints, ' knownpep hints')
                logfn.write(EG+' '+str(num_knownpep_hints)+' knownpep hints'+'\n')
    
            #refseq entries
            num_ncbi_hints = 0
            for entry in refseq_gff_entries[chromosome]:
                hint = list(entry[0:])
                if from_loc<int(hint[3])<to_loc and from_loc<int(hint[4])<to_loc:
                    hint_len = int(hint[4]) - int(hint[3]) +1
                    hint_start = int(hint[3]) - from_loc +1 #1 index
                    hint_stop = (hint_start + hint_len)-1 #1 index
                    hint[3] = str(hint_start)
                    hint[4] = str(hint_stop)
                    all_hints.append('\t'.join(hint))
                    num_ncbi_hints+=1        
            if num_ncbi_hints!=0:
                #print(num_ncbi_hints, ' refseq hints')
                logfn.write(EG+' '+str(num_ncbi_hints)+' refseq hints'+'\n')
    
            #ensembl entries
            num_ens_hints = 0
            for entry in ensembl_gff_entries[chromosome]:
                hint = list(entry[0:])
                if from_loc<int(hint[3])<to_loc and from_loc<int(hint[4])<to_loc:
                    hint_len = int(hint[4]) - int(hint[3]) +1
                    hint_start = int(hint[3]) - from_loc +1 #1 index
                    hint_stop = (hint_start + hint_len)-1 #1 index
                    hint[3] = str(hint_start)
                    hint[4] = str(hint_stop)
                    all_hints.append('\t'.join(hint))
                    num_ens_hints+=1        
            if num_ens_hints!=0:
                #print(num_ens_hints, ' ensembl hints')
                logfn.write(EG+' '+str(num_ens_hints)+' ensembl hints'+'\n')
    
            #write hints to file
            with open(hints_fn,'w') as hintfile:
                hintfile.write('\n'.join(all_hints)+'\n')

        utr_param = 'on'
        if UTR==False:
            utr_param = 'off'
        augustus_cmd = [augustus_path,
                '--progress=true',
                '--codingseq=on',
                '--AUGUSTUS_CONFIG_PATH='+os.path.join(augusuts_parent_dir,'config'),
                '--species='+species,
                '--hintsfile='+hints_fn,
                '--UTR='+utr_param,
                '--alternatives-from-evidence=true',
                '--allow_hinted_splicesites=atac',
                '--extrinsicCfgFile='+extrinsicCfgFile,
                '--softmasking=on',
                genome_fa]

        with open(augustus_sh_fn,'w') as outfile:
            outfile.write(' '.join(augustus_cmd))     
            
        if run_augustus==1:
            print(EG)
            result = cmdrun(augustus_cmd, augustus_output)
            print(result)

    return None

print('generate Augustus cmds...')
for chromosome in chrname_ncbi_map:
    print(chromosome)
    if chromosome not in enosi_hints:
        continue
    #closest species for trained model, no UTR
    species = 'zebrafish'
    UTR = False
    GenerateAugustusInput(chromosome)

if write_hints==True:
    logfn.close()

print('\n\nComplete ('+str(time.time()-script_start_time)+' seconds)')