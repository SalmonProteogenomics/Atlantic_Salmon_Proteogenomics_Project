# -*- coding: utf-8 -*-
import os
import time
import urllib.request
import requests
import argparse
from subprocess import Popen, PIPE, CalledProcessError
from parse_reference_files import readFasta

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-ouput", "--ouput", dest = "ouput",\
                    default = 'data/reference/Actinopterygii/Actinopterygii_refseq',\
                    help="")
parser.add_argument("-t", "--taxon_list", dest = "taxon_list",\
                    default = 'data/reference/Actinopterygii/ancestor7898_proteomes.tsv',\
                    help="list of taxon ids with ancestor:7898 from UniProt")
parser.add_argument("-email", "--email", dest = "email",\
                    default = '',\
                    help="email to use for eutils")
parser.add_argument("-makeblastdb", "--makeblastdb", dest = "makeblastdb", \
                    default = 'src/ncbi-blast-2.14.0+/bin/makeblastdb', \
                    help="makeblastdb path")
args = parser.parse_args()
 

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))


blast_path = os.path.join(base_dir, os.path.normpath(args.makeblastdb))
if not os.path.isfile(blast_path):
    err_msg = 'blast path DNE'
    raise ValueError(err_msg)

eutils_email = str(args.email)
print(eutils_email)
if not eutils_email:
    raise ValueError('no eutils email provided')

uniprot_taxlist_fn = os.path.join(base_dir, os.path.normpath(args.taxon_list))
if not os.path.isfile(uniprot_taxlist_fn):
    raise ValueError('DNE: ', uniprot_taxlist_fn)

ouput_path = os.path.join(base_dir, os.path.normpath(args.ouput))
fasta_directory = os.path.join(ouput_path, 'fasta')
blast_db_dir = os.path.join(ouput_path, 'blast_db')
for directory in [fasta_directory, blast_db_dir]:
    if not os.path.isdir(directory):
        os.makedirs(directory)
 

eutils_base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"                
esearch = "esearch.fcgi?"                                                    
efetch = "efetch.fcgi?"                                                      

with open(uniprot_taxlist_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            taxon = sp[colidx['Taxon Id']]
            output_fasta = os.path.join(fasta_directory, taxon+'_refseq.fasta')
            if not os.path.isfile(output_fasta):
                print(taxon)

                query = 'db=protein&term=txid'+taxon+'[Organism]'
                
                url = eutils_base+ esearch + query +"&usehistory=y"
                data = requests.get(url)
                #print(url)
                content = data.content.decode('utf8')
                
                query_key = content.split('<QueryKey>')[1].split('</QueryKey>')[0]
                webenv = content.split('<WebEnv>')[1].split('</WebEnv>')[0]
                
                url = eutils_base +\
                      efetch +\
                      "db=protein" \
                      + "&query_key=" + query_key + "&WebEnv=" + webenv + \
                      "&rettype=fasta" + "&retmode=text"\
                      +"&email="+eutils_email
                
                #print(url)                
                urllib.request.urlretrieve(url, output_fasta)
                
print('prepare files for blast db creation...')
combined_fn = os.path.join(blast_db_dir, 'Actinopterygii.fsa')
map_fn = os.path.join(blast_db_dir, 'Actinopterygii_map.txt')
with open(combined_fn,'w') as outfile, open(map_fn,'w') as outfile2:
    for fn in os.listdir(fasta_directory):
        taxon = str(fn).replace('_refseq.fasta','')
        fasta_fn = readFasta(os.path.join(fasta_directory, fn))
        for accession in fasta_fn:
            protseq = fasta_fn[accession]['seq']
            sequence_split = [protseq[i:i+60] for i in range(0, len(protseq), 60)]
            outfile.write('>'+accession+'\n'+'\n'.join(sequence_split)+'\n')
            outfile2.write(accession+' '+taxon+'\n')

#create blast db
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

def cmdrun(cmd_listform):
    start_time = time.time()
    process = Popen(cmd_listform, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    # check the output to determine if a silent error occurred
    if process.returncode == 0:
        if (has_bad_exception(stdout) or has_bad_exception(stderr)):
            raise CalledProcessError(process.returncode, cmd_listform, output=stdout)
        else:
            msg = "complete: "+str(time.time()-start_time)+"s"
            return msg
    else:
        raise CalledProcessError(process.returncode, cmd_listform, output=stdout)


cmd = [blast_path,
       '-in',
       combined_fn,
       '-parse_seqids',
       '-blastdb_version',
       '5',
       '-taxid_map',
       map_fn,
       '-title',
       'Actinopterygii_refseq',
       '-dbtype',
       'prot']

result = cmdrun(cmd)
print(result)
