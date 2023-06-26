'''
//Title:          buildSpliceGraph.py
//Authors:        Sunghee Woo, Seong Won Cha, Vineet Bafna
//Created:        2013
// Copyright 2007,2008,2009 The Regents of the University of California
// All Rights Reserved
//
// Permission to use, copy, modify and distribute any part of this
// program for educational, research and non-profit purposes, by non-profit
// institutions only, without fee, and without a written agreement is hereby
// granted, provided that the above copyright notice, this paragraph and
// the following three paragraphs appear in all copies.
//
// Those desiring to incorporate this work into commercial
// products or use for commercial purposes should contact the Technology
// Transfer & Intellectual Property Services, University of California,
// San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
// Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
//
// IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
// FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
// INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.
//
// THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
// OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
// ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
// REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
// EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
// THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

revised by msl043 6.12.2019

'''
import warnings
import os,argparse
import re
from subprocess import Popen, PIPE
from convert_ms2db_to_fasta import Convert_ms2db
import shutil
import time

def arguments():
    base_dir = os.path.abspath(os.path.join(__file__ ,"../../.."))
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-mem','--RAMmemory', type = str, help='Xmx memory') #3G
    parser.add_argument('-d', '--dataset', default = 'heart')
    parser.add_argument('-CoorMinReads','--min_reads_per_coordinate',\
                        default = 3)
    parser.add_argument('-DNAdir', '--dna_fasta_dir',\
                        default = os.path.join(base_dir, 'data', 'reference',\
                        'RefSeq','GCF_000233375.1_ICSASG_v2_genomic.fna'), type = str)
    
    return parser.parse_args()


def has_bad_exception(content):
    lines = content.splitlines()
    for i, line in enumerate(lines):
        lower_case_line = line.lower()
        if "[error]" in lower_case_line:
            return True
        elif "exception" in lower_case_line:
            if "javax.xml.stream.XMLStreamException: ParseError" in line and len(lines) > (i+1) and "xsi:nil" in lines[i+1]:
                continue
            elif 'com.ctc.wstx.exc.WstxParsingException: Undeclared namespace prefix "xsi" (for attribute "nil")' in line:
                return False
            else:
                return True
    return False

def cmdrun(cmd_listform, outfn):
    process = Popen(cmd_listform, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    # check the output to determine if a silent error occurred
    if process.returncode == 0:
        if (has_bad_exception(stdout) or has_bad_exception(stderr)):
            print ' '.join(cmd_listform)
            raise ValueError("Found exception, exiting with status [1].")
        else:
            outfile = open(outfn, 'w')
            lines = stdout.splitlines()
            for line in lines:
                outfile.write(line.strip()+'\n')
            outfile.close()
    else:
        print ' '.join(cmd_listform)
        raise ValueError(process.returncode)


def PrepDB(Format,SourceFileName):
    import struct, string

    class FASTACompressor:
        """
        Convert a protein database into concatenated format.
        Processes FASTA format.
        """
        def __init__(self, SourceFileName, SquishedFileName, IndexFileName):
            self.SourceFile = open(SourceFileName,"rb")
            self.SquishedFile = open(SquishedFileName,"wb")
            self.IndexFile = open(IndexFileName,"wb")
            self.SquishedFileName = SquishedFileName
            self.IndexFileName = IndexFileName
        def Compress(self):
            RecordNumber = 0
            LineNumber = 0
            FirstRecord = 1
            LineNumberWarnings = 0
            DummyTable = string.maketrans("", "")
            while (1):
                LineNumber += 1
                SourceFilePos = self.SourceFile.tell()
                FileLine  = self.SourceFile.readline()
                if not FileLine:
                    break # end o' file!
                FileLine = FileLine.strip()
                if not FileLine:
                    continue # empty lines (whitespace only) are skipped
                if FileLine[0] == ">":
                    RecordNumber += 1
                    if not FirstRecord:
                        self.SquishedFile.write("*")                
                    ID = FileLine[1:81].strip()
                    # Fix weird characters in the ID:
                    ID = ID.replace("\t", " ")
                    # Note: Important to call tell() *after* writing the asterisk!  (Fixed a bug 1/20/5)
                    SquishedFilePos = self.SquishedFile.tell() 
                    Str = struct.pack("<qi80s", SourceFilePos, SquishedFilePos, ID)
                    self.IndexFile.write(Str)
                    FirstRecord = 0
                else:
                    WarnFlag = 0
                    FileLine = string.translate(FileLine, DummyTable, " \r\n\t*")
                    FileLine = FileLine.upper()
                    Str = ""
                    for Char in FileLine:
                        if Char not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                            WarnFlag = 1
                        else:
                            Str += Char
                    #FileLine = FileLine.replace("*","")
                    if WarnFlag and LineNumberWarnings < 10:
                        print "* Warning: line %s contains non-amino-acid characters:"%LineNumber
                        print FileLine
                        LineNumberWarnings += 1
                        if LineNumberWarnings >= 10:
                            print "(omitting further warnings)"
                    self.SquishedFile.write(Str)
            print "Converted %s protein sequences (%s lines) to .trie format."%(RecordNumber, LineNumber)
            print "Created database file '%s'"%self.SquishedFileName
            
    if Format == "fasta":
        CompressorClass = FASTACompressor
    else:
        errmsg= "Unknown source database format '%s'"%Format
        raise ValueError(errmsg)
    
    SquishedFileName = "%s.trie"%os.path.splitext(SourceFileName)[0]
    IndexFileName = "%s.index"%os.path.splitext(SourceFileName)[0]

    # Use FASTACompressor for FASTA format, Compressor for the weird swiss-prot format
    Squasher = CompressorClass(SourceFileName, SquishedFileName, IndexFileName)
    Squasher.Compress()
    
    return None



def BuildSpliceGraph(args):

    ###===Functions===###
    def baseconvert(n, base):
        """convert positive decimal integer n to equivalent in another base (2-36)"""
        digits = "0123456789abcdefghijklmnopqrstuvwxyz"
        try:
            n = int(n)
            base = int(base)
        except:
            return ""
    
        if n < 0 or base < 2 or base > 36:
            return ""
    
        s = ""
        while 1:
            r = n % base
            s = digits[r] + s
            n = n / base
            if n == 0:
                break
    
        return s
    
    def checkParentID(part1, part2):
        if part1[0] == -1 or part2[0] == -1:
            return False
        tmp1 = part1[-1].split(";")[0].split("=")[1]
        tmp2 = part2[-1].split(";")[0].split("=")[1]
        if tmp1 == tmp2:
            return True
        else:
            return False
    
    def decomp(CIGAR,CIGAR_begin,sign,strand):
        del_location = []
        seq_location = []
        pattern = re.compile('S|N|M|D|I|[0-9]+')
        CIGAR_info = [(m.group()) for m in pattern.finditer(CIGAR)]
        CIGAR_interval = [int(n) for i,n in enumerate(CIGAR_info) if i%2 == 0]
        index_D = [i for i,n in enumerate(CIGAR_info) if n == sign]
        for i in index_D:
            del_interval = int(CIGAR_info[i-1])
            del_start = CIGAR_begin
            for j in range((i-1)/2):
                if CIGAR_info[2*j+1] != 'I' and  CIGAR_info[2*j+1] != 'S':
                    del_start += CIGAR_interval[j]
            #del_start = CIGAR_begin + sum(CIGAR_interval[0:(i-1)/2])
            del_end = del_start + del_interval
            del_location.append((del_start,del_end))
            seq_finder = 0
            for j in range((i-1)/2):
                if CIGAR_info[2*j+1] != 'D' and CIGAR_info[2*j+1] != 'N':
                    seq_finder += CIGAR_interval[j]
            seq_location.append((seq_finder,seq_finder+CIGAR_interval[(i-1)/2]))
        return [del_location,seq_location]
    
    def addSAM(SourceFile):
        check_sample = [0,1000] #0: check off, 1: check on
        id = 1
        for line in SourceFile:
            if line == "":
                continue
            if line.startswith(">"):
                continue
            if line.startswith("@"):
                continue
            if line.startswith("track"):
                continue
            data = line.strip()
            data = data.split()
            QNAME = data[0]
            FLAG = int(data[1])
            RNAME = data[2]
            POS = int(data[3])-1
            CIGAR = data[5]
            TLEN = int(data[8])
    	    
            if check_sample[0] == 1:
                if 'N' in CIGAR:
                    id += 1
                    if id > check_sample[1]:
                        break
    	    
            tmp_flag = baseconvert(FLAG,2)
            if tmp_flag == "0" or tmp_flag == 0:
                strand = "+"
            elif len(tmp_flag)<5:
                strand = "+"
            elif int(tmp_flag[-5]) == 0:
                strand = "+"
            elif int(tmp_flag[-5]) == 1:
                strand = "-"
            else:
                print QNAME, CIGAR, POS, TLEN, "can't parse strand from: ", FLAG
                continue
    	    
            chr_name = RNAME 
        		    
            if chr_name not in chr_ref:
                print chr_name+' not in chr_ref'
                continue
    	    
            if chr_name == -1:
                continue
    	    #############################################
            CIGAR_begin = POS
    	    #############################################
            if 'S' in CIGAR:
                continue
    	
            if 'N' in CIGAR:   
                sign = 'N'
                [data_location, seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
                data = full_data.get(chr_name)
                for coord_info in data_location:  
                    if coord_info[0] in data:
                        check = 0
                        for i,check_end in enumerate(data[coord_info[0]]):
                            if check_end[0] == coord_info[1] and strand == check_end[2]:
                                if SAM_file_name in check_end[1]:
                                    data[coord_info[0]][i][1][SAM_file_name].append(QNAME)
                                else:
                                    data[coord_info[0]][i][1][SAM_file_name] = [QNAME]  	    					
                                check = 1
                        if check == 0:
                            data[coord_info[0]].append([coord_info[1],{SAM_file_name:[QNAME]},strand])
                    else:
                        data[coord_info[0]] = [[coord_info[1],{SAM_file_name:[QNAME]},strand]]
                full_data[chr_name] = data
    
    ###=== ===###    
    base_dir = os.path.abspath(os.path.join(__file__ ,"../../.."))
    RAMgb = args.RAMmemory
 
    dataset_dir = os.path.join(base_dir, 'data', str(args.dataset))
    SAM_directory = os.path.join(dataset_dir, 'SAM_files')
    SPLICEdb_directory = os.path.join(dataset_dir, 'SpliceDB')
    
    path_to_dna_fasta_files = os.path.normpath(args.dna_fasta_dir)
    min_read_counts = int(args.min_reads_per_coordinate)

    RAMgb = args.RAMmemory
    construct_splice_graph_from_tmp = os.path.join(base_dir, 'src',\
                        'SpliceDBcreation', 'ConstructSpliceGraphFromTmp.jar') #ConstructSpliceGraphFromTmp.jar

    
    print '\nRead genomic fasta files...'
    chr_ref = []
    dna_TRIE_file_name_per_chr = {}
    Files = os.listdir(path_to_dna_fasta_files)    
    
    for F in Files:
        file_path_and_name = os.path.join(path_to_dna_fasta_files,F)
        (Stub, Extension) = os.path.splitext(file_path_and_name)
    
        chr_ref_tmp = F.split("_")[0]
        chr_ref.append(chr_ref_tmp)
			
        if chr_ref_tmp not in dna_TRIE_file_name_per_chr:
            dna_TRIE_file_name_per_chr[chr_ref_tmp] = Stub + ".trie"
        
        if os.path.exists(Stub + ".trie") == False:
            PrepDB("fasta", file_path_and_name)
            
    chr_ref = list(set(chr_ref))
  
    print '\nConstructing SpliceDB...'
    constructed_splice_databases = []
    for sub_dir in os.listdir(SAM_directory):
        sub_dir_start_time = time.time()
        print(sub_dir)
        sub_dir_path = os.path.join(SAM_directory, sub_dir)
        SAM_file_list = sorted(os.listdir(sub_dir_path))
        if not SAM_file_list:
            errmsg = 'no RNAseq files detected in '+sub_dir+'...'
            print errmsg
            warnings.warn(errmsg)
            continue
    
        # Using chr representation, we initiate the data structure 
        full_data = {}
        in_data = {}
        del_data = {}
        mu_data = {}
        for i in range(0,len(chr_ref)):
            full_data[chr_ref[i]] = {}
            in_data[chr_ref[i]] = {}
            del_data[chr_ref[i]] = {}
            mu_data[chr_ref[i]] = {}
            
        #output directories
        intermediate_out_dir = os.path.join(SPLICEdb_directory,\
                                 sub_dir, 'intermediate_files')
        outputdir = os.path.join(SPLICEdb_directory,\
                                 sub_dir, 'final_output')
        for directory in [intermediate_out_dir, outputdir]:
            if not os.path.isdir(directory):
                print 'creating directory: ', directory
                os.makedirs(directory)

        ##Retreive folder
        print '\nParsing SAM files...'
        file_name_info = {}
        for F in SAM_file_list:
        	
            if F in file_name_info:
                print 'Current file already in memory :', F
                continue
            
            SAM_file_name = len(file_name_info)
            
            file_name_info[F] = SAM_file_name
            file_path_and_name = os.path.join(sub_dir_path,F)
            
            print "Processing .sam file :",F
    
            SourceFile = open(file_path_and_name,'r')
            addSAM(SourceFile)
    
        sorted_chr_list = sorted(full_data.keys())
        
        print '\nConstructing SpliceDB per chromosome...'
        for chromosome in sorted_chr_list:
            print chromosome
            Stub = "Merged_BAM_SAM_SPL"
            Output_TMP_file_name_per_chr = os.path.join(intermediate_out_dir, Stub+"_"+chromosome+".spl")
            spl_with_readQNAME = os.path.join(intermediate_out_dir, Stub+"_"+chromosome+".QNAME.spl")
            Output_ms2db_file_name_per_chr = os.path.join(intermediate_out_dir, Stub+"_"+chromosome+".ms2db")
            Output_fasta_file_name_per_chr = os.path.join(intermediate_out_dir, Stub+"_"+chromosome+".fa")
            ms2db_logfn = os.path.join(intermediate_out_dir, Stub+"_"+chromosome+".spl.log")
        
            print "Writing output .spl file :", os.path.basename(Output_TMP_file_name_per_chr)
            s = open(Output_TMP_file_name_per_chr,'w')
            s2 = open(spl_with_readQNAME,'w')
        
            s.write('#')
            s2.write('#')
            for i in file_name_info:
                s.write(i+':'+str(file_name_info[i]))
                s2.write(i+':'+str(file_name_info[i]))
                if i != file_name_info.keys()[-1]:
                    s.write(',')
                    s2.write(',')
            s.write('\n')
            s2.write('\n')
        
            s.write('#Splice'+'\n')
            s2.write('#Splice'+'\n')
            splice_data = sorted(full_data[chromosome].items(), key = lambda x: float(x[0]))
            for key,value in splice_data:
                #data[coord_info[0]].append([coord_info[1],{SAM_file_name:[QNAME]},strand])
                #key == coord_info[0] (start coordinate?)
                #value is a list of end coordinate, dictionary of sam filename and read information, strand
                for entry in value:
                    #end_coor = entry[0]
                    #Reads_list = entry[1]
                    #strand = entry[2]
                    readcounts = sum([len(QNAMES) for QNAMES in entry[1].values()])
                    if readcounts < min_read_counts:
                        continue
                    s.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    s2.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    for i in entry[1]:
                        #SAM_file_name == i
                        s.write(str(i)+':'+str(len(entry[1][i])))
                        s2.write(str(i)+':'+'@'.join(entry[1][i]))
                        if i != entry[1].keys()[len(entry[1])-1]:
                            s.write(',')
                            s2.write(',')
                    s.write('\t'+entry[2]+'\n')
                    s2.write('\t'+entry[2]+'\n')
        	
            s.write('#Deletion'+'\n')
            s2.write('#Deletion'+'\n')
            deletion = sorted(del_data[chromosome].items(), key = lambda x: float(x[0]))
            for key,value in deletion:
                for entry in value:
                    readcounts = sum([len(QNAMES) for QNAMES in entry[1].values()])
                    if readcounts < min_read_counts:
                        continue
                    s.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    s2.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    for i in entry[1]:
                        s.write(str(i)+':'+str(len(entry[1][i])))
                        s2.write(str(i)+':'+'@'.join(entry[1][i]))
                        if i != entry[1].keys()[len(entry[1])-1]:
                            s.write(',')
                            s2.write(',')
                    s.write('\t'+entry[2]+'\n')		
                    s2.write('\t'+entry[2]+'\n')	
              		
            s.write('#Insertion'+'\n')
            s2.write('#Insertion'+'\n')
            insertion = sorted(in_data[chromosome].items(), key = lambda x: float(x[0]))
            for key,value in insertion:
                for entry in value:
                    readcounts = sum([len(QNAMES) for QNAMES in entry[1].values()])
                    if readcounts < min_read_counts:
                        continue
                    s.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    s2.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    for i in entry[1]:
                        s.write(str(i)+':'+str(len(entry[1][i])))
                        s2.write(str(i)+':'+'@'.join(entry[1][i]))
                        if i != entry[1].keys()[len(entry[1])-1]:
                            s.write(',')
                            s2.write(',')
                    s.write('\t'+entry[2]+'\n')
                    s2.write('\t'+entry[2]+'\n')
       
            s.write('#Mutation'+'\n')
            s2.write('#Mutation'+'\n')
            Mutation = sorted(mu_data[chromosome].items(), key = lambda x: float(x[0]))
            for key,value in Mutation:
                for entry in value:
                    readcounts = sum([len(QNAMES) for QNAMES in entry[1].values()])
                    if readcounts < min_read_counts:
                        continue
                    s.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    s2.write(chromosome+'\t'+str(key)+'\t'+str(entry[0])+'\t')
                    for i in entry[1]:
                        s.write(str(i)+':'+str(len(entry[1][i])))
                        s2.write(str(i)+':'+'@'.join(entry[1][i]))
                        if i != entry[1].keys()[len(entry[1])-1]:
                            s.write(',')
                            s2.write(',')
                    s.write('\t'+entry[2]+'\n')
                    s2.write('\t'+entry[2]+'\n')
    
            s.close()
            s2.close()
        
            print "Constructing splice graph from .spl file: ", os.path.basename(Output_ms2db_file_name_per_chr)
            cmd = ['java',
                   "-Xmx"+RAMgb+"G",
                   "-jar",
                   construct_splice_graph_from_tmp,
                   "-p",
                   Output_TMP_file_name_per_chr,
                   "-w",
                   Output_ms2db_file_name_per_chr,
                   "-s",
                   dna_TRIE_file_name_per_chr[chromosome]]
    
            cmdrun(cmd, ms2db_logfn)
        	
            print "Converting splice graph(.spl) to .fa file: ", os.path.basename(Output_fasta_file_name_per_chr)
            print Convert_ms2db(Output_ms2db_file_name_per_chr, Output_fasta_file_name_per_chr, 30)
    
        print 'Merging SpliceDB files...'
        with open(os.path.join(outputdir, 'merged_SPLICEdatabase.fa'),'w') as outfile:
            for fn in os.listdir(intermediate_out_dir):
                if fn.endswith('.fa'):
                    with open(os.path.join(intermediate_out_dir, fn),'r') as infile:
                        for line in infile:
                            outfile.write(line.strip()+'\n')
        
        #remove intermediate directory 
        try:
           shutil.rmtree(intermediate_out_dir)
           print "directory removed successfully"
        except OSError as x:
           print "Error occured: %s : %s" % (intermediate_out_dir, x.strerror)
    
        print '\nsub_dir SpliceDB: '+str(time.time()-sub_dir_start_time)+' seconds...\n'

        constructed_splice_databases.append(sub_dir)

    if len(constructed_splice_databases)==len(os.listdir(SAM_directory)):
        return True
    else:
        message = 'constructed '+ str(len(constructed_splice_databases))+' databases out of '+ str(len(os.listdir(SAM_directory))) + ' input directories provided...'
        return message


#### MAIN ####
def main():
    start_time = time.time()
    args = arguments()
    results = BuildSpliceGraph(args)
    if results==True:
        print '\n\nConstruction of SpliceDB Complete ('+str(time.time()-start_time)+' seconds)'
    else:
        print results
        warnings.warn(results)

if __name__ == "__main__":
    main() 