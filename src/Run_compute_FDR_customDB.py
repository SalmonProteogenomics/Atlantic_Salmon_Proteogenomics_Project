# -*- coding: utf-8 -*-
'''
Author: Miin S. Lin
Created: June 26, 2023
'''

import argparse, os
from ComputeFDR import Compute_FDR
from collections import OrderedDict
import shutil

def arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-d", "--dataset", dest = "dataset",\
                        default = '2019_tissues_n23', help="")
    parser.add_argument('-wo', '--wo_output', dest = "wo_output",\
                        default = '2_RefSeq_Ensembl_Search')
    return parser.parse_args()

def readFasta_add_decoy(fasta_list):
    
    decoy_prefix = 'XXX_'
    
    fasta_proteins = {}
    for fastafn in fasta_list:
        sequence = []
        accession = ''
        with open(fastafn,'rb') as infile:
            for line in infile:
                if line[0]=='>':
                    if sequence:
                        fasta_proteins[accession] = ''.join(sequence)
                        #add decoy
                        fasta_proteins[decoy_prefix+accession] = ''.join(sequence)[::-1] 
                        
                    accession = line.strip()[1:].split()[0]
                    sequence = []
                else:
                    sequence.append(line.strip())
            if sequence:
                fasta_proteins[accession] = ''.join(sequence)
                #add decoy
                fasta_proteins[decoy_prefix+accession] = ''.join(sequence)[::-1]
    return fasta_proteins

def compute_fdr(args):
    base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))
    ##reference files
    reference_dir = os.path.join(base_dir, 'data', 'reference')
    
    ##data directory
    data_dir = os.path.join(base_dir, 'data', args.dataset)
    workflow_dir = os.path.join(data_dir, args.wo_output, 'workflow_output')
    
    db_fn = ''
    if args.wo_output == '2_RefSeq_Ensembl_Search':
        db_fn = os.path.join(data_dir, '1_RefSeq_SpliceDB_Search',\
                             'combined_Enosi_Output', 'GCF_000233375.1_ICSASG_v2_protein+Ensembl.fasta')
    if args.wo_output == '5_RefSeq_Ensembl_Augustus_Search':
        db_fn = os.path.join(data_dir, args.wo_output,\
                             'fasta_db', 'RefSeq+Ensembl+Augustus.fasta')
    
    if not os.path.isfile(db_fn):
        raise ValueError('database path dow not exist...')
    
    print('Reading database...')    
    Database = readFasta_add_decoy([db_fn])    
    #GPM contaminants database
    print('Reading GPM contaminants fasta file...')
    gpm_fn = os.path.join(reference_dir, 'contaminant_db','gpm_cRAP_012915.fasta')
    contaminant_proteins = readFasta_add_decoy([gpm_fn])
    
    for ac in contaminant_proteins:
        Database[ac] = str(contaminant_proteins[ac])
    
    #MSGF peptidematch CMD output directories 
    EM_dir = os.path.join(workflow_dir, 'ExactMatch_output')
    
    #create MSGF_fdrdir
    fdr_dir = os.path.join(workflow_dir, 'MSGF_fdrdir')
    if os.path.isdir(fdr_dir):
        shutil.rmtree(fdr_dir)
    os.makedirs(fdr_dir)
    
    all_combinedfn = os.path.join(workflow_dir, 'MSGF_combinedtsv','MSGF_tsvdir.tsv')
    
    fdr = float(0.01)
    PeptideLevel_score_FDR = os.path.join(fdr_dir, 'PeptideLevel_score_FDR_values.txt')
    PSMlevel_txt = os.path.join(fdr_dir, os.path.basename(all_combinedfn)[:-4]+'_PSMLevelFDR_targets.fdr')
    PSMlevel_decoys = os.path.join(fdr_dir, os.path.basename(all_combinedfn)[:-4]+'_PSMLevelFDR_decoys.fdr')
    Precursorlevel_txt = os.path.join(fdr_dir, os.path.basename(all_combinedfn)[:-4]+'_PrecursorLevelFDR_targets.fdr')
    Precursorlevel_decoys = os.path.join(fdr_dir, os.path.basename(all_combinedfn)[:-4]+'_PrecursorLevelFDR_decoys.fdr')
    
    decoy_prefix = 'XXX_'
    peptide_protein = {}
    for exact_match_fn in os.listdir(EM_dir):
        with open(os.path.join(EM_dir, exact_match_fn), 'rb') as infile:
            for line in infile:
                sp = line.strip().split('\t')
                if line[0:2]=='##':
                    col = {x.replace('#',''):xi for xi,x in enumerate(sp)}
                elif line[0]!='#':
                    peptide = sp[col['Query']].replace('I','L')
                    protein = sp[col['Subject']]
                    if protein=='No match':
                        continue
                    if peptide not in peptide_protein:
                        peptide_protein[peptide] = set()
                    peptide_protein[peptide].add(protein)
    
    Peptide_TD = {}
    proteins_peptides = {}
    for pep in peptide_protein:
        Peptide_TD[pep] = {}
        entry = peptide_protein[pep]
        decoy_matches = [P for P in entry if P[:4]==decoy_prefix]
        target_matches = entry-set(decoy_matches)
        
        if not decoy_matches:
            Peptide_TD[pep]['td'] = '_'
            Peptide_TD[pep]['Proteins'] = entry
        else:
            if not target_matches:
                Peptide_TD[pep]['td'] = decoy_prefix
                Peptide_TD[pep]['Proteins'] = entry
            else:
                Peptide_TD[pep]['td'] = '_'
                Peptide_TD[pep]['Proteins'] = target_matches

        for prot in Peptide_TD[pep]['Proteins']:
            if prot not in proteins_peptides:
                proteins_peptides[prot] = set()
            proteins_peptides[prot].add(pep)
        
    ############################################################################
    MSGF_tsvdir_A_dir = os.path.join(workflow_dir, 'MSGF_tsvdir_A')
    header = ''
    PSMs_all= []
    
    if os.path.isfile(all_combinedfn):
        header = ''
        matched_spectra = OrderedDict()
        print('Reading MSGF_tsvdir.tsv...')
        with open(all_combinedfn,'rb') as infile:
            for i,line in enumerate(infile):
                sp = line.strip().split('\t')
                if i==0:
                    col = {x:xi for xi,x in enumerate(sp)}
                    header = sp[0:]
                else:
                    specf = sp[col['#SpecFile']].replace('.mgf','')
                    scannum = sp[col['ScanNum']]
                    if specf not in matched_spectra:
                        matched_spectra[specf] = OrderedDict()
                    matched_spectra[specf][scannum] = sp[0:]
    
    if not os.path.isfile(all_combinedfn):
        if not os.path.isdir(MSGF_tsvdir_A_dir):
            raise ValueError('MSGF_tsvdir_A directory DNE...')
        print('MSGF_tsvdir.tsv DNE...parsing tsv directory...')
        header = ['#SpecFile', 'SpecID', 'ScanNum', 'Title', 'FragMethod', 'Precursor', 'IsotopeError', 'PrecursorError(ppm)', 'Charge', 'Peptide', 'Protein', 'DeNovoScore', 'MSGFScore', 'SpecEValue', 'EValue', 'peptide_nonmodseq', 'precursor_pepcharge', 'target_decoy']
        tsv_list = []
        tsv_to_parent_dir = {}
        for x in os.listdir(MSGF_tsvdir_A_dir):
            if x.endswith('.tsv') and os.path.isfile(os.path.join(MSGF_tsvdir_A_dir, x)):
                tsv_to_parent_dir[x] = MSGF_tsvdir_A_dir
                tsv_list.append(x)

        tsv_list = sorted(tsv_list)
        tsv_list = [os.path.join(tsv_to_parent_dir[x], x) for x in tsv_list]
        print('Parsing ',len(tsv_list), ' tsv files...')
 
        matched_spectra = {}
        for fn in tsv_list:
            with open(fn, 'rb') as infile:
                for line in infile:
                    sp = line.strip().split('\t')
                    if line[0]=='#':
                        sp.append('peptide_nonmodseq')
                        sp.append('precursor_pepcharge')
                        sp.append('target_decoy')
                        col = {x:xi for xi,x in enumerate(sp)}
                         
                    else:
                        SpecFile = os.path.splitext(sp[col['#SpecFile']])[0]
                        if SpecFile not in matched_spectra:
                            matched_spectra[SpecFile] = {}                
                        Charge = sp[col['Charge']]
                        Title = sp[col['Title']]
                        if 'File:""' in Title:
                            #use <RunID> as filename?...
                            RunID = '.'.join(Title.split('File')[0].strip().split('.')[:-3])
                            Title = Title.replace('File:""','File:"'+RunID+'.EXT'+'"')
                            sp[col['Title']] = Title
                        msgfscore = float(sp[col['MSGFScore']])
                        specevalue = float(sp[col['SpecEValue']])                       
                        scannum = Title.strip().split('scan=')[1].split('"')[0]
                        sp[col['SpecID']] = 'index='+scannum
                        sp[col['ScanNum']] = scannum
                        Peptide = sp[col['Peptide']][2:-2]
                        Charge = sp[col['Charge']]                
                        nonmod_peptide = ''.join([aa for aa in Peptide if aa.isalpha()])
                        nonmmod_peptide_IL = nonmod_peptide.replace('I','L')
                        td_type = Peptide_TD[nonmmod_peptide_IL]['td']
                        proteins = Peptide_TD[nonmmod_peptide_IL]['Proteins']
                        sp[col['Protein']] = ';'.join(proteins)
                        Precursor = Peptide.replace('I','L')+'|'+Charge
                        sp.append(nonmod_peptide)
                        sp.append(Precursor)
                        sp.append(td_type)
                        
                        if scannum in matched_spectra[SpecFile]:
                            #larger msgf+ score the better
                            PSM_score_diff = msgfscore-float(matched_spectra[SpecFile][scannum][col['MSGFScore']])
                            #if MSGFscore the same
                            if PSM_score_diff==0:
                                #if currently stored is decoy
                                if matched_spectra[SpecFile][scannum][col['target_decoy']]==decoy_prefix:
                                    # if current psm is target, replace decoy with target
                                    if td_type!=decoy_prefix:
                                        matched_spectra[SpecFile][scannum] = sp[0:]
                                    # else if current psm is also decoy, choose decoy with lower specevalue
                                    else:
                                        if specevalue<float(matched_spectra[SpecFile][scannum][col['SpecEValue']]):
                                            matched_spectra[SpecFile][scannum] = sp[0:]                                      
                                # if currently stored is target
                                else:
                                    # if current psm is also target
                                    if td_type!=decoy_prefix:
                                        # if current specevalue less than stored specevalue, replace
                                        if specevalue<float(matched_spectra[SpecFile][scannum][col['SpecEValue']]):
                                            matched_spectra[SpecFile][scannum] = sp[0:]
                                continue
                            #larger msgf+ score used to choose peptide assigned to spectrum
                            if PSM_score_diff/abs(PSM_score_diff) == 1:
                                matched_spectra[SpecFile][scannum] = sp[0:]
                        #initialize        
                        elif scannum not in matched_spectra[SpecFile]:
                            matched_spectra[SpecFile][scannum] = sp[0:]

    #PSM level FDR (pooled)
    PSMs_all = []
    for specf in matched_spectra:
        entries = matched_spectra[specf].values()
        for sp in entries:
            PSMs_all.append(sp)
    
    ############################################################################
    #FDR computation
    FDR_calc_score_direction = -1
    FDRcompute_targetpass_col = 0
    FDRcompute_decoypass_col = 1
    FDRcompute_pvalue_col = 2
    FDRcompute_tdnum_col = 3
    FDRcompute_FDR_col = 4
    FDRcompute_FDR_values = 5

    print('\nComputing peptide level FDR...')
    #Peptide Level FDR (pooled)
    Peptides_all = {}
    for specf in matched_spectra:
        entries = matched_spectra[specf].values()
        for sp in entries:
            nonmodpep = sp[col['peptide_nonmodseq']].replace('I','L')
            PSMscore = float(sp[col['SpecEValue']])
            td_type = sp[col['target_decoy']]
            
            # choose best psm for pepgroup
            if nonmodpep in Peptides_all:
                PSMovw = PSMscore-float(Peptides_all[nonmodpep][col['SpecEValue']])
                #if same score
                if PSMovw==0:
                    #if currently stored is decoy, if current psm is target, replace decoy with targets
                    if Peptides_all[nonmodpep][col['target_decoy']]==decoy_prefix and td_type!=decoy_prefix:
                        Peptides_all[nonmodpep] = sp[0:]
                    continue
                #if better score, replace
                if PSMovw/abs(PSMovw) == FDR_calc_score_direction:
                    Peptides_all[nonmodpep] = sp[0:]
            #initialize
            else:
                Peptides_all[nonmodpep] = sp[0:]
    
    #MS-GF:PepQValue
    #Peptide-level Q-value estimated using the target-decoy approach.
    #Reported only if "-tda 1" is specified.
    #If multiple spectra are matched to the same peptide, only the best scoring PSM (lowest SpecProb) is retained.
    #After that, MS-GF:PepQValue is calculated as #DecoyPSMs>s / #TargetPSMs>s among the retained PSMs.
    #This approximates the Q-value of the set of unique peptides.
    #In the MS-GF+ output, the same PepQValue value is given to all PSMs sharing the peptide.
    #Thus, even a low-quality PSM may get a low PepQValue (if it has a high-quality "sibling" PSM sharing the peptide).
    #Note that this should not be used to count the number of identified PSMs.

    Peptides_all_list = Peptides_all.values()
    del Peptides_all
    PeptidelevelFDR = Compute_FDR(Peptides_all_list, 0.90, col['SpecEValue'], FDR_calc_score_direction, col['target_decoy'], decoy_prefix)
    if PeptidelevelFDR == 'cannot reach specified FDR':
        raise ValueError('Cannot reach specified FDR.')
    
    #(s, pass_decoy_scores, pass_target_scores, currentrate)
    with open(PeptideLevel_score_FDR,'wb') as outfile:
        peplevel_fdr_header = header[0:]
        peplevel_fdr_header.extend(['pass_decoy_scores','pass_target_scores','PepQValue'])
        outfile.write('\t'.join(peplevel_fdr_header)+'\n')
        
        PeptideLevel_FDRvalues = {}
        for s in PeptidelevelFDR[FDRcompute_FDR_values]:
            score = str(float(s[0][0]))
            td_type = s[0][1]
            pass_decoy_scores = str(s[1])
            pass_target_scores = str(s[2])
            PeptideLevel_FDRvalues[score] = {'td_type':td_type,
                                              'pass_decoy_scores':pass_decoy_scores,\
                                              'pass_target_scores':pass_target_scores}
        #one pqvalue score for each peptide
        for sp in PeptidelevelFDR[FDRcompute_targetpass_col]:
            specevalue = str(float(sp[col['SpecEValue']]))
            FDR_value = PeptideLevel_FDRvalues[specevalue]
            pass_decoy_scores = float(FDR_value['pass_decoy_scores'])
            pass_target_scores = float(FDR_value['pass_target_scores'])
            PepQValue = pass_decoy_scores/pass_target_scores
            entry = sp[0:]
            entry.extend([str(FDR_value['pass_decoy_scores']), str(FDR_value['pass_target_scores']), str(PepQValue)])
            outfile.write('\t'.join(entry)+'\n')
    
    print('### Peptide-level ####')
    print('p-value threshold: ', PeptidelevelFDR[FDRcompute_pvalue_col])
    print('decoy|target: ', PeptidelevelFDR[FDRcompute_tdnum_col])
    print('FDR: ', PeptidelevelFDR[FDRcompute_FDR_col])


    print('\nComputing PSM level FDR...')       
    PSMlevelFDR = Compute_FDR(PSMs_all, fdr, col['SpecEValue'], FDR_calc_score_direction, col['target_decoy'], decoy_prefix)
    if PSMlevelFDR == 'cannot reach specified FDR':
        raise ValueError('Cannot reach specified 1% level FDR.')

    target_PSMs_all = PSMlevelFDR[FDRcompute_targetpass_col]
    with open(PSMlevel_txt,'wb') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for sp in target_PSMs_all:
            outfile.write('\t'.join(sp)+'\n')
    
    decoy_PSMs_all = PSMlevelFDR[FDRcompute_decoypass_col]
    with open(PSMlevel_decoys,'wb') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for sp in decoy_PSMs_all:
            outfile.write('\t'.join(sp)+'\n')
    
    print('### 1% PSM-level FDR ####')
    print('p-value threshold: ', PSMlevelFDR[FDRcompute_pvalue_col])
    print('decoy|target: ', PSMlevelFDR[FDRcompute_tdnum_col])
    print('FDR: ', PSMlevelFDR[FDRcompute_FDR_col])
    del PSMlevelFDR
    
    print('\nComputing precursor level FDR...')
    #Precursor level FDR (pooled)
    Precursors_all = {}
    for specf in matched_spectra:
        entries = matched_spectra[specf].values()
        for sp in entries:
            precursor = sp[col['precursor_pepcharge']]
            PSMscore = float(sp[col['SpecEValue']])
            td_type = sp[col['target_decoy']]
            
            # choose best psm for pepgroup
            if precursor in Precursors_all:
                PSMovw = PSMscore-float(Precursors_all[precursor][col['SpecEValue']])
                #if same score
                if PSMovw==0:
                    #if currently stored is decoy, if current psm is target, replace decoy with targets
                    if Precursors_all[precursor][col['target_decoy']]==decoy_prefix and td_type!=decoy_prefix:
                        Precursors_all[precursor] = sp[0:]
                    continue
                #if better score, replace
                if PSMovw/abs(PSMovw) == FDR_calc_score_direction:
                    Precursors_all[precursor] = sp[0:]
            #initialize
            else:
                Precursors_all[precursor] = sp[0:]
    
    Precursors_all_list = Precursors_all.values()
    del Precursors_all
    PrecursorlevelFDR = Compute_FDR(Precursors_all_list, fdr, col['SpecEValue'], FDR_calc_score_direction, col['target_decoy'], decoy_prefix)
    if PrecursorlevelFDR == 'cannot reach specified FDR':
        raise ValueError('Cannot reach specified 1% level FDR.')
    
    target_Precursors_all = PrecursorlevelFDR[FDRcompute_targetpass_col]
    with open(Precursorlevel_txt,'wb') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for sp in target_Precursors_all:
            outfile.write('\t'.join(sp)+'\n')
    
    decoy_Precursors_all = PrecursorlevelFDR[FDRcompute_decoypass_col]
    with open(Precursorlevel_decoys,'wb') as outfile:
        outfile.write('\t'.join(header)+'\n')
        for sp in decoy_Precursors_all:
            outfile.write('\t'.join(sp)+'\n')
    
    print('### 1% Precursor-level FDR ####')
    print('p-value threshold: ', PrecursorlevelFDR[FDRcompute_pvalue_col])
    print('decoy|target: ', PrecursorlevelFDR[FDRcompute_tdnum_col])
    print('FDR: ', PrecursorlevelFDR[FDRcompute_FDR_col])
    del PrecursorlevelFDR
    
    ###########################################################################
    print('\nComputing protein level FDR...')
    
    #Protein level FDR (pooled)
    protein_level_fdr = 0.01
    if protein_level_fdr=='NA':
        print('skipping protein level computation...')
        return None
    
    from operator import itemgetter
    import math, time
    
    protein_level_fdr = float(protein_level_fdr)
    Proteinlevel_PSMs_txt = os.path.join(fdr_dir, 'Proteinlevel_'+str(protein_level_fdr)+'_FDR_PSMs.fdr')
    Proteinlevel_proteins_txt = os.path.join(fdr_dir, 'Proteinlevel_'+str(protein_level_fdr)+'_FDR_proteins.txt')

    ############### Precursors ###############
    Precursors_scores = {}
    peptide_precursor = {}
    Precursor_list = {}
    with open(Precursorlevel_txt,'rb') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                colidx = {x:xi for xi,x in enumerate(sp)}
            else:
                pep = sp[colidx['peptide_nonmodseq']].replace('I','L')
                SpecEValue = float(sp[colidx['SpecEValue']])
                precursor = sp[colidx['precursor_pepcharge']]
                Precursors_scores[precursor] = -1*math.log10(SpecEValue)
                if pep not in peptide_precursor:
                    peptide_precursor[pep] = set()
                peptide_precursor[pep].add(precursor)
                Precursor_list[precursor] = sp[0:]
    with open(Precursorlevel_decoys,'rb') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            if line[0]=='#':
                colidx = {x:xi for xi,x in enumerate(sp)}
            else:
                pep = sp[colidx['peptide_nonmodseq']].replace('I','L')
                SpecEValue = float(sp[colidx['SpecEValue']])
                precursor = sp[colidx['precursor_pepcharge']]
                Precursors_scores[precursor] = -1*math.log10(SpecEValue)
                if pep not in peptide_precursor:
                    peptide_precursor[pep] = set()
                peptide_precursor[pep].add(precursor)
                Precursor_list[precursor] = sp[0:]
    
    ############### Proteins ############### 
    all_matched_proteins = list(proteins_peptides.keys())
    Protein_sequences = {}
    for prot in all_matched_proteins:
        #only include peptide matches if peptide belongs to precursor that passed 1% Precursor-level FDR
        peps = set([x for x in proteins_peptides[prot] if x in peptide_precursor])
        if peps:
            proteins_peptides[prot] = peps
        else:
            del proteins_peptides[prot]
            continue
        protseq = Database[prot]
        if protseq not in Protein_sequences:
            Protein_sequences[protseq] = set()
        Protein_sequences[protseq].add(prot)
    del all_matched_proteins
    #print(len(proteins_peptides), 'proteins')
    Protein_sequences = {'P'+str(xi):Protein_sequences[x] for xi,x in enumerate(sorted(Protein_sequences.keys()))}
    #print(len(Protein_sequences), 'protein sequences')
    
    #Basically, the protein FDR is calculated by first remapping peptides to all proteins in the database. If a peptide matches to a target, all decoy proteins are removed. PSMs that 
    #failed PSM-level or precursor-level FDR are filtered out according to the FDR thresholds in the input form. The top protein per precursor and score per protein are then calculated 
    #according to this algorithm:
    #1) Set the score of each protein to the sum of -log10 EValues for all (not necessarily uniquely-mapped) precursors where score of precursor is the best PSM of a precursor.
    #2) Pick the highest scoring protein, remove the precursors mapping to that protein and recalculate scores.
    #3) Iterate until there are no precursors left.
    Protein_scores = []
    Protein_precursors = {}
    Protein_peptides = {}
    Protein_td = {}
    for pgroup in Protein_sequences:    
        score = 0
        precursors = set()
        prot = list(Protein_sequences[pgroup])[0]
        for pep in proteins_peptides[prot]:
            precursors.update(peptide_precursor[pep])
            precursor_scores = [Precursors_scores[x] for x in peptide_precursor[pep]]
            score+=sum(precursor_scores)
        Protein_precursors[pgroup] = precursors
        Protein_scores.append((score, pgroup))
        Protein_peptides[pgroup] = proteins_peptides[prot]
    
        td = set([decoy_prefix if x[:4]==decoy_prefix else '_' for x in Protein_sequences[pgroup]])
        if len(td)>1:
            td = '_'
            print(td)
        else:
            td = list(td)[0]
        Protein_td[pgroup] = td
    
    all_proteins = set(Protein_sequences.keys())
    all_precursors = set(Precursor_list.keys())
    covered_precursors = set()
    Protein_List = {}
    t1 = time.time()
    c=0
    while covered_precursors != all_precursors:
        
        if not Protein_scores:
            print('no Protein_scores left...')
            break
        
        c+=1
        top_scoring = max(Protein_scores, key=itemgetter(0))
        max_score = top_scoring[0]
        top_protg = top_scoring[1]    
        Protein_List[top_protg] = (top_protg, max_score, Protein_td[top_protg], Protein_precursors[top_protg]-covered_precursors)

        covered_precursors|=Protein_precursors[top_protg]

        Protein_scores = []
        for prot in all_proteins:
            if prot not in Protein_List:
                score = sum([Precursors_scores[x] for x in Protein_precursors[prot]-covered_precursors])
                if score!=0:
                    Protein_scores.append((score, prot))

        if len(covered_precursors)%5000==0:
            print(c, time.time() - t1, len(covered_precursors))

    print('all_precursors: ', len(all_precursors))
    print('covered_precursors: ', len(covered_precursors))
    
    
    ####### Protein-level FDR ###############
    FDR_calc_score_direction = 1
    FDRcompute_targetpass_col = 0
    FDRcompute_decoypass_col = 1
    FDRcompute_pvalue_col = 2
    FDRcompute_tdnum_col = 3
    FDRcompute_FDR_col = 4
    FDRcompute_FDR_values = 5
    
    scored_protein_list = Protein_List.values()
    #Compute_FDR(list, fdr, col['SpecEValue'], FDR_calc_score_direction, col['target_decoy'], decoy_prefix)
    ProteinlevelFDR = Compute_FDR(scored_protein_list, protein_level_fdr, 1, FDR_calc_score_direction, 2, decoy_prefix)
    if ProteinlevelFDR == 'cannot reach specified FDR':
        raise ValueError('Cannot reach specified '+str(int(protein_level_fdr*100))+'% level FDR.')
    
    print('### '+str(int(protein_level_fdr*100))+'% Protein-level FDR ####')
    print('p-value threshold: ', ProteinlevelFDR[FDRcompute_pvalue_col])
    print('decoy|target: ', ProteinlevelFDR[FDRcompute_tdnum_col])
    print('FDR: ', ProteinlevelFDR[FDRcompute_FDR_col])
    
    #PSMs
    with open(Proteinlevel_PSMs_txt,'wb') as outfile, open(Proteinlevel_proteins_txt,'wb') as outfile2:
        target_proteins = set()
        outfile2.write('\t'.join(['#Pg', 'Score', 'Proteins', 'Peptides', 'Prot_sequence'])+'\n')
        for entry in ProteinlevelFDR[FDRcompute_targetpass_col]:
            pg = entry[0]
            prot = list(Protein_sequences[pg])[0]
            score = str(entry[1])
            outfile2.write('\t'.join([pg, score, ';'.join(Protein_sequences[pg]), ';'.join(Protein_peptides[pg]), Database[prot]])+'\n')
            target_proteins.update(Protein_sequences[pg])
        print('target_proteins: ', len(target_proteins))
        
        header = ''
        PSM_t = 0
        with open(PSMlevel_txt,'rb') as infile:
            for line in infile:
                sp = line.strip().split('\t')
                if line[0]=='#':
                    header = sp[0:]
                    outfile.write('\t'.join(header)+'\n')
                    colidx = {x:xi for xi,x in enumerate(sp)}
                else:
                    precursor = sp[colidx['precursor_pepcharge']]
                    if precursor in Precursor_list:
                        proteins = set(sp[colidx['Protein']].split(';'))
                        if proteins-target_proteins != proteins:
                            PSM_t+=1
                            outfile.write('\t'.join(sp)+'\n')
        print('target PSMs', PSM_t)

    return None

def main():
    args = arguments()
    compute_fdr(args)
    
if __name__ == "__main__":
    main()