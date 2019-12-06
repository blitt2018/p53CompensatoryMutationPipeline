#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 13:55:57 2019

@author: benlitterer
"""
#This file uses a tsv file with mutations to look for as well as the full pfam P53 set of sequences. It outputs a tsv file named after the first command line argument containing the mutations found. 

from reindexForHMMERAlignment import translate_indeces
import pandas as pd 
import sys
from Bio import SeqIO

pd.set_option('display.max_columns', 25)
pd.set_option('display.max_rows', 20)

experimentalCompensMuts = "/home/benlitterer/Academic/Research/ProjectP53/compareWCompensatoryDataset/experimentallyValidatedCompensMuts.tsv"
pfamP53File = "/home/benlitterer/Academic/Research/ProjectP53/knownCompensMuts/PF00870_full.faa"

#Need to find which indeces used to refer to mutation positions correspond with the indeces on our aligned proteins in pfam 
FULL_REF_STR = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK\
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE\
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS\
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP\
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG\
GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"

#the amount of the referance sequence that is aligned in pfam 
REF_STR = "SQKTYQGSYGFRLGFLHSGTAK\
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE\
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS\
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL"

ACTUAL_STR = "............................................................\
....s-QKTYQGSYGFRLGFLH-.....SGTAKSVTCT......................\
................................................YSP..AL..NKM\
FCQ.LA.K.TCPVQLWVDSTP.........................PPGTRVRAMA.I.Y\
KQS........Q.HM.......T.......E.......VVR......R.....CPHHERC\
SDSDGL..............APPQH....L.IRVEGN...lRVEYL.-DDRNT--..-FR\
HSV....................VVPYEPPE...............vGSD......CTTI\
..............HY............................................\
....NYMCNSS..C..M...G....G..MNR.R.P.ILT.IITL.EDS...SG.......\
.........N.LLGRNSF.EVRVCACPGRDRRTEEEN-l.....................\
.......................................".upper()

translated_dict = translate_indeces(REF_STR, ACTUAL_STR)

#this shows that the translation was done successfully, but we need to remember to subtract one from what is in the dict (since 1 based indexing is used in publications etc...)
for key, val in translated_dict.items(): 
    if not REF_STR[key-1] == ACTUAL_STR[translated_dict[key]-1]: 
        print("translation unsuccessful")

#the amount that we need to shift our indices for the referance string due to not having an alignment for the full sequence (Ref starts at residue indexed at 98)
SHIFT_FACTOR = 99

#get the mutations that we want to find in pfam sequences 
sur_df = pd.read_csv(experimentalCompensMuts, "\t")
pfamP53Alignments = SeqIO.parse(pfamP53File, "fasta")


#print(REF_STR[239-SHIFT_FACTOR])
#print(ACTUAL_STR[translated_dict[239-99]])
    
indeces_in_question = [227, 228, 229, 230, 232, 233, 234, 235, 236, 239, 240]
        
DEL_MUTS = []
COMPENS_MUTS = []

for index, row in sur_df.iterrows(): 
    del_mut = row[0].strip("\n ")
    
    #only want unique values in our set (will be used as columns in dataframe eventually)
    if del_mut not in DEL_MUTS and del_mut != "NA" and del_mut != "nan": 
        DEL_MUTS.append(del_mut)
    compens_muts = [item.strip("\n ") for item in row[1].split(",")]
    for mut in compens_muts: 
        if mut not in COMPENS_MUTS and del_mut != "NA" and del_mut != "nan": 
            COMPENS_MUTS.append(mut)

print(DEL_MUTS)
print(COMPENS_MUTS)

#this takes a sequence and returns whether there is a delitirious mutation at each index of interest 
def check_del_muts(seq):
    del_out = []
    for del_str in DEL_MUTS: 
        del_mut = del_str[-1]
        del_index = int(del_str[1:-1])
        
        #here we want to find whether the residue at the location where we want to find a delitirious mutation is actually delitirious 
        if seq[translated_dict[del_index-SHIFT_FACTOR]] == del_mut: 
            del_out.append(True) 
        else: 
            del_out.append(False)
    return del_out 

#this takes a sequence and returns whether there is a compensatory mutation at each index of interest 
def check_compens_muts(seq): 
    compens_out = []
    for compens_str in COMPENS_MUTS: 
        compens_mut = compens_str[-1]
        compens_index = int(compens_str[1:-1])
        
        #here we want to find whether the residue at the location where we want to find a compensatory mutation is actually compensatory
        if seq[translated_dict[compens_index-SHIFT_FACTOR]] == compens_mut: 
            compens_out.append(True) 
        else: 
            compens_out.append(False)
    return compens_out

#column names for our output dataframe. Includes empty column to seperate delitirious mutations from compensatory ones
names = ["prot_id"] + DEL_MUTS + ["<<<<delitirious   compensatory>>>>"] + COMPENS_MUTS + ["sequence_similarity", "residue_similarity"]
out_nested_list = []

#this gives a similarity score based on whether indeces have the same character 
def get_similarity(seqOne, seqTwo):
    seqOne = seqOne.upper()
    seqTwo = seqTwo.upper()
    sameCount = 0
    longerSeqLen = max([len(seqOne), len(seqTwo)])
    for resNum in range(longerSeqLen):
        firstRes = seqOne[resNum]
        secRes = seqTwo[resNum]
        if firstRes == secRes or ((firstRes == "-" and secRes == ".") or (firstRes == "." and secRes == "-")): 
            sameCount += 1
    percentSim = sameCount/longerSeqLen
    return(percentSim)
    
#this gives a similarity score based on whether indeces where at least one sequence actually has a residue are the same 
#the function above inflates similarity by counting indeces where there are gaps in both sequences being compared 
def get_residue_similarity(seqOne, seqTwo):
    seqOne = seqOne.upper()
    seqTwo = seqTwo.upper()
    sameCount = 0
    hasResidues = 0
    longerSeqLen = max([len(seqOne), len(seqTwo)])
    for resNum in range(longerSeqLen):
        firstRes = seqOne[resNum]
        secRes = seqTwo[resNum]
        if firstRes == secRes:
            #we only care about similarity when there is a residue in play, since we know they're the same, we only need to check if the first residue is a gap character 
            if firstRes != "-" and firstRes != ".":  
                sameCount += 1
                
        #if either the first res or the second has a residue that isn't a gap, we can add to our total amount of positions with residues 
        if (firstRes != "-" and firstRes != ".") or (secRes != "-" and secRes != "."): 
            hasResidues += 1
                
    percentSim = sameCount/hasResidues
    return(percentSim)
    
#for each sequence record in pfam P53 profile 
for seqRec in pfamP53Alignments: 
    pfamSeq = seqRec.seq
    pfamName = seqRec.name
    
    #get compensatory and delitirous boolean information for each index in question 
    del_list = check_del_muts(pfamSeq) 
    compens_list = check_compens_muts(pfamSeq)
    
    #only use prots that have compensatory muts in the output dataframe. Tacks on sequence similarity to the end
    if True in compens_list: 
        out_list = [pfamName] + del_list + ["------"] + compens_list + [get_similarity(pfamSeq, ACTUAL_STR), get_residue_similarity(pfamSeq, ACTUAL_STR)]
        out_nested_list.append(out_list)
 

out_df = pd.DataFrame(out_nested_list, columns=names)

#the first argument is used as the output location 
out_df.to_csv(sys.argv[1], "\t")