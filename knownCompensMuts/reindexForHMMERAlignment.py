#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:52:30 2019

@author: benlitterer
"""
from Bio import SeqIO
import sys 
import itertools
import pandas as pd

"TODO: problem with counting, not sure what it is" 

"""
HMMER_file = open(sys.argv[1],"r") 
HMMERParserObject = SeqIO.parse(HMMER_file, "fasta")

all_org_file = open(sys.argv[2], "r")
OrgParserObject = SeqIO.parse(all_org_file, "fasta")
"""

def translate_indeces(ref_seq, new_seq): 
    ref_res_counter = 0
    new_res_counter = 0
    translation_dict = {}
    
    #consider that we can have dash or period
    #check while loop conditions? 
    while ref_res_counter < len(ref_seq) and new_res_counter < len(new_seq): 
        ref_res = ref_seq[ref_res_counter]
        new_res = new_seq[new_res_counter]
        
        #print("ref_res vals:")
        while (ref_res != "-" and ref_res != "." ) and (new_res == "-" or new_res == "."): 
            new_res_counter += 1
            #print(new_res_counter)
            new_res = new_seq[new_res_counter]
         
        #print("new_res vals:")
        while (ref_res == "-" or ref_res == ".") and (new_res != "-" and new_res != "."):
            ref_res_counter += 1
            #print(ref_res_counter)
            ref_res = ref_seq[ref_res_counter]
        
        
        
        """
        for checking that it works
        print(ref_res_counter)
        print(new_res_counter)
        print("-----------------------------")
        """
        if ref_res != new_res:
            
            #this says: if it is not the case that both ref_res and new_res contain gap characters
            if not ((ref_res == "." or ref_res == "-") and (new_res == "." or new_res == "-")): 
                raise Exception("Index Translation has failed")
            
        translation_dict[ref_res_counter+1] = new_res_counter+1
        new_res_counter += 1
        ref_res_counter += 1
    return translation_dict


ref_prac = "..SKYPGPFNFQLQFIDDS.....GGKPVPRATPFT......................................................................YSP..QA...QRLFVQ.MN.H.KVPFRFSISDRV..............\
...........PNGTWL.RIR..T.R.FTQPE.FR.......N...........................................................................................................................\
......................................................................................................................................................................\
..............................................." 

experi_string = "...........................................................................................................................\
......................................................................................................................................................................\
..............................................."

sexperi_string = "............................................................................................................................\
......................................................................................................................................................................\
....................................."

new_prac = "..........................................................................................................................................................\
............S---KYPGPFNFQLQFIDDS..GGKPVPRATPFT......................................................................YSPQA...QRLFVQMN.H.KVPFRFSISDRV...................\
......PNGTWL.RIR..T.R.FTQPE.FR.......-.......-.......---......--------------..............----..-....-.--.----....-----..--------..------....................---..----\
-................---......----..............--................................................----.---..-..-...-....-..---.....-.-.---.-..---.---..--.................\
........-----.....---.-------------------N............................................................................................................................\
......................................................................................................................................................................\
....................................."
"""
simp_ref = "..SKYPGPF..-NF-.QLQF-..-"
simp_prac = "..SKY..--PG-PF..-N..F-.QLQF-..-"
"""
output = translate_indeces(ref_prac, new_prac)
            
print(output)
