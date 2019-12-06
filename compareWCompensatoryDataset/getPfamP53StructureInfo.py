#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:35:26 2019

@author: benlitterer
"""

import requests 
import pandas as pd
import sys 
import numpy as np

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

infile = sys.argv[1]

mutationDf = pd.read_csv(infile, "\t")

#retrieves the NCBI protein accesssion id's from the input file (dataframe)
#THIS ASSUMES THAT THE ACCESSION IDS COLUMN IS 3 FROM THE RIGHT 
prot_accessions = [prot.split("/")[0] for prot in mutationDf["prot_id"]]


#simply the weirdest syntax known to man 
id_query_str = "+OR+".join(prot_accessions)
print(id_query_str)

# example: https://www.uniprot.org/uniprot/?query=A0A0L7RJB6+OR+A0A182Q5M9+OR+Q80ZA1+OR+A0A1I7VBH3&sort=score
#uses entrez etools api to submit a batch query requesting xml sequence output 
response = requests.get("https://www.uniprot.org/uniprot/?query=" + id_query_str + "&format=tab&columns=id,database(PDBe),database(SMR),database(ModBase),database(PDBe-KB),3d")

print(response.status_code)
if response.status_code == 200: 
    print("entrez api call successful")
else: 
    print("entrez api call unsuccessful")
    
#split the output into lines and don't include the header 
lines = response.text.split("\n")[1:]


#create dictionary to hold key value pairs where the key is the protein id and the value is the database value
pfamOrgDict = {}

for line in lines: 
    if line.strip() != "": 
        prot_id = line.split("\t")[0]
        db_hits = line.split("\t")[1:]
        pfamOrgDict[prot_id] = db_hits

orgDf = pd.DataFrame.from_dict(pfamOrgDict, columns=["PDBe", "SMR", "ModBase", "PDBe-KB","3d"], orient="index")
orgDf = orgDf.replace(r'^\s*$', np.nan, regex=True)
orgDf = orgDf.dropna(how="all")
orgDf = orgDf.dropna(axis="columns")
orgDf["prot_id_stripped"] = orgDf.index 
orgDf = pd.merge(orgDf, mutationDf, on="prot_id_stripped", how="left")

orgDf.to_csv(infile.split(".")[0] + ".WithStructs.tsv", sep="\t", header=True, index=None)
