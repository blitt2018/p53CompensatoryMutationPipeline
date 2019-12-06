#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:37:53 2019

@author: benlitterer
"""

import requests 
import pandas as pd
import sys 

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

infile = sys.argv[1]

mutationDf = pd.read_csv(infile, "\t", index_col=0)

#retrieves the NCBI protein accesssion id's from the input file (dataframe)
#THIS ASSUMES THAT THE ACCESSION IDS COLUMN IS 3 FROM THE RIGHT 
prot_accessions = [prot.split("/")[0] for prot in mutationDf["prot_id"]]

#simply the weirdest syntax known to man 
id_query_str = "+OR+".join(prot_accessions)
print(id_query_str)

# example: https://www.uniprot.org/uniprot/?query=A0A0L7RJB6+OR+A0A182Q5M9+OR+Q80ZA1+OR+A0A1I7VBH3&sort=score
#uses entrez etools api to submit a batch query requesting xml sequence output 
response = requests.get("https://www.uniprot.org/uniprot/?query=" + id_query_str + "&format=tab&columns=id,organism")

if response.status_code == 200: 
    print("entrez api call successful")
else: 
    print("entrez api call unsuccessful")
    
#split the output into lines and don't include the header 
lines = response.text.split("\n")[1:]

#create dictionary to hold key value pairs where the key is the protein id and the value is the organism name 
pfamOrgDict = {}

for line in lines: 
    if line.strip() != "": 
        prot_id = line.split("\t")[0]
        organism = line.split("\t")[1]
        pfamOrgDict[prot_id] = organism 

orgDf = pd.DataFrame.from_dict(pfamOrgDict, columns=["organism"], orient="index")

#we need to merge the dataframe based on a column, not the index, so make the index a column 
orgDf["prot_id_stripped"] = orgDf.index

#uniprot sends back a stripped version of the protein, we need to match this in the original dataset 
mutationDf["prot_id_stripped"] = [prot_id.split("_")[0] for prot_id in mutationDf["prot_id"]]

#add the organism information onto the input dataframe 
outDf = pd.merge(mutationDf, orgDf, on="prot_id_stripped")


outfile_path = infile.split(".")[0] + ".WithOrg.tsv" 

outDf.to_csv(outfile_path, sep = "\t", header=True, index=None)