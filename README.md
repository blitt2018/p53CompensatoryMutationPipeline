# p53CompensatoryMutationPipeline

This repository contains the pipeline used to find compensatory mutations in the pfam P53 profile. Finding these mutations is useful in validating computational methods for compensatory mutation discovery. Ideally, organisms containing compensatory mutations could be used as a model to understand these same mutations in humans. The pipeline is composed of the following steps: 

## Gathering compensatory mutations 

Mutations were manually collected from a variety of papers and put into a spreadsheet containing the residue location, type of compensation, and placement in a higher order grouping of mutations. This file is not provided in the repository. 

## Finding compensatory mutations in the pfam P53 profile 

 - Pfam is a database that uses Hidden Markov Models to aggregate protein sequences with similar characteristics. The full P53 protein family was used in this case.
 - In order to find mutations in sequences other than the one they were originally located on, the corresponding residue locations on the aligned pfam sequences were found using a function in reindexForHMMERAlignment.py called translate_indeces. 
 - Next, each protein sequence was searched for both compensatory and delitirious mutations. Sequence similarity data was computed as well. See findCompensInPfam.py. 
 - Finally, Organism and Structural Data was added using the Uniprot REST API. This was done at seperate times and is broken into the files addOrgsToPfamCompensMutProts.py, and getPfamP53StructureInfo.py
 
## Analysis of compensatory mutations 
 - One question is whether or not certain biological groups have similar compensating mutations, and if so, what mutations occur in which groups? This question was answered visually in createCompensMutsByGenusHist.R. 
 - Additionally, it is useful to know if delitirious mutations in the human sequence occur in wild type P53 sequences for different organisms. This question was answered visually in plotDelMutsForOrgCompensMuts.R. 
 
## Result: 
 - An experimentally validated human compensatory mutation was found in a wild-type mouse sequence. This sequence has very high sequence identity to the human referance sequence. Work for this step was done in pymol and evidence was shown in mouse_compensatory_mutation_visualized.png. The mutation found is validated to stabilize the human p53 structure in experimental studies. 