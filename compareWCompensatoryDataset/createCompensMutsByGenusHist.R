library(ggplot2)
library(dplyr)
library(tidyr)
pfamCompensMuts = read.csv("/home/benlitterer/Academic/Research/ProjectP53/compareWCompensatoryDataset/pfamCompensMutsWSeqSim.WithOrg.tsv", sep="\t",na.strings=c(""," ","NA"))
pfamCompensMuts = separate(pfamCompensMuts, "organism", c("org_first_word"), remove=FALSE, extra="drop")

#add a frequency column to the table 
org_freqs = data.frame(sort(table(pfamCompensMuts$org_first_word), decreasing = TRUE))
names(org_freqs) = c("org_first_word", "freq")
just_over_ten = pfamCompensMuts %>% 
  inner_join(org_freqs, by=c("org_first_word")) %>% 
  arrange(desc(freq)) %>%
  filter(freq>5)
  
just_over_ten = select(just_over_ten, "H168R","T284R","N239Y","N268D","T123A","T123P","S240N","N239F","H233Y","N235K","N239W","S240R","sequence_similarity","residue_similarity", "org_first_word", "prot_id")
just_over_ten$org_first_word = as.factor(just_over_ten$org_first_word)
ggplot(just_over_ten, aes(org_first_word)) + geom_bar()

#get the relevant cols to be logical vectors 
just_over_ten[,1:12] = lapply(just_over_ten[,1:12], as.logical.factor)


just_over_ten = gather(just_over_ten, "mutation", "mut_bool", 1:12)
just_over_ten = filter(just_over_ten, mut_bool==TRUE)

ggplot(just_over_ten, aes(mutation, fill=org_first_word)) + geom_bar() + facet_wrap(~org_first_word) + theme(axis.text.x = element_text(angle = 90)) + labs(title="compensatory mutations by organsim",fill="genus")

#this outputs the proteins used to create the graphic as strings with a space inbetween 
#easy to copy and use for other things
cat(as.character(arrange(just_over_ten, desc(residue_similarity))))