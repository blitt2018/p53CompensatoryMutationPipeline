library(ggplot2)
library(dplyr)
library(tidyr)
pfamCompensMuts = read.csv("/home/benlitterer/Academic/Research/ProjectP53/compareWCompensatoryDataset/pfamCompensMutsWSeqSim.WithOrg.tsv", sep="\t",na.strings=c(""," ","NA"))
pfamCompensMuts = separate(pfamCompensMuts, "organism", c("org_first_word"), remove=FALSE, extra="drop")

pfamCompensMuts$org_first_word = as.factor(pfamCompensMuts$org_first_word)

justDels = pfamCompensMuts[,c(2:13, 31)] 

#get the relevant cols to be logical vectors 
justDels[,1:12] = lapply(justDels[,1:12], as.logical.factor)

#for verification of the next step 
#get the sum of the columns 
#justDels[,1:12] = lapply(justDels[,1:12], sum)
justDels = gather(justDels, "mutation", "mut_bool", 1:12)
justDels = filter(justDels, mut_bool == TRUE)

ggplot(justDels, aes(org_first_word)) + geom_bar()

compensWithDels = pfamCompensMuts[,c(2:13,15:26, 31)] 

compensWithDels[,1:24] = lapply(compensWithDels[,1:24], as.logical.factor)

#vector of delitirious mutations from original table 
delMuts = names(pfamCompensMuts)[2:13]
compensWithDels$delCount = rowSums(compensWithDels[,2:13])

#shows that there are actually counts in the delCount column 
#sum(compensWithDels$delCount)

#for verification of the next step 
#get the sum of the columns 
#compensWithDels[,1:24] = lapply(compensWithDels[,1:24], sum)

compensWithDels = gather(compensWithDels, "mutation", "mut_bool", 1:24)
compensWithDels = mutate(compensWithDels, isDel = ifelse(mutation %in% delMuts, TRUE, FALSE))
compensWithDels = filter(compensWithDels, delCount > 0 & mut_bool == TRUE)

ggplot(compensWithDels, aes(mutation, fill=isDel)) + geom_bar() + facet_wrap(~org_first_word) + theme(axis.text.x = element_text(angle = 90)) + labs(title="compensatory mutations by organism genus",fill="delitirious?") + scale_fill_manual(values=c("#95d5e6", "#f25757"))
