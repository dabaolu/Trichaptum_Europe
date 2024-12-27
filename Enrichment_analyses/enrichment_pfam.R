#R script to 
#1) extract pfam of genes within a given genomic window
#2) use fischers exact test to see if the pfam occurence in a given subset (SNPs or genomic windows)
#are significantly different from that of the whole genome

setwd("/Users/dabaosl/UiO Dropbox/Dabao Lu/UiO/Phd/Data_analysis_global")
setwd("/Users/dabaosl/UiO Dropbox/Dabao Lu/UiO/Phd/Europe_study/Enrichment_analyses")

library(tidyverse)
library(stringr)
library(stringi)
library(stats)

genome <- read.table("TA10106M1_final_MAKER.putative_function_domains.gff", sep="\t")
str(genome) #223110 obs, a lot of these are duplicated
genome.genes <- filter(genome, V3=="gene") 
str(genome.genes) #9298 obs, should represent unique genes
head(genome.genes)
pfam <- grep("Pfam", genome.genes$V9, fixed = TRUE) 
str(pfam) #4422 genes with pfam
head(genome.genes, 15)
genome.genes.pfam <- genome.genes[pfam,] #get genes with pfam annotation
colnames(genome.genes.pfam) <- c("seqid","source","type","start","end","score", "strand","phase","attributes") 
str(genome.genes.pfam)

#some genes have multiple pfams, extract these and put in new columns:
str_extract_all(genome.genes.pfam[1,9],"Pfam:PF[0-9]{5}")
pfams <- str_extract_all(genome.genes.pfam$attributes,"Pfam:PF[0-9]{5}") #extracts pfams with regular expression
length(pfams) #matches number of genes with pfams (4422)
lengths(pfams) #number of pfams in each gene
max(lengths(pfams)) #13
pfam.df <- as.data.frame(t(stri_list2matrix(pfams))) #convert to dataframe, each pfam for each gene is a column
str(pfam.df) 
genome.genes.pfam.split <- cbind(genome.genes.pfam, pfam.df)
str(genome.genes.pfam.split)
dim(genome.genes.pfam.split)

#gather and tally pfams on second half of scaffold 5:
genome.genes.pfam.split.scf5.2nd <- filter(genome.genes.pfam.split, seqid=="Scaffold05" & end >= 3213128)
str(genome.genes.pfam.split.scf5.2nd) #390 genes with pfams, some of these have multiple pfams


#make a count of pfams:
#1) all genes (with pfam):
n_distinct(genome.genes.pfam.split$V1) #1675 unique pfams
sapply(genome.genes.pfam.split, function(x) n_distinct(x))
genome.pfams <- c(genome.genes.pfam.split$V1,
                  genome.genes.pfam.split$V2,
                  genome.genes.pfam.split$V3,
                  genome.genes.pfam.split$V4,
                  genome.genes.pfam.split$V5,
                  genome.genes.pfam.split$V6,
                  genome.genes.pfam.split$V7,
                  genome.genes.pfam.split$V8,
                  genome.genes.pfam.split$V9,
                  genome.genes.pfam.split$V10,
                  genome.genes.pfam.split$V11,
                  genome.genes.pfam.split$V12,
                  genome.genes.pfam.split$V13)
str(genome.pfams) 
length(genome.pfams) #57486 pfams in genome
genome.pfams.counts <- table(genome.pfams)
str(genome.pfams.counts)
genome.pfams.counts <- as.data.frame(genome.pfams.counts)
colnames(genome.pfams.counts) <- c("pfam","freq_all")
str(genome.pfams.counts) #2680 unique pfams in total 4422 genes

#2) scaffold5 2nd half genes with pfam:
selected.genes.pfams <- c(genome.genes.pfam.split.scf5.2nd$V1,
                          genome.genes.pfam.split.scf5.2nd$V2,
                          genome.genes.pfam.split.scf5.2nd$V3,
                          genome.genes.pfam.split.scf5.2nd$V4,
                          genome.genes.pfam.split.scf5.2nd$V5,
                          genome.genes.pfam.split.scf5.2nd$V6,
                          genome.genes.pfam.split.scf5.2nd$V7,
                          genome.genes.pfam.split.scf5.2nd$V8,
                          genome.genes.pfam.split.scf5.2nd$V9,
                          genome.genes.pfam.split.scf5.2nd$V10,
                          genome.genes.pfam.split.scf5.2nd$V11)
length(selected.genes.pfams) #4290 pfams in this region of the scaffold
selected.genes.pfams.counts <- table(selected.genes.pfams)
selected.genes.pfams.counts <- as.data.frame(selected.genes.pfams.counts)
colnames(selected.genes.pfams.counts) <- c("pfam", "freq_selected")
str(selected.genes.pfams.counts) #331 pfams in 390 pfam genes on scaffold5
head(selected.genes.pfams.counts)

#combine 1) and 2):
combined.table.all <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam", all.x = TRUE)
dim(combined.table.all) #2680 pfams
head(combined.table.all,20)
combined.table.all[is.na(combined.table.all)] = 0 #replace NA with 0
head(combined.table.all,20)
fisher.table <- combined.table.all[,2:3]
rownames(fisher.table) <- combined.table.all$pfam
head(fisher.table)
fisher.table <- t(fisher.table)
dim(fisher.table)
test <- fisher.test(fisher.table, simulate.p.value=TRUE) #must use simulate to not exceed memory limitation
test #p-value=1, i.e. no statistically significant difference between two sets, that is whole genome and 2nd half of scaffold 5

###Do fischers exact test for individual pfams:
#two classes: 1)pfam all genes and 2)pfam genes with snps
#occurences of each pfam in each class

#Actually classes should be mutually exclusive so they can add up to total, so what you want is:
#1) all pfam in snp genes VERSUS all pfam not in snp genes
#2) all pfam outside snp genes VERSUS all pfam not outside snp genes

combined.2ndscf5.table <- merge(genome.pfams.counts, selected.genes.pfams.counts, by="pfam")
dim(combined.2ndscf5.table) #331: i.e. omits pfams not present in 2nd half of scaffold 5
head(combined.2ndscf5.table)
not.2ndscf5.table <- subset(combined.table.all, !(combined.table.all$pfam %in% combined.2ndscf5.table$pfam))
dim(not.2ndscf5.table) #2349 pfams not in 2nd half scaffold 5, i.e. only found outside this region
#pfam not found outside 2nd hald scf5:
only.2ndscf5.table <- filter(combined.2ndscf5.table, freq_all == freq_selected)
dim(only.2ndscf5.table) #72 pfams only in 2nd half scaffold 5, i.e. 72 pfams absent outside 2nd half scaffold 5.
table(genome.genes.pfam.split$V1 %in% only.2ndscf5.table$pfam) #38, does this equal to the number of genes?
table(genome.genes.pfam.split$V2 %in% only.2ndscf5.table$pfam) #19
table(genome.genes.pfam.split$V3 %in% only.2ndscf5.table$pfam) #10
table(genome.genes.pfam.split$V4 %in% only.2ndscf5.table$pfam) #6
table(genome.genes.pfam.split$V5 %in% only.2ndscf5.table$pfam) #2
table(genome.genes.pfam.split$V6 %in% only.2ndscf5.table$pfam) #0
table(genome.genes.pfam.split$V7 %in% only.2ndscf5.table$pfam) #0
table(genome.genes.pfam.split$V8 %in% only.2ndscf5.table$pfam) #0

table(not.2ndscf5.table$pfam %in% only.2ndscf5.table$pfam) #appears to not be overlapping
table(only.2ndscf5.table$pfam  %in% not.2ndscf5.table$pfam) #appears to not be overlapping

#Try to test pfam by pfam with counts, and then correct for multiple testing:
head(combined.table.all,3)
#counts pfam present outside 2ndscf5 genes: freq_all - freq_selected
#counts pfam present inside 2ndscf5 genes: freq_selected
#counts pfams absent outside 2ndscf5 genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside 2ndscf5 genes: 4422 - freq_selected

#counts pfam present outside snp genes: freq_all - freq_selected
#counts pfam present inside snp genes: freq_selected
#counts pfams absent outside snp genes: 4422 - (freq_all - freq_selected)
#counts pfams absent inside snp genes: 4422 - freq_selected

pvalues <- c() #loop through each pfam and do a Fischers exact test:
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within <- freq_selected
  present_outside <- freq_all-freq_selected
  absent_within <- 390-freq_selected
  absent_outside <- 4422-390-present_outside
  genes_pfam_present <- c(present_within, present_outside)
  genes_pfam_absent <- c(absent_within, absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present, genes_pfam_absent))
  #cont_table.list <- c(cont_table, cont_table.list)
  rownames(cont_table) <- c("within_gene","outside_gene")
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues) #2680 (=number of unique pfams)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues) #combine pvalues with pfam information
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),] #sort with smalles p-value first
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction

#test instead subset of snp genes versus all genes (but does this violate the assumption of exclusivity of categories?)
#replace "outside_snp_genes" where pfam is absent with all genes (4422)
# but shouldn´t "outside_snp_genes" where pfam is present also be replaced with all genes?
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within <- freq_selected
  present_outside <- freq_all-freq_selected
  absent_within <- 390-freq_selected
  absent_outside <- 4422-present_outside
  genes_pfam_present <- c(present_within,present_outside)
  genes_pfam_absent <- c(absent_within,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 pfams significant

#go all the way by also replacing "outside_snp_gene" for pfam present with counts of total:
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within <- freq_selected
  present_all <- freq_all
  absent_within <- 390-freq_selected
  absent_all <- 4422-freq_all
  genes_pfam_present <- c(present_within, present_all)
  genes_pfam_absent <- c(absent_within, absent_all)
  cont_table <- as.data.frame(cbind(genes_pfam_present, genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","all_genes")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}

length(pvalues)
cont_table
combined.table.all$freq_outside <- combined.table.all$freq_all - combined.table.all$freq_selected
final.table <- cbind(combined.table.all, pvalues)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH") #make adjustment for multiple testing
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0 significant after correction

#exclude pfams that are not present in 2nd half of scf5 genes? (will make BH correction softer as fewer tests are conducted)
#but this is probably "unfair"...#well, chatgpt seems to think that it is ok..

pfam.genes <- filter(combined.table.all, freq_selected > 0)
dim(pfam.genes) #331 pfams present in 185 genes
head(pfam.genes, 5)
a <- filter(pfam.genes, freq_selected == freq_all)
dim(a) #72 pfams present only in scf5 2nd half genes
a <- filter(pfam.genes, freq_selected < freq_all)
dim(a) #259 pfams shared scf5 2nd half and genes outside scf5 2nd half
str(genome.genes.pfam.split)
b <- filter(genome.genes.pfam.split, genome.genes.pfam.split$V1 %in% a$pfam)
dim(b)#2047 genes with shared pfam occuring in total, i.e. 2047-390=1657 genes outside 2nd half of scf5 carrying shared pfams with that region?

pvalues <- c()
for(i in 1:nrow(pfam.genes)){
  freq_all <- pfam.genes[i,2]
  freq_selected <- pfam.genes[i,3]
  present_within <- freq_selected
  present_outside <- freq_all-freq_selected
  absent_within <- 390-freq_selected
  absent_outside <- 1657-present_outside
  genes_pfam_present <- c(present_within,present_outside)
  genes_pfam_absent <- c(absent_within,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_delta_gene","outside_delta_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
  print(cont_table)
}

length(pvalues)
cont_table
pfam.genes$freq_outside <- pfam.genes$freq_all - pfam.genes$freq_selected
final.table <- cbind(pfam.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #1 pfams significant after correction: PF13881-ubiquitin domain found in ubiquitin-like protein

#try with pfam counts instead of gene counts:
sum(combined.table.all$freq_all) #7167 total counts of pfams (including duplicates, i.e. this are not necessarily unique)
sum(combined.table.all$freq_selected) #583 total counts of pfams in scf5 2nd half (including duplicates, i.e. these are not necessarily unique)
pvalues <- c()
for(i in 1:nrow(combined.table.all)){
  freq_all <- combined.table.all[i,2]
  freq_selected <- combined.table.all[i,3]
  present_within <- freq_selected
  present_outside <- freq_all-freq_selected
  absent_within <- 583-freq_selected
  absent_outside <- 7167-583-present_within-present_outside
  genes_pfam_present <- c(present_within,present_outside)
  genes_pfam_absent <- c(absent_within,absent_outside)
  cont_table <- as.data.frame(cbind(genes_pfam_present,genes_pfam_absent))
  rownames(cont_table) <- c("within_snp_gene","outside_snp_gene")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.genes$freq_outside <- pfam.genes$freq_all - pfam.genes$freq_selected
final.table <- cbind(combined.table.all, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #0


#try counts of pfams that only occur in 2nd half of scf5 (i.e. ecluding pfams that are not present in 2nd half of scf5)
sum(pfam.genes$freq_all) #2913 counts of pfams that are present both in background genome and scaffold 5 2nd half
sum(pfam.genes$freq_selected) #583 total counts of pfams in scf5 2nd half
pvalues <- c()
for(i in 1:nrow(pfam.genes)){
  freq_all <- pfam.genes[i,2]
  freq_selected <- pfam.genes[i,3]
  present_within <- freq_selected
  present_outside <- freq_all-freq_selected
  absent_within <- 583-freq_selected
  absent_outside <- 2913-583-present_within-present_outside
  pfam_present <- c(present_within, present_outside)
  pfam_absent <- c(absent_within, absent_outside)
  cont_table <- as.data.frame(cbind(pfam_present, pfam_absent))
  rownames(cont_table) <- c("within_pfam","outside_pfam")
  #cont_table.list <- c(cont_table, cont_table.list)
  test <- fisher.test(cont_table)
  p <- test$p.value
  pvalues <- c(p, pvalues)
}
length(pvalues)
cont_table
pfam.genes$freq_outside <- pfam.genes$freq_all - pfam.genes$freq_selected
final.table <- cbind(pfam.genes, pvalues)
head(final.table)
sorted.final.table <- final.table[order(final.table$pvalues),]
head(sorted.final.table,20)
sorted.final.table$BH <- p.adjust(sorted.final.table$pvalues, method = "BH")
head(sorted.final.table,20)
filter(sorted.final.table, BH < 0.05) #1 pfams significant after correction: PF13881-ubiquitin domain found in ubiquitin-like protein
