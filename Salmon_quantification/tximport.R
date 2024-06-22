''' R v3.4'''
#Tximport to DEseq to GO analysis 
#Tximport is a step to import quant.sf file to counts 
#following the tutorial of Sambomics https://github.com/mousepixels/sanbomics_scripts/blob/main/salmon_to_deseq.Rmd

rm(list=ls()) # clear all variables and graphics
graphics.off()

####download and install required packages###

BiocManager::install("tximport")
library(tximport)

#######tximport#####
#setwd
setwd('/Users/baiprongfah/SALMON/AtRTD')

#obtain the genome annotation file GTF AtRTD2_19April2016.gtf from https://ics.hutton.ac.uk/atRTD/
#make sure to have the gtf file in the working directory
gtf <- rtracklayer::import('AtRTD2_19April2016.gtf')
gtf = as.data.frame(gtf) 
head(gtf)

#we need 2 dim vector of transcrip id and gene id 
tx2gene <- gtf[, c('transcript_id','gene_id')]


#####extract quant files#####
#1 assign the directory
quant_dir <- "sal_quant_all/" #sal_quant_all contains the quant.sf files of all samples in one folder
#2 create vector with the quant files 
quant_file <- list.files(quant_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
quant_file
#3 vector of sample names
sample_names <- c('col0_1', 'col0_2', 'sr45_1', 'sr45_2')
#4 combine 2 vectors
names(quant_file) <- sample_names
quant_file


##tximport
txi <- tximport(quant_file, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = TRUE)
txim <- as.data.frame(txi)
#save file
write.csv(txi, 'salmon_quant_count.csv', row.names=TRUE)
head(txim)

