''''R 4.4

rm(list=ls()) # clear all variables and graphics
graphics.off()

##load the package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)


packageVersion('IsoformSwitchAnalyzeR')
#>[1] ‘2.4.0’

setwd('/Users/baiprongfah/SALMON/AtRTD')

###import data
salmonQuant <- importIsoformExpression(parentDir = "sal_quant_all_2/")

#normalised abundance and count using edgeR
head(salmonQuant$abundance, 2)
head(salmonQuant$counts, 2)

myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)
myDesign

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "AtRTDv2_QUASI_19April2016.gtf" ,            
  isoformNtFasta       = "AtRTD2_19April2016.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE
)
summary(aSwitchList)

###run the first part of the isoform switch analysis workflow 
##0: add ORF from GTF file
a <- addORFfromGTF(aSwitchList, "AtRTD2_19April2016.gtf")
##1:filters for non-expressed genes/isoforms
filtered <- preFilter(
  switchAnalyzeRlist = a,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)
#>The filtering removed 22512 ( 38.2% of ) transcripts. There is now 36417 isoforms left

##2: Identify isoforn switches using DEXSeq
dexseq <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = filtered,
  reduceToSwitchingGenes=TRUE)
#>Step 1 of 2: Testing each pairwise comparisons with DEXSeq (this might be a bit slow)...
#>Estimated run time is: 15 min
#>Step 2 of 2: Integrating result into switchAnalyzeRlist...
#>Isoform switch analysis was performed for 10757 gene comparisons (100%).
#>Total runtime: 1.45 min
#>Done

extractSwitchSummary(dexseq)
#>Comparison nrIsoforms nrSwitches nrGenes
#>1 col0 vs sr45         38         42      29
#>
#>
"orfAnalysis" %in% names(a)
#>TRUE

##3: extract nucleotide and amino acid seq
nu <- extractSequence(dexseq, 
  pathToOutput = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out',
  writeToFile=FALSE)

nu
##3.1 Add CPC2 analysis
switch <- analyzeCPC2(
  switchAnalyzeRlist   = nu,
  pathToCPC2resultFile = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/result_cpc2.txt',
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)
#>Added coding potential to 133 (100%) transcripts
#>
##3.2 Add PFAM analysis
switch <- analyzePFAM(
  switchAnalyzeRlist   = switch,
  pathToPFAMresultFile = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/result_pfam.txt',
  showProgress= FALSE
)
?analyzePFAM
##3.3 Add SignalP analysis
switch <- analyzeSignalP(
  switchAnalyzeRlist       = switch,
  pathToSignalPresultFile  = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/result_signal_IP.txt'
)
#>Added signal peptide information to 1 (0.75%) transcripts 
#>
##3.4 
switch <- analyzeIUPred2A(
  switchAnalyzeRlist        = switch,
  pathToIUPred2AresultFile = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/IUPred2A.txt',
  showProgress = FALSE
)
#>Step 1 of 4 : Reading results into R...
#>Step 2 of 4 : Extracting regions of interest...
#>Step 3 of 4 : Integrating IDR with binding site predictions...
#>Step 4 of 4 : Converting AA coordinats to transcript and genomic coordinats...
#>Added IDR information to 15 (11.28%) transcripts
#>
##3.5 Add DeepLoc2 analysis
switch <- analyzeDeepLoc2(
  switchAnalyzeRlist = switch,
  pathToDeepLoc2resultFile = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/deeploc_results.csv',
  quiet = FALSE
)
#>Added subcellular information to 63 (47.37%) transcripts

##3.6 Add DeepTMHMM analysis
switch <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = switch,
  pathToDeepTMHMMresultFile = '/Users/baiprongfah/Master-ICL/MPMS-project2/isoformswitchanalyzer/switch_list_out/deeptmr_results.gff3',
  showProgress=FALSE)
#>Step 1 of 2: Reading results into R...
#>Step 2 of 2: Converting AA coordinats to transcript and genomic coordinats...
#>Added topology information to 63 transcripts (47.37%).

switch
##4: predicting alternative splicing 
splice <- analyzeAlternativeSplicing(switchAnalyzeRlist = switch,
                                     quiet=TRUE)
table(splice$AlternativeSplicingAnalysis$IR )
#4.1: Global splicing analysis
extractSplicingSummary(splice,
                       splicingToAnalyze = c('IR','A3','A5','ES'))
extractSplicingEnrichment(splice)
extractSplicingGenomeWide(splice,
                          featureToExtract = 'isoformUsage',
                          alpha=0.05,
                          dIFcutoff = 0.1,
                          log2FCcutoff = 1,)

###5:analysis of consequences
# the consequences highlighted in the text above
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

consequence <- analyzeSwitchConsequences(
  splice,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.1, 
  showProgress=FALSE
)

extractConsequenceSummary(consequence,
                          consequencesToAnalyze='all',
                          includeCombined=FALSE,
                          asFractionTotal=FALSE,
                          alpha=0.05,
                          dIFcutoff=0.1,
                          plot=FALSE,
                          plotGenes=FALSE,
                          simplifyLocation = TRUE,
                          removeEmptyConsequences = FALSE,
                          localTheme=theme_bw(),
                          returnResult = TRUE)
extractConsequenceEnrichment(consequence)
extractConsequenceEnrichmentComparison(consequence)
extractConsequenceGenomeWide(consequence)

### Extract top switching genes (by q-value)
switchingIso <- extractTopSwitches( 
  consequence, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = FALSE
)
write.csv(switchingIso, 'switching_iso.csv')
subset(switchingIso, gene_id == 'AT3G25540')
switchPlot(consequence, gene = 'AT4G11120')

switchingIso$gene_id
