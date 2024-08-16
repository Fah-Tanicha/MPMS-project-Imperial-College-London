rm(list=ls()) # clear all variables and graphics
graphics.off()

library(ggplot2)

resLFC <- read.csv('/Users/baiprongfah/Master-ICL/MPMS-project2/alba_original_count/deseq2/resLFC.csv')
resFLC
#label the genes
resLFC$diffexpressed <- 'NO'
resLFC$diffexpressed[resLFC$logFC > 1 & resLFC$P.Value < 0.05] <- 'UP'
resLFC$diffexpressed[resLFC$logFC < -1 & resLFC$P.Value < 0.05] <- 'DOWN'

resLFC$delabel <- NA

#see how many genes are up or down
table(unlist(strsplit(tolower(resLFC$diffexpressed), " ")))

ggplot(data = resLFC, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.9, 0.9), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential expression', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Differential gene expression of wildtype vs sr45-1 mutant by DESeq2') # Plot title
