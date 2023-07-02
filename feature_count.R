# running R commands to get feature counts of the bam files

# featureCounts is part of the Rsubread suite of tools
library(Rsubread)

# now I want to loop through the bam files and use each bam in the feature counts tool
# first i want to create a txt file of all bam names and then use readlines to look through
# the names on each line in the bam txt file


# featureCounts also takes an annotation file, so I will use the same gtf file from 
# ther cut&run pipeline

annotation_gtf_file = '../hg38.knownGene.gtf.gz'

bam_files = readLines(paste("./store_bam_files", "/bam_list.txt", sep = ""))

for (x in bam_files) {

file_name = paste0("./store_bam_files/", x)
output_name = paste0(x, "_counts.txt")

fc = featureCounts( files = file_name,
annot.ext = annotation_gtf_file,
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id",
isPairedEnd = TRUE,
countReadPairs = TRUE
) 


write.table(fc$counts, file = output_name, sep = "\t", quote = FALSE)

}


library(tidyverse)

# a list of all the count files in this directory
files = list.files(path = ".", pattern = "*counts.txt")

total_files = length(files)

list_names = list()

for (x in c(1:total_files)) { file_name = files[x]
var_name = paste0("x", gsub("-", "_", file_name)) 
assign(var_name, read.csv(files[x], sep = "\t", header = FALSE, skip = 1))



}


# the names will be stored as strings. we dont want that so use get

df_list = lapply(list_names, get)


first_column = df_list[[1]][,1]

second_columns = lapply(df_list, function(df) df[, 2])

combined_df = do.call(cbind, c(list(first_column), second_columns))

rownames(combined_df) = combined_df[,1]

combined_df = combined_df[,-1]

colnames(combined_df) = list_names


# now we make the matrix 
cts = as.matrix(combined_df)
mode(cts) = "numeric"

coldata = data.frame( condition = c("control", 
"control","control", "control", "control", 
"control","knockdown","knockdown","knockdown","knockdown","knockdown", "knockdown"), row.names = 
list_names)

coldata$condition = as.factor(coldata$condition)


library("DESeq2")

dds = DESeqDataSetFromMatrix( countData = cts, colData = coldata, design = ~ condition)

dds = DESeq(dds)

res = results(dds)


# below I use lfcShrink for log fold change shrinkage for visulization and ranking
# helps with ranking of genes
# remember to site that we used apeglm as the shrinkage type. find in deseq2 package manuel
resLFC = lfcShrink(dds, coef="condition_knockdown_vs_control", type= "apeglm")


# now I want to order the results by pvalue. most significant to least significant

resOrdered = res[order(res$pvalue),]

# do the same for the lfc also
resLFCordered = resLFC[order(resLFC$pvalue),]

# if i want the adjusted pvalue to be set to 0.05, then I need to change alpha from defualt at 0.1

res05 = results(dds, alpha = 0.05)


# plotting plotma png 
png(file = "ma_plot_lfc.png", height = 500, width = 500)
ma_plot_lfc = plotMA(resLFC, ylim = c(-2,2))
dev.off()


# now we take a subset of the genes that pass this value, then save as a table

resLFC_0.05 = subset(resLFCordered, padj <= 0.05)

write.table(as.data.frame(resLFC_0.05), file = "condition_kd_vs_ctr.csv")


# hopefully installing this and choosing all is no problem
BiocManager::install('EnhancedVolcano')
a
library(EnhancedVolcano)

v_plot = EnhancedVolcano(resLFC, lab = rownames(resLFC), x = 'log2FoldChange', y = 'padj', pCutoff 
=10e-5, FCcutoff = 0.5)
png(file = 'v_plot.png', height = 500, width = 500)
v_plot
dev.off()
