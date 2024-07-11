library(ggplot2)
library(limma)
library(maSigPro)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(reshape2)
library(DESeq2)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(singscore)
library(preprocessCore)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(rtracklayer)
library(VennDiagram)
library(GenomicRanges)
library(genomation)
library(ChIPpeakAnno)
library(ChIPseeker)
library(bedr)
library(enrichR)
library(EnhancedVolcano)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT")
setEnrichrSite("Enrichr")
database <- listEnrichrDbs()
dbs=c("Reactome_2022","WikiPathway_2023_Human","KEGG_2021_Human","GO_Molecular_Function_2023",
      "MSigDB_Hallmark_2020")
# LNCaP ADT comprehensive analysis
count=read.table("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/RNAseq/star_salmon/salmon.merged.gene_counts.tsv",sep = '\t', header = T)
tpm=read.table("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/RNAseq/star_salmon/salmon.merged.gene_tpm.tsv",sep = '\t', header = T)
anno=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/LNCaP_RANseq_samples.csv")
rownames(count)=count$gene_id
rownames(tpm)=tpm$gene_id
count=count[,-c(1:2)]
tpm=tpm[,-c(1:2)]
rownames(anno)=anno$ID
anno$Batch=gsub("[\t].*","",anno$Batch)
# will not include H660 samples in this study
anno=subset(anno,anno$Treatment_refine != "H660")
anno=subset(anno,anno$Treatment_refine != "LNCaP-abl")
list=intersect(colnames(count),anno$ID)
count=count[,list]
tpm=tpm[,list]
anno=anno[list,]
anno$Time[which(grepl("DHT",anno$Treatment_refine))]="0D"
anno$Time[which(grepl("R1881",anno$Treatment_refine))]="0D"
anno$Time[which(grepl("ENZ",anno$Treatment_refine))]=anno$Treatment_refine[which(grepl("ENZ",anno$Treatment_refine))]
anno$Time=gsub("[_]","",anno$Time)
anno$Time=factor(anno$Time,levels = c("0D","0.5D","1D","2D","2.5D","2.7D","3D","3.5D","4D","ENZ2D","ENZ4D","ENZ5D"))
anno$class="CS-FBS"
anno$class[which(anno$Time=="0D")]="Control"
anno$class[which(grepl("ENZ",anno$Time))]="ENZ"
anno$class=factor(anno$class,levels = c("Control","CS-FBS","ENZ"))
# now convert TPM to log2TPM
tpm=log2(tpm+1)

pca=prcomp(t(tpm))
pca=as.data.frame(pca$x[,1:2])
identical(rownames(pca),rownames(anno))
pca$ID=anno$ID
pca$batch=anno$Batch
pca$group=anno$Treatment_refine

ggplot(pca,aes(x=PC1,y=PC2,color=batch))+geom_point()+theme_minimal()+theme(legend.position = "none")
# these.TPM data showed strong batch effect
# now remove batch effect on TPM
batch=anno$Batch
design_batch=model.matrix(~ batch)
tpm_cor=removeBatchEffect(tpm,batch = batch)

pca=prcomp(t(tpm_cor))
pca=as.data.frame(pca$x[,1:2])
identical(rownames(pca),rownames(anno))
pca$ID=anno$ID
pca$batch=anno$Batch
pca$group=anno$Treatment_refine
ggplot(pca,aes(x=PC1,y=PC2,color=batch))+geom_point()+geom_text_repel(aes(label = batch)) + theme_minimal()

# now remove outliar batches
anno=subset(anno,anno$Batch!="GSE122923")
anno=subset(anno,anno$Batch!="GSE106560")
anno=subset(anno,anno$Batch!="GSE117305")
anno=subset(anno,anno$Batch!="GSE130534")
anno=subset(anno,anno$Batch!="GSE56512")
anno=subset(anno,anno$Batch!="GSE99795")
anno=subset(anno,anno$Batch!="GSE137833")
anno=subset(anno,anno$Batch!="GSE79301")
anno=subset(anno,anno$Batch!="GSE152254")
anno=subset(anno,anno$Batch!="GSE110903")
tpm=tpm[,rownames(anno)]
count=count[,rownames(anno)]

# redo batch removal 
# now remove batch effect on TPM
batch=anno$Batch
design_batch=model.matrix(~ batch)
tpm_cor=removeBatchEffect(tpm,batch = batch)
pca=prcomp(t(tpm_cor))
pca=as.data.frame(pca$x[,1:2])
identical(rownames(pca),rownames(anno))
pca$ID=anno$ID
pca$batch=anno$Batch
pca$group=anno$Treatment_refine
pca$time=anno$Time
pca$class=anno$class
ggplot(pca,aes(x=PC1,y=PC2,color=group))+geom_point()+geom_text_repel(aes(label = group)) + 
  theme_classic()+theme(legend.position = "none")+coord_cartesian(xlim = c(-50, 50), ylim = c(-50, 50))+
  scale_color_hue(l=40)
ggsave("PCA_group.pdf",width = 20,height = 20,units = "cm")
ggplot(pca,aes(x=PC1,y=PC2,color=batch))+geom_point()+geom_text_repel(aes(label = batch)) + 
  theme_classic()+theme(legend.position = "none")+coord_cartesian(xlim = c(-50, 50), ylim = c(-50, 50))+
  scale_color_hue(l=40)
ggsave("PCA_batch.pdf",width = 20,height = 20,units = "cm")
ggplot(pca,aes(x=PC1,y=PC2,color=time))+geom_point()+geom_text_repel(aes(label = time)) + 
  theme_classic()+coord_cartesian(xlim = c(-50, 50), ylim = c(-50, 50))+
  scale_color_hue(l=40)
ggsave("PCA_time.pdf",width = 20,height = 20,units = "cm")
ggplot(pca,aes(x=PC1,y=PC2,color=class))+geom_point()+geom_text_repel(aes(label = class)) + 
  theme_classic()+coord_cartesian(xlim = c(-50, 50), ylim = c(-50, 50))+
  scale_color_hue(l=40)
ggsave("PCA_class.pdf",width = 20,height = 20,units = "cm")
# the result is good
# visulize genes first

tpm_cor = normalize.quantiles(tpm_cor,keep.names = T)
saveRDS(tpm_cor,"Normalized_log2_TPM.rds")
tpm_cor=round(tpm_cor,digits = 2)
tpm_cor=as.data.frame(t(tpm_cor))
tpm_cor$ID=rownames(tpm_cor)
identical(rownames(tpm_cor),rownames(anno))
tpm_cor$time=anno$Time
tpm_cor$class=anno$class
tpm_cor=melt(tpm_cor,id=c("ID","time","class"))
colnames(tpm_cor)=c("ID","time","class","Gene","Log2TPM")
saveRDS(tpm_cor,"log2TPM.rds")
saveRDS(anno,"LNCaP_ADT_anno.rds")

gene="FOXA1"
tem=subset(tpm_cor,tpm_cor$Gene==gene)
ggplot(tem, aes(x = class, group = class, y = Log2TPM, fill = class)) +
  geom_boxplot(outlier.colour = NA, lwd = 2) +
  geom_jitter(alpha = 0.5) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  scale_fill_npg() +
  stat_compare_means(comparisons = list(c("CS-FBS", "Control"), c("ENZ", "Control")),
                     method = "t.test", label = "p.signif",
                     size = 8, # Adjust line width
                     tip.length = 0.1, # Adjust length of the tip of the bracket
                     vjust = 0.2, # Adjust vertical justification of text
                     label.x.npc = "right", # Position of label
                     label.y.npc = c(0.9, 0.8),
                     bracket.size = 1) # Position of label



# now perform DE analysis
count=round(count,digits = 0)
anno_cs=subset(anno,anno$class != "ENZ")
colData=anno_cs[,c(5,11)]
count_cs=count[,rownames(anno_cs)]
tpm_cor_matrix=readRDS("Normalized_log2_TPM.rds")
tpm_table=as.data.frame(tpm_cor_matrix)
tpm_cs <- tpm_table[rownames(anno_cs), ]

dds <- DESeqDataSetFromMatrix(countData = count_cs,
                              colData = colData, 
                              design = ~ class + Batch)
dds <- DESeq(dds)
resShrink <- lfcShrink(dds, coef="class_CS.FBS_vs_Control", type="apeglm")
de_result=as.data.frame(resShrink)
EnhancedVolcano(de_result,lab = rownames(de_result),x="log2FoldChange",y="padj",xlim = c(-8,8),
                boxedLabels = T,widthConnectors = 0.5,drawConnectors = T,colAlpha = 1,pointSize = 0.8,
                selectLab = c("KLK3","KLK2","FKBP5","AMIGO2","PLD1","ROR2"))+theme_classic()
ggsave("volcanoplot.pdf",width = 18,height = 22,units = "cm")

de_result=subset(de_result,(padj<0.05)&(baseMean>20))
saveRDS(de_result,"DE_table.rds")



# do enrichments
# this time I will use the batch correct TPM 
net <- get_progeny(organism = 'human', top = 500)
sample_acts <- run_mlm(mat=tpm_cor_matrix, net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

sample_acts_mat=as.data.frame(sample_acts_mat)
sample_acts_mat$ID=rownames(sample_acts_mat)
sample_acts_mat=merge(sample_acts_mat,anno[,c(1,11)],by="ID")
sample_acts_mat=melt(sample_acts_mat,id=c("ID","class"))
colnames(sample_acts_mat)=c("ID","group","pathway","score")

pathways <- unique(sample_acts_mat$pathway)

# Create a directory to save the plots
dir.create("signaling_activity", showWarnings = FALSE)

# Loop over each pathway and create a plot
for (pthway in pathways) {
  # Subset data for the current pathway
  subset_data <- subset(sample_acts_mat, pathway == pthway)
  
  # Create ggplot
  p <- ggplot(subset_data, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.colour = NA, lwd = 2) +
    geom_jitter(alpha = 0.5) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none",
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black")) +
    scale_fill_npg() +
    stat_compare_means(comparisons = list(c("Control", "CS-FBS"), c("Control", "ENZ")),
                       method = "t.test", label = "p.signif",
                       size = 8, tip.length = 0.1,
                       vjust = 0.2, label.x.npc = "right",
                       label.y.npc = c(0.9, 0.8),
                       bracket.size = 1) +
    labs(title = paste("Pathway:", pthway))
  # Save the plot as a PDF
  ggsave(paste0("signaling_activity/plot_", pthway, ".pdf"), p)
}


# now TF activity
net <- get_collectri(organism='human', split_complexes=FALSE)
sample_acts <- run_ulm(mat=tpm_cor_matrix, net=net, .source='source', .target='target',
                       .mor='mor', minsize = 5)

sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

sample_acts_mat=as.data.frame(sample_acts_mat)
sample_acts_mat$ID=rownames(sample_acts_mat)
sample_acts_mat=merge(sample_acts_mat,anno[,c(1,11)],by="ID")
sample_acts_mat=melt(sample_acts_mat,id=c("ID","class"))
colnames(sample_acts_mat)=c("ID","group","pathway","score")

pathways <- unique(sample_acts_mat$pathway)

# Create a directory to save the plots
dir.create("TF_activity", showWarnings = FALSE)

# Loop over each pathway and create a plot
for (pthway in pathways) {
  # Subset data for the current pathway
  subset_data <- subset(sample_acts_mat, pathway == pthway)
  
  # Create ggplot
  p <- ggplot(subset_data, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.colour = NA, lwd = 2) +
    geom_jitter(alpha = 0.5) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none",
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black")) +
    scale_fill_npg() +
    stat_compare_means(comparisons = list(c("Control", "CS-FBS"), c("Control", "ENZ")),
                       method = "t.test", label = "p.signif",
                       size = 8, tip.length = 0.1,
                       vjust = 0.2, label.x.npc = "right",
                       label.y.npc = c(0.9, 0.8),
                       bracket.size = 1) +
    labs(title = paste("Pathway:", pthway))
  # Save the plot as a PDF
  ggsave(paste0("TF_activity/plot_", pthway, ".pdf"), p)
}

# Functional enrichment using EnrichR
de_result=readRDS("DE_table.rds")
up=subset(de_result,de_result$log2FoldChange>0.5)
up=rownames(up)
websiteLive <- TRUE
if (websiteLive) {
  enriched <- enrichr(up, dbs)
}
tem=enriched[[1]]
for (i in 2:length(enriched)) {
  tem=rbind(tem,enriched[[i]])
}
tem=subset(tem,(Adjusted.P.value<0.05)&(tem$Combined.Score>5))
enrich_up=tem

down=subset(de_result,de_result$log2FoldChange<(-0.5))
down=rownames(down)
websiteLive <- TRUE
if (websiteLive) {
  enriched <- enrichr(down, dbs)
}
tem=enriched[[1]]
for (i in 2:length(enriched)) {
  tem=rbind(tem,enriched[[i]])
}
tem=subset(tem,(Adjusted.P.value<0.05)&(tem$Combined.Score>5))
enrich_down=tem


write.csv(enrich_up,"Up_signaling.csv")
write.csv(enrich_down,"Dowm_signaling.csv")

# the result is okay but can be better if only keep the protein coding genes
gtf=readGFF("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/analysis_for_other_labs/Mu/IGV/bigwig/gencode.gtf")
gtf=subset(gtf,gtf$type=="gene")
gtf=subset(gtf,gtf$gene_type=="protein_coding")
protein_coding=unique(na.omit(gtf$gene_name))

de_result=de_result[rownames(de_result) %in% protein_coding,]


up=subset(de_result,de_result$log2FoldChange>0.5)
up=rownames(up)
websiteLive <- TRUE
if (websiteLive) {
  enriched <- enrichr(up, dbs)
}
tem=enriched[[1]]
for (i in 2:length(enriched)) {
  tem=rbind(tem,enriched[[i]])
}
tem=subset(tem,(Adjusted.P.value<0.05)&(tem$Combined.Score>5))
enrich_up=tem

down=subset(de_result,de_result$log2FoldChange<(-0.5))
down=rownames(down)
websiteLive <- TRUE
if (websiteLive) {
  enriched <- enrichr(down, dbs)
}
tem=enriched[[1]]
for (i in 2:length(enriched)) {
  tem=rbind(tem,enriched[[i]])
}
tem=subset(tem,(Adjusted.P.value<0.05)&(tem$Combined.Score>5))
enrich_down=tem
write.csv(enrich_up,"Up_signaling.csv")
write.csv(enrich_down,"Dowm_signaling.csv")



### now run all selected signatures
# first assemble all related RNAseq data
tpm_cor_matrix=readRDS("Normalized_log2_TPM.rds")
tpm_table=as.data.frame(tpm_cor_matrix)
anno=readRDS("LNCaP_ADT_anno.rds")
anno2=anno[order(anno$class,anno$Time),]
tpm_table=tpm_table[,rownames(anno2)]
# only keep the protein coding gene
tpm_table=tpm_table[rownames(tpm_table) %in% protein_coding,]


# now load signatures
sig=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/selected_signatures.csv")
all_genes=unique(na.omit(unlist(sig)))
sig_list = setNames(
  lapply(sig, function(column) {
    # Remove empty strings
    clean_column <- column[column != ""]
    as.list(clean_column)
  }), 
  names(sig)
)
sig_list2=lapply(sig_list,unique)
gene_sets <- lapply(names(sig_list2), function(name) {
  GeneSet(setName = name, geneIds = unlist(sig_list2[[name]]))
})
gene_set_collection <- GeneSetCollection(gene_sets)

rankData=rankGenes(tpm_table)
score=multiScore(rankData = rankData,upSetColc = gene_set_collection)
score=as.data.frame(score$Scores)

### calculate the median value of enrichment score for better visulization
score_long <- score %>%
  rownames_to_column("Pathway") %>%
  gather(key = "Sample", value = "Score", -Pathway)

# Merge score_long with anno2 based on sample names
merged_data <- score_long %>%
  left_join(anno2, by = c("Sample" = "ID"))

# Subset the data to include only the "CS-FBS" and "Control" groups
subset_data <- merged_data %>%
  filter(class %in% c("CS-FBS", "Control"))

# Function to perform double-tailed t-test for each pathway
perform_t_test <- function(data) {
  t_test_result <- t.test(Score ~ class, data = data)
  return(t_test_result$p.value)
}

# Apply the t-test function to each pathway
t_test_results <- subset_data %>%
  group_by(Pathway) %>%
  summarise(p_value = perform_t_test(pick(everything())))
t_test_results=as.data.frame(t_test_results)
rownames(t_test_results)=t_test_results$Pathway
t_test_results=t_test_results[rownames(score),]






pdf("ADT_enrichment_heatmap.pdf",width = 30,height = 20)
pheatmap(score,cluster_cols = F,scale = "row",gaps_col = c(106, 134),gaps_row = c(1,23),
         color = colorRampPalette(c("blue", "white","red"))(length(seq(-2,2,0.1))-1),
         breaks = seq(-2,2,0.1),border_color = NA,
         cellheight = 15,cellwidth = 10,cluster_rows = F,show_colnames = F)
dev.off()


# Scale the data by row (each pathway) before calculating fold-change
scaled_score <- t(apply(score, 1, scale))
scaled_score <- as.data.frame(scaled_score)
rownames(scaled_score) <- rownames(score)
colnames(scaled_score) <- colnames(score)

# Convert scaled scores to long format
scaled_score_long <- scaled_score %>%
  rownames_to_column("Pathway") %>%
  gather(key = "Sample", value = "Scaled_Score", -Pathway)

# Merge scaled_score_long with anno2 based on sample names
scaled_merged_data <- scaled_score_long %>%
  left_join(anno2, by = c("Sample" = "ID"))

# Subset the scaled data to include only the "CS-FBS" and "Control" groups
scaled_subset_data <- scaled_merged_data %>%
  filter(class %in% c("CS-FBS", "Control"))

# Calculate the mean of scaled scores in CS-FBS vs Control for each pathway
mean_score_results <- scaled_subset_data %>%
  group_by(Pathway, class) %>%
  summarise(mean_score = mean(Scaled_Score))

# Perform power 2 transformation
mean_score_results <- mean_score_results %>%
  mutate(mean_score_power2 = 2^mean_score)
mean_score_results$class=factor(mean_score_results$class,levels = c("Control","CS-FBS"))
mean_score_results=mean_score_results[,-3]
# Calculate the fold-change of transformed mean scores in CS-FBS vs Control for each pathway
fold_change_results <- mean_score_results %>%
  spread(key = class, value = mean_score_power2) %>%
  mutate(fold_change = `CS-FBS` / Control)

# Perform log2 transformation on the fold-change
fold_change_results <- fold_change_results %>%
  mutate(log2_fold_change = log2(fold_change))

# Convert to data frame and reorder according to t_test_results
fold_change_results <- as.data.frame(fold_change_results)
rownames(fold_change_results) <- fold_change_results$Pathway
fold_change_results <- fold_change_results[rownames(t_test_results),]

# Reorder the fold_change_results according to t_test_results
fold_change_results$Pathway <- factor(fold_change_results$Pathway, levels = unique(rownames(t_test_results)))

# Plot the log2 fold-change as a bar plot
ggplot(fold_change_results, aes(x = Pathway, y = log2_fold_change)) +
  geom_bar(stat = "identity", fill = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(title = "Log2 Fold-Change of Scaled Scores in CS-FBS vs Control", y = "Log2 Fold-Change", x = "Pathway")
ggsave("ssGSEA_enrichment_barplot.pdf",width = 18,height = 30)





# now visulize the DE genes by heatmap

hm=as.data.frame(tpm_cor_matrix[c(up,down),])

anno2=anno[order(anno$class,anno$Time),]
hm=hm[,rownames(anno2)]
anno3=anno2[,c(10,11),drop=F]


breaks <- seq(-2, 2, length.out = 100)
colors <- colorRampPalette(c("navy", "black", "firebrick3"))(length(breaks)-1)
pdf("DE_expression_heatmap.pdf",width = 18,height = 18)
pheatmap(hm, 
         scale = "row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = anno3, 
         gaps_row = 1004, 
         gaps_col = c(106, 134), 
         breaks = breaks, 
         color = colors,
         show_rownames = F,
         show_colnames = F)
dev.off()
# now try to generate heatmap on TPM without batch removal
# the result is not good

#hm=as.data.frame(tpm[c(up,down),])
#anno2=anno[order(anno$class,anno$Time),]
#hm=hm[,rownames(anno2)]
#anno3=anno2[,c(10,11),drop=F]
#breaks <- seq(-2, 2, length.out = 100)
#colors <- colorRampPalette(c("navy", "black", "firebrick3"))(length(breaks)-1)
#pheatmap(hm, 
#         scale = "row", 
#         cluster_rows = FALSE, 
#         cluster_cols = FALSE, 
#         annotation_col = anno3, 
#         gaps_row = 1004, 
#         gaps_col = c(106, 134), 
#         breaks = breaks, 
#         color = colors,
#         show_rownames = F,
#         show_colnames = F)




# now generate bed files to characterize other omics
de=readRDS("DE_table.rds")
tpm_cor_matrix=readRDS("Normalized_log2_TPM.rds")
up=rownames(subset(de,de$log2FoldChange>0.5))
down=rownames(subset(de,de$log2FoldChange<(-0.5)))
ave=apply(tpm_cor_matrix,1,mean)
high=rownames(as.data.frame(tpm_cor_matrix))[which(ave>2)]
high=setdiff(high,unique(c(up,down)))
low=rownames(as.data.frame(tpm_cor_matrix))[which(ave<=2)]
low=setdiff(low,unique(c(up,down)))
# now try another GRCh38 GTF which I used for ChIPseq alignment
gtf=readGFF("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/analysis_for_other_labs/Mu/IGV/bigwig/gencode.gtf")
gtf=subset(gtf,gtf$type=="gene")
bed=gtf[,c(1,4,5,9,11)]
write.table(bed, file = "all.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% up, ]
write.table(subset_bed_data, file = "up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% down, ]
write.table(subset_bed_data, file = "down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% high, ]
write.table(subset_bed_data, file = "high.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% low, ]
write.table(subset_bed_data, file = "low.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

### identifhy AR binding regions
edb=TxDb.Hsapiens.UCSC.hg38.knownGene
ar1=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E17_AR_peaks.narrowPeak")
ar1 <- annotatePeak(ar1, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar1=as.data.frame(ar1)


ar2=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E18_AR_peaks.narrowPeak")
ar2 <- annotatePeak(ar2, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar2=as.data.frame(ar2)

ar3=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E20_DHT_AR_peaks.narrowPeak")
ar3 <- annotatePeak(ar3, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar3=as.data.frame(ar3)

ar4=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E22_DHT_AR_peaks.narrowPeak")
ar4 <- annotatePeak(ar4, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar4=as.data.frame(ar4)

ar5=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E25_AR_peaks.narrowPeak")
ar5 <- annotatePeak(ar5, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar5=as.data.frame(ar5)

ar6=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E28_AR_peaks.narrowPeak")
ar6 <- annotatePeak(ar6, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")
ar6=as.data.frame(ar6)


ar_all=rbind(ar1,ar2,ar3,ar4,ar5,ar6)
ar_all=subset(ar_all,ar_all$score>100)

ar_bed=ar_all[,c(1,2,3,6)]
write.table(ar_bed, file = "ARBS.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



### 
edb=TxDb.Hsapiens.UCSC.hg38.knownGene
ar1=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E17_AR_peaks.narrowPeak")
ar2=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E18_AR_peaks.narrowPeak")
ar3=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E20_DHT_AR_peaks.narrowPeak")
ar4=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E22_DHT_AR_peaks.narrowPeak")
ar5=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E25_AR_peaks.narrowPeak")
ar6=readNarrowPeak("/Volumes/SC1/LNCaP_all_ChIP_TFs/star/mergedLibrary/macs2/narrowPeak/E28_AR_peaks.narrowPeak")

tem=findOverlapsOfPeaks(ar1,ar2,ar3,ar4,ar5)
tem2=tem$mergedPeaks

tem3=annotatePeak(tem2, tssRegion=c(-6000, 6000),
                  TxDb=edb, annoDb="org.Hs.eg.db")
tem4=as.data.frame(tem3)
tem4$merged_name=sapply(tem4$peakNames, function(x) paste(x, collapse = "_"))
ar_bed_overlap=tem4[,c(1,2,3,19)]
write.table(ar_bed_overlap, file = "ARBS_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



ar1 <- annotatePeak(ar1, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")



ar2 <- annotatePeak(ar2, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")


ar3 <- annotatePeak(ar3, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")


ar4 <- annotatePeak(ar4, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")


ar5 <- annotatePeak(ar5, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")


ar6 <- annotatePeak(ar6, tssRegion=c(-6000, 6000),
                    TxDb=edb, annoDb="org.Hs.eg.db")


### identify LNCaP and H660 specific ATAC regions
lncap_2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_3=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_4=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_5=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_6=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_3_REP1.mLb.clN_peaks.narrowPeak")
lncap_7=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_4_REP1.mLb.clN_peaks.narrowPeak")


# find overlapping lncap peaks

lncap=findOverlapsOfPeaks(lncap_2,lncap_4,lncap_7)
lncap_overlap=lncap$mergedPeaks

H660=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/H660_E3_1_REP1.mLb.clN_peaks.narrowPeak")

overlap=findOverlapsOfPeaks(lncap_overlap,H660)
edb=TxDb.Hsapiens.UCSC.hg38.knownGene

lncap_unique=overlap$peaklist$lncap_overlap
lncap_unique=annotatePeak(lncap_unique, tssRegion=c(-6000, 6000),
                  TxDb=edb, annoDb="org.Hs.eg.db")
lncap_unique=as.data.frame(lncap_unique)
lncap_unique$merged_name=sapply(lncap_unique$peakNames, function(x) paste(x, collapse = "_"))
lncap_unique=lncap_unique[,c(1,2,3,19)]
write.table(lncap_unique, file = "lncap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



H660_unique=overlap$peaklist$H660
H660_unique=annotatePeak(H660_unique, tssRegion=c(-6000, 6000),
                          TxDb=edb, annoDb="org.Hs.eg.db")
H660_unique=as.data.frame(H660_unique)
H660_unique$merged_name=sapply(H660_unique$peakNames, function(x) paste(x, collapse = "_"))
H660_unique=H660_unique[,c(1,2,3,19)]
write.table(H660_unique, file = "h660.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


common=overlap$peaklist$`lncap_overlap///H660`
common=annotatePeak(common, tssRegion=c(-6000, 6000),
                          TxDb=edb, annoDb="org.Hs.eg.db")
common=as.data.frame(common)
common$merged_name=sapply(common$peakNames, function(x) paste(x, collapse = "_"))
common=common[,c(1,2,3,19)]
write.table(common, file = "common.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# identify the H660 and DU145 specific regions
H660=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/H660_E3_1_REP1.mLb.clN_peaks.annotatePeaks.txt")
H660=subset(H660,grepl("TSS",H660$Annotation))
H660=subset(H660,H660$Peak.Score>500)
H660=subset(H660,abs(H660$Distance.to.TSS)<3000)
du145=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/DU145_E2_1_REP1.mLb.clN_peaks.annotatePeaks.txt")
du145=subset(du145,grepl("TSS",du145$Annotation))
du145=subset(du145,du145$Peak.Score>500)
du145=subset(du145,abs(du145$Distance.to.TSS)<3000)
lncap=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_4_REP1.mLb.clN_peaks.annotatePeaks.txt")
#lncap=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_2_REP1.mLb.clN_peaks.annotatePeaks.txt")
lncap=subset(lncap,grepl("TSS",lncap$Annotation))
lncap=subset(lncap,lncap$Peak.Score>500)
lncap=subset(lncap,abs(lncap$Distance.to.TSS)<3000)

all=intersect(intersect(lncap$Gene.Name,du145$Gene.Name),H660$Gene.Name)
unique_h660=setdiff(H660$Gene.Name,unique(c(lncap$Gene.Name,du145$Gene.Name)))
unique=setdiff(c(H660$Gene.Name,du145$Gene.Name),lncap$Gene.Name)
unique_lncap=setdiff(lncap$Gene.Name,c(H660$Gene.Name,du145$Gene.Name))
negative=setdiff(rownames(tpm_cor_matrix),unique(c(lncap$Gene.Name,du145$Gene.Name,H660$Gene.Name)))

venn.diagram(list(up=up,down=down,all=all,unique=unique,lncap=unique_lncap),"venn.tiff")


subset_bed_data=bed[bed$gene_name %in% all, ]
write.table(subset_bed_data, file = "positive.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% unique, ]
write.table(subset_bed_data, file = "lp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% unique_lncap, ]
write.table(subset_bed_data, file = "lncap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$gene_name %in% negative, ]
write.table(subset_bed_data, file = "negative.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)














#lncap_1=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E13_LNCaP_DHT5h_REP1.mLb.clN_peaks.narrowPeak")
lncap_2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_3=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_4=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_5=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_6=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_3_REP1.mLb.clN_peaks.narrowPeak")
lncap_7=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_4_REP1.mLb.clN_peaks.narrowPeak")

edb=TxDb.Hsapiens.UCSC.hg38.knownGene

lncap_2 <- annotatePeak(lncap_2, tssRegion=c(-3000, 3000),
                             TxDb=edb, annoDb="org.Hs.eg.db")
lncap_2=as.data.frame(lncap_2)
lncap_2=subset(lncap_2,(lncap_2$score>150)&grepl("Promoter",lncap_2$annotation))
lncap_3 <- annotatePeak(lncap_3, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_3=as.data.frame(lncap_3)
lncap_3=subset(lncap_3,(lncap_3$score>150)&grepl("Promoter",lncap_3$annotation))
lncap_4 <- annotatePeak(lncap_4, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_4=as.data.frame(lncap_4)
lncap_4=subset(lncap_4,(lncap_4$score>150)&grepl("Promoter",lncap_4$annotation))
lncap_5 <- annotatePeak(lncap_5, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_5=as.data.frame(lncap_5)
lncap_5=subset(lncap_5,(lncap_5$score>150)&grepl("Promoter",lncap_5$annotation))
lncap_6 <- annotatePeak(lncap_6, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_6=as.data.frame(lncap_6)
lncap_6=subset(lncap_6,(lncap_6$score>150)&grepl("Promoter",lncap_6$annotation))
lncap_7 <- annotatePeak(lncap_7, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_7=as.data.frame(lncap_7)
lncap_7=subset(lncap_7,(lncap_7$score>150)&grepl("Promoter",lncap_7$annotation))


lncap_E16=unique(intersect(lncap_2$transcriptId,lncap_3$transcriptId))
lncap_E22=unique(intersect(intersect(intersect(lncap_4$transcriptId,lncap_5$transcriptId),lncap_6$transcriptId),lncap_7$transcriptId))
lncap=intersect(lncap_E16,lncap_E22)

lncap_adt1=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_CX_16h_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_adt2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_ENZ168h_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_adt3=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_ENZ168h_2_REP1.mLb.clN_peaks.narrowPeak")

lncap_adt1 <- annotatePeak(lncap_adt1, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
lncap_adt1=as.data.frame(lncap_adt1)
lncap_adt1=subset(lncap_adt1,(lncap_adt1$score>150)&grepl("Promoter",lncap_adt1$annotation))

lncap_adt2 <- annotatePeak(lncap_adt2, tssRegion=c(-3000, 3000),
                           TxDb=edb, annoDb="org.Hs.eg.db")
lncap_adt2=as.data.frame(lncap_adt2)
lncap_adt2=subset(lncap_adt2,(lncap_adt2$score>150)&grepl("Promoter",lncap_adt2$annotation))

lncap_adt3 <- annotatePeak(lncap_adt3, tssRegion=c(-3000, 3000),
                           TxDb=edb, annoDb="org.Hs.eg.db")
lncap_adt3=as.data.frame(lncap_adt3)
lncap_adt3=subset(lncap_adt3,(lncap_adt3$score>150)&grepl("Promoter",lncap_adt3$annotation))


lncap_adt=intersect(intersect(lncap_adt1$transcriptId,lncap_adt2$transcriptId),lncap_adt3$transcriptId)

du145_1=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/DU145_E1_1_REP1.mLb.clN_peaks.narrowPeak")
du145_2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/DU145_E2_1_REP1.mLb.clN_peaks.narrowPeak")
du145_1 <- annotatePeak(du145_1, tssRegion=c(-3000, 3000),
                           TxDb=edb, annoDb="org.Hs.eg.db")
du145_1=as.data.frame(du145_1)
du145_1=subset(du145_1,(du145_1$score>150)&grepl("Promoter",du145_1$annotation))

du145_2 <- annotatePeak(du145_2, tssRegion=c(-3000, 3000),
                           TxDb=edb, annoDb="org.Hs.eg.db")
du145_2=as.data.frame(du145_2)
du145_2=subset(du145_2,(du145_2$score>150)&grepl("Promoter",du145_2$annotation))


du145=intersect(du145_1$transcriptId,du145_2$transcriptId)

h660=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/H660_E3_1_REP1.mLb.clN_peaks.narrowPeak")
h660 <- annotatePeak(h660, tssRegion=c(-3000, 3000),
                        TxDb=edb, annoDb="org.Hs.eg.db")
h660=as.data.frame(h660)
h660=subset(h660,(h660$score>150)&grepl("Promoter",h660$annotation))
h660=unique(h660$transcriptId)


high=unique(intersect(lncap,c(du145,h660)))
lncap_uni=setdiff(lncap,c(du145,h660))
lp_unique=setdiff(unique(c(du145,h660)),lncap)
low=setdiff(bed$transcript_id,c(high,lncap_uni,lp_unique))

gtf=readGFF("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/hg38.knownGene.gtf")
gtf=subset(gtf,gtf$type=="transcript")
bed=gtf[,c(1,4,5,9,10)]
write.table(bed, file = "all.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% high, ]
write.table(subset_bed_data, file = "positive.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% lp_unique, ]
write.table(subset_bed_data, file = "lp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% lncap_uni, ]
write.table(subset_bed_data, file = "lncap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% low, ]
write.table(subset_bed_data, file = "negative.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)






h660=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/H660_E3_1_REP1.mLb.clN_peaks.narrowPeak")
h660 <- annotatePeak(h660, tssRegion=c(-3000, 3000),
                     TxDb=edb, annoDb="org.Hs.eg.db")
h660=as.data.frame(h660)
h660=subset(h660,(h660$score>150)&grepl("Promoter",h660$annotation))


subset_h660 <- h660[h660$transcriptId %in% lp_unique, ]






































##### not used
gene=readBed("all.bed")
gene_bed_flank <- flank(gene, width = 3000,both = T)
lncap_2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_2=subset(lncap_2,lncap_2$score>150)
overlaps <- findOverlaps(lncap_2,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
lncap_promoter=gene_subset
lncap_promoter_name2=lncap_promoter$score

gene_bed_flank <- flank(gene, width = 3000,both = T)
lncap_3=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_ATAC_result2/star/merged_library/macs2/narrow_peak/E16_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_3=subset(lncap_3,lncap_3$score>150)
overlaps <- findOverlaps(lncap_3,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
lncap_promoter=gene_subset
lncap_promoter_name3=lncap_promoter$score

gene_bed_flank <- flank(gene, width = 3000,both = T)
lncap_4=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_1_REP1.mLb.clN_peaks.narrowPeak")
lncap_4=subset(lncap_4,lncap_4$score>150)
overlaps <- findOverlaps(lncap_4,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
lncap_promoter=gene_subset
lncap_promoter_name4=lncap_promoter$score


gene_bed_flank <- flank(gene, width = 3000,both = T)
lncap_5=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/LNCaP_atac_result3/star/merged_library/macs2/narrow_peak/E22_LNCaP_2_REP1.mLb.clN_peaks.narrowPeak")
lncap_5=subset(lncap_5,lncap_5$score>150)
overlaps <- findOverlaps(lncap_5,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
lncap_promoter=gene_subset
lncap_promoter_name5=lncap_promoter$score


lncap_promoter_name=unique(intersect(lncap_promoter_name2,lncap_promoter_name3))
#lncap_promoter_name=unique(c(intersect(lncap_promoter_name2,lncap_promoter_name3),intersect(lncap_promoter_name4,lncap_promoter_name5)))

du145_1=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/DU145_E1_1_REP1.mLb.clN_peaks.narrowPeak")
du145_1=subset(du145_1,du145_1$score>150)
overlaps <- findOverlaps(du145_1,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
du145_promoter=gene_subset
du145_promoter_name1=du145_promoter$score

du145_2=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/DU145_E2_1_REP1.mLb.clN_peaks.narrowPeak")
du145_2=subset(du145_2,du145_2$score>150)
overlaps <- findOverlaps(du145_2,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
du145_promoter=gene_subset
du145_promoter_name2=du145_promoter$score
du145_promoter_name=intersect(du145_promoter_name1,du145_promoter_name2)


h660=readNarrowPeak("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/LNCaP_ADT/ATACseq/H660_ATAC_result/star/merged_library/macs2/narrow_peak/H660_E3_1_REP1.mLb.clN_peaks.narrowPeak")
h660=subset(h660,h660$score>150)
overlaps <- findOverlaps(h660,gene_bed_flank)
gene_subset <- gene_bed_flank[queryHits(overlaps)]
index=duplicated(gene_subset$score)
gene_subset=gene_subset[!index,]
h660_promoter=gene_subset
h660_promoter_name=h660_promoter$score


all_name=unique(gene$score)

lncaP_unique=setdiff(lncap_promoter_name,unique(c(du145_promoter_name,h660_promoter_name)))
common=intersect(lncap_promoter_name,unique(c(du145_promoter_name,h660_promoter_name)))
lp=setdiff(unique(c(du145_promoter_name,h660_promoter_name)),lncap_promoter_name)
rest=setdiff(all_name,unique(c(lncap_promoter_name,du145_promoter_name,h660_promoter_name)))

subset_bed_data=bed[bed$transcript_id %in% common, ]
write.table(subset_bed_data, file = "common.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% lp, ]
write.table(subset_bed_data, file = "lp_unique.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% lncaP_unique, ]
write.table(subset_bed_data, file = "lncap_unique.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
subset_bed_data=bed[bed$transcript_id %in% rest, ]
write.table(subset_bed_data, file = "rest.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

