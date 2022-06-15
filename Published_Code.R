
## Update gencodev22 gene ids to the latest version of HGNC symbols
# Select a dataset from BioMart database
library(biomaRt)
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl')
# Create a vector of gene_ids of ensemble gene annotations from gencodev22
gene_ids <-gsub("\\..*","", gencodev22$gene_id)
# Create a biomaRt query
BM_new <-getBM(attributes = c("ensembl_gene_id","hgnc_symbol",'hgnc_id',"entrezgene_id"),
               filters = "ensembl_gene_id",
               values = gene_ids,
               mart = ensembl)
# Create a dataframe merging gencodev22 annotation and new annotations from biomart
gencode_BM <-merge(gencodev22, BM_new, all.x=T, by=1)
# Replace ids of genes without corresponding hgnc_symbols with gene_ids
index <- which(is.na(gencode_BM$hgnc_symbol))
ids <- gencode_BM$gene_name[index]
gencode_BM$hgnc_symbol <- replace(gencode_BM$hgnc_symbol, index, ids)

# List SCLC samples
EGAD_SCLC <- list.files(".")
a=1
p <- list()
for (i in EGAD_SCLC){
  p[[a]] <- read.table(i, header = F)
  names(p)[a]<-i
  a=a+1
}
# Merge all the elements in the list by the variable "V1" of Ensemble gene names
merged <- Reduce(function(x, y) merge(x, y, by = "V1"), p)
# Rename the variables
colnames(merged)<-c("gene", EGAD_SCLC)
# Remove chr strings and rename samples
merged_new <- merged[-(1:5),]
colnames(merged_new) <- gsub(pattern = "*.count.txt", replacement =  '', x=colnames(merged_new))
colnames(merged_new) <- gsub(pattern = "_1", replacement =  '', x=colnames(merged_new))
# Merge between ensemble annotations of count matrix and the hgnc symbols
merged_new$gene<-gsub("\\..*","", merged_new$gene)
SCLC_count_matrix <- merge(gencode_BM, merged_new, by=1)
# Remove duplicates
SCLC_count_matrix <- unique(SCLC_count_matrix)
# Remove gene_name e gene_id variables
SCLC_count_matrix <- SCLC_count_matrix[, 3:length(colnames(SCLC_count_matrix))]
# Aggregation of counts
SCLC_count_matrix <- aggregate(SCLC_count_matrix[,-1], by = list(SCLC_count_matrix$hgnc_symbol), FUN = sum)
# Gene names as rownames and remove Gene names variable
rownames(SCLC_count_matrix) <- SCLC_count_matrix$Gene
SCLC_count_matrix <- SCLC_count_matrix[,-1]

## Differential Expression Analysis

# Create coldata table
coldata_SCLC <- as.data.frame("samples" = colnames(SCLC_count_matrix), "type"= info_SCLC$type, "subtype"= "SCLC")
coldata_SCLC <- column_to_rownames(coldata_SCLC, "samples")
coldata_SCLC$type <- factor(coldata_SCLC$type)
coldata_SCLC$subtype <- factor(coldata_SCLC$subtype)

# Check rownames coldata and colnames count matrix
all(rownames(coldata_SCLC) == colnames(SCLC_count_matrix))

# DESEQ Dataset
library(DESeq2)
dds_SCLC <- DESeqDataSetFromMatrix(countData = SCLC_count_matrix,
                                              colData = coldata_SCLC,
                                              design = ~ type)

# Pre-filtering
keep <- rowSums(counts(dds_SCLC)) >= 5
dds_SCLC <- dds_SCLC[keep,]

# DESEQ
dds_SCLC <- DESeq(dds_SCLC)

#--------------------------------- FigS1
#PCA plot
vsd <- vst(dds_SCLC, blind=FALSE)
plotPCA(vsd, intgroup=c("type"))





#---------------------------------------
# Get counts
dds_SCLC_Exp <- counts(dds_SCLC)

# Filtering of tumour samples
dds_SCLC_Exp <- dds_SCLC[, colnames(dds_SCLC) %in% tumor_samples]

## Differential Expression Analysis for intergenic TEs

# coldata table for tumour and matched normal samples
#coldata_matched is a table with tumour and normal annotations for each sample
colnames(coldata_matched)[1] <- "Tumour"
colnames(coldata_matched)[2] <- "Normal"
coldata_matched <- gather(coldata_matched, "Type", "Samples", 1:2)
coldata_matched[,"Subtype"] <- "SCLC"
colnames(coldata_matched)[1] <- "Patient_ID"
coldata_matched <- column_to_rownames(coldata_matched, "Samples")
coldata_matched$Type <- factor(coldata_matched$Type)
coldata_matched$Patient_ID <- factor(coldata_matched$Patient_ID)

# Import intergenic RE raw counts
Intergenic <-readRDS(file= "RE_intergenic_1_raw_counts.RDS")

# Get counts
Intergenic <- Intergenic$counts

# Round counts
Intergenic <- round(Intergenic)

# Rename samples
colnames(Intergenic) <- gsub(pattern = "_1", replacement =  '', x=colnames(Intergenic))

# Check rownames coldata and colnames intergenic Re count matrix
all(rownames(coldata_matched) == colnames(REall_raw_counts))

# DeSeq DataSet
INTERGENIC_dds_TEs_matched <- DESeqDataSetFromMatrix(countData = Intergenic,
                                                   colData = coldata_matched,
                                                   design = ~ Patient_ID + Type)
# Pre-filtering
keep <- rowSums(counts(INTERGENIC_dds_TEs_matched)) >= 5
INTERGENIC_dds_TEs_matched <- INTERGENIC_dds_TEs_matched[keep,]

# DESEQ
INTERGENIC_dds_TEs_matched <- DESeq(INTERGENIC_dds_TEs_matched)
INTERGENIC_ddsresults_TEs_matched <- results(INTERGENIC_dds_TEs_matched,
                                             lfcThreshold = 0.58, altHypothesis="greaterAbs", alpha = 0.05 ) %>%
                                             as.data.frame()
# Shrinkage
INTERGENIC_TEs_shrink <- lfcShrink(INTERGENIC_dds_TEs_matched, coef = "Type_Tumour_vs_Normal",
                                   type = "apeglm", lfcThreshold = 0.58,
                                   res = INTERGENIC_ddsresults_TEs_matched) %>%
                                   as.data.frame()

# Create dataframe selecting lfcMLE, padj e lfcMAP
INTERGENIC_ddsresults_TEs_matched_table <- INTERGENIC_ddsresults_TEs_matched[,c("log2FoldChange","padj")]
INTERGENIC_ddsresults_TEs_matched_table <- rownames_to_column(INTERGENIC_ddsresults_TEs_matched_table, "TEs")
INTERGENIC_TEs_shrink_table <- rownames_to_column(INTERGENIC_TEs_shrink, "TEs")
INTERGENIC_results <- merge(INTERGENIC_TEs_shrink_table[,c("log2FoldChange", "TEs")],
                              INTERGENIC_ddsresults_TEs_batch_table,
                              by="TEs")
colnames(INTERGENIC_results)[c(2,3)] <- c("lfcMAP", "lfcMLE")
# Import annotation table of TEs
rep_annotation <- read_tsv(file="rollup_annotation/repName_repFamily_repClass_map.tsv")
# Filtering for Classes of interest
rep_annotation_filtered <- rep_annotation[rep_annotation$repClass %in% c("DNA", "Retroposon",
                                                                         "LINE", "SINE", "LTR"),]
# Filtering of intergenic counts for TEs belonging to classes of interest
INTERGENIC_results_filtered <- INTERGENIC_results[INTERGENIC_results$TEs %in% rep_annotation_filtered$repName,]
# Filtering for padj threshold
INTERGENIC_results_filtered <- INTERGENIC_results_filtered[INTERGENIC_results_filtered$padj <= 0.05,]
# Remove NAs
INTERGENIC_results_filtered <- na.omit(INTERGENIC_results_filtered) # 53

#----------------- Fig1-------------------------------------------------------------------------------------------
# Volcano plot
library(EnhancedVolcano)

labels_volcano <- INTERGENIC_results[INTERGENIC_results$TEs %in% rep_annotation_filtered$repName, ] %>% pull(TEs)

EnhancedVolcano(INTERGENIC_results[INTERGENIC_results$TEs %in% rep_annotation_filtered$repName,],
                lab = labels_volcano,
                selectLab = INTERGENIC_results_filtered$TEs,
                x = 'lfcMAP',
                y = 'padj',
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 1,
                labSize = 3.0,
                cutoffLineWidth = 0.2,
                drawConnectors = TRUE)



## PCA plot for intergenic TEs
vsd <- vst(Intergenic_allsamples_dds, blind=FALSE)
plotPCA(vsd, intgroup=c("type"))

# Barplot downregulated
# LTR
classe_barplot_LTR <- rep(c("LTR"), 4)
family_barplot_LTR <- c("ERV1", "ERVL", "ERVK", "Gypsy")
numero_barplot_LTR <- c(10, 1, 1, 3 )
LTR_barplot <- data.frame(classe_barplot_LTR, family_barplot_LTR, numero_barplot_LTR)
# DNA
classe_barplot_DNA <- rep(c("DNA"), 5)
family_barplot_DNA <- c("MULE-MuDR", "hAT?", "TcMar-Tigger", "hAT-Tip100", "hAT-Charlie")
numero_barplot_DNA <- c(1, 1, 4, 2, 2 )
DNA_barplot <- data.frame(classe_barplot_DNA, family_barplot_DNA, numero_barplot_DNA)
# LINE
classe_barplot_LINE <- c("LINE")
family_barplot_LINE <- c(0)
numero_barplot_LINE <- c(0)
LINE_barplot <- data.frame(classe_barplot_LINE, family_barplot_LINE, numero_barplot_LINE)
# SINE
classe_barplot_SINE <- c("SINE")
family_barplot_SINE <- c("Alu")
numero_barplot_SINE <- c(1)
SINE_barplot <- data.frame(classe_barplot_SINE, family_barplot_SINE, numero_barplot_SINE)
# All
nomi_tab <-list(LTR_barplot, DNA_barplot, LINE_barplot, SINE_barplot)
tab_rename<-lapply(nomi_tab, function(x){
  colnames(x)<-c("classe_barplot","family_barplot","numero_barplot")
  return(x)
})
tab_rename<-data.table::rbindlist(tab_rename)

# Stacked BARPLOT
ColourPalleteMulti <- function(df, group, subgroup){
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}
# Create data
df <- tab_rename
df$group <- paste0(df$class_barplot, "-", df$family_barplot, sep = "")

# Build the colour pallete
colours <-ColourPalleteMulti(df, "class_barplot", "family_barplot")
# Plot results
df %>% mutate(class_barplot =fct_relevel(class_barplot, "LINE", "LTR", "DNA", "SINE")) %>%
  ggplot(aes(x=class_barplot)) +
  geom_bar(aes(fill = group, y = numero_barplot),stat ='identity', position = "stack", colour = "grey27") +
  scale_fill_manual("Class_Family", values=colours) + ylab("N. of TEs")+ xlab("Class")+
  theme_minimal()



# ----------------------------- FigS1 -----------------------------------------------------------
# Heatmap T/N ratio
INTERGENICExp_TEs_normal <- INTERGENICExp_TEs[rownames(INTERGENICExp_TEs) %in% normal_samples,]
INTERGENICExp_TEs_tumour <- INTERGENICExp_TEs[rownames(INTERGENICExp_TEs) %in% tumour_samples,]

# Check colnames
all(colnames(INTERGENICExp_TEs_normal) == colnames(INTERGENICExp_TEs_tumour))

# Log2 Transformation
INTERGENICExp_TEs_normal <- INTERGENICExp_TEs_normal +1
INTERGENICExp_TEs_tumour <- INTERGENICExp_TEs_tumour +1
INTERGENICExp_ratio <- log2(INTERGENICExp_TEs_tumour/INTERGENICExp_TEs_normal)

### Heatmap annotations
# Tissue type annotation
Type_anno <- coldata_matched$Type
names(Type_anno) <- rownames(Type_anno)

anno_row_ratio = rowAnnotation(Type = Type_anno[rownames(INTERGENICExp_ratio)],
                               col = list(Type = c(Tumour= "darkred4", Normal="darkgreen")),
                               simple_anno_size = unit(2, "mm"))


# LFC annotation
LogFC_anno <- INTERGENIC_results_filtered$lfcMAP
names(LogFC_anno) <- INTERGENIC_results_filtered$TEs

anno_col_ratio <- HeatmapAnnotation(Log2FC = LogFC_anno[colnames(INTERGENICExp_ratio)],
                 col= list(Log2FC=col_fun_logFC),  gp = gpar(col = "grey37"))

# Heatmap
heatmap <-Heatmap(INTERGENICExp_ratio, cluster_rows = F,
                  name = "ratio T/N", row_names_gp = gpar(fontsize = 5), right_annotation = anno_row_ratio,
                  bottom_annotation = anno_col_ratio, col=zscore_col_fun,
                  width = unit(15, "cm"), height = unit(15, "cm"),
                  column_names_gp = gpar(fontsize=5), cell_fun = function(j, i,x,y, width, height, fill) {
                  })
draw(heatmap)

col_fun_logFC = colorRamp2(c(-4, -2, 0, 2, 4), c("royalblue4", "royalblue","white", "red", "red4"))


# Heatmap z-score counts grouped by patient_ID
zscore_counts <- scale(INTERGENICExp_TEs)

# Annotazione media tumour samples
tumour_samples <- coldata_matched[coldata_matched$Type=="Tumour", ] %>% pull(rownames(tumor_samples))
mean_counts_tumour <- INTERGENICExp_TEs[rownames(INTERGENICExp_TEs) %in% tumour_samples,]
mean_counts_tumour <- apply(mean_counts_tumour, 2, function(x){mean(x)})
mean_counts_tumour <- as.matrix(mean_counts_tumour)
mean_counts_tumour <- log2(mean_counts_tumour) %>% round()
colnames(mean_counts_tumour) <- "Mean_tumour"

# Duplicate matrix
zscore_counts_paste <- zscore_counts %>% t()
DETE_rep_annotation$name_paste <- paste(DETE_rep_annotation$repFamily, DETE_rep_annotation$repName, sep = ":")
# Check
all(rownames(DETE_rep_annotation) == rownames(zscore_counts_paste))
# Assign rep_family:rep_name as matrix rownames
rownames(zscore_counts_paste) <- DETE_rep_annotation$name_paste
# Check
all(rownames(mean_counts_tumour)==rownames(mean_counts_normal))
# Bind vectors with mean expression for normal and tumour samples
anno_points_bind <- cbind(mean_counts_normal, mean_counts_tumour)
all(rownames(anno_points_bind) == colnames(zscore_counts))
# Row annotation
INTERGENIC_risultati_filtered_order <- INTERGENIC_risultati_filtered[order(INTERGENIC_risultati_filtered$lfcMAP),]
zscore_row_anno = rowAnnotation(Type = Type_vector_v[colnames(zscore_counts_paste)],
                                col = list(Type = c(Tumour= "darkred", Normal="royalblue3")), simple_anno_size = unit(2, "mm"))
# Column annotation
zscore_col_anno = HeatmapAnnotation(Log2FC = LogFC_vector_int[colnames(zscore_counts)],
                                    log2_Mean = anno_lines(anno_points_bind, gp = gpar(col = c("royalblue3", "darkred")),
                                                           add_points = TRUE, pt_gp = gpar(col = c("royalblue3", "darkred")), pch = c(16, 16)),
                                    height = unit(2, "cm"),
                                    col= list(Log2FC=col_fun_logFC),  gp = gpar(col = "grey37"))

# Cell colour function
zscore_col_fun = colorRamp2(c(-2, 0, 2), c("royalblue4", "white", "red"))

# Heatmap
heatmap <-Heatmap(t(zscore_counts_paste),
                  name = "z-score", row_names_gp = gpar(fontsize = 5), right_annotation = zscore_row_anno,
                  bottom_annotation = zscore_col_anno, column_split = 2,
                  col=zscore_col_fun,
                  width = unit(15, "cm"), height = unit(15, "cm"),
                  row_split = coldata_matched$Patient_ID,
                  border = T,
                  column_names_gp = gpar(fontsize=5), cell_fun = function(j, i,x,y, width, height, fill) {
                  })
draw(heatmap)

## GSVA analysis
# List signatures Msigdb hallmark
m_df_hallmark = msigdbr(species = "Homo sapiens", category = "H")
m_list_hallmark = m_df_hallmark %>% split(x = .$human_gene_symbol, f = .$gs_name)
# List custom innate immune signatures
immune_sig <- list(Response_typeI_all, Paper_Response_CD8Teffector, Paper_Response_IL1B, Paper_Response_NFKB,
                           Paper_Response_typeII)
names(immune_sig) <- c("Response_typeI_all", "Paper_Response_CD8Teffector", "Paper_Response_IL1B", "Paper_Response_NFKB",
                                "Paper_Response_typeII")
# GSVA
Gsva_hallmark = gsva(as.matrix(dds_SCLC_Exp), m_list_hallmark, method="gsva", verbose=TRUE)
Gsva_hallmark <- t(Gsva_hallmark)
Gsva_sigs <- gsva(as.matrix(dds_SCLC_Exp), immune_sig, method="gsva", verbose=TRUE)
Gsva_sigs <- t(Gsva_sigs)

# ------------------- Fig2 ----------------------------------------------------------------------------
# Correlation between gsva enrichment scores and DETEs expression
INTCor_exp_TEs_pathway_allsamples <- apply(Intergenic_Exp_allsamples_dds,2,function(x)
{apply(Gsva_hallmark[rownames(Intergenic_Exp_allsamples_dds),], 2, function(y){cor.test(x, y, method="spearman")["estimate"]})})
INTCor_exp_TEs_pathway_allsamples <- lapply(INTCor_exp_TEs_pathway_allsamples, function(x){unlist(x)})
INTCor_exp_TEs_pathway_allsamples <- as.data.frame(INTCor_exp_TEs_pathway_allsamples)
colnames(INTCor_exp_TEs_pathway_allsamples) <- gsub(pattern = '\\.', replacement =  "-", x=colnames(INTCor_exp_TEs_pathway_allsamples))
INTCor_exp_TEs_pathway_allsamples <- as.matrix(INTCor_exp_TEs_pathway_allsamples)

# Transpose matrix
INTCor_exp_TEs_pathway_allsamples_t <- t(as.matrix(INTCor_exp_TEs_pathway_allsamples))

# Order matrix by lfcMAP
INTCor_exp_TEs_pathway_allsamples_t_ordered <- INTCor_exp_TEs_pathway_allsamples_t[order(match(rownames(
  INTCor_exp_TEs_pathway_allsamples_t), INTERGENIC_results_filtered_order$TEs)), , drop=FALSE ]

# Duplicate correlation matrix
INTCor_exp_TEs_pathway_allsamples_paste <- INTCor_exp_TEs_pathway_allsamples_t_ordered

# Check
all(rownames(DETE_rep_annotation_paste) == rownames(INTCor_exp_TEs_pathway_allsamples_paste))

# Assign rep_family:rep_name as matrix rownames
rownames(INTCor_exp_TEs_pathway_allsamples_paste) <- DETE_rep_annotation_paste$name_paste

# Row annotation
Class_anno <- DETE_rep_annotation$repClass
names(Class_anno) <- rownames(Class_anno)

ha_rowanno = rowAnnotation(Class = Class_anno[rownames(INTCor_exp_TEs_pathway_allsamples_t_ordered)],
                           col = list(Class = c(SINE = "indianred1", LTR = "lightgoldenrod", DNA = "skyblue2", LINE= "mediumseagreen"),
                                      Log2FC = col_fun_logFC),
                           Log2FC = LogFC_vector_int[rownames(INTCor_exp_TEs_pathway_allsamples_t_ordered)],
                           simple_anno_size= unit(0.3, "cm"))

# Get p-values
pval_mat_allsamples <- apply(Intergenic_Exp_allsamples_dds[2:length(colnames(Intergenic_Exp_allsamples_dds))],2,function(x)
{apply(Gsva_hallmark[rownames(Intergenic_Exp_allsamples_dds),], 2, function(y){cor.test(x, y, method="spearman")["p.value"]})})
pval_mat_allsamples <- lapply(pval_mat_allsamples, function(x){unlist(x)})
pval_mat_allsamples <- as.data.frame(pval_mat_allsamples)
colnames(pval_mat_allsamples) <- gsub(pattern = '\\.', replacement =  "-", x=colnames(pval_mat_allsamples))
pval_mat_allsamples <- as.data.frame(pval_mat_allsamples)
# Create matrix with p-value asterisk
pval_mat_asterisk_allsamples <- pval_mat_allsamples
pval_mat_asterisk_allsamples[pval_mat_allsamples > 0.05] <- ""
pval_mat_asterisk_allsamples[pval_mat_allsamples>(1e-3) & pval_mat_allsamples<=(0.05)] <- "*"
pval_mat_asterisk_allsamples[pval_mat_allsamples<=(1e-3) & pval_mat_allsamples>(1e-5)] <- "**"
pval_mat_asterisk_allsamples[pval_mat_allsamples<=(1e-5)] <- "***"
# Traspose matrix
pval_mat_asterisk_allsamples <- t(as.matrix(pval_mat_asterisk_allsamples))

# Check
all(rownames(pval_mat_asterisk_allsamples) == rownames(INTCor_exp_TEs_pathway_allsamples_t_ordered))

# Heatmap
heatmap <-Heatmap(INTCor_exp_TEs_pathway_allsamples_paste,
                  name = "Cor", col=col_fun, row_names_gp = gpar(fontsize = 5), right_annotation = ha_rowanno,
                  width = unit(18, "cm"), height = unit(18, "cm"),
                  row_order = rownames(INTCor_exp_TEs_pathway_allsamples_paste),
                  column_names_gp = gpar(fontsize=5),
                  cell_fun = function(j, i,x,y, width, height, fill) {
                  grid.text(pval_mat_asterischi_allsamples_ordered[i, j],x,y, gp=gpar(fontsize=6))})

draw(heatmap)


#------------------- Fig3 --------------------------------------

## Correlation between DETE expression and innate sensors

# Select innate sensors
sensors_corr <- dds_Exp_SCLCallsamples %>% t() %>% as.data.frame() %>%
  select(c("TLR7", "TLR8", "TLR3", "IFIH1", "STING1", "IFI16", "AIM2", "DDX58", "CGAS")) %>%
  rownames_to_column(var = "samples")
# Rownames to column for GSVA signatures samples
Gsva_sigs <- Gsva_sigs %>% as.data.frame() %>% rownames_to_column( var = "samples")
# Rownames to column for TEs samples
Intergenic_Exp_allsamples_dds <- Intergenic_Exp_allsamples_dds %>% rownames_to_column( var = "samples")
# Merge
scatter_correlation <-  inner_join(Gsva_sigs, sensors_corr, by = "samples") %>%
  inner_join(Intergenic_Exp_allsamples_dds, by = "samples") %>%
  inner_join(Gsva_hall, by="samples")
scatter_correlation[,c(7:72)]<-log2(scatter_correlation[,c(7:72)]+1)
scatter_correlation <- as.data.frame(scatter_correlation)
# Create groups for Response to type I interferon
samples_up <- Gsva_sigs %>% filter(Response_typeI_all > 0.2) %>% pull(samples)
samples_down <- Gsva_sigs %>% filter(Response_typeI_all < -0.2) %>% pull(samples)
samples_up_down <- c(samples_up, samples_down)

# Plot
give.n <- function(x){
  return(c(y = max(x) + 1.7, label = length(x)))
}

library(ggpubr)
fun_corr_box <- function(TE, genes, pathway){
  maxval_TE <- scatter_correlation %>% filter(samples %in% samples_up_down) %>%
    mutate(score = ifelse(samples %in% samples_up, "up", "down")) %>% select(c(TE)) %>% max()
  maxval_genes <- scatter_correlation %>% filter(samples %in% samples_up_down) %>%
    mutate(score = ifelse(samples %in% samples_up, "up", "down")) %>% select(c(genes)) %>% max()
  box_new_TE <- scatter_correlation %>% filter(samples %in% samples_up_down) %>%
    mutate(score = ifelse(samples %in% samples_up, "up", "down")) %>% select(c("samples", TE, "score")) %>%
    ggviolin(y = TE, x = "score", fill = "score", palette = c("red4","cornflowerblue"), orientation ="horizontal", color = "black",
             add = "boxplot",
             add.params = list(fill = "white"), ylab = paste("log2(", TE,"+1) Normalized Counts"))+
    stat_compare_means(label = "p.signif", inherit.aes = T, size = 4, method = "wilcox.test", label.x = 1.3, label.y=(maxval_TE+1.8), label.x.npc = 0.3)+
    stat_summary(fun.data = give.n, aes(y = maxval_TE+1), col="black", size = 3, geom = "text")+
    theme_light() +
    rremove("legend") + rremove("ylab")

  box_new_gene <- scatter_correlation %>% filter(samples %in% samples_up_down) %>%
    mutate(score = ifelse(samples %in% samples_up, "up", "down")) %>% select(c("samples", genes, "score")) %>%
    ggviolin(y = genes, x = "score", fill = "score", palette = c("red4","cornflowerblue"), color = "black", add = "boxplot",
             add.params = list(fill = "white"), ylab = paste("log2(", genes,"+1) Normalized Counts"))+
    stat_compare_means(label = "p.signif", inherit.aes = T, size = 4, method = "wilcox.test", label.y = (maxval_genes+1.8))+
    stat_summary(fun.data = give.n, aes(y = maxval_genes+1), col="black", size = 3, geom = "text")+
    theme_light()+
    rremove("legend")+rremove("xlab")

  scatterPlot <- ggscatter(data = scatter_correlation, y = genes, x = TE, fill = pathway , conf.int = TRUE, cor.coef = TRUE,  add = "reg.line",
                           cor.coeff.args = list(method = "spearman", label.x = 0, label.sep = "\n"), shape = 21, size = 4, color = "gray",
                           add.params = list(color = "black", fill = "lightgray", size = 0.3), xlab = paste("log2(",TE,"+1) Normalized Counts"),
                           ylab = paste("log2(", genes,"+1) Normalized Counts"))+
    gradient_fill(c("cornflowerblue","white", "red4"))+
    theme_light()+
    theme(
      legend.position = c(1.2,1.1),
      legend.margin = margin(1.1,1.1,1.1,1.1), legend.key.size = unit(0.5, 'cm'), legend.direction = "vertical", legend.justification = "bottom" )

  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank() )
  library("gridExtra")
  return(grid.arrange(box_new_TE, blankPlot, scatterPlot,box_new_gene ,
                      ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4)))}
fun_corr_box("LTR22C", "DDX58","Response_typeI_all")


# PCA plots for DETEs

library(ggfortify)

list_te <- list()
a=1

for (i in c("LTR30", "LTR9C", "LTR22C", "MER61F", "MER74C")){
  list_te[[a]] <- tab_exp_pca %>%
    column_to_rownames("samples") %>%
    select("TLR7", "TLR8", "TLR3", "IFIH1", "STING1", "IFI16", "AIM2", "DDX58", "CGAS", i) %>%
    prcomp(scale. = TRUE) %>%
    autoplot(data = tab_exp_pca, fill="Response_typeI_all", loadings = TRUE, size=3,
             loadings.label = TRUE, loadings.label.size = 3, loadings.colour = "gray", shape=21, col="gray",
             loadings.label.repel=T)+
    scale_fill_gradient2( high="red4", mid = "white", low= "cornflowerblue") + theme_bw() +
    labs(title = i, fill="Response to type I Interferon")
  names(lista_te)[a] <- i
  a=a+1
}

do.call(grid.arrange, list_te)

#------------------------------ Fig4 -----------------------------------
library("LOLA")

## Universe set
# BED file for expressed intergenic TEs
Intergenic_expressed <- INTERGENIC_ddsresults_TEs_batch_table[INTERGENIC_ddsresults_TEs_batch_table$TEs %in%
                                                                rep_annotation_filtered$repName, 1]
write.table(x = Intergenic_expressed, file = "Intergenic_expressed", quote = F, sep = " ", row.names = F,
            col.names = F)


# Filter LTRs
REP_LTR <- rep_annotation_filtered %>% filter(repClass == "LTR") %>% pull(repName)
intersect_universe_LTR <- intersect_universe_table %>% filter(X4 %in% REP_LTR)
# Save as bed file
write.table (intersect_universe_LTR, file = "intersect_universe_LTR.bed", col.names = F,
             quote = F, sep ="\t",row.names = F)

# Import coordinates for intergenic TEs
coord_intergenic <- readRDS(file= "~/NEW_TEs/rollup_annotation/rmsk_annotation.RDS")

# Create a list of upregulated LTRs to be pasted into UCSC Table browser to obtain coordinates and strandness
list_allLTR_upreg <- df_up %>% filter(repClass == "LTR") %>% pull(repName)
write.table(x = list_allLTR_upreg, file = "list_allLTR_upreg", quote = F, sep = " ", row.names = F, col.names = F)

# Import BED file from UCSC
list_allLTR_upreg <- read_tsv("/home/administrator/NEW_TEs/BED_from_UCSC/lista_allLTR_upregolatiBED", col_names = F)

# Select coordinates for upregulated LTRs and split idx
coord_intergenic_filtered_upLTR <- coord_intergenic %>%
  filter(selected_feature == "intergenic" & repName %in% list_allLTR_upreg$X4) %>%
  select(-md5) %>% separate_rows(idx, sep = ",", convert = T)

# Uniform start coordinates between UCSC BED and rep_annotation (0-start)
colnames(coord_intergenic_filtered_upLTR)[c(7,8,9)] <- c("chr","genoStart", "genoEnd")
coord_intergenic_filtered_upLTR$genoStart <- coord_intergenic_filtered_upLTR$genoStart -1
coord_intergenic_filtered_upLTR$genoEnd <- coord_intergenic_filtered_upLTR$genoEnd
# Merge with UCSC bed file to add strandness
Merge_UCSC_LTR <- inner_join(coord_intergenic_filtered_upLTR, lista_allLTR_upregolati[c("genoStart", "strand")],
                             by="genoStart")
# Remove duplicates
Merge_UCSC_LTR <- unique(Merge_UCSC_LTR)

### LOLA
# Create a bed file to add score variable and rename chr variable name
Merge_UCSC_LTR <- Merge_UCSC_LTR %>%
  mutate(score = rep(1, nrow(Merge_UCSC_LTR))) %>%
  select("chr", "genoStart", "genoEnd", "repName", "score", "strand")
Merge_UCSC_LTR$chr <- paste0("chr", Merge_UCSC_LTR$chr)

# Save in a .bed format
write.table (Merge_UCSC_LTR, file="Merge_UCSC_LTR.bed", col.names = F, quote = F, sep ="\t", row.names = F)

# Load regions from Lola Database
regionDB = loadRegionDB("LOLACore/hg38/")

# Split Granges list for each LTR
regionSetUP_LTR <- readBed("Merge_UCSC_LTR.bed")
regionSetUP_LTR<- split(regionSetUP_LTR,as.factor(regionSetUP_LTR@ranges@NAMES))

# RUN LOLA
lolaUP_LTR = runLOLA(regionSetUP_LTR, intersect_universe_all_LTR, regionDB, cores=8)

## Plot Cistrome Epigenome
# Filter cistrome epigenome data
lola_up_LTR_epigenome <- lolaUP_LTR %>% filter(collection %in%  "cistrome_epigenome")  %>%
  mutate(marker = ifelse(antibody == "NA", filename, antibody)) %>%
  mutate(perc_support = (support/(support+c)*100)) %>%
  select(c("pValueLog", "userSet", "new", "collection", "perc_support")) %>%
  group_by(userSet) %>% slice_max(pValueLog, n = 10, with_ties = F) %>% ungroup() %>%
  mutate(ranking = rep(1:10, 19)) %>% select(-"collection")

# Marker variable as factor
lola_up_LTR_epigenome$marker<-as.factor(lola_up_LTR_epigenome$marker)

# Plot
ggplot(lola_up_LTR_epigenome, aes())+
  geom_point(aes(x = userSet, y = rev(ranking), size = perc_support,fill= pValueLog),pch=21,color="grey24")+
  geom_label(aes(x = userSet, y = rev(ranking),label = marker, color = marker), size = 3, show.legend = F,
             nudge_y =0.4)+
  scale_fill_gradient2(high = "purple4", low = "darkred", midpoint = 2, limits=c(0,20), oob=squish)+
  scale_size(range = c(1.5, 12))+
  theme_ipsum() + colScale +
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  labs(title = "Cistrome Epigenome", x = "Upregulated LTRs", y= "Ranking", fill= "Log(p-value)", size="% support")+
  theme(axis.text.y = element_blank(), axis.ticks = element_blank())


# GSEA analysis for LTR30 and Ridgeplot for downregulated TEs
library(clusterProfiler)
library(DESeq2)
library(tidyverse)
#####################################################
####in Shell
## Use closest bedtools to get the genes closest to downregulated TEs
# Sort the bed file of downregulated TEs
sort -k1,1 -k2,2n bed_intergenic_down.bed > sorted_bed_intergenic_down.bed
# Sort the bed file of genes
sort -k1,1 -k2,2n bed_genes.bed > sorted_bed_genes.bed
#closest bedtools line command
bedtools closest -s -iu -D a -t first -a bed_intergenic_down.bed -b sorted_bed_genes.bed > output_down.bed
###################################################
# Filter genes downstream up to 20kb from TEs as TERM2GENE
DOWN_intergenic_genes.Select <- output_down %>% mutate(X13 = abs(X13)) %>% filter(X13 < 20000) %>%
  select(c(4,10))
colnames(DOWN_intergenic_genes.Select)<-c("Tes","genes")
# List of genes with stat
Genes_SCLC <- results(dds_SCLC_matched) %>% as.data.frame %>% rownames_to_column() %>%
  dplyr::select(rowname, stat) %>%
  na.omit() %>% distinct() %>% deframe()

# GSEA
GSEA.DOWN_intergenic<- GSEA(sort(Geni_SCLC, decreasing = T),minGSSize = 0,
                            maxGSSize = Inf, TERM2GENE = DOWN_intergenic_genes.Select,
                            pvalueCutoff = 1,verbose = T,by = "fgsea",eps = 0,seed = 1)
# Select LTR30
DOWN_intergenic_genes.Select[DOWN_intergenic_genes.Select$gs_name=="LTR30",]

# Gsea plot
gseaplot(GSEA.DOWN_intergenic,"LTR30",title = "LTR30 proximal genes")
# Ridgeplot plot
ridgeplot(GSEA.DOWN_intergenic,fill="pvalue",decreasing=T)


#--------------Fig5 ----------------------------

Survival_untreated <- Survival %>% filter(`previous therapeutic treatment for SCLC` == "untreated")%>%
  as_tibble() %>% dplyr::rename(sample_ID ="Sample-ID",
                                treatment = "previous therapeutic treatment for SCLC",
                                chemotherapy = "chemotherapy (yes/no)",
                                prog_free_surv_months = "progression-free_survival (months)",
                                overall_surv_months = "overall_survival (months)",
                                radiation = "radiation (yes/no)",
                                Survival_status = "Status (at time of last follow-up)",
                                tissue_sampling = "tissue sampling")

# Add variable of LTR30 exp
ALL_SAMPLES_LTR30 <- Intergenic_Exp_allsamples_dds[, c("samples", "LTR30")]
ALL_SAMPLES_LTR30 <- ALL_SAMPLES_LTR30 %>% dplyr::rename(sample_ID = "samples")
# Add variable of exp levels
ALL_SAMPLES_LTR30[ALL_SAMPLES_LTR30$LTR30 <= 4.507529 , "exp"] <- "Low"
ALL_SAMPLES_LTR30[ALL_SAMPLES_LTR30$LTR30 >= 5.70731, "exp"] <- "High"


LTR30_table_survival <- Survival_untreated %>% as.tibble() %>% # select untreated samples
  inner_join(ALL_SAMPLES_LTR30, by = "sample_ID") %>% # add LTR30 expression variable
  mutate(LTR30 = log2(LTR30+1)) %>% # transform LTR30 expression
  filter(chemotherapy != "NA") %>%
  mutate(chemotherapy = factor(chemotherapy)) %>% # chemotherapy as factor
  filter( Survival_status != "NA") %>%
  mutate(Survival_status = if_else(Survival_status == "alive", "0", "1")) %>% # change survival annotation
  filter(overall_surv_months != "NA") %>%
  mutate(Survival_status = as.numeric(Survival_status), # survival as numeric
         overall_surv_months = as.numeric(overall_surv_months)) %>%
 filter(LTR30 <= quantile(LTR30,0.25)| LTR30 >= quantile(LTR30, 0.75) ) # filter LTR30 samples based on exp levels

# Survival plot

library("survival")
library("survminer")

fit_LTR30_new <- survfit(Surv(overall_surv_months, Survival_status) ~ exp + chemotherapy, data=LTR30_table_survival)

ggsurvplot(fit_LTR30_new,
           pval = TRUE, conf.int = FALSE,
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           color = "strata")

# ------------- FigS7 ---------------------------------------------------------

# Paired boxplot

SCLCExp_matched_box <- SCLCExp_matched %>% t() %>% as.data.frame
SCLCExp_matched_box <- SCLCExp_matched_box %>% rownames_to_column("samples")
coldata_matched_box <- coldata_matched %>% rownames_to_column("samples")
merge <- merge(SCLCExp_matched_box, coldata_matched_box[,c("Type", "samples")], by= "samples")

genes_boxplot <- function(gene){
  tab_all <- merge[, c("Type", gene)]
  p <- ggplot(data = tab_all, aes(x=tabella_all[,1], y= log2(tab_all[,2]+1)))+
    geom_boxplot(aes(fill=Type), colour= "black", width = 0.5)
  p <- p + ylab(paste0(gene, " ", "log2(Normalized Counts)")) + xlab("SCLC")
  p <- p + scale_fill_manual(values = c("cornflowerblue", "darkred"))
  p <- p + theme_minimal ()
  p <- p + theme (text = element_text(size = 12))

  return(p)

}

box_plot<-c("ESR1", "KDM1A", "GATA6", "STAG1", "RAD21", "TRIM24")
box_plot<-lapply(box_plot, function(x){
  p<-genes_boxplot(x)
  return(p)
})
do.call(grid.arrange,box_plot)

# Paired t-test

library(PairedData)

t_test_genes <- function(gene){
  tab_genes <- merge[,c("Type", gene)]
  Tumour_vector <- subset(tab_genes, Type == "Tumour", gene, drop=TRUE)
  Normal_vector <-  subset(tab_genes, Type == "Normal", gene, drop=TRUE)
  t_test_res <- t.test(Tumour_vector, Normal_vector, paired = TRUE, alternative = "two.sided")
  return(t_test_res)
}

t_test_genes ("ESR1")
