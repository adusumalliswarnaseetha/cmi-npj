#################################################################################
### STutility for all timepoints (1, 4 and 12 weeks)
library(STutility)
library(ggplot2)
library(tidyverse)
library(biomaRt)
library("org.Hs.eg.db")
library("org.Ss.eg.db")
library("clusterProfiler")
library(clustree)

load("./Chronic_All_Analysis.Rdata")
infoTable <- read.table("infotable.txt",
                        header=T, sep="\t", stringsAsFactors = FALSE)
############################ QUALITY CONTROL ############################
test <- infoTable[c(3,4,7,8,11,12) ,]
se <- InputFromTable(infotable = test, platform =  "Visium", disable.subset = TRUE)
dim(as.data.frame(se@assays$RNA@counts)) #[32180 genes, 19334 spots]
############################ PER SPOT STATS ############################
col_attr <- data.frame(nUMI = Matrix::colSums(se@assays$RNA@counts), 
                       ngene = Matrix::colSums(se@assays$RNA@counts > 0))
quantile(col_attr$ngene, probs = c(0.05))
quantile(col_attr$nUMI, probs = c(0.05))
###########################  PER GENE STATS ############################
gene_attr <- data.frame(nUMI = Matrix::rowSums(se@assays$RNA@counts), 
                        nSpots = Matrix::rowSums(se@assays$RNA@counts > 0))
quantile(gene_attr$nSpots, probs = c(0.1))
quantile(gene_attr$nUMI, probs = c(0.1))
######################################
se <- InputFromTable(infotable = test,
                     minSpotsPerGene =20,
                     minUMICountsPerGene=4,
                     minGenesPerSpot=303,
                     minUMICountsPerSpot=434,
                     platform =  "Visium")
st.object <- se@tools$Staffli
dim(st.object[[]])
dim(as.data.frame(se@assays$RNA@counts))   

############################  Mitochrondrial genes ############################ 
mt.genes <- grep(pattern = ("^GRCh38-MT-|^GRCh38-NDUF-|^GRCh38-COX-|^Ssus11-MT-|^Ssus11-ND-|^Ssus11-COX-|^Ssus11-NDUF-"), x = rownames(se), value = TRUE)
se$percent.mito <- (Matrix::colSums(se@assays$RNA@counts[mt.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100
length(which(se$percent.mito >10))
ST.FeaturePlot(se, features = "percent.mito", dark.theme = F, cols = c("dark blue", "cyan", "yellow", "red", "dark red"), ncol=4)

rp.genes <- grep(pattern = "^GRCh38-RPS|^GRCh38-MRPS|^Ssus11-RPS|^Ssus11-MRPS", x = rownames(se), value = TRUE)
se$percent.ribo <- (Matrix::colSums(se@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100
rp.genes <- grep(pattern = "^GRCh38-RPL|^GRCh38-MRPL|^Ssus11-RPL|^Ssus11-MRPL", x = rownames(se), value = TRUE)
se$percent.ribo <- (Matrix::colSums(se@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

ST.FeaturePlot(se, features = "percent.ribo", dark.theme = F, cols = c("dark blue", "cyan", "yellow", "red", "dark red"), ncol=4)
se@meta.data$slide_id2 <- factor(se@meta.data$slide_id, levels = unique(se@meta.data$slide_id))
VlnPlot(se, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mito"), group.by = "slide_id2", pt.size=0, ncol=2)
#VlnPlot(se, features = c("GRCh38-ACTN2", "GRCh38-GJA1"), group.by = "slide_id", pt.size=0, assay = 'scale.data')
se <- SubsetSTData(se, expression = percent.mito < 10)
cat("Spots: ", ncol(se), "\n")
st.object <- se@tools$Staffli
dim(st.object[[]])
se.section1 <- se
dim(as.data.frame(se.section1@assays$RNA@counts))
ST.FeaturePlot(se.section1, features = "percent.mito", dark.theme = F, cols = c("dark blue", "cyan", "yellow", "red", "dark red"), ncol=2)
ST.FeaturePlot(se.section1, features = "percent.ribo", dark.theme = F, cols = c("dark blue", "cyan", "yellow", "red", "dark red"), ncol=2)
VlnPlot(se, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mito"), group.by = "slide_id2", pt.size=0, ncol=2)

###################################### Normalization ######################################
se.section1 <- LoadImages(se.section1, time.resolve = F, verbose = T)
rm(se)
ImagePlot(se.section1, method = "raster", darken = F, type = "raw", ncols = 2)

se.section1 <- SCTransform(se.section1,  variable.features.n = 3000)
se.section1$section <- paste0("section_", GetStaffli(se.section1)[[, "sample", drop = T]])

########################################## NMF ##########################################
se.section1 <- RunNMF(se.section1, nfactors = 14, n.cores=10)
cscale <- c("darkblue", "cyan", "yellow", "red", "darkred")

ST.DimPlot(se.section1, dims = 1:14, ncol = 6,  grid.ncol = 2,
           reduction = "NMF", dark.theme = F,
           pt.size = 1, center.zero = F, cols = cscale)

FactorGeneLoadingPlot(se.section1, factor = 14, topn = 30, dark.theme = F)
###################################  Clustering ##########################################
cluster_resolutions=seq(0.2,0.9,by=0.1)
se.section1 <- FindNeighbors(object = se.section1, verbose = FALSE, reduction = "NMF", dims = 1:14)
se.section1 <- FindClusters(object = se.section1, verbose = FALSE, resolution=cluster_resolutions)

clustree(se.section1, prefix= "SCT_snn_res.")
se.section1 <- FindClusters(object = se.section1, verbose = FALSE, resolution = 0.4)

#cols = col_vector, 
col_vector = c("#B7CEEC", "#BEAED3","#FDC085","maroon","#FFC0CB", "#FFE4C4","#BEAED4", "magenta","chocolate","#F7E7CE", 
               "#7C9D8E", "#C0C0C0", "blue", "#E1D9D1",  "#D6DBDF")#, "linen"
ST.FeaturePlot(object = se.section1, features = "seurat_clusters", dark.theme = F, pt.size = 1, cols = col_vector,ncol = 2,show.sb = F)# min.cutoff='q1'
ST.FeaturePlot(object = se.section1, features = "seurat_clusters", dark.theme = F, pt.size = 1, split.labels = T, 
               indices = 6, show.sb = FALSE, ncol=4)

de.markers_se.section1 <- FindAllMarkers(se.section1)
de.markers_se.section1 <- de.markers_se.section1[de.markers_se.section1$p_val_adj<0.05 ,]
write.table(de.markers_se.section1,"Chronic_deMarkers_ALL_Clusters.txt",sep = "\t",quote = F)

tot_genes <- as.numeric(table(de.markers_se.section1$cluster))
human_genes <- de.markers_se.section1 %>% filter(grepl("GRCh38-", gene)) %>%
  group_by(cluster) %>% summarise(nHuman= n())
de.markers_se.section1_HUMAN <- de.markers_se.section1 %>% filter(grepl("GRCh38-", gene)) %>%
  group_by(cluster)
write.table(de.markers_se.section1_HUMAN,"Chronic_deMarkers_HUMAN_Clusters.txt",sep = "\t",quote = F)
human_genes <- (human_genes)
human_genes$totGenes <- tot_genes
human_genes$humPer <- round((human_genes$nHuman/human_genes$totGenes)*100,1)
write.table(human_genes,"Chronic_deMarkers_HUMAN_PERCENT.txt",sep = "\t",quote = F)

############################   DEG across days
de_1vs4 <- FindMarkers(se.section1,ident.1 = c(3), ident.2 = 12)
de_1vs4_sig <- de_1vs4[de_1vs4$p_val_adj<0.05 ,]
de_1vs4_sig$gene <- row.names(de_1vs4_sig)
de_1vs4_sig_human <- de_1vs4_sig %>% filter(grepl("GRCh38-", gene))
pos <- de_1vs4_sig_human[de_1vs4_sig_human$avg_log2FC >= 0 ,]
neg <- de_1vs4_sig_human[de_1vs4_sig_human$avg_log2FC < 0 ,]
write.table(de_1vs4_sig,"DE_4_vs_1_sig.txt",sep = "\t",quote = F)
write.table(pos,"DE_4_vs_1_sig_pos.txt",sep = "\t",quote = F)
write.table(neg,"DE_4_vs_1_sig_neg.txt",sep = "\t",quote = F)

de_4vs12 <- FindMarkers(se.section1,ident.1 = c(7,8), ident.2 = 3)
de_4vs12_sig <- de_4vs12[de_4vs12$p_val_adj<0.05 ,]
de_4vs12_sig$gene <- row.names(de_4vs12_sig)
de_4vs12_sig_human <- de_4vs12_sig %>% filter(grepl("GRCh38-", gene))
pos_4vs12 <- de_4vs12_sig_human[de_4vs12_sig_human$avg_log2FC >= 0 ,]
neg_4vs12 <- de_4vs12_sig_human[de_4vs12_sig_human$avg_log2FC < 0 ,]
write.table(de_4vs12_sig,"DE_12_vs_4_sig.txt",sep = "\t",quote = F)
write.table(pos_4vs12,"DE_12_vs_4_sig_pos.txt",sep = "\t",quote = F)
write.table(neg_4vs12,"DE_12_vs_4_sig_neg.txt",sep = "\t",quote = F)

de_1vs12 <- FindMarkers(se.section1,ident.1 = c(7,8), ident.2 = 12)
de_1vs12_sig <- de_1vs12[de_1vs12$p_val_adj<0.05 ,]
de_1vs12_sig$gene <- row.names(de_1vs12_sig)
de_1vs12_sig_human <- de_1vs12_sig %>% filter(grepl("GRCh38-", gene))
pos_1vs12 <- de_1vs12_sig_human[de_1vs12_sig_human$avg_log2FC >= 0 ,]
neg_1vs12 <- de_1vs12_sig_human[de_1vs12_sig_human$avg_log2FC < 0 ,]
write.table(de_1vs12_sig,"DE_12_vs_1_sig.txt",sep = "\t",quote = F)
write.table(pos_1vs12,"DE_12_vs_1_sig_pos.txt",sep = "\t",quote = F)
write.table(neg_1vs12,"DE_12_vs_1_sig_neg.txt",sep = "\t",quote = F)
############################  MARKERS  ############################
setwd("~/Documents/nBox/CM_ST_ALL_Chronic/FINAL_RESULTS/Factors_14/Markers/")
FeatureOverlay(se.section1, 
               features = c("GRCh38-NOS3"), 
               pt.size = 1.5,
               cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
               dark.theme = F, 
               type = "raw",
               sampleids = c(1,2,3,4,5,6),
               ncols = 2,
               min.cutoff = 'q1',
               show.sb=F
)
se.section2 <- subset(x = se.section1, idents = c(12,3,7))
se.section2@meta.data$seurat_clusters <- factor(se.section2@meta.data$seurat_clusters, levels = c(12,3,7))
VlnPlot(se.section2, features =c("GRCh38-CACNA1C","GRCh38-CACNB2","GRCh38-RYR2","GRCh38-KCNH2","GRCh38-KCNJ8","GRCh38-KCNMB4",
                                 "GRCh38-CD36","GRCh38-GJA1","GRCh38-ESRRA"), group.by = "seurat_clusters", ncol=3)
VlnPlot(se.section2, features =c("GRCh38-ACTC1","GRCh38-ACTN2","GRCh38-MYL2","GRCh38-NKX2-5","GRCh38-TNNT2","GRCh38-TPM1"), group.by = "seurat_clusters", ncol=3)
VlnPlot(se.section2, features =c("GRCh38-ACTC1","GRCh38-ACTN2","GRCh38-MYL2","GRCh38-NKX2-5","GRCh38-TNNT2","GRCh38-TPM1"), group.by = "seurat_clusters", ncol=3, pt.size=0)
VlnPlot(se.section2, features =c("Ssus11-LAMA1","Ssus11-LAMA2","Ssus11-LAMA3","Ssus11-LAMA4","Ssus11-LAMA5",
                                 "Ssus11-LAMB1","Ssus11-LAMB2","Ssus11-LAMB3","Ssus11-LAMB4",
                                 "Ssus11-LAMC1","Ssus11-LAMC2","Ssus11-LAMC3"), group.by = "seurat_clusters", ncol=4)
VlnPlot(se.section2, features =c("GRCh38-LAMA1","GRCh38-LAMA2","GRCh38-LAMA3","GRCh38-LAMA4","GRCh38-LAMA5",
                                 "GRCh38-LAMB1","GRCh38-LAMB2","GRCh38-LAMB3","GRCh38-LAMB4",
                                 "GRCh38-LAMC1","GRCh38-LAMC2","GRCh38-LAMC3"), group.by = "seurat_clusters", ncol=4)
VlnPlot(se.section2, features =c("GRCh38-MDK","GRCh38-SDC2","GRCh38-SDC4","GRCh38-LRP1","GRCh38-NCL","GRCh38-VEGFA","GRCh38-VEGFB","GRCh38-VEGFR1"), group.by = "seurat_clusters", ncol=4, pt.size=0)
VlnPlot(se.section2, features =c("GRCh38-POSTN","GRCh38-ITGV","GRCh38-ITGB5"), group.by = "seurat_clusters", ncol=2)
VlnPlot(se.section2, features =c("GRCh38-MYH6","GRCh38-MYH7","GRCh38-MYL7","GRCh38-MYL2","GRCh38-TNNI1","GRCh38-TNNI3"), group.by = "seurat_clusters", ncol=2)
VlnPlot(se.section1, features =c("Ssus11-PECAM1","GRCh38-VWF"), group.by = "slide_id2", ncol=2)

markers <- c(
  # "Ssus11-ACTC1","Ssus11-ACTN2",
  # "Ssus11-MYH7","Ssus11-MYL2",#"Ssus11-MYH6",
  # "Ssus11-MYL4",
  # "Ssus11-GJA1","Ssus11-NKX2-5",
  # "Ssus11-TNNI1","Ssus11-TNNI3","Ssus11-TNNC1",
  # "Ssus11-TNNT2","Ssus11-TPM1"#,"Ssus11-TTN"
  # "GRCh38-ACTC1","GRCh38-ACTN2",
  # "GRCh38-MYH6","GRCh38-MYH7","GRCh38-MYL2",
  # "GRCh38-MYL4",
  # "GRCh38-GJA1","GRCh38-NKX2-5",
  # "GRCh38-TNNI1","GRCh38-TNNI3","GRCh38-TNNC1",
  # "GRCh38-TNNT2","GRCh38-TPM1","GRCh38-TTN"  
  # "GRCh38-LAMA1","GRCh38-LAMA2","GRCh38-LAMA3","GRCh38-LAMA4","GRCh38-LAMA5",
  # "GRCh38-LAMB1","GRCh38-LAMB2","GRCh38-LAMB3",#"GRCh38-LAMB4",
  # "GRCh38-LAMC1","GRCh38-LAMC2","GRCh38-LAMC3"
  "Ssus11-LAMA1","Ssus11-LAMA3","Ssus11-LAMA4",#"Ssus11-LAMA2","Ssus11-LAMA5",
  "Ssus11-LAMB1","Ssus11-LAMB2","Ssus11-LAMB3",#"Ssus11-LAMB4",
  "Ssus11-LAMC1","Ssus11-LAMC2","Ssus11-LAMC3"
)
#path_saveoutputs <- "~/CM_ST_ALL_Chronic/FINAL_RESULTS/Factors_14/Markers/"
path_saveoutputs <- "~/CM_ST_ALL_Chronic/FINAL_RESULTS/Laminins/"

for (IDplotvar in 1:length(markers)){
  name_file <- paste(markers[IDplotvar],"_12wk_rep2",sep="")
  test1 <- FeatureOverlay(se.section1, features = markers[IDplotvar],
                          sampleids = 6,
                          cols = c("darkblue", "cyan", "yellow", "red", "darkred"),
                          #pt.size = 1.5,
                          #pt.alpha = 0.5,
                          #ncols.samples = 1,
                          dark.theme = T,
                          #min.cutoff='q1',
                          #ncol=2
  )
  ggplot2::ggsave(filename=file.path(path_saveoutputs,paste0(name_file,".png")),test1, height=6,width=6)
  ggplot2::ggsave(filename=file.path(path_saveoutputs,paste0(name_file,".eps")),test1, height=6,width=6)
  
}
############################ Quantification ###########################
se.section3_meta <- as.data.frame(se.section1@meta.data)
se.section3_meta$Barcode <- row.names(se.section3_meta)
cells1_eng <- se.section3_meta[se.section3_meta$SCT_snn_res.0.4==c(12) ,] ##12-1wk, 3-4wk, 7,8-12wk
data <- as.data.frame(se.section1@assays$RNA@data)
cells1_Eng_data <- data[, cells1_eng$Barcode]
cells1_Eng_data$gene <- row.names(cells1_Eng_data)
# write.table(cells1_Eng_data, "cells_engrafted_liy.tsv", sep="\t")
markers <- as.data.frame(c("GRCh38-MYH7", "Ssus11-MYH7"))

# markers <- as.data.frame(c("GRCh38-ACTC1","GRCh38-ACTN2","GRCh38-NKX2-5","GRCh38-MYL2","GRCh38-TNNT2","GRCh38-TPM1", "Ssus11-ACTC1","Ssus11-ACTN2","Ssus11-NKX2-5","Ssus11-MYL2","Ssus11-TNNT2","Ssus11-TPM1"
# ))

colnames(markers) <- "gene"
cells_MarkersExp <- markers %>% left_join(cells1_Eng_data, by="gene")
row.names(cells_MarkersExp) <- cells_MarkersExp$gene
#cells_MarkersExp <- cells_MarkersExp[, 1:190] 
#rowMeans(cells_MarkersExp)
apply(cells_MarkersExp, 1, median)
engrafted <- cells_MarkersExp
engrafted$gene <- row.names(engrafted)
engrafted_df <- gather(engrafted, key="cell", value = "Exp", -gene)
engrafted_df <- engrafted_df %>% arrange(gene)
engrafted_df$log2 <-log2(engrafted_df$Exp+1)
library(ggpubr)

my_cols <- c('#ea5353','#E69F00')
engrafted_df2 <- engrafted_df %>% separate(gene,sep="-",into=c("Species","Gene")) %>% dplyr::select(-c(cell,log2))
ggbarplot(engrafted_df2, x="Gene", y="Exp", xlab = NULL, ylab = "Normalized Expression",
          add=c("mean_se"), color ="Species", palette = my_cols,
          fill="Species", 
          position = position_dodge(0.72), error.plot = "upper_errorbar")

engrafted_df2 %>% group_by(Gene, Species) %>% 
  summarise(value=list(Exp)) %>% spread(Species, value) %>% 
  group_by(Gene) %>% 
  mutate(p_value = t.test(unlist(GRCh38), unlist(Ssus11))$p.value, paired = FALSE, t_value = t.test(unlist(GRCh38), unlist(Ssus11))$statistic)
############################ Functional Analysis for human cluster ###########################
setwd("~/Documents/nBox/CM_ST_ALL_Chronic/FINAL_RESULTS/Factors_14/FunctionalAna")
human_clus <- de.markers_se.section1 %>% dplyr::filter(cluster %in% c(7,8)) %>%
  filter(grepl("GRCh38-", gene))
human_clus$GeneSym <- sapply(strsplit(row.names(human_clus), "-"), "[[", 2)
human_clus$avg_log2FC <- as.numeric(human_clus$avg_log2FC)
pos <- human_clus[human_clus$avg_log2FC >= 0 ,]
#neg <- human_clus[human_clus$avg_log2FC < 0 ,]
pos$GeneSym2 <- sapply(strsplit((pos$gene), "-"), "[[", 2)
#neg$GeneSym2 <- sapply(strsplit((neg$gene), "-"), "[[", 2)
clust1_human_id <- bitr(pos$GeneSym2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
clust1_human_kegg <- enrichKEGG(gene = clust1_human_id$ENTREZID, organism     = 'hsa', pvalueCutoff = 0.05)
clust1_human_kegg  <- setReadable(clust1_human_kegg , OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
dotplot(clust1_human_kegg, showCategory = 25)
write.table(clust1_human_kegg, "Chronic_Human_12wk_kegg.txt", sep="\t")
#======================================. GSEA ==========================================
common_human_fc <- pos %>%  dplyr::select(GeneSym2, avg_log2FC) %>% arrange(-avg_log2FC)
common_human_fc_id <- bitr(common_human_fc$GeneSym2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
common_human_fc_id <- common_human_fc_id %>% left_join(common_human_fc, by=c("SYMBOL"="GeneSym2" ))
geneList = common_human_fc_id[,3]
names(geneList) = as.character(common_human_fc_id[,2])
geneList = sort(geneList, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
common_human_gseKEGG <- setReadable(kk2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
#dotplot(common_human_gseKEGG, showCategory=6, font.size=12)
common_human_gseKEGG  <- as.data.frame(common_human_gseKEGG@result)
write.table(common_human_gseKEGG, "Chronic_Human_12wk_KeggGSEA.txt", sep="\t")

#### plotting
pathdf <- common_human_gseKEGG[,c("Description","NES", "setSize")]
pathdf <- pathdf %>% arrange(NES)
pathdf$Description <- factor(pathdf$Description, levels=pathdf$Description)
pathdf$xlabel <- 'NES'
ggplot(pathdf,aes(x=xlabel, y=Description, fill=NES, label=setSize)) +
  geom_tile() + scale_fill_gradient(low="#deebf7", high = "#2171b5") +
  geom_text() + ggtitle("Chronic Human Genes Enriched Pathways")+ theme_classic()  +
  theme(axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold", color="black",   size=12),axis.line = element_blank())


VlnPlot(se.section1, features =c("GRCh38-ACTC1","GRCh38-ACTN2","GRCh38-NKX2-5","GRCh38-MYL2","GRCh38-TNNT2","GRCh38-TPM1", "Ssus11-ACTC1","Ssus11-ACTN2","Ssus11-NKX2-5","Ssus11-MYL2","Ssus11-TNNT2","Ssus11-TPM1"), group.by = "slide_id2")
 VlnPlot(se.section1, features =c("GRCh38-ACTC1","GRCh38-ACTN2","GRCh38-NKX2-5","GRCh38-MYL2","GRCh38-TNNT2","GRCh38-TPM1"), group.by = "slide_id2")
############################ SAVE/READ RDS ############################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
save(list = ls(), file = "Chronic_All_Analysis.Rdata")
#load("~/Documents/nBox/CM_ST_ALL_Acute/RESULTS_FINAL/Acute_All_Analysis.Rdata")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



