library(tidyverse)
library(ggpubr)
###load chronic all analysis RData and opn the r script
load("Chronic_All_Analysis.Rdata")
se.section4_meta <- as.data.frame(se.section1@meta.data)
se.section4_meta$Barcode <- row.names(se.section4_meta)

data2 <- as.data.frame(se.section1@assays$RNA@data)
genelist <- c("Ssus11-COL1A1","Ssus11-COL1A2","Ssus11-FN1","Ssus11-LUM")

#genelist <- c("GRCh38-MDK")
#genelist <- grep(pattern = "^GRCh38-CACN|^GRCh38-KCN", x = rownames(se.section1), value = TRUE)

meta_tp <- se.section4_meta[se.section4_meta$SCT_snn_res.0.4==c(12,3,7,8 ), ]
meta_tp$timepoint <- 'NA'
meta_tp$timepoint[which(meta_tp$SCT_snn_res.0.4%in%c(12))] <- '1wk'
meta_tp$timepoint[which(meta_tp$SCT_snn_res.0.4%in%c(3))] <- '4wk'
meta_tp$timepoint[which(meta_tp$SCT_snn_res.0.4%in%c(7,8))] <- '12wk'
meta_tp <- meta_tp %>% select(Barcode, timepoint)
sel_genes <- c(genelist)
exp_data_ori <- data2[row.names(data2)%in%sel_genes, meta_tp$Barcode]
exp_data_ori$gene <- row.names(exp_data_ori)
exp_data <- gather(exp_data_ori, key="cell", value="Exp", -gene)
exp_data <- exp_data %>% inner_join(meta_tp, by = c("cell"="Barcode")) %>% arrange(gene)
exp_data$log2 <-log2(exp_data$Exp+1)
exp_data <- exp_data %>% separate(gene,sep="-",into=c("Species","Gene")) %>% dplyr::select(-c(cell,log2))
exp_data$timepoint <- factor(exp_data$timepoint, levels = c("1wk","4wk","12wk"))
exp_data$Species <- factor(exp_data$Species, levels=c("Ssus11", "GRCh38"))


#exp_data %>% group_by(Species,timepoint) %>% summarise(value=list(log2)) %>% spread(Species, value) %>% group_by(timepoint) %>% mutate(p_value = t.test(unlist(GRCh38), unlist(Ssus11))$p.value)
my_cols <- c('#ea5353','#E69F00', '#7C9D8E')
head(exp_data)
# exp_data$Gene <- factor(exp_data$Gene, levels = c("MDK","VEGFA",
#                                                   "VEGFB","VEGFC"))
# 
# exp_data$Gene <- factor(exp_data$Gene, levels = c("MDK","VEGFA",
 #                                                 "VEGFB","VEGFC"))
exp_datal <- exp_data %>% filter(Exp>0) %>% mutate(Exp_l = log1p(Exp))
ggbarplot(exp_datal,
          x="Gene", y="Exp_l", xlab = NULL,
          ylab = "Normalized Expression", add=c("jitter","mean_sd"),
          add.params = list(size = 0.5),
          color ="timepoint", palette = my_cols,size=0.5,
          position = position_dodge(0.72), error.plot = "upper_errorbar") +
  #geom_jitter(width = 0.5, height = 0.5)+
  font("xy.text", size = 9) +
  font("xlab", size = 10)+
  font("ylab", size = 10)

subsetted <- subset(se.section1, subset = SCT_snn_res.0.4==c(12,3,7))
VlnPlot(subsetted, features = c("GRCh38-ACTN2", "Ssus11-ACTN2"))
# ggpubr::ggbarplot(exp_datal %>% filter(timepoint=='1wk'),
#                   x="Gene", y="Exp_l", color = "Species", 
#                   #fill="Species", 
#                   ylab = "Normalized Expression",
#                   palette = c('#E69F00','#ea5353'),
#                   add=c("jitter","mean_sd"),
#                   position = position_dodge(0.9),
#                   error.plot = "upper_errorbar")


cmi_markers_pval <- as.data.frame(exp_data %>% group_by(Gene, timepoint) %>% 
  summarise(value=list(Exp)) %>% spread(timepoint, value) %>% 
  group_by(Gene) %>% 
  mutate(wk1_vs_wk4 = t.test(unlist(`1wk`), unlist(`4wk`))$p.value, paired = FALSE, 
         wk4_vs_wk12 = t.test(unlist(`4wk`), unlist(`12wk`))$p.value, paired = FALSE))

pvals <- select(cmi_markers_pval, c('Gene','wk1_vs_wk4','wk4_vs_wk12'))
write.table(pvals, "CMI_markers_pvalues.txt", sep="\t")

exp_data %>% 
  group_by(Gene,timepoint) %>% 
  summarize(value=list(Exp)) %>% 
  spread(timepoint, value) %>% 
  mutate(p_value_1v4 = t.test(unlist(`1wk`), unlist(`4wk`))$p.value, 
         p_value_4v12 = t.test(unlist(`4wk`), unlist(`12wk`))$p.value, 
         ) 
