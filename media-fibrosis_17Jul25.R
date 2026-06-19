library(tidyverse)
library(ggpubr)
###load chronic all analysis RData and opn the r script
load("Chronic_All_Analysis.Rdata")

se.section3_meta <- as.data.frame(se.section1@meta.data)
se.section3_meta$Barcode <- row.names(se.section3_meta)
a_wk1_inf <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('4') ,] ##12-1wk, 3-4wk, 7,8-12wk
a_wk2_inf <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('12') ,] ##12-1wk, 3-4wk, 7,8-12wk
a_wk1_ni <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('15') ,] ##12-1wk, 3-4wk, 7,8-12wk
a_wk2_ni <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('11') ,] ##12-1wk, 3-4wk, 7,8-12wk


c_wk1_inf <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('2') ,] ##12-1wk, 3-4wk, 7,8-12wk
c_wk4_inf <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('1','13') ,] ##12-1wk, 3-4wk, 7,8-12wk
c_wk12_inf <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('13','14','8') ,] ##12-1wk, 3-4wk, 7,8-12wk

c_wk1_ni <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('0') ,] ##12-1wk, 3-4wk, 7,8-12wk
c_wk4_ni <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('10','5') ,] ##12-1wk, 3-4wk, 7,8-12wk
c_wk12_ni <- se.section3_meta[se.section3_meta$SCT_snn_res.0.6==c('3') ,] ##12-1wk, 3-4wk, 7,8-12wk

data <- as.data.frame(se.section1@assays$RNA@data)
############################################
feature <- c("Ssus11-COL1A1", "Ssus11-COL1A2","Ssus11-FN1","Ssus11-LUM")
markers <- as.data.frame(feature)

colnames(markers) <- "Gene"


############################################
a_wk1_inf_data <- data[, a_wk1_inf$Barcode]
a_wk2_inf_data <- data[, a_wk2_inf$Barcode]
a_wk1_ni_data <- data[, a_wk1_ni$Barcode]
a_wk2_ni_data <- data[, a_wk2_ni$Barcode]

c_wk1_inf_data <- data[, c_wk1_inf$Barcode]
c_wk4_inf_data <- data[, c_wk4_inf$Barcode]
c_wk12_inf_data <- data[, c_wk12_inf$Barcode]

c_wk1_ni_data <- data[, c_wk1_ni$Barcode]
c_wk4_ni_data <- data[, c_wk4_ni$Barcode]
c_wk12_ni_data <- data[, c_wk12_ni$Barcode]

colnames(a_wk1_inf_data) <- paste("Inf_a_wk1", colnames(a_wk1_inf_data), sep="_")
colnames(a_wk2_inf_data) <- paste("Inf_a_wk2", colnames(a_wk2_inf_data), sep="_")
colnames(a_wk1_ni_data) <- paste("NI_a_wk1", colnames(a_wk1_ni_data), sep="_")
colnames(a_wk2_ni_data) <- paste("NI_a_wk2", colnames(a_wk2_ni_data), sep="_")

colnames(c_wk1_inf_data) <- paste("Inf_c_wk1", colnames(c_wk1_inf_data), sep="_")
colnames(c_wk4_inf_data) <- paste("Inf_c_wk4", colnames(c_wk4_inf_data), sep="_")
colnames(c_wk12_inf_data) <- paste("Inf_c_wk12", colnames(c_wk12_inf_data), sep="_")
colnames(c_wk1_ni_data) <- paste("NI_c_wk1", colnames(c_wk1_ni_data), sep="_")
colnames(c_wk4_ni_data) <- paste("NI_c_wk4", colnames(c_wk4_ni_data), sep="_")
colnames(c_wk12_ni_data) <- paste("NI_c_wk12", colnames(c_wk12_ni_data), sep="_")

a_wk1_inf_data$gene <- row.names(a_wk1_inf_data)
a_wk2_inf_data$gene <- row.names(a_wk2_inf_data)
a_wk1_ni_data$gene <- row.names(a_wk1_ni_data)
a_wk2_ni_data$gene <- row.names(a_wk2_ni_data)
c_wk1_inf_data$gene <- row.names(c_wk1_inf_data)
c_wk4_inf_data$gene <- row.names(c_wk4_inf_data)
c_wk12_inf_data$gene <- row.names(c_wk12_inf_data)
c_wk1_ni_data$gene <- row.names(c_wk1_ni_data)
c_wk4_ni_data$gene <- row.names(c_wk4_ni_data)
c_wk12_ni_data$gene <- row.names(c_wk12_ni_data)

merged2 <- a_wk1_inf_data %>% 
  inner_join(a_wk2_inf_data,by='gene') %>%
  inner_join(a_wk1_ni_data, by = 'gene')%>%
  inner_join(a_wk2_ni_data, by = 'gene')%>%
  
  inner_join(c_wk1_inf_data, by = 'gene')%>%
  inner_join(c_wk4_inf_data, by = 'gene')%>%
  inner_join(c_wk12_inf_data, by = 'gene')%>%
  inner_join(c_wk1_ni_data, by = 'gene')%>%
  inner_join(c_wk4_ni_data, by = 'gene')%>%
  inner_join(c_wk12_ni_data, by = 'gene')

merged <- merged2[which(merged2$gene%in%markers$Gene),]
dim(merged[merged$gene%in%markers$Gene,])

merged_t <- merged %>% dplyr::select(-gene) %>% t()
colnames(merged_t) <- merged$gene
merged_t <- as.data.frame(merged_t)
head(merged_t)
merged_t$Group <- as.character(unlist(lapply(row.names(merged_t), function(x) { paste0(unlist(strsplit(x, "_"))[1:3],collapse = "_")})))
merged_t$Group <- as.character(unlist(lapply(row.names(merged_t), function(x) { paste0(unlist(strsplit(x, "_"))[1:3],collapse = "_")})))

merged_t <- merged_t %>% gather("Gene","Exp", -Group)
merged_t$Group <- factor(merged_t$Group, levels = c("Inf_a_wk1", "Inf_a_wk2", "NI_a_wk1", "NI_a_wk2",
                                                    "Inf_c_wk1","Inf_c_wk4", "Inf_c_wk12",
                                                    "NI_c_wk1",   "NI_c_wk4",   "NI_c_wk12" ))


###### PLOTS
library(ggpubr)
my_cols <- c('#ea5353','#E69F00', '#7C9D8E')
merged_t$Gene <- factor(merged_t$Gene, levels = c("Ssus11-COL1A1", "Ssus11-COL1A2", "Ssus11-FN1", "Ssus11-LUM"))

# merged_data <- merged_t %>% extract(Group, into=c("GraftType","Timepoint"), 
#                                     regex = "([[:alnum:]]+)_([[:alnum:]]+)")
# 
# merged_data$Timepoint <- factor(merged_data$Timepoint, levels = c("1wk", "4wk", "12wk"))
# merged_data$GraftType <- factor(merged_data$GraftType, levels = c("Inf","eng"))

#infracted_plot <- 
                  ggbarplot(merged_t %>% mutate(Exp_l = log1p(Exp)) %>% 
                              filter(Group!="NI_a_wk1" & Group!="NI_a_wk2" & Group!="NI_c_wk1" & Group!="NI_c_wk4"& Group!="NI_c_wk12" & Group!="Inf_a_wk1"& Group!="Inf_a_wk2" ),
                            x="Gene", y="Exp_l", xlab = NULL, #fill= "Group", 
                            ylab = "Normalized Expression", add=c("jitter","mean_sd"), color ="Group", palette = my_cols,
                            position = position_dodge(0.72), error.plot = "upper_errorbar", title = "Infarcted - pig genes")+font("xy.text", size = 9) +font("xlab", size = 10)+font("ylab", size = 10)

engrafted_plot <- ggbarplot(merged_data %>% filter(GraftType=="eng") %>% mutate(Exp_l = (Exp)),
                            x="Gene", y="Exp_l", xlab = NULL, fill= "Timepoint", 
                            ylab = "Normalized Expression", add=c("mean_se"), color ="Timepoint", palette = my_cols,
                            position = position_dodge(0.72), error.plot = "upper_errorbar", title = "Engrafted - pig genes")

#############################################
media_infarcted_pval <- merged_t %>% group_by(Gene, Group) %>%
  summarise(value = list(Exp)) %>% spread(Group, value) %>%
  group_by(Gene) %>%
  mutate(
    #Inf_acute_1vs2_t = t.test(unlist(`Inf_a_wk1`), unlist(`Inf_a_wk2`), paired = FALSE)$p.value,
    Inf_chronic_1vs4_t = t.test(unlist(`Inf_c_wk1`), unlist(`Inf_c_wk4`), paired = FALSE)$p.value,
    Inf_chronic_4vs12_t = t.test(unlist(`Inf_c_wk4`), unlist(`Inf_c_wk12`), paired = FALSE)$p.value ) %>% dplyr::select(contains("_t"))

write.table(media_infarcted_pval,'media_infarcted_pval.tsv', sep = "\t", col.names =T, row.names=F)

NUM_SIM = 1000
N_SAMPLE = 200

res_list = list()
for (i_sim in 1:NUM_SIM) {
  cat("working for iteration:", i_sim,"\n")
  merged_r = merged_t %>% group_by(Group) %>% sample_n(N_SAMPLE, replace = FALSE)
  test_res <- merged_t %>% group_by(Gene, Group) %>%
    summarise(value = list(Exp)) %>% spread(Group, value) %>%
    group_by(Gene) %>%
    mutate(
      #Inf_acute_1vs2_t = t.test(unlist(`Inf_a_wk1`), unlist(`Inf_a_wk2`), paired = FALSE)$p.value,
      Inf_chronic_1vs4_t = t.test(unlist(`Inf_c_wk1`), unlist(`Inf_c_wk4`), paired = FALSE)$p.value,
      Inf_chronic_4vs12_t = t.test(unlist(`Inf_c_wk4`), unlist(`Inf_c_wk12`), paired = FALSE)$p.value 
     )%>% dplyr::select(contains("_t"))
  
  test_res$iter <- i_sim
  res_list[[i_sim]] = test_res %>% gather(key = "test_comp", value = "pvalue",-Gene,-iter)
}

sim_res = do.call(rbind, res_list)
sim_res_summary2 <- sim_res %>% 
  arrange(iter, Gene) %>% 
  group_by(Gene, test_comp) %>%
  reframe(padj= p.adjust(pvalue, method="BH")) %>%
  group_by(Gene, test_comp) %>%
  summarise(n_sig = sum(.data[["padj"]]<0.05), n_fail = sum(.data[["padj"]]>=0.05))  

sim_res_summary2_table <- sim_res_summary2 %>% dplyr::select(-n_fail) %>% spread(test_comp, n_sig)
sim_res_summary2_table
############################
FeatureOverlay(se.section1, 
               features = c("Ssus11-TGFB1"), 
               pt.size = 1.8,
               cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
               dark.theme = F, 
               type = "raw",
               sampleids = c(3,4,5),
               ncols = 1,
               #min.cutoff = 'q1',
               show.sb=F
)

############################ SAVE/READ RDS ############################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
save(list = ls(), file = "Acute_chronic_Fibrosis_Analysis.Rdata")
