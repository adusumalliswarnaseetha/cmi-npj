# SCT_snn_res.0.5
se.section1

RNA_data <- as.data.frame(se.section1@assays$RNA@data)

#clusno=0
for(clusno in unique(se.section1@meta.data$SCT_snn_res.0.5)){
  spot_ids <- se.section1@meta.data %>% filter(SCT_snn_res.0.5==clusno) %>% row.names()
  cat("cluster: ", clusno, "nSpots:", length(spot_ids),"\n")
  rna_data <- RNA_data[, spot_ids]
  res <- apply(rna_data, 2, function(xdf) { cal_percent(xdf) })
  number.df <- do.call("rbind", res)
  number.df$genes_diff = number.df$nHuman - number.df$nPig
  number.df$exp_diff = number.df$avgExpHuman - number.df$avgExpPig
  fname <- paste0("Rep2_D_wk_Cluster_counts_0.5_clusno", clusno,".tsv")
  write.table(number.df, fname, sep="\t")
}

cal_percent <- function(ldf){
  #ldf2 = ldf[ldf>0]
  names(ldf) <- row.names(rna_data)
  gnames <- names(ldf)
  gnames <- gnames[ldf>20]
  ldf2 <- ldf[ldf>20]
  human_var <- ldf2[grepl("^GRCh38", gnames)]
  pig_var <- ldf2[grepl("Ssus11", gnames)]
  nHuman <- sum(grepl("^GRCh38", gnames))
  nPig <- sum(grepl("Ssus11", gnames))
  avgExp_Human <- mean(human_var)
  avgExp_Pig <- mean(pig_var)
  human_genes <- paste(gnames[grepl("^GRCh38", gnames)],collapse = ",")
  pig_genes <- paste(gnames[grepl("Ssus11", gnames)],collapse = ",")
  #cat("total genes:", length(ldf), "fil genes:", sum(ldf>10), "nHuman:", nHuman, "nPig:", nPig ,"\n")
  tmpdf <- data.frame(nHuman=nHuman, nPig = nPig, avgExpHuman = avgExp_Human, avgExpPig = avgExp_Pig,
                      genesHuman = human_genes, genesPig=pig_genes)
  return(tmpdf)
}



#se.section1 <- se
# se.section1 <- SCTransform(se.section1)
# se. <- SCTransform(se)
# spots_human <- read.table("/Users/swarna/Documents/nBox/CM_Spatial_2_12wk/Acute_2w/results/TEST_Cluster_counts_rep1_0.5_clusno2.tsv", sep="\t")
# 
# human_top25 <- spots_human %>% arrange(desc(nHuman)) %>% head(200) %>% row.names
# 
# se.section1@meta.data$human_top <- 0
# se.section1@meta.data$human_top[row.names(se.section1@meta.data)%in%human_top25] <- 1
# 
# se.section1 <- LoadImages(se.section1, time.resolve = F, verbose = T)
# ImagePlot(se.section1, method = "raster", darken = TRUE, type = "raw")
# ST.FeaturePlot(object = se.section1, 
#                features = c("human_top"), dark.theme = F,
#                cols = c("violet", "red")
# )


FeatureOverlay(se.section1, 
               features = c("GRCh38-TNNC1"), 
               pt.size = 1.5,
               cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
               dark.theme = T, 
               type = "raw",
               sampleids = 1,min.cutoff="q1")

FeatureOverlay(se.section1, 
               features = c("Ssus11-TNNC1"), 
               pt.size = 1.5,
               cols = c("dark blue", "cyan", "yellow", "red", "dark red"), 
               dark.theme = T, 
               type = "raw",
               sampleids = 1,min.cutoff="q1")


# spot_info <- as.data.frame(t(apply(RNA_data, 2, function(x) { cal_proportion(x, row.names(RNA_data))})))
# colnames(spot_info) <- c("numGenes_Human", "numGenes_Monkey", "avgExp_Human", "avgExp_Monkey","topGenes_Human","topGenes_Monkey")
# spot_info$numGenes_Human <- as.numeric(spot_info$numGenes_Human)
# spot_info$numGenes_Monkey <- as.numeric(spot_info$numGenes_Monkey)
# spot_info$avgExp_Human <- as.numeric(spot_info$avgExp_Human)
# spot_info$avgExp_Monkey <- as.numeric(spot_info$avgExp_Monkey)
# 
# spot_info$diff_GeneCounts <- spot_info$numGenes_Human - spot_info$numGenes_Monkey
# spot_info$diff_avgExp <- spot_info$avgExp_Human - spot_info$avgExp_Monkey
# write.table(spot_info, "D1_combined_SpotInfo.txt", sep="\t")

##########################################
###### per cluster - spot info summarized

files <- list.files(pattern = "*.tsv")

process_file <- function(x){
  df = read.table(x, sep="\t", header=T)
  nonzero_human <- sum(df$nHuman>0)
  nonzero_pig <- sum(df$nPig>0)
  nHuman_avg <- mean(df$nHuman, na.rm=TRUE)
  nPig_avg <- mean(df$nPig, na.rm=TRUE)
  avgExpHuman_mean <- mean(df$avgExpHuman, na.rm=TRUE)
  avgExpPig_mean <- mean(df$avgExpPig, na.rm=TRUE)
  resdf <- data.frame(file=x, nonzero_human=nonzero_human, nonzero_pig=nonzero_pig, 
                      nHuman_avg=nHuman_avg, nPig_avg=nPig_avg,
                      avgExpHuman_mean=avgExpHuman_mean, avgExpPig_mean=avgExpPig_mean)
  return(resdf)
}

res_stats <- do.call("rbind", lapply(files, process_file))
