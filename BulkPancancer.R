rm(list = ls())
######

TCGAproject <- TCGAbiolinks::getGDCprojects()$project_id %>%
  .[grepl("TCGA-",.)] 

source("R/downloadTCGA.R")
rna_list <- lapply(TCGAproject,function(i) getrna(i,type = "tpm"))
names(rna_list) <- TCGAproject

FBXO2_list <- lapply(rna_list,function(i){
  exp <- i
  colnames(exp) <- substr(colnames(exp),1,16)
  exp <- exp[,!duplicated(colnames(exp))]
  exp <- log2(exp + 1)
  exp <- exp["FBXO2",]
  return(exp)
})
names(FBXO2_list) <- TCGAproject

FBXO2_tcga <- do.call(cbind,FBXO2_list) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(project = str_split(sample_id,pattern = ".TCGA",simplify = T)[,1],
         sample_id = str_split(sample_id,pattern = "\\.",simplify = T)[,2]) %>%
  mutate(type = ifelse(as.numeric(substr(sample_id,14,15)) > 9,"Paracancer","Cancer"),
         cancer = gsub("TCGA-","",project)) %>%
  mutate(type = factor(type,level =c("Paracancer","Cancer")))
  #mutate(type = ifelse(as.numeric(substr(sample_id,14,15)) > 9,"normal",
  #                    ifelse(as.numeric(substr(sample_id,14,15)) < 10,"tumor","other"))) %>%
  #dplyr::filter(type != "other")

p2 <- ggplot2::ggplot(FBXO2_tcga,aes(x= cancer,y=FBXO2,color=type,fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(3,8)) +
  #introdataviz::geom_split_violin(alpha = .4, trim = FALSE,width = 2) +
  labs(y="FBXO2 mRNA expression level") +
  scale_color_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[6])) +
  scale_fill_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[6])) +
  ggpubr::stat_compare_means(aes(group=type),label = 'p.signif',method = "wilcox.test",
                             symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, Inf), 
                                                symbols = c("***", "**", "*", " "))) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1)) 


p2
#rm(rna_list)

tcgapairexp <- lapply(TCGAproject, function(i) getpairrna(i,genetype = 'mRNA',type = "tpm"))
names(tcgapairexp) <- gsub("TCGA-","",TCGAproject)
tcgapairexp <- tcgapairexp[unlist(lapply(names(tcgapairexp), function(i){!is.null(tcgapairexp[[i]])}))]

FBXO2_pair <- tcgapairexp[[1]] 
FBXO2_pair <- FBXO2_pair["FBXO2",] %>%
  t() %>% as.data.frame() %>%
  mutate(sample_id = rownames(.))

normalsample <- rownames(FBXO2_pair)[as.numeric(substr(rownames(FBXO2_pair),14,15)) > 9]
tumorsample <- rownames(FBXO2_pair)[as.numeric(substr(rownames(FBXO2_pair),14,15)) < 10]
index <- match(substr(normalsample,1,12),substr(tumorsample,1,12))
FBXO2_pair <- FBXO2_pair[c(normalsample[index],tumorsample[index]),]

FBXO2_pair <- lapply(tcgapairexp,function(i){
  exp <- i
  exp <- exp["FBXO2",]
  return(exp)
})

FBXO2_pair <- do.call(cbind,FBXO2_pair) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(FBXO2 = log2(FBXO2 + 1),
         project = str_split(sample_id,pattern = ".TCGA",simplify = T)[,1],
         sample_id = str_split(sample_id,pattern = "\\.",simplify = T)[,2]) %>%
  mutate(type = ifelse(as.numeric(substr(sample_id,14,15)) > 9,"Paracancer","Cancer"),
         sample = as.factor(substr(sample_id,1,12))) %>%
  mutate(type = factor(type,level =c("Paracancer","Cancer")))

df_p_val <- FBXO2_pair %>% group_by(project) %>%
  rstatix::wilcox_test(FBXO2  ~ type,paired = T) %>%
  #rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p",cutpoints = c(0,0.001, 0.01, 0.05,Inf),
                            symbols = c("***", "**", "*","ns")) %>%
  rstatix::add_x_position(x = 'type')

p3 <- ggplot(FBXO2_pair,aes(x = type,y=FBXO2)) +
  geom_point(aes(fill = type,group = sample),size = 0.8,color = "grey70") +
  geom_boxplot(aes(fill = type,color = type),alpha = 0.7) +
  facet_grid(~project,scales = "free_x",switch = "x") +
  labs(x=NULL,y="FBXO2 mRNA expression level") +
  scale_y_continuous(expand = c(0,0),limits = c(0,9)) +
  geom_line(aes(group = sample),size = 0.5,color = "grey70",alpha=0.7) +
  scale_color_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[6])) +
  scale_fill_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[6])) +
  ggpubr::stat_pvalue_manual(df_p_val,y.position = 8.5,label = "p.signif",hjust = 1.2,
                             label.size = 6,hide.ns = T,remove.bracket = T) +
  theme_grey() +
  theme(axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.text.x = element_blank()) 

p3

library(cowplot)

Figure1 <- plot_grid(p2,p3,ncol = 1,labels = c("A","B"))

save_plot("R/组蛋白中介/image/Figure1.pdf",Figure1,base_height = 12,base_width = 18)



####甲基化
library(ChAMPdata)
data(probe.features)
all_probe <-  probe.features[probe.features$gene!='',] %>%
  tibble::rownames_to_column("probe")%>%
  dplyr::select(probe,feature,gene)%>%
  dplyr::arrange(gene) %>%
  dplyr::filter(gene == "FBXO2") 


load("ProcessData/TCGA-BRCA/methy.rda")

meth <- SummarizedExperiment::assay(data) %>%
  na.omit() %>%
  .[rownames(.) %in% all_probe$probe,] %>%
  t(.) %>%
  as.data.frame() %>%
  mutate(mean = rowMeans(.),
         sample_id = substr(rownames(.),1,16),
         type = ifelse(as.numeric(substr(sample_id,14,15)) > 9,"Paracancer","Cancer")) %>%
  mutate(type = factor(type,level =c("Paracancer","Cancer")))


p4 <- ggplot(meth,aes(x=type,y=mean,fill=type,color=type)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar",width=0.3) +
  scale_color_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[2])) +
  scale_fill_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[2])) +
  theme_bw(base_size = 14) +
  labs(y="Promoterr methylation level of FBXO2",x=NULL) +
  theme(legend.position = c(0.1,0.1),legend.title = element_blank(),
        #axis.text.x = element_text(angle = 45,hjust = 1)
        ) +
  ggpubr::stat_compare_means(label = "p.signif",label.x.npc = 'center',
                             comparisons = list(c("Paracancer","Cancer")),
                             hide.ns = F,color = 'red',
                             symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, Inf), 
                                                symbols = c("***", "**", "*", " ")),
                             label.y = 0.68,size=5)+
  theme(legend.position = "none") 

p4
############
FBXO2_BRCA <- FBXO2_tcga[FBXO2_tcga$cancer=="BRCA",] %>%
  dplyr::filter(type == "Cancer") 

BRCA_clinical <- getclinical('TCGA-BRCA')

BRCA_clinical <- inner_join(FBXO2_BRCA[,1:2],BRCA_clinical)

p5 <- ggplot(BRCA_clinical[BRCA_clinical$stage %in% c("I","II",'III','IV'),],aes(x=stage,y=FBXO2,fill=stage,color =stage)) + 
  geom_boxplot(outlier.shape = NA,width=0.5) +
  stat_boxplot(geom = "errorbar",width=0.05) +
  geom_violin(trim = F,alpha=0.5) +
  ggsci::scale_fill_aaas(alpha = 0.7) +
  ggsci::scale_color_aaas(alpha = 0.7) +
  ggpubr::stat_compare_means(label.y = 9.4) +
  theme_bw(base_size = 14) +
  labs(x=NULL,y="FBXO2 mRNA expression level") +
  theme(legend.position = "none")
p5

########蛋白
CPTAC <- data.table::fread("RAW/TCGA/TCGA-BRCA/CPTAC2_Breast_Proteome_tmt10.txt") %>% 
  as.data.frame() %>%
  column_to_rownames(var = "Gene") %>%
  .[-c(1:3),] %>%
  .[,grepl("Unshared.Log.Ratio",colnames(.))]

colnames(CPTAC) <- gsub(" Unshared Log Ratio","",colnames(CPTAC))
colnames(CPTAC) <- unlist(lapply(colnames(CPTAC),function(i){strsplit(i,split = "\\.")[[1]][1]}))

CPTAC <- CPTAC[,!duplicated(colnames(CPTAC))]

CPTAC_clinial <- read.csv("RAW/TCGA/TCGA-BRCA/PDC_study_biospecimen_10262023_143312.csv") %>%
  dplyr::filter(`Case.Submitter.ID` !="Internal Reference - Pooled Sample") %>%
  dplyr::filter(`Sample.Type` != "Not Reported")


CPTAC <- CPTAC[CPTAC_clinial$Aliquot.Submitter.ID] 
CPTAC <- impute::impute.knn(as.matrix(CPTAC))

CPTAC <- CPTAC$data[rownames(CPTAC$data) == "FBXO2",] %>%
  as.data.frame() %>%
  mutate(class = CPTAC_clinial$Sample.Type) %>%
  mutate(class = ifelse(class == "Solid Tissue Normal","Normal","Tumor"))

CPTAC_melt <- data.table::melt(CPTAC) 

head(CPTAC_melt)

p6 <- ggplot(CPTAC_melt,aes(x= class,y=value,fill=class,color=class)) +
  geom_boxplot(outlier.color = "white",width = 0.5) +
  stat_boxplot(geom="errorbar",width=0.08) +
  geom_violin(trim = F,alpha=0.3) +
  scale_y_continuous(limits = c(-2,4)) +
  scale_color_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[2])) +
  scale_fill_manual(values = c(ggsci::pal_aaas(alpha = 0.7)(6)[1],ggsci::pal_aaas(alpha = 0.7)(7)[2])) +
  theme_bw(base_size = 14) +
  labs(y = "Protein Expression of FBXO2 (Z value)",x=NULL) +
  theme(legend.position = "none",
       # axis.text.x = element_text(angle = 40,hjust = 1)
       ) +
  ggpubr::stat_compare_means(comparisons = list(c("Normal","Tumor")),
                             #method = "t.test",
                             label = "p.signif",
                             symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, Inf), 
                                                symbols = c("***", "**", "*", " ")),
                            color ='red',
                             hide.ns = F,label.y = 3.4,size=5)

library(cowplot)

Figure3 <- plot_grid(p5,p4,p6,ncol = 3,labels = c("A","B","C"),
                     rel_widths = c(1,0.5,0.5),label_size = 18,align = "hv")
Figure3 
save_plot("R/组蛋白中介/image/Figure3.pdf",Figure3,base_height = 4.5,base_width = 8.5)





