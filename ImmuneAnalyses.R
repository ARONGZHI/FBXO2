rm(list = ls())
library(IOBR)
library(tidyverse)

exp_list <- lapply(rna_list, function(i){
  exp <- i
  colnames(exp) <- substr(colnames(exp),1,16)
  exp <- exp[,!duplicated(colnames(exp))]
  #normalsample <- colnames(exp)[as.numeric(substr(colnames(exp),14,15)) > 9]
  #tumorsample <- colnames(exp)[as.numeric(substr(colnames(exp),14,15)) == 1]
  #exp <- exp[,c(normalsample,primmalysample)]
  exp <- exp[rowMeans(exp)>0,]
  exp <- log2(exp + 1)
  return(exp)
})


cibersort_list <- lapply(exp_list, function(i){
  cibersort <- deconvo_tme(eset = i,method = "cibersort",arrays = FALSE,perm = 1000)
  return(cibersort)
})

##load("ProcessData/TCGAPancan/cibersort.rda")

cor_res <- lapply(TCGAproject, function(i){
  
  immune_exp <- FBXO2_list[[i]] %>% t() %>%
    as.data.frame() %>% rownames_to_column(var = "ID") %>%
    inner_join(.,rownames_to_column(exp_list[[i]][signature_collection$Immune_Checkpoint,] %>% t() %>% as.data.frame(),var="ID")) %>%
    inner_join(.,cibersort_list[[i]][,1:23])
  
  
  cor_res <- lapply(colnames(immune_exp)[3:ncol(immune_exp)], function(x){
    corT <- cor.test(immune_exp[,"FBXO2"],immune_exp[,x],method="pearson")
    cor <- corT$estimate
    pValue <- corT$p.value
    res <- data.frame(cor = cor,pValue = pValue)
  })
  names(cor_res) <- colnames(immune_exp)[3:ncol(immune_exp)]
  cor_res <- do.call(rbind,cor_res) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'cell')
  
  cor_res$project <-  i
  
  return(cor_res)
})

corimmune <- do.call(rbind,cor_res) %>%
  mutate(psig = ifelse(.$pValue < 0.05,
                       ifelse(.$pValue < 0.01,
                              ifelse(.$pValue < 0.001,"***", "**"), "*"), ""),
         Cancer = gsub("TCGA-","",project),
         cor = round(cor,2)) %>%
  mutate(celltype = gsub("_CIBERSORT","",cell)) %>%
  mutate(celltype = gsub("_CIBERSORT","",celltype)) %>%
  mutate(type = ifelse(celltype %in% signature_collection$Immune_Checkpoint,"Immune Checkpoint","CIBERSORT"))

corimmune[,c(6,7,8)] <- lapply(corimmune[,c(6,7,8)],as.factor)
corimmune$type <- factor(corimmune$type,levels = c("Immune Checkpoint","CIBERSORT"))

p3 <- ggplot(corimmune,aes(x=Cancer,y=celltype)) +
  facet_wrap(~type,scales = 'free_y',ncol = 1) +
  geom_tile(aes(fill = cor), colour = "grey",size=0.1) +
  scale_fill_gradient2(low = "#3B4992FF",mid = "white",high = "#EE0000FF") +
  geom_text(aes(label=psig),col ="black",size = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5)) +
  labs(y=NULL,x=NULL)

save_plot("R/组蛋白中介/image/Figure6.pdf",p3,base_height = 8,base_width = 11)
