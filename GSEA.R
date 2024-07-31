rm(list = ls())
source("R/downloadTCGA.R")

exp <- getrna("TCGA-BRCA",type = 'tpm')
colnames(exp) <- substr(colnames(exp),1,16)
exp <- exp[,!duplicated(colnames(exp))] 

primarysample <- colnames(exp)[as.numeric(substr(colnames(exp),14,15)) == 1]

expprimary <- exp[,colnames(exp)[as.numeric(substr(colnames(exp),14,15)) == 1]] 

expprimary <- log2(expprimary + 1)

highFBXO2 <- colnames(expprimary[,expprimary['FBXO2',] > as.vector(quantile(as.numeric(expprimary['FBXO2',]))[4])])
lowFBXO2 <- colnames(expprimary[,expprimary['FBXO2',] <  as.vector(quantile(as.numeric(expprimary['FBXO2',]))[2])])

expprimary <- expprimary[rownames(expprimary) != "FBXO2",] 

expdiff <- expprimary[,c(highFBXO2,lowFBXO2)] 

group <- factor(c(rep("highFBXO2",length(highFBXO2)),rep('lowFBXO2',length(lowFBXO2))),
                levels = c("highFBXO2",'lowFBXO2'))

source("R/all_diff.R")

diff_res2 <- diff_analysis(expdiff,is_count = F,group = group)

library(org.Hs.eg.db)
library(clusterProfiler)

map <- bitr(diff_res2$deg_limma$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)

diffexp <- merge(diff_res2$deg_limma,map,by.x = "genesymbol",by.y = "SYMBOL")

geneList <- diffexp$logFC
names(geneList) <- diffexp$ENTREZID
geneList <- sort(geneList, decreasing = T)
geneList <- geneList[geneList != 0]


reskegg <- gseKEGG(geneList, organism = "hsa", 
                   pvalueCutoff = 0.05)

keggexp <- as.data.frame(reskegg@result) %>%
  dplyr::filter(!grepl("disease",Description)) %>%
  dplyr::filter(!grepl("lupus erythematosus",Description)) %>%
  dplyr::filter(!grepl("lateral sclerosis",Description)) %>%
  mutate(Description = factor(Description))


p2 <- ggplot(data = keggexp,aes(x = -log10(p.adjust), y = reorder(Description,-log10(p.adjust)),
                          fill = Description)) +
  geom_bar(width = 0.8, stat = 'identity') +
  scale_x_continuous(expand = c(0,0)) +
  ggsci::scale_fill_igv(alpha = 0.7) +
  labs(x = "-Log10 (adjust P value)",y=NULL) +
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  geom_text(data = keggexp,
            aes(x = 0.1, 
                y = Description,
                label = Description),
            size = 4,
            hjust = 0) +
  theme(legend.position = "none")


resgo <- gseGO(geneList, 'org.Hs.eg.db', 
               keyType = "ENTREZID", ont="ALL", pvalueCutoff=0.05) 

goexp <- as.data.frame(resgo@result) %>%
  dplyr::filter(!grepl("disease",Description))

#########
tmp <- split(goexp,list(resgo$ONTOLOGY))

tmp <- lapply(tmp, function(i){
  ep <- i
  i <- i[order(i$p.adjust, decreasing = FALSE),] %>%
    .[1:10,]
})
tmp$BP$Description[5] <- "adaptive immune response based on somatic recombination of immune receptors built from ISD"



p <- lapply(tmp, function(i){
  tmp1 <- ggplot(i,aes(y = reorder(Description,-log10(p.adjust)), x = -log10(p.adjust),fill = Description)) +
    geom_bar(width = 0.8, stat = 'identity') +
    scale_x_continuous(expand = c(0,0)) +
    ggsci::scale_fill_igv(alpha = 0.7) +
    labs(x = "-Log10 (adjust P value)") +
    theme_classic() +
    theme(axis.text.y = element_blank()) +
    geom_text(data = i,
              aes(x = 0.1, 
                  y = Description,
                  label = Description),
              size = 4.5,
              hjust = 0) +
    theme(legend.position = "none",
          axis.title.y = element_blank())
})

cowplot::plot_grid(p$BP+labs(title = "BP"),p$CC+labs(title = "CC"),p$MF+labs(title = "MF"),ncol = 3)

Figure5 <- cowplot::plot_grid(p2+labs(title = "KEGG"),p$BP+labs(title = "BP"),p$CC+labs(title = "CC"),p$MF+labs(title = "MF"),
                   ncol = 1,align = 'hv',rel_heights = c(1,0.8,0.8,0.8))

save_plot(filename = "R/组蛋白中介/image/Figure5.pdf",Figure5,base_height = 12,base_width = 8)

save_plot(filename = "R/组蛋白中介/image/Figure5A.pdf",p2+labs(title = "KEGG"),base_height = 4,base_width = 8)
save_plot(filename = "R/组蛋白中介/image/Figure5B.pdf",p$BP+labs(title = "BP"),base_height = 4,base_width = 8)
save_plot(filename = "R/组蛋白中介/image/Figure5C.pdf",p$CC+labs(title = "CC"),base_height = 4,base_width = 8)
save_plot(filename = "R/组蛋白中介/image/Figure5D.pdf",p$MF+labs(title = "MF"),base_height = 4,base_width = 8)

        
        