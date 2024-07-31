rm(list = ls())

library(TwoSampleMR)
library(tidyverse)

###########
exposure <- extract_instruments(outcomes = "eqtl-a-ENSG00000116661")

outcome <- extract_outcome_data(snps = exposure$SNP,outcomes = "ieu-a-1130",proxies = F)

mrdat13 <- harmonise_data(exposure,outcome,action = 2) %>%
  mutate(R2 = (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure))/
           (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure) + 
              2 * samplesize.exposure*eaf.exposure * (1 - eaf.exposure) * se.exposure^2)) %>%
  mutate(Fvalue = R2 * (samplesize.exposure - 2) / (1 - R2)) %>%
  dplyr::filter(Fvalue > 10) %>%
  dplyr::filter(pval.outcome >= 1e-5)

#load("R/组蛋白中介/mrdat13.rda")

mr_res13 <- mr(mrdat13)  %>%
  generate_odds_ratios()

pleiotropy13 <- mr_pleiotropy_test(mrdat13)
heterogeneity13 <- mr_heterogeneity(mrdat13)

singlesnp_res13 <- mr_singlesnp(mrdat13)


#########
exposure <- extract_instruments(outcomes = "ieu-a-1130")

outcome <- extract_outcome_data(snps = exposure$SNP,outcomes = "eqtl-a-ENSG00000116661",proxies = F)

mrdat31 <- harmonise_data(exposure,outcome,action = 2) %>%
  mutate(R2 = (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure))/
           (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure) + 
              2 * samplesize.exposure*eaf.exposure * (1 - eaf.exposure) * se.exposure^2)) %>%
  mutate(Fvalue = R2 * (samplesize.exposure - 2) / (1 - R2)) %>%
  dplyr::filter(Fvalue > 10) %>%
  dplyr::filter(pval.outcome >= 1e-5)

#load("R/组蛋白中介/mrdat31.rda")

mr_res31 <- mr(mrdat31)  %>%
  generate_odds_ratios()

pleiotropy31 <- mr_pleiotropy_test(mrdat31)
heterogeneity31 <- mr_heterogeneity(mrdat31)

singlesnp_res31 <- mr_singlesnp(mrdat31)

########

###################
library(forestploter)

forestdat <- data.frame(Exposure = c("FBXO2","","","Breast Cancer","",""),
                        Outcome = c("Breast Cancer","","","FBXO2","",''),
                        Method = c(mr_res13$method[1:3],mr_res31$method[1:3]),
                        Nsnp = c(4,"","",9,"",""),
                        OR = c(round(mr_res13$or[1:3],3),round(mr_res31$or[1:3],3)),
                        lor = c(round(mr_res13$or_lci95[1:3],3),round(mr_res31$or_lci95[1:3],3)),
                        uor = c(round(mr_res13$or_uci95[1:3],3),round(mr_res31$or_uci95[1:3],3))) %>%
  mutate( CI = paste0(OR," (",lor,
                      " - ",uor,")"),
          Index = paste(rep(" ",12),collapse = " "),
          #P = format(mr_res$pval,scientific=T,digits=3)
          P = c(round(mr_res13$pval[1:3],3),round(mr_res31$pval[1:3],3))) %>%
  mutate("P value" = ifelse(P < 0.05,paste0(P,"*"),P))

colnames(forestdat) <- c("Exposure","Outcome","Method","N variants",
                         "OR","lor","uor","OR (95%CI)"," ","P","P Value")

tm <- forest_theme(ci_col = "#008B45FF",refline_col = "#EE0000FF")

pdf(file = "R/组蛋白中介/image/FigureS1.pdf",width = 10,height = 5)
p1 <- forest(forestdat[,c(1,2,3,4,9,8,11)],est = forestdat$OR,lower = forestdat$lor,
       upper = forestdat$uor,ci_column = 5,xlim = c(0.7,1.3),
       ticks_at = c(0.7,1,1.3),ref_line = 1,theme = tm)

#dev.off()

edit_plot(p1,row = which(forestdat$P <0.05),
                gp = grid::gpar(col = "red", fontface = "bold"))
dev.off()
