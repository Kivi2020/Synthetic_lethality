TCGA.BRCA.MUT <- read.maf('C:/Users/kirthana/Downloads/TCGA.BRCA.muse.cd0d1636-ae95-4141-94b3-8218eb3b8b25.DR-10.0.protected.maf')
TCGA.STAGE1 <- read.tsv(denseDataOnlyDownload)

TCGA.BRCA.STAGE <- cbind(denseDataOnlyDownload[[1]], denseDataOnlyDownload[[3]])
TCGA.BRCA.STAGE <- data.frame(TCGA.BRCA.STAGE)

TCGA.BRCA.DATA1 <- data.frame(TCGA.BRCA.MUT@data)
View(TCGA.BRCA.DATA1)

library(dplyr)
library(tidyr)
Tumor_sample_barcode <- data.frame(TCGA.BRCA.DATA1$Tumor_Sample_Barcode)
y <- Tumor_sample_barcode %>% separate(TCGA.BRCA.DATA1.Tumor_Sample_Barcode, c('Tumor', 'Sample', 'barcode', 'patient_number'))
cols <- c('Tumor', 'Sample', 'barcode', 'patient_number')
Tumor_sample_barcode <- apply( y[ , cols ] , 1 , paste , collapse = "-" )

TCGA.BRCA.DATA1 <- data.frame(Tumor_sample_barcode, TCGA.BRCA.DATA1)
View(TCGA.BRCA.DATA1)

TCGA.STAGE2 <- TCGA.BRCA.STAGE[(match(TCGA.BRCA.DATA1[ ,1], TCGA.BRCA.STAGE[ ,1] )), ]
View(TCGA.STAGE2)

TCGA.MUT.STAGE <- data.frame(TCGA.STAGE2, TCGA.BRCA.DATA1)

TCGA.MUT.STAGE <- TCGA.MUT.STAGE %>% select(c(-Tumor_sample_barcode))
View(TCGA.MUT.STAGE)

TCGA.MUT.STAGE.maf <- read.maf(TCGA.MUT.STAGE)
final_BRCA <- somaticInteractions(TCGA.MUT.STAGE.maf, pvalue = c(0.01))
View(final_BRCA)

cv <- data.frame(TCGA.MUT.STAGE$X2, TCGA.MUT.STAGE)
stage1 <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage i"))
stage1a <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage ia"))
stage1b <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage ib"))
stage4 <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iv"))
stage2 <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage ii"))
stage2a <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iia"))
stage2b <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iib"))
stage3 <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iii"))
stage3a <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iiia"))
stage3b <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iiib"))
stage3c <- cbind(subset(cv, cv$TCGA.MUT.STAGE.X2 == "stage iiic"))


stage1.maf <- read.maf(stage1)
stage1.maf.SI <- somaticInteractions(stage1.maf, pvalue = c(0.01))
stage1a.maf <- read.maf(stage1a)
stage1a.maf.SI <- somaticInteractions(stage1a.maf, pvalue = c(0.01))
stage1b.maf <- read.maf(stage1b)
stage1b.maf.SI <- somaticInteractions(stage1b.maf, pvalue = c(0.01))

stage4.maf <- read.maf(stage4)
stage4.maf.SI <- somaticInteractions(stage4.maf, pvalue = c(0.01))

stage1.maf.1 <- data.frame(stage1.maf.SI$Event, stage1.maf.SI)
stage1_BRCA_co <- cbind(subset(stage1.maf.1, stage1.maf.1$stage1.maf.SI.Event == "Co_Occurence"))
stage1_BRCA_ME <- cbind(subset(stage1.maf.1, stage1.maf.1$stage1.maf.SI.Event == "Mutually_Exclusive"))

stage1a.maf.1 <- data.frame(stage1a.maf.SI$Event, stage1a.maf.SI)
stage1a_BRCA_co <- cbind(subset(stage1a.maf.1, stage1a.maf.1$stage1a.maf.SI.Event == "Co_Occurence"))
stage1a_BRCA_ME <- cbind(subset(stage1a.maf.1, stage1a.maf.1$stage1a.maf.SI.Event == "Mutually_Exclusive"))

stage1b.maf.1 <- data.frame(stage1b.maf.SI$Event, stage1b.maf.SI)
stage1b_BRCA_co <- cbind(subset(stage1b.maf.1, stage1b.maf.1$stage1b.maf.SI.Event == "Co_Occurence"))
stage1b_BRCA_ME <- cbind(subset(stage1b.maf.1, stage1b.maf.1$stage1b.maf.SI.Event == "Mutually_Exclusive"))

stage4.maf.1 <- data.frame(stage4.maf.SI$Event, stage4.maf.SI)
stage4_BRCA_co <- cbind(subset(stage4.maf.1, stage4.maf.1$stage4.maf.SI.Event == "Co_Occurence"))
stage4_BRCA_ME <- cbind(subset(stage4.maf.1, stage4.maf.1$stage4.maf.SI.Event == "Mutually_Exclusive"))

SL_to_SV_1_4 <- cbind(intersect(stage1_BRCA_ME$pair, stage4_BRCA_co$pair))
SL_to_SV_1a_4 <- cbind(intersect(stage1a_BRCA_ME$pair, stage4_BRCA_co$pair))
SL_to_SV_1b_4 <- cbind(intersect(stage1b_BRCA_ME$pair, stage4_BRCA_co$pair))

SL_to_SV_1_to_4 <- rbind(SL_to_SV_1_4, SL_to_SV_1a_4, SL_to_SV_1b_4)
View(SL_to_SV_1_to_4)
SL_to_SV_1_to_4 <- unique(SL_to_SV_1_to_4)






