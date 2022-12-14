
SL.data <- data.frame(All_synthetic_lethality_gene_pairs)
SL.data2 <- data.frame(negative_intr_biogrid)
SL.data1 <-subset(SL.data, SL.data$V1!="our prediction" & SL.data$V1!="SLDB" & SL.data$V1!="Michael Boettcher")

SL.cgidb <- cbind(SL.data1[[3]], SL.data1[[4]])
SL.biogrid <- cbind(SL.data2[[1]], SL.data2[[2]])

SL.cgidb <- data.frame(SL.cgidb)
SL.biogrid <- data.frame(SL.biogrid)

SLG1.c <- c(SL.cgidb[[1]])
SLG2.c <- c(SL.cgidb[[2]])
SLG1.b <- c(SL.biogrid[[1]])
SLG2.b <- c(SL.biogrid[[2]])

SLG1 <- rbind(SL.cgidb, SL.biogrid)
SL_data <- data.frame(SLG1)
rm(SLG1)

SL_data <- unique(SL_data[c("X1", "X2")])

TSGs <- data.frame(Human_TSGs)
Oncogenes <- data.frame(ongene_human)
ongene_human <- read.delim("D:/Internship 2022 (Dr. Manjari)/ongene_human.txt")
View(ongene_human)
Human_TSGs <- read.delim("D:/Internship 2022 (Dr. Manjari)/Human_TSGs.txt")
View(Human_TSGs)
TSGs <- data.frame(Human_TSGs)
Oncogenes <- data.frame(ongene_human)

O1 <- cbind(intersect(Oncogenes$OncogeneID, SL_data$X1))
O2 <- cbind(intersect(Oncogenes$OncogeneID, SL_data$X2))
T1 <- cbind(intersect(TSGs$GeneID, SL_data$X1))
T2 <- cbind(intersect(TSGs$GeneID, SL_data$X2))
O1.SL1 <- cbind(subset(SL_data$X1, !SL_data$X1 %in% O1))
O2.SL2 <- cbind(subset(SL_data$X2, !SL_data$X2 %in% O2))
T1.SL1 <- cbind(subset(SL_data$X1, !SL_data$X1 %in% T1))
T2.SV2 <- cbind(subset(SL_data$X2, !SL_data$X2 %in% T2))
SLG1 <- data.frame(SL_data$X1)
SLG2 <- data.frame(SL_data$X2)
V1 <- c(O1.SL1)
SL1_NA_O1 <- sapply(SLG1$SL_data.X1,
                    function(x) replace(x, x %in% V1, NA))
V2 <- c(O2.SL2)
SL2_NA_O2 <- sapply(SLG2$SL_data.X2,
                    function(x) replace(x, x %in% V2, NA))
V3 <- c(T1.SL1)
SL1_NA_T1 <- sapply(SLG1$SL_data.X1,
                    function(x) replace(x, x %in% V3, NA))
V4 <- c(T2.SL2)
SL2_NA_T2 <- sapply(SLG2$SL_data.X2,
                    function(x) replace(x, x %in% V4, NA))
T2.SL2 <- cbind(subset(SL_data$X2, !SL_data$X2 %in% T2))
V4 <- c(T2.SL2)
SL2_NA_T2 <- sapply(SLG2$SL_data.X2,
                    function(x) replace(x, x %in% V4, NA))
SL_NA_pairs <- cbind(SL1_NA_O1, SL2_NA_O2, SL1_NA_T1, SL2_NA_T2)
View(SL_NA_pairs)
SL_NA_pairs <- data.frame(SL_NA_pairs)
SL.O1.O2 <- cbind(SL_NA_pairs$SL1_NA_O1, SL_NA_pairs$SL2_NA_O2)
SL.O1.T2 <- cbind(SL_NA_pairs$SL1_NA_O1, SL_NA_pairs$SL2_NA_T2)
SL.O2.T1 <- cbind(SL_NA_pairs$SL2_NA_O2, SL_NA_pairs$SL1_NA_T1)
SL.T1.T2 <- cbind(SL_NA_pairs$SL1_NA_T1, SL_NA_pairs$SL2_NA_T2)
SL.O1.O2 <- data.frame(SL.O1.O2)
SL.O1.T2 <- data.frame(SL.O1.T2)
SL.O2.T1 <- data.frame(SL.O2.T1)
SL.T1.T2 <- data.frame(SL.T1.T2)
View(SL.O1.O2)
O1.O2 <- unique(SL.O1.O2[c("X1", "X2")])
O1.T2 <- unique(SL.O1.T2[c("X1", "X2")])
O2.T1 <- unique(SL.O2.T1[c("X1", "X2")])
T1.T2 <- unique(SL.T1.T2[c("X1", "X2")])
View(O1.O2)
sum(is.na(O1.O2$X1))
sum(is.na(O1.O2$X2))
View(O1.T2)
sum(is.na(O1.T2$X2))
sum(is.na(O1.T2$X1))
View(O2.SL2)
View(O2.T1)
sum(is.na(O2.T1$X2))
sum(is.na(O2.T1$X1))
View(T1.T2)
sum(is.na(T1.T2$X1))
sum(is.na(T1.T2$X2))
x <- c('O-O', 'O-NA', 'O-T', 'T-NA', 'T-T')
y <- c(214, 801, 589, 1029, 482)
Plot_SL <- data.frame(x, y)
barplot(Plot_SL$y, names.arg=Plot_SL$x, ylim=c(0, 1100), xlab = "Type of interaction", ylab = "Number of SL pairs ", main = "Plot of SL pairs of CGIDB and Biogrid data")
