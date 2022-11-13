SL.data <- data.frame(All_synthetic_lethality_gene_pairs)
View(SL.data)

SL.data1 <-subset(SL.data, SL.data$V1!="our prediction" & SL.data$V1!="SLDB" & SL.data$V1!="Michael Boettcher")
SL.data1 <- cbind(SL.data1$V3, SL.data1$V4)
View(SL.data1)
SL.data1 <- unique(SL.data1)
SL.data1 <- data.frame(SL.data1)

Oncogene1 <- cbind(intersect(ongene_human$OncogeneID, SL.data1$X1))
Oncogene2 <- cbind(intersect(ongene_human$OncogeneID, SL.data1$X2))
TSG1 <- cbind(intersect(Human_TSGs$GeneID, SL.data1$X1))
TSG2 <- cbind(intersect(Human_TSGs$GeneID, SL.data1$X2))

SLG1 <- data.frame(SL.data1$X1)
OG1 <- data.frame(Oncogene1)

rows.in.SLG1.that.are.not.in.OG1  <- function(SLG1,OG1)
{
  SLG1.vec <- apply(SLG1, 1, paste, collapse = "")
  OG1.vec <- apply(OG1, 1, paste, collapse = "")
  SLG1.without.OG1.rows<- SLG1[!SLG1.vec %in% OG1.vec,]
  return(SLG1.without.OG1.rows)
}

OG.SLG1 <- rows.in.SLG1.that.are.not.in.OG1(SLG1, OG1)
OG.SLG1 <- cbind(OG.SLG1)
View(OG.SLG1)

SLG2 <- data.frame(SL.data1$X2)
OG2 <- data.frame(Oncogene2)

OG.SLG2 <- rows.in.SLG2.that.are.not.in.OG2  <- function(SLG2,OG2)
{
  SLG2.vec <- apply(SLG2, 1, paste, collapse = "")
  OG2.vec <- apply(OG2, 1, paste, collapse = "")
  SLG2.without.OG2.rows<- SLG2[!SLG2 %in% OG2.vec,]
  return(SLG2.without.OG2.rows)
}
OG.SLG2 <- rows.in.SLG2.that.are.not.in.OG2(SLG2, OG2)
OG.SLG2 <- cbind(OG.SLG2)
View(OG.SLG2)

SLG1_new_OG <- replace(SLG1, OG.SLG1 == SLG1 , "NA")
SLG2_new_OG <- replace(SLG2, OG.SLG2 == SLG2, "NA")
sum(SLG1_new_OG$SL.data1.X1 == "NA") 
sum(SLG2_new_OG == "NA") 


         