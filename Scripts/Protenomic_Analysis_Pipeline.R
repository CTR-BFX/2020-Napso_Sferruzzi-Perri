#!/usr/local/bin/Rscript
# R 3.6.2
#---------------------------------------------------------------------------------
# mass spectrum protenomic data analysis
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2020-Napso_Sferruzi-Perri
#
#
# Analysis Performed by Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


message("+---------------------------------------------------------------------------------------")
message("+                   Install useful packages and setup working directory                 ")
message("+---------------------------------------------------------------------------------------")

suppressPackageStartupMessages({
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("cowplot")
  library("Seurat")
  library("Matrix")
  library("useful")
  library("reshape2")
  library("biomaRt")
  library("scran")
  library("scater")
  library("SingleCellExperiment")
  library("bigmemory")
  library("mltools")
  library("rhdf5")
  library("monocle3")
  library("recommenderlab")
  library("readxl")
  library("GEOquery")
  library("UniProt.ws")
  library("rgdal")
  library("AnnotationDbi")
  library("org.Mm.eg.db")
  library("readxl")
})
options(future.globals.maxSize = 4000 * 1024^2)

Base.dir <- "/Users/xz289/Documents/CTR_ans48_0003"
Data.dir <- paste0(Base.dir, "/Original_Data")
Project <- "CTR_ans48_0003"
Out.dir <- paste0(Base.dir, "/Tina_Paper_Data")

message("+-----            General settings for three data sets           ----------------+")

ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)          
head(ensEMBL2id)
subEnsembl <- unique(ensEMBL2id[,2:3])
colnames(subEnsembl) <- c("MouseGeneName", "description")

## Check each data sets by using Uniprot retrive to get the proper UniprotKB entry name 

Revised.names <- paste0(c("F8VQ40", "F8WJE0", "E9PW66", "ES1", "CTGF", "SYHC", "CYR61",
                          "NIBL1","GRP78","GPR98","GPX41","YB039", "CS043", "CQ062", "GPX42",
                          "SK2L2", "USMG5", "CN166", "NIBAN", "K1468", "WISP1",
                          "F6VW30", "E9PW66","A0A0R4J026"), "_MOUSE")
NewRevi.names <- paste0(c("LAMA1", "SAMH1", "NP1L1", "GAL3A", "CCN2", "HARS1", "CCN1", 
                          "NIBA2", "BIP", "AGRV2","GPX4","MTLN","TRIR", "CYBC1","GPX4",
                          "MTREX", "ATPMD", "RTRAF", "NIBA1", "RELCH", "CCN4",
                          "1433T", "NP1L1", "FST"), "_MOUSE")
Rev.ProteinID <- c("P19137","Q60710","E9PW66","Q9D172","P29268","Q61035","P18406","Q8R1F1","P20029","Q8VHN7",
                   "O70325","Q8BT35","Q9D735","Q3TYS2","O70325","Q9CZU3","Q78IK2","Q9CQE8","Q3UW53","Q148V7",
                   "O54775","P68254", "P28656","P47931")

message("+--------------------------------------------------------------------------------------+")
message("+ Cultured_Primary: Read in different Raw Proteomic Data sets                          +")
message("+--------------------------------------------------------------------------------------+") 

PData3 <- read.csv(paste0(Data.dir, "/Cultured_Primary_Trophoblast_Data.csv"), header=T) ## 2131 proteins
PData3.sel <- PData3[, c(3,4,17,18:23,25)]
PData3.sel$ProteinID <-  unlist(lapply(as.character(PData3.sel$Accession), 
                                       function(x)  strsplit(x, split="[|]")[[1]][1]))
PData3.sel$Accession.Number <- unlist(lapply(as.character(PData3.sel$Accession), 
                                             function(x)  strsplit(x, split="[|]")[[1]][2]))
PData3.sel$Accession.New <- PData3.sel$Accession.Number
PData3.sel$MouseGeneName <- unlist(lapply(unlist(lapply(as.character(PData3.sel$Description),
                                                        function(x) strsplit(x, split="GN=")[[1]][2])),function(x) strsplit(x, split=" ")[[1]][1]))

for(i in 1:length(Revised.names)){
  rev.ind <- which(PData3.sel$Accession.Number==Revised.names[i])
  print(rev.ind)
  if(length(rev.ind)!=0){
    PData3.sel$Accession.New[rev.ind] <- NewRevi.names[i]
    PData3.sel$ProteinID[rev.ind] <- Rev.ProteinID[i]
  }
}

PData3.sel.QC <- subset(PData3.sel, X.Unique >=2 & X.10lgP >=12.7) 
## 1534, Protein-1533(GPX41,GPX42), Genes 1531 (Gnas, Naca, Gpx4), (Peptides <2, 597)
## Q9D735 (Trir) and Q3TYS2 (Cybc1) which has no gene names in the data,
## but found the gene names in uniport.

miss.ind <- which(is.na(PData3.sel.QC$MouseGeneName))
miss.ind.protein <- PData3.sel.QC[miss.ind,]

## Quality control, peptide >=2, 10lgP
PData3.sel.QC$MouseGeneName[miss.ind] <- c("Trir", "Cybc1")

## save the PData3.sel.QC for overlap use later
Protein.dat3c <- PData3.sel.QC[,c(13,14,11,10,4:8,2,3,1,12)]
colnames(Protein.dat3c) <- c("Accession.New", "MouseGeneName", "ProteinID", "Description", "Sample1", "Sample2", 
                             "Sample3","Sample4", "Sample5", "10lgP", "UniquePeptide", "Accession.Number.ori","Accession.Number")

## filter the proteins which expressed less than 3 samples
QC.submat <- PData3.sel.QC[,4:8]
QC.submatind <- ifelse(QC.submat==0,0,1)
QC.rind <- rowSums(QC.submatind)
PData3.sel.QC.GN.S4 <- PData3.sel.QC[QC.rind>=4,] ## 1208 (1206 unique genes)

PData3.final <- PData3.sel.QC.GN.S4[,c(13,14,11,10,4:8,2,3,1,12)]

colnames(PData3.final) <- colnames(Protein.dat3c) 
write.csv(PData3.final, file = paste0(Data.dir, "/Cultured_Trophoblast_N1208_G1206_Gnas_Naca_Filtered_Step2_Data.csv"), row.names=F, quote=T)


test3<- PData3.final[PData3.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
length(unique(test.sub$Accession.New)) ## 67
write.csv(test3, file = paste0(Data.dir, "/Cultured_Trophoblast_N67_G67_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)

message("+--------------------------------------------------------------------------------------+")
message("+ Conditonal_Medium:Read in different Raw Proteomic Data sets                          +")
message("+--------------------------------------------------------------------------------------+")

PData1 <- read.csv(paste0(Data.dir, "/Conditional_Medium_Data.csv"), header = T)
##  1467 x 14
PData1$Accession.Number.ori <- PData1$Accession.Number
Accession.Number <- PData1$Accession.Number
Protein.dat1 <- PData1[grepl("*_MOUSE*", Accession.Number)==T, ]
## 1445
Protein.dat1$Accession.Number <- unlist(lapply(as.character(Protein.dat1$Accession.Number.ori), function(x) strsplit(x, split=" ")[[1]][1]))
spindex <- which(grepl("^sp*", Protein.dat1$Accession.Number)==T)
Protein.dat1$Accession.Number[spindex] <- unlist(lapply(as.character(Protein.dat1$Accession.Number[spindex]), 
                                                        function(x) strsplit(x, split="[|]")[[1]][3]))

## Another two protein accession number with -DECOY
decoyindex <- which(grepl("*DECOY", Protein.dat1$Accession.Number)==T)
Protein.dat1$Accession.Number[decoyindex] <- unlist(lapply(as.character(Protein.dat1$Accession.Number[decoyindex]), 
                                                           function(x) strsplit(x, split="[-]")[[1]][1]))

Protein.dat1$MouseGeneName <- unlist(lapply(unlist(lapply(as.character(Protein.dat1$Identified.Proteins..1467.),
                                                          function(x) strsplit(x, split="GN=")[[1]][2])),function(x) strsplit(x, split=" ")[[1]][1]))
## 1445
## F6VW30, E9PW66,A0A0R4J026,F6QL70, D3Z536,A0A140T8L3, deleted by UniProt, april, 2020
## Ywhaq, Nap1l1,Fst, Gm17669,Gm8225,Rpl7a-ps5
## 1433T_MOUSE, NP1L1_MOUSE, FST_MOUSE

Protein.dat1$Accession.New <- Protein.dat1$Accession.Number

for(i in 1:length(Revised.names)){
  rev.ind <- which(Protein.dat1$Accession.Number==Revised.names[i])
  print(rev.ind)
  if(length(rev.ind)!=0){
    Protein.dat1$Accession.New[rev.ind] <- NewRevi.names[i]
  }
}


Protein.dat1c <- Protein.dat1[,c(18,17,4,10:14,16,5)]
colnames(Protein.dat1c) <- c("Accession.New", "MouseGeneName",  "Description", "Sample1", "Sample2", 
                             "Sample3","Sample4", "Sample5", "Accession.Number.ori","Accession.Number")

## 1445,--1443(KPYM, FBLN1),1441 (Psg16, Pkm, Fbln1,Tpm1)

missGindex <- which(is.na(Protein.dat1c$MouseGeneName))
## The protein TEST2  has no gene name corresponding.

Protein.dat1.newc <- Protein.dat1[-missGindex,c(18,17,4,10:14,16,5)] ## Remove TEST2 and othe useless columns.
## 1444
colnames(Protein.dat1.newc) <- colnames(Protein.dat1c)


## Filter the 4 out of 5 samples
QC.submat1 <- Protein.dat1.newc[,4:8]
QC.submatind1 <- ifelse(QC.submat1==0,0,1)
QC.rind1 <- rowSums(QC.submatind1)
Protein.dat1.new <- Protein.dat1.newc[QC.rind1>=4,]
## KPYM_MOUSE, sp|P52480|KPYM_MOUSE and sp|P52480-2|KPYM_MOUSE; D3YXQ6_MOUSE and D0VY58_MOUSE with gene Psg16
## 924
## Uniprot Retrive to get the UniProt ID.

## check
check.genes <- c("Ywhaq", "Nap1l1", "Fst", "Gm17669","Gm8225","Rpl7a-ps5")
Protein.dat1.new[unlist(lapply(check.genes, function(x) which(Protein.dat1.new$MouseGeneName==x))),]

Protein.dat1.uni <- read.csv(paste0(Data.dir, "/Conditional_Medium_UniProt_N924_filtered.csv"), header=T)[,1:2]
Protein.dat1.final <- merge(Protein.dat1.new,Protein.dat1.uni, by="Accession.New",all.x=T)
Protein.dat1.final[unlist(lapply(check.genes, function(x) which(Protein.dat1.final$MouseGeneName==x))),]
Protein.dat1.final[,11] <- as.character(Protein.dat1.final[,11])
Protein.dat1.final[6,11] <- "P68254"
## 924-protein 923 KPYM Gene 922 Psg16,Pkm

write.csv(Protein.dat1.final, file = paste0(Data.dir, "/Conditional_Medium_D924_G922_Psg16_Pkm_Filtered_Step2_Data.csv"),row.names =F)

test1 <- Protein.dat1.final[Protein.dat1.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
write.csv(test1, file = paste0(Data.dir, "/Conditional_Medium_Trophoblast_N29_G29_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)

message("+--------------------------------------------------------------------------------------+")
message("+ Sorted Trophoblast cell:Read in different Raw Proteomic Data sets                    +")
message("+--------------------------------------------------------------------------------------+")
## Remove the supernatant data
#Protein.dat2 <- read_excel(paste0(Data.dir, "/Sorted_Trophblast_Cells_Data.xlsx"))
#colnames(Protein.dat2) <- c("Accession.Number", "S3", "S4", "S5", "P3", "P4", "P5")
#Protein.dat2 <- Protein.dat2[,-c(2:4)]
#write.csv(Protein.dat2, file = paste0(Data.dir, "/Sorted_Trophblast_Cells_Data.csv"), row.names=F)
Protein.dat2 <- read.csv(paste0(Data.dir, "/Sorted_Trophblast_Cells_Data.csv"), header=T)

Revised.names <- paste0(c("F8VQ40", "F8WJE0", "E9PW66", "ES1", "CTGF", "SYHC", "CYR61", "NIBL1","GRP78","GRP98","GPX41","YB039"), "_MOUSE")
NewRevi.names <- paste0(c("LAMA1", "SAMH1", "NP1L1", "GAL3A", "CCN2", "HARS1", "CCN1", "NIBA2", "BIP", "AGRV2","GPX4","MTLN"), "_MOUSE")

for(i in 1:length(Revised.names)){
  rev.ind <- which(Protein.dat2$Accession.Number==Revised.names[i])
  print(rev.ind)
  Protein.dat2$Accession.Number[rev.ind] <- NewRevi.names[i]
  Protein.dat2
}

Protein.dat2c <- Protein.dat2
## 1142

## Filter the proteins at least exist in 2 samples
QC.mat <- Protein.dat2c[,2:4]
QC.mat.ind <- ifelse(QC.mat==0,0,1)
QC.mat.ind1 <- which(rowSums(QC.mat.ind)>=2)

Protein.dat2c.new <- Protein.dat2c[QC.mat.ind1,]
Protein.dat2c.new$Accession.New <- Protein.dat2c.new$Accession.Number
## Go to the Uniprot revive website to convert the entryname to Uniprot ID and GeneName
## Connect with the MouseGeneName using UniProt.ws library.
Protein.dat2.uni <- read.csv(paste0(Data.dir, "/Sorted_Trophoblast_UniProt_N682_filtered.csv"), header=T)[,1:4]

Protein.dat2.final <- merge(Protein.dat2c.new, Protein.dat2.uni, by = "Accession.New")
Protein.dat2.final$MouseGeneName <- as.character(Protein.dat2.final$MouseGeneName)
Protein.dat2.final[which(Protein.dat2.final$Accession.New=="H2B1K_MOUSE"),8] <- "Hist1h2bk"
## 654 protein, and  gene 653 Naca
colnames(Protein.dat2.final) <- c("Accession.New", "Accession.Number.ori", "P3", "P4", "P5", 
                                  "Accession.Number.uni", "ProteinID", "MouseGeneName")

write.csv(Protein.dat2.final, file = paste0(Data.dir, "/Sorted_TrophoblastP_N682_G681_Naca_Filtered_Step2_Data.csv"),row.names =F)

test2 <- Protein.dat2.final[Protein.dat2.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
write.csv(test2, file = paste0(Data.dir, "/Sorted_cell_TrophoblastP_N41_G41_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)
## 40
message("+---- Check the common proteins among three data sets with the same gene names ---------+")

Protein.dat3.final <- PData3.final
colnames(Protein.dat2.final) <- c("Accession.New", "Accession.Number",  "P3", "P4", "P5",
                                  "Accession.Number.uni", "ProteinID", "MouseGeneName")

over12 <- merge(Protein.dat1.final, Protein.dat2.final, by = "Accession.Number") ## 321
length(unique(over12$MouseGeneName.x)) ## 320
over12g <- merge(Protein.dat1.final, Protein.dat2.final, by = "MouseGeneName") ## 358
length(unique(over12g$MouseGeneName)) ## 356

nonover12g.12 <- over12g[over12g$MouseGeneName%in%over12$MouseGeneName.x==F,]
## these are having different protein name but with the same gene name.
write.csv(nonover12g.12, file = paste0(Data.dir, "/ConditionM_sorted_overlapG_noverlapPro_N48.csv"), row.names=F, quote=T) 


over13 <- merge(Protein.dat1.final, Protein.dat3.final, by = "Accession.Number")  ## 476
length(unique(over13$MouseGeneName.x)) ## 475
over13g <- merge(Protein.dat1.final, Protein.dat3.final, by = "MouseGeneName") ## 569
length(unique(over13g$MouseGeneName)) ## 567

nonover13g.13 <- over13g[over13g$MouseGeneName%in%over13$MouseGeneName.x==F,] 
write.csv(nonover13g.13, file = paste0(Data.dir, "/ConditionM_Culture_overlapG_noverlapPro_N93.csv"), row.names=F, quote=T) 

## 
over32 <- merge(Protein.dat3.final, Protein.dat2.final, by = "Accession.Number")  ## 527
length(unique(over32$MouseGeneName.x)) ## 526
over32g <- merge(Protein.dat3.final, Protein.dat2.final, by = "MouseGeneName") ## 507
length(unique(over32g$MouseGeneName)) ## 503

nonover32.32g <- over32[over32$MouseGeneName.y%in%over32g$MouseGeneName==F,]  
## These are the same protein but the gene names difference, need to change the gene names for sorted data and cultured data. make them consistent
write.csv(nonover32.32g, file = paste0(Data.dir, "/Culture_Sorted_overlapG_noverlapPro_N24.csv"), row.names=F, quote=T) 

checkMouseGeneNames <- unique(c(as.character(nonover32.32g$MouseGeneName.x), as.character(nonover32.32g$MouseGeneName.y)))
checkMouseGeneNames <- checkMouseGeneNames[order(checkMouseGeneNames)]
testens <- lapply(checkMouseGeneNames, function(x) which(ensEMBL2id$external_gene_name==x))
testens.len <- unlist(lapply(testens, length))
replace.genes <- checkMouseGeneNames[which(testens.len==0)]
replace.genes.new <- c("Aars", "Atp5a1", "Atp5b", "Cars","Dars","Ecpas","Eprs","Gars","H13", "H14","H13", "H14","H2bc14",
                       "Iars","Kars", "Lars","Mars1","Nars","Qars","Rars","Sars","Septin11","Septin2","Septin7","Septin9",
                       "Wars")
replaceGmat <- cbind(replace.genes,replace.genes.new)
gind.len1 <- NULL
Protein.dat1.final.new <- Protein.dat1.final
for(i in 1:length(replace.genes.new)){
  gind <- which(Protein.dat1.final.new$MouseGeneName==replace.genes[i])
  gind.lensub <- length(gind)
  gind.len1 <- c(gind.len1, gind)
  gind.len1
}
Protein.dat1.final.new$MouseGeneName <- as.character(Protein.dat1.final.new$MouseGeneName)
Protein.dat1.final.new$MouseGeneName[gind.len1] <- c("Atp5b","Septin11","Septin2","Septin7","Septin9")
##
gind.len2 <- NULL
Protein.dat2.final.new <- Protein.dat2.final
for(i in 1:length(replace.genes.new)){
  gind <- which(Protein.dat2.final.new$MouseGeneName==replace.genes[i])
  gind.lensub <- length(gind)
  gind.len2 <- c(gind.len2, gind)
  gind.len2
}
Protein.dat2.final.new$MouseGeneName <- as.character(Protein.dat2.final.new$MouseGeneName)
Protein.dat2.final.new$MouseGeneName[gind.len2] <- c("Aars", "Atp5a1", "Atp5b", "Cars","Dars","Eprs","Gars","H13", "H14",
                                                     "Iars","Kars", "Lars","Nars","Qars","Rars","Sars","Wars")

##
gind.len3 <- NULL
Protein.dat3.final.new <- Protein.dat3.final
for(i in 1:length(replace.genes.new)){
  gind <- which(Protein.dat3.final.new$MouseGeneName==replace.genes[i])
  gind.lensub <- length(gind)
  gind.len3 <- c(gind.len3, gind)
  gind.len3
}
Protein.dat3.final.new$MouseGeneName <- as.character(Protein.dat3.final.new$MouseGeneName)
Protein.dat3.final.new$MouseGeneName[gind.len3] <- c("Ecpas","H13", "H14","H2bc14","Mars1","Septin11","Septin2","Septin7","Septin9")


write.csv(Protein.dat1.final.new, file = paste0(Data.dir, "/Conditional_Medium_D924_G922_Psg16_Pkm_Filtered_Step2_correctGeNames_Data.csv"),row.names =F,quote=T)
write.csv(Protein.dat2.final.new, file = paste0(Data.dir, "/Sorted_TrophoblastP_N654_G653_Naca_Filtered_Step2_correctGeNames_Data.csv"),row.names =F, quote=T)
write.csv(Protein.dat3.final.new, file = paste0(Data.dir, "/Cultured_Trophoblast_N1208_G1206_Gnas_Naca_Filtered_Step2_correctGeNames_Data.csv"), row.names=F, quote=T)



message("+----                 Perform Overlap with the GEO public Data           --------+")

MousePublic <- read.csv(paste0(Data.dir, "/GEO_Control_Mouse_D3_N47936_GeneName_List.csv"), header=T)
colnames(MousePublic) <- "MouseGeneName"

Dat1.overlap1 <- Protein.dat1.final.new[Protein.dat1.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 908-907--906, Psg16, Pkm
Dat2.overlap1 <- Protein.dat2.final.new[Protein.dat2.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 621-621--620, Naca
Dat3.overlap1 <- Protein.dat3.final.new[Protein.dat3.final.new$MouseGeneName%in%MousePublic[,1]==T,]
## 1180-1180-1178, Gnas, Naca, 28 genes are not found in ensEMBLE and also not found in public Mouse.


message("+----                 Perform Homolog overlap           --------+")

HumanPublic <- read.csv(paste0(Data.dir, "/GEO_Control_Human_D8_N36552_GeneName_List.csv"), header=T)
HomoData <- read.csv(paste0(Data.dir, "/", Project, "-Ensemble_MGI_NCBI_human_mouse_Homolog_april_2020.csv"), header=T)
## 29524
HomoData.overPub <- HomoData[as.character(HomoData$HumanGeneName)%in%as.character(HumanPublic[,1])==T,]
## 26532
Filter.proteins <- read.csv(paste0(Data.dir, "/Filter_Proteins_List.csv"), header=T)

message("+----                 Perform Homolog overlap with Conditioned Medium          --------+")

Dat1.overlap2 <- merge(Protein.dat1.final.new, HomoData.overPub, by = "MouseGeneName")
## make HomoData.overPub unique for mouse and Human matching and remove all of the 1-more 

MusMultinames.dat1 <- unique(Dat1.overlap2$MouseGeneName[duplicated(Dat1.overlap2$MouseGeneName)==T])
MusMultinames.dat1.dupind <- lapply(MusMultinames.dat1, function(x) which(Dat1.overlap2$MouseGeneName==x))
MusMultinames.dat1.nodup <- Dat1.overlap2[-unlist(MusMultinames.dat1.dupind),]
MusMultinames.dat1.dup <-Dat1.overlap2[unlist(MusMultinames.dat1.dupind),]
MusMultinames.dat1.dup$order <- c(1:150)
## manually checking the one-to-one ortholog. 
MusMultinames.dat1.dup[,c(1,12,13)]
rm.homoHum.ind1 <- c(1,4,5,7,9,11,13,15,18,20,21,24,26,27,29,31,32:35,37,38,40,43,45,46,49,51,52,55,56,
                    58,59,60,62,64,67,68,70,72,74,76,78,79,80,82,85,87,89,91,93,94,96,97,99,100,102,
                    103,105,107,110,111,113,114,117,119,123,124,125,128,130,132,133,135,138,140,141,
                    144,145,148,150)
MusMultinames.dat1.dupnew <- MusMultinames.dat1.dup[-rm.homoHum.ind1,-13]

Dat1.overlap2.new <- rbind(MusMultinames.dat1.nodup, MusMultinames.dat1.dupnew)
## 876-875, Pkm 
Dat1.overlap2.new.filter <- Dat1.overlap2.new[Dat1.overlap2.new$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
##859-858, Pkm

Dat1.MHoverlap <- unique(merge(Dat1.overlap1, Dat1.overlap2.new, by = "MouseGeneName"))[-c(563:564),]
## 875-874, Pkm

Dat1.Munique <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.MHoverlap$MouseGeneName==F,]

message("+----                 Perform Homolog overlap with Sorted cells          --------+")

Dat2.overlap2 <- merge(Protein.dat2.final.new, HomoData.overPub, by = "MouseGeneName")
MusMultinames.dat2 <- unique(Dat2.overlap2$MouseGeneName[duplicated(Dat2.overlap2$MouseGeneName)==T])
MusMultinames.dat2.dupind <- lapply(MusMultinames.dat2, function(x) which(Dat2.overlap2$MouseGeneName==x))
MusMultinames.dat2.nodup <- Dat2.overlap2[-unlist(MusMultinames.dat2.dupind),]
MusMultinames.dat2.dup <-Dat2.overlap2[unlist(MusMultinames.dat2.dupind),]
MusMultinames.dat2.dup$order <- c(1:93)
## manually checking the one-to-one ortholog. 
MusMultinames.dat2.dup[,c(1,9,10)]

rm.homoHum.ind2 <- c(1,3,6,7,10,12,14,16,18,19,20,23,24,27,29,30,32,35,36,39,41,42,44,46,49,50,51,54,55,
                    57:60,63,64,68,70,73,75,76,78,81,82,85,87,89,91,93)
MusMultinames.dat2.dupnew <- MusMultinames.dat2.dup[-rm.homoHum.ind2,-10]

Dat2.overlap2.new <- rbind(MusMultinames.dat2.nodup, MusMultinames.dat2.dupnew)
## 631-630, Naca 
Dat2.overlap2.new.filter <- Dat2.overlap2.new[Dat2.overlap2.new$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
##627-626, Naca

Dat2.MHoverlap <- unique(merge(Dat2.overlap1, Dat2.overlap2.new, by = "MouseGeneName"))[-c(316:317),]
## 615-614, Naca

Dat2.Munique <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.MHoverlap$MouseGeneName==F,]

message("+----                 Perform Homolog overlap with Sorted cells cultured         --------+")

Dat3.overlap2 <- merge(Protein.dat3.final.new, HomoData.overPub, by = "MouseGeneName")
MusMultinames.dat3 <- unique(Dat3.overlap2$MouseGeneName[duplicated(Dat3.overlap2$MouseGeneName)==T])
MusMultinames.dat3.dupind <- lapply(MusMultinames.dat3, function(x) which(Dat3.overlap2$MouseGeneName==x))
MusMultinames.dat3.nodup <- Dat3.overlap2[-unlist(MusMultinames.dat3.dupind),]
MusMultinames.dat3.dup <-Dat3.overlap2[unlist(MusMultinames.dat3.dupind),]
MusMultinames.dat3.dup$order <- c(1:168)
## manually checking the one-to-one ortholog. 
MusMultinames.dat3.dup[,c(1,14,15)]
rm.homoHum.ind3 <- c(1,3,6,7,9,11, 14,12,13,16,18,19,21,23,24,25,28,29,30,32,34,36,38,39,42,43,45,48,50,51,
                     54,56,58,59,61,65,68,69,71,72,74:76,79,81,82,85,87,89,90,93,94:96,98,99,101:103,105,
                     107,110,112,115,117,119,121,122,123,126,128,129:131,133:136,137,140,141,144,145,148,
                     150,151,153,156,158,159,162,164,166,167)
MusMultinames.dat3.dupnew <- MusMultinames.dat3.dup[-rm.homoHum.ind3,-15]

Dat3.overlap2.new <- rbind(MusMultinames.dat3.nodup, MusMultinames.dat3.dupnew)
## 1170-1168, Naca, Gnas
Dat3.overlap2.new.filter <- Dat3.overlap2.new[Dat3.overlap2.new$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
##1167-1165, Naca, Gnas

Dat3.MHoverlap <- unique(merge(Dat3.overlap1, Dat3.overlap2.new, by = "MouseGeneName"))[-c(413:414,647:648),]
## 1170-1168, Naca, Gnas
Dat3.Munique <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.MHoverlap$MouseGeneName==F,]


message("+---- Filter out the Prl family proteins ------------+")

# Tina_MHlist <- read_excel(paste0(Out.dir,"/Secreted_Final_List_Tina.xlsx"), sheet=4)
# Tina_MUlist <- read_excel(paste0(Out.dir,"/Secreted_Final_List_Tina.xlsx"), sheet=7)
Dat1.MHoverlap.filter <- Dat1.MHoverlap[Dat1.MHoverlap$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
Dat1.Munique.filter <- Dat1.Munique[Dat1.Munique$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]

Dat2.MHoverlap.filter <- Dat2.MHoverlap[Dat2.MHoverlap$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
Dat2.Munique.filter <- Dat2.Munique[Dat2.Munique$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
Dat2.Munique.filter.new <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.MHoverlap.filter$MouseGeneName==F,]

Dat3.MHoverlap.filter <- Dat3.MHoverlap[Dat3.MHoverlap$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]
Dat3.Munique.filter <- Dat3.Munique[Dat3.Munique$MouseGeneName%in%Filter.proteins$MouseGeneName==F,]


message("+----  Suggestions from Tina to remove some genes, due to multiple homologus with human       ----+")


test1.mh.over <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==T,]
test1.mh.nover <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==F,]

test2.mh.over <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==T,]
test2.mh.nover <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==F,]

test3.mh.over <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==T,]
test3.mh.nover <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==F,]

message("+---write out each step proteins list for flowchart Data1, 2,3 ----------------------+")

write.csv(PData3.sel.QC, file = paste0(Data.dir, "/Cultured_Trophoblast_N1534_G1531_Gnas_Naca_Gpx4_Filtered_pep2_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap1, file = paste0(Data.dir, "/Cultured_Trophoblast_N1180_G1178_Gnas_Naca_Filtered_pep2_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap2.new, file = paste0(Data.dir, "/Cultured_Trophoblast_N1170_G1168_Gnas_Naca_Filtered_pep2_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.Munique,file = paste0(Data.dir, "/Cultured_Trophoblast_N10_G10_Filtered_pep2_Step4_MusUni_Data.csv"), row.names=F, quote=T )

write.csv(Protein.dat2c, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N1142_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.overlap1, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N621_G620_Naca_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.MHoverlap.filter, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N611_G610_Naca_Filtered_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.Munique.filter.new, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N10_G10_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)


write.csv(Protein.dat1c, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N1445_G1441_Psg16_Pkm_Tpm1_Fbln1_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat1.overlap1, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N908_G906_Psg16_Pkm_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat1.MHoverlap, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N876_G875_Psg16_Filtered_Step3_HumanPub_Data.csv"), row.names=F, quote=T)
write.csv(Dat1.Munique, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N33_G32_Psg16_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)


message("+------- Secreted process: signalP + GO pathway Uniprot keyword  ---------------------------+")

sect.dir <- "/storage/CTR-Projects/CTR_ans48/CTR_ans48_0003/SignalP/"
setwd(sect.dir)
library("UniProt.ws")
availableUniprotSpecies("musculus")
up <- UniProt.ws(taxId=10090)
## List the possible values for columns
columns(up)## List the possible values for keytypes
keytypes(up)## get some values back
## list all possible keys of type entrez gene ID.
## (this process is not instantaneous)
if(interactive()){
  keys <- head(keys(up, keytype="UNIPROTKB"))
  keys}


message("+-------            generate fastq file for proteins.           -----------------+")

## Data1 

Data1_MHP_common <- read.csv("Conditioned_Medium_Trophoblast_N881_G879_Psg16_Pkm_Filtered_Step3_HumPub_Data.csv", header=T)
Data1_MHP_common_fasta<- select(up, keys=Data1_MHP_common$ProteinID, 
                                columns=c("UNIPROTKB", "SEQUENCE"),
                                keytype="UNIPROTKB")

colnames(Data1_MHP_common_fasta) <- c("ProteinID", "SEQUENCE")

Data1_MHPC_mer <- merge(Data1_MHP_common_fasta, Data1_MHP_common, by = "ProteinID")
Data1_MHPC_mer <- Data1_MHPC_mer [,c(1,2,3,11,12)]

Data1_fasta_names_1 <- apply(Data1_MHPC_mer, 1, function(x) 
  paste0(">UNIPROT:", x[3], " ", x[1], " ", x[4]))
Data1_fasta_names_2 <- cbind(Data1_fasta_names_1, Data1_MHPC_mer[,2])
Data1_fasta_names_3 <- apply(Data1_fasta_names_2, 1, function(x) paste(x[1], x[2], sep="\n", collapse= "" ))

write.table(Data1_fasta_names_3, file = "Data1_Mouse_Human_Public_Common.fasta", row.names=F, col.names=F, quote=F) 

## Data2 

Data2_MHP_common <- read.csv("Sorted_cell_Trophoblast_N634_G633_Naca_Filtered_Step3_HumPub_Data.csv", header=T)
Data2_MHP_common_fasta<- select(up, keys=Data2_MHP_common$ProteinID, 
                                columns=c("UNIPROTKB", "SEQUENCE"),
                                keytype="UNIPROTKB")

colnames(Data2_MHP_common_fasta) <- c("ProteinID", "SEQUENCE")

Data2_MHPC_mer <- merge(Data2_MHP_common_fasta, Data2_MHP_common, by = "ProteinID")
Data2_MHPC_mer <- Data2_MHPC_mer [,c(1,2,3,4,11)]

Data2_fasta_names_1 <- apply(Data2_MHPC_mer, 1, function(x) 
  paste0(">UNIPROT:", x[3], " ", x[1], " ", x[4]))
Data2_fasta_names_2 <- cbind(Data2_fasta_names_1, Data2_MHPC_mer[,2])
Data2_fasta_names_3 <- apply(Data2_fasta_names_2, 1, function(x) paste(x[1], x[2], sep="\n", collapse= "" ))

write.table(Data2_fasta_names_3, file = "Data2_Mouse_Human_Public_Common.fasta", row.names=F, col.names=F, quote=F) 

##  Data3

Data3_MHP_common <- read.csv("Cultured_Trophoblast_N1172_G1170_Gnas_Naca_Filtered_pep2_Step3_HumPub_Data.csv", header=T)
Data3_MHP_common_fasta<- select(up, keys=Data3_MHP_common$ProteinID, 
                                columns=c("UNIPROTKB", "SEQUENCE"),
                                keytype="UNIPROTKB")
colnames(Data3_MHP_common_fasta) <- c("ProteinID", "SEQUENCE")

Data3_MHPC_mer <- merge(Data3_MHP_common_fasta, Data3_MHP_common, by = "ProteinID")
Data3_MHPC_mer <- Data3_MHPC_mer[,c(1,2,3,13,14)]

Data3_fasta_names_1 <- apply(Data3_MHPC_mer, 1, function(x) 
  paste0(">UNIPROT:", x[3], " ", x[1], " ", x[4]))
Data3_fasta_names_2 <- cbind(Data3_fasta_names_1, Data3_MHPC_mer[,2])
Data3_fasta_names_3 <- apply(Data3_fasta_names_2, 1, function(x) paste(x[1], x[2], sep="\n", collapse= "" ))

write.table(Data3_fasta_names_3, file = "Data3_Mouse_Human_Public_Common.fasta", row.names=F, col.names=F, quote=F) 

message("+------ Run signalP for the above three fasta files   -------------+")

for(names in c("Data1", "Data2", "Data3")){
  
  system(paste0("/storage/Software/packages/signalp-4.1/signalp -t euk -f short ", names, "_Mouse_Human_Public_Common.fasta >",
                names,"_fsa.short_out_SignalP"))
}

message("+------------ Check the UniProt GO pathway with the key word extracellular --------------+")

## Data1 
Dat1.pathway <- select(up, keys=Data1_MHP_common$ProteinID, 
                       columns=c("UNIPROTKB", "GO"),
                       keytype="UNIPROTKB")
selword <- "extracellular"
go.ind1 <- grep(selword, Dat1.pathway[,2])
Dat1.pathway.sub <- Dat1.pathway[go.ind1,]
Dat1.pathway$go.ind <- ifelse(Dat1.pathway$UNIPROTKB%in%Dat1.pathway.sub$UNIPROTKB==T, "Y", "N")

Dat1.pathway.sec <- Dat1.pathway[,c(1,3)]
colnames(Dat1.pathway.sec) <- c("ProteinID", "go.ind")

Dat1.signalP <- read.table("Data1_fsa.short_out_SignalP", header = F)[,c(1,10)]
Dat1.signalP[,1] <- gsub("UNIPROT:", "", Dat1.signalP[,1])
colnames(Dat1.signalP) <- c("Accession.New", "SignalP")

Dat1.signalP.mer <- merge(Dat1.signalP, Data1_MHP_common, by = "Accession.New")
Dat1.signalP.pathway.mer <- merge(Dat1.signalP.mer, Dat1.pathway.sec, by = "ProteinID")

Dat1.signalP.pathway.mer.new <- unique(Dat1.signalP.pathway.mer) ## KPYM with its splice. 881

## Data2 
Dat2.pathway <- select(up, keys=Data2_MHP_common$ProteinID, 
                       columns=c("UNIPROTKB", "GO"),
                       keytype="UNIPROTKB")

go.ind2 <- grep(selword, Dat2.pathway[,2])
Dat2.pathway.sub <- Dat2.pathway[go.ind2,]
Dat2.pathway$go.ind <- ifelse(Dat2.pathway$UNIPROTKB%in%Dat2.pathway.sub$UNIPROTKB==T, "Y", "N")

Dat2.pathway.sec <- Dat2.pathway[,c(1,3)]
colnames(Dat2.pathway.sec) <- c("ProteinID", "go.ind")

Dat2.signalP <- read.table("Data2_fsa.short_out_SignalP", header = F)[,c(1,10)]
Dat2.signalP[,1] <- gsub("UNIPROT:", "", Dat2.signalP[,1])
colnames(Dat2.signalP) <- c("Accession.New", "SignalP")

Dat2.signalP.mer <- merge(Dat2.signalP, Data2_MHP_common, by = "Accession.New")
Dat2.signalP.pathway.mer <- merge(Dat2.signalP.mer, Dat2.pathway.sec, by = "ProteinID")

Dat2.signalP.pathway.mer.new <- unique(Dat2.signalP.pathway.mer) ##  634

## Data3 
Dat3.pathway <- select(up, keys=Data3_MHP_common$ProteinID, 
                       columns=c("UNIPROTKB", "GO"),
                       keytype="UNIPROTKB")

go.ind3 <- grep(selword, Dat3.pathway[,2])
Dat3.pathway.sub <- Dat3.pathway[go.ind3,]
Dat3.pathway$go.ind <- ifelse(Dat3.pathway$UNIPROTKB%in%Dat3.pathway.sub$UNIPROTKB==T, "Y", "N")

Dat3.pathway.sec <- Dat3.pathway[,c(1,3)]
colnames(Dat3.pathway.sec) <- c("ProteinID", "go.ind")

Dat3.signalP <- read.table("Data3_fsa.short_out_SignalP", header = F)[,c(1,10)]
Dat3.signalP[,1] <- gsub("UNIPROT:", "", Dat3.signalP[,1])
colnames(Dat3.signalP) <- c("Accession.New", "SignalP")

Dat3.signalP.mer <- merge(Dat3.signalP, Data3_MHP_common, by = "Accession.New")
Dat3.signalP.pathway.mer <- merge(Dat3.signalP.mer, Dat3.pathway.sec, by = "ProteinID")

Dat3.signalP.pathway.mer.new <- unique(Dat3.signalP.pathway.mer) ##  1172

## Summary of the secreted number for both approaches in each data sets.

Dat1.signalP.pathway.mer.new$secreted <- ifelse(Dat1.signalP.pathway.mer.new$go.ind=="Y", "Y", "N")
Dat2.signalP.pathway.mer.new$secreted <- ifelse(Dat2.signalP.pathway.mer.new$SignalP=="Y"|Dat2.signalP.pathway.mer.new$go.ind=="Y", "Y", "N")
Dat3.signalP.pathway.mer.new$secreted <- ifelse(Dat3.signalP.pathway.mer.new$SignalP=="Y"|Dat3.signalP.pathway.mer.new$go.ind=="Y", "Y", "N")

## 312, 128, 214
Dat1.signalP.pathway.mer.newsub <- subset(Dat1.signalP.pathway.mer.new, secreted == "Y")
Dat2.signalP.pathway.mer.newsub <- subset(Dat2.signalP.pathway.mer.new, secreted == "Y")
Dat3.signalP.pathway.mer.newsub <- subset(Dat3.signalP.pathway.mer.new, secreted == "Y")

message("+-------Tina filtered and double checked  secreted list     -----------+")

Tina_secreted_final <- read_excel("All_secreted_lists_final.xlsx", sheet=1)
PDat1.secr <- Tina_secreted_final[1:257,1:2]
colnames(PDat1.secr) <- c("Accession.New", "MouseGeneName")
PDat2.secr <- Tina_secreted_final[1:105,3:4]
colnames(PDat2.secr) <- c("Accession.New", "MouseGeneName")
PDat3.secr <- Tina_secreted_final[1:158,5:6]
colnames(PDat3.secr) <- c("Accession.New", "MouseGeneName")

## Data1
PDat1.secr.mer <- merge(PDat1.secr, Dat1.signalP.pathway.mer.new, by = "Accession.New")
PDat1.secr.nmer.1 <- PDat1.secr[PDat1.secr$Accession.New%in%Dat1.signalP.pathway.mer.new$Accession.New==F,]
PDat1.secr.nmer.2 <- Dat1.signalP.pathway.mer.newsub[Dat1.signalP.pathway.mer.newsub$Accession.New%in%PDat1.secr$Accession.New==F,]

## manually checking the PDat1.secr.nmer.2 and found the following proteins are also secreted.
add1.sec.names <- paste0(c("1433S", "ROA2", "FRIH", "HS90B", "GELS", "LAMP2", "YBOX1", "NHRF1", "HS105",
                           "ISG15", "Q8CFH0", "ASL", "OLR1", "ERAP1", "INAR2", "PLSL", "FSCN1", "ASC"), "_MOUSE")
all1.sec.names <- c(add1.sec.names, PDat1.secr.mer$Accession.New)
add1.sec.index <- unique(unlist(lapply(all1.sec.names, function(x) which(Dat1.signalP.pathway.mer.new$Accession.New==x)) ))
add1.sec.index <- add1.sec.index[order(add1.sec.index)]

Data1_secreted_final <-  unique(Dat1.signalP.pathway.mer.new[add1.sec.index,]) ## 270
## There are six from Tina_secreted_final which have go.ind and signalP is "N',
## A0A1B0GSG5_MOUSE, GSTM2_MOUSE, GSTO1_MOUSE, KAPCA_MOUSE, PLBL2_MOUSE
## I think PLBL2_MOUSE is not secreted, the rest of with either both Y or at least one of them Y.
rm.ind1 <- which(Data1_secreted_final$Accession.New =="PLBL2_MOUSE")
Data1_secreted_final <- Data1_secreted_final[-rm.ind1,]

write.csv(Data1_secreted_final, file = "Conditioned_Medium_Trophoblast_secretedList_SignalP_GO_Tina_N270_G269_Pkm.csv", row.names=F, quote=T)

## Data2
rm.ind2 <- which(duplicated(PDat2.secr$Accession.New)==T)
## H2B1F_MOUSE   Hist1h2bj    
## H2B1F_MOUSE   Hist1h2bl    
## H2B1F_MOUSE   Hist1h2bn    
## H2B1C_MOUSE   Hist1h2be    
## H2B1C_MOUSE   Hist1h2bg 
PDat2.secr <- PDat2.secr[-rm.ind2,]
PDat2.secr.mer <- merge(PDat2.secr, Dat2.signalP.pathway.mer.new, by = "Accession.New") ## 98
PDat2.secr.nmer.1 <- PDat2.secr[PDat2.secr$Accession.New%in%Dat2.signalP.pathway.mer.new$Accession.New==F,]
PDat2.secr.nmer.2 <- Dat2.signalP.pathway.mer.newsub[Dat2.signalP.pathway.mer.newsub$Accession.New%in%PDat2.secr$Accession.New==F,]

## manually checking the PDat2.secr.nmer.2 and found the following proteins are also secreted.
add2.sec.names <- paste0(c("BIP", "FSCN1", "ITAV", "PLSL"), "_MOUSE")
all2.sec.names <- c(add2.sec.names, PDat2.secr.mer$Accession.New)
add2.sec.index <- unique(unlist(lapply(all2.sec.names, function(x) which(Dat2.signalP.pathway.mer.new$Accession.New==x)) ))
add2.sec.index <- add2.sec.index[order(add2.sec.index)]

Data2_secreted_final <-  unique(Dat2.signalP.pathway.mer.new[add2.sec.index,]) ## 102
## There are three from Tina_secreted_final which have go.ind and signalP is "N',
## "CAPG_MOUSE"  "VIGLN_MOUSE" "RINI_MOUSE" 
## after manually checking they are secreted

write.csv(Data2_secreted_final, file = "Sorted_cell_Trophoblast_secretedList_SignalP_GO_Tina_N102_G102.csv", row.names=F, quote=T)

## Data3
PDat3.secr.mer <- merge(PDat3.secr, Dat3.signalP.pathway.mer.new, by = "Accession.New") ## 158
PDat3.secr.nmer.1 <- PDat3.secr[PDat3.secr$Accession.New%in%Dat3.signalP.pathway.mer.new$Accession.New==F,]
PDat3.secr.nmer.2 <- Dat3.signalP.pathway.mer.newsub[Dat3.signalP.pathway.mer.newsub$Accession.New%in%PDat3.secr$Accession.New==F,]

## manually checking the PDat3.secr.nmer.2 and found the following proteins are also secreted.
add3.sec.names <- paste0(c("ANXA3", "C1QBP", "DHB12", "ROA2", "BIP", "HMGB2", "CD81", "FSCN1", "ITAV", 
                           "PLSL", "ERP29", "LXN", "NHRF1", "AQP1", "PGBM", "DPYL3", "DNJC9", "TINAL", 
                           "HNRPM"), "_MOUSE")
all3.sec.names <- c(add3.sec.names, PDat3.secr.mer$Accession.New) ## 177
add3.sec.index <- unique(unlist(lapply(all3.sec.names, function(x) which(Dat3.signalP.pathway.mer.new$Accession.New==x)) ))
add3.sec.index <- add3.sec.index[order(add3.sec.index)]

Data3_secreted_final <-  unique(Dat3.signalP.pathway.mer.new[add3.sec.index,]) ## 177
## There are three from Tina_secreted_final which have go.ind and signalP is "N',
## "CAPG_MOUSE"  "VIGLN_MOUSE" "RINI_MOUSE" 
## after manually checking they are secreted

write.csv(Data3_secreted_final, file = "Cultured_Trophoblast_secretedList_SignalP_GO_Tina_N177_G177.csv", row.names=F, quote=T)


message("+------  get an unique secreted list for all three data sets -------------+")

Data1_sec_sub <- Data1_secreted_final[,c(4,11)]
Data2_sec_sub <- Data2_secreted_final[,c(12,4)]
Data3_sec_sub <- Data3_secreted_final[,c(4,14)] 

colnames(Data1_sec_sub) <- c("MouseGeneName", "Accession")
colnames(Data2_sec_sub) <- c("MouseGeneName", "Accession")
colnames(Data3_sec_sub) <- c("MouseGeneName", "Accession")

secmer.all <- unique(rbind(Data1_sec_sub, Data2_sec_sub, Data3_sec_sub)) ## 368
secmer.all[308,1] <- "Hist1h2bk"
dupGenes <- secmer.all[which(duplicated(secmer.all[,1])==T),1]
checkdup <- lapply(dupGenes, function(x) secmer.all[which(secmer.all[,1]==x),])
replace.ass <- c("NUCB2_MOUSE (+1)", "HA1B_MOUSE (+1)", "BIP_MOUSE", "GRN_MOUSE (+1)" , "ACBP_MOUSE (+1)",
                 "LEG3_MOUSE (+1)", "LAMA1_MOUSE")
replace.ass.new <- c("NUCB2_MOUSE", "HA1B_MOUSE", "GRP78_MOUSE", "GRN_MOUSE",  "ACBP_MOUSE", "LEG3_MOUSE",
                     "F8VQ40_MOUSE")
reptest <- lapply(replace.ass, function(x) which(secmer.all[,2]==x))
secmer.all[unlist(reptest),2] <- replace.ass.new
secmer.all.final <- unique(secmer.all) ## 362
write.csv(secmer.all.final, file = "Cultured_ConditionM_SortedC_Trophoblast_secretedList_N360_G333.csv", row.names = F, quote=T)

message("+----------          SingleCell Checking and Heatmap      -----------------------+")

suppressPackageStartupMessages({
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("cowplot")
  library("Seurat")
  library("Matrix")
  library("useful")
  library("reshape2")
  library("biomaRt")
  library("scran")
  library("scater")
  library("SingleCellExperiment")
  library("bigmemory")
  library("mltools")
  library("rhdf5")
  library("monocle3")
  library("recommenderlab")
  library("readxl")
})
options(future.globals.maxSize = 4000 * 1024^2)

Base.dir <- "/storage/CTR-Projects/CTR_ans48/CTR_ans48_0003/scRNASeq_CS_Wang_2018"
Project <- "CTR_ans48_0003"

setwd(Base.dir)

message("+--------------------------------------------------------------------------------------------------------+")
message("+-----------------                   Seurat DEG Analysis, res0.6                            -------------+")
message("+--------------------------------------------------------------------------------------------------------+")

load("GSE89497_matrix_res0.6_april_2020.RData")
Mouse_Human_HomList <- read.csv("CTR_ans48_0003-Ensemble_MGI_NCBI_human_mouse_Homolog.csv", header=T)

Final.secret <- secmer.all.final
Secret.merE <- merge(Final.secret, Mouse_Human_HomList , by = "MouseGeneName", all.x=T) ## 409

Secret.merE.new <- unique(Secret.merE[,-2])  ## 382

## Sort out the duplication Human Genes
Secret.merE.dupIndex <- which(duplicated(Secret.merE.new$MouseGeneName)==T)
Secret.merE.dupSub <- Secret.merE.new[c(Secret.merE.dupIndex),]
dupGeneName <- Secret.merE.dupSub$MouseGeneName
dupSet <- lapply(dupGeneName, function(x) Secret.merE[Secret.merE$MouseGeneName==x,])
dropInd <- c(13, 30,45,52,88,97,100,140,146,152,156,158,164,184,199,221:224,243,244,280,283:284,293,295,
             310,312,313,315,319:337,345)
## Psg22 dropped due to the homolog has multiple in human. which I can only found the homolog in ensemble not in MGI and NCBI.

Secret.merE.new.1 <- unique(Secret.merE.new[-dropInd, ]) ## 332, Mouse Gene 294, Human 325
## only need the  positive ones to check the overlap genes cell type enrichment.

GSE89497_matrix.su@meta.data$FINAL_CLUSTERS <- GSE89497_matrix.su@meta.data$RNA_snn_res.0.6

new.cluster.ids <- c("EVT_8w_2", "CTB_8w_3", "EVT_8w_1", "EVT_24w_1", "Macro_1", "CTB_8w_2", "Macro_2",
                     "EVT_8w_12", "Blood_cell", "STB_8w", "EVT_8w_3", "Mes_2", "EVT_24w_2", "Mes_1", "CTB_8w_1")
names(new.cluster.ids) <- levels(GSE89497_matrix.su)
GSE89497_matrix.su <- RenameIdents(GSE89497_matrix.su , new.cluster.ids)

message("+--------------------------------------------------------------------------------------------------------+")
message("+-----------------                   Heatmap                                                -------------+")
message("+--------------------------------------------------------------------------------------------------------+")

library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(lattice)

selected.genes <- unique(Secret.merE.new.1[,2]) ## 325

selected_genes <- as.character(selected.genes[selected.genes %in%  rownames(GSE89497_matrix.su@assays[["RNA"]]@scale.data) ]) 
## 301 overlap
selected_metadata <- GSE89497_matrix.su@meta.data
scale_data <- GetAssayData(GSE89497_matrix.su, assay = "RNA", slot = "scale.data")
matrix <- as.matrix(scale_data[rownames(scale_data) %in% selected_genes,])
matrix <- as.data.frame(t(matrix))
tmatrix <- t(matrix)

matrix$Split <- selected_metadata$Type[match(rownames(selected_metadata), rownames(matrix))]
split.List <- unique(as.character(matrix$Split))
testList <- lapply(split.List, function(x) matrix[matrix$Split==x,-ncol(matrix)])
testList.meanmat <- lapply(testList, function(x) colMeans(apply(x, 2, as.numeric)))

scheat.mat <- cbind(testList.meanmat[[1]], testList.meanmat[[2]], testList.meanmat[[3]], testList.meanmat[[4]],
                    testList.meanmat[[5]])
colnames(scheat.mat) <- c("HE24W_EVT", "HE8W_CTB",  "HE8W_EVT",  "HE8W_STB",  "HE8W_STR")
colnames(scheat.mat) <- c("EVT-24w", "CTB-8w",  "EVT-8w",  "STB-8w",  "STR-8w")
#write.csv(scheat.mer, file = "scRNA_overlap_geneList_expression.csv")


## Heatmap for 301 genes, clustering

pdf("Secreted_N325_Psg22Drop_SelN301_scRNA_Heatmap_June_2020_addEVT24_newColor.pdf", width = 3, height = 40)
col_fn2 <- colorRamp2(c(-4,-2, 0, 2,4),  c("yellow4", "yellow", "white", "blue", "blue4"))
#col_fn2 <- colorRamp2(c(-4,-2, 0, 2,4),  c("blue", "cyan", "pink", "orange", "red"))

Heatmap(as.matrix(scheat.mat[,c(4,2,3,1)]), col = col_fn2, 
        cluster_columns = F, cluster_rows = T, row_title=NULL, 
        show_row_dend = F,    show_row_names = TRUE, 
        show_column_names = TRUE, row_dend_reorder = TRUE,
        row_km = 10, row_names_gp = gpar(col = rep("black",10), fontsize = rep(5,10)),
        #row_split = 8, row_names_gp = gpar(col = rep("black",8), fontsize = rep(5,8)),
        heatmap_legend_param = list(title = "Expression", title_gp = gpar(col = "black", fontsize = 6),
                                    legend_height = unit(2, "cm"), title_position = "topleft", 
                                    grid_width = unit(0.2, "cm")),  
        column_dend_height = unit(2, "cm"), column_names_centered = TRUE, column_names_rot = 45,
        row_dend_width = unit(4, "cm") , row_title_rot = 0) 
dev.off()


message("+---- Heatmap for selected significant genes in the groups ---------+")

scheat.dat <- as.data.frame(scheat.mat)
scheat.dat$HumanGeneName <- rownames(scheat.dat)
scheat.mer <- merge(scheat.dat, Secret.merE.new, by = "HumanGeneName")
## select a subset of genes to plot which enriched in STB, CTB, EVT 8 wks.
scheat.sub <- subset(scheat.mer, scheat.mer[,3]>=0.6 | scheat.mer[,4]>=0.6 | scheat.mer[,5] >= 1 )  ## 64
scheat.submat <- subset(scheat.mat, scheat.mat[,2]>=0.6 | scheat.mat[,3]>=0.6 | scheat.mat[,4] >= 1 ) ## 63
scheat.submat1 <- subset(scheat.mat, scheat.mat[,2]>=1 | scheat.mat[,3]>=1 | scheat.mat[,4] >= 1 | scheat.mat[,5] >= 1) ## 61
## 38
colnames(scheat.submat1) <- c("EVT-24w", "CTB-8w", "EVT-8w", "STB-8w", "STR-8w")

pdf("Secreted_Protien_N325_Psg22Drop_SelN301_scRNA_Heatmap_N301_selGenes38_8w_STB1_CTB1_EVT1_EVT21_June_2020_new.pdf", width=3, height=6)
col_fn2 <- colorRamp2(c(-4,-2, 0, 2,4),  c( "darkblue", "blue", "white", "yellow", "yellow2"))

Heatmap(as.matrix(scheat.submat1[,c(4,2,3,1)]), col = col_fn2, 
        cluster_columns = F, cluster_rows = TRUE, row_title=NULL, 
        show_row_dend = F,    show_row_names = TRUE, 
        show_column_names = T, row_dend_reorder = TRUE,
        row_km = 3, row_names_gp = gpar(col = rep("black",3), fontsize = rep(6,3)),
        column_dend_height = unit(2, "cm"), column_names_centered = TRUE, 
        row_dend_width = unit(4, "cm") , row_title_rot = 0, column_names_rot = 0,
        heatmap_legend_param = list(title = "log1 (TPM)", title_gp = gpar(col = "black", fontsize = 5),
                                    legend_height = unit(2, "cm"), title_position = "topleft", 
                                    grid_width = unit(0.3, "cm")),
        column_names_gp=gpar(fontsize=6)) 
dev.off()


tiff("Fig3E.tiff", units="in", width=3, height=6, res=300)
Heatmap(as.matrix(scheat.submat1[,c(4,2,3,1)]), col = col_fn2, 
        cluster_columns = F, cluster_rows = TRUE, row_title=NULL, 
        show_row_dend = F,    show_row_names = TRUE, 
        show_column_names = T, row_dend_reorder = TRUE,
        row_km = 3, row_names_gp = gpar(col = rep("black",3), fontsize = rep(6,3)),
        column_dend_height = unit(2, "cm"), column_names_centered = TRUE, 
        row_dend_width = unit(4, "cm") , row_title_rot = 0, column_names_rot = 0,
        heatmap_legend_param = list(title = "log1 (TPM)", title_gp = gpar(col = "black", fontsize = 5),
                                    legend_height = unit(2, "cm"), title_position = "topleft", 
                                    grid_width = unit(0.3, "cm"),
                                    labels_gp = gpar(fontsize = 4)),
        column_names_gp=gpar(fontsize=6)) 
dev.off()

message("+--- Check Transcript Factor list within scRNASeq with expression for Cell type -------+")

TF.dat <- read_excel("TF_check.xlsx", sheet = 1)

scale_data <- GetAssayData(GSE89497_matrix.su, assay = "RNA", slot = "scale.data")
matrix <- as.matrix(scale_data)
matrix <- as.data.frame(t(matrix))
tmatrix <- t(matrix)

matrix$Split <- selected_metadata$Type[match(rownames(selected_metadata), rownames(matrix))]
split.List <- unique(as.character(matrix$Split))
testList <- lapply(split.List, function(x) matrix[matrix$Split==x,-ncol(matrix)])
testList.meanmat <- lapply(testList, function(x) colMeans(apply(x, 2, as.numeric)))

scheat.mat <- cbind(testList.meanmat[[1]], testList.meanmat[[2]], testList.meanmat[[3]], testList.meanmat[[4]],
                    testList.meanmat[[5]])
colnames(scheat.mat) <- c("HE24W_EVT", "HE8W_CTB",  "HE8W_EVT",  "HE8W_STB",  "HE8W_STR")
colnames(scheat.mat) <- c("EVT-24w", "CTB-8w",  "EVT-8w",  "STB-8w",  "STR-8w")

TF.list1oversc <- scheat.mat[rownames(scheat.mat)%in%(as.matrix(TF.dat[1:23,1]))==T,] ## 23
TF.list2oversc <- scheat.mat[rownames(scheat.mat)%in%(as.matrix(TF.dat[1:77,2]))==T,] ## 74
TF.list3oversc <- scheat.mat[rownames(scheat.mat)%in%(as.matrix(TF.dat[1:70,3]))==T,] ## 68

write.csv(TF.list1oversc, file = "TF_list1_oversc_N23.csv", row.names=T, quote=T)
write.csv(TF.list2oversc, file = "TF_list2_oversc_N74_Drop3_EVX1_FOXA2_HNF1A.csv", row.names=T, quote=T)
write.csv(TF.list3oversc, file = "TF_list3_oversc_N74_Drop2_FOXA2_HNF1A.csv", row.names=T, quote=T)

message("+----------Find all markers for seurat object and check whether the TF's are in specific group of cell types----+")

GSE89497_matrix.su <- RunUMAP(GSE89497_matrix.su, dims = 1:20)
pdf("Dimplot_scRNA.pdf", width=10, height=4)
plot1<- DimPlot(GSE89497_matrix.su, reduction = "umap")
plot2 <- DimPlot(GSE89497_matrix.su, reduction = "umap", group.by="Celltype")
plot1+plot2
dev.off()
sig.markers <- FindAllMarkers(GSE89497_matrix.su, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
test1 <- sig.markers[sig.markers$gene%in%(as.matrix(TF.dat[1:23,1]))==T,]
test2 <- sig.markers[sig.markers$gene%in%(as.matrix(TF.dat[1:77,2]))==T,]
test3 <- sig.markers[sig.markers$gene%in%(as.matrix(TF.dat[1:70,2]))==T,]

test1.sub <- subset(test1, cluster!="Macro_1"&cluster!="Blood_cell"&cluster!="Mes_1"&cluster!="Mes_2" & cluster!="Macro_2")
test2.sub <- subset(test2, cluster!="Macro_1"&cluster!="Blood_cell"&cluster!="Mes_1"&cluster!="Mes_2" & cluster!="Macro_2")
test3.sub <- subset(test3, cluster!="Macro_1"&cluster!="Blood_cell"&cluster!="Mes_1"&cluster!="Mes_2" & cluster!="Macro_2")

write.csv(test1.sub, file = "TF_list1_oversc_sigMarker_P23_N14.csv", row.names=T, quote=T)
write.csv(test2.sub, file = "TF_list2_oversc_sigMarker_P77_N43.csv", row.names=T, quote=T)
write.csv(test3.sub, file = "TF_list3_oversc_sigMarker_P70_N40.csv", row.names=T, quote=T)

message("+----                    Tissue enrichment Analysis               ---------------+")

library("TissueEnrich")

genes1 <- c("TFPI2", "SERPINE2", "PAPPA2", "IGF2", "FLT1")
genes2 <- c("TNFRSF9", "TNFRSF11B", "TFPI2", "TFPI", "SERPINE1", "SERPINB9G", "SERPINB9E", "SERPINB9B",
            "SCT", "SBSN", "PGA5", "PAPPA2", "LAMA1", "INHBA", "GZMF", "GZMC", "GKN1", "FLT1", "FBLN7",
            "CTSK", "CSF1R", "CREG1", "CEACAM5", "C1QTNF1", "A2M")
genes3 <- rownames(scheat.mat)

###### Heatmap for enrich human placenta Genes 4 plot

gs1 <- GeneSet(geneIds=genes1, organism="Homo Sapiens", geneIdType=SymbolIdentifier())
output1 <- teEnrichment(inputGenes = gs1, rnaSeqDataset = 1)
seEnrichmentOutput1 <- output1[[1]]
enrichmentOutput1 <- setNames(data.frame(assay(seEnrichmentOutput1),row.names = rowData(seEnrichmentOutput1)[,1]), colData(seEnrichmentOutput1)[,1])
enrichmentOutput1$Tissue <- row.names(enrichmentOutput1)


library(tidyr)

seExp1 <- output1[[2]][["Placenta"]]
exp1 <- setNames(data.frame(assay(seExp1), row.names = rowData(seExp1)[,1]), colData(seExp1)[,1])
exp1$Gene <- row.names(exp1)
exp1 <- exp1%>% gather(key = "Tissue", value = "expression",1:(ncol(exp1)-1))


pdf("TissueEnrich_Heatmap_G5.pdf", width=4, height=1)
ggplot(exp1, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                            colour = "white") + scale_fill_gradient(low = "grey87",
                                                                                    high = "red2")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title = element_text(size=2),
        legend.title=element_text(size=4), legend.text=element_text(size=3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size= 3),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank(),
        axis.text.y = element_text(vjust = 1, hjust = 1, size= 4),
        legend.key.size = unit(0.5,"line"),
        panel.border = element_blank(),
        axis.ticks = element_blank()) 
dev.off()

###### Heatmap for enrich mouse placenta Genes 25 plot

gs2 <- GeneSet(geneIds=genes2, organism="Mus Musculus", geneIdType=SymbolIdentifier())
output2 <- teEnrichment(inputGenes = gs2, rnaSeqDataset = 3)
seEnrichmentOutput2 <- output2[[1]]
enrichmentOutput2 <- setNames(data.frame(assay(seEnrichmentOutput2),row.names = rowData(seEnrichmentOutput2)[,1]), colData(seEnrichmentOutput2)[,1])
enrichmentOutput2$Tissue <- row.names(enrichmentOutput2)


seExp2 <- output2[[2]][["E14.5-Placenta"]]
exp2 <- setNames(data.frame(assay(seExp2), row.names = rowData(seExp2)[,1]), colData(seExp2)[,1])
exp2$Gene <- row.names(exp2)
exp2 <- exp2%>% gather(key = "Tissue", value = "expression",1:(ncol(exp2)-1))


pdf("TissueEnrich_Heatmap_G25.pdf", width=4, height=3)
ggplot(exp2, aes(Tissue, Gene)) + geom_tile(aes(fill = expression),
                                            colour = "white") + scale_fill_gradient(low = "grey87",
                                                                                    high = "red2")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 5),axis.title = element_text(size=2),
        legend.title=element_text(size=5), legend.text=element_text(size=5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size= 4),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank(),
        axis.text.y = element_text(vjust = 1, hjust = 1, size= 4),
        legend.key.size = unit(0.5,"line"),
        panel.border = element_blank(),
        axis.ticks = element_blank()) 
dev.off()

message("+---- Featureplot for selected proteins in TF and mouse unique --------+")

markers1 <- c("FLT1", "TFPI2", "ANGPT2")
markers2 <- c("ARNT2", "ELF3", "PLAG1", "SP2", "MEF2D", "MYCN", "FOS", "NFYC", "CREB1", "IRF3")
markers <- c(markers1, markers2)
pdf("Dimplot_scRNA_Liu.pdf")

DimPlot(GSE89497_matrix.su, group.by="Celltype", pt.size = 0.8, label=T)

dev.off()


plt <- list()
for( i in 1:13){
  plt[[i]] <- FeaturePlot(GSE89497_matrix.su, features=markers[i], pt.size = 0.5)
  plt
}

pdf("FeaturePlot_scRNA_UniMarker_Liu.pdf", width= 9, height= 2)
plot_grid(plt[[1]], plt[[2]], plt[[3]], nrow = 1, ncol = 3)
dev.off()

pdf("FeaturePlot_scRNA_TFMarker_Liu.pdf", width= 8, height= 15 )
plot_grid(plt[[4]], plt[[5]], plt[[6]], plt[[7]], plt[[8]], plt[[9]],plt[[10]],
          plt[[11]], plt[[12]],plt[[13]], nrow = 5, ncol=2)
dev.off()



message("+-----------              END of Script        ------------------------+")
