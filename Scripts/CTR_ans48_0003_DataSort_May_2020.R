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

Protein.dat2 <- read_excel(paste0(Data.dir, "/Sorted_Trophblast_Cells_Data.xlsx"))
colnames(Protein.dat2) <- c("Accession.Number", "S3", "S4", "S5", "P3", "P4", "P5")


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
QC.mat <- Protein.dat2c[,2:4]+Protein.dat2c[,5:7]
QC.mat.ind <- ifelse(QC.mat==0,0,1)
QC.mat.ind1 <- which(rowSums(QC.mat.ind)>=2)

Protein.dat2c.new <- Protein.dat2c[QC.mat.ind1,]
Protein.dat2c.new$Accession.New <- Protein.dat2c.new$Accession.Number
## Go to the Uniprot revive website to convert the entryname to Uniprot ID and GeneName
## Connect with the MouseGeneName using UniProt.ws library.
Protein.dat2.uni <- read.csv(paste0(Data.dir, "/Sorted_Trophoblast_UniProt_N682_filtered.csv"), header=T)[,1:4]

Protein.dat2.final <- merge(Protein.dat2c.new, Protein.dat2.uni, by = "Accession.New")
Protein.dat2.final$MouseGeneName <- as.character(Protein.dat2.final$MouseGeneName)
Protein.dat2.final[which(Protein.dat2.final$Accession.New=="H2B1K_MOUSE"),11] <- "Hist1h2bk"
colnames(Protein.dat2.final) <- c("Accession.New", "Accession.Number.ori", "S3", "S4", "S5", "P3", "P4", "P5", 
                                  "Accession.Number.uni", "ProteinID", "MouseGeneName")

write.csv(Protein.dat2.final, file = paste0(Data.dir, "/Sorted_Trophoblast_N682_G681_Naca_Filtered_Step2_Data.csv"),row.names =F)

test2 <- Protein.dat2.final[Protein.dat2.final$MouseGeneName%in%ensEMBL2id$external_gene_name==F,]
write.csv(test2, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N41_G41_NotEnsem_Filtered_Step2_Data.csv"), row.names=F, quote=T)


message("+----                 Perform Overlap with the GEO public Data           --------+")

MousePublic <- read.csv(paste0(Data.dir, "/GEO_Control_Mouse_D3_N47747_GeneName_List.csv"), header=T)
colnames(MousePublic) <- "MouseGeneName"

Dat1.overlap1 <- Protein.dat1.final[Protein.dat1.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 908-907--906, Psg16, Pkm
Dat2.overlap1 <- Protein.dat2.final[Protein.dat2.final$MouseGeneName%in%MousePublic[,1]==T,] 
## 642-642--641, Naca
Dat3.overlap1 <- PData3.final[PData3.final$MouseGeneName%in%MousePublic[,1]==T,]
## 1180-1180-1178, Gnas, Naca, 28 genes are not found in ensEMBLE and also not found in public Mouse.


message("+----                 Perform Homolog overlap           --------+")

HumanPublic <- read.csv(paste0(Data.dir, "/GEO_Control_Human_D8_N34673_GeneName_List.csv"), header=T)
HomoData <- read.csv(paste0(Data.dir, "/", Project, "-Ensemble_MGI_NCBI_human_mouse_Homolog_april_2020.csv"), header=T)

HomoData.overPub <- HomoData[as.character(HomoData$HumanGeneName)%in%as.character(HumanPublic[,1])==T,]
## make HomoData.overPub unique for mouse and Human matching and remove all of the 1-more or more-1 genes
Mmulit <- HomoData.overPub[duplicated(HomoData.overPub$MouseGeneName)==T,]
## 5362
Hmulit <- HomoData.overPub[duplicated(HomoData.overPub$HumanGeneName)==T,]
## 6679

HomoData.overPub.Muni <- unique(as.character(HomoData.overPub$MouseGeneName)) 
colnames(HomoData.overPub.Muni) <- "MouseGeneName"

Dat1.overlap2 <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 898-898--896, Psg16, Pkm
Dat2.overlap2 <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 638-638--637, Naca
Dat3.overlap2 <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%HomoData.overPub.Muni==T,] 
## 1175-1175-1173, Gnas, Naca

message("+----  Suggestions from Tina to remove some genes, due to multiple homologus with human       ----+")

Filter.proteins <- read.csv(paste0(Data.dir, "/Filter_Proteins_List.csv"), header=T)

Dat1.overlap2.new <- Dat1.overlap2[Dat1.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 881-881--879, Psg16, Pkm
Dat2.overlap2.new <- Dat2.overlap2[Dat2.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 634-634--633, Naca
Dat3.overlap2.new <- Dat3.overlap2[Dat3.overlap2$MouseGeneName%in%Filter.proteins$MouseGeneName==F,] 
## 1172-1172-1170, Gnas, Naca

test1.mh.over <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==T,]
test1.mh.nover <- Dat1.overlap1[Dat1.overlap1$MouseGeneName%in%Dat1.overlap2.new$MouseGeneName==F,]

test2.mh.over <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==T,]
test2.mh.nover <- Dat2.overlap1[Dat2.overlap1$MouseGeneName%in%Dat2.overlap2.new$MouseGeneName==F,]

test3.mh.over <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==T,]
test3.mh.nover <- Dat3.overlap1[Dat3.overlap1$MouseGeneName%in%Dat3.overlap2.new$MouseGeneName==F,]

message("+---write out each step proteins list for flowchart Data1, 2,3 ----------------------+")

write.csv(PData3.sel.QC, file = paste0(Data.dir, "/Cultured_Trophoblast_N1534_G1531_Gnas_Naca_Gpx4_Filtered_pep2_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap1, file = paste0(Data.dir, "/Cultured_Trophoblast_N1180_G1178_Gnas_Naca_Filtered_pep2_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat3.overlap2.new, file = paste0(Data.dir, "/Cultured_Trophoblast_N1172_G1170_Gnas_Naca_Filtered_pep2_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(test3.mh.nover,file = paste0(Data.dir, "/Cultured_Trophoblast_N8_G8_Filtered_pep2_Step4_MusUni_Data.csv"), row.names=F, quote=T )

write.csv(Protein.dat2c, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N1142_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.overlap1, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N642_G641_Naca_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(Dat2.overlap2.new, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N634_G633_Naca_Filtered_Step3_HumPub_Data.csv"), row.names=F, quote=T)
write.csv(test2.mh.nover, file = paste0(Data.dir, "/Sorted_cell_Trophoblast_N8_G8_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)


write.csv(Protein.dat1c, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N1445_G1441_Psg16_Pkm_Tpm1_Fbln1_Step1_Data.csv"), row.names=F, quote=T)
write.csv(Dat1.overlap1, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N881_G879_Psg16_Pkm_Filtered_Step3_MousePub_Data.csv"), row.names=F, quote=T)
write.csv(test1.mh.nover, file = paste0(Data.dir, "/Conditioned_Medium_Trophoblast_N27_G27_Filtered_Step4_MusUni_Data.csv"), row.names=F, quote=T)


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

pdf("Secreted_N325_Psg22Drop_SelN301_scRNA_Heatmap_June_2020_addEVT24.pdf", width = 3, height = 40)
col_fn2 <- colorRamp2(c(-4,-2, 0, 2,4),  c("darkorchid4", "darkviolet", "white", "darkolivegreen", "darkgreen"))
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
col_fn2 <- colorRamp2(c(-4,-2, 0, 2,4),  c("darkorchid4", "darkviolet", "white", "darkolivegreen", "darkgreen"))

Heatmap(as.matrix(scheat.submat1[,c(4,2,3,1)]), col = col_fn2, 
        cluster_columns = F, cluster_rows = TRUE, row_title=NULL, 
        show_row_dend = F,    show_row_names = TRUE, 
        show_column_names = T, row_dend_reorder = TRUE,
        row_km = 3, row_names_gp = gpar(col = rep("black",3), fontsize = rep(6,3)),
        column_dend_height = unit(2, "cm"), column_names_centered = TRUE, 
        row_dend_width = unit(4, "cm") , row_title_rot = 0, column_names_rot = 45,
        heatmap_legend_param = list(title = "Expression", title_gp = gpar(col = "black", fontsize = 7),
                                    legend_height = unit(2, "cm"), title_position = "topleft", 
                                    grid_width = unit(0.2, "cm")) ) 

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

message("+-----------              END of Script        ------------------------+")
