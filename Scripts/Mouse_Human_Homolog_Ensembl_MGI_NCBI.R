## !/usr/local/bin/Rscript
##
## Mouse Human Homolog list from Ensembl, MGI and NCBI
##
## CopyRight: Xiaohui Zhao (xz289@cam.ac.uk)
##


message("+---------------------------------------------------------------------------------------+")
message("+                   Install useful packages and setup working directory                 +")
message("+---------------------------------------------------------------------------------------+")

library("biomaRt")
library("GEOquery")
library("UniProt.ws")
library("rgdal")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("readxl")

Base.dir <- "/Users/xz289/Documents/CTR_ans48_0003"
Data.dir <- paste0(Base.dir, "/P147_1_20170707")
Project <- "CTR_ans48_0003"
Out.dir <- paste0(Base.dir, "/Tina_Paper_Data")


message("+---------------------------------------------------------------------------------------+")
message("+   get list of mouse/human genes, getLDS homolog in ensemble, MGI and NCBI             +")
message("+---------------------------------------------------------------------------------------+")
## Ensemble getLDS homolog list calling
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")


ensemblM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", 
                                 "uniprot_gn_symbol", "uniprot_gn_id"),
                  mart = mouse)
ensemblH <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", 
                                 "uniprot_gn_symbol", "uniprot_gn_id"),
                  mart = human)

ensemblMG <- getBM(attributes = c("external_gene_name"), mart = mouse)
ensemblHG <- getBM(attributes = c("external_gene_name"), mart = human)

Ensembl_MHhomolog <- getLDS(attributes = c("external_gene_name"),
                            filters = "external_gene_name", values= ensemblMG$external_gene_name,
                            mart = mouse,
                            attributesL = c("external_gene_name"), 
                            martL = human) 

## 26640 unique mouse genes is 20682; unique human genes is 19673
colnames(Ensembl_MHhomolog) <- c("MouseGeneName", "HumanGeneName")

write.csv(Ensembl_MHhomolog, file = paste0(Out.dir, "/", Project, "-Ensemble_human_mouse_Homolog_april_2020.csv"), row.names =F)

Ensembl_MHhomolog_test <- getLDS(attributes = c("external_gene_name", "uniprot_gn_id"),
                                 filters = "external_gene_name", values= ensemblM$external_gene_name,
                                 mart = mouse,
                                 attributesL = c("external_gene_name","uniprot_gn_id"), 
                                 martL = human) 

message("+---------------------------------------------------------------------------------------+")
message("+ MGI homolog list calling:                                                             +")
message("+ Complete List of Human and Mouse Homologs with phenotype annotations (plain text)     +")
message("+---------------------------------------------------------------------------------------+")

MGI.MHhomolog.all <- read.xlsx(paste0(Out.dir, "/", Project, "-MGI_human_mouse_Homolog_Phenotype.xlsx"), 
                               sheetIndex=1) 
## 18748 in total, no 678; yes 18070
MGI_MHhomolog <- subset(MGI.MHhomolog.all, MGI.MHhomolog.all$Homo=="yes")[,c(1,5)]

#### Test uniquely from MGI and Ensemble homolog genes list unique and common
testMGI.ensembl <- merge(Ensembl_MHhomolog, MGI_MHhomolog, by = "MouseGeneName")
length(unique(testMGI.ensembl[,1]))
## common mouse genes number is 16977
MGI.ensembl.FalseM <- MGI_MHhomolog[MGI_MHhomolog$MouseGeneName%in%Ensembl_MHhomolog$MouseGeneName==F,]
length(unique(MGI.ensembl.FalseM[,1]))
## 866
ensembl.MGI.FalseM <- Ensembl_MHhomolog[Ensembl_MHhomolog$MouseGeneName%in%MGI_MHhomolog$MouseGeneName==F,]
length(unique(ensembl.MGI.FalseM[,1]))
## 3705

message("+---------------------------------------------------------------------------------------+")
message("+ NCBI homolog list calling:                                                            +")                 
message("+---------------------------------------------------------------------------------------+")

NCBI_Human <- read.xlsx(paste0(Out.dir, "/", Project, "-NCBI_human_mouse_Homolog.xlsx"), sheetIndex=2)
NCBI_Mouse <- read.xlsx(paste0(Out.dir, "/", Project, "-NCBI_human_mouse_Homolog.xlsx"), sheetIndex=1)
NCBI_mer <- merge(NCBI_Mouse, NCBI_Human, by = "ID")

NCBI_MHhomolog <- NCBI_mer[,-1]
## 17355

message("+---------------------------------------------------------------------------------------+")
message("+            UpsetR plot to show the overlap and unique of three Lists homolog          +")               
message("+---------------------------------------------------------------------------------------+")

Ensembl_MGI_NCBI <- unique(c(as.character(Ensembl_MHhomolog$MouseGeneName), 
                             as.character(MGI_MHhomolog$MouseGeneName),
                             as.character(NCBI_MHhomolog$MouseGeneName))) 
## 22791
EMN_mat <- matrix(0, ncol=3, nrow = length(Ensembl_MGI_NCBI))
EMN_mat[,1] <- ifelse(Ensembl_MGI_NCBI%in%Ensembl_MHhomolog$MouseGeneName==T, 1, 0)
EMN_mat[,2] <- ifelse(Ensembl_MGI_NCBI%in%MGI_MHhomolog$MouseGeneName==T, 1, 0)
EMN_mat[,3] <- ifelse(Ensembl_MGI_NCBI%in%NCBI_MHhomolog$MouseGeneName==T, 1, 0)
rownames(EMN_mat) <- Ensembl_MGI_NCBI
colnames(EMN_mat) <- c("Ensembl", "MGI", "NCBI")
EMN_matplot <- as.data.frame(EMN_mat)

library(UpSetR)
pdf(paste0(Out.dir, "/", Project, "-Ensembl_MGI_NCBI_human_mouse_upsetR_plot_april_2020.pdf"))
upset(EMN_matplot, sets = c("Ensembl", "MGI", "NCBI"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "off")
dev.off()


message("+---------------------------------------------------------------------------------------+")
message("+            get final list of homolog genes/proteins for mouse and human               +")
message("+---------------------------------------------------------------------------------------+")

MGI_MHhomolog <- MGI_MHhomolog[,c(2,1)]
Final_MHhomolog <- unique(rbind(Ensembl_MHhomolog, MGI_MHhomolog, NCBI_MHhomolog)) 
## 29254, 22791unique mouse; 21549 unique human.
write.csv(Final_MHhomolog, paste0(Out.dir, "/", Project, "-Ensemble_MGI_NCBI_human_mouse_Homolog_april_2020.csv"), row.names =F)

message("+------------------------------------FIN-----------------------------------------------+")
