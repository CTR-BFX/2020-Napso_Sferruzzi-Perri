## !/usr/local/bin/R
##
## Rscript for downloading and sort the public Mouse and Human data set
##
## Copyright: Xiaohui Zhao (xz289@cam.ac.uk)
##

message("+-------------------------------------------------------------------------------+")
message("+    Install the packages and basic settings                                    +")
message("+-------------------------------------------------------------------------------+")

library("GEOquery")
library("biomaRt")
library("preprocessCore") 
## normalise.quantile
Base.dir <- "/Users/xz289/Documents/CTR_ans48_0003"
data.dir <- "/Users/xz289/Documents/CTR_ans48_0003/Temp_GEOData"
out.dir <- "/Users/xz289/Documents/CTR_ans48_0003/Original_Data"
setwd(data.dir)

message("+-------------------------------------------------------------------------------+")
message("+      L2_Human_noraml_2ndtrimester_Affymetrix_GSE9984(GSM252353-252356)        +")
message("+ 2 columns, affy id and RMA normalised signal log2 value ind. sample           +")
message("+                 RNA-seq Data, 2nd trimester is in cols 6-9                    +")
message("+-------------------------------------------------------------------------------+")
#
# Call the human ensemble set with affymetrix id to match the RNA-seq matrix
#
ensemblH    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
# Use [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array 
ensEMBL2affidH <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 
                                     'entrezgene_id', 'affy_hg_u133_plus_2'), mart = ensemblH)
ensEMBL2affidH.new <- ensEMBL2affidH[-which(ensEMBL2affidH$affy_hg_u133_plus_2==""),]
colnames(ensEMBL2affidH.new) <- c("EnsembID.Human", "GeneName.Human", "EntrezGene.Human", "AffyID.Human")


DatH_normal_Aff <- getGEO("GSE9984",  destdir = data.dir,
                          GSEMatrix = TRUE, AnnotGPL = T, getGPL = TRUE)

Human_2ndTrim_RNAmat <- matrix(scan(gzfile("GSE9984_series_matrix.txt.gz"), 
                                    what="character", skip = 62, nlines = 54676), ncol = 13, 
                               byrow = T)[,c(1,5:8)]
colnames(Human_2ndTrim_RNAmat) <- c("AffyID.Human", Human_2ndTrim_RNAmat[1,c(2:5)])
Human_Dat1 <- as.data.frame(Human_2ndTrim_RNAmat[-1,])

# normalised by quantile normalisation
dat1 <- normalize.quantiles((apply(Human_Dat1[,-1],2, as.numeric)))
ndat1 <- data.frame(Human_Dat1[,1], dat1)
colnames(ndat1) <- colnames(Human_Dat1)
Merge.dat11 <- merge(ndat1, ensEMBL2affidH.new, by = "AffyID.Human")

# dim 2881, 918 unique human genes, 922 unique mouse genes (3 genes missed)

# This data set we lose Nccrp1, Hist1h4i, Hnrnpa1 mouse genes
# simplified the data with GeneName.Mouse, UniProtID.Mouse, GeneName.Human, UniProtID.Human, Mean value of expression
Hdat1 <- Merge.dat11

write.table(Hdat1, file = paste0(Base.dir, "/Level1_output/Human_2ndTrim_RNA_GSE9984_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)


message("+-------------------------------------------------------------------------------+")
message("+      L2_Human_noraml_thirdtrimester_Microarrays_GSE28551(GSM707067-707087)    +")
message("+            Placenta normal third trimester rep1.- 21.                         +")
message("+      L2_Human_noraml_firsttrimester_Microarrays_GSE28551(GSM707051-707066)    +")
message("+            Placenta normal third trimester rep1.- 16.                         +")
message("+           Table with ID_REF and VALUE (Normalized signal intensity)           +")
message("+           ABI Human Genome Survey Micoarray version 2, GPL2986                +")
message("+-------------------------------------------------------------------------------+")

DatH_microarray_trimester <- getGEO("GSE28551",  destdir = data.dir,
                                    GSEMatrix = TRUE, AnnotGPL = T, getGPL = TRUE)
Human_3rdTrim_Micomat <- matrix(scan(gzfile("GSE28551_series_matrix.txt.gz"), 
                                     what="character", skip = 66, nlines = 17045), ncol = 38, 
                                byrow = T)[,c(1, 18:38)]
colnames(Human_3rdTrim_Micomat) <- c("ABI.ID", Human_3rdTrim_Micomat[1,2:22])
Human_DatM <- Human_3rdTrim_Micomat[-1,]
#
# Log2 then qunatile normalized
#
Human_DatM_Norm <- normalize.quantiles(as.matrix(apply(Human_DatM[,-1], 2, as.numeric)))
GSE28551_Normal <- apply(Human_DatM_Norm, 1, mean)
GPL2986 <- matrix(scan("GPL2986_annot.txt", what="character", fill = T, sep="\t"), ncol = 3, byrow=T)
GPL2986.annot <- GPL2986[-which(GPL2986[,2]==""),]
colnames(GPL2986.annot) <- c("ABI.ID", "GeneName.Human", "EntrezGene.Human")
# 32878---17916
# 
# Merge GPL and Matrix
Human_DatM_new <- data.frame(Human_DatM[,1], GSE28551_Normal)
colnames(Human_DatM_new) <- c("ABI.ID", "GSE28551_Normal")
Merge.dat21 <- merge(Human_DatM_new, GPL2986.annot, by ="ABI.ID")
Hdat2 <- Merge.dat21

write.table(Hdat2, file = paste0(Base.dir, "/Level1_output/Human_3rdTrim_Array_GSE28551_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)



message("+--------------------------------------------------------------------------------------+")
message("+Mouse Placenta_micorarry E15 is chosen, all cell type of the placenta (fetal origin)  +")
message("+      Placenta/Decidua (e8.5,9.0,10.5,12.0,13.5,15.0,17.0,19.0; p0)                   +")
message("+ This superseries is composed of the following subseries, GSE11220 and GSE11222       +")
message(" GSE11220 Timecourse of developing mouse placenta with placental and decidual   
        tissues profiled separately, series 1                                                   +")
message(" GSE11222 Placental and decidual timecourse samples normalized and modeled with
        an undissected e17 sample, series 2                                                     +")
message("+Placenta RNA-seq series 1 with 3 replicates,GSM282776-282778,series 2 GSM282825-282827+")
message("+Decidue RNA-seq series 1 with 2 replicates,GSM282796-282797,series 2 GSM282845-282846 +")
message("+        GPL1261 [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array                    +")
message("+      https://github.com/uiuc-cgm/placenta-viviparous.git                             +")
message("+dCHIP absolute expression value.                                                      +")
message("+--------------------------------------------------------------------------------------+")

DatM_RNA_E15 <- getGEO("GSE11224",  destdir = data.dir,
                       GSEMatrix = TRUE, AnnotGPL = T, getGPL = TRUE)
Mouse_E15_RNAmat <- matrix(scan(gzfile("GSE11224_series_matrix.txt.gz"), 
                                what="character", skip = 63, nlines = 45102), ncol = 87, 
                           byrow = T)[,c(1,18:20, 58:60, 38:39,78:79)]

Mouse_RNADat_E15 <- Mouse_E15_RNAmat[-1,]
colnames(Mouse_RNADat_E15) <- c("Affy.ID", Mouse_E15_RNAmat[1,2:11])
Mpla1 <- apply(Mouse_RNADat_E15[,c(2:7)], 2, as.numeric)
Mpla2 <- apply(Mpla1, 1, mean)
GSE11224_P <- normalize.quantiles(as.matrix(log2(Mpla2)))

Mdec1 <- apply(Mouse_RNADat_E15[,c(8:11)], 2, as.numeric)
Mdec2 <- apply(Mdec1, 1, mean)
GSE11224_D <- normalize.quantiles(as.matrix(log2(Mdec2)))

MouseP <- data.frame(Mouse_RNADat_E15[,1], GSE11224_P)
MouseD <- data.frame(Mouse_RNADat_E15[,1], GSE11224_D)
colnames(MouseP) <- c("Affy.ID", "GSE11224_P")
colnames(MouseD) <- c("Affy.ID", "GSE11224_D")
#
# The GPL similar as the Human ones, not proper formatted, will use ensembl affy_mouse430_2
#
ensemblM    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2affidM <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 
                                     'entrezgene_id', 'affy_mouse430_2'), mart = ensemblM)
ensEMBL2affidM.new <- ensEMBL2affidM[-which(ensEMBL2affidM$affy_mouse430_2==""),]
colnames(ensEMBL2affidM.new) <- c("EnsembID.Mouse", "GeneName.Mouse", "EntrezGene.Mouse", "Affy.ID")
# 83563----40931
#
# Merge three data sets, Tina's, ensemble affymetrix and public data to get common genes mat.
Merge.ME <- merge(MouseP, ensEMBL2affidM.new, by = "Affy.ID")
MDat1 <- Merge.ME

Merge.MP <- merge(MouseD, ensEMBL2affidM.new, by = "Affy.ID")
MDat2 <- Merge.MP

write.table(MDat1, file = paste0(Base.dir, "/Level1_output/Mouse_E15_placenta_decidue_RNA_GSE11224_Placenta_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)
write.table(MDat2, file = paste0(Base.dir, "/Level1_output/Mouse_E15_placenta_decidue_RNA_GSE11224_Decidue_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)

message("+-------------------------------------------------------------------------------+")
message("+      Species8_placentas_RNASeq_GSE79121(GSM2086260 Mouse)                     +")
message("+      https://github.com/uiuc-cgm/placenta-viviparous.git                      +")
message("+ Dasypus novemcinctus; Ateles fusciceps; Pan paniscus; Canis lupus familiaris; 
        Loxodonta africana; Bos taurus; Mus musculus; Monodelphis domestica              +")
message("+                  GPL9250, fpkm normalisation                                  +")
message("+ Expression Profiling of Term Placenta in Viviparous Mammals by RNA-Seq        +")
message("+   SRX1629206(SRR3222431), GSM2086260_mus_musculus_genes.fpkm_tracking.gz      +")
message("+   GSE66622 human placental villus parenchyma, 	GPL10999")
message("+-------------------------------------------------------------------------------+")

S8_mouse_RNA_FPKM <- matrix(scan(gzfile("GSM2086260_mus_musculus_genes.fpkm_tracking.gz"),
                                 what = "character", fill = T, sep = "\t"), ncol = 13, byrow=T) 
S8_mouse_RNA_FPKM_Dat <- S8_mouse_RNA_FPKM[,c(4, 5, 10)]
S8mouse_RNAFPKM_Dat <- S8_mouse_RNA_FPKM_Dat[-1,]
colnames(S8mouse_RNAFPKM_Dat) <- c("EnsembID.Mouse", "GeneName.Mouse", "FPKM")
GSE79121_Normal <- normalize.quantiles(as.matrix(log2(as.numeric(S8mouse_RNAFPKM_Dat[,3])+2)))
# 45390 in total, 36281 with 0 FPKM, 9109 with values.
#
# merge mouse RNA FPKM with Tina_Data
MDat3 <- data.frame(S8mouse_RNAFPKM_Dat[,1:2], GSE79121_Normal)
colnames(MDat3) <- c("EnsembID.Mouse", "GeneName.Mouse", "GSE79121_Normal")


# missed Aprt,Galns,Ftl1-ps1,Plpbp,Elob,Rack1,Mesd,Rab1a,Ptpa,Nectin2 
write.table(MDat3, file = paste0(Base.dir, "/Level1_output/Mouse_S8RNA_FPKM_Term_placenta_allcellT_GSE79121_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)



message("+----------------------------------------------------------------------------------------------------+")
message("+ Compared 7 GEO datasets for normal vs PE pregnancy--- we use the normal                            +")
message("+ L2_Human_placentas_26normal_pregnancy_GSE10588(GSM225470-225480,225483-225497) Array               +")
message("+ GPL2986 ABI Human Genome Survey Microarray Version 2                                               +")
message("+ Quantile normalized value of raw signal count. The raw signal value is the corrected,              +") 
message("+ background subtracted measurement of chemiluminescent signal, as outputted by the ABI 1700 software+")
message("+----------------------------------------------------------------------------------------------------+")

for( name in c("GSE10588", "GSE43942", "GSE4707", "GSE25906","GSE24129", "GSE30186", "GSE44711")){
  getGEO(name, destdir = data.dir,
         GSEMatrix = TRUE, AnnotGPL = T, getGPL = TRUE)
}


# call the data and match the Annotation

L2_mata <- matrix(scan(gzfile("GSE10588_series_matrix.txt.gz"), skip = 71, what = "character",
                       nlines = 32879, sep="\t", fill = T), ncol = 44, byrow = T)[,c(1,2:12, 14:28)]
colnames(L2_mata) <- c("ID", L2_mata[1,2:27])
L2_mata <- L2_mata[-1,]
L2_mat1 <- apply(L2_mata[,-1], 2, as.numeric)
L2_mat2 <- log2(L2_mat1)
GSE10588_Normal <- apply(L2_mat2, 1, mean)
L2_mat3 <- data.frame(L2_mata[,1], GSE10588_Normal)
colnames(L2_mat3) <- c("ID", "GSE10588_Normal")

L2_mata_GPL <- read.table("GPL2986.annot.txt", header = F, fill = T)
L2_mata_GPL.1 <- L2_mata_GPL[-which(L2_mata_GPL[,2]==""),] # 43932----17916
colnames(L2_mata_GPL.1) <- c("ID", "GeneName.Human", "EntrezID.Human")

# merge data with Entrez ID
Merge_L2mata1 <- merge(L2_mat3, L2_mata_GPL.1, by = "ID")

entrezID <-ensEMBL2affidH[,c(1:3)]
colnames(entrezID) <- c("EnsembID.Human", "GeneName.Human", "EntrezID.Human")

Merge_TEL2a <- merge(Merge_L2mata1, entrezID, by = "GeneName.Human")

Hdat3 <- Merge_TEL2a

write.table(Hdat3, file = paste0(Base.dir, "/Level1_output/Human_placentas_26normal_pregnancy_GSE10588_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)


message("+-------------------------------------------------------------------------------------------------+")
message("+ L2_Human_4Normal_GSE4707(GSM106249,106250,106274,106275) Array                                  +")
message("+ GPL1708 Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Feature Number version)      +")
message("+ Log10 Ratio of normalized Cy5 (red) signal to normalized Cy3 (green) signal                     +")
message("+-------------------------------------------------------------------------------------------------+")

L2_matc <- matrix(scan(gzfile("GSE4707_series_matrix.txt.gz"), skip = 62, what = "character",
                       nlines = 43932, sep="\t", fill = T), ncol = 15, byrow = T)[,c(1,2:5)]
colnames(L2_matc) <- c("ID", L2_matc[1,c(2:5)])
L2_matc_new <- L2_matc[-1,]
L2_matc_new1 <- apply(L2_matc_new[,-1], 2, as.numeric)
GSE4707_mean <- apply(L2_matc_new1, 1, mean)
GSE4707_Normal <- normalize.quantiles(as.matrix(log2(10^(GSE4707_mean)+8)))
L2_matc_new2 <- data.frame(L2_matc_new[,1], GSE4707_Normal)
colnames(L2_matc_new2) <- c("ID", "GSE4707_Normal")

L2_matc_GPL <- read.table("GPL1708.annot.txt", header = F, fill = T)
colnames(L2_matc_GPL) <- c("ID", "GeneName.Human", "EntrezID.Human")

# merge Data
Merge_matCGPL <- merge(L2_matc_new2, L2_matc_GPL, by = "ID")
Hdat5 <- Merge_matCGPL

write.table(Hdat5, file = paste0(Base.dir, "/Level1_output/Human_placentas_4normal_pregnancy_GSE4707_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)

message("+------------------------------------------------------------------------------------------------+")
message("+ L2_Human_37Control_GSE25906(GSM635927-635963) Array                                            +")
message("+ GPL6102 Illumina human-6 v2.0 expression beadchip                                              +")
message("+ log2 transformed quantile normalized signal intensity                                          +")
message("+------------------------------------------------------------------------------------------------+")

L2_matd <- matrix(scan(gzfile("GSE25906_series_matrix.txt.gz"), skip = 67, what = "character",
                       nlines = 48702, sep="\t", fill = T), ncol = 61, byrow = T)[,c(1,25:61)]
colnames(L2_matd) <- L2_matd[1,]
L2_matd_new <- L2_matd[-1, ]
L2_matd_new1 <- apply(L2_matd_new[,-1], 2, as.numeric)
GSE25906_Normal <- apply(L2_matd_new1, 1, mean)
L2_matd_new2 <- data.frame(L2_matd_new[,1], GSE25906_Normal)
colnames(L2_matd_new2) <- c("ID_REF", "GSE25906_Normal")

L2_matd_GPL <- read.table("GPL6102.annot.txt", header = T, fill = T)
colnames(L2_matd_GPL) <- c("ID_REF", "GeneName.Human", "EntrezID.Human")

Merged1 <- merge(L2_matd_new2, L2_matd_GPL, by = "ID_REF")[,c(1:6)]

Hdat4 <- Merged1

write.table(Hdat4, file = paste0(Base.dir, "/Level1_output/Human_placentas_37normal_pregnancy_GSE25906_april_new.txt"), 
            row.names=F, col.names = T, quote = F)

message("+-------------------------------------------------------------------------------------------------+")
message("+ L2_Human_8Control_GSE24129(GSM594041-594048)                                                    +")
message("+ GPL6244 [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]          +")
message("+ Quantile normalized gene level expression values from GeneBASE                                  +")
message("+-------------------------------------------------------------------------------------------------+")

L2_mate <- matrix(scan(gzfile("GSE24129_series_matrix.txt.gz"), skip = 66, what = "character",
                       nlines = 24806, sep="\t", fill = T), ncol = 25, byrow = T)[,c(1:9)]
colnames(L2_mate) <- c("affy_hugene_1_0_st_v1", L2_mate[1,c(2:9)])
L2_mate1 <- L2_mate[-1,]
L2_mate2 <- apply(L2_mate1[,-1], 2, as.numeric)
L2_mate3 <- apply(L2_mate2, 1, mean)
GSE24129_Normal <- log2(L2_mate3)
L2_mate4 <- data.frame(L2_mate1[,1], GSE24129_Normal)
colnames(L2_mate4) <- c("affy_hugene_1_0_st_v1", "GSE24129_Normal")
# L2_mate_GPL <- read.table("GPL6244.annot.txt", header = F, fill =T)
ens_mate <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 
                               'entrezgene_id', 'affy_hugene_1_0_st_v1'), mart = ensemblH)
colnames(ens_mate) <- c("EnsembID.Human", "GeneName.Human", "EntrezID.Human", "affy_hugene_1_0_st_v1")
Mergee1 <- merge(L2_mate4, ens_mate, by = "affy_hugene_1_0_st_v1")
Hdat7 <- Mergee1

write.table(Hdat7, file = paste0(Base.dir, "/Level1_output/Human_placentas_8normal_pregnancy_GSE24129_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)

message("+-------------------------------------------------------------------------------------------------+")
message("+ L2_Human_6Control_GSE30186(GSM747304 - 747309) Array, quantile normalized                       +")
message("+ GPL10558 Illumina HumanHT-12 V4.0 expression beadchip 	                                         +")
message("+-------------------------------------------------------------------------------------------------+")

L2_matf <- matrix(scan(gzfile("GSE30186_series_matrix.txt.gz"), skip = 63, what = "character",
                       nlines = 47324, sep="\t", fill = T), ncol = 13, byrow = T)[,c(1, 8:13)]
colnames(L2_matf) <- L2_matf[1,]
L2_matf_new <- L2_matf[-1,]
L2_matf_new1 <- apply(L2_matf_new[,-1], 2, as.numeric)
L2_matf_new2 <- apply(L2_matf_new1, 1, mean)
GSE30186_Normal <- log2(L2_matf_new2+87.82)
L2_matf_new3 <- data.frame(L2_matf_new[,1], GSE30186_Normal)
colnames(L2_matf_new3) <- c("ID_REF", "GSE30186_Normal")
L2_matfg_GPL <- matrix(scan("GPL10558.annot.txt", what="character", fill = T, sep ="\t"), ncol = 3, byrow=T)
colnames(L2_matfg_GPL) <- c("ID_REF", "GeneName.Human", "EntrezID.Human")
Mergef1 <- merge(L2_matf_new3, L2_matfg_GPL, by = "ID_REF")

Hdat6 <- Mergef1
write.table(Hdat6, file = paste0(Base.dir, "/Level1_output/Human_placentas_6normal_pregnancy_GSE30186_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)

message("+-------------------------------------------------------------------------------------------------+")
message("+ L2_Human_8Control_GSE44711(GSM1089237 - 1089244) Array                                          +")
message("+ GPL10558 Illumina HumanHT-12 V4.0 expression beadchip 	                                         +")
message("+ Quantile normalized value of raw signal count. The raw signal value is the corrected, 
        background subtracted measurement of chemiluminescent signal, as outputted by the ABI 1700 software.+")
message("+--------------------------------------------------------------------------------------------------+")

L2_matg <- matrix(scan(gzfile("GSE44711_series_matrix.txt.gz"), skip = 69, what = "character",
                       nlines = 47228, sep="\t", fill = T), ncol = 17, byrow = T)[,c(1, 10:17)]
colnames(L2_matg) <- L2_matg[1,]
L2_matg <- L2_matg[-1,]
L2_matg1 <- apply(L2_matg[,-1], 2, as.numeric)
L2_matg2 <- apply(L2_matg1, 1, mean)
GSE44711_Normal <- log2(L2_matg2 + 12)
L2_matg3 <- data.frame(L2_matg[,1], GSE44711_Normal)
colnames(L2_matg3) <- c("ID_REF", "GSE44711_Normal")

Mergeg1 <- merge(L2_matg3, L2_matfg_GPL, by = "ID_REF")

Hdat8 <- Mergeg1

write.table(Hdat8, file = paste0(Base.dir, "/Level1_output/Human_placentas_8normal_pregnancy_GSE44711_april_2020.txt"), 
            row.names=F, col.names = T, quote = F)


message("+-----------------------------------------------------------------------------------+")
message("+ Save all of the human and mouse data                                              +")
message("+-----------------------------------------------------------------------------------+")

datalist = list(Hdat1, Hdat2, Hdat3, Hdat4, Hdat5, Hdat6, Hdat7, Hdat8, MDat1, MDat2, MDat3)

save(datalist, file = paste0(Base.dir, "/Level1_output/Mergedat_all_Normalised_Public_Control_Mouse_Human_april_2020.RData"))

message("+            write the human/mouse public data genes                                +")

HumanGeneNames <- unique(c(as.character(datalist[[1]]$GeneName.Human),as.character(datalist[[2]]$GeneName.Human),
                           as.character(datalist[[3]]$GeneName.Human),as.character(datalist[[4]]$GeneName.Human),
                           as.character(datalist[[5]]$GeneName.Human),as.character(datalist[[6]]$GeneName.Human),
                           as.character(datalist[[7]]$GeneName.Human),as.character(datalist[[8]]$GeneName.Human)))
MouseGeneNames <- unique(c(as.character(datalist[[9]]$GeneName.Mouse), as.character(datalist[[10]]$GeneName.Mouse), 
                           as.character(datalist[[11]]$GeneName.Mouse)))
print(length(HumanGeneNames)); print(length(MouseGeneNames))
write.csv(HumanGeneNames, file = paste0(out.dir,"/GEO_Control_Human_D8_N36552_GeneName_List.csv"), row.names=F)
write.csv(MouseGeneNames, file = paste0(out.dir,"/GEO_Control_Mouse_D3_N47936_GeneName_List.csv"), row.names=F)

message("+--- Old version of ensembl to generate the public data sets -------+")
## load("/Users/xz289/Documents/CTR_ans48_0003/Homolog_Protein_output/CTR_ans48_0003-GEOcontrol_Human_Mouse_MerData_All.RData")
## FinalH.GeneNames = unique(c(as.character(unique(GSE9984_Hmat.mer[,6])), as.character(unique(GSE28551_Hmat.mer[,23])),
##                            as.character(unique(GSE10588_Hmat.mer[,28])), as.character(unique(GSE4707_Hmat.mer[,6])),
##                            as.character(unique(GSE25906_Hmat.mer[,39])), as.character(unique(GSE30186_Hmat.mer[,8])),
##                            as.character(unique(GSE24129_Hmat.mer[,11])), as.character(unique(GSE44711_Hmat.mer[,10]))))

## FinalM.GeneNames <- unique(c(as.character(unique(GSE79121_Mmat[,2])), as.character(unique(GSE11224_MmatP.mer[,9])),
##                             as.character(unique(GSE11224_MmatD.mer[,7]))))

## write.csv(FinalH.GeneNames, file = paste0(out.dir,"/GEO_Control_Human_D8_N33804_GeneName_List.csv"), row.names=F)
## write.csv(FinalM.GeneNames, file = paste0(out.dir,"/GEO_Control_Mouse_D3_N47466_GeneName_List.csv"), row.names=F)


message("+-------------                        FIN                   ---------------------+")














