######################Figures 3, 4, 5., Supplementary Figures 3, 5, 6, 8.##############
#GSE132771 processing for cell signatures
#processed fron 10X output provided by authors
#see https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132771&format=file
.libPaths(c(.libPaths(), "/sc/wo/home/karmajx/R/x86_64-pc-linux-gnu-library/3.4", "/sc/wo/home/karmajx/R/x86_64-pc-linux-gnu-library/3.5"))
options(stringsAsFactors = F)
library(sctransform)
library(Seurat, lib.loc = "/sc/wo/home/karmajx/R/library")
library(httpuv, lib.loc = "/sc/wo/home/karmajx/R/library")
library(cowplot)
library(harmony)
library(magrittr)
library(foreach)
library(future)
library(dplyr)
gc()
for (mydir in list.dirs(recursive = F)){
  assign(paste0(gsub('./', '', mydir, fixed = T), '_10X'), Read10X(mydir,gene.column = 1))
}
gsmlist <- mget(ls(pattern='GSM'))
rm(list=ls(pattern = 'GSM[0-9]'))
for (i in 1:length(gsmlist)){
  assign(paste0(names(gsmlist)[i], '_seuratobject'), CreateSeuratObject(gsmlist[[i]], project = names(gsmlist)[i], min.cells = 5))
}

seuratlist <- mget(ls(pattern = 'seuratobject'))
for (i in 1:24){
  seuratlist[[i]] <- PercentageFeatureSet(seuratlist[[i]], pattern = "^MT-", col.name = "percent.mt")
  seuratlist[[i]] <- subset(seuratlist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
  seuratlist[[i]] <- NormalizeData(seuratlist[[i]], verbose = TRUE)
  seuratlist[[i]] <- FindVariableFeatures(seuratlist[[i]], selection.method = "vst", nfeatures = 2000)
  seuratlist[[i]] <- ScaleData(seuratlist[[i]], verbose = TRUE)
  seuratlist[[i]] <- RunPCA(seuratlist[[i]], npcs = 30, verbose = TRUE)
  seuratlist[[i]] <- RunTSNE(seuratlist[[i]], reduction = "pca", dims = 1:20)
  seuratlist[[i]] <- FindNeighbors(seuratlist[[i]], reduction = "pca", dims = 1:20)
  seuratlist[[i]] <- FindClusters(seuratlist[[i]], resolution = 0.6)
  seuratlist[[i]] <- RenameCells(seuratlist[[i]], add.cell.id=names(seuratlist)[i])
}

#GSE132771 Total lung cell suspension data
humanseurat_all <- seuratlist[c(seq(10, 24, by=2))]
ipf_all_anchors <- FindIntegrationAnchors(object.list = humanseurat_all, dims = 1:20)
ipf_all_combined <- IntegrateData(anchorset = ipf_all_anchors, dims = 1:20)
DefaultAssay(ipf_all_combined) <- "integrated"
ipf_all_combined <- ScaleData(ipf_all_combined, verbose = TRUE)
ipf_all_combined <- RunPCA(ipf_all_combined, npcs = 30, verbose = TRUE)
ipf_all_combined <- JackStraw(ipf_all_combined, dims = 30)
ipf_all_combined <- RunTSNE(ipf_all_combined, reduction = "pca", dims = 1:30)
ipf_all_combined <- FindNeighbors(ipf_all_combined, reduction = "pca", dims = 1:30)
ipf_all_combined <- FindClusters(ipf_all_combined, resolution = 0.8)
ipf_all_combined@meta.data$celltype <- plyr::mapvalues(ipf_all_combined@active.ident, from = sort(unique(ipf_all_combined@active.ident)), 
                                                       to = c('SPP1_monocytes_0', 'Infl_monocytes_1', 'ACKR1pos_endo_2', 'ACKR1neg_endo_3',
                                                              'Fibroblasts_4', 'AT2_5', 'Th_6', 'Pericytes_7', 'HLAhigh_mac_8',
                                                              'Sm_9', 'HLAhigh_mac_10', 'Bcells_11', 'Tc_12', 'AT1_13', 'PC_14', 'Endo_15',
                                                              'Ciliated_16', 'Monocytes_17', 'Monocytes_18', 'Cluster_19', 'Cluster_20',
                                                              'Bcells_21', 'Pericytes_22', 'AT2_23', 'Endo_24'))
#Supplementary Figure 3A:
Idents(ipf_all_combined) <- 'celltype'
DimPlot(ipf_all_combined, label = T, label.size = 6) + NoLegend()

#GSE132771 Lin- data (CD45-Epcam-CD235a-)
humanseurat_lin <- seuratlist[c(9, 11, 13, 15, 17, 19, 21, 23)]
ipf_lin_anchors <- FindIntegrationAnchors(object.list = humanseurat_lin, dims = 1:20)
ipf_lin_combined <- IntegrateData(anchorset = ipf_lin_anchors, dims = 1:20)
DefaultAssay(ipf_lin_combined) <- "integrated"
ipf_lin_combined <- ScaleData(ipf_lin_combined, verbose = TRUE)
ipf_lin_combined <- RunPCA(ipf_lin_combined, npcs = 30, verbose = TRUE)
ipf_lin_combined <- JackStraw(ipf_lin_combined, dims = 30)
ipf_lin_combined <- RunTSNE(ipf_lin_combined, reduction = "pca", dims = 1:30)
ipf_lin_combined <- FindNeighbors(ipf_lin_combined, reduction = "pca", dims = 1:30)
ipf_lin_combined <- FindClusters(ipf_lin_combined, resolution = 0.4)
ipf_lin_combined@meta.data$celltype <- plyr::mapvalues(ipf_lin_combined@active.ident, from = sort(unique(ipf_lin_combined@active.ident)), 
                                                       to = c('THY1high_alv_fib_0', 'THY1pos_sm_1', 'THY1neg_sm_2', 'CTHRC1pos_3', 'Adventitial_4',
                                                              'THY1neg_alv_fib_5', 'Pericytes_6', 'Peribronchial_7', 'Sm_8', 'Alveolar_9',
                                                              'Alveolar_10', 'Epi_11', 'Hematopoietic_12', 'Sm_13', 'Hematopoietic_14'))
Idents(ipf_lin_combined) <- 'celltype'
#Supplementary Figure 3B:
DimPlot(ipf_lin_combined, label = T, label.size = 6) + NoLegend()
FeaturePlot(ipf_lin_combined, 'THY1', pt.size = 1)
save(ipf_lin_combined, file = "/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_lin_noig_combined_symbols_fromraw_062620.RData")

#GSE135893 processing
setwd("/sc/wo/tri_data/jk/2020_05_Hs_IPF_GSE135893")
gse135893 <- readRDS('GSE135893_ILD_annotated_fullsize.rds') # RDS file deposited by authors of GSE135893
#see https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893_ILD_annotated_fullsize.rds.gz
#removing IG genes so they are not used as plasma cell signature
notigh <- row.names(gse135893@assays$RNA@counts)[setdiff(1:nrow(gse135893@assays$RNA@counts), grep('^IGH', row.names(gse135893@assays$RNA@counts)))]
notigl <- setdiff(notigh, notigh[grep('^IGL', notigh)])
notigk <- setdiff(notigl, notigl[grep('^IGK', notigl)])
notigk[grep('^IG', notigk)]

gse135893_ipf_control <- subset(gse135893, features = notigk, cells = row.names(subset(gse135893@meta.data, gse135893@meta.data$Diagnosis %in% c('IPF', 'Control'))))
DefaultAssay(gse135893_ipf_control) <- 'RNA'
gse135893_ipf_control <- NormalizeData(gse135893_ipf_control, verbose = TRUE)
gse135893_ipf_control <- FindVariableFeatures(gse135893_ipf_control, selection.method = "vst", nfeatures = 2000)
gse135893_ipf_control <- ScaleData(gse135893_ipf_control, verbose = TRUE)
gse135893_ipf_control <- RunPCA(gse135893_ipf_control, npcs = 30, verbose = TRUE)
#gse135893_ipf_control <- JackStraw(gse135893_ipf_control, num.replicate = 100)
#gse135893_ipf_control <- ScoreJackStraw(gse135893_ipf_control, dims = 1:20)
gse135893_ipf_control <- RunUMAP(gse135893_ipf_control, reduction = "pca", dims = 1:20)
gse135893_ipf_control <- FindNeighbors(gse135893_ipf_control, reduction = "pca", dims = 1:20)
gse135893_ipf_control <- FindClusters(gse135893_ipf_control, resolution = 1)
DimPlot(gse135893_ipf_control, label = T, label.size = 6, split.by = 'Diagnosis') + NoLegend()
DimPlot(gse135893_ipf_control, label = T, label.size = 4) + NoLegend()
FeaturePlot(gse135893_ipf_control, features = c('C1QA', 'SPP1', 'HLA-DBP1', 'ACKR1', 'GZMB'), pt.size = 1)
FeaturePlot(gse135893_ipf_control, features = c('COL1A1', 'CTHRC1', 'HAS1', 'HAS2'), pt.size = 1)

save(gse135893_ipf_control, file = 'gse135893_ipf_control_seuratnorm_062820.RData')
FeaturePlot(gse135893_ipf_control, features = c('SFTPC', 'SFTPA1', 'SFTPA2', 'MUC5B'), pt.size = 1)
FeaturePlot(gse135893_ipf_control, features = c('SOX4', 'ELF3', 'TSC22D1', 'NUPR1', 'JUND'), pt.size = 1, split.by = 'Diagnosis')
FeaturePlot(gse135893_ipf_control, features = c('GSN'), pt.size = 1, split.by = 'Diagnosis')
gse135893_ipf_control@meta.data$celltypenew <- plyr::mapvalues(gse135893_ipf_control@active.ident, 
                                                               from = sort(unique(gse135893_ipf_control@active.ident)), 
                                                               to = c('Ciliated_0', 'Ciliated_1', 'AT2_2',
                                                                      'SPP1_mac_3', 'C1QA_mac_4', 'C1QA_mac_5',
                                                                      'cDC_6', 'Mono_7', 'Tc_8',
                                                                      'C1QA_mac_9', 'Th_10', 'AT1_11', 'C1QA_mac_12',
                                                                      'AT2_13', 'ACKR1_pos_endo_14',
                                                                      'MUC5Bpos_AT1_15', 'ACKR1_neg_endo_16',
                                                                      'Basal_AT1_17', 'Diff_cil_18', 
                                                                      'Fibroblasts_19', 'ACKR1_neg_endo_20',
                                                                      'Monocytes_21', 'Prolif_mac_22', 'Fibroblasts_23', 
                                                                      'Ly_endo_24', 'Bcells_25', 'Sm_26', 'MC_27', 'PC_28',
                                                                      'AT2_29', 'AT2_30', 'Mesothelial_31', 'Mac_32'))
Idents(gse135893_ipf_control) <- 'celltypenew'

#Supplementary Figure 5A:
DimPlot(gse135893_ipf_control, label = T, label.size = 4) + NoLegend()

#calculate cell type signatures. Repeated for each dataset (GSE132771 total lung suspension, GSE132771 Lin- data, GSE135893).
#i.e. total 3 gene signature sets are developed.
#example of GSE132771 total lung shown.

data <- as.matrix(ipf_all_combined@assays$RNA@data)
celltypes <- as.character(ipf_all_combined@active.ident)
metadata <- ipf_all_combined@meta.data
cells <- unique(metadata[,"celltype"])
group <- unique(metadata[,"celltype"])
optimizeGeneList <- function(geneList,data,label){
  glist <- data.frame(geneIn=0,AUC=0,stringsAsFactors = F)
  geneList <- geneList[order(-geneList[,2]),]
  n <- ifelse(nrow(geneList)>100,100,nrow(geneList))
  for(j in c(1:n)){
    ge <- rownames(geneList)[1:j]
    geD2 <- apply(data[ge,,drop=FALSE],2,mean)
    gROC2 <- roc(label,geD2)
    gROC_auc2 <- as.numeric(auc(gROC2))
    
    glist[j,1] <- j
    glist[j,2] <- gROC_auc2
  }
  increaseThr <- 0.005
  maxAUC <- max(glist[,2])
  glist[,3] <- maxAUC-glist[,2]
  mu <- min(which(glist[,3]<increaseThr))
  mu <- ifelse(mu<10,10,mu)
  
  g <- rownames(geneList)[1:mu]
  return(g)
}
library(foreach)
library(doParallel)
library(pROC)
cl <- makeCluster(8)
registerDoParallel(cl)
betweenGroupCompare <- foreach(i = 1:length(unique(celltypes)), .packages=c("Seurat","pROC")) %dopar% {
  .libPaths(c(.libPaths(), "/sc/wo/home/karmajx/R/x86_64-pc-linux-gnu-library/3.4", "/sc/wo/home/karmajx/R/x86_64-pc-linux-gnu-library/3.5"))
  library(MAST)
  c <- unique(celltypes)[i]
  id <- Idents(ipf_all_combined)
  label <- ifelse(id==c,1,0)
  
  geneList <- FindMarkers(object=ipf_all_combined,ident.1=c,test.use="MAST",min.pct=0.05,only.pos = TRUE)
  g <- optimizeGeneList(geneList,data,label)
  
  r <- list(cellgroup=c,result=geneList,signature=g)
  return(r)
}
stopCluster(cl)

fcthr <- 0.38
fdrthr <- 0.05
cl <- makeCluster(length(betweenGroupCompare))
registerDoParallel(cl)
x <- foreach(i = 1:length(betweenGroupCompare), .packages=c("Seurat","pROC")) %dopar% {
  
  g <- betweenGroupCompare[[i]]$cellgroup
  sigR <- betweenGroupCompare[[i]]$result
  sig <- rownames(sigR)[sigR[,2]>fcthr & sigR[,5]<fdrthr]
  
  cellT <- as.character(unique(metadata[metadata[,"celltype"]==g,"celltype"]))
  
  cellSig <- list()
  ci <- 1
  if(length(cellT)>1){
    subm <- metadata[metadata[,"celltype"]==g,]
    subd <- data[sig,rownames(subm)]
    subo <- CreateSeuratObject(subd)
    subo <- AddMetaData(object=subo,metadata=subm)
    subo <- SetIdent(object=subo,value=subm[,"celltype"])
    for(j in c(1:length(cellT))){
      id <- Idents(subo)
      label <- ifelse(id==cellT[j],1,0)
      genelist <- FindMarkers(object=subo,ident.1=cellT[j],test.use="MAST",min.pct=0.05,only.pos = TRUE)
      g <- optimizeGeneList(genelist,subd,label)
      
      r <- list(celltype=cellT[j],genelist=genelist,signature=g)
      cellSig[[ci]] <- r
      ci <- ci+1
    }
  }else{
    r <- list(celltype=g,genelist=sigR,signature=betweenGroupCompare[[i]]$signature)
    cellSig[[ci]] <- r
    ci <- ci+1
  }
  return(cellSig)
}
stopCluster(cl)

gene_sig <- data.frame(celltype=NA,signature=NA,stringsAsFactors = F)
gi <- 1
for(i in c(1:length(x))){
  a <- x[[i]]
  if(length(a)==1){
    l <- a[[1]]$celltype
    b <- a[[1]]$signature
    gene_sig[(gi:(gi+length(b)-1)),1] <- rep(l,length(b))
    gene_sig[(gi:(gi+length(b)-1)),2] <- b
    gi <- gi+length(b)
  }else{
    for(j in c(1:length(a))){
      l <- a[[j]]$celltype
      b <- a[[j]]$signature
      gene_sig[(gi:(gi+length(b)-1)),1] <- rep(l,length(b))
      gene_sig[(gi:(gi+length(b)-1)),2] <- b
      gi <- gi+length(b)
    }
  }
}
write.table(gene_sig,file="ipf_all_combined_signature.txt",row.names=F,col.names = T,sep="\t",quote=F)

#########plot heatmap######
#gene_sig <- gene_sig[!is.na(gene_sig[,2]),]
gene_sig <- cellsignatures
c <- unique(gene_sig[,1])
g <- unique(gene_sig[,2])
cgMatrix <- array(0,dim=c(length(c),length(g)))
rownames(cgMatrix) <- c
colnames(cgMatrix) <- g

for(i in c(1:nrow(gene_sig))){
  cgMatrix[gene_sig[i,1],gene_sig[i,2]] <- 1
}

library(pheatmap)

data <- gse135893_ipf_control@assays$RNA@data
metadata <- gse135893_ipf_control@meta.data
data <- as.matrix(data)

filter_d <- data[g,]
filter_d <- 2^filter_d-1
c <- unique(gene_sig[,1])
finalSig <- array(0,dim=c(length(g),length(c)))
rownames(finalSig) <- g
colnames(finalSig) <- c
for(i in c(1:length(c))){
  k <- c[i]
  l <- rownames(metadata)[metadata[,"celltypenew"]==k]
  b <- filter_d[,l]
  b <- apply(b,1,mean)
  if(length(which(is.na(b)))>0){
    cat("I:",i,"\n")
  }
  finalSig[,i] <- b
}

a <- (finalSig-apply(finalSig,1,mean))/apply(finalSig,1,sd)
a[a<(-2)] <- (-2)
a[a>2] <- 2

library(pheatmap)
png("IPF_gse135893_lin_manuscript_heatmap.png",width=1000,height=1000,res=100)

red <- colorRampPalette(c("blue","white","red"))(255)
mybreak <- seq(-2,2,length.out=256)

pheatmap(a,color=red,breaks=mybreak,scale="none",cluster_rows=F,cluster_cols=F,show_rownames=FALSE,show_colnames=T)

dev.off()

write.table(finalSig,file="IPF_scRNAseq_17cells_signature_MAST_referenceMatrix.txt",row.names=T,col.names=T,sep="\t",quote=F)


#calculate signature scores in GSE47460. Repeated for each signature set (GSE132771 total lung suspension, GSE132771 Lin- data, GSE135893).
#example of GSE132771 total lung shown.

#correct Sheppard lin
cellsignatures_lin <- read.table('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_lin_scRNAseq_cells_signature_MAST_062620.txt', header = T)
#correct Sheppard all
cellsignatures <- read.table('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_all_res08_scRNAseq_cells_signature_MAST_062620.txt', header = T)
#correct GSE135893
cellsignatures_gse135893 <- read.table('/sc/wo/tri_data/jk/2020_05_Hs_IPF_GSE135893/gse135893_ipf_control_res1_noig_seuratnorm_signature_062720.txt', header = T)

cellsignatures$celltype <- gsub('-', '_', cellsignatures$celltype, fixed = T)
celltypes <- unique(cellsignatures$celltype)
cellsignatures_genes <- lapply(celltypes, function(k) subset(cellsignatures, celltype==k)$signature)
names(cellsignatures_genes) <- celltypes
sort(unlist(lapply(cellsignatures_genes, function(k) length(intersect(row.names(gse47460_matrix), k)))))

cellsignatures_lin$celltype <- gsub('-', '_', cellsignatures_lin$celltype, fixed = T)
celltypes_lin <- unique(cellsignatures_lin$celltype)
cellsignatures_lin_genes <- lapply(celltypes_lin, function(k) subset(cellsignatures_lin, celltype==k)$signature)
names(cellsignatures_lin_genes) <- celltypes_lin
sort(unlist(lapply(cellsignatures_lin_genes, function(k) length(intersect(row.names(gse47460_matrix), k)))))

cellsignatures_gse135893$celltype <- gsub('-', '_', cellsignatures_gse135893$celltype, fixed = T)
celltypes_gse135893 <- unique(cellsignatures_gse135893$celltype)
cellsignatures_gse135893_genes <- lapply(celltypes_gse135893, function(k) subset(cellsignatures_gse135893, celltype==k)$signature)
names(cellsignatures_gse135893_genes) <- celltypes_gse135893
sort(unlist(lapply(cellsignatures_gse135893_genes, function(k) length(intersect(row.names(gse47460_matrix), k)))))

#removed 'Cluster_19 in GSE135893' (mitochondrial genes high) as it has no overlap with genes in GSE47460
cellsignatures_genes <- cellsignatures_genes[c(1:15, 17:25)]
names(cellsignatures_genes)

gse47460_celltypegsva <- lapply(cellsignatures_genes, function(k) gsva(as.matrix(gse47460_matrix), gset.idx.list = list(k)))
names(gse47460_celltypegsva) <- names(cellsignatures_genes)
gse47460_celltypegsva_table <- as.data.frame(do.call(rbind, gse47460_celltypegsva))
rownames(gse47460_celltypegsva_table) <- names(gse47460_celltypegsva)

gse47460_celltypegsva_lin <- lapply(cellsignatures_lin_genes, function(k) gsva(as.matrix(gse47460_matrix), gset.idx.list = list(k)))
names(gse47460_celltypegsva_lin) <- names(cellsignatures_lin_genes)
gse47460_celltypegsva_lin_table <- as.data.frame(do.call(rbind, gse47460_celltypegsva_lin))
rownames(gse47460_celltypegsva_lin_table) <- names(gse47460_celltypegsva_lin)

gse47460_celltypegsva_gse135893 <- lapply(cellsignatures_gse135893_genes, function(k) gsva(as.matrix(gse47460_matrix), gset.idx.list = list(k)))
names(gse47460_celltypegsva_gse135893) <- names(cellsignatures_gse135893_genes)
gse47460_celltypegsva_gse135893_table <- as.data.frame(do.call(rbind, gse47460_celltypegsva_gse135893))
rownames(gse47460_celltypegsva_gse135893_table) <- names(gse47460_celltypegsva_gse135893)

#Supplementary Figures 3C, 3D, 5B.
#correlation matrix from GSVA scores
#Supplementary Figure 3C.
corr <- round(cor(t(gse47460_celltypegsva_table)), 3)
p.mat <- cor_pmat(t(gse47460_celltypegsva_table))
ggcorrplot(corr, hc.order = TRUE, outline.col = "white")
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           outline.col = "white",
           #ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), lab = T)

#Supplementary Figure 3D.
corr_lin <- round(cor(t(gse47460_celltypegsva_lin_table)), 3)
p.mat_lin <- cor_pmat(t(gse47460_celltypegsva_lin_table))
ggcorrplot(corr_lin, hc.order = TRUE, outline.col = "white")
ggcorrplot(corr_lin, hc.order = TRUE, type = "lower",
           outline.col = "white",
           #ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), lab = T)

#Supplementary Figure 5B.
corr_gse135893 <- round(cor(t(gse47460_celltypegsva_gse135893_table)), 3)
p.mat_gse135893 <- cor_pmat(t(gse47460_celltypegsva_gse135893_table))
ggcorrplot(corr_gse135893, hc.order = TRUE, outline.col = "white")
ggcorrplot(corr_gse135893, hc.order = TRUE, type = "lower",
           outline.col = "white",
           #ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), lab = T)

gse47460_celltypegsva_consensus <- merge(consensusdf2_gse47460, t(gse47460_celltypegsva_table), by.x = 0, by.y = 0)
rownames(gse47460_celltypegsva_consensus) <- gse47460_celltypegsva_consensus$Row.names
gse47460_celltypegsva_consensus$Row.names <- NULL

gse47460_celltypegsva_test <- apply(gse47460_celltypegsva_consensus[,2:ncol(gse47460_celltypegsva_consensus)], 2, 
                                    function(k) dunnTest(k~consensusclass, data = gse47460_celltypegsva_consensus))

gse47460_celltypegsva_lin_consensus <- merge(consensusdf2_gse47460, t(gse47460_celltypegsva_lin_table), by.x = 0, by.y = 0)
rownames(gse47460_celltypegsva_lin_consensus) <- gse47460_celltypegsva_lin_consensus$Row.names
gse47460_celltypegsva_lin_consensus$Row.names <- NULL

gse47460_celltypegsva_lin_test <- apply(gse47460_celltypegsva_lin_consensus[,2:ncol(gse47460_celltypegsva_lin_consensus)], 2, 
                                        function(k) dunnTest(k~consensusclass, data = gse47460_celltypegsva_lin_consensus))

gse47460_celltypegsva_gse135893_consensus <- merge(consensusdf2_gse47460, t(gse47460_celltypegsva_gse135893_table), by.x = 0, by.y = 0)
rownames(gse47460_celltypegsva_gse135893_consensus) <- gse47460_celltypegsva_gse135893_consensus$Row.names
gse47460_celltypegsva_gse135893_consensus$Row.names <- NULL

gse47460_celltypegsva_gse135893_test <- apply(gse47460_celltypegsva_gse135893_consensus[,2:ncol(gse47460_celltypegsva_gse135893_consensus)], 2, 
                                              function(k) dunnTest(k~consensusclass, data = gse47460_celltypegsva_gse135893_consensus))

#Figures 3 and 4.
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
rm(target_gg_gse47460)
for (i in colnames(gse47460_celltypegsva_consensus)[2:ncol(gse47460_celltypegsva_consensus)]){
  assign(paste0('gg_celltypes2_gse47460_', i), ggplot(gse47460_celltypegsva_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse47460 <- mget(ls(pattern = 'gg_celltypes2_gse47460_'))
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
do.call(plot_grid, target_gg_gse47460)

#Figure 5.
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
rm(target_gg_gse47460)
for (i in colnames(gse47460_celltypegsva_lin_consensus)[2:ncol(gse47460_celltypegsva_lin_consensus)]){
  assign(paste0('gg_celltypes2_gse47460_', i), ggplot(gse47460_celltypegsva_lin_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse47460 <- mget(ls(pattern = 'gg_celltypes2_gse47460_'))
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
do.call(plot_grid, target_gg_gse47460)

#Supplementary Figure 6.
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
rm(target_gg_gse47460)
for (i in colnames(gse47460_celltypegsva_gse135893_consensus)[2:ncol(gse47460_celltypegsva_gse135893_consensus)]){
  assign(paste0('gg_celltypes2_gse47460_', i), ggplot(gse47460_celltypegsva_gse135893_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse47460 <- mget(ls(pattern = 'gg_celltypes2_gse47460_'))
rm(list=ls(pattern = 'gg_celltypes2_gse47460_'))
do.call(plot_grid, target_gg_gse47460)