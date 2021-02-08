#GSE135893 signatures
cellsignatures_gse135893 <- read.table('/sc/wo/tri_data/jk/2020_05_Hs_IPF_GSE135893/gse135893_ipf_control_res1_noig_seuratnorm_signature_062720.txt', header = T)

cellsignatures_gse135893$celltype <- gsub('-', '_', cellsignatures_gse135893$celltype, fixed = T)
celltypes_gse135893 <- unique(cellsignatures_gse135893$celltype)
cellsignatures_gse135893_genes <- lapply(celltypes_gse135893, function(k) subset(cellsignatures_gse135893, celltype==k)$signature)
names(cellsignatures_gse135893_genes) <- celltypes_gse135893
sort(unlist(lapply(cellsignatures_gse135893_genes, function(k) length(intersect(row.names(gse134692), k)))))

cellsignatures_gse135893_genes_ensembl <- lapply(cellsignatures_gse135893_genes, function(k) row.names(subset(gse134692_annotation, GeneName %in% k)))
sort(unlist(lapply(cellsignatures_gse135893_genes_ensembl, function(k) length(intersect(k, row.names(gse134692))))))

gse134692_celltypegsva_gse135893 <- lapply(cellsignatures_gse135893_genes_ensembl, function(k) gsva(as.matrix(gse134692), gset.idx.list = list(k)))
names(gse134692_celltypegsva_gse135893) <- names(cellsignatures_gse135893_genes_ensembl)
gse134692_celltypegsva_gse135893_table <- do.call(rbind, gse134692_celltypegsva_gse135893)
rownames(gse134692_celltypegsva_gse135893_table) <- names(gse134692_celltypegsva_gse135893)
gse134692_celltypegsva_gse135893_consensus <- merge(consensusdf2_gse134692, t(gse134692_celltypegsva_gse135893_table), by.x = 0, by.y = 0)
rownames(gse134692_celltypegsva_gse135893_consensus) <- gse134692_celltypegsva_gse135893_consensus$Row.names
gse134692_celltypegsva_gse135893_consensus$Row.names <- NULL

#Supplementary Figure 7.
rm(list=ls(pattern = 'gg_celltypes2_gse134692_gse135893_signatures_'))
rm(target_gg_gse134692_gse135893_signatures)
for (i in colnames(gse134692_celltypegsva_gse135893_consensus)[3:ncol(gse134692_celltypegsva_gse135893_consensus)]){
  assign(paste0('gg_celltypes2_gse134692_gse135893_signatures_', i), ggplot(gse134692_celltypegsva_gse135893_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse134692_gse135893_signatures <- mget(ls(pattern = 'gg_celltypes2_gse134692_gse135893_signatures_'))
rm(list=ls(pattern = 'gg_celltypes2_gse134692_gse135893_signatures_'))
do.call(plot_grid, target_gg_gse134692_gse135893_signatures)

gse134692_gse135893_signatures_test <- apply(gse134692_celltypegsva_gse135893_consensus[,3:ncol(gse134692_celltypegsva_gse135893_consensus)], 2, 
                                             function(k) dunnTest(k~reverseclass, data = gse134692_celltypegsva_gse135893_consensus))

####################sheppard signatures all cells###################
#correct Sheppard lin
cellsignatures_lin <- read.table('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_lin_scRNAseq_cells_signature_MAST_062620.txt', header = T)
#correct Sheppard all
cellsignatures <- read.table('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_all_res08_scRNAseq_cells_signature_MAST_062620.txt', header = T)

cellsignatures$celltype <- gsub('-', '_', cellsignatures$celltype, fixed = T)
celltypes <- unique(cellsignatures$celltype)
cellsignatures_genes <- lapply(celltypes, function(k) subset(cellsignatures, celltype==k)$signature)
names(cellsignatures_genes) <- celltypes
sort(unlist(lapply(cellsignatures_genes, function(k) length(intersect(row.names(gse134692), k)))))

cellsignatures_genes_ensembl <- lapply(cellsignatures_genes, function(k) row.names(subset(gse134692_annotation, GeneName %in% k)))
sort(unlist(lapply(cellsignatures_genes_ensembl, function(k) length(intersect(k, row.names(gse134692))))))

gse134692_celltypegsva <- lapply(cellsignatures_genes_ensembl, function(k) gsva(as.matrix(gse134692), gset.idx.list = list(k)))
names(gse134692_celltypegsva) <- names(cellsignatures_genes_ensembl)
gse134692_celltypegsva_table <- do.call(rbind, gse134692_celltypegsva)
rownames(gse134692_celltypegsva_table) <- names(gse134692_celltypegsva)
gse134692_celltypegsva_consensus <- merge(consensusdf2_gse134692, t(gse134692_celltypegsva_table), by.x = 0, by.y = 0)
rownames(gse134692_celltypegsva_consensus) <- gse134692_celltypegsva_consensus$Row.names
gse134692_celltypegsva_consensus$Row.names <- NULL

#Supplementary Figure 7.
rm(list=ls(pattern = 'gg_celltypes2_gse134692_signatures_'))
rm(target_gg_gse134692_signatures)
for (i in colnames(gse134692_celltypegsva_consensus)[3:ncol(gse134692_celltypegsva_consensus)]){
  assign(paste0('gg_celltypes2_gse134692_signatures_', i), ggplot(gse134692_celltypegsva_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse134692_signatures <- mget(ls(pattern = 'gg_celltypes2_gse134692_signatures_'))
rm(list=ls(pattern = 'gg_celltypes2_gse134692_signatures_'))
do.call(plot_grid, target_gg_gse134692_signatures)

gse134692_signatures_test <- apply(gse134692_celltypegsva_consensus[,3:ncol(gse134692_celltypegsva_consensus)], 2, 
                                             function(k) dunnTest(k~consensusclass, data = gse134692_celltypegsva_consensus))

####################sheppard signatures lin###################
#correct Sheppard lin
cellsignatures_lin <- read.table('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/ipf_control_lin_scRNAseq_cells_signature_MAST_062620.txt', header = T)

cellsignatures_lin$celltype <- gsub('-', '_', cellsignatures_lin$celltype, fixed = T)
celltypes <- unique(cellsignatures_lin$celltype)
cellsignatures_lin_genes <- lapply(celltypes, function(k) subset(cellsignatures_lin, celltype==k)$signature)
names(cellsignatures_lin_genes) <- celltypes
sort(unlist(lapply(cellsignatures_lin_genes, function(k) length(intersect(row.names(gse134692), k)))))

cellsignatures_lin_genes_ensembl <- lapply(cellsignatures_lin_genes, function(k) row.names(subset(gse134692_annotation, GeneName %in% k)))
sort(unlist(lapply(cellsignatures_lin_genes_ensembl, function(k) length(intersect(k, row.names(gse134692))))))

gse134692_lin_celltypegsva <- lapply(cellsignatures_lin_genes_ensembl, function(k) gsva(as.matrix(gse134692), gset.idx.list = list(k)))
names(gse134692_lin_celltypegsva) <- names(cellsignatures_lin_genes_ensembl)
gse134692_lin_celltypegsva_table <- do.call(rbind, gse134692_lin_celltypegsva)
rownames(gse134692_lin_celltypegsva_table) <- names(gse134692_lin_celltypegsva)
gse134692_lin_celltypegsva_consensus <- merge(consensusdf2_gse134692, t(gse134692_lin_celltypegsva_table), by.x = 0, by.y = 0)
rownames(gse134692_lin_celltypegsva_consensus) <- gse134692_lin_celltypegsva_consensus$Row.names
gse134692_lin_celltypegsva_consensus$Row.names <- NULL

#Supplementary Figure 7.
rm(list=ls(pattern = 'gg_celltypes2_gse134692_lin_signatures_'))
rm(target_gg_gse134692_lin_signatures)
for (i in colnames(gse134692_lin_celltypegsva_consensus)[3:ncol(gse134692_lin_celltypegsva_consensus)]){
  assign(paste0('gg_celltypes2_gse134692_lin_signatures_', i), ggplot(gse134692_lin_celltypegsva_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none") + 
           ylab('Signature score') + ggtitle(i) + theme(plot.title = element_text(size = 12), legend.position = "none")) + xlab('consensusclass')
}
target_gg_gse134692_lin_signatures <- mget(ls(pattern = 'gg_celltypes2_gse134692_lin_signatures_'))
rm(list=ls(pattern = 'gg_celltypes2_gse134692_lin_signatures_'))
do.call(plot_grid, target_gg_gse134692_lin_signatures)

gse134692_lin_signatures_test <- apply(gse134692_lin_celltypegsva_consensus[,3:ncol(gse134692_lin_celltypegsva_consensus)], 2, 
                                   function(k) dunnTest(k~consensusclass, data = gse134692_lin_celltypegsva_consensus))
