######################Figure 6.#########################
#Figure 6A generated using ArrayStudio. We cannot post code for this figure here.
#Figure 6B (code below is repeated for each chemokine):
FeaturePlot(gse135893_ipf_control, features = 'CCL5', pt.size = 1)
ggplot(gse47460_phenom_consensus, aes(consensusclass, CCL5, fill = consensusclass)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none")
dunnTest(CXCL17~consensusclass, data = gse47460_phenom_consensus)

######################Figure 7., Supplementary Figures 8 and 9.####################
#process GSE135893 to separate 'Ciliated_high' and 'Ciliated_low donors
gse135893 <- readRDS('GSE135893_ILD_annotated_fullsize.rds')
gse135893_ipf <- subset(gse135893, cells = row.names(subset(gse135893@meta.data, gse135893@meta.data$Diagnosis=='IPF')))
DefaultAssay(gse135893_ipf) <- 'RNA'
gse135893_ipf <- PercentageFeatureSet(gse135893_ipf, pattern = "^MT-", col.name = "percent.mt")
gse135893_ipf <- NormalizeData(gse135893_ipf, verbose = TRUE)
gse135893_ipf <- FindVariableFeatures(gse135893_ipf, selection.method = "vst", nfeatures = 2000)
gse135893_ipf <- ScaleData(gse135893_ipf, verbose = TRUE, vars.to.regress = 'percent.mt')

gse135893_ipf <- RunPCA(gse135893_ipf, npcs = 30, verbose = TRUE)
gse135893_ipf <- JackStraw(gse135893_ipf, num.replicate = 100)
gse135893_ipf <- ScoreJackStraw(gse135893_ipf, dims = 1:20)
gse135893_ipf <- RunUMAP(gse135893_ipf, reduction = "pca", dims = 1:20)
gse135893_ipf <- FindNeighbors(gse135893_ipf, reduction = "pca", dims = 1:20)
gse135893_ipf <- FindClusters(gse135893_ipf, resolution = 0.6)

gse135893_fibroblasts_ciliatedhighvslow <- FindMarkers(gse135893_ipf, test.use = 'MAST', logfc.threshold = 0,
                                                       subset.ident = 'Fibroblasts_17', group.by = 'ciliated', 
                                                       ident.1 = 'Ciliated_high', ident.2 = 'Ciliated_low')
gse135893_fibroblasts_ciliatedhighvslow$gene <- row.names(gse135893_fibroblasts_ciliatedhighvslow)
gse135893_ipf@meta.data$celltypenew <- plyr::mapvalues(gse135893_ipf@active.ident, 
                                                       from = sort(unique(gse135893_ipf@active.ident)),
                                                       to = c('SPP1pos_macs_0', 'Ciliated_1', 'C1QA_mac_2', 'Ciliated_3', 'AT1_4', 'AT2_5', 'C1QA_mac_6', 
                                                              'ACKR1pos_endo_7', 'Monocytes_8', 'AT1_9', 'Th_10', 'AT1_11', 'Macs_12', 'Tc_13', 'HAS1_fibro_14',
                                                              'Diff_ciliated_15', 'ACKR1neg_endo_16', 'Fibroblasts_17', 'Prolif_macs_18', 'Lymph_endo_19',
                                                              'Sm_20', 'Bcells_21', 'Macs_22', 'PC_23', 'AT2_24', 'MC_25', 'AT1_26', 'Macs_27', 'Ciliated_28', 
                                                              'Fibroblast_29'))
Idents(gse135893_ipf) <- 'celltypenew'
DimPlot(gse135893_ipf, label = T, label.size = 4) + NoLegend()
save(gse135893_ipf, file = 'gse135893_ipfonly_seurat_withseuratnorm_062620.RData')

cellclusters <- function(x){
  abstable <- as.matrix.data.frame(table(x@meta.data$orig.ident, x@meta.data$celltypenew))
  perctable <- round(prop.table(abstable, 1)*100, 2)
  colnames(perctable) <- colnames(table(x@meta.data$orig.ident, x@meta.data$celltypenew))
  rownames(perctable) <- row.names(table(x@meta.data$orig.ident, x@meta.data$celltypenew))
  perctable <- as.data.frame(perctable)
  return(perctable)
}
gse135893_ipf_perctable <- cellclusters(gse135893_ipf)
#celltypebydonor <- as.matrix.data.frame(table(gse135893_ipf@meta.data$orig.ident, gse135893_ipf@meta.data$celltype))
write.csv(gse135893_ipf_perctable, file = 'gse135893_ipfonly_celltypebydonorpercentage_seuratnorm.csv')
gse135893_ipf@meta.data$ciliated <- plyr::mapvalues(gse135893_ipf@meta.data$orig.ident, 
                                                    from = sort(unique(gse135893_ipf@meta.data$orig.ident)),
                                                    to = c(rep('Ciliated_low', 8), rep('Ciliated_high', 4), rep('Ciliated_low', 2), 
                                                           rep('Ciliated_high', 2), 'Ciliated_low', 'Ciliated_high', 'Ciliated_low'))

#Supplementary Figure 8:
DimPlot(gse135893_ipf, label = T, label.size = 4, split.by = 'ciliated') + NoLegend()

#Supplementary Table 3:
gse135893_ipf_perctable$totalciliated <- gse135893_ipf_perctable$Ciliated_1+gse135893_ipf_perctable$Ciliated_2+gse135893_ipf_perctable$Ciliated_3+
  gse135893_ipf_perctable$Ciliated_28+gse135893_ipf_perctable$Diff_ciliated_15
ciliatedvector <- vector()
for (i in 1:nrow(gse135893_ipf_perctable)){
  if (gse135893_ipf_perctable$totalciliated[i]>20){
    ciliatedvector[i] <- 'ciliated_high'
  }else{
    ciliatedvector[i] <- 'ciliated_low'
  }
}

gse135893_ipf_perctable$totalmyeloid <- gse135893_ipf_perctable$SPP1pos_macs_0 + gse135893_ipf_perctable$C1QA_mac_2 + gse135893_ipf_perctable$C1QA_mac_6 +  
  gse135893_ipf_perctable$Monocytes_8 + gse135893_ipf_perctable$Macs_12 + gse135893_ipf_perctable$Prolif_macs_18 + gse135893_ipf_perctable$Macs_22 + gse135893_ipf_perctable$Macs_27

gse135893_ipf_perctable$totalendothelial <- gse135893_ipf_perctable$ACKR1pos_endo_7 + gse135893_ipf_perctable$ACKR1neg_endo_16 +gse135893_ipf_perctable$Lymph_endo_19

gse135893_ipf_perctable$ciliatedphenotype <- ciliatedvector
ciliatedtests <- apply(gse135893_ipf_perctable[,1:31], 2, function(k) t.test(k~ciliatedphenotype, data = gse135893_ipf_perctable))
ciliated_p <- lapply(ciliatedtests, function(k) c(k$p.value, (k$estimate[1]-k$estimate[2])))
ciliatedteststable <- as.data.frame(do.call(rbind, ciliated_p))
colnames(ciliatedteststable) <- c('p_value', 'difference')
ciliatedteststable$padjust <- p.adjust(ciliatedteststable$p_value, method = 'BH')

#Supplementary Figure 8B. Repeated for all three populations.
ggplot(gse135893_ipf_perctable, aes(ciliatedphenotype, totalciliated, fill = ciliatedphenotype)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#999999", "#E69F00")) + theme(legend.position = "none")

#############################ciliated clusterprofiler#####################
#Supplementary Tables 2A and 2B are generated in IPA. We cannot post code for these tables here.

gse135893_ipf_ciliatedhighvshealthy <- FindMarkers(gse135893_ipf_control_noig, test.use = 'MAST', 
                                                   group.by = 'ciliated', ident.1 = 'Ciliated_high', ident.2 = 'Healthy')
gse135893_ipf_ciliatedhighvshealthy$gene <- row.names(gse135893_ipf_ciliatedhighvshealthy)
gse135893_ipf_ciliatedhighvshealthy_sig <- subset(gse135893_ipf_ciliatedhighvshealthy, p_val_adj<0.05)
xlsx::write.xlsx(gse135893_ipf_ciliatedhighvshealthy_sig[,c(6, 2, 5)], file = 'gse135893_ipf_ciliatedhighvshealthy.xlsx', row.names = F)

gse135893_ipf_ciliatedlowvshealthy <- FindMarkers(gse135893_ipf_control_noig, test.use = 'MAST', 
                                                  group.by = 'ciliated', ident.1 = 'Ciliated_low', ident.2 = 'Healthy')
gse135893_ipf_ciliatedlowvshealthy$gene <- row.names(gse135893_ipf_ciliatedlowvshealthy)
gse135893_ipf_ciliatedlowvshealthy_sig <- subset(gse135893_ipf_ciliatedlowvshealthy, p_val_adj<0.05)
xlsx::write.xlsx(gse135893_ipf_ciliatedlowvshealthy_sig[,c(6, 2, 5)], file = 'gse135893_ipf_ciliatedlowvshealthy.xlsx', row.names = F)

ciliatedlowlist <- row.names(gse135893_ipf_ciliatedlowvshealthy_sig)
ciliatedlowentrez <- bitr(ciliatedlowlist, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
ciliatedlowlistentrez <- merge(ciliatedlowentrez, gse135893_ipf_ciliatedlowvshealthy_sig, by.x = 'SYMBOL', by.y='gene')
ciliatedlowlistentrez <- ciliatedlowlistentrez[order(ciliatedlowlistentrez$avg_logFC, decreasing = T),]
ciliatedlowlist_vector <- ciliatedlowlistentrez[,4]
names(ciliatedlowlist_vector) <- ciliatedlowlistentrez$ENTREZID
ciliatedlowlist_go <- gseGO(geneList = ciliatedlowlist_vector, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, 
                            minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.2, verbose = TRUE)
ciliatedlowlist_kegg <- gseKEGG(geneList = ciliatedlowlist_vector, organism = 'hsa', nPerm = 1000, minGSSize = 15, 
                                pvalueCutoff = 0.2, verbose = TRUE)
ciliatedlowlist_goresults <- ciliatedlowlist_go@result
ciliatedlowlist_keggresults <- ciliatedlowlist_kegg@result

ciliatedlowlist_reactome <- gsePathway(ciliatedlowlist_vector, nPerm=10000, minGSSize = 15,
                                       pvalueCutoff=1, pAdjustMethod="BH", verbose=TRUE)
ciliatedlowlist_reactome_results <- ciliatedlowlist_reactome@result
ciliatedlowlist_reactome_results <- ciliatedlowlist_reactome_results[order(ciliatedlowlist_reactome_results$NES, decreasing = T),]
#Supplementary Table 2C
ciliatedlowlist_reactome_results_up <- ciliatedlowlist_reactome_results[which(ciliatedlowlist_reactome_results$NES>0),]
xlsx::write.xlsx(ciliatedlowlist_reactome_results_up, file = 'gse135893_ciliatedlowlist_reactome_results_up.xlsx')
#emapplot(ciliatedlowlist_reactome)

ciliatedhighlist <- row.names(gse135893_ipf_ciliatedhighvshealthy_sig)
ciliatedhighentrez <- bitr(ciliatedhighlist, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
ciliatedhighlistentrez <- merge(ciliatedhighentrez, gse135893_ipf_ciliatedhighvshealthy_sig, by.x = 'SYMBOL', by.y='gene')
ciliatedhighlistentrez <- ciliatedhighlistentrez[order(ciliatedhighlistentrez$avg_logFC, decreasing = T),]
ciliatedhighlist_vector <- ciliatedhighlistentrez[,4]
names(ciliatedhighlist_vector) <- ciliatedhighlistentrez$ENTREZID
ciliatedhighlist_go <- gseGO(geneList = ciliatedhighlist_vector, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, 
                             minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.2, verbose = TRUE)
ciliatedhighlist_kegg <- gseKEGG(geneList = ciliatedhighlist_vector, organism = 'hsa', nPerm = 1000, minGSSize = 15, 
                                 pvalueCutoff = 0.2, verbose = TRUE)
ciliatedhighlist_goresults <- ciliatedhighlist_go@result
ciliatedhighlist_keggresults <- ciliatedhighlist_kegg@result
ciliatedhighlist_reactome <- gsePathway(ciliatedhighlist_vector, nPerm=10000, minGSSize = 15,
                                        pvalueCutoff=1, pAdjustMethod="BH", verbose=TRUE)
ciliatedhighlist_reactome_results <- ciliatedhighlist_reactome@result
ciliatedhighlist_reactome_results <- ciliatedhighlist_reactome_results[order(ciliatedhighlist_reactome_results$NES, decreasing = T),]
#Supplementary Table 2D
ciliatedhighlist_reactome_results_up <- ciliatedhighlist_reactome_results[which(ciliatedhighlist_reactome_results$NES>0),]
xlsx::write.xlsx(ciliatedhighlist_reactome_results_up, file = 'gse135893_ciliatedhighlist_reactome_results_up.xlsx')

#Figure 7A:
hist(gse135893_ipf_perctable$totalciliated, col = 'blue')

#PyMINEr and NicheNet
library(CBDD) # Implementation of PyMINEr by Clarivate Analytics
library(CBDDnetworks)
library(nichenetr)
load('~/Rpackages/nichenetr_start_060420.RData') # compilation of NicheNet ligand receptor files by authors of NicheNet
#https://github.com/saeyslab/nichenetr

gse135893_ipfmatrix <- as.matrix(gse135893_ipf@assays$RNA@data)
gse135893_ipfmatrix <- gse135893_ipfmatrix[which(rowSums(gse135893_ipfmatrix)>0),]
gse135893_ipfidents <- as.data.frame(gse135893_ipf@active.ident)
colnames(gse135893_ipfidents) <- 'celltype'
gse135893_ipfidents$celltype <- as.character(gse135893_ipfidents$celltype)
gse135893_ipfidents$cellnames <- row.names(gse135893_ipfidents)
#ligand-receptor interactions obtained from PMID 26198319, formatted to file named 'receptor_ligand.csv'
reclig <- read.csv('/sc/wo/tri_data/jk/2019_07_Hs_Sheppard_IPF_ILD_scRNAseq/receptor_ligand.csv', header = T)
reclig <- subset(reclig, reclig$Evidence %in% c('literature_supported', 'putative'))
reclig <- reclig[,c(2, 3)]
colnames(reclig) <- c('ligand', 'receptor')

ciliated_highcells <- row.names(subset(gse135893_ipf@meta.data, ciliated=='Ciliated_high'))
ciliated_lowcells <- row.names(subset(gse135893_ipf@meta.data, ciliated=='Ciliated_low'))
ciliated_highmatrix <- gse135893_ipfmatrix[,ciliated_highcells]
ciliated_lowmatrix <- gse135893_ipfmatrix[,ciliated_lowcells]
library(psych)
ciliated_highidents <- gse135893_ipfidents[ciliated_highcells,]
ciliated_lowidents <- gse135893_ipfidents[ciliated_lowcells,]
ciliated_high_miner <- PyMINEr(ciliated_highmatrix, reclig, ciliated_highidents$celltype, 'logcounts')
ciliated_high_miner$meanz <- apply(ciliated_high_miner[,5:6], 1, geometric.mean)
ciliated_low_miner <- PyMINEr(ciliated_lowmatrix, reclig, ciliated_lowidents$celltype, 'logcounts')
ciliated_low_miner$meanz <- apply(ciliated_low_miner[,5:6], 1, geometric.mean)
save(ciliated_high_miner, ciliated_low_miner, file = 'gse135893_pyminer_seuratnorm_ipf_cilitedsubsets.RData')

###########ciliated low mac (condition 1, Figure 7B)################
ciliated_low_miner_mac <- ciliated_low_miner[c(grep('Mac', ciliated_low_miner$celltype1), grep('Mono', ciliated_low_miner$celltype1)),]
ciliated_low_miner_mac <- subset(ciliated_low_miner_mac, z1>quantile(ciliated_low_miner_mac$z1, 0.95))
ciliated_low_miner_mac <- subset(ciliated_low_miner_mac, z2>quantile(ciliated_low_miner_mac$z2, 0.95))

grid_col_ligand =rand_color(length(unique(ciliated_low_miner_mac$ligand)))
names(grid_col_ligand) <- unique(ciliated_low_miner_mac$ligand)
grid_col_target = rand_color(length(unique(ciliated_low_miner_mac$receptor)))
names(grid_col_target) <- unique(ciliated_low_miner_mac$receptor)

grid_col_tbl_ligand = tibble(ligand = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_ligand
grid_col_tbl_target = tibble(receptor = grid_col_target %>% names(), color_target_type = grid_col_target)
grid_col_tbl_target
circos_links <- ciliated_low_miner_mac[,c(1:4, 7)]
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand, by = 'ligand') %>% inner_join(grid_col_tbl_target, by = 'receptor')
links_circle = circos_links %>% select(ligand,receptor, meanz)
#circos_links = circos_links %>% mutate(ligand = paste(ligand," "))

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparent ~ ligand-target potential score
transparency = circos_links %>% mutate(meanz =(meanz-min(meanz))/(max(meanz)-min(meanz))) %>% mutate(transparency = 1-meanz) %>% .$transparency

circos_links$condensedcelltype <- plyr::mapvalues(circos_links$celltype2, from = sort(unique(circos_links$celltype2)),
                                                  to = c('AT1', rep('B', 2), 'Basal', 'Ciliated', 'Club', rep('DC', 3),
                                                         'Fibroblast', 'Goblet', 'Endothelial', rep('Mac', 3), 'Fibroblast',
                                                         'Epithelial', 'NK', 'Sm', 'T', rep('Endothelial', 5)))
receptoroccurence <- as.matrix.data.frame(table(circos_links$condensedcelltype, circos_links$receptor))
rownames(receptoroccurence) <- sort(unique(circos_links$condensedcelltype))
colnames(receptoroccurence) <- sort(unique(circos_links$receptor))
receptoroccurence <- as.data.frame(receptoroccurence)
receptoroccurence_total <- apply(receptoroccurence, 2, function(k) sum(k>0))
condensedcelltype_links <- vector()
for (i in 1:nrow(circos_links)){
  if (receptoroccurence_total[circos_links$receptor[i]]==1){
    condensedcelltype_links[i] <- circos_links$condensedcelltype[i]
  }else{
    condensedcelltype_links[i] <- 'mixed'
  }
}
circos_links$targetcelltype <- condensedcelltype_links
circos_links_mixed <- subset(circos_links, circos_links$targetcelltype=='mixed')
circos_links_unique <- subset(circos_links, circos_links$targetcelltype!='mixed')
circos_links_unique <- circos_links_unique[order(circos_links_unique$targetcelltype),]
circos_links_reorg <- as.data.frame(rbind(circos_links_unique, circos_links_mixed))

target_order = circos_links_reorg$receptor %>% unique()
ligand_order = circos_links_reorg[order(circos_links_reorg$celltype1),]$ligand %>% unique() %>% sort()
order = c(ligand_order,target_order)
#Figure 7B-7C
circos.clear()
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$meanz,annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.075), list(track.height = 0.2)))
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) #
highlight.sector(circos_links_reorg$ligand, track.index = 1, col = "red", font = 2, facing = 'bending.outside', text.vjust = '0.5mm',
                 text = "Ligands produced by macrophages", cex = 0.6, text.col = "white", niceFacing = TRUE)
lapply(unique(circos_links_reorg$targetcelltype), 
       function(k) highlight.sector(subset(circos_links_reorg, targetcelltype==k)$receptor, track.index = 1, col = rand_color(length(unique(circos_links_reorg$targetcelltype))), font = 2, facing = 'bending.inside', text.vjust = '0.5mm',
                                    text = k, cex = 0.6, text.col = "white", niceFacing = TRUE))

##############ciliated high mac (condition 2, Figure 7C)#####################
ciliated_high_miner_mac <- ciliated_high_miner[c(grep('Mac', ciliated_high_miner$celltype1), grep('Mono', ciliated_high_miner$celltype1)),]
ciliated_high_miner_mac <- subset(ciliated_high_miner_mac, z1>quantile(ciliated_high_miner_mac$z1, 0.9))
ciliated_high_miner_mac <- subset(ciliated_high_miner_mac, z2>quantile(ciliated_high_miner_mac$z2, 0.9))

grid_col_ligand =rand_color(length(unique(ciliated_high_miner_mac$ligand)))
names(grid_col_ligand) <- unique(ciliated_high_miner_mac$ligand)
grid_col_target = rand_color(length(unique(ciliated_high_miner_mac$receptor)))
names(grid_col_target) <- unique(ciliated_high_miner_mac$receptor)

grid_col_tbl_ligand = tibble(ligand = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_ligand
grid_col_tbl_target = tibble(receptor = grid_col_target %>% names(), color_target_type = grid_col_target)
grid_col_tbl_target
circos_links <- ciliated_high_miner_mac[,c(1:4, 7)]
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand, by = 'ligand') %>% inner_join(grid_col_tbl_target, by = 'receptor')
links_circle = circos_links %>% select(ligand,receptor, meanz)
#circos_links = circos_links %>% mutate(ligand = paste(ligand," "))

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparent ~ ligand-target potential score
transparency = circos_links %>% mutate(meanz =(meanz-min(meanz))/(max(meanz)-min(meanz))) %>% mutate(transparency = 1-meanz) %>% .$transparency

circos_links$condensedcelltype <- plyr::mapvalues(circos_links$celltype2, from = sort(unique(circos_links$celltype2)),
                                                  to = c('Basal', 'Ciliated', 'Club', 'Fibroblast',
                                                         'Endothelial', rep('Mac', 5), 'Fibroblast',
                                                         'Epithelial', 'NK', 'Pericyte', rep('T', 2), 'Endothelial'))
receptoroccurence <- as.matrix.data.frame(table(circos_links$condensedcelltype, circos_links$receptor))
rownames(receptoroccurence) <- sort(unique(circos_links$condensedcelltype))
colnames(receptoroccurence) <- sort(unique(circos_links$receptor))
receptoroccurence <- as.data.frame(receptoroccurence)
receptoroccurence_total <- apply(receptoroccurence, 2, function(k) sum(k>0))
condensedcelltype_links <- vector()
for (i in 1:nrow(circos_links)){
  if (receptoroccurence_total[circos_links$receptor[i]]==1){
    condensedcelltype_links[i] <- circos_links$condensedcelltype[i]
  }else{
    condensedcelltype_links[i] <- 'mixed'
  }
}
circos_links$targetcelltype <- condensedcelltype_links
circos_links_mixed <- subset(circos_links, circos_links$targetcelltype=='mixed')
circos_links_unique <- subset(circos_links, circos_links$targetcelltype!='mixed')
circos_links_unique <- circos_links_unique[order(circos_links_unique$targetcelltype),]
circos_links_reorg <- as.data.frame(rbind(circos_links_unique, circos_links_mixed))

target_order = circos_links_reorg$receptor %>% unique()
ligand_order = circos_links_reorg[order(circos_links_reorg$celltype1),]$ligand %>% unique() %>% sort()
order = c(ligand_order,target_order)

#Figure 7B-7C
circos.clear()
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$meanz,annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.075), list(track.height = 0.2)))
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) #


###################ciliated_low (condition 3, Figure 7D)##################
ciliated_low_miner_ciliated <- ciliated_low_miner[grep('Ciliated', ciliated_low_miner$celltype1),]
ciliated_low_miner_ciliated <- subset(ciliated_low_miner_ciliated, z1>quantile(ciliated_low_miner_ciliated$z1, 0.95))
ciliated_low_miner_ciliated <- subset(ciliated_low_miner_ciliated, z2>quantile(ciliated_low_miner_ciliated$z2, 0.95))

circos.clear()
circos_links <- ciliated_low_miner_ciliated[,c(1:4, 7)]
grid_col_ligand =rand_color(length(unique(ciliated_low_miner_ciliated$ligand)))
names(grid_col_ligand) <- unique(ciliated_low_miner_ciliated$ligand)
grid_col_ligand
grid_col_target = rand_color(length(unique(ciliated_low_miner_ciliated$receptor)))
names(grid_col_target) <- unique(ciliated_low_miner_ciliated$receptor)

grid_col_tbl_ligand = tibble(ligand = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_ligand
grid_col_tbl_target = tibble(receptor = grid_col_target %>% names(), color_target_type = grid_col_target)
grid_col_tbl_target

circos_links = circos_links %>% inner_join(grid_col_tbl_ligand, by = 'ligand') %>% inner_join(grid_col_tbl_target, by = 'receptor')
links_circle = circos_links %>% select(ligand,receptor, meanz)
#circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparent ~ ligand-target potential score
transparency = circos_links %>% mutate(meanz =(meanz-min(meanz))/(max(meanz)-min(meanz))) %>% mutate(transparency = 1-meanz) %>% .$transparency

circos_links$condensedcelltype <- plyr::mapvalues(circos_links$celltype2, from = sort(unique(circos_links$celltype2)),
                                                  to = c('Basal', 'Ciliated', 'Club', 'Fibroblast', 'Goblet', 'Endothelial', 
                                                         rep('Mac', 2), 'Fibroblast', 'Pericyte', 'Sm', rep('T', 2), rep('Endothelial', 4)))
receptoroccurence <- as.matrix.data.frame(table(circos_links$condensedcelltype, circos_links$receptor))
rownames(receptoroccurence) <- sort(unique(circos_links$condensedcelltype))
colnames(receptoroccurence) <- sort(unique(circos_links$receptor))
receptoroccurence <- as.data.frame(receptoroccurence)
receptoroccurence_total <- apply(receptoroccurence, 2, function(k) sum(k>0))
condensedcelltype_links <- vector()
for (i in 1:nrow(circos_links)){
  if (receptoroccurence_total[circos_links$receptor[i]]==1){
    condensedcelltype_links[i] <- circos_links$condensedcelltype[i]
  }else{
    condensedcelltype_links[i] <- 'mixed'
  }
}
circos_links$targetcelltype <- condensedcelltype_links
circos_links_mixed <- subset(circos_links, circos_links$targetcelltype=='mixed')
circos_links_unique <- subset(circos_links, circos_links$targetcelltype!='mixed')
circos_links_unique <- circos_links_unique[order(circos_links_unique$targetcelltype),]
circos_links_reorg <- as.data.frame(rbind(circos_links_unique, circos_links_mixed))

target_order = circos_links_reorg$receptor %>% unique()
ligand_order = circos_links_reorg[order(circos_links_reorg$celltype1),]$ligand %>% unique() %>% sort()
order = c(ligand_order,target_order)
#Figure 7D-7E
circos.clear()
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$meanz,annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.075), list(track.height = 0.2)))
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) #
highlight.sector(circos_links_reorg$ligand, track.index = 1, col = "red", font = 2, facing = 'bending.outside', text.vjust = '0.5mm',
                 text = "Ligands produced by ciliated epithelial cells", cex = 0.6, text.col = "white", niceFacing = TRUE)
lapply(unique(circos_links_reorg$targetcelltype), 
       function(k) highlight.sector(subset(circos_links_reorg, targetcelltype==k)$receptor, track.index = 1, col = rand_color(length(unique(circos_links_reorg$targetcelltype))), font = 2, facing = 'bending.inside', text.vjust = '0.5mm',
                                    text = k, cex = 0.6, text.col = "white", niceFacing = TRUE))

###############ciliated_high (condition 4, Figure 7E)#############
ciliated_high_miner_ciliated <- ciliated_high_miner[grep('Ciliated', ciliated_high_miner$celltype1),]
ciliated_high_miner_ciliated <- subset(ciliated_high_miner_ciliated, z1>quantile(ciliated_high_miner_ciliated$z1, 0.95))
ciliated_high_miner_ciliated <- subset(ciliated_high_miner_ciliated, z2>quantile(ciliated_high_miner_ciliated$z2, 0.95))

circos.clear()
circos_links <- ciliated_high_miner_ciliated[,c(1:4, 7)]
grid_col_ligand = rand_color(length(unique(ciliated_high_miner_ciliated$ligand)))
names(grid_col_ligand) <- unique(ciliated_high_miner_ciliated$ligand)
grid_col_target = rand_color(length(unique(ciliated_high_miner_ciliated$receptor)))
names(grid_col_target) <- unique(ciliated_high_miner_ciliated$receptor)

grid_col_tbl_ligand = tibble(ligand = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(receptor = grid_col_target %>% names(), color_target_type = grid_col_target)

#circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand, by = 'ligand') %>% inner_join(grid_col_tbl_target, by = 'receptor')
links_circle = circos_links %>% select(ligand,receptor, meanz)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparent ~ ligand-target potential score
transparency = circos_links %>% mutate(meanz =(meanz-min(meanz))/(max(meanz)-min(meanz))) %>% mutate(transparency = 1-meanz) %>% .$transparency

circos_links$condensedcelltype <- plyr::mapvalues(circos_links$celltype2, from = sort(unique(circos_links$celltype2)),
                                                  to = c('Ciliated', 'Club', 'DC', 'Epithelial', 'T', rep('Endothelial', 3)))
receptoroccurence <- as.matrix.data.frame(table(circos_links$condensedcelltype, circos_links$receptor))
rownames(receptoroccurence) <- sort(unique(circos_links$condensedcelltype))
colnames(receptoroccurence) <- sort(unique(circos_links$receptor))
receptoroccurence <- as.data.frame(receptoroccurence)
receptoroccurence_total <- apply(receptoroccurence, 2, function(k) sum(k>0))
condensedcelltype_links <- vector()
for (i in 1:nrow(circos_links)){
  if (receptoroccurence_total[circos_links$receptor[i]]==1){
    condensedcelltype_links[i] <- circos_links$condensedcelltype[i]
  }else{
    condensedcelltype_links[i] <- 'mixed'
  }
}
circos_links$targetcelltype <- condensedcelltype_links
circos_links_mixed <- subset(circos_links, circos_links$targetcelltype=='mixed')
circos_links_unique <- subset(circos_links, circos_links$targetcelltype!='mixed')
circos_links_unique <- circos_links_unique[order(circos_links_unique$targetcelltype),]
circos_links_reorg <- as.data.frame(rbind(circos_links_unique, circos_links_mixed))

target_order = circos_links_reorg$receptor %>% unique()
ligand_order = circos_links_reorg[order(circos_links_reorg$celltype1),]$ligand %>% unique() %>% sort()
order = c(ligand_order,target_order)

#chord diagram
#process below repeated for all four conditions above, producing 4 different diagrams
#Figure 7B-7C
circos.clear()
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$meanz,annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.075), list(track.height = 0.2)))
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) #
highlight.sector(circos_links_reorg$ligand, track.index = 1, col = "red", font = 2, facing = 'bending.outside', text.vjust = '0.5mm',
                 text = "Ligands produced by macrophages", cex = 0.6, text.col = "white", niceFacing = TRUE)
lapply(unique(circos_links_reorg$targetcelltype), 
       function(k) highlight.sector(subset(circos_links_reorg, targetcelltype==k)$receptor, track.index = 1, col = rand_color(length(unique(circos_links_reorg$targetcelltype))), font = 2, facing = 'bending.inside', text.vjust = '0.5mm',
                                    text = k, cex = 0.6, text.col = "white", niceFacing = TRUE))

#Figure 7D-7E
circos.clear()
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$meanz,annotationTrack = "grid", 
             preAllocateTracks = list(list(track.height = 0.075), list(track.height = 0.2)))
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) #
highlight.sector(circos_links_reorg$ligand, track.index = 1, col = "red", font = 2, facing = 'bending.outside', text.vjust = '0.5mm',
                 text = "Ligands produced by ciliated epithelial cells", cex = 0.6, text.col = "white", niceFacing = TRUE)
lapply(unique(circos_links_reorg$targetcelltype), 
       function(k) highlight.sector(subset(circos_links_reorg, targetcelltype==k)$receptor, track.index = 1, col = rand_color(length(unique(circos_links_reorg$targetcelltype))), font = 2, facing = 'bending.inside', text.vjust = '0.5mm',
                                    text = k, cex = 0.6, text.col = "white", niceFacing = TRUE))

#nichenet output
library(nichenetr)

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = gse135893_ipf, 
  receiver = c('SPP1pos_macs_0', 'Ciliated_1', 'C1QA_mac_2', 'Ciliated_3', 'AT1_4', 'AT2_5', 'C1QA_mac_6', 
               'ACKR1pos_endo_7', 'Monocytes_8', 'AT1_9', 'Th_10', 'AT1_11', 'Macs_12', 'Tc_13', 'HAS1_fibro_14',
               'Diff_ciliated_15', 'ACKR1neg_endo_16', 'Fibroblasts_17', 'Prolif_macs_18', 'Lymph_endo_19',
               'Sm_20', 'Bcells_21', 'Macs_22', 'PC_23', 'AT2_24', 'MC_25', 'AT1_26', 'Macs_27', 'Ciliated_28', 
               'Fibroblast_29'), 
  condition_colname = "ciliated", condition_oi = "Ciliated_high", condition_reference = "Ciliated_low", 
  sender = c('SPP1pos_macs_0', 'Ciliated_1', 'C1QA_mac_2', 'Ciliated_3', 'AT1_4', 'AT2_5', 'C1QA_mac_6', 
             'ACKR1pos_endo_7', 'Monocytes_8', 'AT1_9', 'Th_10', 'AT1_11', 'Macs_12', 'Tc_13', 'HAS1_fibro_14',
             'Diff_ciliated_15', 'ACKR1neg_endo_16', 'Fibroblasts_17', 'Prolif_macs_18', 'Lymph_endo_19',
             'Sm_20', 'Bcells_21', 'Macs_22', 'PC_23', 'AT2_24', 'MC_25', 'AT1_26', 'Macs_27', 'Ciliated_28', 
             'Fibroblast_29'), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human",
  filter_top_ligands = FALSE)
#Supplementary Figure 9:
DotPlot(gse135893_ipf, features = nichenet_output$top_ligands[1:30] %>% rev(), cols = "RdYlBu") + RotatedAxis()
save(nichenet_output, file = 'gse135893_nichenet_ipfonly_ciliated_seuratnorm.RData')
