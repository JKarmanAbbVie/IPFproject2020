###############Pathway enrichment tables###################
#Table 2A and 2B are generated in IPA (Qiagen). Code for these tables cannot be posted here.

gene1list <- row.names(tT_gse47460_1vshealthy_sig1)
gene1entrez <- bitr(gene1list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
gene1listentrez <- merge(gene1entrez, tT_gse47460_1vshealthy_sig1, by.x = 'SYMBOL', by.y='gene')
gene1listentrez <- gene1listentrez[order(gene1listentrez$logFC, decreasing = T),]
gene1list_vector <- gene1listentrez[,3]
names(gene1list_vector) <- gene1listentrez$ENTREZID
gene1list_go <- gseGO(geneList = gene1list_vector, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, 
                      minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.2, verbose = TRUE)
gene1list_kegg <- gseKEGG(geneList = gene1list_vector, organism = 'hsa', nPerm = 1000, minGSSize = 15, 
                          pvalueCutoff = 0.2, verbose = TRUE)
gene1list_goresults <- gene1list_go@result
gene1list_keggresults <- gene1list_kegg@result

gene1list_reactome <- gsePathway(gene1list_vector, nPerm=10000, minGSSize = 15,
                                 pvalueCutoff=0.2, pAdjustMethod="BH", verbose=TRUE)
gene1list_reactome_results <- gene1list_reactome@result
#Table 2C
gene1list_reactome_results <- gene1list_reactome_results[order(gene1list_reactome_results$NES, decreasing = T),]
gene1list_reactome_results_up <- gene1list_reactome_results[which(gene1list_reactome_results$NES>0),]
write.xlsx(gene1list_reactome_results_up, file = 'gse47460_subset1vshealthy_reactome_results_up_081720.xlsx')

gene2list <- row.names(tT_gse47460_2vshealthy_sig1)
gene2entrez <- bitr(gene2list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
gene2listentrez <- merge(gene2entrez, tT_gse47460_2vshealthy_sig1, by.x = 'SYMBOL', by.y='gene')
gene2listentrez <- gene2listentrez[order(gene2listentrez$logFC, decreasing = T),]
gene2list_vector <- gene2listentrez[,3]
names(gene2list_vector) <- gene2listentrez$ENTREZID
gene2list_go <- gseGO(geneList = gene2list_vector, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, 
                      minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.2, verbose = TRUE)
gene2list_kegg <- gseKEGG(geneList = gene2list_vector, organism = 'hsa', nPerm = 1000, minGSSize = 15, 
                          pvalueCutoff = 0.2, verbose = TRUE)
gene2list_goresults <- gene2list_go@result
gene2list_keggresults <- gene2list_kegg@result
gene2list_reactome <- gsePathway(gene2list_vector, nPerm=10000, minGSSize = 15,
                                 pvalueCutoff=0.2, pAdjustMethod="BH", verbose=TRUE)
gene2list_reactome_results <- gene2list_reactome@result
gene2list_reactome_results <- gene2list_reactome_results[order(gene2list_reactome_results$NES, decreasing = T),]
#Table 2D
gene2list_reactome_results_up <- gene2list_reactome_results[which(gene2list_reactome_results$NES>0),]
write.xlsx(gene2list_reactome_results_up, file = 'gse47460_subset2vshealthy_reactome_results_up_081720.xlsx')
emapplot(gene2list_reactome)
viewPathway("Extracellular matrix organization", readable=TRUE, foldChange=gene1list_vector)
gseaplot(gene1list_reactome, geneSetID = "R-HSA-1474244")
viewPathway("Extracellular matrix organization", readable=TRUE, foldChange=gene2list_vector)
gseaplot(gene2list_reactome, geneSetID = "R-HSA-1474244")

#Figures 1D and 2A, Supplementary Figures 2 and 4: same code below repeated for each Figure, with appropriate 'genes' vector used
genes <- c('RPGRIP1', 'DNAH6', 'DNAH7', 'DNAI1', 'MUC5B')

rm(list=ls(pattern = 'Figure1_'))
rm(figure1_plots)
for (i in genes){
  assign(paste0('Figure1_', i), ggplot(gse47460_phenom_consensus, aes_string('consensusclass', i, fill = 'consensusclass')) + 
           geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
           scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "none"))
}
figure1_plots <- mget(ls(pattern = 'Figure1_'))
rm(list=ls(pattern = 'Figure1_'))
do.call(plot_grid, figure1_plots)

#Figure 1E:
gse32537_reported <- read.csv('gse32537_reported_difflist.csv', header = T)
gse32537_reported_genes <- unique(gse32537_reported$Gene)
gse32537_reported_bygene <- lapply(gse32537_reported_genes, function(k) subset(gse32537_reported, Gene==k)$logFC)
names(gse32537_reported_bygene) <- gse32537_reported_genes
gse32537_reported_geneavs <- lapply(gse32537_reported_bygene, function(k) mean(k)*(-1))
gse32537_reported2 <- as.data.frame(unlist(gse32537_reported_geneavs))
tT_gse47460_tT_gse32537_reported2 <- merge(tT_gse47460, gse32537_reported2, by.x = 'gene', by.y = 0)
rownames(tT_gse47460_tT_gse32537_reported2) <- tT_gse47460_tT_gse32537_reported2$gene
colnames(tT_gse47460_tT_gse32537_reported2)[c(2,8)] <- c('logFC_GSE47460', 'logFC_GSE32537')
#Figure 1E below:
ggplot(tT_gse47460_tT_gse32537_reported2, aes(x=logFC_GSE47460, y=logFC_GSE32537)) + 
  geom_point(shape=18, color="blue")+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="blue")
cor.test(tT_gse47460_tT_gse32537_reported2$logFC_GSE47460, tT_gse47460_tT_gse32537_reported2$logFC_GSE32537)

#Figure 2B:
gse47460_phenom_consensus_dlco <- subset(gse47460_phenom_consensus, dlco!='no_value')
gse47460_phenom_consensus_dlco$dlco <- as.numeric(gse47460_phenom_consensus_dlco$dlco)
dunnTest(dlco~consensusclass, 
         data = gse47460_phenom_consensus_dlco)
TukeyHSD(aov(dlco~consensusclass, data = gse47460_phenom_consensus_dlco))
gg_dlco <- ggplot(gse47460_phenom_consensus_dlco, aes(consensusclass, dlco, fill = consensusclass)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=nicecolors) + theme(legend.position = "none")
nicecolors <- c("#999999", "#E69F00", "#56B4E9")
gse47460_phenom_consensus_fvc_pre <- subset(gse47460_phenom_consensus, fvc_pre!='no_value')
gse47460_phenom_consensus_fvc_pre$fvc_pre <- as.numeric(gse47460_phenom_consensus_fvc_pre$fvc_pre)
dunnTest(fvc_pre~consensusclass, data = gse47460_phenom_consensus_fvc_pre)
gg_fev <- ggplot(gse47460_phenom_consensus_fvc_pre, aes(consensusclass, fvc_pre, fill = consensusclass)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=nicecolors) + theme(legend.position = "none")

gse47460_phenom_consensus_fev1_pre <- subset(gse47460_phenom_consensus, fev1_pre!='no_value')
gse47460_phenom_consensus_fev1_pre$fev1_pre <- as.numeric(gse47460_phenom_consensus_fev1_pre$fev1_pre)
dunnTest(fev1_pre~consensusclass, data = gse47460_phenom_consensus_fev1_pre)
gg_fvc <- ggplot(gse47460_phenom_consensus_fev1_pre, aes(consensusclass, fev1_pre, fill = consensusclass)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=nicecolors) + theme(legend.position = "none")

plot_grid(gg_dlco, gg_fev, gg_fvc)