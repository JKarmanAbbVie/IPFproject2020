###############Figure 9.########################
pirfenidone_up <- read.csv('pirfenidone_up_signature_Kwapiszewska_2018_paper_p005_lung_homogenates_lfc141.csv', header = T)
gse47460_pirfenidone_up <- gsva(as.matrix(gse47460_matrix), gset.idx.list = list(pirfenidone_up$Gene))
gse47460_pirfenidone_up_consensus <- merge(t(gse47460_pirfenidone_up), consensusdf2_gse47460, by.x = 0, by.y = 0)
rownames(gse47460_pirfenidone_up_consensus) <- gse47460_pirfenidone_up_consensus$Row.names
gse47460_pirfenidone_up_consensus$Row.names <- NULL
colnames(gse47460_pirfenidone_up_consensus)[1] <- 'Pirfenidone_signature_score'

ggplot(gse47460_pirfenidone_up_consensus, aes(consensusclass, Pirfenidone_signature_score, fill = consensusclass)) + 
  geom_boxplot(width=0.5, show.legend = F) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=nicecolors) + 
  ylab('Signature score') + ggtitle('Pirfenidone signature score') + theme(plot.title = element_text(size = 12), legend.position = "none")
dunnTest(Pirfenidone_signature_score~consensusclass, data = gse47460_pirfenidone_up_consensus)