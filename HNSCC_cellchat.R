suppressPackageStartupMessages({
	library(Seurat)
	library(CellChat)
	library(patchwork)
	library(ggalluvial)
	library(NMF)
	library(ggplot2)
	library(future)
})


seurat_ob = readRDS("seurat.rds")
assay = "RNA"

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
PPI.species <- PPI.human
ident2use <- "celltype"
threads <- 10
output_dir <- getwd()

# =================================================================================
# Load data
# =================================================================================
data.input <- GetAssayData(seurat_ob, assay = assay, slot = "data") # normalized data matrix
seurat_ob = SetIdent(seurat_ob, value = ident2use)
labels <- Idents(seurat_ob)
meta <- data.frame(group = labels, row.names = names(labels))
# Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

### Preprocessing
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = threads) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.species)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pdf(file.path(output_dir, 'interaction_number_network.pdf'))
par(mfrow = c(1,1), xpd=T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, arrow.size=1, title.name = "Number of interactions")
dev.off()
interaction_number_data = cbind(rownames(cellchat@net$count), cellchat@net$count)
colnames(interaction_number_data) = c("ligand_cell/receptor_cell",colnames(interaction_number_data)[-1])
write.table(interaction_number_data , file.path(output_dir,  "interaction_number_network.xls") , sep="\t", quote=F, row.names=F)

pdf(file.path(output_dir, 'interaction_strength_network.pdf'))
par(mfrow = c(1,1), xpd=T)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, arrow.size=1, title.name = "Interaction weights/strength")
dev.off()
interaction_strength_data = cbind(rownames(cellchat@net$weight), cellchat@net$weight)
colnames(interaction_strength_data) = c("ligand_cell/receptor_cell",colnames(interaction_strength_data)[-1])
write.table(interaction_strength_data , file.path(output_dir, "interaction_strength_network.xls") , sep="\t", quote=F, row.names=F)

nrow = length(unique(meta$group))/2
mat <- cellchat@net$count
pdf(file.path(output_dir, "interaction_number_network_grouby.pdf"), height = 7 * (nrow-1))
par(mfrow = c(nrow, 2), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

mat <- cellchat@net$weight
pdf(file.path(output_dir, "interaction_strength_network_grouby.pdf"), height = 7 * (nrow-1))
par(mfrow = c(nrow, 2), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

outdir <- file.path(output_dir, "signaling_pathway_visualize")
dir.create(outdir, recursive = T)
for(pathway in cellchat@netP$pathways){
	pathways.show <- pathway
	# Circle plot
	pdf(file.path(outdir, paste0(pathway,'_Circle_plot.pdf')))
	par(mfrow=c(1,1), xpd=TRUE)
	netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
	dev.off()

	# Chord diagram
	pdf(file.path(outdir, paste0(pathway,'_Chord_diagram.pdf')))
	par(mfrow=c(1,1))
	netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
	dev.off()

	# Heatmap
	pdf(file.path(outdir, paste0(pathway,'_Heatmap.pdf')))
	par(mfrow=c(1,1))
	print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
	dev.off()

	# Access all the signaling pathways showing significant communications
	p = netAnalysis_contribution(cellchat, signaling = pathways.show )
	ggsave(file.path(outdir, paste0(pathway,'_L-R_contribution.pdf')), plot=p)

	p = plotGeneExpression(cellchat, signaling = pathways.show )
	ggsave(file.path(outdir, paste0(pathway,'_GeneExpression.pdf')), plot=p) 

	# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
	pdf(file.path(outdir, paste0(pathway,'_signalingRole.pdf')))
	p = netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 10)
	dev.off()

}

p <- netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = T)
p$communication$pval <- gsub("1", "p > 0.05", p$communication$pval)
p$communication$pval <- gsub("2", "0.01 < p < 0.05", p$communication$pval)
p$communication$pval <- gsub("3", "p < 0.01", p$communication$pval)
write.table(p$communication, file.path(output_dir, "communication.xls"), sep = "\t", row.names = F, quote = F)

pp <- netVisual_bubble(cellchat, remove.isolate = FALSE)
ggsave(file.path(output_dir, "significant_interactions_bubble_plot.pdf"), plot = pp, height = length(unique(p$communication$interaction_name)) / 7 + 1, width = length(unique(p$communication$source.target)) / 7 + 1)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(file.path(output_dir, 'Signaling_role.pdf'),plot=gg1) 
# Signaling role analysis on the cell-cell communication networks of interest
pdf(file.path(output_dir, 'outgoing_signaling_role_heatmap.pdf'))
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()
pdf(file.path(output_dir, "incoming_signaling_role_heatmap.pdf"))
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

saveRDS(cellchat, file.path(output_dir, "cellchat_results.rds"))
