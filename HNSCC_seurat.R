# R version 3.6.1
library(Seurat)
library(DoubletFinder)
library(tictoc)
library(presto)
library(future)
library(future.apply)
library(dplyr)
library(tibble)

countMatrixSparse <- Read10X(matrix_path, gene.column = 2)
seurat_ob <- CreateSeuratObject( countMatrixSparse[["Gene Expression"]], names.field = 2, assay = "RNA", names.delim = "-" )

# mito genes
mito.genes <- c("MT-ATP6","MT-ATP8","MT-CO1","","MT-CO2","","MT-CO3","","MT-CYB","","MT-ND1","","MT-ND2","","MT-ND3","","MT-ND4","","MT-ND4L","MT-ND5","","MT-ND6")
raw.counts = GetAssayData(seurat_ob, slot = "counts")
percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
# HB genes
hb.genes <- c("HBB","HBD","HBG1","HBG2","HBE1","HBZ","HBM","HBA2","HBA1","HBQ1")
percent.hb <- Matrix::colSums(raw.counts[hb.genes, ])/Matrix::colSums(raw.counts)
seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.hb, col.name = "percent.hb")
# log10GenesPerUMI
data_ob@meta.data$log10GenesPerUMI <- log10(data_ob@meta.data$nFeature_RNA)/log10(data_ob@meta.data$nCount_RNA)

# Filter out low-quality cells
seurat_by_sample =SplitObject(seurat_ob, split.by = "sampleid" )
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "nFeature_RNA", low.threshold = 200)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "nCount_RNA", low.threshold = 1000)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "percent.mito",
        low.threshold = "-Inf", high.threshold = "0.15")
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "percent.hb",
        low.threshold = "-Inf", high.threshold = "0.05")
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "log10GenesPerUMI", low.threshold = "0.7")
    seurat_by_sample[[idx]] = object
}

merged_seurat = seurat_by_sample[[1]]   
for( idx in 2:length(seurat_by_sample) ){
    merged_seurat = merge(x = merged_seurat, y = seurat_by_sample[[idx]],
    do.scale = F, do.center = F, do.normalize = F)
}

seurat_ob = merged_seurat

print("removing doublet using DoubletFinder")
RemoveDoublets <- function(
  object,
  identity,
  doublet.rate,
  pN=0.25,
  PCs=1:30,
  use.SCT=FALSE,
  num.cores=1,
  quietly=TRUE
){

  tic("get sweep parameters")
  # calculate parameters
  if (quietly==TRUE){
    invisible(capture.output(sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)))
    invisible(capture.output(sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)))
    ff <- tempfile()
    png(filename=ff)
    invisible(capture.output(bcmvn <-
find.pK(sweep.stats)))
    dev.off()
    unlink(ff)
  }else{
    sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)
    sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
    ff <- tempfile()
    png(filename=ff)
    bcmvn <-
find.pK(sweep.stats)
    dev.off()
    unlink(ff)
  }
  toc()

  # choose parameters
  maxBCmetric    <- max(bcmvn$BCmetric, na.rm = TRUE)
  pK <- as.numeric(as.character(bcmvn[bcmvn$BCmetric==maxBCmetric, ]$pK))

  # compute doublet scores
  tic("Removing doublets")
  annotations    <- object@meta.data$identity  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi       <- round(doublet.rate*length(colnames(x = object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
  seu.scored     <- doubletFinder_v3(object, PCs =PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = use.SCT)
  toc()

  # pick out doublets
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")

  # remove doublets
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  return(list(removed=seu.removed, original=seu.scored))
}

rate <- function(num){
        if (num >10000) rate= 0.1 else rate= num*0.0000077
        return(rate)
}

seurat_ob = SetIdent( seurat_ob, value = "sampleid")
obj = SplitObject(seurat_ob,split.by = "sampleid")
obj_rm=list()
for( i in names(obj)){
	obj[[i]] <- NormalizeData(obj[[i]])
	obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
	obj[[i]] <- ScaleData(obj[[i]])
	obj[[i]] <- RunPCA(obj[[i]])
	obj[[i]] <- RunUMAP(obj[[i]], dims = 1:10)
	obj_rm[[i]] = RemoveDoublets(obj[[i]], doublet.rate=rate(dim(obj[[i]])[2]),  num.cores=4)
}
# metadata + doublet pANN DF.classifications
removed= lapply(obj_rm,FUN=function(x) x = x$removed)
if ( length(removed) > 1 ) {
	seurat_ob = merge(removed[[1]],do.call(c,removed[-1]))
} else {
	seurat_ob = removed[[1]]
}



seurat_ob <- NormalizeData(object = seurat_ob,
                            normalization.method = 'LogNormalize',scale.factor = 10000)
seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
                    clip.max = "auto", mean.function = "FastExpMean",
                    dispersion.function = "FastLogVMR", num.bin = 20,
                    nfeature = 2000, binning.method = "equal_width" )

seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
                    vars.to.regress = c("nCount_RNA","percent.mito"), verbose = T )

seurat_ob <- RunPCA(seurat_ob, features = VariableFeatures(seurat_ob),assay = "RNA")
seurat_ob = FindNeighbors( seurat_ob, reduction = "pca", 
                                features = VariableFeatures(seurat_ob),
                               nn.eps = 0, force.recalc = T, verbose = F)
seurat_ob <- FindClusters(object = seurat_ob, resolution = 0.4, algorithm = 1, verbose = F)
seurat_ob = RunUMAP(seurat_ob,dims = 1:10,verbose = F,seed.use=1,reduction = "pca", n.components = 2)
# find markers
global_DEGs = wilcoxauc(seurat_ob, group_by  = 'clusters', assay = "data")
data = GetAssayData(seurat_ob, slot = "counts")
results = future_lapply(unique(Idents(seurat_ob)), function(idx){
	cells.1 = Cells(subset(seurat_ob, idents = idx))
	cells.2 = Cells(subset(seurat_ob, idents = idx, invert = T))
	deg4idx = global_DEGs %>% filter(group == idx)
	pct.1 = round( Matrix::rowSums(data[deg4idx$feature, cells.1] > 0)/length(cells.1), digits = 3)
	pct.2 = round( Matrix::rowSums(data[deg4idx$feature, cells.2] > 0)/length(cells.2), digits = 3)
	data.alpha <- cbind(pct.1, pct.2)
	colnames(x = data.alpha) <- c("pct.1", "pct.2")
	xe = cbind(deg4idx, data.alpha)
	alpha.max <- apply( data.alpha, MARGIN = 1, FUN = max)
	names(x = alpha.max) <- rownames(data.alpha)
	# alpha.diff <- alpha.max - apply(X = data.alpha, MARGIN = 1, FUN = min)
	features <- names( which(alpha.max > 0.25 ) )
	if (length(x = features) == 0) { stop("No features pass min.diff.pct threshold") }
	alpha = xe %>% filter( feature %in% features)
})
global_DEGs =  do.call(rbind, results)
global_DEGs = global_DEGs %>% filter( logFC > 0 ) %>%
				select( feature, group, logFC, pval,padj,pct.1, pct.2, auc) %>%
				rename(gene = feature, cluster = group,avg_logFC = logFC, p_val_adj = padj, p_val = pval )
global_DEGs = global_DEGs %>%
    dplyr::mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>%
    dplyr::select( gene, everything() )
#top10 markers
topn_markers  = global_DEGs %>% group_by(cluster) %>% 
            arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
            top_n(10,gene_diff)
saveRDS(seurat,"seurat.rds")

