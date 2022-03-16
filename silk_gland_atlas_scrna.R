library(Seurat)

countMatrixSparse <- Read10X(matrix_path, gene.column = 2)
seurat_ob <- CreateSeuratObject( countMatrixSparse[["Gene Expression"]], names.field = 2, assay = "RNA", names.delim = "-" )

mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(seurat_ob), value = T,perl=T)
raw.counts = GetAssayData(seurat_ob, slot = "counts")
percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
#remove the outliers cells using the linear model
UMIs_per_cell = Matrix::colSums(raw.counts)
genes_per_cell = Matrix::colSums(raw.counts>0) 
df = data.frame(UMIs_per_cell=UMIs_per_cell, genes_per_cell=genes_per_cell)
df = df[order(df$UMIs_per_cell),] 
m <- rlm(genes_per_cell~UMIs_per_cell,data=df) 
pb <- data.frame(predict(m, interval='prediction',level = 1-1e-3, type="response"))
# identifier outliers as having observed genes_per_cell outside the prediction confidence interval
outliers <- rownames(df)[df$genes_per_cell > pb$upr | df$genes_per_cell < pb$lwr]
seurat_ob = subset(seurat_ob, cells = colnames(raw.counts)[!colnames(raw.counts) %in% outliers])
# Filter out low-quality cells
seurat_by_sample =SplitObject(seurat_ob, split.by = "sampleid" )
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    upper_bound <- 10^(mean(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])) + 2*sd(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])))
    lower_bound <- 10^(mean(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])) - 2*sd(log10(object@meta.data[,"nFeature_RNA"][object@meta.data[,"nFeature_RNA"]>0])))
    object = SubsetData(object, subset.name = "nFeature_RNA",
        low.threshold = lower_threshold, high.threshold = upper_threshold)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    upper_bound <- 10^(mean(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])) + 2*sd(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])))
    lower_bound <- 10^(mean(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])) - 2*sd(log10(object@meta.data[,"nCount_RNA"][object@meta.data[,"nCount_RNA"]>0])))
    object = SubsetData(object, subset.name = "nCount_RNA",
        low.threshold = lower_threshold, high.threshold = upper_threshold)
    seurat_by_sample[[idx]] = object
}
for( idx in 1:length(seurat_by_sample) ){
    object = seurat_by_sample[[idx]]
    object = SubsetData(object, subset.name = "percent.mito",
        low.threshold = "-Inf", high.threshold = "0.3")
    seurat_by_sample[[idx]] = object
}
merged_seurat = seurat_by_sample[[1]]   
for( idx in 2:length(seurat_by_sample) ){
    merged_seurat = merge(x = merged_seurat, y = seurat_by_sample[[idx]],
    do.scale = F, do.center = F, do.normalize = F)
}

seurat_ob = merged_seurat
seurat_ob <- NormalizeData(object = seurat_ob,
                            normalization.method = opt$normmeth,scale.factor = 10000)
seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
                    clip.max = "auto", mean.function = "FastExpMean",
                    dispersion.function = "FastLogVMR", num.bin = 20,
                    nfeature = 4000, binning.method = "equal_width" )

seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
                    vars.to.regress = c("nCount_RNA","percent.mito"), verbose = T )

seurat_ob <- RunMnn.Seurat(seurat_ob, features = VariableFeatures(seurat_ob),assay = "RNA",
                                batch = "batchid", npc.use = 10)
seurat_ob = FindNeighbors( seurat_ob, reduction = "mnn", dims = 1:10,
                                features = VariableFeatures(seurat_ob),
                               nn.eps = 0, force.recalc = T, verbose = F)
seurat_ob <- FindClusters(object = seurat_ob, resolution = 0.4, algorithm = 1, verbose = F)
seurat_ob = RunUMAP(seurat_ob,dims = 1:10,verbose = F,seed.use=1,reduction = "mnn", n.components = 2)
# find markers
global_DEGs = FindAllMarkers(object = seurat_ob,only.pos = T,test.use = "bimod", logfc.threshold = 0, min.pct = 0.25)
global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% select( gene, everything())
#top10 markers
topn_markers  = global_DEGs %>% group_by(cluster) %>% 
            arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
            top_n(10,gene_diff)
saveRDS(seurat,"seurat.rds")

#==============================monocle============================#
seurat_ob = readRDS( "seurat.rds" )
gbm_cds = as.CellDataSet( seurat_ob )
#Estimate size factors and dispersions
gbm_cds = estimateSizeFactors(gbm_cds)
gbm_cds = estimateDispersions(gbm_cds)
#Filtering low-quality cells
gbm_cds = detectGenes(gbm_cds,min_expr = 1)
#keep only the genes expressed in at least 10 cells of the data set.
expressed_genes = row.names(subset(fData(gbm_cds),num_cells_expressed>10))
pData(gbm_cds)$Total_mRNA = Matrix::colSums(exprs(gbm_cds))
#Find ordering genes
clustering_DEGs = differentialGeneTest(gbm_cds,fullModelFormulaStr = ~clusters)
featureData(gbm_cds)@data[rownames(clustering_DEGs),"pval"]=clustering_DEGs$pval
featureData(gbm_cds)@data[rownames(clustering_DEGs),"qval"]=clustering_DEGs$qval
ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
gbm_cds = setOrderingFilter(gbm_cds,ordering_genes = ordering_genes)
gbm_cds = reduceDimension(gbm_cds,max_components = 2,verbose = T,check_duplicates = F)
gbm_cds = orderCells(gbm_cds,reverse = F)
saveRDSMC(gbm_cds,"pseudotime_results.rds")

