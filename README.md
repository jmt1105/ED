N16 <- Read10X(data.dir ="D:/For research/R work package/单细胞测序分析/CON VED 单细胞测序/N16_N17_result.tar/N16_N17_result/N16/outs/filtered_feature_bc_matrix") 
N17 <- Read10X(data.dir ="D:/For research/R work package/单细胞测序分析/CON VED 单细胞测序/N16_N17_result.tar/N16_N17_result/N17/outs/filtered_feature_bc_matrix") 
Con <- CreateSeuratObject(counts = N16, project = "Con",min.cells = 3, min.features = 200)
VED <- CreateSeuratObject(counts = N17, project = "VED",min.cells = 3, min.features = 200)
Con <- NormalizeData(Con)
Con <- FindVariableFeatures(Con, nfeatures = 4000)
VED <- NormalizeData(VED)
VED <- FindVariableFeatures(VED, nfeatures = 4000)
sampleList <- list(Con, VED)
scedata <- FindIntegrationAnchors(object.list = sampleList, dims = 1:20)
scedata <- IntegrateData(anchorset = scedata, dims = 1:20)
scedata <- ScaleData(scedata, vars.to.regress = "nFeature_RNA", verbose = FALSE)
scedata <- RunPCA(scedata, npcs = 20, verbose = FALSE)
scedata <- FindNeighbors(scedata, reduction = "pca", dims = 1:20)
scedata <- FindClusters(scedata, 
                        resolution = 0.5)
scedata <-RunTSNE(scedata, reduction = "pca", dims = 1:20) 
DefaultAssay(scedata) <- "RNA"
all.markers  <- FindAllMarkers(scedata,
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75,
                               grouping.var = "seurat_clusters", verbose = FALSE)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]#根据基因标注细胞
SMC <- subset(scedata,celltype=="SMC")
diff_gene <- FindMarkers(SMC, group.by = "orig.ident", ident.1 = "VED", ident.2 = "Con", verbose = FALSE)
diff_gene <- read.csv(file="diff_gene_SMC.CSV")
rt=diff_gene
genes=as.vector(rt[,1])	
entrezIDs=mget(genes, org.Rn.egSYMBOL2EG, ifnotfound=NA)	
View(entrezIDs)
entrezIDs=as.character(entrezIDs)	
rt=cbind(rt,entrezID=entrezIDs)		
rt=rt[is.na(rt[,"entrezID"])==F,]	
gene=rt$entrezID
GOpathway=enrichGO(gene = gene,
                   OrgDb = org.Rn.eg.db, 
                   pvalueCutoff =1,	
                   qvalueCutoff = 1,	
                   ont="all",	
                   readable =T)	
GOpathway1=as.data.frame(GOpathway)
write.table(GOpathway1,file="GOpathway1.txt",sep="\t",quote=F,row.names = F)
keggpathway1 <- enrichKEGG(rt$entrezID, organism = 'rno', keyType = 'kegg', pvalueCutoff = 0.5,
                           pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 1, use_internal_data = FALSE)        
keggpathway=as.data.frame(keggpathway1)
write.table(keggpathway,"keggpathway.txt") 
