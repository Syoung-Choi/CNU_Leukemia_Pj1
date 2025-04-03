# load packages
library(remotes)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)

#################################
# start using Seurat in R and choose the folder to work on.
library(Seurat)
setwd("./your/dir") 
getwd() 

################################
# load dataasets
AML_List <- list(
  "CR_List" = readRDS("./01_RDS/01.CR_integrated.rds"),
  "ER_List" = readRDS("./01_RDS/02.ER_integrated.rds")
) 


#doing anchors combination of AML_ER and AML_CR
nchors = FindIntegrationAnchors(object.list = AML_List,dims = 1:30,anchor.features=4000)
integrated = IntegrateData(anchorset = anchors, dims = 1:30)
integrated = ScaleData(integrated, verbose = FALSE,features=rownames(integrated))
integrated = RunPCA(integrated, npcs = 50, verbose = FALSE,features=VariableFeatures(object = integrated))
integrated = JackStraw(integrated,num.replicate = 100)
integrated = ScoreJackStraw(integrated,dims=1:20)
integrated = FindNeighbors(integrated, dims = 1:35)
integrated = FindClusters(integrated, resolution = 0.5) 
integrated = RunTSNE(integrated,dims=1:35)
integrated = RunUMAP(integrated,dims=1:38)
integrated$source2 = ifelse(substr(integrated$orig.ident,1,6)=='AML_ER', integrated$orig.ident, 'AML_CR')
integrated$source3 = ifelse(substr(integrated$orig.ident,1,6)=='AML_ER', 'AML_ER', 'AML_CR')


cell_counts <- table(integrated@meta.data$source3)
Control_cells <- cell_counts["AML_CR"]
Treatment_cells <- cell_counts["AML_ER"]
print(paste("Control_cells:",Control_cells))
print(paste("Treatment_cells:",Treatment_cells))

#====================================================================================================

integrated <- readRDS("./RData/03.CRER_intgrated.rds")

DimPlot(integrated, pt.size = 1.5,label = TRUE,label.size = 8) + NoLegend() + ggtitle(label = "Cluster of CR and ER patients")

FeaturePlot(integrated, c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_NCAM1", "rna_FCGR3A", "rna_CD14", "rna_MS4A7", "rna_ITGAX", "rna_CD34", "rna_CD19", "rna_CD99", "rna_CD79A"), pt.size = 0.8, ncol = 3, max.cutoff = 2, min.cutoff = 0)

VlnPlot(integrated, c("rna_HAVCR2"), group.by = "orig.ident", split.by = "source3")

All.markers <- FindAllMarkers(integrated, only.pos = TRUE, test.use = "roc")

All.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1) %>%
  slice_head(n=5) %>%
  ungroup() -> top5 

p1 <- DoHeatmap(integrated, features = top5$gene, slot = "scale.data", group.bar.height = 0.02) + NoLegend()
p1

a1 <- DimPlot(integrated, pt.size = 1.5, label = T, label.size = 8) + NoLegend() + ggtitle(label = "Clusters of CR and ER patients")
a2 <- FeaturePlot(integrated, c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_NCAM1", "rna_FCGR3A", "rna_CD14", "rna_MS4A7", "rna_ITGAX", "rna_CD34", "rna_CD99", "rna_CD19", "rna_CD79A"), pt.size = 0.8, ncol = 3, max.cutoff = 2, min.cutoff = 0)
a1+a2


CRER <- readRDS(file = "./01.RDS/04.CRER_Find_markers.rds")

#======================================================
plot_integrated_Celltype = function(srat){ 
  ##take an integrated Seurat object, plot distributions over org.ident.
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)

count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$source3)
count_mtx <- as.data.frame.matrix(count_table)
count_mtx$cluster <- rownames(count_mtx)
melt_mtx <- melt(count_mtx)
melt_mtx$cluster <- as.factor(melt_mtx$cluster)

cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)), decreasing = T))
cluster_size$cluster <- factor(cluster_size$cluster, levels = sorted_labels)
melt_mtx$cluster <- factor(melt_mtx$cluster, levels = sorted_labels)

colnames(melt_mtx)[2] <- "dataset"

p1 <- ggplot(cluster_size, aes(y = cluster, x = value)) + geom_bar(position = "dodge", stat = "identity", fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cell per cluster, log10 scale") + ylab("")

p2 <- ggplot(melt_mtx,aes(x = cluster, y = value, fill = dataset)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + 
  scale_fill_brewer(palette = "Set2") + ylab("Fraction of cells in each dataset") + xlab("Cell type") + theme(legend.position = "top")

p2+p1 + plot_layout(widths = c(3,1))

}

##==================================================================================================================
plot_integrated_Celltype(CRER)

DimPlot(CRER,pt.size = 1.5, label = T,label.size = 8)+ NoLegend() + ggtitle(label = "Clusters of CR and ER patients")

FeaturePlot(CRER, features = ("rna_HBB"), max.cutoff = 8, pt.size = 1)  

CRER_Tcell <- subset(CRER, rna_CD3D > 1|rna_CD3E > 1|rna_CD3G > 1, ident = c(0,1,2,3,4,5,6,8))
ElbowPlot(CRER)

CRER_Tcell <- FindNeighbors(CRER_Tcell, dims = 1:30)
CRER_Tcell <- FindClusters(CRER_Tcell, resolution = 0.6)
CRER_Tcell <- RunTSNE(CRER_Tcell, dims = 1:30)
CRER_Tcell <- RunUMAP(CRER_Tcell, dims = 1:40)

ER_DNT_MK <- c("rna_CD4", "rna_CD8A", "rna_IKZF2", "rna_CD160", "rna_LYZ", "rna_TRDC", "rna_THEMIS", "rna_KLRF1", "rna_ZBTB16") #MK: marker

a3 <- DimPlot(CRER_Tcell, pt.size = 1.5, label = T, label.size = 8)+ NoLegend()+ ggtitle(label = "Clusters of CR and ER CD3+ Tcell subset")

a4 <- FeaturePlot(CRER_Tcell, features = c("rna_CD4", "rna_CD8A", "rna_IKZF2", "rna_CD160", "rna_TRAC", "rna_TRDC", "rna_THEMIS", "rna_KLRF1", "rna_ZBTB16"), pt.size = 1, min.cutoff = 0, max.cutoff = 2)

a3+a4

saveRDS(object = CRER_Tcell, file = "./RData/05.CRER_Tcell_reclusters.rds")

DimPlot(CRER_Tcell, pt.size = 1.5, label = T, label.size = 8, split.by = "source3")+ NoLegend() + ggtitle(label = "Clusters of CR and ER patients")

FeaturePlot(CRER_Tcell, features = c("rna_CD4", "rna_CD8A", "rna_IKZF2", "rna_CD160"), pt.size = 1, min.cutoff = 0, max.cutoff = 2, ncol = 2)
FeaturePlot(CRER_Tcell, features = c("rna_LYZ", "rna_TRDC", "rna_THEMIS", "rna_KLRF1", "rna_ZBTB16"),  split.by = "source3", pt.size = 1, min.cutoff = 0, max.cutoff = 2, ncol = 4)
FeaturePlot(CRER_Tcell, features = c("rna_TRAC", "rna_TRBC1", "rna_TRDC", "rna_TRGC1"), pt.size = 1, min.cutoff = 0, max.cutoff = 2) # check TCR (a-b-g-d)

VlnPlot(CRER_Tcell, features = c("rna_CD4", "rna_CD8A", "rna_IKZF2", "rna_CD160", "rna_TRAC", "rna_TRDC", "rna_THEMIS", "rna_KLRF1", "rna_ZBTB16", "rna_CD14", "rna_SELL", "rna_CCR7"), pt.size = 0)


new.cluster.ids <- c("DNT cells", "CD8+ Tcells", "CD4+ Tcells", "CD4+ Tcells", "DNT cells", "DNT cells", "CD8+ Tcells", "CD8+ Tcells", "MAIT cells", "MAIT cells")
names(new.cluster.ids) <- levels(CRER_Tcell)
CRER_Tcell <- RenameIdents(CRER_Tcell, new.cluster.ids)
DimPlot(CRER_Tcell, pt.size = 1.5, label = T, label.size = 6)
#levels(CRER_Tcell)

CRER_Tcell <- readRDS("./RData/06.CRER_Tcell_reclusters_annotated.rds")

VlnPlot(CRER_Tcell, features = ("rna_HAVCR2"), split.by = "source3", pt.size = 0.5, cols = c("cyan","goldenrod"))
VlnPlot(CRER_Tcell, features = ("rna_IKZF2"), split.by = "source3", pt.size = 0.5, cols = c("cyan", "goldenrod"))
VlnPlot(CRER_Tcell, features = c("rna_TRAC", "rna_TRBC1", "rna_TRDC", "rna_TRGC1"), pt.size = 0, ncol = 2)

FeaturePlot(CRER_Tcell, features = ("rna_HAVCR2"), split.by = "source3", pt.size = 1)

DotPlot(CRER_Tcell, features = ER_DNT_MK, scale = T, col.min = 0, col.max = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#==============================================================================
#identification of cell clusters
CRER_Allmarkers <- FindAllMarkers(CRER_Tcell, only.pos = T, test.use = "roc")

CRER_Allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1)

CRER_Allmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n=5) %>%
  ungroup() -> top5

DoHeatmap(CRER_Tcell,features = top5$gene, slot = "scale.data", group.bar.height = 0.02) + NoLegend()


#comp_CR vs ER volcano
library(EnhancedVolcano)
library(MAST)

a <- FindMarkers(CRER_Tcell, ident.1 = "DNT cells", subset.ident = "source3", test.use = "MAST")
VlnPlot(CRER_Tcell, features = "rna_GZMB", group.by = "seurat_clusters", split.by = "source3")
CRER_Tcell@meta.data$seurat_clusters

a %>%
  dplyr::filter(avg_log2FC > 1 | avg_log2FC < -1) %>% 
  dplyr::arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 50) %>%
  ungroup() -> Comp_DNT

#a1 <- a%>%
#  arrange(desc(myAUC))%>%
#  slice_head(n = 15) %>%
#  ungroup() -> Comp_DNT using for #a above


Comp_DNT
T1<-subset(CRER_Tcell, ident = "DNT cells")
DoHeatmap(T1, features = c(features = rownames(Comp_DNT)),group.by = "source3",slot = "scale.data") + NoLegend()

title_name <- "DNT CR vs ER" 
head(a,5)
EnhancedVolcano(a,
                lab = rownames(a),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = title_name,
                pCutoff = 10e-100,
                FCcutoff = 2,
                pointSize = 1.5, labCol = "black",
                labSize = 5.0,
                col = c('black', 'black', 'black', "red3"),
                colAlpha = 1)



EnhancedVolcano(a,
                lab = rownames(a),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = title_name,
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.5,labCol = "black",
                labSize = 5.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1) 



keyvals <- ifelse(
  a$avg_log2FC < -2 & a$avg_log2FC < 10e-50, 'blue',
  ifelse(a$avg_log2FC > 2 & a$avg_log2FC < 10e-50, 'red',
         'black')) 



keyvals[is.na(keyvals)] <- 'black' 
names(keyvals)[keyvals == 'red'] <- 'CR_upregulated'
names(keyvals)[keyvals == 'black'] <- 'non'
names(keyvals)[keyvals == 'blue'] <- 'ER_upregulated'

EnhancedVolcano(a,
                lab = rownames(a),
                x = 'avg_log2FC', 
                y = 'p_val_adj',  
               selectLab = rownames(a)[which(names(keyvals) %in% c('ER_upregulated', 'CR_upregulated'))],
               xlab = bquote(~Log[2]~ 'fold change'), 
               title = title_name,
               pCutoff = 10e-50, 
               FCcutoff = 2,
               
               pointSize = 2.5,            
               labSize = 3,                
               colCustom = keyvals,        
               colAlpha = 1,              
               
               legendPosition = 'bottom',
               legendLabSize = 15,
               labFace = 'bold',
               legendIconSize = 5.0,
               
               drawConnectors = T,         
               widthConnectors = 1.0, boxedLabels = T, 
               colConnectors = 'black',
               arrowheads = F,            
               
               gridlines.major = T,        
               gridlines.minor = F,
               
               border = 'partial',         
               borderWidth = 1.5,          
               borderColour = 'black')   



BiocManager::install("msigdbr")
BiocManager::install("fgsea")
remotes::install_github("immunogenomics/presto")

library(msigdbr)    
library(fgsea)      
library(data.table) 
library(presto)      
library(tidyr)      
library(dplyr)      
library(ggplot2)    
library(tibble)     
 
AML_CRER <- subset(CRER_Tcell, idents = c("DNT cells"))
DimPlot(AML_CRER, split.by = "source3")

table(Idents(CRER_Tcell))             
prop.table(table(Idents(CRER_Tcell))) 

Idents(AML_CRER) <- "source3"
table(Idents(AML_CRER))

Assay_title <- "ER vs CR"       
DefaultAssay(AML_CRER) <- "RNA" 
AML_CRER = JoinLayers(AML_CRER) 

Idents(AML_CRER) <- "celltype.condition" #gọi AML_CRER là celltype.condition.
myData.genes<- wilcoxauc(AML_CRER, "source3") 
                                              

#tạo danh sách tập hợp gene phục vụ cho phân tích GSEA, enrichment 
m_df<- msigdbr(species = "Homo sapiens", category = "H") 
                                                         
                                                         
                                                         



fgsea_set<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(m_df)                       
dplyr::count(myData.genes, group) 

myData.genes %>%                        
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n=10)                          

AML_ER.genes <- myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)         

ranks <- deframe(AML_ER.genes)        
head(ranks)
##############################

fgseaRes <- fgsea(fgsea_set, stats = ranks, nperm = 1000) 

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

#NES: thấp lên cao
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()

ggplot(fgseaResTidy %>% 
         dplyr::filter(padj < 0.05) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +  
  geom_col(aes(fill = NES<0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = Assay_title,) +
  theme_minimal()

Pathway <-"HALLMARK_OXIDATIVE_PHOSPHORYLATION" 

plotEnrichment(fgsea_set[[Pathway]],                     
               ranks) + labs(title=Pathway)              
                                                        
class(top5)
class(fgsea_set)

top5[1,]
fgsea_set[Pathway] 
fgsea_set[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]       
                                                        
###########################
AML_SUB= readRDS('./RData/T_NK_subset_V.1.rds')
DNT_ALL <- c("rna_CD160", "rna_CD44", "rna_FCER1G", "rna_ANXA2", "rna_IFNG", "rna_S100A11", "rna_XCL1", "rna_PIM1", "rna_IKZF2", "rna_ATF4")
DotPlot(AML_SUB,features = DNT_ALL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

AML_DNT= readRDS('./RData/T_NK_subset_V.1.rds')

FeaturePlot(AML_DNT, feature = c("rna_HAVCR2", "rna_TIGIT", "rna_CTLA4", "rna_PDCD1", "rna_LAG3"), ncol = 5, split.by = "orig.ident") 

FeaturePlot(AML_DNT, features = c("rna_LGALS9", "rna_PVR", "rna_CD80", "rna_CD274", "rna_CEACAM1"), ncol = 4, split.by ="orig.ident")


DNT_surface <- c("rna_CD160", "rna_CD44", "rna_FCER1G")
DNT_secreted <- c("rna_ANXA2", "rna_IFNG", "rna_S100A11", "rna_XCL1")
DNT_IF <- c("rna_SUB1", "rna_PIM1", "rna_IKZF2", "rna_ATF4")

DoHeatmap(AML_DNT, features = DNT_ALL, group.by = "source3")
DNT_ALL <- gsub("rna_", "", DNT_ALL) 

# 
library(openxlsx)
library(HGNChelper)
library(SeuratData)
library(EnhancedVolcano)

AML_DNT <- subset(CRER_Tcell, idents = "DNT cells")
DimPlot(AML_DNT, reduction = "umap", label = F, split.by = "orig.ident")
#===================================
AML_DNT[["RNA"]]
AML_DNT[["integrated"]]

class(AML_DNT[["RNA"]])
class(AML_DNT[["integrated"]])

names(AML_DNT)


DefaultAssay(AML_DNT) <- "RNA"
AML_DNT_1 = JoinLayers(AML_DNT)

cs <- colSums(GetAssayData(AML_DNT_1,assay="RNA",layer="counts")>0) 
ce <- cs[which(cs>5)]                                               
sf <- cs[names(ce)]                                                 

Idents(AML_DNT)<-"celltype.condition"
myData.genes <-wilcoxauc(AML_DNT_1,"source3") 


m_df<- msigdbr(species = "Homo sapiens", collection = "C2",subcollection = "CP:KEGG_MEDICUS") 

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

head(m_df) 
dplyr::count(myData.genes,group)

myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(logFC),desc(auc)) %>%
  head(n=10)

AML_ER.genes<- myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
rank <- deframe(AML_ER.genes)
head(rank)

################
fgseaRes <-fgsea(fgsea_sets, stats = ranks, nperm = 1000)  

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, nMoreExtreme) %>%
  arrange(padj) %>%
  head()

conflicts_prefer(dplyr::filter) 

ggplot(fgseaResTidy %>% dplyr::filter(padj < 0.05) %>% head(n = 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES<0)) +
  coord_flip() +
  labs(x="pathway", y="Normalized Enrichment Score",
       title = "Hallmark pathways NES from GSEA") + 
  theme_minimal()

a <- "KEGG_LEISHMANIA_INFECTION"
plotEnrichment(fgsea_set[[a]], ranks) + labs(title = a)

VlnPlot(object = CRER_Tcell, features = c("rna_TRDC", "rna_IKZF2", "rna_ZBTB16", "rna_CD28"), pt.size = 0, ncol = 4)
VlnPlot(object = CRER_Tcell, features = ER_DNT_MK, pt.size = 0)

TCR_SET <-c("rna_TRAC", "rna_TRBC1", "rna_TRDC", "rna_TRGC1")                                                                                    #TCRab and TCRgd
Notch <- c("rna_ADAM17", "rna_APH1A", "rna_CREBBP", "rna_CTBP1", "rna_CTBP2", "rna_DTX3L", "rna_DTX4", "rna_EP300", "rna_HDAC1", "rna_HDAC2")    #các gene biểu hiện cho tín hiệu Notch: Tín hiệu Notch quá mức có thể dẫn đến ung thư bạch cầu cấp dòng T (T-ALL), ung thư vú, tuyến tụy.

DotPlot(object = CRER_Tcell, features = Notch, scale = T, split.by = "source3", cols = c("blue", "blue")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("CR vs ER NOTCH signal genes")

VlnPlot(object = CRER_Tcell, features = c("rna_TRAC", "rna_TRBC1", "rna_TRDC", "rna_TRGC1", "rna_TRAC"), pt.size = 0, ncol = 1)
FeaturePlot(CRER_Tcell, features = TCR_SET, pt.size = 1.0, max.cutoff = 2)


############### USING ESCAPE TO PERFORM GENE SET ENRICHMENT ANALYSES ON SINGLE-CELL RNA-SEQ-DATA

library(escape)                            
library(SingleCellExperiment)              
library(scran)                             
library(Seurat)
library(SeuratObject)
library(RColorBrewer)                      
library(ggplot2)
     

# saveRDS(object = CRER_Tcell, file = "./RData/CR_ER_T_SUB.rds")
CRER_Tcell <- readRDS(file = "./RData/06.CRER_Tcell_reclusters_annotated.rds")
View(CRER_Tcell@meta.data)
SUBset <- CRER_Tcell
SUBset <- subset(CRER_Tcell, ident = c("DNT cells"))
DefaultAssay(SUBset) <- "RNA"
SUBset = JoinLayers(SUBset)

DimPlot(SUBset, label = T, label.size = 4)

msigdbr_collections() # Kiểm tra các bộ gene sets hợp lệ

sce.DNT <-as.SingleCellExperiment(SUBset, assay = "RNA")
#GS.db <- getGeneSets(library = "C2", subcategory = "CP:KEGG_MEDICUS")    
GS.db <- getGeneSets(library = "C2",subcategory = "CP:KEGG_LEGACY")


enrichment.scores <- escape.matrix(SUBset, gene.sets = GS.db, groups = 100, min.size = 5)  

table(SUBset$source3) 

################# Tabls 1 total Tcells pop
idents <- Idents(CRER_Tcell)
cell_counts <- table(idents)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("ident", "cell_count")

orig_idents <-CRER_Tcell@meta.data$orig.ident
ident_data <- data.frame(orig_ident = orig_idents, ident = idents)  
cell_counts <- table(ident_data$ident, ident_data$orig_ident)       
cell_counts_df <- as.data.frame(cell_counts)                        
colnames(cell_counts_df) <- c("Cell type", "Annotation", "Cells")   
                                                                                           
write.csv(cell_counts_df, file = "./RData/Table/Main)_individual_P_counts.csv")
view(orig_idents)                                                  

ggplot(data = as.data.frame(enrichment.scores), 
       mapping = aes(enrichment.scores[,1], enrichment.scores[,2])) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.title = element_blank())

GS.db
SUBset <- runEscape(SUBset, method = "ssGSEA", gene.sets = GS.db, groups = 1000, min.size = 5, new.assay.name = "escape.ssGSEA")

sce.DNT <- runEscape(sce.DNT, method = "UCell", gene.sets = GS.db, groups = 1000, min.size = 5, new.assay.name = "escape.UCell")


colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)       
names(GS.db)

list(GS.db[["KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"]])


Target_GS <- "KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"

print(SUBset@assays)
FeaturePlot(SUBset, Target_GS,pt.size = 1.5) + 
  scale_color_gradientn(colors = colorblind_vector) + ggtitle(label = Target_GS)  

########################### Normalization
colorblind_vector <- hcl.colors(n=7, palette = "heat", fixup = TRUE)
colorblind_vector <- colorRampPalette(c("blue", "white", "red"))(100)

FeaturePlot(SUBset, features = Target_GS) + scale_color_gradientn(colors = colorblind_vector) + ggtitle(label = "HALLMARK-INTERFERON-GAMMA-RESPONE")

SUBset <- performNormalization(SUBset,                                    
                               assay = "escape.ssGSEA",
                               gene.sets = GS.db,
                               scale.factor = SUBset$nFeature_RNA)        
                                                                          
################ Plot heat map
P1 <- heatmapEnrichment(SUBset, 
                        group.by = "source3",
                        #gene.set.use = Target_GS,
                        cluster.rows = TRUE,       
                        cluster.columns = TRUE,
                        assay = "escape.ssGSEA") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                                                         axis.text.y = element_text(angle = 180,hjust = 0, size = 8),
                                                         legend.position = "top",
                                                         legend.text = element_text(angle = 90))
heatmapEnrichment(sce.DNT,
                  group.by = "source3",
                  assay = "escape.UCell",
                  scale = TRUE,
                  cluster.rows = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 180,hjust = 0),legend.position = "top",
                                               legend.text = element_text(angle = 90))

heatmapEnrichment(SUBset,
                  assay = "escape.ssGSEA",
                  palette = "Spectral",
                  group.by = "source3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



P2 <- heatmapEnrichment(SUBset,
                        group.by = "source3",
                        cluster.rows = TRUE,
                        cluster.columns = TRUE,
                        gene.set.use = rownames(SUBset@assays$escape.ssGSEA@data)[1:100],   
                        assay = "escape.ssGSEA") + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 1), 
                                                         axis.text.y= element_text(size = 8, face = "bold"),
                                                         legend.position = "top",
                                                         legend.key.size = unit(0.8, "cm")) 

P1
P2

################geyserEnrichment: Generate a ridge plot to examine enrichment distributions
list(GS.db$`KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION`)
Target_GS <- "KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"

geyserEnrichment(SUBset,
                 assay = "escape.ssGSEA",
                 gene.set = Target_GS,
                 order.by = "mean")

geyserEnrichment(SUBset, 
                 assay = "escape.ssGSEA",
                 gene.set = Target_GS,
                 facet.by = "groups")

geyserEnrichment(SUBset, 
                 assay = "escape.ssGSEA",
                 gene.set = Target_GS, order.by = "mean",
                 color.by = Target_GS, group.by = "source3", facet.by = "integrated_snn_res.0.48") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

list(GS.db$`KEGG-OXIDATIVE-PHOSPHORYLATION`)
Target_GS <- "KEGG-OXIDATIVE-PHOSPHORYLATION"

H1 <- geyserEnrichment(SUBset,
                       assay = "escape.ssGSEA",
                       gene.set = Target_GS, order.by = "mean",
                       color.by = Target_GS, group.by = "source3") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 15, face = "bold"),
        legend.position = "top", legend.key.size = unit(0.8, "cm"))
H1

# ggsave(filename = "KEGG_geneset_5.png", plot = H1, path = "./RData/@test", width = 6, height = 8, dpi = 600)

############### 

list(GS.db$'HALLMARK-ALLOGRAFT-REJECTION')
a <- "KEGG-OXIDATIVE-PHOSPHORYLATION"

geyserEnrichment(SUBset, 
                 assay = "escape.ssGSEA",
                 gene.set = a,
                 color.by = a, group.by = "source3")

ridgeEnrichment(sce.DNT, 
                assay = "escape.UCell",
                gene.set = a,
                add.rug = TRUE,                     
                scale = TRUE, group.by = "source3")  

densityEnrichment(sce.DNT,
                  gene.set.use = a,                                                 
                  gene.sets = GS.db, group.by = "source3") + ggtitle(label = a)     
####################
SUBset <- performPCA(SUBset,                        
                     assay = "escape.ssGSEA",
                     n.dim = 1:10)                  

pcaEnrichment(SUBset,                              
              dimRed = "escape.PCA",                
              x.axis = "PC1",                                          
              y.axis = "PC2")

pcaEnrichment(SUBset, facet.by = "source3",
              dimRed = "escape.PCA",
              x.axis = "PC1",
              y.axis = "PC2",
              add.percent.contribution = TRUE,                  
              display.factors = TRUE,                            
              number.of.factors = 15, style = "hex") + theme()   

SUBset <- performNormalization(SUBset,
                               assay = "escape.ssGSEA",
                               gene.sets = GS.db,
                               scale.factor = SUBset$nFeature_RNA)
SUBset2<- SUBset
SUBset2 <- SetIdent(SUBset2, value = SUBset2@meta.data$source3)  

DNT.all.markers <- FindAllMarkers(SUBset2,
                                  assay = "escape.ssGSEA_normalized",
                                  min.pct = 0,                   
                                  logfc.threshold = 0)          
head(DNT.all.markers)

write.csv(DNT.all.markers, file = "./RData/Table/Hallmarkers.csv")

######################################
######################################

library(ggpubr)
VlnPlot(CRER_Tcell, features = ("rna_HAVCR2"), split.by = "source3", idents = c("DNT cells", "MAIT cells"))
cell_counts <- table(CRER_Tcell@meta.data$source3)
view(cell_counts)

#In số lượng tế bào của từng loại
control_cells <- cell_counts["AML_CR"]     
treatment_cells <- cell_counts["AML_ER"]  

print(paste("Control cells:", control_cells))
print(paste("Treatment cells:", treatment_cells))


CRER_Tcell@meta.data$Diagnosis <- rep(c("AML_CR", "AML_ER"), time = c(9168, 10213-9168))   


#=====================
vp_case1 <- function(gene_signature, file_name, test_sign, y_max, nident) {
  plot_case1 <- function(signature, y_max = NULL) {
    p <- VlnPlot(CRER_Tcell, features = signature,
                 pt.size = 0, 
                 group.by = "Diagnosis", 
                 y.max = y_max, idents = c(nident)) + 
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") + 
      stat_compare_means(comparisons = test_sign, label = "p.format") 
    return(p)   
  }
  
  plot_list <- list() 
  y_max_list <- list()  
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]], na.rm = TRUE) 
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1)) 
  }
  
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 6, height = 8)
}

gene_sig <- "rna_NFKB1"
comparisons <- list(c("AML_CR", "AML_ER"))
vp_case1(gene_signature = gene_sig, file_name = "./RData/@test/test",
         test_sign = comparisons, y_max = 7,nident = "DNT cells")

VlnPlot(CRER_Tcell, features = c("rna_GZMB"), split.by = "source3", pt.size = 0)

#############Additional Analysis
DimPlot(SUBset)

DNT_subset <- subset(CRER_Tcell, ident="DNT cells")
DNT_subset <- RunPCA(DNT_subset, npcs = 30, verbose = FALSE)
DNT_subset =  FindNeighbors(DNT_subset, dims = 1:30)
DNT_subset =  FindClusters(DNT_subset, resolution = 0.5, algorithm = 1)
DNT_subset =  RunTSNE(DNT_subset, dims = 1:30)
DNT_subset =  RunUMAP(DNT_subset, dims = 1:10)

DimPlot(DNT_subset, split.by = "integrated_snn_res.0.48")
FeaturePlot(DNT_subset, features = c("rna_TP53"))
VlnPlot(DNT_subset, features = c("rna_TP53"),split.by = "source3")
FeaturePlot(DNT_subset, features = c("rna_GZMB"), max.cutoff = 2)



# sessionInfo()
# 
# 
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Korean_Korea.utf8  LC_CTYPE=Korean_Korea.utf8    LC_MONETARY=Korean_Korea.utf8
# [4] LC_NUMERIC=C                  LC_TIME=Korean_Korea.utf8    
# 
# time zone: Asia/Seoul
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3          scran_1.32.0                scuttle_1.14.0              SingleCellExperiment_1.26.0
# [5] SummarizedExperiment_1.34.0 Biobase_2.64.0              GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
# [9] IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
# [13] matrixStats_1.3.0           escape_2.0.0                EnhancedVolcano_1.22.0      ggrepel_0.9.5              
# [17] SeuratData_0.2.2.9001       HGNChelper_0.8.14           openxlsx_4.2.5.2            presto_1.0.0               
# [21] Rcpp_1.0.12                 data.table_1.15.4           fgsea_1.30.0                msigdbr_7.5.1              
# [25] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
# [29] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [33] tidyverse_2.0.0             patchwork_1.2.0             ggplot2_3.5.1               SeuratDisk_0.0.0.9021      
# [37] Seurat_5.1.0                SeuratObject_5.0.2          sp_2.1-4                   
# 
# loaded via a namespace (and not attached):
#   [1] GSVA_1.52.3               spatstat.sparse_3.1-0     httr_1.4.7                tools_4.4.1              
# [5] sctransform_0.4.1         utf8_1.2.4                R6_2.5.1                  HDF5Array_1.32.0         
# [9] lazyeval_0.2.2            uwot_0.2.2                ggdist_3.3.2              rhdf5filters_1.16.0      
# [13] withr_3.0.0               gridExtra_2.3             progressr_0.14.0          cli_3.6.3                
# [17] spatstat.explore_3.2-7    fastDummies_1.7.3         spatstat.data_3.1-2       ggridges_0.5.6           
# [21] pbapply_1.7-2             R.utils_2.12.3            parallelly_1.37.1         limma_3.60.3             
# [25] rstudioapi_0.16.0         RSQLite_2.3.7             generics_0.1.3            ica_1.0-3                
# [29] spatstat.random_3.2-3     distributional_0.4.0      zip_2.3.1                 Matrix_1.7-0             
# [33] fansi_1.0.6               abind_1.4-5               R.methodsS3_1.8.2         lifecycle_1.0.4          
# [37] edgeR_4.2.0               rhdf5_2.48.0              SparseArray_1.4.8         Rtsne_0.17               
# [41] grid_4.4.1                blob_1.2.4                promises_1.3.0            dqrng_0.4.1              
# [45] crayon_1.5.3              miniUI_0.1.1.1            lattice_0.22-6            beachmat_2.20.0          
# [49] cowplot_1.1.3             annotate_1.82.0           KEGGREST_1.44.1           magick_2.8.3             
# [53] metapod_1.12.0            pillar_1.9.0              rjson_0.2.21              future.apply_1.11.2      
# [57] codetools_0.2-20          fastmatch_1.1-4           leiden_0.4.3.1            glue_1.7.0               
# [61] vctrs_0.6.5               png_0.1-8                 spam_2.10-0               gtable_0.3.5             
# [65] cachem_1.1.0              S4Arrays_1.4.1            mime_0.12                 survival_3.7-0           
# [69] statmod_1.5.0             bluster_1.14.0            fitdistrplus_1.2-1        ROCR_1.0-11              
# [73] nlme_3.1-165              bit64_4.0.5               RcppAnnoy_0.0.22          irlba_2.3.5.1            
# [77] KernSmooth_2.23-24        splitstackshape_1.4.8     colorspace_2.1-0          DBI_1.2.3                
# [81] UCell_2.8.0               tidyselect_1.2.1          bit_4.0.5                 compiler_4.4.1           
# [85] AUCell_1.26.0             graph_1.82.0              BiocNeighbors_1.22.0      hdf5r_1.3.11             
# [89] DelayedArray_0.30.1       plotly_4.10.4             scales_1.3.0              lmtest_0.9-40            
# [93] rappdirs_0.3.3            SpatialExperiment_1.14.0  digest_0.6.36             goftest_1.2-3            
# [97] spatstat.utils_3.0-5      XVector_0.44.0            htmltools_0.5.8.1         pkgconfig_2.0.3          
# [101] sparseMatrixStats_1.16.0  fastmap_1.2.0             rlang_1.1.4               htmlwidgets_1.6.4        
# [105] UCSC.utils_1.0.0          shiny_1.8.1.1             DelayedMatrixStats_1.26.0 zoo_1.8-12               
# [109] jsonlite_1.8.8            BiocParallel_1.38.0       R.oo_1.26.0               BiocSingular_1.20.0      
# [113] magrittr_2.0.3            GenomeInfoDbData_1.2.12   dotCall64_1.1-1           Rhdf5lib_1.26.0          
# [117] munsell_0.5.1             babelgene_22.9            reticulate_1.38.0         stringi_1.8.4            
# [121] zlibbioc_1.50.0           MASS_7.3-61               plyr_1.8.9                parallel_4.4.1           
# [125] listenv_0.9.1             deldir_2.0-4              Biostrings_2.72.1         splines_4.4.1            
# [129] tensor_1.5                hms_1.1.3                 locfit_1.5-9.10           igraph_2.0.3             
# [133] spatstat.geom_3.2-9       RcppHNSW_0.6.0            reshape2_1.4.4            ScaledMatrix_1.12.0      
# [137] XML_3.99-0.17             tzdb_0.4.0                httpuv_1.6.15             RANN_2.6.1               
# [141] polyclip_1.10-6           future_1.33.2             scattermore_1.2           rsvd_1.0.5               
# [145] xtable_1.8-4              RSpectra_0.16-1           later_1.3.2               viridisLite_0.4.2        
# [149] ggpointdensity_0.1.0      memoise_2.0.1             AnnotationDbi_1.66.0      cluster_2.1.6            
# [153] timechange_0.3.0          globals_0.16.3            GSEABase_1.66.0   


########################################################################################################
########################################################################################################