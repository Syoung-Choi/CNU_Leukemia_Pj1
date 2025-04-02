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

table(Idents(CRER_Tcell))             #Đếm số tế bào trong từng nhóm
prop.table(table(Idents(CRER_Tcell))) #Tính tỷ lệ của mỗi nhóm so với tổng số tế bào

Idents(AML_CRER) <- "source3"
table(Idents(AML_CRER))
#Idents(AML_CRER)[Idents(AML_CRER) == "CD4"] <- "source3": Chỉ đổi nhóm tế bào CD4 thành "source3", không ảnh hưởng nhóm khác.

Assay_title <- "ER vs CR"       #Tạo một biến Assay_title có giá trị "ER vs CR DNT".Không ảnh hưởng đến dữ liệu trong AML_CRER
DefaultAssay(AML_CRER) <- "RNA" #Thiết lập "RNA" làm assay mặc định cho đối tượng AML_CRER -> "RNA" (biểu hiện gene gốc)/ "SCT" (dữ liệu đã chuẩn hóa với SCTransform)/ "ADT" (dữ liệu protein từ CITE-seq -> Nếu bạn đang làm việc với dữ liệu RNA và muốn phân tích biểu hiện gene, thì cần đặt "RNA" làm assay mặc định. Việc đặt assay mặc định giúp tiết kiệm thời gian, vì bạn không cần phải chỉ định cụ thể assay mỗi khi gọi một hàm. 
AML_CRER = JoinLayers(AML_CRER) #gộp các layers dữ liệu trong một đối tượng Seurat.Hữu ích khi có dữ liệu từ nhiều lớp khác nhau (RNA, ADT, ATAC-seq).Đảm bảo rằng tất cả layers có cùng kích thước và thứ tự ô để thực hiện phân tích tiếp theo.

Idents(AML_CRER) <- "celltype.condition" #gọi AML_CRER là celltype.condition.
myData.genes<- wilcoxauc(AML_CRER, "source3") #Area Under the Curve (AUC): phân loại để đánh giá hiệu quả của mô hình hoặc sự khác biệt giữa các nhóm
                                              #Kết quả của hàm wilcoxauc() được lưu trong đối tượng myData.genes, có thể chứa các giá trị AUC cho các gene hoặc nhóm tế bào trong AML_CRER liên quan đến điều kiện "source3".

#tạo danh sách tập hợp gene phục vụ cho phân tích GSEA, enrichment 
m_df<- msigdbr(species = "Homo sapiens", category = "H") #trích xuất tập hợp gene từ MSigDB (Molecular Signatures Database).
                                                         #species = "Homo sapiens": Chỉ định rằng chúng ta muốn lấy dữ liệu cho loài người.
                                                         #"H"allmark gene sets: một tập hợp gồm khoảng 50 bộ gene đại diện cho các tín hiệu sinh học đặc trưng của tế bào.giúp giảm nhiễu và tập trung vào các con đường tín hiệu sinh học có ý nghĩa sinh học quan trọng.
                                                         #m_df sẽ là một data frame chứa danh sách các gene thuộc nhóm Hallmark trong MSigDB dành cho loài người.


#Chuyển đổi m_df thành một danh sách (list) chứa các tập hợp gen cho phân tích FGSEA (Fast Gene Set Enrichment Analysis).
fgsea_set<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) #.$gene_symbol: Cột chứa tên gen.
                                                             #.$gs_name: Cột chứa tên tập hợp gen.
                                                             # split(): nhóm các giá trị gene_symbol theo từng nhóm gs_name -> tạo ra một danh sách trong đó mỗi phần tử là một tập hợp các gen thuộc một nhóm.

head(m_df)                        # mặc định là hiện thị 6 dòng đầu tiên trong m_df
dplyr::count(myData.genes, group) # Đếm số lượng phần tử trong mỗi nhóm (group) trong myData.genes

myData.genes %>%                        
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n=10)                          #Lấy 10 gen có logFC và sau đó auc cao nhất trong nhóm "AML_ER" ( sắp xếp giảm dần logFC rồi đến auc).

AML_ER.genes <- myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)         #Lọc các gen thuộc nhóm "AML_ER", sắp xếp theo auc giảm dần và chỉ giữ lại hai cột feature và auc.

ranks <- deframe(AML_ER.genes)        #deframe() giúp chuyển data frame (2 cột) thành named vector -> chuyển từ bảng dọc thành bảng ngang.Rất hữu ích cho GSEA/FGSEA
head(ranks)
##############################

fgseaRes <- fgsea(fgsea_set, stats = ranks, nperm = 1000) #phân tích mức độ phong phú (enrichment) của các gene sets trong một bộ dữ liệu được xếp hạng.
                                                          #stats: truyền vào hàm fgsea là vector xếp hạng của các gen
                                                          #nperm(permutation)= 1000: thực hiện 1000 lần hoán vị để kiểm tra tính ngẫu nhiên của kết quả -> kết quả từ phân tích thực tế (trong dữ liệu của bạn) sẽ được so sánh với phân phối ngẫu nhiên này để tính toán giá trị p-value cho mức độ phong phú của mỗi gene set.nMoreExtreme nhỏ cho thấy rằng sự phong phú của gene set trong dữ liệu thực tế có thể không phải ngẫu nhiên, đáng tin cậy.

#NES: cao xuống thấpthấp
fgseaResTidy <- fgseaRes %>%    
  as_tibble() %>%               #Chuyển fgseaRes thành một tibble để dễ dàng thao tác với dữ liệu.
  arrange(desc(NES))            #Sắp xếp dữ liệu theo cột NES(Normalized Enrichment Score_mức độ phong phú của gene set) giảm dần (tức là gene set có NES cao nhất sẽ đứng đầu).

#NES: thấp lên cao
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%     #-leadingEdge (Thông tin về các gen trong vùng gen quan trọng nhất của gene set), -ES ( Enrichment Score, là điểm số phong phú của gene set), -nMoreExtreme(Số lượng hoán vị có kết quả cực đoan hơn so với kết quả quan sát thực tế) có nghĩa là loại bỏ các cột này khỏi fgseaResTidy. 
  arrange(padj) %>%                                       # padj là adjusted p-value (giá trị p đã điều chỉnh),(multiple testing correction), mặc định là tăng dần.
  head()                                                  # head() sẽ chỉ giữ lại 6 gene sets đầu tiên sau khi đã sắp xếp theo padj.

ggplot(fgseaResTidy %>% 
         dplyr::filter(padj < 0.05) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +   #⚠ Lưu ý: Nếu bạn muốn chắc chắn chọn 50 pathways có NES cao nhất/thấp nhất, bạn nên sắp xếp dữ liệu trước bằng arrange(desc(NES)) hoặc arrange(NES)_ mặc định từ thấp lên cao.
  geom_col(aes(fill = NES<0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = Assay_title,) +
  theme_minimal()

#dplyr::filter (Dùng để lọc hàng (rows) của một data.frame hoặc tibble dựa trên điều kiện cụ thể, xử lý dữ liệu dạng bảng): khác với filter trong stats(khi làm việc với chuỗi số liệu (numeric vector, time series)., thực hiện phép làm mượt dữ liệu bằng bộ lọc tuyến tính)
#dplyr::filter(padj < 0.05): Chỉ định sử dụng hàm filter từ dplyr để lọc dữ liệu mà padj < 0.05.
#head(n = 50): Chọn 50 pathway có NES cao nhất hoặc thấp nhất sau khi lọc. -> chỉ có 35 hallmarks (pathway) đủ điều kiện p_ajust<0.05 
#Xác định các trục aes(): aes(reorder(pathway, NES), NES) → 
                                                         #pathway (trục x): Được sắp xếp theo giá trị NES (reorder(pathway, NES)).
                                                         #NES (trục y): Giá trị làm giàu của từng pathway.
#geom_col() → Vẽ biểu đồ cột (bar plot) với chiều cao của mỗi cột là NES và phân chia cột làm 2 màu.
                                 #Nếu NES > 0 → Một màu (thường là xanh).
                                 #Nếu NES < 0 → Một màu khác (thường là đỏ).
#Lật trục (coord_flip()): Chuyển biểu đồ cột từ dạng dọc sang ngang để dễ đọc hơn.Giờ đây, pathway sẽ nằm trên trục y thay vì trục x.
#theme_minimal() → Xóa bỏ các đường kẻ nền, tạo giao diện đơn giản và dễ nhìn hơn.

Pathway <-"HALLMARK_OXIDATIVE_PHOSPHORYLATION" 

#Đường cong biểu thị enrichment score (ES), điểm cao nhất là điểm làm giàu tối đa.
plotEnrichment(fgsea_set[[Pathway]],                     #fgsea_set[[Pathway]] trích xuất danh sách các gene từ Biến Pathway là tên của một bộ gen 
               ranks) + labs(title=Pathway)              #ranks: Là một vector có thứ hạng của các gene (thường được xếp hạng dựa trên thống kê như log fold-change hoặc p-value).Được sử dụng để đánh giá mức độ làm giàu của tập hợp gene trong đường dẫn Pathway.
                                                        



class(top5)
class(fgsea_set)

top5[1,]
fgsea_set[Pathway] 
fgsea_set[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]       #fgsea_sets là một list trong R. Mỗi phần tử trong list này đại diện cho một pathway, chứa các gene liên quan.
                                                        #fgsea_sets[[Pathway]] trả về nội dung thực sự của phần tử đó (không phải là list mà là vector hoặc object bên trong).
###########################
AML_SUB= readRDS('./RData/T_NK_subset_V.1.rds')
DNT_ALL <- c("rna_CD160", "rna_CD44", "rna_FCER1G", "rna_ANXA2", "rna_IFNG", "rna_S100A11", "rna_XCL1", "rna_PIM1", "rna_IKZF2", "rna_ATF4")
DotPlot(AML_SUB,features = DNT_ALL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #tạo bảng biểu hiện mức độ tập hợp các gene trên mỗi clusters. 
                                                           #axis.text.x: Đây là một thuộc tính của theme(), được dùng để định dạng nhãn trên trục x.
                                                           #element_text(): Hàm dùng để tùy chỉnh văn bản (text), bao gồm màu sắc, kích thước, độ nghiêng, căn chỉnh, v.v.
                                                                            #angle = 90: Xoay văn bản của trục x một góc 90 độ (tức là hiển thị theo chiều dọc).
                                                                            #hjust = 1: (horizontal justification) điều chỉnh căn chỉnh theo chiều ngang.

AML_DNT= readRDS('./RData/T_NK_subset_V.1.rds')

FeaturePlot(AML_DNT, feature = c("rna_HAVCR2", "rna_TIGIT", "rna_CTLA4", "rna_PDCD1", "rna_LAG3"), ncol = 5, split.by = "orig.ident") #ncol = 5: Sắp xếp các biểu đồ theo 5 cột trong bố cục hiển thị.
FeaturePlot(AML_DNT, feature = c("rna_HAVCR2"), ncol = 4, split.by = "orig.ident")
FeaturePlot(AML_DNT, feature = c("rna_IKZF2"), ncol = 4, split.by = "orig.ident")

FeaturePlot(AML_DNT, features = c("rna_LGALS9", "rna_PVR", "rna_CD80", "rna_CD274", "rna_CEACAM1"), ncol = 4, split.by ="orig.ident")
FeaturePlot(AML_DNT, features = c("rna_LGALS9"), ncol = 4, split.by = "orig.ident")

DNT_surface <- c("rna_CD160", "rna_CD44", "rna_FCER1G")
DNT_secreted <- c("rna_ANXA2", "rna_IFNG", "rna_S100A11", "rna_XCL1")
DNT_IF <- c("rna_SUB1", "rna_PIM1", "rna_IKZF2", "rna_ATF4")

DoHeatmap(AML_DNT, features = DNT_ALL, group.by = "source3")
DNT_ALL <- gsub("rna_", "", DNT_ALL) #Các gene trong rownames(AML_DNT) không có prefix "rna_", trong khi danh sách DNT_ALL của bạn lại có "rna_" ở đầu tên gene.Vì vậy, DoHeatmap() không tìm thấy các gene do tên không khớp.
                                     #gsub("rna_", "", DNT_ALL) sẽ loại bỏ "rna_" khỏi tên gene, giúp nó khớp với rownames(AML_DNT --> 


#install.packages("openxlsx", dependencies = TRUE) #dùng để đọc, ghi và xử lý file Excel (.xlsx).dependencies = TRUE giúp tự động cài đặt các package phụ thuộc nếu cần.
#install.packages('HGNChelper')                    #Giúp đảm bảo các gene trong dữ liệu RNA-seq của bạn đúng với danh pháp chuẩn.HGNC (HUGO Gene Nomenclature Committee).
#remotes::install_github("satijalab/seurat-data")  #dùng để cài đặt trực tiếp từ GitHub thay vì CRAN.

# CRAN (Comprehensive R Archive Network): Là kho lưu trữ chính thức của R, nơi chứa các package ổn định, được kiểm duyệt kỹ trước khi xuất bản.
# GitHub Là nền tảng dành cho các nhà phát triển chia sẻ code. cập nhật nhanh, nhiều tính năng mới nhưng có thể chưa ổn định.

library(openxlsx)
library(HGNChelper)
library(SeuratData)
library(EnhancedVolcano)

R.version.string # kiểm tra phiên bản R hiện tại

AML_DNT <- subset(CRER_Tcell, idents = "DNT cells") #gọi DNT cells từ CRER_T cells là AML_DNT.
DimPlot(AML_DNT, reduction = "umap", label = F, split.by = "orig.ident")
#===================================
AML_DNT[["RNA"]]
AML_DNT[["integrated"]] 

class(AML_DNT[["RNA"]])
class(AML_DNT[["integrated"]])

names(AML_DNT)

#Dữ liệu của bạn chứa hai assay: "RNA" và "integrated". Mỗi assay có một tập hợp các đặc trưng (features) khác nhau.
                                 #Assay "RNA" chứa 27,947 đặc trưng (genes) cho 20,277 tế bào.
                                 #Assay "integrated" chứa 2,578 đặc trưng cho 20,277 tế bào, với các đặc trưng biến đổi (variable features) được chọn trong quá trình tích hợp dữ liệu.


DefaultAssay(AML_DNT) <- "RNA"       # Trong một object Seurat có thể có nhiều assay, chẳng hạn như RNA (expression data), ADT (protein expression), SCT (scaled expression data), v.v. Bạn cần xác định rõ assay nào sẽ được sử dụng mặc định cho các phân tích tiếp theo (chẳng hạn như tính toán PCA, clustering, v.v.).
AML_DNT_1 = JoinLayers(AML_DNT)      # để join phần intergrated với RNA --> lỗi do JoinLayers yêu cầu một lớp (class) cụ thể của đối tượng để thực thi. 

cs <- colSums(GetAssayData(AML_DNT_1,assay="RNA",layer="counts")>0) #lấy dữ liệu từ assay "RNA" trong đối tượng Seurat AML_DNT_1, đặc biệt là dữ liệu trong lớp counts. Lớp counts chứa thông tin về số lượng đọc cho mỗi gen trong mỗi tế bào.Nó sẽ tạo ra một ma trận boolean, trong đó mỗi phần tử có giá trị TRUE nếu giá trị ban đầu lớn hơn 0 và FALSE nếu giá trị đó bằng 0
ce <- cs[which(cs>5)]                                               #colSums: đếm số lượng gen (hoặc biểu hiện gen) có giá trị lớn hơn 0 trong mỗi tế bào (cột). Kết quả là một vector, trong đó mỗi phần tử đại diện cho số lượng gen biểu hiện trong một tế bào cụ thể.
sf <- cs[names(ce)]                                                 #lấy được giá trị từ cs chỉ cho các tế bào có ít nhất 6 gen có giá trị lớn hơn 0.

Idents(AML_DNT)<-"celltype.condition"
myData.genes <-wilcoxauc(AML_DNT_1,"source3") #Kết quả của phép kiểm tra Wilcoxon với AUC sẽ được lưu vào biến myData.genes. Kết quả này có thể chứa các gen khác biệt giữa các điều kiện, kèm theo các giá trị AUC tương ứng.

#m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG") 
msigdbr_collections()
m_df<- msigdbr(species = "Homo sapiens", collection = "C2",subcollection = "CP:KEGG_MEDICUS") #KEGG Legacy Pathways chứa các đường dẫn cũ, trong khi KEGG Medicus Pathways chứa các đường dẫn cập nhật với các thông tin y học hiện đại.

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name) # tạo ra file mới lấy data từ m_df mà trong đó gồm ít nhất 2 cột là gene_symbol và gs_name.Hàm split() dùng để chia nhỏ dữ liệu trong cột gene_symbol dựa trên giá trị của cột gs_name (tên bộ gen).
                                                               #x = .$gene_symbol: Đây là phần dữ liệu mà bạn muốn chia (các biểu tượng gen).
                                                               #f(factor) = .$gs_name: Đây là yếu tố phân nhóm (tên bộ gen). Các biểu tượng gen sẽ được nhóm theo giá trị trong cột gs_name.

head(m_df) 
dplyr::count(myData.genes,group)

myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%    #Chỉ giữ lại những dòng có giá trị group bằng "AML_ER".
  arrange(desc(logFC),desc(auc)) %>%      #logFC (log fold-change, thể hiện mức độ thay đổi biểu hiện gene). auc (area under curve, có thể dùng để đánh giá độ quan trọng của gene)
  head(n=10)

AML_ER.genes<- myData.genes %>%
  dplyr::filter(group == "AML_ER") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)            #Chỉ giữ lại hai cột: feature (tên gene) và auc (giá trị AUC).

rank <- deframe(AML_ER.genes)            #chuyển đổi dữ liệu từ AML_ER.genes (data frame/tibble) sang một named vector bằng hàm deframe() của gói tibble.
head(rank)                               #Giúp dễ dàng truy xuất giá trị dựa trên tên.

################
#fgsea() là hàm trong thư viện fgsea, dùng để thực hiện Gene Set Enrichment Analysis (GSEA), một phương pháp thống kê nhằm tìm các tập hợp gene (gene sets) có sự thay đổi đáng kể trong một nghiên cứu phân tích dữ liệu gene expression.
fgseaRes <-fgsea(fgsea_sets, stats = ranks, nperm = 1000)  #stats: Đây là vector chứa các giá trị thống kê được tính từ dữ liệu của bạn, chẳng hạn như logFC (log fold change) hoặc AUC (Area Under Curve) cho mỗi gene.
                                                           #nperm: hoán vị ngẫu nhiên (permutation). thực hiện 1000 lần hoán vị ngẫu nhiên để đánh giá tính đáng tin cậy của kết quả GSEA.
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))                                       #NES (Normalized Enrichment Score)

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, nMoreExtreme) %>%
  arrange(padj) %>%
  head()

conflicts_prefer(dplyr::filter) #ưu tiên sử dụng filter trong hàm dplyr thay vì hàm stat. 

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
#https://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html.

#BiocManager::install("escape")
#BiocManager::install("scran")
#BiocManager::install("pathview")
library(escape)                            #Hỗ trợ ssGSEA (single-sample GSEA).Phân tích hoạt động của signaling pathways trên từng tế bào.
library(SingleCellExperiment)              #Hỗ trợ tiền xử lý, lọc và chuẩn hóa dữ liệu.Kết hợp với các package khác như scran, Seurat.
library(scran)                             #Chuẩn hóa dữ liệu scRNA-seq (deconvolution normalization).Phát hiện các gene có sự biến thiên cao (highly variable genes).Nhóm tế bào bằng graph-based clustering.
library(Seurat)
library(SeuratObject)
library(RColorBrewer)                      #Cung cấp bảng màu đẹp và chuyên nghiệp để trực quan hóa dữ liệu trong R.
library(ggplot2)
base::library(package = "pathview")        #Pathview là một công cụ để vẽ sơ đồ KEGG pathway có tích hợp dữ liệu biểu hiện gene hoặc protein.

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
#GS.db <- getGeneSets(library = "C2", subcategory = "CP:KEGG_MEDICUS")    #getGeneSets(): Đây là hàm dùng để lấy một hoặc nhiều gene sets từ MSigDB. MSigDB có các bộ gene sets chia theo các loại khác nhau như "H", "C5", "C2", v.v.
                                                                          #library = "C2": Chỉ định rằng bạn muốn lấy gene sets từ thư viện "C2". Thư viện C2 chứa các gene sets liên quan đến các đường truyền (pathways), như KEGG, Reactome, v.v.

GS.db <- getGeneSets(library = "C2",subcategory = "CP:KEGG_LEGACY")

#SingleCellExperiment là định dạng chuẩn cho nhiều công cụ phân tích dữ liệu đơn tế bào trong R, chẳng hạn như các gói scran, scater, và Seurat. Bằng cách chuyển đổi dữ liệu thành SingleCellExperiment, bạn có thể tận dụng tất cả các chức năng và công cụ của các gói này một cách mượt mà.

enrichment.scores <- escape.matrix(SUBset, gene.sets = GS.db, groups = 100, min.size = 5)  
#escape.matrix() để tính toán điểm làm giàu (enrichment scores) của một tập hợp gene trong tập dữ liệu SUBset dựa trên cơ sở dữ liệu gene GS.db. 
#groups(n=100): The number of cells to separate the enrichment calculation.

table(SUBset$source3) #kiểm tra số lượng cell trong mỗi nhóm (CR và ER)

################# Tabls 1 total Tcells pop
#idents <- Idents(CRER_Tcell)
#cell_counts <- table(idents)
#cell_counts_df <- as.data.frame(cell_counts)
#colnames(cell_counts_df) <- c("ident", "cell_count")
#write.csv(cell_counts_df, file = "./RData/Table/Main)_Tcells_pop.csv")

orig_idents <-CRER_Tcell@meta.data$orig.ident
ident_data <- data.frame(orig_ident = orig_idents, ident = idents)  #Tạo DataFrame chứa thông tin về tế bào
cell_counts <- table(ident_data$ident, ident_data$orig_ident)       #Tính số lượng tế bào theo từng nhóm
cell_counts_df <- as.data.frame(cell_counts)                        #Chuyển đổi bảng đếm thành DataFrame
colnames(cell_counts_df) <- c("Cell type", "Annotation", "Cells")   #Đổi tên các cột thành:# "Cell type": Loại tế bào (ident).
                                                                                           # "Annotation": Nguồn gốc/điều kiện (orig_ident).
                                                                                           # "Cells": Số lượng tế bào (table() đã đếm).
write.csv(cell_counts_df, file = "./RData/Table/Main)_individual_P_counts.csv")
view(orig_idents)                                                  #ident = active.ident: nhóm phân loại tế bào dựa trên cluster (DNT, CD8,...)
                                                                   #orig.ident: Nguồn gốc của tế bào (bệnh nhân, điều kiện thí nghiệm)_AML_CR,AML_ER.



ggplot(data = as.data.frame(enrichment.scores), 
       mapping = aes(enrichment.scores[,1], enrichment.scores[,2])) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.title = element_blank())
#scatter plot (biểu đồ phân tán) dựa trên dữ liệu enrichment.scores.
#enrichment.scores[a,b]: a là số hàng muốn xét, b là số cột tương ứng --> để trống a là chỉ muốn xét các gene trong cột thứ nhất và thứ 2)
#as.data.frame(enrichment.scores): Nếu enrichment.scores là một ma trận (matrix), as.data.frame() giúp ggplot nhận dạng nó đúng cách.
#aes(enrichment.scores[,1], enrichment.scores[,2]): Chọn cột 1 làm trục X, cột 2 làm trục Y.
#geom_point(): Mỗi điểm (gene) tương ứng với một hàng trong enrichment.scores.

GS.db
SUBset <- runEscape(SUBset, method = "ssGSEA", gene.sets = GS.db, groups = 1000, min.size = 5, new.assay.name = "escape.ssGSEA")
#ssGSEA với runEscape(): để đánh giá mức độ hoạt động của các tập hợp gen trong từng tế bào -> trực quan hóa kết quả bằng biểu đồ và lưu vào Seurat/SinglecellExperimen 
#ssGSEA:calculates the enrichment score using a rank-normalized approach

sce.DNT <- runEscape(sce.DNT, method = "UCell", gene.sets = GS.db, groups = 1000, min.size = 5, new.assay.name = "escape.UCell")
#UCell:calculates a Mann-Whitney U statistic based on the gene rank list

colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)         #tạo một vector màu gồm 7 màu từ bảng màu "inferno" của hàm hcl.colors() trong R. Đây là một bảng màu có độ tương phản cao, giúp hỗ trợ người bị mù màu phân biệt dễ dàng.

names(GS.db) #kiểm tra tên hợp lệ.

list(GS.db[["KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"]])
# KEGG-CELL-CYCLE
# KEGG-DNA-REPLICATION
# KEGG-TGF-BETA-SIGNALING-PATHWAY
# KEGG-APOPTOSIS
# KEGG-PATHWAYS-IN-CANCER
# KEGG-T-CELL-RECEPTOR-SIGNALING-PATHWAY
# KEGG-JAK-STAT-SIGNALING-PATHWAY
# KEGG-OXIDATIVE-PHOSPHORYLATION
# KEGG-ACUTE-MYELOID-LEUKEMIA

Target_GS <- "KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"

print(SUBset@assays)
FeaturePlot(SUBset, Target_GS,pt.size = 1.5) + 
  scale_color_gradientn(colors = colorblind_vector) + ggtitle(label = Target_GS)  #Tạo một gradient màu với nhiều màu sắc.

########################### Normalization
colorblind_vector <- hcl.colors(n=7, palette = "heat", fixup = TRUE)
colorblind_vector <- colorRampPalette(c("blue", "white", "red"))(100)
#hcl.colors sử dụng một palette định nghĩa sẵn (trong trường hợp này là "heat"), trong khi colorRampPalette cho phép bạn tùy chỉnh dải màu bằng cách chỉ định các màu sắc đầu tiên và số lượng màu cần tạo ra.
#hcl.colors chỉ định số màu (7 màu trong ví dụ), còn colorRampPalette cho phép chỉ định số lượng màu (100 màu trong ví dụ).

FeaturePlot(SUBset, features = Target_GS) + scale_color_gradientn(colors = colorblind_vector) + ggtitle(label = "HALLMARK-INTERFERON-GAMMA-RESPONE")

SUBset <- performNormalization(SUBset,                                    #chuẩn hóa dữ liệu trong đối tượng SUBset.
                               assay = "escape.ssGSEA",
                               gene.sets = GS.db,
                               scale.factor = SUBset$nFeature_RNA)        #scale.factor được dùng để xác định yếu tố chuẩn hóa cho dữ liệu. 
                                                                          #Việc sử dụng nFeature_RNA có thể giúp chuẩn hóa dữ liệu biểu hiện gen dựa trên số lượng gen được phát hiện, giúp giảm tác động của sự khác biệt về số lượng gen đo được giữa các tế bào.
################ Plot heat map
P1 <- heatmapEnrichment(SUBset, 
                        group.by = "source3",
                        #gene.set.use = Target_GS,
                        cluster.rows = TRUE,       #cluster.rows = TRUE, cluster.columns = TRUE:cho phép bạn nhóm (cluster) các dòng (rows) và cột (columns) trong heatmap.Việc cluster giúp bạn phát hiện các nhóm tế bào hoặc nhóm bộ gen có sự tương đồng cao về mức độ biểu hiện hoặc các tính chất khác. 
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

#write.table(myData.genes,file = "./RData/MyData.genes.txt") -> export excel/text file

P2 <- heatmapEnrichment(SUBset,
                        group.by = "source3",
                        cluster.rows = TRUE,
                        cluster.columns = TRUE,
                        gene.set.use = rownames(SUBset@assays$escape.ssGSEA@data)[1:100],   #chọn một tập hợp con của các gen để vẽ(Chọn 100 gene đầu tiên từ danh sách gene của tập dữ liệu escape.ssGSEA).
                        assay = "escape.ssGSEA") + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 1), #face = "bold":in đậm 
                                                         axis.text.y= element_text(size = 8, face = "bold"),
                                                         legend.position = "top",
                                                         legend.key.size = unit(0.8, "cm")) #Điều chỉnh kích thước khung của legend thành 0.8 cm.

P1
P2

# ggsave(filename = "KEGG.png", plot = P1, path = "./RData/@test", width = 10, height = 40, dpi = 600)
# ggsave(filename = "KEGG_escape_ssGSEA1.png", plot = P2, path = "./RData/@test", width = 10, height = 12, dpi = 600)

################geyserEnrichment: Generate a ridge plot to examine enrichment distributions
# kiểm tra phân bố mức độ làm giàu của các tập hợp gene (gene sets) trong các nhóm mẫu khác nhau.

list(GS.db$`KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION`)
Target_GS <- "KEGG-ANTIGEN-PROCESSING-AND-PRESENTATION"

geyserEnrichment(SUBset,
                 assay = "escape.ssGSEA",
                 gene.set = Target_GS,
                 order.by = "mean")      #"mean" có nghĩa là sắp xếp theo giá trị trung bình của mức độ làm giàu trong các mẫu.

geyserEnrichment(SUBset, 
                 assay = "escape.ssGSEA",
                 gene.set = Target_GS,
                 facet.by = "groups")


colnames(SUBset@meta.data) #--> don't have "group" inside.

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
#ridgeEnrichment  : sự phân bố của Enrichment Scores theo từng nhóm cell_type, với mỗi nhóm là một đường density.
#densityEnrichment: vẽ biểu đồ mật độ (density plot) nhằm thể hiện sự phân bố enrichment scores của một tập hợp gene (gene.set.use = a) trên dữ liệu Single-Cell Expression (sce.DNT).

list(GS.db$'HALLMARK-ALLOGRAFT-REJECTION')
a <- "KEGG-OXIDATIVE-PHOSPHORYLATION"

geyserEnrichment(SUBset, 
                 assay = "escape.ssGSEA",
                 gene.set = a,
                 color.by = a, group.by = "source3")

ridgeEnrichment(sce.DNT, 
                assay = "escape.UCell",
                gene.set = a,
                add.rug = TRUE,                      # Hiển thị dấu gạch nhỏ dưới trục x
                scale = TRUE, group.by = "source3")  # Chuẩn hóa dữ liệu

densityEnrichment(sce.DNT,
                  gene.set.use = a,                                                 #Xác định tập hợp gene cần phân tích
                  gene.sets = GS.db, group.by = "source3") + ggtitle(label = a)     #cơ sở dữ liệu chứa các bộ gen mà bạn muốn sử dụng để tính toán sự phong phú.

####################
SUBset <- performPCA(SUBset,                        #PCA giúp giảm chiều dữ liệu và xác định các trục biến thiên quan trọng nhất.Điều này có thể hữu ích cho việc nhóm mẫu, phân tích cụm, hoặc tích hợp dữ liệu.
                     assay = "escape.ssGSEA",
                     n.dim = 1:10)                  #n.dim = 1:10: Giữ lại 10 thành phần chính (Principal Components - PCs) đầu tiên.

pcaEnrichment(SUBset,                               #pcaEnrichment(): Đây là một hàm có thể được dùng để đánh giá mức độ phong phú (enrichment) của các tập hợp gene trong không gian PCA.
              dimRed = "escape.PCA",                #dimRed = "escape.PCA": Sử dụng kết quả PCA vừa tính ở bước trước để vẽ biểu đồ
              x.axis = "PC1",                                               #PCA giúp giảm chiều xuống 2D hoặc 3D để dễ phân tích hơn (dimensionality-reduced data) 
              y.axis = "PC2")

pcaEnrichment(SUBset, facet.by = "source3",
              dimRed = "escape.PCA",
              x.axis = "PC1",
              y.axis = "PC2",
              add.percent.contribution = TRUE,                   #Hiển thị phần trăm đóng góp của PC1 và PC2 (VD: "PC1 (35%)", "PC2 (20%)").
              display.factors = TRUE,                            #Hiển thị các yếu tố có đóng góp lớn nhất trong PCA.
              number.of.factors = 15, style = "hex") + theme()   #Hiển thị 15 yếu tố có ảnh hưởng lớn nhất trong PCA.
                                                                 #Sử dụng biểu đồ dạng lưới hình lục giác (hexbin) thay vì scatter plot.
SUBset <- performNormalization(SUBset,
                               assay = "escape.ssGSEA",
                               gene.sets = GS.db,
                               scale.factor = SUBset$nFeature_RNA)
SUBset2<- SUBset
SUBset2 <- SetIdent(SUBset2, value = SUBset2@meta.data$source3)  # khi sử dụng FindALLMarkers cần có 2 nhóm tương ứng để so sánh. Nhưng ban đầu lại gán Subset là "DNT cells" nên không có nhóm tương ứng để đối chứng. 


DNT.all.markers <- FindAllMarkers(SUBset2,
                                  assay = "escape.ssGSEA_normalized",
                                  min.pct = 0,                   # min.pct              : xác định tỷ lệ tối thiểu của tế bào trong mỗi nhóm phải có sự biểu hiện của một gene để gene đó được xem là marker.  min.pct = 0 có nghĩa là không yêu cầu tỷ lệ tối thiểu (tất cả các gene đều có thể được kiểm tra, bất kể tần suất xuất hiện trong các nhóm).
                                  logfc.threshold = 0)           # log fold-change (FC) : những gene có sự thay đổi log FC lớn hơn ngưỡng này mới được coi là marker. logfc.threshold = 0 có nghĩa là không yêu cầu sự thay đổi biểu hiện (tất cả các gene có sự thay đổi nào đó đều sẽ được đưa vào phân tích).
head(DNT.all.markers)

write.csv(DNT.all.markers, file = "./RData/Table/Hallmarkers.csv")

######################################
######################################

#install.packages("ggpubr") -> vẽ biểu đồ khoa học dễ dàng hơn.
library(ggpubr)
VlnPlot(CRER_Tcell, features = ("rna_HAVCR2"), split.by = "source3", idents = c("DNT cells", "MAIT cells"))
cell_counts <- table(CRER_Tcell@meta.data$source3)
view(cell_counts)

#In số lượng tế bào của từng loại
control_cells <- cell_counts["AML_CR"]     #[]: Dùng để truy xuất hoặc chỉ định phần tử trong một đối tượng (như vector, table, data frame, matrix).
treatment_cells <- cell_counts["AML_ER"]   #(): Dùng để gọi hàm và thực thi một hành động trong R.

print(paste("Control cells:", control_cells))
print(paste("Treatment cells:", treatment_cells))

#thêm một cột Diagosis vào metadata, giúp xác định mỗi tế bào thuộc nhóm "AML_CR" hay "AML_ER", dựa trên số lượng đã đếm trước đó.
CRER_Tcell@meta.data$Diagnosis <- rep(c("AML_CR", "AML_ER"), time = c(9168, 10213-9168))    #rep(): Lặp lại giá trị theo số lần chỉ định.
#                                       control, sample           con.cells, max-con"ER"   "AML_CR" sẽ được lặp lại 9168 lần.
#                                                                                          "AML_ER" sẽ được lặp lại (10213 - 9168) = 1045 lần.  
#                                                                                          -> Kết quả là một vector chứa 9168 giá trị "AML_CR" và 1045 giá trị "AML_ER".
#           Dữ liệu cần đúng số lượng tế bào ở từng nhóm để 
#           tránh lỗi khi chạy các phân tích downstream như FindAllMarkers() hoặc vẽ biểu đồ với VlnPlot().


#=====================Vẽ biểu đồ violin (VlnPlot) để hiển thị sự khác biệt biểu hiện của các gene giữa các nhóm chẩn đoán (AML_CR, AML_ER).
vp_case1 <- function(gene_signature, file_name, test_sign, y_max, nident){
  plot_case1 <- function(signature, y_max =NULL){                               #Khai báo hàm con plot_case1  #y_max = NULL, tham số y.max trong VlnPlot() sẽ không bị giới hạn,
    VlnPlot(CRER_Tcell, features = signature,                                   #Chọn gene cụ thể để vẽ.
            pt.size = 0.1, 
            group.by = "Diagnosis",
            y.max = y_max, idents = c(nident)                                   # Lọc chỉ hiển thị dữ liệu của nhóm nident (ví dụ: "DNT cells").
            ) + stat_compare_means(comparisons = test_sign, label = "p.signif") # Tính toán và hiển thị giá trị p-value giữa các nhóm chẩn đoán.
  }
  plot_list <-list()
  y_max_list <- list()                                                          # Lưu giá trị tối đa trên trục Y cho từng gene.
  for(gene in gene_signature){                                                  # VÒNG LẶP ĐỂ VẼ BIỂU ĐỒ CHO TỪNG GENE. 
    plot_list[[gene]] <- plot_case1(gene) 
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]])                   # Truy cập dữ liệu thô ($data[[gene]]) từ biểu đồ ggplot để lấy giá trị lớn nhất của gene này.
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1))     #+1 được thêm vào để đảm bảo không bị che mất giá trị lớn nhất trên trục y của biểu đồ violin (VlnPlot).
  }
  cowplot::plot_grid(plotlist = plot_list)                     # Dùng plot_grid() từ thư viện cowplot để ghép tất cả biểu đồ thành một hình.
  file_name <- paste0(file_name, "_r.png")                     # Ghép thêm chuỗi "r_png" vào tên file.
  ggsave(file_name, width = 6, height = 8)
}

#+++++++++++++++++++++++++++++++++++
#+Tham số đầu vào:

#gene_signature: Danh sách các gene cần vẽ biểu đồ violin.

#file_name: Đường dẫn và tên file ảnh sẽ lưu kết quả.

#test_sign: Danh sách các nhóm so sánh để tính giá trị p-value.

#y_max: Giá trị tối đa của trục Y (giúp hiển thị đầy đủ giá trị p-value).

#nident: Nhóm tế bào cụ thể cần vẽ.
#+++++++++++++++++++++++++++++++++++


vp_case1 <- function(gene_signature, file_name, test_sign, y_max, nident) {
  plot_case1 <- function(signature, y_max = NULL) {
    p <- VlnPlot(CRER_Tcell, features = signature,
                 pt.size = 0, 
                 group.by = "Diagnosis", 
                 y.max = y_max, idents = c(nident)) + # Nếu có giá trị y_max, giới hạn trục y.
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") + # tính toán và vẽ các thống kê tóm tắt từ dữ liệu trong biểu đồ.Vẽ điểm trung bình của dữ liệu trong mỗi nhóm, dùng hình dạng 18 (hình tam giác) và màu đen.
      stat_compare_means(comparisons = test_sign, label = "p.format") 
    return(p)   #Biểu đồ sẽ không được lưu trữ nếu không có return(p), điều này giúp lưu trữ kết quả để sử dụng sau.
  }
  
  plot_list <- list() 
  y_max_list <- list()  #Khởi tạo danh sách rỗng để lưu biểu đồ và giá trị y_max.
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]], na.rm = TRUE) #na.rm = TRUE: bỏ qua giá trị NA.
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1)) #Vẽ lại biểu đồ với y_max mới để làm đẹp biểu đồ (+1 để tránh bị cắt phần trên).
  }
  
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 6, height = 8)
}

#Idents(CRER_Tcell) <- "integrated_snn_res.0.48"

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
########################################################################################################
########################################################################################################
base::library(package = "BiocParallel")         #Optionally multithread some steps
base::library(package = "DT")                   #Display nice table in HTML
base::library(package = "ggplot2")              #Draw graphs and plots
base::library(package = "ggpubr")               #Draw nicer graphs
base::library(package = "rstatix")              #Base R statistics
base::library(package = "knitr")                #Build  file HTML, PDF, Word from R Markdown (.Rmd) hoặc LaTeX (.Rnw).
base::library(package = "dplyr")                #Handle big tables
base::library(package = "Seurat")               #Handle SingCell analyses
base::library(package = "SeuratObject")         #Contains the data structures (object structure) of Seurat (assay, metadata, idents, slot, v.v.).
                                                #helping to store and manage single-cell data.
base::library(package = "SingleCellExperiment") #Handle SingleCell file formats.
base::library(package = "escape")               #Perform Gene set enrichment, pathway analysis
base::library(package = "clusterProfiler")      #Perform Gene ontology enrichment, pathway enrichment, GSEA analysis
base::library(package = "dittoSeq")             #Draw nice plots based on Seurat
base::library(package = "org.Hs.eg.db")         #Human genome annotation
base::library(package = "Cairo")                #Graphs library

CRER_Tcell <- base::readRDS(file = "./RData/06.CRER_Tcell_reclusters_annotated.rds")
DefaultAssay(CRER_Tcell) <- "RNA"
CRER_Tcell <- JoinLayers(CRER_Tcell)
DimPlot(CRER_Tcell, label = T, label.size = 4)

SUBset <- subset(CRER_Tcell, ident= c("DNT cells"))

GS <- getGeneSets(library = "H")                                   # Hiện tại library của escape đã được thay thế bằng collection từ msigdbr. 

GS <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") # đây cú pháp mới được thay thế. 

if (!requireNamespace("irGSEA", quietly = TRUE)) {                 #requireNamespace(): kiểm tra xem package irGSEA đã được cài đặt chưa.
  devtools::install_github("chuiqin/irGSEA", force = T)            #devtools::install_github("chuiqin/irGSEA") tải và cài đặt package từ GitHub thay vì CRAN/Bioconductor.
                                                                   #force = T ép buộc cài đặt lại package ngay cả khi đã có sẵn.
}
library(irGSEA)

DNT_cells <- irGSEA.score(object = SUBset, assay = "RNA", 
                          slot = "data", seeds = 123, ncores = 1,
                          min.cells = 3, min.feature = 0,
                          custom = F, geneset = NULL, msigdb = T, 
                          species = "Homo sapiens", category = "H",  
                          subcategory = NULL, geneid = "symbol",
                          method = c("AUCell", "UCell", "singscore", 
                                     "ssgsea", "JASMINE", "viper"),
                          aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

#Vì doMC không hoạt động trên Windows, bạn có thể dùng doParallel cho AUCell.
cl <- makeCluster(1)
registerDoParallel(cl)







