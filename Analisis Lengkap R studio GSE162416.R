#Modul: Analisis Ekspresi Gen makrofag terhadap Mycobacterium leprae 
#Dataset: GSE162416 (Infected Macrofag vs Normal)
#Platform: Microarray (Illumina GPL10558)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 

#Author: Rhizna Juniarti
##############################################################################
Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(gplots)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 


#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE162416", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]


#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PRE-PROCESSING DATA EKSPRESI
ex <- exprs(gset)

#Quantile untuk interpretasi log fold change
qx <- as.numeric(quantile(ex,c(0,0.25,0.5,0.75,0.99,1), na.rm= TRUE))

#LogTransform adalah variabel logika (TRUE/FALSE)
#Operator logika:
#>  : lebih besar dari 
#|| : OR (atau)
#&& : AND (dan)

LogTransform<- (qx[5]>100)|| (qx[6]-qx[1]> 50 && qx[2]>0)


#IF Statement:
#Jika LogTransform= TRUE, maka lakukan log2
if (LogTransform) {
  #Nilai <= 0 tidak boleh di -log, maka diubah menjadi NA
  ex[ex<= 0]<- NA
  ex <- log2(ex)
}

#PART E. DEFINISI KELOMPOK SAMPEL 

#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]

#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)

#factor():
#Mengubah data kategorik menjadi faktor
#Faktor sangat penting untuk analisis statistik di R
gset$group <- factor(groups)

#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)

#########################################
#DESIGN MATRIX (KERANGKA STATISTIK)
#########################################
#model.matrix():
#membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design<- model.matrix(~0 + gset$group)

#colnames(): memberi nama kolom agar mudah dibasa
colnames(design)<- levels(gset$group)

#menentukan perbandingan biologis
grup_terinfeksi<- nama_grup[2]
grup_normal<- nama_grup[1]

contrast_formula<- paste(grup_terinfeksi, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

#################################
#Analisis Differential Expression (LIMMA)
#################################

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01  -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)

#PART H. ANOTASI NAMA GEN 

#Penting:
#Pada data microarray Illumina, unit analisis awal adalah PROBE,
#bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan
#database resmi Bioconductor.

#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db, 
  keys = probe_ids, 
  columns = c("SYMBOL", "GENENAME"), 
  keytype = "PROBEID"
)

head(probe_ids)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI 

#Boxplot digunakan untuk:
#- Mengecek distribusi nilai ekspresi antar sampel
#- Melihat apakah ada batch effect
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar

#Set warna berdasarkan grup
#Set warna berdasarkan grup
# 1. Definisikan palet warna manual
my_palette <- c("grey", "indianred") 

# 2. Mapping warna ke data
group_colors <- my_palette[as.numeric(gset$group)]

# 3. Gambar Boxplot
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

# 4. Legend 
legend(
  "topright",
  legend = levels(gset$group),
  fill = my_palette,
  cex = 0.8
)
#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )


#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)

#UMAP digunakan untuk:
#- Mereduksi ribuan gen menjadi 2 dimensi
#- Melihat pemisahan sampel secara global
#- Alternatif PCA (lebih sensitif ke struktur lokal)

#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )


#PART J.1 VISUALISASI VOLCANO PLOT 

#Volcano plot menggabungkan:
#- Log fold change (efek biologis)
#- Signifikansi statistik

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Makrofag vs M. leprae")

#PART J.2 VISUALISASI HEATMAP 

#Heatmap digunakan untuk melihat pola ekspresi gen
#antar sampel berdasarkan gen-gen paling signifikan

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)


#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

#PART K. MENJAWAB SOAL 
#1. Gen apa saja yang mengalami upregulation dan downregulation, ditampilkan dalam bentuk volcano plot.
#2. 50 Differentially Expressed Genes (DEGs) teratas, ditampilkan dalam bentuk heatmap.
#3. Analisis enrichment, yang mencakup:
#- Gene Ontology (GO)
#- KEGG Pathway: Hasil analisis enrichment wajib disertai visualisasi/plot yang relevan.

#Gunakan Filter Packages dplyr
#Gen upregulated dan downregulated 

# 1. Filter yang Upregulated (Naik)
# Syarat: logFC > 1 dan adj.P.Val < 0.01
# Ambil daftar gen yang Upregulated (UP)
# Ambil daftar gen yang Upregulated (UP)
gen_up <- volcano_data %>% 
  filter(status == "UP") %>% 
  dplyr::select(Gene, logFC, adj.P.Val) %>% # Menambahkan dplyr:: di depan select
  arrange(desc(logFC))

# Lihat hasilnya (50 gen teratas)
head(gen_up, 50)
#Jumlah gen upregulated
nrow(gen_up)

# Ambil daftar gen yang Downregulated (DOWN)
gen_down <- volcano_data %>% 
  filter(status == "DOWN") %>% 
  dplyr::select(Gene, logFC, adj.P.Val) %>% 
  arrange(logFC)

#50 gen downregulated
head(gen_down, 50)
#jumlah gen downregulated
nrow(gen_down)



#Jawaban no 3
#load dulu packages untuk pathway
library(clusterProfiler)
library(org.Hs.eg.db) # Database untuk Gen Manusia
library(enrichplot)

# Ambil semua gen yang signifikan (tidak termasuk status 'NO')
genes_to_test <- volcano_data$Gene[volcano_data$status != "NO"]

# Jalankan analisis GO
ego <- enrichGO(gene          = genes_to_test,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL', # nama gen (Symbol)
                ont           = "BP",     # BP = Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)

# Visualisasi GO (Dotplot)
dotplot(ego, showCategory=15) + 
  ggtitle("Gene Ontology - Biological Process")

# A. Konversi SYMBOL ke ENTREZID
gene_conv <- bitr(genes_to_test, 
                  fromType = "SYMBOL", 
                  toType   = "ENTREZID", 
                  OrgDb    = org.Hs.eg.db)

# B. Jalankan analisis KEGG
kk <- enrichKEGG(gene         = gene_conv$ENTREZID,
                 organism     = 'hsa', # 'hsa' artinya Homo sapiens (Manusia)
                 pvalueCutoff = 0.05)

# C. Visualisasi KEGG (Barplot)
dotplot(kk, showCategory=15) + 
  ggtitle("KEGG Pathway Analysis")

#PART L. MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_Analisis_GSE162416.csv")

#Menggabungkan data yang sebelumnya belum dilengkapi analisis GO dan KEGG
#Menambahkan data Up regulated & down regulated dalam satu file yang sama namun berbeda sheet

#Load Packages 
install.packages("openxlsx")
library("openxlsx")

# 1. Buat Workbook baru
wb <- createWorkbook()

# Sheet 1: DATA KESELURUHAN ((telah di download sebelumnya))
# Masukkan data asli kamu (misal namanya topTableResults)
addWorksheet(wb, "Hasil_Analisis_GSE162416")
writeData(wb, "Hasil_Analisis_GSE162416", topTableResults)

# Sheet 2: GEN UP REGULATED
addWorksheet(wb, "Gen_UP_Regulated")
writeData(wb, "Gen_UP_Regulated", gen_up)

# Sheet 3: GEN DOWN REGULATED
addWorksheet(wb, "Gen_DOWN_Regulated")
writeData(wb, "Gen_DOWN_Regulated", gen_down)

# Sheet 4: GENE ONTOLOGY (GO)
addWorksheet(wb, "Gene_Ontology_Results")
writeData(wb, "Gene_Ontology_Results", as.data.frame(ego))

# Sheet 5: KEGG PATHWAY
addWorksheet(wb, "KEGG_Pathway_Results")
writeData(wb, "KEGG_Pathway_Results", as.data.frame(kk))

# 2. Simpan ke dalam satu file
saveWorkbook(wb, "Analisis_Lengkap_GSE162416.xlsx", overwrite = TRUE)

message("Semua data (Lengkap, UP, DOWN, GO, KEGG) sudah tersimpan dalam satu file Excel!")

message("Analisis selesai. File hasil telah disimpan.")