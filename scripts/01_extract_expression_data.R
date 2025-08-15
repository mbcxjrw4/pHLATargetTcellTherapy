# 1. load data
# tsv.gz file of gene id and name is brought into R.
UCSC_Path <- "data/Input/UCSC/"
file4 <- "TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv.gz"
gene <- data.table::fread(file = paste0(UCSC_Path, file4), check.names=FALSE, stringsAsFactors=F)

# TSV file of meta data is brought into R.
file5 <-  "TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv"
meta <- data.table::fread(file=paste0(UCSC_Path, file5), check.names=FALSE, stringsAsFactors=F)

# TSV file of integrated RNAseq data are brought into R.
file6 <- "TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv"
rna <- data.table::fread(file=paste0(UCSC_Path, file6), check.names=FALSE, stringsAsFactors=F)

# 2.Data analysis 
# ct antigen
ct_antigen_targets <- membrane_targets[membrane_targets %in% gene$gene]

# UCSC Xena
# remove data from TARGET
meta <- meta[Study!="TARGET"]
meta <- meta[Primary.site!=""]
meta <- meta[Diagnosis!=""]

# Thyroid
meta[Primary.site=="Thyroid Gland", Primary.site := "Thyroid"] 

# meta data
meta_norm <- meta[Diagnosis == "Normal"]
meta_norm[, c("Diagnosis", "Gender", "Study") := NULL]
data.table::setnames(meta_norm, old = c("Sample.ID", "Primary.site"), new = c("sample_id", "tissue"))

meta_tum <- meta[Diagnosis != "Normal"]
meta_tum[, c("Primary.site", "Gender", "Study") := NULL]
data.table::setnames(meta_tum, old = c("Sample.ID", "Diagnosis"), new = c("sample_id", "indication"))

# rna data
rna <- data.table::transpose(rna, keep.names = "Sample.ID")
data.table::setnames(rna, old=colnames(rna), new=c("Sample.ID", gene$gene))

# select genes
feature <- c("Sample.ID", ct_antigen)
rna <- rna[, ..feature]

norm_expr <- rna[Sample.ID %in% meta_norm$sample_id]
norm_expr <- data.table::transpose(norm_expr, keep.names = "Gene")
data.table::setnames(norm_expr, old = colnames(norm_expr), new = as.character(norm_expr[1,]))
norm_expr <- norm_expr[-1,]
rownames(norm_expr) <- norm_expr$Sample.ID
data.table::setnames(norm_expr, old = "Sample.ID", new = "gene")

tum_expr <- rna[Sample.ID %in% meta_tum$sample_id]
tum_expr <- data.table::transpose(tum_expr, keep.names = "Gene")
data.table::setnames(tum_expr, old = colnames(tum_expr), new = as.character(tum_expr[1,]))
tum_expr <- tum_expr[-1,]
rownames(tum_expr) <- tum_expr$Sample.ID
data.table::setnames(tum_expr, old = "Sample.ID", new = "gene")

# 3. Output
data.table::fwrite(meta_norm, file = "data/Processed/meta_norm.csv")
data.table::fwrite(meta_tum, file = "data/Processed/meta_tum.csv")
data.table::fwrite(norm_expr, file = "data/Processed/norm_expr.csv")
data.table::fwrite(tum_expr, file = "data/Processed/tum_expr.csv")
