# Cut-off for each gene: 
# 1. calculate the absoluate max expression and relevant max expression (75% qunatile + 1.5*IQR) for each indication. 
# 2. calculate the 95% quantile of relevant max and absolute max and 
# 3. choose the smaller one as the threshold

# 1. read in the data
norm_expr <- data.table::fread("data/Processed/norm_expr.csv")
meta_norm <- data.table::fread("data/Processed/meta_norm.csv")

# 2. Data analysis
# 2.1 Data integration
norm_expr <- data.table::melt(norm_expr, id.vars = "gene", variable.name = "sample_id", value.name = "TPM")
combined <- merge(meta_norm, norm_expr, by = "sample_id")
combined[, TPM := as.numeric(TPM)]

# 2.2. Cutoff calculation
relevant_max <- function(y){
    if(IQR(y)==0){
        background <- quantile(y, probs = 0.75)
        y <- y[y>background]
    }
    
    if(length(y)>1){
        return((quantile(y, probs = c(0.75)) + 1.5*IQR(y)))
    }else{
        return(background)
    }
}

threshold <- combined[, list(abs_max=max(TPM), rel_max=relevant_max(TPM)), by=list(gene, tissue)]

threshold <- threshold[, list(abs_max_95quantile = quantile(abs_max, probs = 0.95), rel_max_95quantile = quantile(rel_max, probs = 0.95)), by = list(tissue)]

threshold[, cutoff := data.table::fifelse(rel_max_95quantile < abs_max_95quantile, rel_max_95quantile, abs_max_95quantile)]

threshold[, c("rel_max_95quantile", "abs_max_95quantile") := NULL]

# 3. Output
data.table::fwrite(threshold, file = "data/Processed/cutoffs.csv")
