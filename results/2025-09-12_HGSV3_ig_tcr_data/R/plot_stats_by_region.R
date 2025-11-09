results <- "/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data"

#1. COUNT HOW MANY IN EACH REGION PER SAMPLE
plot_contig_counts_by_region <- function() {
  df <- read.csv(file.path(results, "contig_counts_by_region.tsv"), sep='\t')
  df_counts <- df[, 3:ncol(df)]
  
  for (region_name in names(df_counts)[3:ncol(df_counts)]) {
    region <- df_counts[[region_name]]
    region_name <- colnames(df_counts)[i]
    png(paste0("plots/", region_name, "_contigs_barplot.png"), width = 1200, height = 300)
    barplot(region, names.arg=paste(df$sample, df$hap, sep="_"), las = 2, cex.names = 0.6, ylab="# contigs", xlab="sample assembly", main=paste0(region_name, " contigs in subsetted assemblies"))
    dev.off()
  }
}

plot_num_bases_by_region <- function() {
  df <- read.csv(file.path(results, "num_bases_by_region.tsv"), sep='\t')
  
  df_counts <- df[, 3:ncol(df)]
  for (region_name in names(df_counts)) {
    region <- df_counts[[region_name]]
    par(mar = c(10, 4, 4, 2) + 0.1)
    png(paste0("plots/", region_name, "_num_bases_barplot.png"), width = 1200, height = 300)
    barplot(region, names.arg=paste(df$sample, df$hap, sep="_"), las = 2, cex.names = 0.6, ylab="# bases", xlab="sample assembly", main=paste0("# of ", region_name, " bases in subsetted assemblies"))
    dev.off()
  }
}
plot_contig_counts_by_region()
plot_num_bases_by_region()




#GET POTENTIALLY PROBLEMATIC SAMPLES WITH LOW PVALUES
make_contig_counts_pvalues_by_region <- function() {
  df <- read.csv(file.path(results, "contig_counts_by_region.tsv"), sep = '\t')
  
  df_counts <- df[, 3:ncol(df)]
  pvalues <- data.frame(sample = df$sample, hap = df$hap)
  
  for (region_name in names(df_counts)) {
    region <- df_counts[[region_name]]
    mean_val <- mean(region, na.rm = TRUE)
    sd_val <- sd(region, na.rm = TRUE)
    p_col <- if (sd_val == 0) rep(1, length(region)) else 2 * pnorm(-abs((region - mean_val) / sd_val))
    pvalues[[region_name]] <- p_col
  }
  
  write.table(
    pvalues,
    file.path(results, "contig_counts_pvalues_by_region.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
make_contig_counts_pvalues_by_region()

make_num_bases_pvalues_by_region <- function() {
  df <- read.csv(file.path(results, "num_bases_by_region.tsv"), sep = '\t')
  
  df_counts <- df[, 3:ncol(df)]
  pvalues <- data.frame(sample = df$sample, hap = df$hap)
  
  for (region_name in names(df_counts)) {
    region <- df_counts[[region_name]]
    mean_val <- mean(region, na.rm = TRUE)
    sd_val <- sd(region, na.rm = TRUE)
    p_col <- 2 * pnorm(-abs((region - mean_val) / sd_val))
    pvalues[[region_name]] <- p_col
  }
  
  write.table(pvalues, file.path(results, "num_bases_pvalues_by_region.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
make_num_bases_pvalues_by_region()

