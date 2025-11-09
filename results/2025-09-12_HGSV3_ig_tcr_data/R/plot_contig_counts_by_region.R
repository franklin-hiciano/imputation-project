results <- "/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data"
df <- read.csv(file.path(results, "contig_counts_by_region.tsv"), sep='\t')
df_counts <- df[, 3:ncol(df)]
for (i in 1:ncol(df_counts)) {
  region <- df_counts[, i]
  region_name <- colnames(df_counts)[i]
  png(paste0("plots/", region_name, "_contigs_barplot.png"), width = 1200, height = 300)
  barplot(region, names.arg=paste(df$sample, df$hap, sep="_"), las = 2, cex.names = 0.6, ylab="# contigs", xlab="sample assembly", main=paste0(region_name, " contigs in subsetted assemblies")
  dev.off()
}

