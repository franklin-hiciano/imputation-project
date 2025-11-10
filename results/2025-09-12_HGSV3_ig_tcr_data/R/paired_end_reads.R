results <- "/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data"
counts <- read.csv(file.path(results, "paired_end_reads_counts.tsv"), sep='\t', header=TRUE)
#THIS IS INCORRECT COUNTING
if (all(counts$R1==counts$R2)) {
  png(paste0("plots/", region_name, "_contigs_barplot.png"), width = 1200, height = 300)
  barplot(counts$R1, names.arg=counts$sample, las = 2, cex.names = 0.6, xlab="sample", ylab="reads", main="1000G short-read paired-end read counts for 66 samples")
}
