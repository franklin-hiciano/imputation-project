#!/bin/bash
set -e -x
set -euo pipefail

THREADS=24 #same as bsub
data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
hg38_reference_dir=${data}/hg38_reference
franken_reference_dir=${data}/franken_reference
assemblies_dir=${data}/assemblies/assemblies
paf_dir=${data}/paf

function get_assemblies {
	pass	# acquired through globus
}

function sum_contigs_per_sample {
  infile="contig_counts.tsv"
  outfile="contig_counts_summed.tsv"
  echo -e "sample\ttotal_contigs" > ${outfile}

  # Use awk to sum hap1 + hap2 per sample
  awk -F'\t' '
  {
    counts[$1] += $3
  }
  END {
    for (sample in counts)
      printf "%s\t%d\n", sample, counts[sample]
  }' ${infile} | sort -k1,1 >> ${outfile}
}
	
function count_contigs_fai {
	ls ${data}/assemblies/assemblies | while read sample
	do
		
		for hap in hap1 hap2
		do
			num=`cat ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz.fai | wc -l`
			printf "%s\t%s\t%s\n" "${sample}" "${hap}" "${num}"
		done
	done > contig_counts.tsv
	sum_contigs_per_sample
}

function download_hg38 {
	mkdir -p ${hg38_reference_dir}
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
}

function index_hg38_with_minimap2 {
	module load minimap2
	minimap2 -d ${hg38_reference_dir}/GRCh38_no_alt.mmi ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
}


function align_assemblies {
	module load minimap2
	mkdir -p ${paf_dir}
	
	for hap in hap1 hap2
do
	all_assemblies=$(ls -d ${assemblies_dir}/*/*${hap}.fasta)
	while read assembly
	do
		echo "RUNNING ASSEMBLY ${assembly}"
		minimap2 -x asm5 -t "$THREADS" -c --secondary=no "${data}/hg38_reference/GRCh38_no_alt.mmi" "${assembly}" > "${paf_dir}/$(basename ${assembly}).paf"
	done <<< "${all_assemblies}"
done
}

function download_franken_reference {
	
	mkdir ${franken_reference_dir}
	wget http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta -O ${franken_reference_dir}
}

#get_assemblies
#count_contigs_fai
#download_hg38
#index_hg38_with_minimap2
#align_assemblies
download_franken_reference
