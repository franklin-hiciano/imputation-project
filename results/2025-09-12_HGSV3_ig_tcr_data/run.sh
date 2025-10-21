#!/bin/bash
set -e -x
set -euo pipefail

assemblies_dir=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/assemblies/assemblies
sample_assembly_names=$(ls ${assemblies_dir}/)

data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long

function get_assemblies {
	pass	# acquired through globus
}

function count_contigs {
	sample=HG00096
	gunzip -c ${assemblies_dir}/${sample}/${sample}.vrk-ps-sseq.asm-hap1.fasta.gz | grep -c ">"
	gunzip -c ${assemblies_dir}/${sample}/${sample}.vrk-ps-sseq.asm-hap2.fasta.gz | grep -c ">"
}

function count_contigs_multithreaded {
	sample=HG00096
	copy=${assemblies_dir}/${sample}/${sample}.vrk-ps-sseq.asm-hap1_copy.fasta.gz
	f=${assemblies_dir}/${sample}/${sample}.vrk-ps-sseq.asm-hap1_copy.fasta
	cp ${assemblies_dir}/${sample}/${sample}.vrk-ps-sseq.asm-hap1.fasta.gz ${copy}
	unpigz ${copy}
	grep -c ">" ${f}
}

function count_all_contigs_multithreaded {
	counts=${assemblies_dir}/contig_counts.txt
	touch ${counts}
	echo "" > ${counts}
	for hap in hap1 hap2
do
	all_assemblies=$(ls -d /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/assemblies/assemblies/*/*${hap}.fasta.gz)
	while read assembly
	do
		unpigz -kf ${assembly}
		num_contigs=$(grep -c "^>" ${assembly%.gz})
		printf "%s\t%s\n" "${assembly}" "${num_contigs}" >> "$counts"

	done <<< "${all_assemblies}"
done
}

function count_contigs_fai {
	ls ${data}/assemblies/assemblies | while read sample
	do
		for hap in hap1 hap2
		do
			num=`cat ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz.fai | wc -l`
			echo "${sample}\t${hap}\t${num}"
		done
	done > count.txt
}
	

function download_hg38 {
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O ${assemblies_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 

}


count_all_contigs_multithreaded
function count_assemblies {
	ls /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/assemblies/assemblies | wc -l
}


function get_hg38 {
	out=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/hg38.fa.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O ${out} 
	gunzip ${out}
}

function align_assemblies_to_hg38 {
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip hg38.fa.gz	
}

#count_assemblies
#count_contigs
#count_contigs_multithreaded
#count_all_contigs_multithreaded
download_hg38
