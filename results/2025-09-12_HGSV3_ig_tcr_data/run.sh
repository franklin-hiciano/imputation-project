#!/bin/bash
set -e -x
set -euo pipefail

data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
results=/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data
databases=/sc/arion/work/hiciaf01/databases/

mkdir -p ${data}/assemblies/
mkdir -p ${data}/short_reads/
mkdir -p ${data}/long_reads
mkdir -p ${data}/fastq
mkdir -p ${data}/hg38_reference
mkdir -p ${data}/hg38_aligned_to_assemblies
mkdir -p ${data}/franken_reference_dir
mkdir -p ${data}/long_reads_paf

#HUMAN REFERENCE GENOME (HGR)
function download_hg38 {
        wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O /sc/arion/work/hiciaf01/databases/references/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function subset_hg38 {
	module load samtools
	samtools faidx ${databases}/references/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
	samtools faidx ${databases}/references/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > ${databases}/references/GRCh38_reference_genome/hg38_subset.fa
		
}

#REFERENCE GENOMES (RG)
function get_IG_TCR_loci_on_hg38 {
	wget https://raw.githubusercontent.com/franklin-hiciano/imputation-project/refs/heads/main/results/2025-09-12_HGSV3_ig_tcr_data/hg38_ig_and_tcr_coordinates.bed -O ${results}/hg38_ig_and_tcr_coordinates.bed
}

function get_samples {
	wget https://raw.githubusercontent.com/nefte48/Imputation_data/refs/heads/main/results/2025_08_27_1kg_IG_snps/sample_names_shortR_longR.txt?token=GHSAT0AAAAAADMTK43JVGEXL2SHIJLY2KKA2ICKLYQ -O ${results}/sample_names_shortR_longR.txt
}

function make_globus_batch_file_for_assemblies {
  for batch in 20230818_verkko_batch1 20230927_verkko_batch2 20240201_verkko_batch3; do
	echo /1000g/ftp/data_collections/HGSVC3/working/${batch}/assemblies/ ${data#/sc/arion}/assemblies/
  done < ${results}/sample_names_shortR_longR.txt > ${results}/assemblies_batch_transfer.txt
}

function get_assemblies_with_globus {
        module load python
        
	EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2   # EMBL-EBI Public Data
        MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7

        globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC TSV transfer 0001" --batch  ${results}/assemblies_batch_transfer.txt
}

function count_contigs {
    assemblies_dir="${data}/assemblies/assemblies"
    individual_counts_file="contig_counts.tsv"

    mapfile -t fai_files < <(find "${data}/assemblies/assemblies" -maxdepth 2 -mindepth 2 -name '*.vrk-ps-sseq.asm-hap[12].fasta.gz.fai')
    
    echo -e "sample\thap\tcontigs" > "${individual_counts_file}"
    
    for fai_file in "${fai_files[@]}"; do
        base_file=$(basename "${fai_file}")
        sample="${base_file%%.*}"
        hap=$(echo "${base_file}" | grep -o -E 'hap[12]')
        num=$(wc -l < "${fai_file}")
        
        printf "%s\t%s\t%s\n" "${sample}" "${hap}" "${num}"
    done >> "${individual_counts_file}"
}

function sum_contigs_per_sample {
    infile="contig_counts.tsv"
    outfile="contig_counts_summed.tsv"

    echo -e "sample\ttotal_contigs" > "${outfile}"

    awk -F'\t' ' 
    NR > 1 {
      counts[$1] += $3
    }
    END {
      for (sample in counts)
        printf "%s\t%d\n", sample, counts[sample]
    }' "${infile}" | sort -k1,1 >> "${outfile}"
}

function align_assemblies { 
        while IFS=$'\t' read sample hap contigs
do
	job_name=${sample}_${hap}
	bsub_command="module load minimap2 && minimap2 -x asm5 -t 16 -c ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz ${databases}/references/GRCh38_reference_genome/hg38_subset.fa > ${data}/hg38_aligned_to_assemblies/${sample}_${hap}.paf"

    echo "Submitting ${job_name}..."
    bsub -J "${job_name}" \
      -P "acc_oscarlr" \
      -n "16" \
      -R "span[hosts=1]" \
      -R "rusage[mem=8000]" \
      -q express \
      -W 12:00 \
      -o "${data}/hg38_aligned_to_assemblies/${job_name}.out" \
      -e "${data}/hg38_aligned_to_assemblies/${job_name}.err" \
      "${bsub_command}"
done< <(tail -n +2 contig_counts.tsv)
}

function lift_over {
	module load minimap2
	mkdir -p ${data}/lift_over
	while read sample hap contigs
  	do
		paftools.js liftover ${data}/hg38_aligned_to_assemblies/${sample}_${hap}.paf hg38_ig_and_tcr_coordinates.bed > ${data}/lift_over/${sample}_${hap}.bed
		ls ${data}/lift_over/${sample}_${hap}.bed
	done < <(tail -n +2 contig_counts.tsv)
}

function ig_tcr_regions {
	module load samtools
while read sample hap count; do
	samtools faidx "${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz" $(awk '{print $1":"$2+1"-"$3}' "${data}/lift_over/${sample}_${hap}.bed") > "/sc/arion/projects/oscarlr/franklin/imputation/2025-09-12_HGSV3_ig_tcr_data/${sample}_${hap}.fa"
done < <(tail -n +2 contig_counts.tsv)
}

function count_contigs_within_regions {
while read sample hap contigs;
do
	awk -v OFS='\t' '{print $1, 0, $2}' ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz.fai > genome.fai.bed

	bedtools intersect -a hg38_ig_and_tcr_coordinates.bed -b genome.fai.bed -c > ig_tcr_contig_counts_${sample}_${hap}.bed
done < <(tail -n +2 contig_counts.tsv)
}

function download_franken_reference {
	wget http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta -O ${data}/franken_reference_dir/reference.fasta
}

function combine_assemblies_with_franken_reference {
    while read sample hap contigs
    do
        cat /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/franken_reference/reference.fasta <(zcat ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz) > ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta
    done < contig_counts.tsv
}
#SHORT READS

function get_short_reads_index_files {
	for index in https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index; do
		wget ${index} -O ${results}/$(basename ${index})	
	done
}

function short_read_fofn {
	for index in 1000G_698_related_high_coverage.sequence.index 1000G_2504_high_coverage.sequence.index
do
	while read -r sample
	do
		awk -F'\t' -v s="${sample}" '{if ($10 == s) printf "%s\n", $1}' ${results}/${index}
	done < sample_names_shortR_longR.txt
done > short_read_fofn.txt
}	

function make_globus_batch_file_for_short_reads {
while read -r path; do
	echo ${path#ftp://ftp.sra.ebi.ac.uk} ${data#/sc/arion}/short_reads/$(basename ${path})
done < ${results}/short_read_fofn.txt > ${results}/short_reads_batch_transfer.txt
}

function get_short_reads_with_globus {
        module load python

        EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2   # EMBL-EBI Public Data
        MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7

        globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC TSV transfer 0001" --batch  ${results}/short_reads_batch_transfer.txt --verify-checksum --sync-level checksum
}

function convert_cram_to_fastq {
	module load samtools
while read -r path_globus path_minerva; do
	cram=/sc/arion${path_minerva}
	job_name=$(basename ${cram})
	bsub_command="module load samtools && samtools fastq --reference ${data}/hg38_reference/hg38_subset.fa -1 ${data}/fastq/$(basename ${cram})_R1.fastq.gz -2 $(basename ${cram})_R2.fastq.gz -0 /dev/null -s /dev/null -n ${cram}"

        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "1" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/fastq/${job_name}.out" \
            -e "${data}/fastq/${job_name}.err" \
            "${bsub_command}"
done < ${results}/short_reads_batch_transfer.txt
}


function concat_long_reads {
    cat ${results}/sample_names_shortR_longR.txt | while read -r sample
    do
        job_name="${sample}_joined_long_reads"
        bsub_command="awk -v s='${sample}' '\$1 == s {print \$2}' long_read_fofn.txt | awk -F/ '{print \$NF}' | xargs -I {} cat ${data}/long_reads/{} > ${data}/long_reads/${job_name}.fasta.gz"
        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "1" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/long_reads/${job_name}.out" \
            -e "${data}/long_reads/${job_name}.err" \
            "${bsub_command}"
    done
}


#LONG READS
function long_read_fofn {
  # builds: sample<TAB>remote_path
  while read -r sample; do
    awk -F'\t' -v s="$sample" '
      $6 == s { printf "%s\t%s\n", s, $1 }
    ' "${results}/igsr_HGSVC3.tsv.tsv"
  done < sample_names_shortR_longR.txt > long_read_fofn.txt
}

function make_globus_batch_file {
  while IFS=$'\t' read -r sample path; do
	path_1000_genomes_endpoint=${path#ftp.sra.ebi.ac.uk}
	path_minerva=${data}/long_reads/$(basename -- "${path}")
	path_mssm_arion_endpoint=${path_minerva#/sc/arion}
    	printf "%s\t%s\n" "${path_1000_genomes_endpoint}" "${path_mssm_arion_endpoint}"
  done < ${results}/long_read_fofn.txt > ${results}/long_read_batch_transfer.txt
}

function get_small_globus_batch_file {
	head -n 10 ${results}/long_read_batch_transfer.txt > ${results}/long_read_batch_transfer_small.txt 
}

function get_long_reads_with_globus {
	module load python
	
	local SRC="14a0be5f-226c-49fe-b65f-dba083d67fc3"   # source collection UUID
	local DST="6621ca70-103f-4670-a5a7-a7d74d7efbb7"   # Minerva/Arion collection UUID
	
	EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2   # EMBL-EBI Public Data
    	MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7

  	globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC TSV transfer 0001" --batch  ${results}/long_read_batch_transfer_small.txt
}

function concat_long_reads {
    cat ${results}/sample_names_shortR_longR.txt | while read -r sample
    do
        job_name="${sample}_joined_long_reads"
	bsub_command="awk -v s='${sample}' '\$1 == s {print \$2}' long_read_fofn.txt | awk -F/ '{print \$NF}' | xargs -I {} cat ${data}/long_reads/{} > ${data}/long_reads/${job_name}.fasta.gz"
        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "1" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/long_reads/${job_name}.out" \
            -e "${data}/long_reads/${job_name}.err" \
            "${bsub_command}"
    done
}

function align_long_reads_to_new_franken_assemblies {
	while read sample hap contigs
do
        job_name=${sample}_long_aligned_to_${sample}_${hap}
        bsub_command="module load minimap2 && minimap2 -x asm5 -t 16 -c ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta ${data}/long_reads/${sample}_joined_long_reads.fasta.gz > ${data}/long_reads_paf/${job_name}.paf"
	echo "Submitting ${job_name}..."
	bsub -J "${job_name}" \
		-P "acc_oscarlr" \
		-n "16" \
		-R "span[hosts=1]" \
		-R "rusage[mem=8000]" \
		-q express \
		-W 12:00 \
		-o "${data}/paf_using_correct_hg38/${job_name}.out" \
		-e "${data}/paf_using_correct_hg38/${job_name}.err" \
		"${bsub_command}"
done < <(tail -n +2 contig_counts.tsv)
}

#download_hg38
#subset_hg38
#get_IG_TCR_loci_on_hg38
#get_samples
#make_globus_batch_file_for_assemblies
#get_assemblies_with_globus
#merge_assembly_dirs
#get_short_reads_index_files
#short_read_fofn
#make_globus_batch_file_for_short_reads
#get_short_reads_with_globus
#convert_cram_to_fastq

