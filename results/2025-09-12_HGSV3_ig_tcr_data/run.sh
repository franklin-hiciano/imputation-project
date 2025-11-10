#!/bin/bash
set -e -x
set -euo pipefail

data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
results=/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data
databases=/sc/arion/work/hiciaf01/databases/

mkdir -p ${data}/aligned_short_reads_aligned/
mkdir -p ${data}/aligned_short_reads/
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

function make_region_names_file {
	awk '{printf "%s_%s_%s\n", $1,$2,$3}' hg38_ig_and_tcr_coordinates.bed > region_names.txt
}

function count_contigs_within_regions {

    echo -e "sample\thap\t$(cut -f4 hg38_ig_and_tcr_coordinates.bed | paste -s -d '\t')" > contig_counts_by_region.tsv

    while read sample hap contigs

do
        while read region
	do
		echo $(grep ${region} ${data}/lift_over/${sample}_${hap}.bed | wc -l)
	done < region_names.txt > region_counts_${sample}_${hap}.txt

	printf "%s\t%s\t%s\n" "${sample}" "${hap}" "$(paste -s region_counts_${sample}_${hap}.txt)"
done < <(tail -n +2 contig_counts.tsv) >> contig_counts_by_region.tsv

}

function count_bases_within_regions {
	echo -e "sample\thap\t$(cut -f4 hg38_ig_and_tcr_coordinates.bed | paste -s -d '\t')" > num_bases_by_region.tsv
	while read sample hap contigs
do
	while read region
        do
		echo $(awk -v s="${region}" 'index($4, s) > 0 {sum += ($3-$2)} END {print sum}' ${data}/lift_over/${sample}_${hap}.bed)
        done < region_names.txt > num_bases_${sample}_${hap}.txt

	printf "%s\t%s\t%s\n" "${sample}" "${hap}" "$(paste -s num_bases_${sample}_${hap}.txt)"
done < <(tail -n +2 contig_counts.tsv) >> num_bases_by_region.tsv
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
while read -r sample hap counts; do
	cram=${data}/short_reads/${sample}.final.cram
	job_name=$(basename ${cram})
	bsub_command="module load samtools && samtools fastq --reference ${databases}/references/GRCh38_reference_genome/hg38_subset.fa -1 ${data}/fastq/$(basename ${cram})_R1.fastq.gz -2 ${data}/fastq/$(basename ${cram})_R2.fastq.gz -0 /dev/null -s /dev/null -n ${cram}"

        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "16" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/fastq/${job_name}.out" \
            -e "${data}/fastq/${job_name}.err" \
            "${bsub_command}"
done < <(tail -n +2 contig_counts.tsv)
}

function index_combined_assemblies {
while read -r sample hap counts; do
        job_name=${sample}_${hap}_frnk
        bsub_command="module load bwa && bwa index ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta"

        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "16" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/aligned_short_reads/${job_name}.out" \
            -e "${data}/aligned_short_reads/${job_name}.err" \
            "${bsub_command}"
    done < <(tail -n +2 contig_counts.tsv)
}


function align_short_reads {
	   module load samtools
while read -r sample hap counts
do
	for hap in hap1 hap2
	do
		cram=${data}/short_reads/${sample}.final.cram
	        job_name=${sample}_${hap}_SR_aligned
        	bsub_command="module load bwa && bwa mem -t 16 ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta ${data}/fastq/$(basename ${cram})_R1.fastq ${data}/fastq/$(basename ${cram})_R2.fastq > ${data}/short_reads_aligned/${job_name}.sam"
	

	        bsub -J "${job_name}" \
            -P "acc_oscarlr" \
            -n "16" \
            -R "span[hosts=1]" \
            -R "rusage[mem=8000]" \
            -q express \
            -W 12:00 \
            -o "${data}/short_reads_aligned_sam/${job_name}.out" \
            -e "${data}/short_reads_aligned_sam/${job_name}.err" \
            "${bsub_command}"
	done
done < <(tail -n +2 contig_counts.tsv)
}	

function subset_short_reads {
	samtools fastq ${data}/short_reads_aligned/.sam	
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

function make_globus_transfer_file_for_long_reads() {
    local sample="$1"

    while IFS=$'\t' read -r path; do
        path_1000_genomes_endpoint="${path#*ftp.sra.ebi.ac.uk}"
        path_minerva="${data}/long_reads/$(basename -- "$path")"
        path_mssm_arion_endpoint="${path_minerva#/sc/arion}"
        printf '%s\t%s\n' "$path_1000_genomes_endpoint" "$path_mssm_arion_endpoint"
    done < <(awk -v s="$sample" '$1==s {print $2}' "$results/long_read_fofn.txt") > "${results}/${sample}_long_read_transfer.txt"
}



function get_sample_long_reads_with_globus() {
	module load python
	sample=$1
        while IFS=$'\t' read -r sample path
do
        path_1000_genomes_endpoint=${path#ftp.sra.ebi.ac.uk}
        path_minerva=${data}/long_reads/$(basename -- "${path}")
        path_mssm_arion_endpoint=${path_minerva#/sc/arion}
        printf "%s\t%s\n" "${path_1000_genomes_endpoint}" "${path_mssm_arion_endpoint}"
done < <(awk -v s=${sample} '$1 == s {print "$1"\t"$2"}' ${results}/long_read_fofn.txt) > ${results}/${sample}_${hap}_long_read_transfer.txt
		
	EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2
	MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7

  	globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC TSV transfer 0001" --batch ${results}/${sample}_${hap}_long_read_transfer.txt
}

get_sample_long_reads_with_globus() {
  module load python 
  local sample="$1"
  local out="${results}/${sample}_long_read_transfer.txt"

  awk -v s="$sample" '$1==s {print $1 "\t" $2}' "$results/long_read_fofn.txt" > ${sample}_long_read_fofn.txt
  while IFS=$'\t' read -r _ path; do
    printf '%s\t%s\n' \
      "${path#*ftp.sra.ebi.ac.uk}" \
      "${data#/sc/arion}/long_reads/$(basename -- "$path")"
  done < ${sample}_long_read_fofn.txt > "${results}/${sample}_long_read_transfer.txt"

  EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2
  MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7
  globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC long reads: ${sample}" --batch "${results}/${sample}_long_read_transfer.txt"
}


function concat_long_reads() {
  while read -r sample; do
    job_name="${sample}_joined_long_reads"
    xfer_file="${results}/${sample}_long_read_transfer.txt"

    # if the transfer list doesn't exist yet, build it first
    [[ -s "$xfer_file" ]] || get_sample_long_reads_with_globus "$sample"

    bsub -J "$job_name" \
      -P "acc_oscarlr" \
      -n "1" \
      -R "span[hosts=1]" \
      -R "rusage[mem=8000]" \
      -q express \
      -W 12:00 \
      -o "${data}/long_reads/${job_name}.out" \
      -e "${data}/long_reads/${job_name}.err" \
      bash -lc "cut -f2 '$xfer_file' | xargs -r -I{} cat '{}' > '${data}/long_reads/${job_name}.fasta.gz'"
  done < "${results}/sample_names_shortR_longR.txt"
}


function concat_long_reads() {
	sample=$1
	job_name="${sample}_joined_long_reads"
        bsub_command="cat ${sample}_long_read_transfer.txt | cut -f2 | xargs -I {} cat ${data}/long_reads/{} > ${data}/long_reads/${job_name}.fasta.gz"
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

}

function align_long_read_to_combined_assembly() {
	sample=$1
	for hap in hap1 hap2
do
	job_name=${sample}_${hap}
	bsub_command="module load minimap2 && minimap2 -x asm5 -t 16 -c ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta ${data}/long_reads/${sample}_joined_long_reads.fasta.gz > ${data}/long_reads_paf/${job_name}.paf"
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
done
}

function subset_long_reads {
	samtools fastq	
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

#make_region_names_file
#count_contigs_within_regions
#count_bases_within_regions
#index_combined_assemblies
align_short_reads

#make_globus_transfer_file_for_long_reads HG00171
#get_sample_long_reads_with_globus HG00171




#align_long_read_to_combined_assembly HG00096
