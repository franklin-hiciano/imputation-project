#!/bin/bash
set -e -x
set -euo pipefail

data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
results=/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data

mkdir -p ${data}/long_reads
mkdir -p ${data}/fastq1
mkdir -p ${data}/fastq2
mkdir -p ${data}/hg38_reference
mkdir -p ${data}/paf_using_correct_hg38
mkdir -p ${data}/franken_reference_dir

mapfile -t assemblies < <(find "${data}/assemblies/assemblies" -maxdepth 2 -mindepth 2 -name '*.vrk-ps-sseq.asm-hap[12].fasta.gz')

function get_samples {
	wget https://raw.githubusercontent.com/nefte48/Imputation_data/refs/heads/main/results/2025_08_27_1kg_IG_snps/sample_names_shortR_longR.txt?token=GHSAT0AAAAAADMTK43JVGEXL2SHIJLY2KKA2ICKLYQ -O ${results}/sample_names_shortR_longR.txt
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

function download_hg38 {
        wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O /sc/arion/work/hiciaf01/databases/references/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function index_hg38_with_minimap2 {
	module load minimap2
        minimap2 -d ${data}/hg38_reference/GRCh38_full_analysis_set_plus_decoy_hla.mmi ${data}/hg38_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
}
function align_assemblies_oscar {


  cat contig_counts.tsv | grep NA19036 | awk '$2 == "hap1"' | while read sample hap contigs
  do
	fname=$(basename ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz)
	job_name="${sample_${hap}}"
    	paf_output="${data}/paf_using_correct_hg38/${job_name}.paf"

 	bsub_command="module load minimap2 && \
      minimap2 -x asm5 -t 16 -c --secondary=no \
      ${data}/hg38_reference/GRCh38_full_analysis_set_plus_decoy_hla.mmi \"${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz\" \
      > \"${data}/paf_using_correct_hg38/${sample}.vrk-ps-sseq.asm-${hap}.paf\""

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

  done
}

function align_assemblies_local {
	module load minimap2
	while read -r sample hap contigs
do 
	minimap2 -x asm5 -t 16 -c --secondary=no "${data}/hg38_reference/GRCh38_full_analysis_set_plus_decoy_hla.mmi" "${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz" > "${data}/paf_using_correct_hg38/${sample}.vrk-ps-sseq.asm-${hap}.paf"
done < contig_counts.tsv; 
}
	

function combine_assemblies_with_franken_reference {
    while read sample hap contigs
    do
        cat /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/franken_reference/reference.fasta <(zcat ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz) > ${data}/assemblies/assemblies/${sample}.vrk-ps-sseq.asm-${hap}_combined_with_franken.fasta
    done < contig_counts.tsv
}



function lift_over {
	mkdir -p ${data}/lift_over
	cat contig_counts.tsv | while read sample hap contigs
  	do
		paftools.js liftover ${data}/paf_using_correct_hg38/${sample}.vrk-ps-sseq.asm-hap${hap}.paf hg38_ig_and_tcr_coordinates.bed > ${data}/lift_over/${sample}_${hap}.bed
	done
}
	




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
function download_franken_reference {
	wget http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta -O ${data}/franken_reference_dir/reference.fasta
}

function convert_cram_to_fastq {
	module load samtools
	short_read_crams=$(ls ${data}/*.final.cram)
	while read cram
do
	samtools fastq -1 ${data}/fastq1/$(basename ${cram}).fq -2 ${data}/fastq2/$(basename ${cram}).fq -0 /dev/null -s /dev/null -n ${cram}
done <<< "${short_read_crams}"
}

#get_samples
#count_contigs
#sum_contigs_per_sample
#download_hg38
#index_hg38_with_minimap2
#align_assemblies_oscar
#long_read_fofn
#make_globus_batch_file
#get_small_globus_batch_file
#get_long_reads_with_globus
#combine_assemblies_with_franken_reference
#download_hg38
align_assemblies_local
