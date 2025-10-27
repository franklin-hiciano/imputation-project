#!/bin/bash
set -e -x
set -euo pipefail

data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
results=/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data
hg38_reference_dir=${data}/hg38_reference
franken_reference_dir=${data}/franken_reference
assemblies_dir=${data}/assemblies/assemblies
paf_dir=${data}/paf_using_correct_hg38
fastq_read1_dir=${data}/fastq1
fastq_read2_dir=${data}/fastq2
long_reads_index=/sc/arion/work/hiciaf01/projects/imputation/results/2025-09-12_HGSV3_ig_tcr_data/igsr_HGSVC3.tsv.tsv
long_reads_dir=${data}/long_reads


mkdir -p ${long_reads_dir}
mkdir -p ${fastq_read1_dir}
mkdir -p ${fastq_read2_dir}
mkdir -p ${hg38_reference_dir}
mkdir -p ${paf_dir}
mkdir -p ${franken_reference_dir}

mapfile -t assemblies < <(find "${assemblies_dir}" -maxdepth 2 -mindepth 2 -name '*.vrk-ps-sseq.asm-hap[12].fasta.gz')

function count_contigs_fai {
    local assemblies_dir="${data}/assemblies/assemblies"
    local individual_counts_file="contig_counts.tsv"

    mapfile -t fai_files < <(find "${assemblies_dir}" -maxdepth 2 -mindepth 2 -name '*.vrk-ps-sseq.asm-hap[12].fasta.gz.fai')
    
    echo -e "sample\thap\tcontigs" > "${individual_counts_file}"
    
    for fai_file in "${fai_files[@]}"; do
        local base_file=$(basename "${fai_file}")
        local sample="${base_file%%.*}"
        local hap=$(echo "${base_file}" | grep -o -E 'hap[12]')
        local num=$(wc -l < "${fai_file}")
        
        printf "%s\t%s\t%s\n" "${sample}" "${hap}" "${num}"
    done >> "${individual_counts_file}"
}

function sum_contigs_per_sample {
    local infile="contig_counts.tsv"
    local outfile="contig_counts_summed.tsv"

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
	url=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
        wget ${url} -O ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function index_hg38_with_minimap2 {
	module load minimap2
        minimap2 -d ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function align_assemblies {

  local num_cores=16

  for assembly_file in "${assemblies[@]}"; do

    local fname=$(basename "${assembly_file}")
    local job_name="${fname%.fasta.gz}"
    local paf_output="${paf_dir}/${job_name}.paf"

    local bsub_command="module load minimap2 && \
      minimap2 -x asm5 -t ${num_cores} -c --secondary=no \
      ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi \"${assembly_file}\" \
      > \"${paf_output}\""

    echo "Submitting ${job_name}..."
    bsub -J "${job_name}" \
      -P "acc_oscarlr" \
      -n "${num_cores}" \
      -R "span[hosts=1]" \
      -R "rusage[mem=8000]" \
      -q express \
      -W 12:00 \
      -o "${paf_dir}/${job_name}.out" \
      -e "${paf_dir}/${job_name}.err" \
      "${bsub_command}"

  done
}
function long_read_fofn {
  # builds: sample<TAB>remote_path
  while read -r sample; do
    awk -F'\t' -v s="$sample" '
      $6 == s { printf "%s\t%s\n", s, $1 }
    ' "$long_reads_index"
  done < sample_names_shortR_longR.txt > long_read_fofn.txt
}

function make_globus_batch_file {
  while IFS=$'\t' read -r sample path; do
	path_1000_genomes_endpoint=${path#ftp.sra.ebi.ac.uk}
	path_minerva=${long_reads_dir}/$(basename -- "${path}")
	path_mssm_arion_endpoint=${path_minerva#/sc/arion}
    	printf "%s\t%s\n" "${path_1000_genomes_endpoint}" "${path_mssm_arion_endpoint}"
  done < long_read_fofn.txt > long_read_batch_transfer.txt
}
function get_long_reads_with_globus {
	module load python
	
	local SRC="14a0be5f-226c-49fe-b65f-dba083d67fc3"   # source collection UUID
	local DST="6621ca70-103f-4670-a5a7-a7d74d7efbb7"   # Minerva/Arion collection UUID
	
	EMBL_EBI_ENDPOINT=47772002-3e5b-4fd3-b97c-18cee38d6df2   # EMBL-EBI Public Data
    	MINERVA_ARION_ENDPOINT=6621ca70-103f-4670-a5a7-a7d74d7efbb7

  	globus transfer ${EMBL_EBI_ENDPOINT} ${MINERVA_ARION_ENDPOINT} --label "HGSVC TSV transfer 0001" --batch  long_read_batch_transfer.txt

}
function download_franken_reference {
	wget http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta -O ${franken_reference_dir}/reference.fasta
}

function convert_cram_to_fastq {
	module load samtools
	short_read_crams=$(ls ${data}/*.final.cram)
	while read cram
do
	samtools fastq -1 ${fastq_read1_dir}/$(basename ${cram}).fq -2 ${fastq_read2_dir}/$(basename ${cram}).fq -0 /dev/null -s /dev/null -n ${cram}
done <<< "${short_read_crams}"
}



#get_assemblies
#count_contigs_fai
#download_hg38
#index_hg38_with_minimap2
#align_assemblies
#download_franken_reference
#convert_cram_to_fastq
#align_assemblies_normally
#long_read_fofn
#download_with_globus
#make_globus_batch_file
#get_long_reads_with_globus
align_assemblies
