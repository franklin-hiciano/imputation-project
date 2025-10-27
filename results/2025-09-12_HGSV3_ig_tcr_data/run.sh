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
	samtools fastq -1 paired1.fq -2 paired2.fq -0 /dev/null -s /dev/null -n input.bam > single_end.fqls ${data}/assemblies/assemblies | while read sample
	do
		
		for hap in hap1 hap2
		do
			num=`cat ${data}/assemblies/assemblies/${sample}/${sample}.vrk-ps-sseq.asm-${hap}.fasta.gz.fai | wc -l`
			printf "%s\t%s\t%s\n" "${sample}" "${hap}" "${num}"
		done
	done > contig_counts.tsv
	sum_contigs_per_sample
}
#
#function download_hg38 {
#	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
#}
#
#function index_hg38_with_minimap2 {
#	module load minimap2
#	minimap2 -d ${hg38_reference_dir}/GRCh38_no_alt.mmi ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#}
#
function download_hg38 {
	url=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
	wget ${url} -O ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function index_hg38_with_minimap2 {
       module load minimap2
       minimap2 -d ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
}

function split_assemblies_into_smaller_lists_for_parallelization {
	num_parts=24
	find ${assemblies_dir} -name "*hap*.fasta" > ${results}/all_assemblies.txt
	split -n l/${num_parts} --additional-suffix=.txt ${results}/all_assemblies.txt ${results}/assemblies_split_
}
function align_assemblies_normally {
    num_cores=48

    module load minimap2

    while read assembly
    do
		fname=$(basename "${assembly}")

        	paf="${paf_dir}/${fname}.paf"

        	minimap2 -x asm5 -t ${num_cores} -c --secondary=no \
        	${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi "${assembly}" \
        	> "${paf}"
    done < ${results}/all_assemblies.txt
}


function align_assemblies_parallelized {
    # Define the number of cores to use, matching the '-n' in bsub and '-t' in minimap2
    num_cores=24

    # Find the 12 partial file lists
    assemblies_lists=$(find ${results} -maxdepth 1 -name 'assemblies_split_*.txt')

    # Iterate over the list files, launching one job per file.
    while read list
    do
        # Use basename to create a descriptive job name and for output files
        job_name=$(basename "${list}" .txt)

        # 1. DEFINE THE COMMAND TO BE SENT TO BSUB
        # This command is a single shell script that will run on the cluster node.
        # Note the use of escaped variables (\$) for those that need to be evaluated *inside* the job.
        # The list file variable ('$list') is not escaped as it is evaluated *before* bsub submission.
        bsub_command="module load minimap2 && \
            while read assembly; do \
                # Evaluate assembly path (inside the job)
                fname=\$(basename \"\$assembly\"); \
                paf=\"${paf_dir}/\${fname}.paf\"; \
                
                # Run minimap2 using the allocated cores
                minimap2 -x asm5 -t ${num_cores} -c --secondary=no \
                ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi \"\$assembly\" \
                > \"\$paf\"; \
            done < ${list}"
            
        # 2. SUBMIT THE JOB USING THE DEFINED COMMAND VARIABLE
        # The ${num_cores} variable is used for both LSF resource allocation and minimap2 threading.
        bsub -J "${job_name}" \
        -P "acc_oscarlr" \
        -n "${num_cores}" \
        -R "span[hosts=1]" \
        -R "rusage[mem=6000]" \
        -q express \
        -W 12:00 \
        -o "${paf_dir}/${job_name}.out" \
        -e "${paf_dir}/${job_name}.err" \
        "${bsub_command}"
        
    done <<< "${assemblies_lists}"

}
#
#function align_assemblies_parallelized {
#    num_cores=16
#    mem=8000 # Memory in MB
#
#    assemblies_lists=$(find ${results} -maxdepth 1 -name 'assemblies_split_*.txt')
#
#    while read list
#    do
#        job_name=$(basename "${list}" .txt)
#
#        bsub_command="module load minimap2 && \
#            while read assembly; do \
#                # Evaluate assembly path (inside the job)
#                fname=\$(basename \"\$assembly\"); \
#                paf=\"${paf_dir}/\${fname}.paf\"; \
#
#                # Run minimap2 using the allocated cores
#                minimap2 -x asm5 -t ${num_cores} -c --secondary=no \
#                ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi \"\$assembly\" \
#                > \"\$paf\"; \
#            done < ${list}"
#
#        bsub -I \
#        -J "${job_name}" \
#        -P "acc_oscarlr" \
#        -n "${num_cores}" \
#        -R "span[hosts=1]" \
#        -R "rusage[mem=${mem}]" \
#        -q interactive \
#        -W 12:00 \
#        "${bsub_command}"
#
#    done <<< "${assemblies_lists}"
#
#    echo "Launched 12 parallel interactive jobs. Check 'bjobs' for status."
#}
#
function align_assemblies {
	module load minimap2
	
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



#
#function align_assemblies_parallelized {
#	module load minimap2
#	alignment_threads=18
#	num_jobs=1
#
#	MIN_BYTES="${MIN_BYTES:-1000000}"	
#	
#	for hap in hap1 hap2
#do
#	all_assemblies=$(ls -d ${assemblies_dir}/*/*${hap}.fasta)
#	while read assembly
#	do
#		echo "RUNNING ASSEMBLY ${assembly}"
#		job_name="align_${assembly}"
#
#		
#      		while [ "$(bjobs -w 2>/dev/null | grep -c 'align_')" -ge ${num_jobs} ]; do
#        		sleep 10
#      		done
#
#		fname=$(basename ${assembly})
#		paf="${paf_dir}/${fname}.paf"
#		if [ -f "${paf}" ]; then
#		  size=$(wc -c < "${paf}" 2>/dev/null || echo 0)
#		  if [ "${size}" -ge "${MIN_BYTES}" ]; then
#		    echo "[SKIP] ${fname}: existing PAF ${size} bytes ≥ ${MIN_BYTES}"
#		    continue
#		  else
#		    echo "[REDO] ${fname}: tiny PAF ${size} bytes < ${MIN_BYTES} — resubmitting"
#		    rm -f "${paf}"
#		  fi
#		fi
#		
#		out=${paf_dir}/${fname}.out
#		err=${paf_dir}/${fname}.err
#		
#		rm -f ${out} ${err} ${paf}
#		bsub -J "${job_name}" \
#		-P "acc_oscarlr" \
#           	-n "${alignment_threads}" \
#           	-R "span[hosts=1]" \
#		-R "rusage[mem=6000]" \
#           	-q express \
#           	-W 4:00 \
#		-o "${out}" \
#		-e "${err}" \
#          	 "module load minimap2 && minimap2 -x asm5 -t ${alignment_threads} -c --secondary=no \
#        	    ${hg38_reference_dir}/GRCh38_full_analysis_set_plus_decoy_hla.mmi ${assembly} \
#		    > ${paf}"
#	done <<< "${all_assemblies}"
#done
#}
#

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
split_assemblies_into_smaller_lists_for_parallelization
align_assemblies_parallelized
