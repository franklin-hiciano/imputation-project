#!/bin/bash
set -e -x
set -euo pipefail

THREADS=24 #same as bsub
data=/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long
hg38_reference_dir=${data}/hg38_reference
franken_reference_dir=${data}/franken_reference
assemblies_dir=${data}/assemblies/assemblies
paf_dir=${data}/paf
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

function download_hg38 {
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
}

function index_hg38_with_minimap2 {
	module load minimap2
	minimap2 -d ${hg38_reference_dir}/GRCh38_no_alt.mmi ${hg38_reference_dir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
}


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


function align_assemblies_parallelized {
	module load minimap2
	alignment_threads=24
	num_jobs=12

	MIN_BYTES="${MIN_BYTES:-1000000}"	
	
	for hap in hap1 hap2
do
	all_assemblies=$(ls -d ${assemblies_dir}/*/*${hap}.fasta)
	while read assembly
	do
		echo "RUNNING ASSEMBLY ${assembly}"
		job_name="align_${assembly}"

		
      		while [ "$(bjobs -w 2>/dev/null | grep -c 'align_')" -ge ${num_jobs} ]; do
        		sleep 10
      		done

		fname=$(basename ${assembly})
		paf="${paf_dir}/${fname}.paf"
		if [ -f "${paf}" ]; then
		  size=$(wc -c < "${paf}" 2>/dev/null || echo 0)
		  if [ "${size}" -ge "${MIN_BYTES}" ]; then
		    echo "[SKIP] ${fname}: existing PAF ${size} bytes ≥ ${MIN_BYTES}"
		    continue
		  else
		    echo "[REDO] ${fname}: tiny PAF ${size} bytes < ${MIN_BYTES} — resubmitting"
		    rm -f "${paf}"
		  fi
		fi

		bsub -J "${job_name}" \
		-P "acc_oscarlr" \
           	-n "${alignment_threads}" \
           	-R "span[hosts=1]" \
		-R "rusage[mem=6000]" \
           	-q express \
           	-W 12:00 \
		-o "${paf_dir}/${fname}.out" \
		-e "${paf_dir}/${fname}.err" \
          	 "module load minimap2 && minimap2 -x asm5 -t ${alignment_threads} -c --secondary=no \
        	    ${data}/hg38_reference/GRCh38_no_alt.mmi ${assembly} \
		    > ${paf}"
	done <<< "${all_assemblies}"
done
}

function long_read_fofn {
  # builds: sample<TAB>remote_path
  while read -r sample; do
    fn=$(awk -F'\t' -v s="$sample" '$6 == s {print $1}' "$long_reads_index")
    printf "%s\t%s\n" "$sample" "$fn"
  done < sample_names_shortR_longR.txt > long_read_fofn.txt
}



function make_batch_file {
	
	
	$(basename file)
}


function download_with_globus {
  # usage: download_with_globus
  # requires: long_reads_index, long_reads_dir set; endpoints below set
  long_read_fofn

  local SRC="14a0be5f-226c-49fe-b65f-dba083d67fc3"   # source collection UUID
  local DST="6621ca70-103f-4670-a5a7-a7d74d7efbb7"   # Minerva/Arion collection UUID

  module load python
  globus whoami >/dev/null 2>&1 || globus login

  mkdir -p "${long_reads_dir}"

  # read sample and path, skip blanks
  while IFS=$'\t' read -r sample remotepath; do
    [[ -z "$sample" || -z "$remotepath" || "$remotepath" == "#"* ]] && continue

    # per-sample destination dir
    dest="${long_reads_dir%/}/${sample}"
    mkdir -p "$dest"

    # submit transfer (recursive in case paths are dirs)
    echo "Submitting transfer for ${sample}: ${remotepath} -> ${dest}/"
    globus ls ${SRC}:${remotepath}
    globus transfer \
      --label "long_reads_${sample}_$(date +%Y%m%d_%H%M%S)" \
      --recursive \
      --preserve-mtime \
      --verify-checksum \
      "${SRC}:${remotepath}" \
      "${DST}:${dest}/"

  done < long_read_fofn.txt
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
#align_assemblies_parallelized
long_read_fofn
#download_with_globus
