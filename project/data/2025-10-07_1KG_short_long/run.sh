#!/bin/bash
set -e -x
set -euo pipefail
function get_sample_list {
	wget -o sample_names_shortR_longR.txt https://raw.githubusercontent.com/nefte48/Imputation_data/refs/heads/main/results/2025_08_27_1kg_IG_snps/sample_names_shortR_longR.txt
}


function short_read_fofn {
    index1=1000G_2504_high_coverage.sequence.index
    index2=1000G_698_related_high_coverage.sequence.index
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index -O "${index1}"
    wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index -O "${index2}"
    cat sample_names_shortR_longR.txt | while read sample
    do
	    fn=`awk -F'\t' -v s=${sample} '$10 == s {print $0}' ${index1} ${index2} | cut -f1`
	    echo -e "${sample}\t${fn}"
    done > short_read_fofn.txt

}


function long_read_fofn {
    # get file locally, transferred to minerva using scp
    # scp igsr_HGSVC3.tsv odrio10@minerva.hpc.mssm.edu:/
    
    


    index=igsr_HGSVC3.tsv.tsv
    cat sample_names_shortR_longR.txt | while read sample    
    do
	    fn=`awk -F'\t' -v s=${sample} '$6 == s {print $1}' ${index} | cut -f1`
	    echo -e "${sample}\t${fn}"
    done > long_read_fofn.txt
}

function download_short {
    echo -e "\n\033[1;33mRun the following commands:\033[0m\n
--------------------------------------------------
\033[1;36mscreen -S short-reads\033[0m

\033[1;36mbsub -P acc_oscarlr -q interactive -n 6 -W 12:00 -R \"span[hosts=1]\" -Is /bin/bash\033[0m

\033[1;36mmodule load aspera-connect\033[0m

\033[1;36mcut -f2 REPLACE_THIS_WITH_THE_DOWNLOAD_LINKS_FILE_NAME_EITHER_1_OR_2 | sed -E 's|^[^/]+://[^/]+||' | \\
  xargs -n1 -I{} /hpc/packages/minerva-centos7/aspera-connect/4.2.4/aspera/bin/ascp \\
  -i /hpc/packages/minerva-centos7/aspera-connect/3.9.6/etc/asperaweb_id_dsa.openssh \\
  -P33001 -O33001 -v -L- \\
  'era-fasp@fasp.sra.ebi.ac.uk:{}' \\
  '/sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long'\033[0m

\033[1;33mTo check if it's downloading:\033[0m
\033[1;36mls /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long\033[0m
--------------------------------------------------\n"

}


function download_short_wget {
	cut -f2 short_read_fofn.txt | xargs -n1 -P1  wget -c -P /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long {}
	
}

function download_long_wget {
	awk '{print $NF}' long_read_fofn.txt | xargs -n1 -P1  wget -c -P /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long {}
}

function get_franken_reference {
	wget http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta -O /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/reference.fasta
	
	wget https://raw.githubusercontent.com/Watson-IG/immune_receptor_genomics/refs/heads/main/240520/reference.fasta.fai -O /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/reference.fasta.fai
}

function create_franken_minimap2_index_file {
	cd /sc/arion/scratch/hiciaf01/projects/imputation/data/2025-10-07_1KG_short_long/
	module load minimap2
	minimap2 -d franken_reference_minimap2_index.mmi reference.fasta
	
	cd -
}




function align_long_reads_minimap {
	
	
	vim run.sh
}
#short_read_fofn
#download_short_wget 
#long_read_fofn
#download_long_wget
#get_franken_reference
create_franken_minimap2_index_file
