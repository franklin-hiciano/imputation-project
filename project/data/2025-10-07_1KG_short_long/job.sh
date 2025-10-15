#!/bin/bash
#BSUB -P acc_oscarlr
#BSUB -W 03:00
#BSUB -q interactive
#BSUB -n 2
#BSUB -J short-reads-download
#BSUB -o short-reads.%J.out
#BSUB -e short-reads.%J.err
#BSUB -cwd .
#BSUB -R span[hosts=1]

set -euo pipefail

# Run your workload (no nohup needed in batch)
bash run.sh



