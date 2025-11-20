#!/bin/bash --login

#SBATCH --job-name=chaz-track-parser
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96GB
#SBATCH --time=0-48:00:00
#SBATCH --partition=Medium
#SBATCH --mail-type=END,FAIL
#SBATCH --output=jobs/%x/%j/stdout.txt

OUTPUT_DATA_DIR="data"

LOG_DIR=jobs/${SLURM_JOB_NAME}/${SLURM_JOB_ID}/
mkdir -p $LOG_DIR

log-proc-mem $$ > $LOG_DIR/proc-mem.csv &
LOG_PROC_PID=$!

micromamba run --name chaz snakemake --unlock
micromamba run --name chaz \
    snakemake \
        --resources mem_mb=96000 \
        --cores 16 \
        --keep-going \
        --rerun-incomplete \
        $OUTPUT_DATA_DIR/out/CHAZ-normalised-freq

kill $LOG_PROC_PID
plot-mem $LOG_DIR/proc-mem.csv $LOG_DIR/proc-mem.png
