#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=nextflow-job
#SBATCH --output=nextflow_%j.log
module load Nextflow
module load Singularity/3.10.2

nextflow run main.nf -c nextflow.config -ansi-log false -with-dag report_final.html 
