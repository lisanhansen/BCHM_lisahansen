#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=PolII_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lisa.hansen@colorado.edu 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
pwd; hostname; date
echo "Lets go"
module load singularity/3.1.1
nextflow run nf-core/chipseq -r 1.2.1 \
-profile singularity \
--paired_end \
--input design.csv \
--fasta /scratch/Shares/rinnclass/CLASS_2022/data/GRCh38.p13.genome.fa \
--gtf /scratch/Shares/rinnclass/CLASS_2022/data/gencode.v32.annotation.gtf \
--macs_gsize 3.2e9 \
--blacklist /scratch/Shares/rinnclass/CLASS_2022/data/hg38-blacklist.v2.bed \
--email lisa.hansen@colorado.edu \
-resume \
-c nextflow.config
date
