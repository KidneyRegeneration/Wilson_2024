#!/bin/bash

#PBS -N citeseqhash
#PBS -q batch
#PBS -l nodes=1:ppn=16
#PBS -l mem=160GB
#PBS -l walltime=100:00:00
#PBS -A kidn1
#PBS -m abe
#PBS -M sean.wilson@mcri.edu.au

cd $PBS_O_WORKDIR

module load R python

source /group/kidn1/hpc/cite-seq-count/bin/activate

for c in SW_CPT{1..2}; do
mkdir $c/cat
cat $c/*R1*.fastq.gz > $c/cat/all_R1.fastq.gz
cat $c/*R2*.fastq.gz > $c/cat/all_R2.fastq.gz

CITE-seq-Count \
-R1 $c/cat/all_R1.fastq.gz \
-R2 $c/cat/all_R2.fastq.gz \
-t $c/tags.csv \
-cbf 1 \
-cbl 16 \
-umif 17 \
-umil 28 \
-cells 20000 \
-o HASHED/$c \
--threads 8 \
--dense 

done
