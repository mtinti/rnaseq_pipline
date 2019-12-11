#!/bin/bash -l
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -o /dev/null
#$ -e /dev/null


grid_out='/homes/mtinti/RNAseq/viper-test/rna_seq/{experiment}'

exec >$grid_out'/multiqc_out.txt' 2>$grid_out'/multiqc_err.txt'

echo 'redirect grid outputs to:' $grid_out

conda activate ritSeq

echo 'run 15'
multiqc '{experiment}/data/' -o '{experiment}/multi_qc'


script_name='count_all_{experiment}.R'
echo 'running' $script_name
Rscript --vanilla $script_name



