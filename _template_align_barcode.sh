#!/bin/bash -l
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -o /dev/null
#$ -e /dev/null


grid_out='/homes/mtinti/RNAseq/viper-test/rna_seq/{experiment}'

exec >$grid_out'/{base_fastq}.out.txt' 2>$grid_out'/{base_fastq}.err.txt'

echo 'redirect grid outputs to:' $grid_out

conda activate ritSeq

echo 'conda env loaded'

#for the future, if always one fastqfile * experiment
#g_version=$1
#experiment=$2
#base_fastq=$3
#echo $g_version $experiment $base_fastq


echo 'load variable'
##variables
path_genome_index='genomes/{g_version}/{g_version}'
path_fastq_files='{experiment}/data/'
base_fastq='{base_fastq}'
path_out=$path_fastq_files$base_fastq'/'
gtf_file='genomes/{g_version}/{g_version}.gtf'
gff_file='genomes/{g_version}/{g_version}.gff'

echo 'make out dirs'
mkdir $path_out
mkdir $path_out'fastqc'
mkdir $path_out'fastp'
##qc of fastq 
echo 'run 0.1' 

inFileLength=$(echo $(zcat $path_fastq_files$base_fastq'1.fastq.gz'|wc -l)/4|bc)
echo 'inFileLength all: '$inFileLength

echo 'run 0.1.1 filter problematic reads' 

fastq_1=$path_fastq_files$base_fastq'1.fastq.gz'
fastq_2=$path_fastq_files$base_fastq'2.fastq.gz'
out_1=$path_fastq_files$base_fastq'f1.fastq.gz'
out_2=$path_fastq_files$base_fastq'f2.fastq.gz'
fastp -i $fastq_1 -I $fastq_2 -o $out_1 -O $out_2 -h $path_out'fastp/fastp_'$base_fastq'.html' -j $path_out'fastp/'$base_fastq'_fastp.json'


echo 'run 0.2 parse fastq for barcodes' 
python remove_barcode.py '{experiment}/data/{base_fastq}' $inFileLength

fastq_1=$path_fastq_files$base_fastq'1_barcode.fastq.gz'
fastq_2=$path_fastq_files$base_fastq'2_barcode.fastq.gz'

inFileLength=$(echo $(zcat $path_fastq_files$base_fastq'1_barcode.fastq.gz'|wc -l)/4|bc)
echo 'inFileLength barcode: '$inFileLength


echo 'run fastqc' 
fastqc  $fastq_1 -o $path_out'fastqc'
fastqc  $fastq_2 -o $path_out'fastqc'

echo 'run 2.0 align all reads'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 -S $path_out$base_fastq'.sam'

echo 'run 2.1'
samtools view -bS $path_out$base_fastq'.sam' > $path_out$base_fastq'.bam'
rm $path_out$base_fastq'.sam' ##saving space

echo 'run 4'
samtools sort $path_out$base_fastq'.bam' -o $path_out$base_fastq'.sorted.bam'

echo 'run 4'
samtools sort -n $path_out$base_fastq'.bam' -o $path_out$base_fastq'.sorted_name.bam'

echo 'run 5'
samtools index $path_out$base_fastq'.sorted.bam'

echo 'run 6'
rm $path_out$base_fastq'.bam'

echo 'run 7'
bedtools genomecov -ibam $path_out$base_fastq'.sorted.bam' -bg > $path_out$base_fastq'_coverage_bg.bed'
echo 'run 7.0'
bedtools genomecov -ibam $path_out$base_fastq'.sorted.bam' -d > $path_out$base_fastq'_coverage_d.bed'
echo 'run 7.1'
python plot_chr_coverage.py $path_out$base_fastq'_coverage_d.bed' '{experiment}/'$base_fastq'coverage_d.png'

#echo 'run 7.1'
#python mod_bed_graph.py $path_out$base_fastq'_coverage_bg.bed'
#echo 'run 7.2'
#bedmap --mean $gff_file'.bed' $path_out$base_fastq'_coverage_bg.bed' > $path_out$base_fastq'_coverage_reference_mean.bed'

unset DISPLAY

echo 'run 8'
mkdir $path_out'qualimap_results_rnaseq'
qualimap rnaseq -bam $path_out$base_fastq'.sorted.bam'  -gtf $gtf_file -outdir $path_out'qualimap_results_rnaseq' -pe

echo 'run 8'
mkdir $path_out'qualimap_results_bamqc'
qualimap bamqc  -bam $path_out$base_fastq'.sorted.bam' -outdir $path_out'qualimap_results_bamqc'

echo 'run 9'
mkdir $path_out'unmap/'
samtools view -u -f 12 -F 256 $path_out$base_fastq'.sorted.bam' > $path_out'unmap/'$base_fastq'.sorted_unmap.bam'

echo 'run 10'
samtools sort -n $path_out'unmap/'$base_fastq'.sorted_unmap.bam' -o $path_out'unmap/'$base_fastq'.sorted_unmap.bam'

echo 'run 11'
bamToFastq -i $path_out'unmap/'$base_fastq'.sorted_unmap.bam' -fq $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq' -fq2 $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq'


echo 'run 12'

script_name='count_'$base_fastq'.R'
echo 'running' $script_name
Rscript --vanilla $script_name

echo 'run 13'
mv $script_name {experiment}

echo 'run 14'
mkdir '{experiment}/data/{base_fastq}/fastq_screen'
fastq_screen --nohits $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq' $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq' \
--conf 'FastQ_Screen_Genomes/fastq_screen.conf' --outdir '{experiment}/data/{base_fastq}/fastq_screen'

#rm $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq'
#rm $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq'






