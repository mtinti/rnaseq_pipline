#!/bin/env Rscript


library(Rsubread)

count_table <- featureCounts(
files=c({bam_files}),
annot.ext="{gtf_file}",
isGTFAnnotationFile=TRUE,
isPairedEnd=TRUE,
requireBothEndsMapped=FALSE,
nthreads=8,
countMultiMappingReads=TRUE,
allowMultiOverlap=TRUE,
countChimericFragments=TRUE,
verbose=TRUE)


write.table(x=data.frame(count_table$annotation[,c("GeneID")], count_table$counts, stringsAsFactors=FALSE), 
file="{count_file}", quote=FALSE, sep="\t", row.names=FALSE)

save(count_table,file="{count_file}.Rda")
