giraffe estimate --input fastq.list --plot --cpu 4
giraffe observe --input fastq.list --plot --cpu 4 --ref Read/ecoli_chrom.fa
giraffe gcbias --input bam.list --plot --ref Read/ecoli_chrom.fa
giraffe modbin --input bed.list --cpu 4 --plot --pos Methylation/zf_promoter.db
