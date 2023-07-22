# Test data

### estimate, observe, and gcbias

To test the functions of `estimate`, `observe`, and `gcbias`, we used the Nanopore R10.4.1 data from *S. aureus*.

```shell
# the code used to run the test datapip
giraffe observe --input Saur.fastq --ref Saur_genome.fasta --cpu 5 --plot
giraffe estimate --input Saur.fastq --cpu 5 --plot
giraffe gcbias --input results/observed_quality/tmp.sort.bam --ref Saur_genome.fasta --plot
```



### modibin

To test the function of **`modibin`**,  we used Nanopore R9.4.1 data from zebrafish kidney marrow. Here, We calculated the 5mC proportion at promoter level in chromosome 1.

```shell
# the code used to run the test datapip
giraffe modibin --input zebraFish_5mC.bed --ref zebraFish_promoter.csv --plot --cpu 5
```