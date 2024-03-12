wget https://figshare.com/ndownloader/files/44967445 -O fastq.list
wget https://figshare.com/ndownloader/files/44967442 -O bed.list
wget https://figshare.com/ndownloader/files/44967499 -O bam.list

# The reference and ONT reads (R10.4.1 and R9.4.1) of E.coli
wget https://figshare.com/ndownloader/files/44967436 -O Read.tar.gz

# The 5mC methylation files of zebrafish blood and kidney samples.
# The bed file is the gene promoter position in chromosome 1. 
wget https://figshare.com/ndownloader/files/44967427 -O Methylation.tar.gz

tar -xzvf Read.tar.gz
tar -xzvf Methylation.tar.gz
rm Read.tar.gz Methylation.tar.gz