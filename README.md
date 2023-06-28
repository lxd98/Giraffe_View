**Giraffe_View** is designed to help assess and visualize the accuracy of a sequencing dataset, specifically for Oxford Nanopore Technologies (ONT) long-read sequencing.



**Workflow:**

```mermaid
graph TD
raw_data --> |Quality control| clean_data
raw_data --> |Basecall| modification_file
modification_file --> modification_distribution
clean_data --> Estimated_accuracy
clean_data --> |Reference| aligned_file
aligned_file --> Homopolymer_analysis
aligned_file --> GC_bias 
aligned_file --> Observed_accuracy
Observed_accuracy --> Read_level
Observed_accuracy --> Contigs/Chromsomes_level
```

The update is available [here](https://github.com/lrslab/Giraffe_View).
