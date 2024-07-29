import pandas as pd
import multiprocessing

def run_regional_methylation(input_methyl, input_target, sample_ID, num_processes):
    # Read methylation and target data from the input files
    methyl = pd.read_csv(input_methyl, sep='\t', header=None, names=["CHROM", "START", "END", "VALUE"])
    target = pd.read_csv(input_target, sep='\t', header=None, names=["CHROM", "START", "END", "ID"])

    # Determine if VALUE is in range 0-1 or 0-100
    max_value = methyl['VALUE'].max()
    value_scale = 100 if max_value > 1 else 1

    # Get the unique chromosomes from the target data
    unique_chromosomes = set(target['CHROM'])

    # Create a multiprocessing pool with the specified number of processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        jobs = []
        for chromosome in unique_chromosomes:
            # Filter the methylation and target data for the current chromosome
            sub_methyl = methyl[methyl["CHROM"] == chromosome][["START", "END", "VALUE"]]
            sub_target = target[target["CHROM"] == chromosome][["START", "END", "ID"]]
            # Create an asynchronous job for processing the data
            jobs.append(pool.apply_async(regional_methylation_bed_worker, (sub_methyl, sub_target, str(sample_ID), chromosome, value_scale)))

        # Wait for all jobs to complete
        for job in jobs:
            job.get()

def regional_methylation_bed_worker(input_methyl, input_target, sample_ID, chromosome, value_scale):
    # Open the output file for writing
    with open(f"Giraffe_Results/4_Regional_modification/Temp_methy_{sample_ID}_{chromosome}.txt", "w") as ff:
        # Iterate over each row in the target data
        for row in input_target.itertuples(index=True, name='Pandas'):
            start = row.START
            end = row.END
            target_ID = row.ID

            # Filter the methylation data for the current region
            target_data = input_methyl[(start <= input_methyl["START"]) & (input_methyl["END"] <= end)]
            
            # Calculate the mean methylation value and scale it
            mean_methylation = target_data["VALUE"].mean() / value_scale

            # Write the result to the output file
            ff.write(f"{target_ID}\t{mean_methylation:.4f}\t{sample_ID}\n")
