import pandas as pd
import os
from tqdm import tqdm
from multiprocessing import Process, Queue, current_process

def bin_methylation(input_str, data):
    """Calculate average methylation level for a target region."""
    input_str = input_str.split(",")
    target_id = input_str[3]
    chrom_data = data[data["ref_ID"] == input_str[0]]
    target_data = chrom_data[(int(input_str[1]) <= chrom_data["start"]) & (chrom_data["end"] <= int(input_str[2]))]
    mean_methylation = target_data["modified_percentage"].mean()
    return (target_id, mean_methylation)


def worker(inqueue, outqueue, data):
    """Process input data and put output data into output queue."""
    for frame in iter(inqueue.get, "STOP"):
        outqueue.put(bin_methylation(frame, data))
    outqueue.put(f"{current_process().name}: BYE!")


def manager(process_num, input_file, input_ref):
    """Distribute work to worker processes and write output to file."""
    data = pd.read_csv(input_file, header=None, names=["ref_ID", "start", "end", "modified_percentage"], sep="\t")
    targets = pd.read_csv(input_ref, header=None, names=["chr", "start", "end", "target_id"])
    total = len(targets)
    inqueue = Queue()
    outqueue = Queue()

    for i in range(process_num):
        Process(target=worker, args=(inqueue, outqueue, data)).start()

    with open(input_ref) as ff:
        for l in ff:
            inqueue.put(l.strip())

    for i in range(process_num):
        inqueue.put("STOP")

    stop_count = 0
    with tqdm(total=total, desc=input_file) as pbar, open(f"results/regional_modification/{input_file.replace('.bed', '')}_{input_ref.replace('.csv', '')}.bed", "w") as of:
        while stop_count < process_num:
            pbar.update()
            result = outqueue.get()
            if len(result) == 2:
                of.write(f"{result[0]}\t{result[1]}\n")
            elif result[-4:] == "BYE!":
                stop_count += 1

def methylation_calculation(input_file, input_ref, process_num):
    """Calculate average methylation level for target regions in input_ref using methylation data in input_file."""
    manager(process_num, input_file, input_ref)