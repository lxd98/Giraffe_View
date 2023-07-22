import numpy as np
import math
import multiprocessing

def calculate_read_accuracy(string):
    """
    Calculate the estimated accuracy for a given string.
    Parameters:
    string (str): A string of bases in ASCII format.
    Returns:
    list: A list containing the estimated accuracy, error proportion, and Q-value.
    """
    error_list = []
    for base_value in string:
        ascii_value = ord(base_value) - 33
        error_proporation = math.pow(10, (-1) * int(ascii_value) / 10)
        error_list.append(error_proporation)
    error_mean = np.mean(error_list)
    return [1 - error_mean, error_mean, (-10) * math.log10(error_mean)]

def process_chunk(chunk):
    """
    Process a chunk of reads.
    Parameters:
    chunk (list): A list of reads.
    Returns:
    list: A list containing the estimated accuracy, error proportion, and Q-value for each read in the chunk.
    """
    results = []
    for line in chunk:
        read_id, quality = line
        quality = calculate_read_accuracy(quality)
        results.append([read_id, quality[0], quality[1], quality[2]])
    return results

def calculate_estimated_accuracy(input_file, num_processes, chunk_size=1000):
    """
    Calculate the estimated accuracy for each read in an input file and write the results to an output file.
    Parameters:
    input_file (str): The path to the input file.
    output_file (str): The path to the output file.
    num_processes (int): The number of processes to use for multiprocessing.
    chunk_size (int): The size of the chunks to split the input file into.
    """
    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    with open(input_file, "r") as input_file:
        count = 1
        chunk = []
        for line in input_file:
            # Get the read ID
            if count % 4 == 1:
                line = line.replace("\n", "")
                line = line.split(" ")
                read_id = line[0]
                count += 1
            elif count % 4 == 0:
                line = line.replace("\n", "")
                chunk.append((read_id, line))
                count += 1
            else:
                count += 1
            if len(chunk) == chunk_size:
                results.append(pool.apply_async(process_chunk, (chunk,)))
                chunk = []
        if len(chunk) > 0:
            results.append(pool.apply_async(process_chunk, (chunk,)))
    pool.close()
    pool.join()

    output_file = open("results/estimated_quality/final_estimated_accuracy.txt", "w")
    # with open("result/estimated/estimated_accuracy.txt", "w") as output_file:
    output_file.write("ID\tacc\terror\tQ_value\n")
    for result in results:
        for line in result.get():
            message = f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}"
            output_file.write(message + "\n")