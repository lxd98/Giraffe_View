from os import system 
import numpy as np
import math
import multiprocessing

def GC_content(string):
    read = str(string).upper()
    length = len(read)
    c = read.count("C")
    g = read.count("G")
    GC = (c+g)/length
    return[length, GC]

def Qvalue_to_accuracy(string):
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
        read_id, sequence, quality = line
        GC = GC_content(sequence)
        quality = Qvalue_to_accuracy(quality)
        results.append([read_id, quality[0], quality[1], quality[2], GC[0], GC[1]])
    return results

def calculate_estimated_accuracy(input_type, input_file, num_processes, chunk_size=1000):
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
            
            # get the read ID
            if count % 4 == 1:
                line = line.replace("\n", "")
                line = line.split(" ")
                read_id = line[0]
                count += 1
            
            # get the read sequence
            elif count % 4 == 2:
                sequence = line.replace("\n", "")
                count += 1

            # get the quality string
            elif count % 4 == 0:
                line = line.replace("\n", "")
                chunk.append((read_id, sequence, line))
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
    input_file.close()

    file = "Giraffe_Results/1_Estimated_quality/" + str(input_type) + ".tmp"
    with open(file, "a") as output_file:
            for result in results:
                for line in result.get():
                    message = f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}"
                    message += f"\t{line[4]}\t{line[5]}\t{input_type}"
                    output_file.write(message + "\n")
    output_file.close()
    output_file.close()


def merge_results():
    with open("Giraffe_Results/1_Estimated_quality/header", "a") as ff:
        ff.write("ReadID\tAccuracy\tError\tQ_value\tLength\tGC_content\tGroup\n")
    ff.close()
    system("cat Giraffe_Results/1_Estimated_quality/header \
        Giraffe_Results/1_Estimated_quality/*tmp > \
        Giraffe_Results/1_Estimated_quality/Estimated_information.txt")
    system("rm Giraffe_Results/1_Estimated_quality/*tmp Giraffe_Results/1_Estimated_quality/header")