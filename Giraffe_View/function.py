import subprocess
import re
import os
from termcolor import colored
from subprocess import Popen, PIPE
import pandas as pd

def print_with_color(input_string):
    print(colored(input_string, "green"))

def error_with_color(input_string):
    print(colored(input_string, "red"))

def loading_dataset(input_file):
    dataset = {}
    with open(input_file) as ff:
        for l in ff:
            l = l.replace("\n", "")
            l = l.split()

            # check the file
            if os.path.exists(l[2]):
                dataset[l[0]] = {}
                dataset[l[0]]["type"] = l[1]
                dataset[l[0]]["path"] = l[2]

            else:
                error_with_color("Please check the path of " + str(l[0]) + "!")

    return dataset

def cmd_shell(cammands, string):
    process = Popen(cammands.split(' '), stdout=subprocess.DEVNULL, universal_newlines=True)
    process.wait()
    err = process.communicate()

    if process.returncode == 0:
        # print('{} SUCCESS'.format(string))
        pass
    else:
        # print('{} FAILED'.format(string))
        error_with_color(err)

def mkdir_d(input_name):
    mes = "Giraffe_Results/" + str(input_name)
    cmd = ["mkdir", "-p", str(mes)]
    subprocess.run(cmd, check=True)

def count_indel_and_snv(str):
    dict = {}
    for i in str:
        dict[i] = dict.get(i, 0) + 1
    return dict


def bam2fastq(input_bam, CPU):
    with open("bam2fq.sh", "w") as ff:
        ff.write("samtools fastq " + str(input_bam) + " -@ " + str(CPU) + " > giraffe_tmp.fastq")
    ff.close()

#remove the insertion (I) in the tail of string
def remove_I(string):
    while string[-1] == "I":
        string = string[:-1]
    return(string)

# remove soft (S) and hard (H) clip in CIGAR and return the matched pairs
def remove_clip_list(input_cigar, input_pairs, input_ID):
    remove_cigarstring = re.findall(r"\d+[S, H]+", input_cigar)
    #HH & 0H & H0 & 00
    if ((len(remove_cigarstring) == 2) and (remove_cigarstring[0][-1] == remove_cigarstring[1][-1] == "H")) or ((len(remove_cigarstring) == 1) and (remove_cigarstring[-1] == "H")) or (len(remove_cigarstring) == 0):
        valid_pairs = input_pairs
    #SS
    elif (len(remove_cigarstring) == 2) and (remove_cigarstring[0][-1] == remove_cigarstring[1][-1] == "S"):
        remove_start_site = int(remove_cigarstring[0][:-1])
        tmp_pairs = input_pairs[remove_start_site:]
        remove_end_site = int(remove_cigarstring[1][:-1])
        valid_pairs = tmp_pairs[:len(tmp_pairs)-remove_end_site]
    # 0S & HS
    elif ((len(remove_cigarstring) == 1) and (input_cigar[-1] == "S")) or (len(remove_cigarstring) == 2) and (remove_cigarstring[0][-1] == "H") and ((remove_cigarstring[1][-1] == "S")):
        remove_end_site = int(remove_cigarstring[-1][:-1])
        valid_pairs = input_pairs[:len(input_pairs)-remove_end_site]
    # S0 & SH
    elif (len(remove_cigarstring) == 1) and (input_cigar[-1] != "S") or ((len(remove_cigarstring) == 2) and (remove_cigarstring[0][-1] == "S") and (remove_cigarstring[1][-1] == "H")):
        remove_start_site = int(remove_cigarstring[0][:-1])
        valid_pairs = input_pairs[remove_start_site:]
    else:
        print(str(input_ID) + ", please recheck this CIGAR and MD!")
    return(valid_pairs)

"""
only for base A T G C
(read_position, ref_position, "ref_base")
none    √   √   Deletion(D)
√   none    none Insertion(I)
√   √   N(A,T,G,C)  Match(M)
√   √   n(a,t,g,c)  Substitution(S)
"""
def get_base_alignment(input_list): 
    map_list = ["A", "T", "G", "C"]
    result = ""
    if input_list[0] == None:
        result = "D" # D = deletion
    else:
        if input_list[1] == None:
            result = "I" # I = insertion
        else:
            if input_list[2] in map_list:
                result = "M" # M = match
            else:
                result =  "S" # S = substitution
    return result

def process_in_chunks(file_path, chunk_size=10000):
    chunks = pd.read_csv(file_path, chunksize=chunk_size, sep="\t")
    results = []
    for chunk in chunks:
        results.append(chunk)
    return pd.concat(results)