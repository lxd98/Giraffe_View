import subprocess
import re
import os
from termcolor import colored
from subprocess import Popen, PIPE

def print_with_color(input_string):
    print(colored(input_string, "green"))

def error_with_color(input_string):
    print(colored(input_string, "red"))

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

def Data_process(read, ref, threads=10):
    # Define the commands as a list of strings to avoid issues with spaces
    # in file names or command arguments
    # path =  os.getcwd()
    cmd0 = ["mkdir", "-p", "results/observed_quality"]
    cmd1 = ["seqkit", "seq", read, "-m", "200", "-Q", "7", "-g", "-j", str(threads), "-o", "results/observed_quality/clean.fastq"]
    cmd2 = ["minimap2", "-ax", "map-ont", "-o", "results/observed_quality/tmp.sam", "--MD", "--secondary=no", "-L", "-t", str(threads), ref, "results/observed_quality/clean.fastq"]
    cmd3 = ["samtools", "view", "-bS", "-F4", "-@", str(threads), "-o", "results/observed_quality/tmp.bam", "results/observed_quality/tmp.sam"]
    cmd4 = ["samtools", "sort", "-@", str(threads), "-o", "results/observed_quality/tmp.sort.bam", "results/observed_quality/tmp.bam"]
    cmd5 = ["samtools", "index", "-@", str(threads), "results/observed_quality/tmp.sort.bam"]
    cmd6 = ["rm", "-rf", "results/observed_quality/tmp.sam", "results/observed_quality/tmp.bam"]

    # Run each command and check the return code
    for i, cmd in enumerate([cmd0, cmd1, cmd2, cmd3, cmd4, cmd5, cmd6]):
        try:
            subprocess.run(cmd, check=True)
            # print("Command {} succeeded".format(i + 1))
        except subprocess.CalledProcessError as e:
            print("Command {} failed with error code {}".format(i + 1, e.returncode))
            print(e.output)
            # Raise an exception to indicate that processing failed
            raise Exception("Data processing failed")

def mkdir_d(input_name):
    mes = "results/" + str(input_name)
    cmd = ["mkdir", "-p", str(mes)]
    subprocess.run(cmd, check=True)


def count_indel_and_snv(str):
    dict = {}
    for i in str:
        dict[i] = dict.get(i, 0) + 1
    return dict

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