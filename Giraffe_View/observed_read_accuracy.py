import pysam
import re
from os import system 
from Giraffe_View.function import *

def data_process(sample_ID, data_type, data_path, ref, threads=10):
    # Define the commands as a list of strings to avoid issues with spaces
    # in file names or command arguments
    # path =  os.getcwd()

    # cmd0 = ["mkdir", "-p", "results/observed_quality"]
    # cmd1 = ["seqkit", "seq", read, "-m", "200", "-Q", "7", "-g", "-j", str(threads), "-o", "results/observed_quality/clean.fastq"]
    output = "Giraffe_Results/2_Observed_quality/" + str(sample_ID) + ".bam"
    
    if data_type == "ONT":
        cmd1 = ["minimap2", "-ax", "map-ont", "-o", "Giraffe_Results/2_Observed_quality/tmp.sam", "--MD", \
        "--secondary=no", "-L", "-t", str(threads), ref, data_path]

    elif data_type == "Pacbio":
        cmd1 = ["minimap2", "-ax", "map-hifi", "-o", "Giraffe_Results/2_Observed_quality/tmp.sam", "--MD", \
        "--secondary=no", "-L", "-t", str(threads), ref, data_path]

    cmd2 = ["samtools", "view", "-bS", "-F4", "-@", str(threads), "-o", "Giraffe_Results/2_Observed_quality/tmp.bam", "Giraffe_Results/2_Observed_quality/tmp.sam"]
    cmd3 = ["samtools", "sort", "-@", str(threads), "-o", output, "Giraffe_Results/2_Observed_quality/tmp.bam"]
    cmd4 = ["samtools", "index", "-@", str(threads), output]
    cmd5 = ["rm", "-rf", "Giraffe_Results/2_Observed_quality/tmp.bam", "Giraffe_Results/2_Observed_quality/tmp.sam"]

    # Run each command and check the return code
    for i, cmd in enumerate([cmd1, cmd2, cmd3, cmd4, cmd5]):
        try:
            subprocess.run(cmd, check=True)
            # print("Command {} succeeded".format(i + 1))
        except subprocess.CalledProcessError as e:
            print("Command {} failed with error code {}".format(i + 1, e.returncode))
            print(e.output)
            # Raise an exception to indicate that processing failed
            raise Exception("Data processing failed")

def identify_match(cigar):
    """
    Identifies the number of matching bases in a read from its CIGAR string.
    """
    cigar_mat = re.findall(r"\d+M", cigar)
    base_num_mat = sum(int(i[:-1]) for i in cigar_mat)
    return base_num_mat

def identify_insertion(cigar):
    """
    Identifies the number of inserted bases in a read from its CIGAR string.
    """
    cigar_ins = re.findall(r"\d+I", cigar)
    base_num_ins = sum(int(i[:-1]) for i in cigar_ins)
    return base_num_ins

def identify_deletion(cigar):
    """
    Identifies the number of deleted bases in a read from its CIGAR string.
    """
    cigar_del = re.findall(r"\d+D", cigar)
    base_num_del = sum(int(i[:-1]) for i in cigar_del)
    return base_num_del

def identify_substitution(md):
    """
    Identifies the number of substitutions in a read from its MD tag.
    """
    return len(re.findall(r"\d+[ATCG]", md))

def observed_accuracy(bam_file, sample_ID):
    """
    Calculates the observed accuracy, insertion, deletion, substitution, match,
    and identity rates for each read in a BAM file.
    """
    quality = {}
    with open("Giraffe_Results/2_Observed_quality/"+sample_ID+".acc_tmp", "w") as outfile:  
        bamfile= pysam.AlignmentFile(bam_file, 'rb')
        for read in bamfile:
                # Skip reads with unmapped or secondary alignments
                # if int(read.flag) not in [0, 16]:
                #     continue
            read_ID = read.query_name
            read_cigar = read.cigarstring
            read_md = read.get_tag("MD")

            Ins = identify_insertion(read_cigar)
            Del = identify_deletion(read_cigar)
            Sub = identify_substitution(read_md)
            Mat = identify_match(read_cigar) - Sub

            if read_ID not in quality.keys():
                quality[str(read_ID)] = {}
                quality[str(read_ID)]["Ins"] = int(Ins)
                quality[str(read_ID)]["Del"] = int(Del)
                quality[str(read_ID)]["Sub"] = int(Sub)
                quality[str(read_ID)]["Mat"] = int(Mat)

            else:
                quality[str(read_ID)]["Ins"] += int(Ins)
                quality[str(read_ID)]["Del"] += int(Del)
                quality[str(read_ID)]["Sub"] += int(Sub)
                quality[str(read_ID)]["Mat"] += int(Mat)

        bamfile.close()

        for read_key in quality.keys():
            read_ID = read_key
            Ins = quality[str(read_key)]["Ins"]
            Del = quality[str(read_key)]["Del"]
            Sub = quality[str(read_key)]["Sub"]
            Mat = quality[str(read_key)]["Mat"]

            total = Ins + Del + Sub + Mat
            Acc = Mat / total if total > 0 else 0
            Iden = Mat / (Mat + Sub) if (Mat + Sub) > 0 else 0
            outfile.write(f"{read_ID}\t{Ins}\t{Del}\t{Sub}\t{Mat}\t{Iden:.4f}\t{Acc:.4f}\t{sample_ID}\n")
    outfile.close()

def merge_results_observed_acc():
    with open("Giraffe_Results/2_Observed_quality/header", "a") as ff:
        ff.write("ID\tIns\tDel\tSub\tMat\tIden\tAcc\tGroup\n")
    ff.close()

    system("cat Giraffe_Results/2_Observed_quality/header \
        Giraffe_Results/2_Observed_quality/*.acc_tmp > \
        Giraffe_Results/2_Observed_quality/Observed_information.txt")
    system("rm Giraffe_Results/2_Observed_quality/*.acc_tmp \
        Giraffe_Results/2_Observed_quality/header")