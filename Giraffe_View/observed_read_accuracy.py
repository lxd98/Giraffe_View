import pysam
import re

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

def observed_accuracy(bam_file):
    """
    Calculates the observed accuracy, insertion, deletion, substitution, match,
    and identity rates for each read in a BAM file.
    """
    quality = {}
    with open("results/observed_quality/final_observed_accuracy.txt", "w") as outfile:  
        outfile.write("ID\tIns\tDel\tSub\tMat\tIden\tAcc\n")

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
            outfile.write(f"{read_ID}\t{Ins}\t{Del}\t{Sub}\t{Mat}\t{Iden:.4f}\t{Acc:.4f}\n")
            
    outfile.close()