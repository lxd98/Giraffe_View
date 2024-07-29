from os import system 
import pandas as pd
import pysam
import re
from Giraffe_View.function import *
import multiprocessing

def homopolymer_summary_1(input_file, sample_ID, chromosome):
	data = {}
	with open(input_file) as ff:
		for line in ff:
			line = line.replace("\n", "")
			line = line.split("\t")
			position = line[0] + "_" + line[1] + "_" + line[2]

			if position not in data.keys():
				data[position] = {}
				data[position]["type"] = str(line[3]) + line[4]
				data[position]["depth"] = 1
				data[position]["mat"] = 0

				if int(line[3]) == int(line[5]):
					data[position]["mat"] += 1

			elif position in data.keys():
				data[position]["depth"] += 1
				if int(line[3]) == int(line[5]):
					data[position]["mat"] += 1
	ff.close()

	output_1 = f"Giraffe_Results/2_Observed_quality/{sample_ID}_homopolymer_in_reference_{chromosome}.txt"
	with open(output_1, "w") as ff:
		# ff.write("pos\tnum_of_mat\tdepth\ttype\tGroup\n")
		for i in data.keys():
			mes = str(i) + "\t" + str(data[i]["mat"]) + "\t" + str(data[i]["depth"]) + "\t" + data[i]["type"]
			ff.write(mes + "\t" + str(sample_ID) + "\n")
		ff.close()

def homopolymer_summary_2(sample_ID):
	input_file = "Giraffe_Results/2_Observed_quality/" + str(sample_ID) + ".homopolymer_in_reference.txt"
	output_file = "Giraffe_Results/2_Observed_quality/" + str(sample_ID) + ".homo_tmp"
	data = pd.read_table(input_file, sep="\t")	
	valid = data[data["depth"] >= 3].copy()

	if len(valid) != 0:
		ff = open(output_file, "w")
		valid["rate"] = valid["num_of_mat"] / valid["depth"]
		
		def Abase(x):
			if re.search(".*A", x): 
				return(True)
			else: 
				return(False)
		
		def Tbase(x):
			if re.search(".*T", x): 
				return(True)
			else: 
				return(False)
		
		def Gbase(x):
			if re.search(".*G", x): 
				return(True)
			else: 
				return(False)

		def Cbase(x):
			if re.search(".*C", x): 
				return(True)
			else: 
				return(False)
        
		T_acc = valid[valid["type"].apply(Tbase)]["rate"].mean()
		G_acc = valid[valid["type"].apply(Gbase)]["rate"].mean()
		C_acc = valid[valid["type"].apply(Cbase)]["rate"].mean()
		A_acc = valid[valid["type"].apply(Abase)]["rate"].mean()

		ff.write(f"A\t{A_acc:.4f}\t{sample_ID}\n")
		ff.write(f"T\t{T_acc:.4f}\t{sample_ID}\n")
		ff.write(f"C\t{C_acc:.4f}\t{sample_ID}\n")
		ff.write(f"G\t{G_acc:.4f}\t{sample_ID}\n")
		ff.close()
	else:
		error_with_color("The read coverage of data was too shallow to conduct the homopolymer analysis!!!") 

def merge_results_observed_homopolymer():
    with open("Giraffe_Results/2_Observed_quality/header", "a") as ff:
        ff.write("Base\tAccuracy\tGroup\n")
    ff.close()

    system("cat Giraffe_Results/2_Observed_quality/header \
        Giraffe_Results/2_Observed_quality/*.homo_tmp > \
        Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt")

    system("rm Giraffe_Results/2_Observed_quality/*homo_tmp \
        Giraffe_Results/2_Observed_quality/header"
        )

def homopolymer_from_bam_worker(input_bamfile, sample_ID, chromosome):
    bamfile = pysam.AlignmentFile(input_bamfile, "rb")
    output = f"Giraffe_Results/2_Observed_quality/{sample_ID}_homopolymer_detail_{chromosome}.txt"
    
    with open(output, "w") as ff:
        for read in bamfile.fetch(chromosome):
            read_ID = read.query_name
            read_pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            read_cigar = read.cigarstring
            read_ref_id = read.reference_name
            read_valid_pair = remove_clip_list(read_cigar, read_pair, read_ID)

            homopolymer_ref = ""
            homopolymer_read = ""
            homopolymer_ref_pos = []
            count = 1

            for base in read_valid_pair:
                base_alignment = get_base_alignment(base)
                if homopolymer_ref == "":
                    if base_alignment != "I":
                        homopolymer_ref = str(base[2]).upper()
                        homopolymer_read = str(base_alignment)
                        homopolymer_ref_pos.append(base[1])
                else:
                    if base[2] is None:
                        homopolymer_read += str(base_alignment)
                    else:
                        base_ref = str(base[2]).upper()
                        if base_ref == homopolymer_ref[0]:
                            homopolymer_ref += base_ref
                            homopolymer_read += str(base_alignment)
                            homopolymer_ref_pos.append(base[1])
                        else:
                            if len(homopolymer_ref) >= 4:
                                homopolymer_ref_pos = [
                                    str(len(homopolymer_ref)) + homopolymer_ref[0],
                                    str(remove_I(homopolymer_read)),
                                    str(read_ref_id),
                                    *homopolymer_ref_pos
                                ]
                                stat_info = count_indel_and_snv(homopolymer_ref_pos[1])
                                stat_info = {k: stat_info.get(k, 0) for k in ['M', 'D', 'S', 'I']}

                                mes = (
                                    f"{homopolymer_ref_pos[2]}\t{homopolymer_ref_pos[3]}\t{homopolymer_ref_pos[-1]}\t"
                                    f"{homopolymer_ref_pos[0][:-1]}\t{homopolymer_ref_pos[0][-1]}\t"
                                    f"{stat_info['M']}\t{stat_info['D']}\t{stat_info['I']}\t{stat_info['S']}\t"
                                    f"{read_ID}\t{sample_ID}")

                                ff.write(mes + "\n")
                                count += 1

                            homopolymer_ref = base_ref
                            homopolymer_read = str(base_alignment)
                            homopolymer_ref_pos = [base[1]]
    bamfile.close()
    ff.close()
    homopolymer_summary_1(output, sample_ID, chromosome)

def run_homopolymer_from_bam(input_bamfile, sample_ID, num_processes=10):
    bamfile = pysam.AlignmentFile(input_bamfile, "rb")
    chromosomes = bamfile.references

    with multiprocessing.Pool(processes=num_processes) as pool:
        jobs = []
        for chromosome in chromosomes:
            jobs.append(pool.apply_async(homopolymer_from_bam_worker, (input_bamfile, sample_ID, chromosome)))

        for job in jobs:
            job.get()
