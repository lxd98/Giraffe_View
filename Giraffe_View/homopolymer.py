import pandas as pd
import pysam
import re
import multiprocessing
from Giraffe_View.function import *

def homopolymer_from_bam(input_bamfile):
	bamfile = pysam.AlignmentFile(input_bamfile, threads=4)
	dir_polymer = {}

	for read in bamfile:
		read_ID = read.query_name
		read_pair = read.get_aligned_pairs(matches_only=False, with_seq=True) #show (read_position, ref_position, "ref_base")
		read_cigar = read.cigarstring # get the CIGAR for each read
		read_ref_id = read.reference_name # the aligned reference chromsome id e.g chr1

		read_valid_pair = remove_clip_list(read_cigar, read_pair, read_ID) # remove the cliping, only the mapped bases were remained.
		
		homoploymer_ref = "" 
		homoploymer_read = "" 
		homoploymer_ref_pos = [] 
		dir_polymer[read_ID] = {} 
		count = 1 

		for base in read_valid_pair:
			#start
			if homoploymer_ref == "":
				if get_base_alignment(base) != "I":
					homoploymer_ref = str(base[2]).upper()
					homoploymer_read = str(get_base_alignment(base))
					homoploymer_ref_pos.append(base[1])
			else:
				if base[2] == None:
					homoploymer_read += str(get_base_alignment(base))
				else:
					if str(base[2]).upper() == homoploymer_ref[0]:
						homoploymer_ref += str(base[2]).upper()
						homoploymer_read += str(get_base_alignment(base))
						homoploymer_ref_pos.append(base[1])
					elif str(base[2]).upper() != homoploymer_ref[0]:
						if len(homoploymer_ref) >= 4:
							homoploymer_ref_pos.insert(0, str(read_ref_id))
							homoploymer_ref_pos.insert(0, str(remove_I(homoploymer_read)))
							homoploymer_ref_pos.insert(0, str(len(homoploymer_ref)) + homoploymer_ref[0])
							dir_polymer[read_ID][str(count)] = homoploymer_ref_pos.copy()
							count += 1
						del homoploymer_ref_pos[:]
						homoploymer_ref = str(base[2]).upper()
						homoploymer_read = str(get_base_alignment(base))
						homoploymer_ref_pos.append(base[1])
	bamfile.close()

	ff = open("results/observed_quality/homo.txt", "w")
	for read in dir_polymer.keys():
		for n in dir_polymer[read].keys():
			stat_info = count_indel_and_snv(str(dir_polymer[read][n][1]))

			if "M" not in stat_info.keys():
				stat_info["M"] = 0
			if "D" not in stat_info.keys():
				stat_info["D"] = 0
			if "S" not in stat_info.keys():
				stat_info["S"] = 0
			if "I" not in stat_info.keys():
				stat_info["I"] = 0

			mes = str(dir_polymer[read][n][2]) + "\t" +  str(dir_polymer[read][n][3]) + "\t" + str(dir_polymer[read][n][-1]) + "\t" # chr start end
			mes += str(dir_polymer[read][n][0][:-1]) + "\t" + str(dir_polymer[read][n][0][-1]) + "\t"  #  4 A
			mes += str(stat_info["M"]) + "\t" + str(stat_info["D"]) + "\t" + str(stat_info["I"]) + "\t" + str(stat_info["S"])  + "\t"  # map del ins sub
			mes += str(read)
			ff.write(mes + "\n")
	ff.close()

def homopolymer_summary_1(input_file):
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

	ff = open("results/observed_quality/homo_tmp.txt", "w")
	ff.write("pos\tnum_of_mat\tdepth\ttype\n")   

	for i in data.keys():
		mes = str(i) + "\t" + str(data[i]["mat"]) + "\t" + str(data[i]["depth"]) + "\t" + data[i]["type"]
		ff.write(mes + "\n")
	ff.close()

def homopolymer_summary_2():
	data = pd.read_table("results/observed_quality/homo_tmp.txt", sep="\t")	
	valid = data[data["depth"] >= 3].copy()

	if len(valid) != 0:
		ff = open("results/observed_quality/final_homo_summary.txt", "w")
		ff.write("base\tvalue\n")

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
		
		ff.write("T\t" + str(valid[valid["type"].apply(Tbase)]["rate"].mean()) + "\n")
		ff.write("G\t" + str(valid[valid["type"].apply(Gbase)]["rate"].mean()) + "\n")
		ff.write("A\t" + str(valid[valid["type"].apply(Abase)]["rate"].mean()) + "\n")
		ff.write("C\t" + str(valid[valid["type"].apply(Cbase)]["rate"].mean()) + "\n")
		ff.close()
	else:
		print("The read coverage of data was shallow, homopolymer analysis failed!!!")