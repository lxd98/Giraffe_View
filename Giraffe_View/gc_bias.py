import os
from os import system
import pandas as pd
from Giraffe_View.function import cmd_shell
import multiprocessing

def classify_by_chromosome(input_file):
    classified_lines = {}

    # Read the input file and classify the lines
    with open(input_file, 'r') as file:
        for line in file:
            # Split the line by tab
            fields = line.strip().split('\t')
            first_field = fields[0]

            # Add the line to the appropriate list in the dictionary
            if first_field not in classified_lines:
                classified_lines[first_field] = []
            classified_lines[first_field].append(line)

    # Write the classified lines to separate output files
    for key, lines in classified_lines.items():
        output_file = f"Giraffe_Results/3_GC_bias/{key}_gcbias_bin.bed"
        with open(output_file, 'w') as file:
            file.writelines(lines)

def get_bin_bed(input_reference, input_binsize):
	if not os.path.exists(f"{input_reference}.fai"):
		system(f"samtools faidx {input_reference}")
	system(f"bedtools makewindows -g {input_reference}.fai -w {input_binsize} > Giraffe_Results/3_GC_bias/bin.bed")
	classify_by_chromosome("Giraffe_Results/3_GC_bias/bin.bed")

def get_bin_GC(args):
	input_reference, input_chromosome, path = args
	system(f"bedtools nuc -fi {input_reference} -bed {path}/{input_chromosome}_gcbias_bin.bed > {path}/{input_chromosome}_bin_GC.tmp")

	input_file = f"{path}/{input_chromosome}_bin_GC.tmp"
	output = f"{path}/{input_chromosome}_bin_GC.txt"
	with open(input_file, "r") as ff:
		with open(output, "w") as of:
			for bin in ff:
				if bin[0] != "#":
					bin = bin.replace("\n", "")
					bin = bin.split()
					bin_chrom = bin[0]
					bin_start = bin[1]
					bin_end = bin[2]
					bin_gc = bin[4]
					mes = str(bin_chrom) + "\t" + str(bin_start) + "\t"
					mes += str(bin_end) + "\t" + str(bin_gc) + "\n"
					of.write(mes)
	system(f"rm {path}/{input_chromosome}_bin_GC.tmp")

def manager_GC_content(input_reference, num_cpus):
	processes = []
	chromosomes = []

	with open(f"{input_reference}.fai", "r") as ff:
		for l in ff.readlines():
			l = l.replace("\n","").split()
			chromosomes.append(l[0])
	ff.close()

	args = [(input_reference, chrom, "Giraffe_Results/3_GC_bias") for chrom in chromosomes]
	with multiprocessing.Pool(processes=num_cpus) as pool:
		pool.map(get_bin_GC, args)

	path = "Giraffe_Results/3_GC_bias"
	system(f"cat {path}/*_bin_GC.txt > {path}/bin_GC.txt")
	system(f"rm {path}/*_bin_GC.txt")

def get_bin_depth(args):
	input_sample_ID, input_bam, input_chromosome, path = args
	system(f"samtools bedcov {path}/{input_chromosome}_gcbias_bin.bed {input_bam} > {path}/{input_sample_ID}_{input_chromosome}_bin_depth.txt")

def manager_bin_depth(input_reference,sample_ID, bamfile, num_cpus):
	processes = []
	chromosomes = []

	with open(f"{input_reference}.fai", "r") as ff:
		for l in ff.readlines():
			l = l.replace("\n","").split()
			chromosomes.append(l[0])
	ff.close()

	args = [(sample_ID, bamfile, chrom, "Giraffe_Results/3_GC_bias") for chrom in chromosomes]
	with multiprocessing.Pool(processes=num_cpus) as pool:
		pool.map(get_bin_depth, args)

	path = "Giraffe_Results/3_GC_bias"
	system(f"cat {path}/*_bin_depth.txt > {path}/{sample_ID}.bin_depth.txt")
	system(f"rm {path}/*_bin_depth.txt")

def compute_GC_bias(ref, bamfile, binsize, sample_ID, num_cpu):
	path="Giraffe_Results/3_GC_bias"
	if os.path.exists(f"{path}/bin.bed") and os.path.exists(f"{path}/bin_GC.txt"):
		manager_bin_depth(ref, sample_ID, bamfile, num_cpu)
	else:
		get_bin_bed(ref, binsize)
		manager_GC_content(ref, num_cpu)
		manager_bin_depth(ref, sample_ID, bamfile, num_cpu)

	system(f"rm {path}/bin.bed")
	system(f"rm {path}/*_bin.bed")

def merge_GC_content_and_depth(binsize, sample_ID):
	data = {}
	input_depth = "Giraffe_Results/3_GC_bias/" + str(sample_ID) + ".bin_depth.txt"
	with open(input_depth) as f1:
		for bins in f1:
			bins = bins.replace("\n", "")
			bins = bins.split("\t")
			if bins[-1] != 0:
				KEY = bins[0] + "_" + bins[1] + "_" + bins[2]
				data[KEY]= {}
				data[KEY]["dp"] = int(bins[3]) /  int(binsize)
	f1.close()

	with open("Giraffe_Results/3_GC_bias/bin_GC.txt") as f2:
		for bins in f2:
			bins = bins.replace("\n", "")
			bins = bins.split("\t")		
			KEY = bins[0] + "_" + bins[1] + "_" + bins[2]
			if KEY in data.keys():
				data[KEY]["GC"] = float(bins[3]) * 100
			else:
				continue
	f2.close()

	merged_data = {}
	merged_data["dp"] = []
	merged_data["GC"] = []
	for i in data.keys():
		tmp_dp = data[i]["dp"]
		tmp_gc = data[i]["GC"]
		merged_data["dp"].append(tmp_dp)
		merged_data["GC"].append(tmp_gc)
	merged_data = pd.DataFrame.from_dict(merged_data)

	output_file = "Giraffe_Results/3_GC_bias/" + str(sample_ID) + "_relationship_raw.txt"
	ff = open(output_file, "w")
	ff.write("GC_content\tDepth\tNumber\tGroup\n")
	for i in range(0,101):
		tmp = merged_data[(i-0.5 <= merged_data["GC"]) & (merged_data["GC"] < i+0.5)].copy()
		if len(tmp) != 0:
			ave_dp = tmp["dp"].mean()
		else:
			ave_dp = 0.0
		ff.write(str(i) + "\t" + str(ave_dp) + "\t" + str(len(tmp)) + "\t" + str(sample_ID) + "\n")
	ff.close()

	# get the 95% data for downstream normalization
	df = pd.read_csv(output_file, sep=r'\s+')
	# df = pd.read_csv(output_file, delim_whitespace=True)
	max_number = df["Number"].max()
	total_number = df["Number"].sum()
	porportion = 0.95
	tmp = df[df["Number"] == max_number].copy()
	
	if len(tmp) == 1:
		for i in tmp["GC_content"]:
			start = i
			end = i
		
		for i in range(1,51):
			t1 = df[(start-1 <=df["GC_content"]) & (df["GC_content"] <= end+1)].copy()
			if t1["Number"].sum() / total_number >= porportion:
				nor_df = t1
				break
			else:
				start -= 1
				end += 1
				continue

	# normalization
	ave_dp = nor_df["Depth"].mean()
	nor_df["Normalized_depth"] = nor_df.apply(lambda row: row["Depth"]/ave_dp, axis=1)
	output_file_1 = "Giraffe_Results/3_GC_bias/" + str(sample_ID) + "_relationship_tmp.txt"
	nor_df.to_csv(output_file_1, sep="\t", index=False, header=False)

	system("rm Giraffe_Results/3_GC_bias/*bin_depth.txt")

def merge_files():
	with open("header", "w") as ff:
		ff.write("GC_content\tDepth\tNumber\tGroup\tNormalized_depth\n")
	ff.close()
	system("cat header Giraffe_Results/3_GC_bias/*_relationship_tmp.txt \
			> Giraffe_Results/3_GC_bias/Relationship_normalization.txt")
	system("rm header Giraffe_Results/3_GC_bias/*_relationship_tmp.txt")

def get_bin_number_within_GC_content():
	df = pd.read_table("Giraffe_Results/3_GC_bias/bin_GC.txt", header=None)
	with open("Giraffe_Results/3_GC_bias/Bin_distribution.txt", "w") as ff:
		ff.write("GC_content\tNumber\n")
		df[3] = df[3] * 100
		for i in range(0,101):
			tmp = df[(i-0.5 <= df[3]) & (df[3] < i+0.5)].copy()
			ff.write(str(i) + "\t" + str(len(tmp)) + "\n")
	ff.close()
	system("rm Giraffe_Results/3_GC_bias/bin_GC.txt")
