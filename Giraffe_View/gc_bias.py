import os
from os import system
import pandas as pd
from Giraffe_View.function import cmd_shell

def compute_GC_bias(ref, bamfile, binsize, sample_ID):
	def get_bin_bed():
		with open("tmp.sh","w") as bash_file:
			cmd1 = "samtools faidx " + str(ref)
			cmd2 = "bedtools makewindows -g " + str(ref) + ".fai -w " + str(binsize) + " > Giraffe_Results/3_GC_bias/bin.bed"
			bash_file.write(cmd1 + "\n")
			bash_file.write(cmd2 + "\n")
		bash_file.close()
		run = "bash tmp.sh"
		cmd_shell(str(run), "Bin size calculating")

	def get_bin_GC_content():
		with open("tmp_1.sh","w") as bash_file:
			cmd = 'bedtools nuc -fi ' + str(ref) + ' -bed Giraffe_Results/3_GC_bias/bin.bed'
			cmd += '> Giraffe_Results/3_GC_bias/bin.tmp'
			bash_file.write(str(cmd))
		bash_file.close()
		run = "bash tmp_1.sh"
		cmd_shell(str(run), "GC content calculating")

		output = open("Giraffe_Results/3_GC_bias/bin_GC_content.txt", "w")
		with open("Giraffe_Results/3_GC_bias/bin.tmp") as ff:
			for bins in ff:
				if bins[0] != "#":
					bins = bins.replace("\n", "")
					bins = bins.split("\t")
					output.write(bins[0] + "\t" + bins[1] + "\t" + bins[2] + "\t" + bins[4] + "\n")
		ff.close()
		output.close()

	def get_bin_depth():
		ff = open("tmp_2.sh", "w")
		mes = "samtools bedcov Giraffe_Results/3_GC_bias/bin.bed " + str(bamfile)
		mes += " > Giraffe_Results/3_GC_bias/" + str(sample_ID) + "_bin_depth.txt"
		ff.write(mes)
		ff.close()
		run = "bash tmp_2.sh"
		cmd_shell(str(run), "Depth calculating")

	if os.path.exists("Giraffe_Results/3_GC_bias/bin.bed") and os.path.exists("Giraffe_Results/3_GC_bias/bin_GC_content.txt"):
		get_bin_depth()
	else:
		get_bin_bed()
		get_bin_GC_content()
		get_bin_depth()
		system("rm ./Giraffe_Results/3_GC_bias/bin.tmp")


	system("rm tmp*sh")
	system("rm Giraffe_Results/3_GC_bias/bin.bed")

def merge_GC_content_and_depth(binsize, sample_ID):
	data = {}
	input_depth = "Giraffe_Results/3_GC_bias/" + str(sample_ID) + "_bin_depth.txt"
	with open(input_depth) as f1:
		for bins in f1:
			bins = bins.replace("\n", "")
			bins = bins.split("\t")
			if bins[-1] != 0:
				KEY = bins[0] + "_" + bins[1] + "_" + bins[2]
				data[KEY]= {}
				data[KEY]["dp"] = int(bins[3]) /  int(binsize)
	f1.close()

	with open("Giraffe_Results/3_GC_bias/bin_GC_content.txt") as f2:
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
	df = pd.read_csv(output_file, delim_whitespace=True)
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
				start -= 2
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
	df = pd.read_table("Giraffe_Results/3_GC_bias/bin_GC_content.txt", header=None)
	with open("Giraffe_Results/3_GC_bias/Bin_distribution.txt", "w") as ff:
		ff.write("GC_content\tNumber\n")
		df[3] = df[3] * 100
		for i in range(0,101):
			tmp = df[(i-0.5 <= df[3]) & (df[3] < i+0.5)].copy()
			ff.write(str(i) + "\t" + str(len(tmp)) + "\n")
	ff.close()
	system("rm Giraffe_Results/3_GC_bias/bin_GC_content.txt")