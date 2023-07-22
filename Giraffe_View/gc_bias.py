import pandas as pd
from Giraffe_View.function import cmd_shell

def compute_GC_bias(ref, bamfile, binsize):
	def get_GC_content():
		output = open("results/GC_bias/GC_content.txt", "w")
		with open("results/GC_bias/BIN.txt") as ff:
			for bins in ff:
				if bins[0] != "#":
					bins = bins.replace("\n", "")
					bins = bins.split("\t")
					output.write(bins[0] + "\t" + bins[1] + "\t" + bins[2] + "\t" + bins[4] + "\n")
		ff.close()
		output.close()		

	def dp():
		ff = open("GC.sh", "w")
		mes = "samtools bedcov results/GC_bias/BIN.bed " + str(bamfile) + " > results/GC_bias/BIN.dp"
		ff.write(mes)
		ff.close()
		run = "bash GC.sh"
		clean = "rm GC.sh"
		cmd_shell(str(run), "Depth calculating")
		cmd_shell(str(clean), "Clean")

	ff = open("GC.sh", "w")
	cmd1 = "samtools faidx " + str(ref)
	cmd2 = "bedtools makewindows -g " + str(ref) + ".fai -w " + str(binsize) + " > results/GC_bias/BIN.bed" 
	cmd3 = 'bedtools nuc -fi ' + str(ref) + ' -bed results/GC_bias/BIN.bed > results/GC_bias/BIN.txt'
	run = "bash GC.sh"
	clean = "rm results/GC_bias/BIN.txt " + str(ref) + ".fai GC.sh"

	ff.write(cmd1 + "\n")
	ff.write(cmd2 + "\n")
	ff.write(cmd3 + "\n")
	ff.close()

	cmd_shell(str(run), "GC content calculating")
	get_GC_content()
	# get_GC_db()
	cmd_shell(str(clean), "Clean")
	dp()


def merge_CG_content_and_depth(binsize):
	data = {}

	with open("results/GC_bias/BIN.dp") as f1:
		for bins in f1:
			bins = bins.replace("\n", "")
			bins = bins.split("\t")
			if bins[-1] != 0:
				KEY = bins[0] + "_" + bins[1] + "_" + bins[2]
				data[KEY]= {}
				data[KEY]["dp"] = int(bins[3]) /  int(binsize)
	f1.close()

	with open("results/GC_bias/GC_content.txt") as f2:
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

	ff = open("results/GC_bias/GC_bias_raw.txt", "w")
	ff.write("GC_content\tdp\tnumber\n")
	for i in range(0,101):
		tmp = merged_data[(i-0.5 <= merged_data["GC"]) & (merged_data["GC"] < i+0.5)].copy()
		if len(tmp) != 0:
			ave_dp = tmp["dp"].mean()
		else:
			ave_dp = 0.0
		ff.write(str(i) + "\t" + str(ave_dp) + "\t" + str(len(tmp)) + "\n")
	ff.close()

	#get the 95% data for downstream normalization
	df = pd.read_csv("results/GC_bias/GC_bias_raw.txt", delim_whitespace=True)
	max_number = df["number"].max()
	total_number = df["number"].sum()
	porportion = 0.90
	tmp = df[df["number"] == max_number].copy()
	
	if len(tmp) == 1:
		for i in tmp["GC_content"]:
			start = i
			end = i
		
		for i in range(1,51):
			t1 = df[(start-1 <=df["GC_content"]) & (df["GC_content"] <= end+1)].copy()
			if t1["number"].sum() / total_number >= porportion:
				nor_df = t1
				break
			else:
				start -= 1
				end += 1
				continue

	# normalization
	ave_dp = nor_df["dp"].mean()
	nor_df["nor_dp"] = nor_df.apply(lambda row: row["dp"]/ave_dp, axis=1)
	nor_df.to_csv("results/GC_bias/final_GC_bias_nor.txt", sep="\t", index=False)