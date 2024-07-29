import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import os
from Giraffe_View.function import process_in_chunks

warnings.filterwarnings('ignore')

def plot_estimate(format='svg', path='Giraffe_Results/1_Estimated_quality'):
	df = process_in_chunks("Giraffe_Results/1_Estimated_quality/Estimated_information.txt")
	df = pd.DataFrame(df)
	df["Accuracy"] = df["Accuracy"] * 100
	df["GC_content"] = df["GC_content"] * 100
	df["Length"] = df["Length"] / 1000

	min_1 = df["Accuracy"].min()
	if min_1 >= 95:
		acc_scale = [95, 100]
		acc_breaks = [i for i in range(95, 101, 0.5)]
	elif min_1 >= 90:
		acc_scale = [90, 100]
		acc_breaks = [i for i in range(90, 101, 1)]
	elif min_1 >= 80:
		acc_scale = [80, 100]
		acc_breaks = [i for i in range(80, 101, 2)]
	elif min_1 >= 70:
		acc_scale = [70, 100]
		acc_breaks = [i for i in range(70, 101, 5)]
	elif min_1 >= 60:
		acc_scale = [60, 100]
		acc_breaks = [i for i in range(60, 101, 5)]
	elif min_1 >= 50:
		acc_scale = [50, 100]
		acc_breaks = [i for i in range(50, 101, 5)]
	elif min_1 >= 40:
		acc_scale = [40, 100]
		acc_breaks = [i for i in range(40, 101, 10)]
	elif min_1 >= 30:
		acc_scale = [30, 100]
		acc_breaks = [i for i in range(30, 101, 10)]
	elif min_1 >= 20:
		acc_scale = [20, 100]
		acc_breaks = [i for i in range(20, 101, 10)]
	elif min_1 >= 10:
		acc_scale = [10, 100]
		acc_breaks = [i for i in range(10, 101, 10)]
	else:
		acc_scale = [0, 100]
		acc_breaks = [i for i in range(0, 101, 5)]

	# plot
	plt.figure(figsize=(8, 6))
	ax = sns.kdeplot(data=df, x="Accuracy", hue="Group", fill=True, 
		alpha=0.6, palette = "Set2", common_norm=False)
	sns.move_legend(ax, "upper left")	
	ax
	plt.xlabel("Estimated read accuracy (%)")
	plt.ylabel("Probability Density Function")
	plt.xlim(acc_scale)
	plt.xticks(acc_breaks)

	plt.tight_layout()
	plt.savefig(f"{path}/1_Read_estimate_accuracy.{format}", format=format, dpi=300)
	plt.close()

	plt.figure(figsize=(8, 4))
	sns.boxplot(data=df, y="Group", x="GC_content", hue="Group", 
		palette="Set2", dodge=False, showfliers=False)
	
	plt.xlabel("GC content (%)")
	plt.xlim(0, 101)
	plt.xticks(range(0, 101, 10))
	
	plt.yticks()	
	plt.legend([],[], frameon=False)  # Hide the legend
	plt.tight_layout()
	plt.savefig(f"{path}/2_Read_GC_content.{format}", format=format, dpi=300)
	plt.close()

	ave = df["Length"].mean()
	if ave <= 1:
		len_scale = [0, 5]
		len_breaks = [i for i in range(0, 6, 0.5)]
	elif ave <= 5:
		len_scale = [0, 10]
		len_breaks = [i for i in range(0, 11, 1)]
	elif ave <= 10:
		len_scale = [0, 20]
		len_breaks = [i for i in range(0, 21, 2)]
	elif ave <= 20:
		len_scale = [0, 30]
		len_breaks = [i for i in range(0, 31, 5)]
	elif ave <= 30:
		len_scale = [0, 50]
		len_breaks = [i for i in range(0, 51, 5)]
	else:
		len_scale = [0, 100]
		len_breaks = [i for i in range(0, 101, 10)]

	plt.figure(figsize=(8, 6))
	sns.kdeplot(data=df, x="Length", hue="Group", 
		fill=True, common_norm=False, alpha=0.6,
		palette="Set2")

	plt.xlabel("Read length (Kb)")
	plt.ylabel("Probability Density Function")
	plt.xlim(len_scale)
	plt.xticks(len_breaks)
	plt.tight_layout()
	plt.savefig(f"{path}/3_Read_length.{format}", format=format, dpi=300)
	plt.close()

def plot_observe_acc(format='svg', path='Giraffe_Results/2_Observed_quality'):
	# color_set = "Set2"
	df = process_in_chunks("Giraffe_Results/2_Observed_quality/Observed_information.txt")
	df = pd.DataFrame(df)

	df["Acc"] = df["Acc"] * 100
	min_1 = df["Acc"].min()

	min_1 = df["Acc"].min()

	if min_1 >= 95:
		acc_scale = [95, 100]
		acc_breaks = [i for i in range(95, 101, 0.5)]
	elif min_1 >= 90:
		acc_scale = [90, 100]
		acc_breaks = [i for i in range(90, 101, 1)]
	elif min_1 >= 80:
		acc_scale = [80, 100]
		acc_breaks = [i for i in range(80, 101, 2)]
	elif min_1 >= 70:
		acc_scale = [70, 100]
		acc_breaks = [i for i in range(70, 101, 5)]
	elif min_1 >= 60:
		acc_scale = [60, 100]
		acc_breaks = [i for i in range(60, 101, 5)]
	elif min_1 >= 50:
		acc_scale = [50, 100]
		acc_breaks = [i for i in range(50, 101, 5)]
	elif min_1 >= 40:
		acc_scale = [40, 100]
		acc_breaks = [i for i in range(40, 101, 10)]
	elif min_1 >= 30:
		acc_scale = [30, 100]
		acc_breaks = [i for i in range(30, 101, 10)]
	elif min_1 >= 20:
		acc_scale = [20, 100]
		acc_breaks = [i for i in range(20, 101, 10)]
	elif min_1 >= 10:
		acc_scale = [10, 100]
		acc_breaks = [i for i in range(10, 101, 10)]
	else:
		acc_scale = [0, 100]
		acc_breaks = [i for i in range(0, 101, 5)]

	# Plot density plot for observed read accuracy
	# sns.set(style='darkgrid')
	plt.figure(figsize=(8, 6))

	ax = sns.kdeplot(data=df, x="Acc", hue="Group", fill=True, 
		common_norm=False, alpha=0.6, palette = "Set2")
	sns.move_legend(ax, "upper left")	
	ax
	
	plt.xlabel("Observed read accuracy (%)")
	plt.ylabel("Probability Density Function")
	plt.xlim(acc_scale)
	plt.xticks(acc_breaks)

	plt.tight_layout()
	plt.savefig(f"{path}/1_Observed_read_accuracy.{format}", format=format, dpi=300)
	plt.close()

	# Compute mismatch proportions
	df["p_ins"] = 100 * df["Ins"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_del"] = 100 * df["Del"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_sub"] = 100 * df["Sub"] / (df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])

	# Melt the dataframe for mismatch proportions
	df1 = pd.melt(df, id_vars=['Group'], value_vars=['p_ins', 'p_del', 'p_sub'])
	max_1 = df["p_ins"].max()
	max_2 = df["p_del"].max()
	max_3 = df["p_sub"].max()
	max_4 = max(max_1, max_2, max_3)

	if max_4 <= 5:
		mis_scale = [0, 5]
		mis_breaks = [i for i in range(0, 6, 0.5)]
	elif max_4 <= 10:
		mis_scale = [0, 10]
		mis_breaks = [i for i in range(0, 11, 1)]
	elif max_4 <= 20:
		mis_scale = [0, 20]
		mis_breaks = [i for i in range(0, 21, 2)]
	elif max_4 <= 30:
		mis_scale = [0, 30]
		mis_breaks = [i for i in range(0, 31, 5)]
	elif max_4 <= 40:
		mis_scale = [0, 40]
		mis_breaks = [i for i in range(0, 41, 5)]
	elif max_4 <= 50:
		mis_scale = [0, 50]
		mis_breaks = [i for i in range(0, 51, 5)]
	elif max_4 <= 60:
		mis_scale = [0, 60]
		mis_breaks = [i for i in range(0, 61, 10)]
	elif max_4 <= 70:
		mis_scale = [0, 70]
		mis_breaks = [i for i in range(0, 71, 10)]
	elif max_4 <= 80:
		mis_scale = [0, 80]
		mis_breaks = [i for i in range(0, 81, 10)]
	elif max_4 <= 90:
		mis_scale = [0, 90]
		mis_breaks = [i for i in range(0, 91, 10)]
	else:
		mis_scale = [0, 100]
		mis_breaks = [i for i in range(0, 101, 10)]

	# Plot boxplot for mismatch proportions
	plt.figure(figsize=(8, 6))
	sns.boxplot(data=df1, x="variable", y="value", hue="Group", 
		showfliers=False, width=0.5, gap=0.1, saturation=0.6, 
		palette = "Set2", linecolor="black")

	plt.ylabel("Mismatch proportion (%)")
	plt.ylim(mis_scale)
	plt.yticks(mis_breaks)
	plt.xticks(ticks=[0, 1, 2], labels=["Deletion", "Insertion", "Substitution"])
	plt.xlabel("")

	# Ensure the legend is created correctly
	handles, labels = plt.gca().get_legend_handles_labels()
	if not handles:
		for group in df["Group"].unique():
			handle = plt.Line2D([0], [0], color=sns.color_palette("pastel")[0], lw=2)
			handles.append(handle)
			labels.append(group)
	
	plt.legend(handles=handles, labels=labels, title='Group')
	plt.tight_layout()
	plt.savefig(f"{path}/2_Observed_mismatch_proportion.{format}", format=format, dpi=300)
	plt.close()

def plot_observe_homo(format='svg', path='Giraffe_Results/2_Observed_quality'):
	# Load data
	df = process_in_chunks("Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt")
	df["Accuracy"] = df["Accuracy"] * 100

	# Determine scale and breaks for y-axis based on minimum accuracy
	min_acc = df["Accuracy"].min()
	max_acc = df["Accuracy"].max()

	min_value = int((min_acc//10) * 10)
	max_value = int((max_acc//10) * 10 + 10)
	homo_scale = [min_value, max_value]

	dif = max_value - min_value

	if dif <= 10:
		homo_breaks = [i for i in range(min_value, max_value+1, 1)]
	elif dif <= 20:
		homo_breaks = [i for i in range(min_value, max_value+1, 2)]
	elif dif <= 30:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 40:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 50:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 60:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 70:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 80:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	elif dif <= 90:
		homo_breaks = [i for i in range(min_value, max_value+1, 5)]
	else:
		homo_breaks = [i for i in range(min_value, max_value+1, 10)]

	# Create the plot
	plt.figure(figsize=(8, 6))
	sns.lineplot(data=df, x='Base', y='Accuracy', hue='Group', linewidth=1.5,
		markers=True, dashes=False, palette = "Set2", alpha=0.6, legend=False)
	sns.scatterplot(data=df, x='Base',y='Accuracy', hue='Group', 
		palette = "Set2", s=50, edgecolor="black")
	
	# Customize plot
	plt.ylim(homo_scale)
	plt.yticks(homo_breaks)
	plt.ylabel('Accuracy of homopolymer identification (%)')
	plt.xlabel('Base')

	# Save plot
	output_path = f"{path}/3_Homoploymer_summary.{format}"
	plt.savefig(output_path, format=format, dpi=300, bbox_inches='tight')
	plt.close()

def plot_GC_bias(input_binsize, format='svg', path='Giraffe_Results/3_GC_bias'):
	# sns.set_style("whitegrid")
	# Load the first dataset
	df = pd.read_csv("Giraffe_Results/3_GC_bias/Bin_distribution.txt", sep="\t")
	accuracy_scale = [0, 100]
	accuracy_breaks = [i for i in range(0, 101, 10)]


	# Plot distribution length
	plt.figure(figsize=(8, 5))

	sns.lineplot(data=df, x="GC_content", y="Number", color="#96D1E8", linewidth=1.5, alpha=0.3)
	sns.scatterplot(data=df, x="GC_content", y="Number", color="#96D1E8", edgecolor="black", s=15)

	plt.xlim(accuracy_scale)
	plt.xticks(accuracy_breaks)
	plt.xlabel("GC content (%)")
	plt.ylabel(f"Number of bins (bin size = {input_binsize} bp)")
	plt.grid(False)
	plt.savefig(f"{path}/1_Bin_distribution.{format}", dpi=300)
	plt.close()

	# Load the second dataset
	df1 = pd.read_csv("Giraffe_Results/3_GC_bias/Relationship_normalization.txt", sep="\t")

	# Plot GC bias
	plt.figure(figsize=(8, 5))
	sns.lineplot(data=df1, x="GC_content", y="Normalized_depth", hue="Group", 
	palette = "Set2", linewidth=1.5, alpha=.6)
	sns.scatterplot(data=df1, x="GC_content", y="Normalized_depth", hue="Group",
	palette = "Set2", edgecolor="black", s=20, legend=False)

	plt.axhline(1, color="grey", linestyle="dotted")
	plt.ylim(0, 2)

	depth_breaks = (0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
	plt.yticks(depth_breaks)

	plt.xlabel("GC content (%)")
	plt.ylabel("Normalized depth")
	plt.grid(False)
	plt.savefig(f"{path}/2_Relationship_normalization.{format}", dpi=300)
	plt.close()

def plot_modi_bin(format="svg", path="Giraffe_Results/4_Regional_modification"):

    df = pd.read_csv("Giraffe_Results/4_Regional_modification/Regional_methylation_proportion.txt", sep="\t", names=["ID", "Value", "Group"])
    df["Value"] = df["Value"] * 100

    # sns.set(style="whitegrid")
    plt.figure(figsize=(20, 5))

    # Violin plot
    sns.violinplot(data=df, y="Group", x="Value",
     width=0.5, alpha=0.7, inner="box", split=True,
     inner_kws=dict(box_width=8, whis_width=2, color=".8"))

    methyl_scale = (0,100)
    methyl_breaks = [i for i in range(0, 101, 10)]

    plt.xlim(methyl_scale)
    plt.xticks(methyl_breaks)
    plt.xlabel("Methylation proportion (%)")
    plt.ylabel("")
    plt.yticks(fontsize=12, color='black')
    plt.legend([],[], frameon=False)  # Remove legend

    plt.savefig(f"{path}/1_Regional_modification.{format}", dpi=300)
    plt.close()