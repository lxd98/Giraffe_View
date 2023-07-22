from plotnine import * 
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

def plot_observe_acc():
	df=pd.read_csv("./results/observed_quality/final_observed_accuracy.txt", sep="\t")
	df=pd.DataFrame(df)

	df["Acc"] = df["Acc"] * 100
	min_1 = df["Acc"].min()
	if min_1 >= 75:
		acc_scale = [75,100]
		acc_breaks = [i for i in range(75, 101, 5)]
	elif min_1 >= 50:
		acc_scale = [50,100]
		acc_breaks = [i for i in range(50, 101, 10)]
	elif min_1 >= 20:
		acc_scale = [20,100]
		acc_breaks = [i for i in range(20, 101, 10)]
	else:
		acc_scale = [0,100]
		acc_breaks = [i for i in range(0, 101, 10)]

	Acc = (
		ggplot(df, aes(x="Acc")) + 
			geom_density(adjust=10, size=1.5, color="#96D1E8") +
	    	scale_x_continuous(name="Observed read accuracy", 
	    		limits=acc_scale, breaks=acc_breaks) +
	    	scale_y_continuous(name="Probability Density Function (PDF)") +
	    	theme_bw() + 
	   		# geom_vline(aes(xintercept=90), color="grey", linetype="dotted") +
	    	# geom_vline(aes(xintercept=99), color="grey", linetype="dotted") +
	    	theme(axis_text=element_text(size=12, color="black"),
	          axis_title=element_text(size=15, color="black"))
	    )

	Acc.save(filename = "Observed_read_accuracy.pdf", 
		width=8, height=6, dpi=300, path="results/observed_quality/")

	df["p_ins"] = 100 * df["Ins"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_del"] = 100 * df["Del"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_sub"] = 100 * df["Sub"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df1 = df[["p_ins", "p_del", "p_sub"]].melt()

	max_1 = df["p_ins"].max()
	max_2 = df["p_del"].max()
	max_3 = df["p_sub"].max()
	max_4 = max(max_1, max_2, max_3)

	if max_4 <= 20:
		mis_breaks = [i for i in range(0, 101, 1)]
	elif max_4 <= 50:
		mis_breaks = [i for i in range(0, 101, 2)]
	else:
		mis_breaks = [i for i in range(0, 101, 5)]		

	Mis = (
		ggplot(df1, aes(x="variable", y="value")) + 
			geom_boxplot(width=0.25, outlier_shape='', fill="#96D1E8", alpha=0.3) +
			theme_bw() +
	    	scale_y_continuous(name="Mismatch proportion (%)",breaks=mis_breaks) +
	    	scale_x_discrete(name="Mismatch type", 
	    		labels=["Deletion", "Insertion", "Substitution"]) +
	    	theme(axis_text=element_text(size=12, color="black"),
	        	  axis_title=element_text(size=15, color="black"))
		)

	Mis.save(filename = "Observed_mismatch_proportion.pdf", 
		width=8, height=6, dpi=300, path="results/observed_quality/")

def plot_observe_homo():
	df2=pd.read_csv("./results/observed_quality/final_homo_summary.txt", sep="\t")
	df2["value"] = df2["value"] * 100
	df2["group"] = "group"

	min_2 = df2["value"].min()
	if min_2 >= 75:
		homo_scale = [75,100]
		homo_breaks = [i for i in range(75, 101, 5)]
	elif min_2 >= 50:
		homo_scale = [50,100]
		homo_breaks = [i for i in range(50, 101, 5)]
	elif min_2 >= 25:
		homo_scale = [25,100]
		homo_breaks = [i for i in range(25, 101, 5)]
	else:
		homo_scale = [0,100]
		homo_breaks = [i for i in range(0, 101, 5)]

	Homo = (
		ggplot(df2, aes(x="base",y="value", group="group"))  + 
			geom_line(size=1, color="#96D1E8", alpha=.7, linetype="dashed") + 
			geom_point(size=3, fill="#96D1E8", color="grey") + 
		    
		    scale_y_continuous(name="Accuracy of homopolymer identification (%)", 
		    	limits=homo_scale, breaks =homo_breaks) +
		    xlab("Base") + theme_bw() + 
		    theme(axis_text=element_text(size=12, color="black"),
		          axis_title=element_text(size=15, color="black"))
		)

	Homo.save(filename = "Observed_homopolymer_identification.pdf", 
		width=8, height=6, dpi=300, path="results/observed_quality/")

def plot_estimate():
	df=pd.read_csv("results/estimated_quality/final_estimated_accuracy.txt", sep="\t")
	df=pd.DataFrame(df)
	df["acc"] = df["acc"] * 100

	min_1 = df["acc"].min()
	if min_1 >= 75:
		acc_scale = [75,100]
		acc_breaks = [i for i in range(75, 101, 5)]
	elif min_1 >= 50:
		acc_scale = [50,100]
		acc_breaks = [i for i in range(50, 101, 10)]
	elif min_1 >= 20:
		acc_scale = [20,100]
		acc_breaks = [i for i in range(20, 101, 10)]
	else:
		acc_scale = [0,100]
		acc_breaks = [i for i in range(0, 101, 10)]

	Acc = (
		ggplot(df, aes(x="acc")) + 
			geom_density(adjust=10, size=1.5, color="#96D1E8") +
            scale_x_continuous(name="Estimated read accuracy (%)",
            	limits=acc_scale, breaks=acc_breaks) + 
      #      	geom_vline(aes(xintercept=90), color="grey", linetype="dotted") +
	    	# geom_vline(aes(xintercept=99), color="grey", linetype="dotted") +
            scale_y_continuous(name="Probability Density Function (PDF)") + 
            theme_bw() +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=15, color="black"))
		)

	Acc.save(filename = "estimated_read_accuracy.pdf", 
		width=8, height=6, dpi=300, path="results/estimated_quality")

def plot_GC_bias():
	df=pd.read_csv("results/GC_bias/GC_bias_raw.txt", sep="\t")
	df=pd.DataFrame(df)

	accurac_scale = [0,100]
	accurac_break = [i for i in range(0, 101, 10)]

	distribution_len = (
		ggplot(df, aes(x="GC_content", y="number")) + 
			geom_line(size=2, color="#96D1E8", alpha=.3) + 
			geom_point(size=2, fill="#96D1E8", color="grey") +
            scale_x_continuous(name = "GC content (%)", 
            	limits=accurac_scale, breaks =accurac_break) +
            scale_y_continuous(name="Number of bins") + theme_bw() +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=15, color="black"))
		)

	df1=pd.read_csv("results/GC_bias/final_GC_bias_nor.txt", sep="\t")
	df1=pd.DataFrame(df1)

	dif = df1["GC_content"].max() - df1["GC_content"].min()
	if dif <= 15:
		gc_breaks = [i for i in range(0, 101, 1)]
	elif dif <= 30:
		gc_breaks = [i for i in range(0, 101, 2)]
	elif dif <= 50:
		gc_breaks = [i for i in range(0, 101, 5)]
	else:
		gc_breaks = [i for i in range(0, 101, 10)]

	GC_bias = (
		ggplot(df1, aes(x="GC_content", y="nor_dp")) + 
			geom_hline(aes(yintercept=1), color="grey", linetype="dotted") + 
			geom_line(size=2, color="#96D1E8", alpha=.3) + 
			geom_point(size=2, fill="#96D1E8", color="grey") +
			scale_y_continuous(name="Normalized depth", 
				limits=[0, 2], breaks=np.arange(0, 2.1, 0.2)) +
			theme_bw() + 
			scale_x_continuous(name="GC content (%)", breaks=gc_breaks) +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=15, color="black"))
		)

	distribution_len.save(filename = "bin_distribution.pdf", 
		width=8, height=6, dpi=300, path="results/GC_bias/")
	GC_bias.save(filename = "GC_bias.pdf", 
		width=8, height=6, dpi=300, path="results/GC_bias/")


def plot_modi_bin(INPUT, REF):
	filename = INPUT.replace('.bed', '') + "_" + REF.replace('.csv', '') + ".bed"
	df=pd.read_csv("results/regional_modification/" + str(filename), sep="\t", 
		names=["ID", "Value"])
	df=df.dropna()
	df=pd.DataFrame(df)
	df["group"] = "Bin"

	diff = df["Value"].max() - df["Value"].min()

	if diff <= 20:
		Breaks = [i for i in range(0, 101, 2)]
	elif diff <= 40:
		Breaks = [i for i in range(0, 101, 5)]
	else:
		Breaks = [i for i in range(0, 101, 10)]

	Bin_box = (
		ggplot(df, aes(x="group", y="Value")) + 
			geom_boxplot(width=0.15, fill="#96D1E8", alpha=0.2, outlier_shape="") +
			scale_y_continuous(name="Modification proportion (%)",
				breaks=Breaks) + xlab("") + theme_bw() +
			theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=15, color="black"))
		)

	Bin_violin = (
		ggplot(df, aes(x="group", y="Value")) + 
			geom_violin(width=0.15, fill="#96D1E8", alpha=0.3) + 
			theme_bw() +
			scale_y_continuous(name="Modification proportion (%)",
				breaks=Breaks) + xlab("") + 
			theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=15, color="black"))
		)

	Bin_box.save(filename = "motif_bin_box.pdf", 
		width=8, height=6, dpi=300, path="results/regional_modification/")

	Bin_violin.save(filename = "motif_bin_violin.pdf", 
		width=8, height=6, dpi=300, path="results/regional_modification/")