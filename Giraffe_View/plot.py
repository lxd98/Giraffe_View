from plotnine import * 
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

def plot_observe():
	df=pd.read_csv("./results/quality/final_observed_accuracy.txt", sep="\t")
	df=pd.DataFrame(df)

	df["Acc"] = df["Acc"] * 100
	accurac_scale = [0,100]
	accurac_break = [i for i in range(0, 101, 10)]
	mismatch_scale = [0,21]
	mismatch_break = [i for i in range(0, 21, 5)]
	
	Acc = (
		ggplot(df, aes(x="Acc")) + geom_density(adjust=10, size=1) +
	    	scale_x_continuous(name="Observed read accuracy", limits=accurac_scale, breaks=accurac_break) +
	    	scale_y_continuous(name="Probability Density Function (PDF)") + theme_bw() + 
	   		# geom_vline(aes(xintercept=90), color="#990000", linetype="dashed") +
	    	# geom_vline(aes(xintercept=99), color="#990000", linetype="dashed") +
	    	theme(axis_text=element_text(size=12, color="black"),
	          axis_title=element_text(size=12, color="black"))
	    )

	df["p_ins"] = 100 * df["Ins"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_del"] = 100 * df["Del"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_sub"] = 100 * df["Sub"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df1 = df[["p_ins", "p_del", "p_sub"]].melt()

	Mis = (
		ggplot(df1, aes(x="variable", y="value")) + geom_boxplot(outlier_size=0.5) + theme_bw() +
	    	scale_y_continuous(name="Mismatch proportion (%)", limits=mismatch_scale, breaks=mismatch_break) +
	    	scale_x_discrete(name="Mismatch type", labels=["Deletion", "Insertion", "Substitution"]) +
	    	theme(axis_text=element_text(size=12, color="black"),
	        	  axis_title=element_text(size=12, color="black"))
		)

	df2=pd.read_csv("./results/quality/final_homo_summary.txt", sep="\t")
	df2["value"] = df2["value"] * 100
	df2["group"] = "group"

	Homo = (
		ggplot(df2, aes(x="base",y="value", group="group"))  + 
		    geom_line(size=1) + geom_point(size=2) +
		    scale_y_continuous(name="Accuracy of homopolymer identification (%)", 
		    	limits=accurac_scale, breaks =accurac_break) +
		    xlab("Base") + theme_bw() + 
		    theme(axis_text=element_text(size=12, color="black"),
		          axis_title=element_text(size=12, color="black"))
		)

	Acc.save(filename = "Observed_read_accuracy.pdf", width=8, height=6, dpi=300, path="results/quality/")
	Mis.save(filename = "Observed_mismatch_proportion.pdf", width=8, height=6, dpi=300, path="results/quality/")
	Homo.save(filename = "Observed_homopolymer_identification.pdf", width=8, height=6, dpi=300, path="results/quality/")

def plot_estimate():
	df=pd.read_csv("results/estimated/final_estimated_accuracy.txt", sep="\t")
	df=pd.DataFrame(df)
	df["acc"] = df["acc"] * 100
	
	accurac_scale = [0,100]
	accurac_break = [i for i in range(0, 101, 10)]

	Acc = (
		ggplot(df, aes(x="acc")) + geom_density(adjust=10, size=1) +
            scale_x_continuous(name="Estimated read accuracy (%)", limits=accurac_scale, breaks=accurac_break) + 
            scale_y_continuous(name="Probability Density Function (PDF)") + theme_bw() +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"))
		)

	Acc.save(filename = "estimated_read_accuracy.pdf", width=8, height=6, dpi=300, path="results/estimated")

def plot_GC_bias():
	df=pd.read_csv("results/GC_bias/GC_bias_raw.txt", sep="\t")
	df=pd.DataFrame(df)

	accurac_scale = [0,100]
	accurac_break = [i for i in range(0, 101, 10)]

	distribution_len = (
		ggplot(df, aes(x="GC_content", y="number")) + geom_point(size=2) + geom_line(size=1) +
            scale_x_continuous(name = "GC content (%)", limits=accurac_scale, breaks =accurac_break) +
            scale_y_continuous(name="Number of bins") + theme_bw() +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"))
		)

	df1=pd.read_csv("results/GC_bias/final_GC_bias_nor.txt", sep="\t")
	df1=pd.DataFrame(df1)

	GC_bias = (
		ggplot(df1, aes(x="GC_content", y="nor_dp")) + geom_point(size=2) + geom_line(size=1) +
			scale_y_continuous(name="Normalized depth", 
				limits=[0, 2], breaks=np.arange(0, 2.1, 0.2)) +
			theme_bw() + 
			scale_x_continuous(name="GC content (%)", breaks=[i for i in range(0, 101, 1)]) +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"))
		)

	distribution_len.save(filename = "bin_distribution.pdf", width=8, height=6, dpi=300, path="results/GC_bias/")
	GC_bias.save(filename = "GC_bias.pdf", width=8, height=6, dpi=300, path="results/GC_bias/")