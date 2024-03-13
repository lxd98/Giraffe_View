from os import system
from plotnine import * 
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

def plot_estimate():
	df=pd.read_csv("Giraffe_Results/1_Estimated_quality/Estimated_information.txt", sep="\t")
	df=pd.DataFrame(df)
	df["Accuracy"] = df["Accuracy"] * 100
	df["GC_content"] = df["GC_content"] * 100
	df["Length"] = df["Length"] / 1000

	min_1 = df["Accuracy"].min()
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
		ggplot(df, aes(x="Accuracy", fill="Group")) + 
			geom_density(adjust=1, size=0.5, alpha=0.5) +
            scale_x_continuous(name="Estimated read accuracy (%)",
            	limits=acc_scale, breaks=acc_breaks) + 
            scale_y_continuous(name="Probability Density Function (PDF)") + 
            theme_classic() + xlab("") +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"),
                  legend_text = element_text(size=12, color="black"),
                  legend_title = element_blank(),
                  legend_position = "right")
            )

	GC_content=(
		ggplot(df, aes(x="Group", y="GC_content", fill="Group")) + 
			geom_boxplot(alpha=0.5, width=0.2, outlier_shape='') + 
			 # outlier_size=0.35) +
            scale_y_continuous(name="GC content (%)", 
            	limits=(0,101), breaks=[i for i in range(0,101,10)]) + 
            theme_classic() +
            theme(axis_text_x=element_text(size=12, color="black", angle=45),
            	axis_text_y=element_text(size=12, color="black"),
                	axis_title_x=element_blank(),
 			axis_title_y=element_text(size=12, color="black"),
                	legend_position = 'none')
            )
	
	ave = df["Length"].mean()

	if ave <= 1:
		len_scale = [0,5]
		len_breaks = [i for i in range(0, 11, 0.5)]
	elif ave <= 5:
		len_scale = [0,10]
		len_breaks = [i for i in range(0, 11, 1)]
	elif ave <= 10:
		len_scale = [0,20]
		len_breaks = [i for i in range(0, 21, 2)]
	elif ave <= 20:
		len_scale = [0,30]
		len_breaks = [i for i in range(0, 31, 5)]
	elif ave <= 30:
		len_scale = [0,50]
		len_breaks = [i for i in range(0, 51, 5)]
	else:
		len_scale = [0,100]
		len_breaks = [i for i in range(0, 101, 10)]

	Length = (
		ggplot(df, aes(x="Length", fill="Group")) + 
			geom_density(adjust=1, size=0.5, alpha=0.5) +
            scale_x_continuous(name="Read length (Kb)", 
            	limits=len_scale, breaks=len_breaks) + 
            scale_y_continuous(name="Probability Density Function (PDF)") + 
            theme_classic() + xlab("") +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"),
                  legend_text = element_text(size=12, color="black"),
                  legend_title = element_blank(),
                  legend_position = "right")
            )

	Acc.save(filename = "1_Read_accuracy.pdf", 
		width=8, height=6, dpi=300, path="Giraffe_Results/1_Estimated_quality")
	Length.save(filename = "2_Read_length.pdf", 
		width=8, height=6, dpi=300, path="Giraffe_Results/1_Estimated_quality")
	GC_content.save(filename = "3_Read_GC_content.pdf", 
		width=4, height=6, dpi=300, path="Giraffe_Results/1_Estimated_quality")

def plot_observe_acc():
	df=pd.read_csv("Giraffe_Results/2_Observed_quality/Observed_information.txt", sep="\t")
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
		ggplot(df, aes(x="Acc", fill="Group")) + 
			geom_density(adjust=1, size=0.5, alpha=0.5) +
	    	scale_x_continuous(name="Observed read accuracy (%)", 
	    		limits=acc_scale, breaks=acc_breaks) +
	    	scale_y_continuous(name="Probability Density Function (PDF)") +
	    	theme_classic() + 
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"),
                  legend_text = element_text(size=12, color="black"),
                  legend_title = element_blank(),
                  legend_position = "right")
            )

	Acc.save(filename = "1_Observed_read_accuracy.pdf", 
		width=8, height=6, dpi=300, path="Giraffe_Results/2_Observed_quality/")

	df["p_ins"] = 100 * df["Ins"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_del"] = 100 * df["Del"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	df["p_sub"] = 100 * df["Sub"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	
	df1 = pd.melt(df, id_vars=['Group'], value_vars=['p_ins', 'p_del', 'p_sub'])
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
		ggplot(df1, aes(x="variable", y="value", fill="Group")) + 
			geom_boxplot(alpha=0.5, width=0.2, outlier_shape='') +
				# outlier_size=0.35) +
			theme_classic() +
	    	scale_y_continuous(name="Mismatch proportion (%)",breaks=mis_breaks) +
	    	scale_x_discrete(name="Mismatch type", 
	    		labels=["Deletion", "Insertion", "Substitution"]) +
	    	theme(axis_text=element_text(size=12, color="black"),
	        	axis_title_y=element_text(size=12, color="black"),
	        	axis_title_x=element_blank(),
	        	legend_text = element_text(size=12, color="black"),
	        	legend_title = element_blank(),
	        	legend_position = "right")
		)

	Mis.save(filename = "2_Observed_mismatch_proportion.pdf",
		width=8, height=6, dpi=300, path="Giraffe_Results/2_Observed_quality/")

def plot_observe_homo():
	df=pd.read_csv("./Giraffe_Results/2_Observed_quality/Homoploymer_summary.txt", sep="\t")
	df["Accuracy"] = df["Accuracy"] * 100

	min_2 = df["Accuracy"].min()
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
		ggplot(df, aes(x="Base",y="Accuracy", fill="Group", color="Group", group="Group")) + 
			geom_line(size=1, alpha=.8) + 
			geom_point(size=3, color="black", alpha=.8) +
		    scale_y_continuous(name="Accuracy of homopolymer identification (%)", 
		    	limits=homo_scale, breaks =homo_breaks) +
		    xlab("Base") + theme_classic() + 
		    theme(axis_text=element_text(size=12, color="black"),
		          axis_title_y=element_text(size=12, color="black"),
		          axis_title_x=element_blank(),
		          legend_text = element_text(size=12, color="black"),
		          legend_title = element_blank(),
		          legend_position = "right")
		    )

	Homo.save(filename = "3_Homoploymer_summary.pdf", 
		width=8, height=6, dpi=300, path="Giraffe_Results/2_Observed_quality/")

def plot_GC_bias():
	df=pd.read_csv("Giraffe_Results/3_GC_bias/Bin_distribution.txt", sep="\t")
	accurac_scale = [0,100]
	accurac_break = [i for i in range(0, 101, 10)]

	distribution_len=(
		ggplot(df, aes(x="GC_content", y="Number")) + 
			geom_line(size=1.5, color="#96D1E8", alpha=.3) + 
			geom_point(size=1.5, fill="#96D1E8", color="black") +
            scale_x_continuous(name = "GC content (%)", 
            	limits=accurac_scale, breaks =accurac_break) +
            scale_y_continuous(name="Number of bins") + theme_classic() +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"))
            )

	df1=pd.read_csv("Giraffe_Results/3_GC_bias/Relationship_normalization.txt", sep="\t")

	dif = df1["GC_content"].max() - df1["GC_content"].min()
	if dif <= 15:
		gc_breaks = [i for i in range(0, 101, 1)]
	elif dif <= 30:
		gc_breaks = [i for i in range(0, 101, 2)]
	elif dif <= 50:
		gc_breaks = [i for i in range(0, 101, 5)]
	else:
		gc_breaks = [i for i in range(0, 101, 10)]

	GC_bias=(
		ggplot(df1, aes(x="GC_content", y="Normalized_depth", 
			group="Group", fill="Group", color="Group")) + 
			geom_hline(aes(yintercept=1), color="grey", linetype="dotted") + 
			geom_line(size=1.5, alpha=.3) + 
			geom_point(size=1.5,color="black") +
			scale_y_continuous(name="Normalized depth", 
				limits=[0, 2], breaks=np.arange(0, 2.1, 0.2)) +
			theme_classic() + 
			scale_x_continuous(name="GC content (%)", breaks=gc_breaks) +
            theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"),
                  legend_title = element_blank(),
                  legend_text = element_text(size=12, color="black"),
                  legend_position = "bottom"
                  )
            )

	distribution_len.save(filename = "1_Bin_distribution.pdf", 
		width=8, height=3, dpi=300, path="Giraffe_Results/3_GC_bias")
	GC_bias.save(filename = "2_Relationship_normalization.pdf", 
		width=8, height=5, dpi=300, path="Giraffe_Results/3_GC_bias")

def plot_modi_bin():
	system('cat Giraffe_Results/4_Regional_modification/*bed | grep -v -w "nan" \
		> Giraffe_Results/4_Regional_modification/plot.txt')

	df=pd.read_csv("Giraffe_Results/4_Regional_modification/plot.txt", 
		sep="\t", names=["ID", "Value", "Group"])
	
	df=df.dropna()
	df["Value"] = df["Value"] * 100
	diff = df["Value"].max() - df["Value"].min()

	if diff <= 20:
		Breaks = [i for i in range(0, 101, 2)]
	elif diff <= 50:
		Breaks = [i for i in range(0, 101, 5)]
	else:
		Breaks = [i for i in range(0, 101, 10)]

	Bin_box = (
		ggplot(df, aes(y="Value", x="Group", fill="Group")) + 
			geom_point(position = position_jitter(width=0.15), size = 0.5,alpha = 0.1) + 
			geom_boxplot(width=0.25,  alpha=0.5, outlier_shape="") +
			geom_violin(position = position_nudge(x = .2), alpha = 0.7, adjust = 0.5, style="right") + 
			coord_flip() +
			scale_y_continuous(name="Modification proportion (%)",
				breaks=Breaks) + xlab("") + theme_classic() +
			theme(axis_text=element_text(size=12, color="black"),
                  axis_title=element_text(size=12, color="black"),
                  legend_position = 'none')
		)


	Bin_box.save(filename = "1_Regional_modification.pdf", 
		width=8, height=6, dpi=300, path="Giraffe_Results/4_Regional_modification")

	system("rm Giraffe_Results/4_Regional_modification/plot.txt")