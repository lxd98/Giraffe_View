import rpy2.robjects as robjects

def plot_observe():
	rscript='''
	library(ggplot2)
	library(reshape2)
	library(patchwork)

	homo <- function(input_file){
	  df <- read.csv(input_file, header = T, sep = "\t")
	  df["group"] <- "test"
	  p <- ggplot(df, aes(x=base,y=value*100, group=group))  + 
	    geom_line(alpha=0.5) + geom_point(size=2) +
	    scale_y_continuous(name="Accuracy of homopolymer identification (%)", 
	                       limits = c(0,100), breaks = seq(0,100,10)) +
	    xlab("Base") + theme_bw() + 
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"),
	          legend.title = element_blank(), legend.position = "bottom")
	  return(p)
	}

	acc <- function(input_file) {
	  df <- read.csv(input_file, header = T, sep = "\t") 
	  p <- ggplot(df, aes(x=Acc*100)) + geom_density(adjust=10, alpha=0.8) +
	    scale_x_continuous(name="Observed read accuracy (%)", 
	                       limits = c(0,100), breaks = seq(0,100,10)) + 
	    scale_y_continuous(name="Probability Density Function (PDF)") + theme_bw() +
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"))
	  return(p)
	}

	mismatch <- function(input_file){
	  df <- read.csv(input_file, header = T, sep = "\t")
	  df["group"] <- "test"
	  df["p_ins"] <- 100*df["Ins"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	  df["p_del"] <- 100*df["Del"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	  df["p_sub"] <- 100*df["Sub"]/(df["Ins"] + df["Del"] + df["Sub"] + df["Mat"])
	  df1 <- df[, 8:11]
	  df2 <- melt(df1)
	  df2$variable <- factor(df2$variable, levels = c("p_ins", "p_del", "p_sub"))
	  p <- ggplot(df2, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA) + theme_bw() +
	    scale_y_continuous(name="Mismatch proportion (%)", 
	                       limits = c(0,100), breaks = seq(0,100,10)) +
	    scale_x_discrete(name="Type", labels=c("Ins", "Del", "Sub")) +
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"))
	  return(p)
	}

	A <- acc("results/quality/final_observed_accuracy.txt") + 
		mismatch("results/quality/final_observed_accuracy.txt") + 
		homo("results/quality/final_homo_summary.txt")
	
	ggsave("observe.png", A, path="./results/quality/", dpi=300, width=12, height=5)
	'''
	robjects.r(rscript)

def plot_estimate():
	rscript = '''
	library(ggplot2)
	library(reshape2)

	estimated <- function(input_file){
	  df <- read.csv(input_file, header = T, sep = "\t")
	  p <- ggplot(df, aes(x=acc*100)) + geom_density(adjust=10, alpha=0.8) +
	    scale_x_continuous(name="Estimated read accuracy (%)", 
	                       limits = c(0,100), breaks = seq(0,100,10)) + 
	    scale_y_continuous(name="Probability Density Function (PDF)") + theme_bw() +
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"))
	  return(p)
	}

	A <- estimated("results/estimated/final_estimated_accuracy.txt")
	ggsave("estimate.png", A, path="./results/estimated/", dpi=300, width=8, height=5)
	'''
	robjects.r(rscript)

def plot_GC_bias():
	rscript = '''
	library(ggplot2)
	library(patchwork)

	bin_count <- function(input_file){
	  df <- read.csv(input_file, sep = "\t", header = T)
	  p <- ggplot(df, aes(x=GC_content, y=number)) + geom_point(size=1) + geom_line() +
	    scale_x_continuous(name = "GC content (%)", 
	                       limits = c(0,100), breaks = seq(0,100,10)) +
	    scale_y_continuous(name="Number of bins") + theme_bw() +
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"),
	          legend.title = element_blank(), legend.position = "bottom")
	  return(p)
	}

	GC_bias <- function(input_file){
	  df <- read.csv(input_file, sep = "\t", header = T)
	  p <- ggplot(df, aes(x=GC_content, y=nor_dp)) + geom_point(size=1) + geom_line() +
	    scale_y_continuous(name="Normalized depth", limits = c(0.5,1.5), 
	                       breaks = seq(0.5,1.5,0.1)) + theme_bw() + 
	    scale_x_continuous(name="GC content (%)", breaks = seq(0,100,1)) +
	    theme(axis.text=element_text(size=12, family = "Arial", color="black"),
	          axis.title=element_text(size=12, family = "Arial", color="black"),
	          legend.text = element_text(size=12, family = "Arial", color="black"),
	          legend.title = element_blank(), legend.position = "bottom")
	  return(p)
	}

	A <- bin_count("results/GC_bias/GC_bias_raw.txt")/GC_bias("results/GC_bias/final_GC_bias_nor.txt")
	ggsave("GC_bias.png", A, path="./results/GC_bias/", dpi=300, width=8, height=5)
	'''
	robjects.r(rscript)