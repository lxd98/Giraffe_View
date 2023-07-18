import sys
import argparse
from function import *
from homopolymer import *
from observed_read_accuracy import *
from GC_bias import *
from estimated_read_accuracy import *
from regional_modification import *
from plot import *

def observed(args):
    Data_process(args.input, args.ref, args.cpu)

    print("Start observed accuracy analysis!")
    observed_accuracy("results/quality/tmp.sort.bam")
    
    print("Start homopolymer identification!")
    homopolymer_from_bam("results/quality/tmp.sort.bam")
    homopolymer_summary_1("results/quality/homo.txt")
    homopolymer_summary_2()
    # Rplot("quality.R")

def methylation(args):
    mkdir_d("regional_modification")
    methylation_calculation(args.input, args.ref, args.cpu)

def GC_bias(args):
    mkdir_d("GC_bias")
    compute_GC_bias(args.ref, args.input, args.binsize)
    merge_CG_content_and_depth(args.binsize)

def estimated(args):
    mkdir_d("estimated")
    calculate_estimated_accuracy(args.input, args.cpu)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Giraffe_view", description="A tool to help you assess quality of your ONT data.")

    # Define subparsers
    subparsers = parser.add_subparsers(dest='function', help="")
    observed_parser = subparsers.add_parser('observe', help='Observed quality in accuracy, mismatch, and homopolymer')
    observed_parser.add_argument("--input", type=str, metavar="<fastq>", required=True, help="input reads")
    observed_parser.add_argument("--ref", type=str, metavar="<reference>", required=True, help="input reference")
    observed_parser.add_argument("--cpu", type=int, metavar="<number>", required=False, help="number of cpu (default:10)", default=10)

    methylation_parser = subparsers.add_parser('modi', help='Average modification proportion of regions')
    methylation_parser.add_argument("--input", type=str, metavar="<bed>", required=True, help="input bed file")
    methylation_parser.add_argument("--ref", type=str, metavar="<reference>", required=True, help="input reference")
    methylation_parser.add_argument("--cpu", type=int, metavar="<number>", required=False, help="number of cpu (default:10)", default=10)

    GC_bias_parser = subparsers.add_parser('GC_bias', help='Relationship between GC content and depth')
    GC_bias_parser.add_argument("--ref", type=str, metavar="<reference>", required=True, help="input reference file")
    GC_bias_parser.add_argument("--input", type=str, metavar="<sam/bam>", required=True, help="input bam/sam file")
    GC_bias_parser.add_argument("--binsize", type=int, metavar="", required=False, help="input bin size (default:1000)", default=1000)

    estimated_parser = subparsers.add_parser('estimate', help='Estimated read accuracy')
    estimated_parser.add_argument("--input", type=str, metavar="<fastq>", required=True, help="input reads")
    estimated_parser.add_argument("--cpu", type=int, metavar="<number>", required=False, help="number of cpu (default:10)", default=10)

    # Add function to print help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Call the appropriate function based on the subparser used
    if args.function == "observe":
        observed(args)
        plot_observe()

    elif args.function == "modi":
        methylation(args)
    
    elif args.function == "GC_bias":
        GC_bias(args)
        plot_GC_bias()
    
    elif args.function == "estimate":
        estimated(args)
        plot_estimate()