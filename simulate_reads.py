from SCRaMbLE_DNA_simulation import DNA_extraction_coverage
from get_segments import cut_reference_in_segments
from get_loxp_segment import cut_reference_in_LUs
from subpath_to_DNA1 import path_to_DNA
import sys
import argparse
import random

# This function takes as input a genome or chromosome in fasta format and output simulated reads in fasta format.
def SIM_reads(filename="IXR_BACnewseq.fa", ID="", synthetic=True, coverage=25, reads_len=8, sigma=3, circular=False):
    if ID == "":
        # Take the name before the dot. find() finds the first occurrence. rfind() finds the last occurrence.
        ID = filename[:filename.rfind(".")]
    # Cut the DNA in segments or LU and extract the LUs sequence.
    if synthetic:
        genome = cut_reference_in_LUs(filename=filename, ID=ID, starting_LU=1)
    #else:
    #    cut_reference_in_segments(filename=filename, ID=ID, starting_LU=starting_LU, segment_size=segment_size, min_segment_size=min_segment_size, mapping_breakpoints=mapping_breakpoints, RANDOM=RANDOM, random_cuts=random_cuts, essential=[])

    # Cuts the input genome
    reads = DNA_extraction_coverage(syn_chr=genome, coverage=coverage, reads_len=reads_len, sigma=sigma, circular=circular)

    # Convert the subpaths into DNA and save them in a fasta format.
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    path_to_DNA(path=reads, LU_fasta=ID+".loxpreg.fa", filename=ID + "_SIM_reads_C" + str(coverage) + "_R" + str(random_seed), ID="")

    return None

# test the code
if __name__ == '__main__':
    #test the script with the synIXR chromosome
    #SIM_reads(filename="IXR_BACnewseq.fa", synthetic=True, coverage=25)

    parser = argparse.ArgumentParser(description="This script simulates long sequencing reads from a genome.")

    parser.add_argument("-ref", "--reference", type=str, required=True, help="The reference fasta file containing a synthetic chromosome")
    parser.add_argument("-ID", "--ID", type=str, required=False, default="", help="The name of the output file")
    parser.add_argument("-synthetic", "--synthetic", type=bool, required=False, default=True, help="Is the reference genome synthetic? If so, the script will only cuts into loxPsym sites.")
    parser.add_argument("-coverage", "--coverage", type=int, required=False, default=25, help="The coverage of the sequencing reads")
    parser.add_argument("-reads_len", "--reads_len", type=int, required=False, default=8, help="The mean of the read length distribution")
    parser.add_argument("-sigma", "--sigma", type=int, required=False, default=3, help="The sigma of the read length distribution")
    parser.add_argument("-circular", "--circular", type=bool, required=False, default=False, help="Is the reference genome circular?")

    args = parser.parse_args()

    SIM_reads(filename=args.reference, ID=args.ID, synthetic=args.synthetic, coverage=args.coverage, reads_len=args.reads_len, sigma=args.sigma, circular=args.circular)
