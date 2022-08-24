from get_segments import cut_reference_in_segments
from get_loxp_segment import cut_reference_in_LUs
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from SCRaMbLE_simulation_3 import SCRaMbLE_muliple_chrs
from subpath_to_DNA1 import path_to_DNA
import argparse


# SCRaMbLE chromosomes sequences as DNA and output SCRaMbLEd chromosomes DNA.
# This is not finished
def SCRaMbLE_DNA(filename="IXR_BACnewseq.fa", ID="", starting_LU=1, segment_size=1000, min_segment_size=100, mapping_breakpoints={}, RANDOM=False, random_cuts=100):
    # Cut the DNA in segments
    cut_reference_in_segments(filename=filename, ID=ID, starting_LU=starting_LU, segment_size=segment_size, min_segment_size=min_segment_size, mapping_breakpoints=mapping_breakpoints, RANDOM=RANDOM, random_cuts=random_cuts)
    #cut_reference_in_LUs(filename=filename, ID=ID, starting_LU=starting_LU)
    # SCRaMbLE the path
    #SCRaMbLEd_chrs = SCRaMbLE_muliple_chrs(list_chr=[1,2,3], Number_events=1, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], Ptra=0.05)

    # Convert the SCRaMbLEd paths into DNA
    #path_to_DNA(path=[1, 2, 3], LU_fasta="IXR_BACnewseq.loxpreg.fa", filename="A_test", ID="path")
    return None


if __name__ == '__main__':

    #"""
    parser = argparse.ArgumentParser(description="Cut the reference chromosomes in segments of a define length.")
    # metavar=""
    parser.add_argument("-filename", "--filename", type=str, required=True, help="The reference fasta file that you want to cut in pieces")
    parser.add_argument("-ID", "--ID", type=str, required=False, default="", help="The name of the output file")
    parser.add_argument("-segment_size", "--segment_size", type=int, required=False, default=1000, help="The size of each segment. Default=1000")
    parser.add_argument("-min_segment_size", "--min_segment_size", type=int, required=False, default=100, help="The minimal size of a segment. Default=100")

    args = parser.parse_args()

    SCRaMbLE_DNA(filename=args.filename, ID=args.ID, starting_LU=1, segment_size=args.segment_size, min_segment_size=args.min_segment_size)
    #"""
