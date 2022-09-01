from get_segments import cut_reference_in_segments
from get_loxp_segment import cut_reference_in_LUs
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from SCRaMbLE_simulation_3 import SCRaMbLE_muliple_chrs
from subpath_to_DNA1 import path_to_DNA
import argparse


# SCRaMbLE chromosomes sequences as DNA and output SCRaMbLEd chromosomes DNA.
# This is not finished
def SCRaMbLE_DNA(filename="IXR_BACnewseq.fa", ID="", synthetic=True, starting_LU=1,segment_size=1000, min_segment_size=100,
                 mapping_breakpoints={}, RANDOM=False, random_cuts=100, Number_SCRaMbLE_events=15):
    if ID == "":
        # Take the name before the dot. find() finds the first occurrence. rfind() finds the last occurrence.
        ID = filename[:filename.rfind(".")]
    # Cut the DNA in segments or LU
    if synthetic:
        cut_reference_in_LUs(filename=filename, ID=ID, starting_LU=starting_LU)
    else:
        cut_reference_in_segments(filename=filename, ID=ID, starting_LU=starting_LU, segment_size=segment_size, min_segment_size=min_segment_size, mapping_breakpoints=mapping_breakpoints, RANDOM=RANDOM, random_cuts=random_cuts, essential=[])

    # Convert the genome in a sequence of segments
    genome = []

    # SCRaMbLE the chromosomes paths
    SCRaMbLEd_chrs = SCRaMbLE_muliple_chrs(list_chr=genome, Number_events=Number_SCRaMbLE_events, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], Ptra=0.05)

    # Convert the SCRaMbLEd paths into DNA
    path_to_DNA(path=SCRaMbLEd_chrs, LU_fasta=ID+".loxpreg.fa", filename=filename+"_SCRaMbLEd_"+str(Number_SCRaMbLE_events)+".fasta", ID=ID+"_SCRaMbLEd")
    return None

# test the code
if __name__ == '__main__':

    #"""
    parser = argparse.ArgumentParser(description="SCRaMbLE chromosomes sequences and output SCRaMbLEd chromosomes fasta files.")
    # metavar=""
    parser.add_argument("-filename", "--filename", type=str, required=True, help="The reference fasta file that you want to cut in pieces")
    parser.add_argument("-ID", "--ID", type=str, required=False, default="", help="The name of the output file")
    parser.add_argument("-segment_size", "--segment_size", type=int, required=False, default=1000, help="The size of each segment. Default=1000")
    parser.add_argument("-min_segment_size", "--min_segment_size", type=int, required=False, default=100, help="The minimal size of a segment. Default=100")

    args = parser.parse_args()

    SCRaMbLE_DNA(filename=args.filename, ID=args.ID, starting_LU=1, segment_size=args.segment_size, min_segment_size=args.min_segment_size)
    #"""
