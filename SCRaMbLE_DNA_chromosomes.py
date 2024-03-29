from get_segments import cut_reference_in_segments
from get_loxp_segment import cut_reference_in_LUs
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from SCRaMbLE_simulation_3 import SCRaMbLE_muliple_chrs
from subpath_to_DNA1 import path_to_DNA
import argparse
import random


# SCRaMbLE chromosomes sequences as DNA and output SCRaMbLEd chromosomes DNA.
def SCRaMbLE_DNA(filename="IXR_BACnewseq.fa", ID="", synthetic=False, Number_SCRaMbLE_events=15, starting_LU=1, segment_size=1000, min_segment_size=100,
                 mapping_breakpoints={}, RANDOM=False, random_cuts=100, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], Ptra=0.05):
    if ID == "":
        # Take the name before the dot. find() finds the first occurrence. rfind() finds the last occurrence.
        ID = filename[:filename.rfind(".")]
    # Cut the DNA in segments or LU and extract the LUs sequence.
    if synthetic:
        genome = cut_reference_in_LUs(filename=filename, ID=ID, starting_LU=starting_LU)
    else:
        genome = cut_reference_in_segments(filename=filename, ID=ID, starting_LU=starting_LU, segment_size=segment_size, min_segment_size=min_segment_size, mapping_breakpoints=mapping_breakpoints, RANDOM=RANDOM, random_cuts=random_cuts)
    print("genome =", genome)
    if isinstance(genome[0], list) and len(genome) == 1:
        # There is only one chromosome
        genome = genome[0]
    # SCRaMbLE the chromosomes paths
    SCRaMbLEd_chrs = SCRaMbLE_muliple_chrs(list_chr=genome, Number_events=Number_SCRaMbLE_events, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability, Ptra=Ptra)
    print("SCRaMbLEd_chrs =", SCRaMbLEd_chrs)
    # Convert the SCRaMbLEd paths into DNA
    if synthetic:
        LU_fasta_ID = ID + ".loxpreg.fa"
    else:
        LU_fasta_ID = ID + ".segments.fa"
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    path_to_DNA(path=SCRaMbLEd_chrs, LU_fasta=LU_fasta_ID, filename=ID+"_SCRaMbLEd_SE"+str(Number_SCRaMbLE_events) + "_R" + str(random_seed), ID=ID+"_SCRaMbLEd")
    return None

# test the code
if __name__ == '__main__':
    #test the script with the synIXR chromosome
    #SCRaMbLE_DNA(filename="IXR_BACnewseq.fa", synthetic=True, Number_SCRaMbLE_events=20)
    #SCRaMbLE_DNA(filename="BY4741_chr09.fa", ID="chr09", synthetic=False, Number_SCRaMbLE_events=20)
    #SCRaMbLE_DNA(filename="BY4741_chr02_03.fasta", ID="chr02_03", synthetic=False, Number_SCRaMbLE_events=20)

    parser = argparse.ArgumentParser(description="SCRaMbLE chromosomes sequences and output SCRaMbLEd chromosomes fasta files.")

    parser.add_argument("-filename", "--filename", type=str, required=True, help="The reference fasta file that you want to cut in pieces")
    parser.add_argument("-ID", "--ID", type=str, required=False, default="", help="The name of the output file. If empty, the script will take the reference's filename before the dot (.)")
    parser.add_argument("-synthetic", "--synthetic", type=bool, required=False, default=False, help="Is the reference genome synthetic? Does it contain loxPsym sites?")
    parser.add_argument("-num_SE", "--Number_SCRaMbLE_events", type=int, required=True, default=15, help="How many SCRaMbLE events do you want to simulate?")
    parser.add_argument("-segment_size", "--segment_size", type=int, required=False, default=1000, help="The size of each segment. Default=1000")
    parser.add_argument("-min_segment_size", "--min_segment_size", type=int, required=False, default=100, help="The minimal size of a segment. Default=100")
    parser.add_argument("-circular", "--circular", type=bool, required=False, default=False, help="Use this flag if the chromosome is circular. Not having this flag will evaluate to False")

    args = parser.parse_args()
    # Note: --synthetic False will evaluate to True! removing --synthetic will evaluate to False

    SCRaMbLE_DNA(filename=args.filename, ID=args.ID, synthetic=args.synthetic, Number_SCRaMbLE_events=args.Number_SCRaMbLE_events , starting_LU=1, segment_size=args.segment_size, min_segment_size=args.min_segment_size, circular=args.circular)
