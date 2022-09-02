from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import pandas as pd
import random

# Cuts the reference genome in segments of the same size. If you have some mapping_breakpoints the program will be force to cut here as well.
def cut_reference_in_segments(filename="IXR_BACnewseq.fa", ID="", starting_LU=1, segment_size=1000, min_segment_size=100, mapping_breakpoints={}, RANDOM=False, random_cuts=100):
    if ID == "":
        # Take the name before the dot. find() finds the first occurrence. rfind() finds the last occurrence.
        ID = filename[:filename.rfind(".")]
    num_chrs = 0    # this is a counter for the number of chromosomes
    # Open the reference fasta file
    with open(filename) as handle:
        # Store info on the LU on lists
        LU_start_list, LU_end_list, LU_size_list, LU_ID_list, LU_DNA_list, LU_list, Chr_list = [], [], [], [], [], [], []
        for record in SeqIO.parse(handle, "fasta"):
            print(record.id)
            num_chrs = num_chrs + 1
            LU_list.append([])
            record_len = len(record.seq)
            # Use range() to define the positions were to cut. [1:] at the end removes the zero at the beginning
            pos_loxP = range(0, record_len, segment_size)[1:]
            # mapping_breakpoints = {chr_ID1: [pos1, pos2, pos3], chr_ID2: [pos1, pos2], ...}
            if mapping_breakpoints != {}:
                # Add the mapping breakpoints
                pos_loxP = list(pos_loxP) + mapping_breakpoints[record.id]
                pos_loxP.sort()
                # Remove cut points that generate segments too small
                for break_point in mapping_breakpoints[record.id]:
                    index = pos_loxP.index(break_point)
                    distance1 = pos_loxP[index] - pos_loxP[index-1]
                    distance2 = pos_loxP[index+1] - pos_loxP[index]
                    #print("distance1, distance2 =", distance1, distance2)
                    if distance2 < min_segment_size:
                        pos_loxP.pop(index + 1)
                    if distance1 < min_segment_size:
                        pos_loxP.pop(index - 1)
            if RANDOM:
                # Choose at random the points where to cut. The number of cuts is decided with the variable random_cuts.
                pos_loxP = sorted(random.sample(range(1, record_len, 1), k=random_cuts))
            print("pos_loxP =", list(pos_loxP))
            # Loop through the loxP site positions and store the LU start, end position and its size
            start_LU = 0
            skip_last_LU = False
            for i in range(len(pos_loxP)):
                if record_len - pos_loxP[i] < min_segment_size:
                    # The last LU is too small and need to be added to the previous one.
                    end_LU = record_len
                    skip_last_LU = True
                else:
                    end_LU = pos_loxP[i]
                LU = record.seq[start_LU:end_LU]
                LU_DNA_list.append(LU)     # Seq(LU) or LU ?
                LU_len = len(LU)
                # the ID of each LU is the ID of the chromosome_start_position_end_position_LUID
                LU_ID = record.id + "_" + str(start_LU+1) + "_" + str(end_LU) + "_" + str(starting_LU)
                LU_ID_list.append(LU_ID)
                LU_start_list.append(start_LU+1)
                LU_end_list.append(end_LU)
                LU_size_list.append(LU_len)
                LU_list[num_chrs-1].append(starting_LU)
                Chr_list.append(record.id)
                start_LU = pos_loxP[i]
                starting_LU = starting_LU + 1
                print(LU_ID, LU_len)
                #print(LU)
            if skip_last_LU:
                # Skip the last LU because it is too small. This was added to the previous LU
                continue
            end_LU = len(record.seq)
            last_LU = record.seq[start_LU:]
            LU_DNA_list.append(last_LU)
            LU_len = len(last_LU)
            LU_ID = record.id + "_" + str(start_LU+1) + "_" + str(end_LU) + "_" + str(starting_LU)
            LU_ID_list.append(LU_ID)
            LU_start_list.append(start_LU+1)
            LU_end_list.append(end_LU)
            LU_size_list.append(LU_len)
            LU_list[num_chrs-1].append(starting_LU)
            Chr_list.append(record.id)
            print(LU_ID, LU_len)
            #print(last_LU)
            #print(LU_start_list, LU_end_list, LU_size_list)

            # Save all the LUs into a single fasta file
            records = []
            LU_list_one_list = []
            for Chr in LU_list:
                LU_list_one_list = LU_list_one_list + Chr
            for i in range(len(LU_ID_list)):        # Guarda i+1!!!
                Description = "Chromosome=" + str(record.id) + ", LU=" + str(LU_list_one_list[i]) + ", start=" + str(LU_start_list[i]) + ", end=" + str(LU_end_list[i]) + ", size=" + str(LU_size_list[i])
                records.append(SeqRecord(LU_DNA_list[i], LU_ID_list[i], description=Description))
            SeqIO.write(records, ID + ".segments.fa", "fasta")

            # Save the info about each LU
            df_LU_position = pd.DataFrame()
            # Add the data in different columns
            df_LU_position["Chromosome"] = Chr_list
            df_LU_position["LU"] = LU_list_one_list
            df_LU_position["LU start"] = LU_start_list
            df_LU_position["LU end"] = LU_end_list
            df_LU_position["LU size"] = LU_size_list
            df_LU_position["LU DNA"] = LU_DNA_list
            # Save the DataFrame in an excel file
            df_LU_position.to_csv(ID + "_segments_info.csv")
    return LU_list

if __name__ == '__main__':
    #cut_reference_in_segments(filename=str(sys.argv[1]))
    #genome = cut_reference_in_segments(filename="IXR_BACnewseq.fa", starting_LU=1, segment_size=1000, min_segment_size=300, mapping_breakpoints={"IXR_BAC_LU44_del": [2500, 3900]})
    #genome = cut_reference_in_segments(filename="BY4741_chr09.fa", starting_LU=1, segment_size=1000, min_segment_size=300)
    #print(genome)
    # IXR_BACnewseq_test.fa has two synIXR chromosomes
    #cut_reference_in_segments(filename="IXR_BACnewseq_test.fa", segment_size=3000, min_segment_size=500)
    #cut_reference_in_segments(filename="SynIII-KC880027.1.fasta")

    parser = argparse.ArgumentParser(description="Cut the reference chromosomes in segments of a define length.")

    parser.add_argument("-filename", "--filename", type=str, required=True, help="The reference fasta file that you want to cut in pieces")
    parser.add_argument("-ID", "--ID", type=str, required=False, default="", help="The name of the output file")
    parser.add_argument("-segment_size", "--segment_size", type=int, required=False, default=1000, help="The size of each segment. Default=1000")
    parser.add_argument("-min_segment_size", "--min_segment_size", type=int, required=False, default=100, help="The minimal size of a segment. Default=100")

    args = parser.parse_args()

    cut_reference_in_segments(filename=args.filename, ID=args.ID, starting_LU=1, segment_size=args.segment_size, min_segment_size=args.min_segment_size)
