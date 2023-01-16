from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import pandas as pd

# This script divides a reference genome (fasta) into a single fasta file containing the LoxPsym Units (LUs) sequences.
# Basically it cuts the reference genome where it finds a loxPsym site. Therefore, at each LU is assigned an identifier integer.

def cut_reference_in_LUs(filename="IXR_BACnewseq.fa", ID="", starting_LU=1):
    if ID == "":
        # Take the name before the dot. find() finds the first occurrence. rfind() finds the last occurrence.
        ID = filename[:filename.rfind(".")]
    # Define the sequence used for the cutting of the reference (restriction site)
    loxpseq = Seq("ATAACTTCGTATAATGTACATTATACGAAGTTAT")
    #print(loxpseq, len(loxpseq))
    #print(loxpseq[0:17], loxpseq[17:34])
    num_chrs = 0    # this is a counter for the number of chromosomes
    # Open the reference fasta file
    with open(filename) as handle:
        # Store info on the LU on lists
        LU_start_list, LU_end_list, LU_size_list, LU_ID_list, LU_DNA_list, LU_list, Chr_list = [], [], [], [], [], [], []
        for record in SeqIO.parse(handle, "fasta"):
            print(record.id)
            num_chrs = num_chrs + 1
            LU_list.append([])
            # Search the positions of loxP sites in the reference and store this in a list
            pos_loxP = SeqUtils.nt_search(str(record.seq), loxpseq)[1:]
            #print("pos_loxP =", pos_loxP)
            # Loop through the loxP site positions and store the LU start, end position and its size
            start_LU = 0
            for i in range(len(pos_loxP)):
                end_LU = pos_loxP[i] + 17
                #print(record.seq[pos_loxP[i]:pos_loxP[i]+34])   # this will print the loxP site
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
                start_LU = pos_loxP[i] + 17
                starting_LU = starting_LU + 1
                print(LU_ID, LU_len)
                #print(LU)
            # add last LU. This will be from the position of the last loxP site until the end of the chromosome
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
        for i in range(len(LU_ID_list)):        # Check i+1!!!
            Description = "Chromosome=" + str(record.id) + ", LU=" + str(LU_list_one_list[i]) + ", start=" + str(LU_start_list[i]) + ", end=" + str(LU_end_list[i]) + ", size=" + str(LU_size_list[i])
            records.append(SeqRecord(LU_DNA_list[i], LU_ID_list[i], description=Description))
        SeqIO.write(records, ID + ".loxpreg.fa", "fasta")

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
        df_LU_position.to_csv(ID + "_LU_info.csv")
    return LU_list

if __name__ == '__main__':
    #cut_reference_in_LUs(filename=str(sys.argv[1]))
    #extracted_LUs_seq = cut_reference_in_LUs(filename="IXR_BACnewseq.fa")
    #extracted_LUs_seq = cut_reference_in_LUs(filename="synII_III.fasta")
    #print(extracted_LUs_seq)
    # IXR_BACnewseq_test.fa has two synIXR chromosomes
    #cut_reference_in_LUs(filename="IXR_BACnewseq_test.fa")
    #cut_reference_in_LUs(filename="SynIII-KC880027.1.fasta")

    parser = argparse.ArgumentParser(description="Divide the reference synthetic chromosome fasta file in a single fasta file containing the LU sequences.")

    parser.add_argument("-ref", "--reference", type=str, required=True, help="The reference fasta file containing a synthetic chromosome")
    parser.add_argument("-ID", "--filename", type=str, required=False, default="", help="The name of the output file")

    args = parser.parse_args()

    cut_reference_in_LUs(filename=args.reference, ID=args.filename, starting_LU=1)
