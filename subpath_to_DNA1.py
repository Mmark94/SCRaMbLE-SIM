from comparison_sol import str_path
from comparison_sol import path_str
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import pandas as pd

# path_str is from comparison_sol.py
def path_str(x: list):
    to_list = [str(i) for i in x]
    to_string = ",".join(to_list)
    return to_string

def extract_LU_len(LU_fasta="IXR_BACnewseq.loxpreg.fa", SEQ=False):
    LU_ID_len_seq = {}
    # Open the reference fasta file
    with open(LU_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Extract the LU integer from the LU_ID. Usually the LU_ID is chrname_start_end_LU. Therefore we need the part after the last "_". This can be access by splitting the LU_ID into a list and by taking the  last element [-1].
            LU_ID = int(record.id.split("_")[-1])
            #print(record.id, LU_ID)
            LU_seq = record.seq
            if SEQ:
                LU_ID_len_seq[LU_ID] = LU_seq
            else:
                # Only save the LU length
                LU_ID_len_seq[LU_ID] = len(LU_seq)
    return LU_ID_len_seq

# If the input is a list of paths, it will save each path in a different fasta file
def path_to_DNA_separate(path=[1,2,3], LU_fasta="IXR_BACnewseq.loxpreg.fa", filename="A_test"):
    if isinstance(path, str):
        # For example:  "1,2,-3" => [1, 2, -3]
        path = str_path(path)
    # Extract the LU sequences. LU_ID_seq = {LU_ID: LU_seq}
    LU_ID_seq = extract_LU_len(LU_fasta=LU_fasta, SEQ=True)
    if isinstance(path[0], list):   # the input is a list of subpaths
        for i in range(len(path)):
            path_to_DNA(path=path[i], LU_fasta="IXR_BACnewseq.loxpreg.fa", filename=filename + "_" + str(i))
        return None
    else:
        # Convert the path into a DNA sequence
        path_DNA = Seq("")
        for LU in path:
            if abs(LU) not in LU_ID_seq.keys():     # Note this assumes that all LUs in LU_ID_seq.keys() are positive
                print("ERROR! The LU", abs(LU), "in the path is not present in the fasta file", LU_fasta)
                print("Check your path or your fasta file!")
                #continue
                return None
            if LU >= 0:
                path_DNA = path_DNA + LU_ID_seq[LU]
            else:
                # The LU is negative. Therefore it is in the reverse orientation.
                path_DNA = path_DNA + LU_ID_seq[abs(LU)].reverse_complement()

        # Save the path DNA into a fasta file
        if isinstance(path[0], int):  # There is only one sequence in "subpaths"
            path = [path]
            path_DNA = [path_DNA]
        records = []
        for i in range(len(path)):
            path_string = path_str(path[i])
            Description = "path=" + path_string + ", path_DNA_size=" + str(len(path_DNA[i]))
            records.append(SeqRecord(path_DNA[i], "path_" + path_string, description=Description))
        SeqIO.write(records, filename + ".sol.fa", "fasta")
        return path_DNA

# Use this function to convert paths or subpaths into DNA sequences in fasta files.
# If the input is a list of paths, it will save each path in same fasta file in different lines
def path_to_DNA(path=[1,2,3], LU_fasta="IXR_BACnewseq.loxpreg.fa", filename="A_test", ID="path"):
    if isinstance(path, str):
        # For example:  "1,2,-3" => [1, 2, -3]
        path = str_path(path)
    # Extract the LU sequences. LU_ID_seq = {LU_ID: LU_seq}
    LU_ID_seq = extract_LU_len(LU_fasta=LU_fasta, SEQ=True)

    if isinstance(path[0], list):   # The input is a list of subpaths
        # Convert the paths into DNA sequences
        path_DNA = []
        for i in range(len(path)):
            path_DNA_temp = Seq("")
            for LU in path[i]:
                if abs(LU) not in LU_ID_seq.keys():     # Note this assumes that all LUs in LU_ID_seq.keys() are positive
                    print("ERROR! The LU", abs(LU), "in the path is not present in the fasta file", LU_fasta)
                    print("Check your path or your fasta file!")
                    # continue
                    return None
                if LU >= 0:
                    path_DNA_temp = path_DNA_temp + LU_ID_seq[LU]
                else:
                    # The LU is negative. Therefore it is in the reverse orientation.
                    path_DNA_temp = path_DNA_temp + LU_ID_seq[abs(LU)].reverse_complement()
            path_DNA.append(path_DNA_temp)
    else:
        # Convert the path into a DNA sequence
        path_DNA = Seq("")
        for LU in path:
            if abs(LU) not in LU_ID_seq.keys():     # Note this assumes that all LUs in LU_ID_seq.keys() are positive
                print("ERROR! The LU", abs(LU), "in the path is not present in the fasta file", LU_fasta)
                print("Check your path or your fasta file!")
                #continue
                return None
            if LU >= 0:
                path_DNA = path_DNA + LU_ID_seq[LU]
            else:
                # The LU is negative. Therefore it is in the reverse orientation.
                path_DNA = path_DNA + LU_ID_seq[abs(LU)].reverse_complement()

    # Save the path DNA into a fasta file
    if isinstance(path[0], int):  # There is only one sequence in "subpaths"
        path = [path]
        path_DNA = [path_DNA]
    records = []
    for i in range(len(path)):
        path_string = path_str(path[i])
        Description = "path=" + path_string + ", path_DNA_size=" + str(len(path_DNA[i]))
        records.append(SeqRecord(path_DNA[i], ID + "_" + str(i+1), description=Description))   # + "_" + path_string
    SeqIO.write(records, filename + ".sol.fa", "fasta")
    #return path_DNA
    return None

if __name__ == '__main__':
    LU_fasta = "IXR_BACnewseq.loxpreg.fa"
    filename = "A_test"
    path = [1, 2, 3, -4]
    path = [[3, -4], [4,-5,-6], [7,-8,-9]]
    #print(extract_LU_len(LU_fasta=LU_fasta))
    #path_to_DNA(path=path, LU_fasta=LU_fasta, filename=filename, ID="path")


    parser = argparse.ArgumentParser(description="Convert the path or subpath into DNA sequences.")

    parser.add_argument("-path", "--path", type=str, metavar="", required=True, help="The path/s that you want to convert in DNA sequences")
    parser.add_argument("-LU", "--LU_fasta", type=str, metavar="", required=True, help="The fasta file that contain the sequences of the LUs")
    parser.add_argument("-name", "--filename", type=str, metavar="", required=False, default="A_test", help="The name of the output file")
    parser.add_argument("-ID", "--seq_ID", type=str, metavar="", required=False, default="path", help="The ID of the sequence inside the fasta file")

    args = parser.parse_args()

    path_to_DNA(path=args.path, LU_fasta=args.LU_fasta, filename=args.filename, ID=args.seq_ID)
