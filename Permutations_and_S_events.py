import itertools
from comparison_sol import invert
from comparison_sol import half_pos2
from comparison_sol import count_solutions
from comparison_sol import str_path
from comparison_sol import remove_duplicate
from comparison_sol import descescent_LU_num
import pandas as pd


# how long are the chromosomes?
chr_L_range = range(1, 7)

# how many different LUs are there?
num_LU_range = range(1, 7)

# Number of permutations and SCRaMbLE events (SE)
num_Permutations = []
num_SE = []

Circular = False
Brute_Force = False

# Plot all the data. Chromosome length, number of solutions and time for the assembly
ALL_chr = []         # Store the path
ALL_chr_len = []     # Store the path length

for chr_L in chr_L_range:
    for num_LU in num_LU_range:
        print(chr_L, num_LU)
        if chr_L >= num_LU:
            # Calculate the number of permutations
            if chr_L % 2 != 0:  # The chromosome has odd length
                number_permutations = (2 ** (chr_L - 1)) * (num_LU ** chr_L)
            else:
                # For even length chromosomes I need to add the number of palindromic paths. These are the permutations of the half chromosome length
                number_permutations = (2 ** (chr_L - 1)) * (num_LU ** chr_L) + 2 ** ((chr_L / 2) - 1) * num_LU ** (chr_L / 2)

            # Calculate the number of SCRaMbLE events
            events = (3/2) * (chr_L + 1) * chr_L - 1   # -1 is the identity element and comes from the inversion of the whole chromosome

            if Circular:
                number_permutations = number_permutations / chr_L
                events = (3 / 2) * chr_L * (chr_L - 1) - 1
            print("number_permutations =", number_permutations)
            print("events =", events)
            if chr_L == num_LU:     # Only take permutations and SE if the chromosome length is equal to the number of num_LU
                num_Permutations.append(int(number_permutations))
                num_SE.append(int(events))

            if Brute_Force == True:
                # we create all the possible combinations of chromosomes with different LUs
                all_LUs = range(1, num_LU + 1)
                all_LUs_neg = [-i for i in all_LUs]
                all_LUs = list(all_LUs) + list(all_LUs_neg)
                #print("all_LUs =", list(all_LUs))
                CHR = list(itertools.combinations_with_replacement(all_LUs, chr_L))
                # convert tuples in lists
                CHR = [list(i) for i in CHR]
                #print("CHR =", CHR)
                # now with the chromosome with the LUs defined we need to find all the possible permutations of the LUs
                Permutations = []
                for i in CHR:
                    Permutations = Permutations + list(itertools.permutations(i))
                # convert tuples in lists
                Permutations = [list(i) for i in Permutations]
                # now the duplicated chromosomes need to be removed
                Permutations = remove_duplicate(Permutations, list_of_list=True)
                Permutations.sort()
                L_permutations = len(Permutations)
                #print("Permutations =", Permutations)
                print("How many?           =", L_permutations)
        print("--------------------------------------------")

print("num_Permutations =", num_Permutations)
print("num_SE =", num_SE)

"""
ALL_chr_len = [len(x) for x in ALL_chr]    # Store the path length
# Put all the data from lists in a DataFrame
df = pd.DataFrame(columns=["path", "path_len"])
df["path"] = ALL_chr
df["path_len"] = ALL_chr_len

# Save the DataFrame in an excel file
max_chr_L = str(max(ALL_chr_len))
df.to_excel("Permutations_chr_L" + max_chr_L + ".xlsx")


Retrieve = False
# Retrieve the data from excel
if Retrieve:
    df_permutations = pd.read_excel("Permutations_chr_L5.xlsx")
    print(df_permutations.head())
    ALL_chr = list(df_permutations["path"])
    ALL_chr_len = list(df_permutations["path_len"])
"""

"""
# Store the data
# Linear chromosomes
num_Permutations =  [1, 10, 108, 2080, 50000, 1493856, 52706752, 2147516416, 99179645184, 5120001600000]
num_SE =            [2, 8, 17, 29, 44, 62, 83, 107, 134, 164]
# Circular chromosomes
num_Permutations =  [1, 5, 36, 520, 10000, 248976, 7529536, 268439552, 11019960576, 512000160000]
num_SE =            [-1, 2, 8, 17, 29, 44, 62, 83, 107, 134]
"""
