import pandas as pd
from comparison_sol import path_str
from comparison_sol import str_path
from Mapping_coverage_MM import plot_LU_CN
from Mapping_coverage_MM import plot_LU_CN_percentage


# Calculate the LU copy number of synIXR data

#import excel file with solutions
ALL_sol = pd.read_excel("synIXR_solutions_ALL2.xlsx", engine="openpyxl")
print(ALL_sol.head())

# In synIXR_solutions_ALL2.xlsx I removed strain JS603 because it has two synIXR. I also manipulated strain JS606 to fit the eccDNA inside the synIXR.

ID = list(ALL_sol["ID"])
solution = list(ALL_sol["solution"])
solution = [str_path(x) for x in solution]

syn_chr = list(range(1, 45, 1))
essential = [2, 7, 9, 10, 12, 19, 20, 24]   # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.

# Note, very important!
# LU 14 contains MET28 and LU 32 contains LYS1. After SCRaMbLE 33 strains were selected to be -LYS, 20 to be -MET and 10 to be -LYS -MET.

#plot_LU_CN(solution, Plot="boxplot", essential=essential, CEN=[2])  # "histogram", "boxplot", "violinplot"
#plot_LU_CN(solution, Plot="histogram", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
#plot_LU_CN(solution, Plot="violinplot", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
plot_LU_CN_percentage(solution, max_CN=5, essential=essential, CEN=[2], filename="synIXR", SE="")    # SE="6.2"
#plot_LU_CN_percentage(solution, max_CN=10, essential=essential, CEN=[2])

import matplotlib.pyplot as plt
from Mapping_coverage_MM import count_essential_LU_CN
# Plot the essential and non-essential LU CN
MAX_LU = max([abs(x) for x in syn_chr])
essential_no_centromere = essential[:]
essential_no_centromere.remove(2)
essential_LUs, non_essential_LUs = count_essential_LU_CN(solution, essential=essential_no_centromere, Max_LU=MAX_LU)
print("essential_LUs     =", essential_LUs)
print("non_essential_LUs =", non_essential_LUs)
plt.figure(figsize=(7, 3.5), dpi=200)
plt.ylabel("LU CN")
plt.title("essential vs non-essential LU CN")
Labels = ["essential LU CN", "non-essential LU CN"]
plt.axhline(y=0, color="grey", linestyle="-", alpha=0.3)
plt.boxplot([essential_LUs, non_essential_LUs], labels=Labels, showfliers=True)
# plt.violinplot([essential_LUs, non_essential_LUs])
plt.show()
plt.close()


# Statistics
from comparison_sol import NG50_calculator
import statistics
import numpy as np

# Statistics on the chromosomes
chr_len_L = []
len_unique_L = []
NG50_L = []
LG50_L = []
duplication_rate = []
for CHR in solution:
    #if isinstance(CHR[0], list):
    #    CHR = CHR[0]
    CHR_stats = NG50_calculator(CHR)
    chr_len_L.append(CHR_stats[0])
    len_unique_L.append(CHR_stats[1])
    NG50_L.append(CHR_stats[2])
    LG50_L.append(CHR_stats[3])
    duplication_rate.append(CHR_stats[5])
print()
print("chr_len_L =", statistics.mean(chr_len_L), "+-", statistics.stdev(chr_len_L), ". Percentiles =", np.percentile(chr_len_L, [25, 50, 75]))
print("len_unique_L =", statistics.mean(len_unique_L), "+-", statistics.stdev(len_unique_L), ". Percentiles =", np.percentile(len_unique_L, [25, 50, 75]))
print("NG50_L =", statistics.mean(NG50_L), "+-", statistics.stdev(NG50_L), ". Percentiles =", np.percentile(NG50_L, [25, 50, 75]))
print("LG50_L =", statistics.mean(LG50_L), "+-", statistics.stdev(LG50_L), ". Percentiles =", np.percentile(LG50_L, [25, 50, 75]))
print("duplication_rate =", statistics.mean(duplication_rate), "+-", statistics.stdev(duplication_rate), ". Percentiles =", np.percentile(duplication_rate, [25, 50, 75]))

"""
ALL synIXR
chr_len_L = 49.17460317460318 +- 28.765761279990645 . Percentiles = [36.  40.  50.5]
len_unique_L = 35.698412698412696 +- 4.726844820381264 . Percentiles = [34. 37. 39.]
NG50_L = 1.5873015873015872 +- 1.4439520213068542 . Percentiles = [1. 1. 2.]
LG50_L = 15.206349206349206 +- 3.4085221949559643 . Percentiles = [13. 17. 17.]
duplication_rate = 1.4019840821385428 +- 0.9355965777733104 . Percentiles = [1.         1.         1.51923077]

Only synIXR unsolved by short reads
chr_len_L = 64.48148148148148 +- 38.94312457375437 . Percentiles = [41.5 54.  70. ]
len_unique_L = 34.592592592592595 +- 5.930749461390961 . Percentiles = [31.5 36.  39. ]
NG50_L = 2.3333333333333335 +- 1.9806758753205742 . Percentiles = [1. 2. 2.]
LG50_L = 12.925925925925926 +- 3.6047610176406972 . Percentiles = [10. 13. 16.]
duplication_rate = 1.8922559183954804 +- 1.2766803020611441 . Percentiles = [1.21111111 1.55882353 1.90541872]
"""
