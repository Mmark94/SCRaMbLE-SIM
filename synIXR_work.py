import pandas as pd
from comparison_sol import path_str
from comparison_sol import str_path
from Mapping_coverage_MM import plot_LU_CN
from Mapping_coverage_MM import plot_LU_CN_percentage


# Calculate the LU copy number of synIXR data

#import excel file with solutions
ALL_sol = pd.read_excel("synIXR_solutions_ALL2.xlsx", engine="openpyxl")
print(ALL_sol.head())

# In synIXR_solutions_ALL2.xlsx I removed strain JS603 because its has two synIXR. I also manipulated strain JS606 to fit the eccDNA inside the synIXR.

ID = list(ALL_sol["ID"])
solution = list(ALL_sol["solution"])
solution = [str_path(x) for x in solution]

syn_chr = list(range(1, 45, 1))
essential = [2, 7, 9, 10, 12, 19, 20, 24]   # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.

# Note, very important!
# LU 14 contains MET28 and LU 32 contains LYS1. After SCRaMbLE 33 strains were selected to be -LYS, 20 to be -MET and 10 to be -LYS -MET.

plot_LU_CN(solution, Plot="boxplot", essential=essential, CEN=[2])  # "histogram", "boxplot", "violinplot"
#plot_LU_CN(solution, Plot="histogram", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
#plot_LU_CN(solution, Plot="violinplot", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
plot_LU_CN_percentage(solution, max_CN=5, essential=essential, CEN=[2])
plot_LU_CN_percentage(solution, max_CN=10, essential=essential, CEN=[2])

