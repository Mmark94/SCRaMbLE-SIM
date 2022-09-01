import random
import matplotlib.pyplot as plt
from SCRaMbLE_simulation_3_circular import SCRaMbLE4_circular
from Mapping_coverage_MM import plot_LU_CN
from Mapping_coverage_MM import plot_LU_CN_percentage
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from SCRaMbLE_simulation_3 import SCRaMbLE4_lin_cir


# test the code
if __name__ == "__main__":
    segments = 44  #number of loxP segments
    syn_chr = list(range(1, segments+1, 1))
    essential = [2, 7, 9, 10, 12, 19, 20, 24]  # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.
    print("syn_chr =", syn_chr)


    #for i in range(10):
    #    print(SCRaMbLE4_lin_cir(syn_chr, 10, essential, False, probability=[2,1,1,2]))
    #    print(SCRaMbLE4_lin_cir(syn_chr, 10, essential, True, probability=[2,1,1,2]))
    #    print(force_SCRaMLE_lin_cir(syn_chr, 10, essential, False, probability=[2,1,1,2]))
    #    print()

    CHR_list = []
    for i in range(1000):
        print(i)
        CHR_list.append(force_SCRaMLE_lin_cir(syn_chr, 50, essential, circular=True, CEN=[2], probability=[0,2,2,1]))
    #plot_LU_CN(CHR_list, Plot="boxplot", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
    #plot_LU_CN(CHR_list, Plot="histogram", essential=essential, CEN=[2])      #"histogram", "boxplot", "violinplot"
    plot_LU_CN_percentage(CHR_list, max_CN=5, essential=essential, CEN=[2])
    # I use the following website to create the giff: https://ezgif.com/maker