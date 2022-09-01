from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from SCRaMbLE_simulation_3 import SCRaMbLE4
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from comparison_sol import NG50_calculator
from work_parallel_functions import sum_within_list
from work_parallel_functions import standard_deviation_from_points_parallel
from work_parallel_functions import divide_by_sim
from work_parallel_functions import sum_two_lists
from SCRaMbLE_simulation_chr_len_events_store2 import plot_chr_len
import matplotlib.pyplot as plt
import random

simulations = 20
Chr_len0 = 51
max_SCRaMbLEd_events = 1001
SCRaMbLEd_events_range = 10
#num_essential = round(len(syn_chr) / 6)
num_essential = 10

syn_chr = list(range(1, Chr_len0, 1))
SCRaMbLEd_events = list(range(0, max_SCRaMbLEd_events, SCRaMbLEd_events_range))

#essential = [2,7,9,10,12,20]
#essential = sorted(random.sample(syn_chr, k=num_essential))

def simulation_chr_len(simulations):
    print(simulations)
    Chr_len0 = 51
    max_SCRaMbLEd_events = 1001
    SCRaMbLEd_events_range = 10
    # num_essential = round(len(syn_chr) / 6)
    num_essential = 10

    syn_chr = list(range(1, Chr_len0, 1))
    SCRaMbLEd_events = list(range(0, max_SCRaMbLEd_events, SCRaMbLEd_events_range))

    essential = sorted(random.sample(syn_chr, k=num_essential))

    # store the chr length in a list
    chr_len = [0 for i in range(len(SCRaMbLEd_events))]
    L_unique = [0 for i in range(len(SCRaMbLEd_events))]
    LG50 = [0 for i in range(len(SCRaMbLEd_events))]

    for events in range(len(SCRaMbLEd_events)):
        S = force_SCRaMLE_lin_cir(syn_chr=syn_chr, Number_events=SCRaMbLEd_events[events], essential=essential, CEN=[2], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])
        NG50 = NG50_calculator(S)
        L_unique[events] = L_unique[events] + NG50[1]
        LG50[events] = LG50[events] + NG50[3]
        chr_len[events] = chr_len[events] + len(S)
    return chr_len, L_unique, LG50

chr_len_L = []
L_unique_L = []
LG50_L = []

if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=6) as executor:
        for result in executor.map(simulation_chr_len, range(0, simulations)):
            chr_len_L.append(result[0])
            L_unique_L.append(result[1])
            LG50_L.append(result[2])
        # Use this if you want to plot individual chr len
        #plot_chr_len(SCRaMbLEd_events=SCRaMbLEd_events, chr_len=chr_len_L, L_unique=L_unique_L,num_essential=num_essential, simulations=1)

        # calculate standard deviation
        chr_len_dv = standard_deviation_from_points_parallel(chr_len_L)
        L_unique_dv = standard_deviation_from_points_parallel(L_unique_L)
        LG50_dv = standard_deviation_from_points_parallel(LG50_L)

        # Put together all the simulations
        chr_len_L = sum_within_list(chr_len_L)
        L_unique_L = sum_within_list(L_unique_L)
        LG50_L = sum_within_list(LG50_L)

        # This is if you have results from previous simulations
        previous_simulations = 0
        chr_len = [0 for i in range(len(SCRaMbLEd_events))]
        L_unique = [0 for i in range(len(SCRaMbLEd_events))]
        LG50 = [0 for i in range(len(SCRaMbLEd_events))]
        """
        previous_simulations = 100
        chr_len = [5000, 5925, 6566, 6753, 7017, 6907, 7531, 7555, 7446, 7260, 7237, 7774, 7030, 7369, 8042, 7178, 7728, 7579, 7415, 8092, 7224, 7340, 6845, 6770, 7222, 6907, 6840, 6948, 7382, 7220, 6379, 6751, 6668, 7279, 7129, 7043, 7106, 6257, 6902, 6856, 6763, 6281, 6470, 6469, 6426, 6193, 6565, 6120, 5927, 6290, 6465, 6312, 5945, 6831, 6040, 6055, 5429, 5863, 6213, 6078, 6189, 6147, 5514, 5825, 6144, 6191, 5927, 6030, 5201, 5960, 5673, 5643, 5263, 5627, 5184, 5247, 5597, 5561, 5280, 5420, 5253, 5721, 5417, 5484, 5502, 5193, 5626, 5259, 5469, 5338, 5806, 5002, 5663, 5734, 5398, 5456, 5367, 5268, 5484, 5395, 5448]
        L_unique = [5000, 4585, 4274, 4051, 3799, 3580, 3450, 3277, 3207, 3099, 2860, 2813, 2770, 2597, 2593, 2476, 2473, 2377, 2347, 2373, 2231, 2185, 2106, 1986, 2020, 1888, 1920, 1863, 1849, 1852, 1724, 1789, 1705, 1733, 1738, 1696, 1678, 1590, 1596, 1518, 1559, 1551, 1555, 1548, 1460, 1483, 1445, 1455, 1433, 1410, 1330, 1377, 1356, 1350, 1346, 1285, 1319, 1300, 1325, 1305, 1286, 1273, 1281, 1283, 1278, 1279, 1239, 1219, 1230, 1186, 1250, 1222, 1189, 1169, 1201, 1206, 1174, 1189, 1159, 1149, 1157, 1156, 1121, 1164, 1149, 1114, 1134, 1131, 1133, 1108, 1146, 1104, 1114, 1147, 1137, 1116, 1088, 1095, 1105, 1107, 1088]
        LG50 = [2400, 1668, 1345, 1171, 1038, 980, 883, 801, 772, 736, 650, 628, 630, 575, 536, 516, 514, 507, 467, 488, 442, 446, 419, 388, 384, 357, 379, 360, 342, 341, 305, 331, 295, 315, 310, 312, 294, 285, 274, 262, 270, 265, 259, 260, 245, 239, 223, 249, 245, 238, 206, 225, 214, 211, 215, 195, 214, 198, 217, 199, 205, 200, 206, 210, 190, 203, 182, 188, 184, 180, 190, 182, 187, 174, 194, 198, 183, 185, 176, 179, 177, 175, 155, 182, 176, 171, 172, 168, 164, 154, 176, 167, 178, 155, 174, 167, 172, 161, 161, 157, 170]
        """
        # Sum the previous simulations to the new
        simulations = simulations + previous_simulations
        chr_len_L = sum_two_lists(chr_len_L,chr_len)
        L_unique_L = sum_two_lists(L_unique_L,L_unique)
        LG50_L = sum_two_lists(LG50_L,LG50)

        print("----------")
        print("syn_chr =", len(syn_chr))
        print("num_essential = ", num_essential)
        print("simulation =", simulations)
        print("chr_len =", chr_len_L)
        print("L_unique =", L_unique_L)
        print("LG50 =", LG50_L)
        print("chr_len_dv =", chr_len_dv)
        print("L_unique_dv =", L_unique_dv)
        print("LG50_dv =", LG50_dv)
        #Divide for the number of simulations. You can use also the function divide_by_sim in work_parallel_functions
        chr_len_L = [x/simulations for x in chr_len_L]
        L_unique_L = [x/simulations for x in L_unique_L]
        LG50_L = [x/simulations for x in LG50_L]

        print("----------")
        print("simulations =", simulations)
        print("chr_len_L =", chr_len_L)
        print("L_unique_L =", L_unique_L)
        print("LG50_L =", LG50_L)

        chr_len = chr_len_L[:]
        L_unique = L_unique_L[:]
        LG50 = LG50_L[:]
        # plot the results
        plot_chr_len(SCRaMbLEd_events=SCRaMbLEd_events, chr_len=chr_len, L_unique=L_unique, LG50=LG50, num_essential=num_essential, simulations=simulations)
        """
        # Use this if you want histograms
        # plt.bar(SCRaMbLEd_events, chr_len, align='center')
        # plt.xticks(SCRaMbLEd_events, SCRaMbLEd_events)
        # Use this if you want lines
        plt.plot(SCRaMbLEd_events, chr_len, label="Chr len", color ="tab:blue")
        plt.plot(SCRaMbLEd_events, L_unique, label="Len unique", color ="tab:orange")
        plt.plot(SCRaMbLEd_events, LG50, label="LG50", color ="tab:green")
        #plt.errorbar(SCRaMbLEd_events, chr_len, xerr=1, yerr=chr_len_dv)
        #plt.errorbar(SCRaMbLEd_events, L_unique, xerr=1, yerr=L_unique_dv)
        #plt.errorbar(SCRaMbLEd_events, LG50, xerr=1, yerr=LG50_dv)
        plt.ylim([0, max(chr_len) + 5])
        plt.ylabel("Chr length mean over " + str(simulations) + " simulations")
        plt.xlabel("SCRaMbLEd events")
        plt.title("Chr length over SCRaMbLE events")
        plt.hlines(num_essential, 0, SCRaMbLEd_events[-1], colors="red", linestyle=":", label="Essential LU")
        plt.text(SCRaMbLEd_events[-1] * -0.16, len(syn_chr) * -0.14, "Simulations = " + str(simulations) + "\n" + "Initial chr length = " + str(len(syn_chr)) + "\n" + "Max SCRaMbLEd events = " + str(SCRaMbLEd_events[-1]) + "\n" + "Essential LU = " + str(num_essential))
        # plt.legend()
        # plt.savefig('chr_length_events.png')
        # I should add the error bar to the histograms
        plt.show()
        plt.close()
        """