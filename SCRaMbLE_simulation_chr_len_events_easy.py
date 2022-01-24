from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from comparison_sol import LoxP_unit_count_Dict_list
from comparison_sol import NG50_calculator
from SCRaMbLE_simulation_chr_len_events_store2 import plot_chr_len
from work_parallel_functions import sum_within_list
from work_parallel_functions import standard_deviation_from_points_parallel
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import statistics

# SCRaMbLEs a synthetic chromosome and outputs chr_len, L_unique, LG50
def SCRaMbLE_SIM_length(syn_chr, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1]):
    chr_len = [[] for _ in range(simulations)]
    L_unique = [[] for _ in range(simulations)]
    LG50 = [[] for _ in range(simulations)]
    chr_len_SD = [0]    # The first value is always 0
    L_unique_SD = [0]
    LG50_SD = [0]

    for s in range(simulations):
        print(s)
        SCRaMbLEd = syn_chr[:]
        for _ in range(events):
            # Perform SCRaMbLE on the synthetic chromosome
            SCRaMbLEd = force_SCRaMLE_lin_cir(SCRaMbLEd, 1, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
            # Get some statistics about the SCRaMbLEd chromosome
            NG50 = NG50_calculator(SCRaMbLEd)
            chr_len[s].append(len(SCRaMbLEd))
            L_unique[s].append(NG50[1])
            LG50[s].append(NG50[3])
    # Calculate the standard deviation
    for s in range(events):
        CHR_L = [sim_L[s] for sim_L in chr_len]
        unique_L = [sim_L[s] for sim_L in L_unique]
        LG50_L = [sim_L[s] for sim_L in LG50]
        chr_len_SD.append(statistics.stdev(CHR_L))
        L_unique_SD.append(statistics.stdev(unique_L))
        LG50_SD.append(statistics.stdev(LG50_L))
    #print("chr_len =", len(chr_len), chr_len)
    #print("L_unique =", L_unique)
    #print("LG50 =", LG50)
    #print("chr_len_SD =", len(chr_len_SD), chr_len_SD)
    #print("L_unique_SD =", L_unique_SD)
    #print("LG50_SD =", LG50_SD)
    return chr_len, L_unique, LG50, chr_len_SD, L_unique_SD, LG50_SD

def plot_SCRaMbLE_chr_len(syn_chr, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name="", SD=False):
    # SCRaMbLE the chromosome and generate chr_len, L_unique, LG50
    S = SCRaMbLE_SIM_length(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=probability)
    chr_len = S[0]
    L_unique = S[1]
    LG50 = S[2]
    if SD == True:
        chr_len_SD = S[3]
        L_unique_SD = S[4]
        LG50_SD = S[5]
    else:
        chr_len_SD = []
        L_unique_SD = []
        LG50_SD = []

    # List of time points
    SCRaMbLEd_events = list(range(0, events+1, 1))  # +1 because we record also the initial value

    # If you want to plot all the SCRaMbLE chromosomes
    #plot_chr_len(SCRaMbLEd_events, chr_len, L_unique=L_unique, LG50=LG50, num_essential=len(essential), simulations=simulations, file_name="")

    # Put together all the simulations
    chr_len_L = sum_within_list(chr_len)
    L_unique_L = sum_within_list(L_unique)
    LG50_L = sum_within_list(LG50)

    # Add the first value
    NG50 = NG50_calculator(syn_chr)
    chr_len_L.insert(0, len(syn_chr)*simulations)
    L_unique_L.insert(0, NG50[0]*simulations)
    LG50_L.insert(0, NG50[3]*simulations)

    # Divide for the number of simulations. You can use also the function divide_by_sim in work_parallel_functions
    chr_len_L = [x / simulations for x in chr_len_L]
    L_unique_L = [x / simulations for x in L_unique_L]
    LG50_L = [x / simulations for x in LG50_L]
    #print("chr_len_L =", chr_len_L)
    #print("L_unique_L =", L_unique_L)
    #print("LG50_L =", LG50_L)

    # Plot the results
    plot_chr_len(SCRaMbLEd_events, chr_len_L, L_unique=L_unique_L, LG50=LG50_L, num_essential=len(essential), simulations=simulations, circular=circular, probability=probability, file_name=file_name, SD_chr_len=chr_len_SD, SD_L_unique=L_unique_SD, SD_LG50=LG50_SD)
    return None

# I split the function plot_SCRaMbLE_chr_len in two pieces
def mean_SCRaMbLE_chr_len(syn_chr, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name=""):
    # Generate the chr len
    S = SCRaMbLE_SIM_length(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=probability)
    chr_len = S[0]
    L_unique = S[1]
    LG50 = S[2]

    # If you want to plot all the SCRaMbLE chromosomes
    #plot_chr_len(SCRaMbLEd_events, chr_len, L_unique=L_unique, LG50=LG50, num_essential=len(essential), simulations=simulations, file_name="")

    # Put together all the simulations
    chr_len_L = sum_within_list(chr_len)
    L_unique_L = sum_within_list(L_unique)
    LG50_L = sum_within_list(LG50)

    # Add the first value
    NG50 = NG50_calculator(syn_chr)
    chr_len_L.insert(0, len(syn_chr)*simulations)
    L_unique_L.insert(0, NG50[0]*simulations)
    LG50_L.insert(0, NG50[3]*simulations)

    # Divide for the number of simulations. You can use also the function divide_by_sim in work_parallel_functions
    chr_len_L = [x / simulations for x in chr_len_L]
    L_unique_L = [x / simulations for x in L_unique_L]
    LG50_L = [x / simulations for x in LG50_L]

    return chr_len_L, L_unique_L, LG50_L


def plot_mean_SCRaMbLE_chr_len(syn_chr, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name=""):
    A_result = mean_SCRaMbLE_chr_len(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=probability)

    # List of time points
    SCRaMbLEd_events = list(range(0, events + 1, 1))  # +1 because we record also the initial value

    # Plot the results
    plot_chr_len(SCRaMbLEd_events, chr_len=A_result[0], L_unique=A_result[1], LG50=A_result[2], num_essential=len(essential), simulations=simulations, circular=circular, probability=probability, file_name=file_name)
    return None

def chr_len_range_SCRaMbLE(events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name=""):

    #chr_len_range = list(range(25, 201, 25))
    chr_len_range = list(range(25, 226, 50))
    # chr_len_range = [25, 50, 75, 100, 125, 150, 175, 200]
    # chr_len_range = [25, 75, 125, 175, 225]
    print("chr_len_range =", chr_len_range)

    chr_len_multi = []
    L_unique_multi = []
    LG50_multi = []

    for i in chr_len_range:
        print("Chr len =", i)
        # Create the chromosomes
        SYN_CHR = list(range(1, i, 1))
        A_result = mean_SCRaMbLE_chr_len(SYN_CHR, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=probability, file_name=file_name)
        chr_len_multi.append(A_result[0])
        L_unique_multi.append(A_result[1])
        LG50_multi.append(A_result[2])

    # List of time points
    SCRaMbLEd_events = list(range(0, events + 1, 1))  # +1 because we record also the initial value

    # plot the results
    plt.figure(figsize=(18, 9))

    # Font size
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 24
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.ylabel("Chr length mean over " + str(simulations) + " simulations")
    plt.xlabel("SCRaMbLEd events")
    plt.title("Chr length over SCRaMbLE events")

    #plt.plot(0, 0, label="Chr len", color="tab:blue")
    plt.plot(0, 0, label="Len unique", color="cornflowerblue")
    # plt.plot(0,0, label="LG50", color="tab:green")
    #Labels = ["Chr = 25", "Chr = 50", "Chr = 75", "Chr = 100", "Chr = 125", "Chr = 150", "Chr = 175", "Chr = 200"]
    Labels = ["Chr = 25", "Chr = 75", "Chr = 125", "Chr = 175", "Chr = 225"]

    for i in range(len(chr_len_multi)):
        plt.plot(SCRaMbLEd_events, chr_len_multi[i], label=Labels[i])
        # plt.errorbar(SCRaMbLEd_events, chr_len[i], yerr=SD_chr_len, elinewidth=0.1)
        if L_unique_multi != [[]]:
            plt.plot(SCRaMbLEd_events, L_unique_multi[i], color="cornflowerblue")
        #if LG50 != [[]]:
        #    plt.plot(SCRaMbLEd_events, LG50[i], label="LG50", color="tab:green")

    plt.hlines(num_essential, 0, SCRaMbLEd_events[-1], colors="red", linestyle=":", label="Essential LU")
    plt.legend()

    # These are some information to add to the saved files.
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    if circular:  # record if the chromosome is linear or circular
        lin_cir = "c"
    else:
        lin_cir = "l"
    probability_str = ""  # record the probabilities of each event
    for i in probability:
        probability_str = probability_str + str(i)

    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".png", dpi=200)
    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".svg", format="svg", dpi=200)
    plt.show()
    plt.close()
    return None

def chr_len_essential_range_SCRaMbLE(syn_chr=50, events=15, simulations=100, CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name=""):

    essential_range = list(range(0, 51, 10))
    # essential_range = [0, 10, 20, 30, 40, 50]
    print("essential_range =", essential_range)

    chr_len_multi = []
    L_unique_multi = []
    LG50_multi = []

    for ESSE in essential_range:
        print("num essential =", ESSE)
        if ESSE == 0:
            essential = []
        else:
            essential = sorted(random.sample(syn_chr, k=ESSE))
        A_result = mean_SCRaMbLE_chr_len(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=probability, file_name=file_name)
        chr_len_multi.append(A_result[0])
        L_unique_multi.append(A_result[1])
        LG50_multi.append(A_result[2])

    # List of time points
    SCRaMbLEd_events = list(range(0, events + 1, 1))  # +1 because we record also the initial value

    # plot the results
    plt.figure(figsize=(18, 9))

    # Font size
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 24
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.ylabel("Chr length mean over " + str(simulations) + " simulations")
    plt.xlabel("SCRaMbLEd events")
    plt.title("Chr length over SCRaMbLE events")

    Labels = ["# essential = 0", "# essential = 10", "# essential = 20", "# essential = 30", "# essential = 40", "# essential = 50"]

    for i in range(len(chr_len_multi)):
        plt.plot(SCRaMbLEd_events, chr_len_multi[i], label=Labels[i])
        # if L_unique != [[]]:
        #    plt.plot(SCRaMbLEd_events, L_unique[i], color="cornflowerblue")
    plt.legend()
    # These are some information to add to the saved files.
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    if circular:  # record if the chromosome is linear or circular
        lin_cir = "c"
    else:
        lin_cir = "l"
    probability_str = ""  # record the probabilities of each event
    for i in probability:
        probability_str = probability_str + str(i)

    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".png", dpi=200)
    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".svg", format="svg", dpi=200)
    plt.show()
    plt.close()
    return None

def chr_len_probabilities_range_SCRaMbLE(syn_chr=50, events=15, simulations=100, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], file_name=""):

    probabilities_range = [[0,4,4,1], [0,4,4,2], [0,4,4,4], [0,4,4,8], [0,4,4,16]]

    chr_len_multi = []
    L_unique_multi = []
    LG50_multi = []

    for PROB in probabilities_range:
        print("Probability =", PROB)
        A_result = mean_SCRaMbLE_chr_len(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=circular, mu=mu, sigma=sigma, force=force, probability=PROB, file_name=file_name)
        chr_len_multi.append(A_result[0])
        L_unique_multi.append(A_result[1])
        LG50_multi.append(A_result[2])

    # List of time points
    SCRaMbLEd_events = list(range(0, events + 1, 1))  # +1 because we record also the initial value

    # plot the results
    plt.figure(figsize=(18, 9))

    # Font size
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 24
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.ylabel("Chr length mean over " + str(simulations) + " simulations")
    plt.xlabel("SCRaMbLEd events")
    plt.title("Chr length over SCRaMbLE events")

    Labels = ["P(Del)/P(Dup) = 4", "P(Del)/P(Dup) = 2", "P(Del)/P(Dup) = 1", "P(Del)/P(Dup) = 1/2", "P(Del)/P(Dup) = 1/4"]

    for i in range(len(chr_len_multi)):
        plt.plot(SCRaMbLEd_events, chr_len_multi[i], label=Labels[i])
        # if L_unique != [[]]:
        #    plt.plot(SCRaMbLEd_events, L_unique[i], color="cornflowerblue")
    plt.hlines(num_essential, 0, SCRaMbLEd_events[-1], colors="red", linestyle=":", label="Essential LU")
    plt.legend()
    # These are some information to add to the saved files.
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    if circular:  # record if the chromosome is linear or circular
        lin_cir = "c"
    else:
        lin_cir = "l"
    probability_str = ""  # record the probabilities of each event
    for i in probability:
        probability_str = probability_str + str(i)

    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".png", dpi=200)
    plt.savefig("SCRaMbLE_chr_len/chr_length_" + lin_cir + "_events_" + str(SCRaMbLEd_events[-1]) + "_sim_" + str(simulations) + "_P" + probability_str + "_" + file_name + random_seed + ".svg", format="svg", dpi=200)
    plt.show()
    plt.close()
    return None

# test the code
if __name__ == "__main__":

    syn_chr = list(range(1, 45, 1))
    essential = [2, 7, 9, 10, 12, 20]
    #num_essential = 10
    #essential = sorted(random.sample(syn_chr, k=num_essential))
    CEN = [2]
    events = 1000
    simulations = 200

    plot_SCRaMbLE_chr_len(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=True, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], SD=True)

    #plot_mean_SCRaMbLE_chr_len(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])

    # chr_len_range_SCRaMbLE
    syn_chr = list(range(1, 25, 1))
    num_essential = 10
    essential = sorted(random.sample(syn_chr, k=num_essential))
    #chr_len_range_SCRaMbLE(events=events, simulations=simulations, essential=essential, CEN=CEN, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])

    # chr_len_essential_range_SCRaMbLE
    syn_chr = list(range(1, 51, 1))
    #chr_len_essential_range_SCRaMbLE(syn_chr, events=events, simulations=simulations, CEN=CEN, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])

    # chr_len_probabilities_range_SCRaMbLE
    syn_chr = list(range(1, 51, 1))
    num_essential = 10
    events = 250
    simulations = 50
    essential = sorted(random.sample(syn_chr, k=num_essential))
    #chr_len_probabilities_range_SCRaMbLE(syn_chr, events=events, simulations=simulations, essential=essential, CEN=CEN, circular=False, mu=0, sigma=10, force=True)
