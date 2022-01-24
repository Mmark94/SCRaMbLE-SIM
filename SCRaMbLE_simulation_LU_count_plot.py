from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
from comparison_sol import LoxP_unit_count_Dict_list
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random


def SCRaMbLE_evolution_LU_count(syn_chr, SCRaMbLE_events, essential, CEN, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1]):
    # LU_count_time stores the number of each LU through time
    LU_count_time = {}
    # At the beginning, all the LU are in one copy
    for LU in syn_chr:
        LU_count_time[LU] = [1]

    SCRaMbLEd = syn_chr[:]
    for i in range(SCRaMbLE_events):
        # Perform SCRaMbLE on the synthetic chromosome
        SCRaMbLEd = force_SCRaMLE_lin_cir(SCRaMbLEd, 1, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
        #print(SCRaMbLEd)
        # Count the LU in the chromosome after SCRaMbLE
        LU_count = LoxP_unit_count_Dict_list(SCRaMbLEd, syn_chr)
        # print(LU_count)

        # Add the LU count to the LU_count_time dictionary. k is the LU, v is the copy number of the LU k.
        for k, v in LU_count.items():
            LU_count_time[k].append(v)
    print("LU_count_time =", LU_count_time)

    # Plot the data

    # get discrete colormap
    # You can choose different gradients here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    colors_non_esse = cm.Blues(np.linspace(0, 1, len(syn_chr) - len(essential)))
    colors_esse = cm.Reds(np.linspace(0, 1, len(essential)))
    ALL_colour = np.zeros(shape=(len(syn_chr), 4))
    counter_esse = 0
    counter_non_esse = 0
    for LU in syn_chr:
        if abs(LU) in essential:
            ALL_colour[LU - 1] = colors_esse[counter_esse]
            counter_esse += 1
        else:
            ALL_colour[LU - 1] = colors_non_esse[counter_non_esse]
            counter_non_esse += 1
    # print(ALL_colour)

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

    # Start creating the plot
    fig, ax = plt.subplots(figsize=(20, 10))

    SCRaMbLEd_events_list = list(range(0, SCRaMbLE_events + 1))
    ax.stackplot(SCRaMbLEd_events_list, LU_count_time.values(), labels=LU_count_time.keys(), colors=ALL_colour)
    ax.legend(loc='upper left', prop={'size': 8})
    ax.set_title('LoxP Units Evolution')
    ax.set_xlabel('SCRaMbLEd_events')
    ax.set_ylabel('Number of LoxP Units')
    # ax.text(45, 60, 'Essential ' + str(essential), fontsize=10)

    # These are some information to add to the saved files.
    # Create a random seed to save the image
    random_seed = str(random.random())[2:6]
    if circular:    # record if the chromosome is linear or circular
        lin_cir = "c"
    else:
        lin_cir = "l"
    probability_str = ""    # record the probabilities of each event
    for i in probability:
        probability_str = probability_str + str(i)
    plt.savefig("SCRaMbLE_evolution_LU/SCRaMbLE_evolution_LU_count_chr_" + lin_cir + "_L" + str(len(syn_chr)) + "_Ev" + str(SCRaMbLE_events) + "_P" + probability_str + "_s" + random_seed + ".png", dpi=200)
    plt.savefig("SCRaMbLE_evolution_LU/SCRaMbLE_evolution_LU_count_chr_" + lin_cir + "_L" + str(len(syn_chr)) + "_Ev" + str(SCRaMbLE_events) + "_P" + probability_str + "_s" + random_seed + ".svg", format="svg", dpi=200)
    #plt.show()
    return None

# test the code
if __name__ == "__main__":

    # Syn9R SCRaMbLE
    syn_chr = list(range(1, 45, 1))
    essential = [2, 7, 9, 10, 12, 20]
    CEN = [2]

    #SCRaMbLE_evolution_LU_count(syn_chr, 100, essential=essential, CEN=CEN, circular=False)

    for i in range(5):
        print(i)
        SCRaMbLE_evolution_LU_count(syn_chr, 500, essential=essential, CEN=CEN, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])

    """
    # Syn3 SCRaMbLE
    syn3 = list(range(1, 100, 1))
    essential_3 = [7, 9, 13, 20, 28, 36, 38, 48, 70, 73, 74, 84, 97]
    CEN3 = [38]

    for i in range(10):
        print(i)
        SCRaMbLE_evolution_LU_count(syn3, 1000, essential=essential_3, CEN=CEN3, circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1])
    """
