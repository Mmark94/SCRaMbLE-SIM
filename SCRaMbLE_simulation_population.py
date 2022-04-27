# SCRaMbLE simulation population
import random
import numpy as np
import matplotlib.pyplot as plt
from SCRaMbLE_simulation_3 import SCRaMbLE4
from SCRaMbLE_simulation_3_circular import SCRaMbLE4_circular


# This add one to all the elements of the list (cells)
def list_plus_one(cells: list):
    new_cells = [x+1 for x in cells]
    return new_cells

# This calculate the frequency of each cell in the population
def dict_percentage(cells: list):
    D_count = {}
    D_percentage = {}
    unique = list(dict.fromkeys(cells))
    L_cells = len(cells)
    for cell_u in unique:
        D_count[cell_u] = cells.count(cell_u)
        D_percentage[cell_u] = round(100*(cells.count(cell_u) / L_cells), 2)
    return D_count, D_percentage

# This calculate the amount of SCRaMbLEd cells and how many SCRaMbLE events there are in a population. The Cre recombinase is under a daughter-specific promoter.
# One replication in yeast takes around 90 min. Therefore, SCRaMbLE for 12 hours means that the cells undergo to replication 8 times.
def daughter_promoter(initial_cells=1, number_replication=8, senescence=20):
    cells = [0 for _ in range(initial_cells)]
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        new_cells = list_plus_one(cells)
        # Remove old cells that entered in senescence
        #new_cells = [x for x in new_cells if x != senescence]

        # The cells double at each replication
        cells = cells + new_cells
    cells.sort()
    Dict_count = dict_percentage(cells)[0]
    Dict_percentage = dict_percentage(cells)[1]
    return Dict_count, Dict_percentage


# This function simulates a SCRaMbLE population, with many different genotypes, some of them SCRaMbLEd others unSCRaMbLEd
def SCRaMbLE_population(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=7, sigma=7, CEN=[]):
    cells = [syn_chr for _ in range(initial_cells)]
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        new_cells = [SCRaMbLE4(Chr, events_for_replication, essential=essential, mu=mu, sigma=sigma, CEN=CEN) for Chr in cells]
        # The cells double at each replication
        cells = cells + new_cells
    return cells


# This function simulates a SCRaMbLE population, with many different genotypes, some of them SCRaMbLEd others unSCRaMbLEd
# In this function I want to exclude the cells that deleted an essential gene.
def SCRaMbLE_population2(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=7, sigma=7, CEN=[]):
    cells = [syn_chr for _ in range(initial_cells)]
    new_cells = []
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        for Chr in cells:
            SCRaMbLEd_cell = SCRaMbLE4(Chr, events_for_replication, essential=essential, mu=mu, sigma=sigma, CEN=CEN)
            # If an essential LU is deleted, the function SCRaMbLE4 will out the same chromosome. Therefore, if the chromosome is the same, I want to remove it.
            if SCRaMbLEd_cell != Chr:
                new_cells.append(SCRaMbLEd_cell)
        # The cells double at each replication
        cells = cells + new_cells
    return cells

def death_rate_by_SCRaMbLE(population, survivor_rate):
    # death_rate = 1 - survivor_rate
    survivors = {}
    for i in population:
        if i == 0:  # this excludes the error: division by zero
            survivors[i] = population[i]
        else:
            survivors[i] = round(population[i] * (survivor_rate ** i), 4)
    # count how many cells survived
    survivors_percentage = 0
    for i in survivors:
        survivors_percentage = survivors_percentage + survivors[i]
    print("survivors_percentage with survivor_rate", survivor_rate, " =", survivors_percentage)
    return survivors

def plot_Dic2(myDictionary, file_name="SCRaMbLEd_cells"):
    # Plot
    plt.figure(figsize=(7, 3.5), dpi=300)

    # Font size
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Plot
    plt.bar(myDictionary.keys(), myDictionary.values(), align='center')
    plt.xticks(range(max(myDictionary.keys()) + 1), range(max(myDictionary.keys()) + 1))
    plt.ylabel("Number cells")
    plt.xlabel("SCRaMbLE events")
    plt.title(file_name)
    plt.savefig("SCRaMbLE_population/" + file_name + ".png", dpi=300, bbox_inches='tight')
    plt.savefig("SCRaMbLE_population/" + file_name + ".svg", format='svg', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    return None

def plot_Dic_survivor(myDictionary, myDictionary2, file_name="survivor_SCRaMbLEd_cells", survivor_rate = 0.5):
    # myDictionary = All cells, myDictionary2 = survivor cells
    # Plot
    plt.figure(figsize=(7, 3.5), dpi=300)

    # Font size
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Plot
    plt.bar(myDictionary.keys(), myDictionary.values(), align='center', color='indianred', label="Dead cells")
    plt.bar(myDictionary2.keys(), myDictionary2.values(), align='center', color='limegreen', label="Alive cells")
    plt.xticks(range(max(myDictionary.keys()) + 1), range(max(myDictionary.keys()) + 1))
    plt.ylabel("Number of cells")
    plt.xlabel("Number of SCRaMbLE events")
    plt.title(file_name)
    plt.text(-0.5, max(myDictionary.values()) * 0.95, "SCRaMbLE event survivor rate = " + str(survivor_rate))
    plt.savefig("SCRaMbLE_population/" + file_name + "_SR_" + str(survivor_rate) + ".png", dpi=300, bbox_inches='tight')
    plt.savefig("SCRaMbLE_population/" + file_name + "_SR_" + str(survivor_rate) + ".svg", format='svg', dpi=300, bbox_inches='tight')
    plt.legend()
    #plt.show()
    plt.close()
    return None
# I use the following website to create giff: https://ezgif.com/maker

def death_rate_by_SCRaMbLE_plot(initial_cells=1, number_replication=8, survivor_rate=0.5):
    cells = daughter_promoter(initial_cells, number_replication)[1]
    survivor = death_rate_by_SCRaMbLE(cells, survivor_rate)
    plot_Dic_survivor(cells, survivor, "survivor_SCRaMbLEd_cells", survivor_rate)
    return None

# test the code
if __name__ == "__main__":

    segments = 44  #number of loxP segments
    syn_chr = list(range(1, segments+1, 1))
    essential = [2, 7, 9, 10, 12, 19, 20, 24]  # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.
    print("syn_chr =", syn_chr)

    Cell = daughter_promoter(3, 8)
    #print(Cell[0])
    print(Cell[1])
    for i in np.arange(0.0, 1.1, 0.1):
        i = round(i, 1)
        print(i)
        survivor = death_rate_by_SCRaMbLE(Cell[1], i)
        print(survivor)

    #plot_Dic_survivor(Cell[1], survivor, "survivor")
    #death_rate_by_SCRaMbLE_plot(1, 8, 0.7)
    range_list = [round(x * 0.1, 1) for x in range(0, 11)]
    for i in range_list:
        death_rate_by_SCRaMbLE_plot(1, 8, i)
    plot_Dic2(Cell[1], "SCRaMbLEd_cells")
    plot_Dic2(survivor, "survivor")

    pop = SCRaMbLE_population(syn_chr, 1, 4, 5)
    for i in pop:
        print(i)
