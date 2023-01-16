# SCRaMbLE simulation population
import random
import statistics
import numpy as np
import matplotlib.pyplot as plt
from SCRaMbLE_simulation_3 import SCRaMbLE4
from SCRaMbLE_simulation_3_circular import SCRaMbLE4_circular

# This script can be used to simulate SCRaMbLE in a yeast population.
# So far, we have modelled and simulated evolution through a single genome/chromosome that keeps accumulating SCRaMbLE events.
# However, a different approach consists of modelling and simulating the dynamics of an entire SCRaMbLE population where each cell has a different number of SEs and cells share common ancestors and, therefore, some SEs.
# As a proof of concept, we modelled and simulated one of the most common SCRaMbLE protocols where the Cre recombinase gene is placed under a daughter-specific promoter, and it is expressed only in the daughter cell but not in the mother cell (Lindstrom & Gottschling, 2009).
# Therefore, at each generation, the cells replicate, and only the daughter cells undergo an x number of SCRaMbLE events (for our simulations, we assumed M = 1 SE for replication).

# This add one to all the elements of the list (cells)
def list_plus_one(cells: list):
    new_cells = [x+1 for x in cells]
    return new_cells

# This calculate the frequency of each cell in the population.
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
def SCRaMbLE_population(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    cells = [syn_chr for _ in range(initial_cells)]
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        new_cells = [SCRaMbLE4(Chr, events_for_replication, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability) for Chr in cells]
        # The cells double at each replication
        cells = cells + new_cells
    return cells

# This function simulates a SCRaMbLE population, with many different genotypes, some of them SCRaMbLEd others unSCRaMbLEd
def SCRaMbLE_population_dict(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    cells = [syn_chr for _ in range(initial_cells)]
    # Record the number of SCRaMbLE events in each cell
    num_SE = [0 for _ in range(initial_cells)]
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        new_cells = [SCRaMbLE4(Chr, events_for_replication, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability) for Chr in cells]
        # Record the number of SCRaMbLE events in each cell
        new_cells_num_SE = [SE + events_for_replication for SE in num_SE]
        num_SE = num_SE + new_cells_num_SE
        # The cells double at each replication
        cells = cells + new_cells
    #print("num_SE =", num_SE)

    pop_dict = {}
    for i in range(0, number_replication+1):
        pop_dict[i * events_for_replication] = []
    #print("pop_dict =", pop_dict)
    for i in range(len(cells)):
        pop_dict[num_SE[i]].append(cells[i])
    #print("pop_dict =", pop_dict)
    return pop_dict

# This function simulates a SCRaMbLE population, with many different genotypes, some of them SCRaMbLEd others unSCRaMbLEd
# In this function I want to exclude the cells that deleted an essential gene.
def SCRaMbLE_population2(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    cells = [syn_chr for _ in range(initial_cells)]
    new_cells = []
    for i in range(number_replication):
        # The daughter cells (half of them) will have plus one SCRaMbLE event
        for Chr in cells:
            SCRaMbLEd_cell = SCRaMbLE4(Chr, events_for_replication, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            # If an essential LU is deleted, the function SCRaMbLE4 will out the same chromosome. Therefore, if the chromosome is the same, I want to remove it.
            if SCRaMbLEd_cell != Chr:
                new_cells.append(SCRaMbLEd_cell)
        # The cells double at each replication
        cells = cells + new_cells
    return cells

# Population structure is a dictionary with the keys = number of SEs and the values the percentage of the cells in the population with that many SEs.
def death_rate_by_SCRaMbLE(population, survivor_rate):
    # death_rate = 1 - survivor_rate
    survivors = {}
    for i in population:
        if i == 0:  # this excludes the error: division by zero
            survivors[i] = population[i]
        else:
            # The fraction of cells population[i] is multiply by the survivor_rate elevated to the power of how many SEs it received.
            survivors[i] = round(population[i] * (survivor_rate ** i), 4)
    # Count how many cells survived
    survivors_percentage = 0
    for i in survivors:
        survivors_percentage = survivors_percentage + survivors[i]
    print("Survival rate of a single SCRaMbLE event =", survivor_rate)
    print("population =", population)
    print("survivors  =", survivors)
    print("Fraction of the cells surviving =", survivors_percentage)
    print()
    return survivors

# Correlation between SCRaMbLE events mortality rate and the amount of cell surviving wihin a SCRaMbLE population.
def plot_SE_mortality_vs_cell_survival(population: dict):
    mortality_range = np.arange(0.0, 1.05, 0.05)
    cell_survival_list = []
    all_cells = sum(population.values())
    for mortality in mortality_range:
        survivor_rate = 1 - mortality
        cell_survival_dict = {}
        for SE, cells in population.items():
            cell_survival_dict[SE] = cells * (survivor_rate ** SE)
        cell_survival_list.append(sum(cell_survival_dict.values()) / all_cells)
    print("mortality_range =", mortality_range)
    print("cell_survival_list =", cell_survival_list)
    # Plot SE_mortality vs cell_survival
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
    plt.scatter(mortality_range, cell_survival_list, s=12)
    plt.ylabel("Percentage of cells (c) surviving (A)")
    plt.xlabel("SCRaMbLE events (E) mortality rate (m)")
    plt.title("SCRaMbLE events mortality vs percentage of cells surviving", fontsize=MEDIUM_SIZE)
    #plt.xlim(0, 1)
    #plt.ylim(0, 1)
    plt.savefig("SCRaMbLE_population/SE_mortality_vs_cell_survival2.png", dpi=300, bbox_inches='tight')
    plt.savefig("SCRaMbLE_population/SE_mortality_vs_cell_survival2.svg", format='svg', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    return None

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

def plot_Dic_survivor(myDictionary, myDictionary2, file_name="survivor_SCRaMbLEd_cells", mortality_rate=0.5, Show=False):
    mortality_rate = round(mortality_rate, 2)
    # myDictionary = All cells, myDictionary2 = survivor cells
    # Count how many cells survived
    survivors_percentage = sum(myDictionary2.values())
    all_cells = sum(myDictionary.values())
    survivors_percentage = survivors_percentage / all_cells
    print("Cells          =", myDictionary)
    print("survivor cells =", myDictionary2)
    print("survivor cells =", survivors_percentage)
    print("All cells      =", all_cells)

    #for i in myDictionary2.keys():
    #    myDictionary[i] = myDictionary[i] / all_cells
    #    myDictionary2[i] = myDictionary2[i] / all_cells

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
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

    # Plot
    plt.bar(myDictionary.keys(), myDictionary.values(), align='center', color='indianred', label="Dead cells")
    plt.bar(myDictionary2.keys(), myDictionary2.values(), align='center', color='limegreen', label="Alive cells")
    plt.xticks(range(max(myDictionary.keys()) + 1), range(max(myDictionary.keys()) + 1))
    plt.ylabel("Percentage of cells")
    #plt.ylabel("Number of cells")
    plt.xlabel("Number of SCRaMbLE events")
    #plt.title("Survival of a SCRaMbLEd population")    # file_name
    #plt.text(-0.6, max(myDictionary.values()) * 0.92, "SCRaMbLE event mortality rate = " + str(mortality_rate) + "\n" + "Percentage of the cells surviving = " + str(round(survivors_percentage, 2)))
    plt.text(-0.6, max(myDictionary.values()) * 0.92, "SE mortality = " + str(round(100*mortality_rate)) + " %" + "\n" + "Cells surviving = " + str(round(100*survivors_percentage, 1)) + " %")
    #plt.text(-0.6, max(myDictionary.values()) * 0.92, "Percentage of the cells surviving = " + str(100*round(survivors_percentage, 3)) + " %")
    plt.savefig("SCRaMbLE_population/" + file_name + "_MR_" + str(mortality_rate) + ".png", dpi=300, bbox_inches="tight")
    plt.savefig("SCRaMbLE_population/" + file_name + "_MR_" + str(mortality_rate) + ".svg", format='svg', dpi=300, bbox_inches="tight")
    plt.legend(fontsize=MEDIUM_SIZE)
    if Show:
        plt.show()
    plt.close()
    return None
# I use the following website to create giff: https://ezgif.com/maker

def death_rate_by_SCRaMbLE_plot(initial_cells=1, number_replication=8, survivor_rate=0.5):
    cells = daughter_promoter(initial_cells, number_replication)[1]
    survivor = death_rate_by_SCRaMbLE(cells, survivor_rate)
    mortality_rate = 1 - survivor_rate
    plot_Dic_survivor(cells, survivor, file_name="Survival rate of SCRaMbLEd population", mortality_rate=mortality_rate)
    return None

# This output True / False
def check_essential_LUs_are_in_chr(Chr, essential=[]):
    new_chr_abs = [abs(ele) for ele in Chr]
    # It checks if all the essential segments are in the chromosome
    essential_LUs_present = all(elem in new_chr_abs for elem in essential)
    # print(new_chr, essential_LUs_present)
    return essential_LUs_present

# Calculate the probability that SCRaMbLE deletes an essential LU.
def percentage_of_chrs_with_essential_LUs_deleted(Chr, simulations=1000, Number_events=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    if Chr == []:
        return 0
    if isinstance(Chr[0], int): # There is only one chromosome
        essential_deleted = 0
        for _ in range(simulations):
            # SCRaMbLE the chromosome
            # You need to set essential=[] in SCRaMbLE4 so essential LUs can be deleted
            new_chr = SCRaMbLE4(Chr, Number_events=Number_events, essential=[], mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            # It checks if all the essential segments are in the synthetic chromosome after the SCRaMbLE events happened.
            essential_LUs_present = check_essential_LUs_are_in_chr(Chr=new_chr, essential=essential)
            #print(new_chr, essential_LUs_present)
            if essential_LUs_present is False:
                # An essential LU has been deleted
                essential_deleted = essential_deleted + 1
        percentage_esse_deleted = essential_deleted / simulations
        #print("Percentage of chromosomes with essential LUs deleted =", percentage_esse_deleted, ". Simulations =", simulations)
        return percentage_esse_deleted
    else:   # There multiple chromosomes, take the average
        percentage_esse_deleted_list = []
        for single_chr in Chr:
            percentage_esse_deleted_list.append(percentage_of_chrs_with_essential_LUs_deleted(single_chr, simulations=simulations, Number_events=Number_events, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability))
        percentage_esse_deleted_average = sum(percentage_esse_deleted_list) / len(percentage_esse_deleted_list)
        print("percentage_esse_deleted_list =", percentage_esse_deleted_list)
        print("percentage_esse_deleted_average =", percentage_esse_deleted_average, ". n =", len(percentage_esse_deleted_list))
        return percentage_esse_deleted_average

# Simulate the SCRaMbLE population. Cre under the daughter-specific promoter.
def simulate_SCRaMbLE_pop_check_survival_rate(syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    # You need to set essential=[] in SCRaMbLE4 so essential LUs can be deleted
    pop_dict = SCRaMbLE_population_dict(syn_chr=syn_chr, initial_cells=initial_cells, number_replication=number_replication, events_for_replication=events_for_replication, essential=[], mu=mu, sigma=sigma, CEN=CEN, probability=probability)
    num_cells_dict = {}
    dead_cells_dict = {}
    survivor_cells_dict = {}
    mortality_list = []
    # Count how many cells survived
    for SE in pop_dict.keys():
        print("SE =", SE)
        num_cells = len(pop_dict[SE])
        num_cells_dict[SE] = num_cells
        dead_cells = 0
        alive_cells_list = []
        for cell in pop_dict[SE]:
            essential_LUs_present = check_essential_LUs_are_in_chr(Chr=cell, essential=essential)
            if essential_LUs_present is False:
                dead_cells = dead_cells + 1
            else:
                alive_cells_list.append(cell)
        dead_cells_dict[SE] = dead_cells
        survivor_cells_dict[SE] = num_cells - dead_cells
        # Find out mortality rate. You can use pop_dict[SE] but there are cells already dead. Better to use alive_cells_list
        mortality = percentage_of_chrs_with_essential_LUs_deleted(alive_cells_list, simulations=100, Number_events=1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
        mortality_list.append(mortality)
        print("mortality =", mortality)
    print("mortality_list =", mortality_list)
    print("mortality_Q1_Q2_Q3 =", len(mortality_list), "=", np.percentile(mortality_list, 25), np.percentile(mortality_list, 50), np.percentile(mortality_list, 75))
    print("mortality_mean =", np.average(mortality_list), "±", np.std(mortality_list))
    return num_cells_dict, survivor_cells_dict

# test the code
if __name__ == "__main__":

    segments = 44  #number of loxP segments
    syn_chr = list(range(1, segments+1, 1))
    essential = [2, 7, 9, 10, 12, 19, 20, 24]  # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.
    print("syn_chr =", syn_chr)
    #"""
    print(percentage_of_chrs_with_essential_LUs_deleted(syn_chr, Number_events=1, simulations=1000, essential=essential, mu=0, sigma=10, CEN=[2], probability=[3, 2, 2, 1]))
    
    Cell = daughter_promoter(3, 8)
    print(Cell[0])
    print(Cell[1])
    for i in np.arange(0.0, 1.1, 0.1):
        i = round(i, 1)
        print(i)
        survivor = death_rate_by_SCRaMbLE(Cell[1], i)
        print(survivor)
    
    #plot_Dic_survivor(Cell[1], survivor, "survivor")
    #death_rate_by_SCRaMbLE_plot(1, 8, 0.7)
    for i in np.arange(0.0, 1.05, 0.1):
        death_rate_by_SCRaMbLE_plot(initial_cells=1, number_replication=8, survivor_rate=i)
    #plot_Dic2(Cell[1], "SCRaMbLEd_cells")
    #plot_Dic2(survivor, "survivor")
    #"""
    #pop = SCRaMbLE_population(syn_chr=syn_chr, initial_cells=1, number_replication=4, events_for_replication=5)
    #for i in pop:
    #    print(i)
    #syn_chr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    #essential = [2, 7, 9, 10, 12, 19, 20, 24, 28, 30, 35, 40]
    #print("survivor_rate =", 1 - percentage_of_chrs_with_essential_LUs_deleted(syn_chr, Number_events=1, simulations=1000, essential=essential, mu=0, sigma=10, CEN=[2], probability=[0, 2, 2, 1]))
    #pop_dict = SCRaMbLE_population_dict(syn_chr=syn_chr, initial_cells=1, number_replication=8, events_for_replication=1, essential=essential, mu=0, sigma=10, CEN=[], probability=[0, 2, 2, 1])
    #print("pop_dict =", pop_dict)

    pop_survival = simulate_SCRaMbLE_pop_check_survival_rate(syn_chr=syn_chr, initial_cells=10, number_replication=8, events_for_replication=1, essential=essential, mu=0, sigma=10, CEN=[], probability=[0, 2, 2, 1])
    print("pop_survival =", pop_survival)
    plot_Dic_survivor(pop_survival[0], pop_survival[1], "population_SIM", mortality_rate=0.206, Show=True)

    Cell = daughter_promoter(3, 8)
    plot_SE_mortality_vs_cell_survival(Cell[1])

    """
    # 2022-05-07
    chr len = 44 LUs
    No. essential LUs = 6
    PDEL/PDUP = 2
    No. cells = 100*2^8 = 25,600
    cells surviving = 38.9 %
    mortality_list = [0.21944, 0.21455, 0.21054, 0.20676, 0.20019, 0.19477, 0.19298, 0.19006, 0.22433]
    mortality_Q1_Q2_Q3 = 0.19477 0.20676 0.21455
    mortality_mean = 0.20596 ± 0.01152
    """
