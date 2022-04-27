# Mapping
from comparison_sol import LoxP_unit_count_Dict_list
from comparison_sol import invert
from SCRaMbLE_simulation_3 import force_SCRaMLE_lin_cir
import matplotlib.pyplot as plt
import statistics
import numpy as np
import matplotlib.cm as cm
import random

def mapping(solution, path):
    if solution == [] or path == []:
        return False
    L_path = len(path)
    c=0
    while c + L_path< len(solution) + 1:
        k_mer = solution[c:L_path+c]
        if path == k_mer or invert(path) == k_mer:
            return [c, L_path+c-1, True]
        c = c + 1
    return False

def paths_in_sol(solution, paths):
    if solution == [] or paths == [] or paths == [[]]:
        return False
    for path in paths:
        if path == []:
            continue
        if not(mapping(solution, path)):
            return False
    return True

def extract_unmapped_reads(solution, paths):
    if solution == [] or paths == [] or paths == [[]]:
        return False
    unmapped_reads = []
    for path in paths:
        if path == []:
            continue
        if not(mapping(solution, path)):
            unmapped_reads.append(path)
    return unmapped_reads

def remove_duplicate(paths:list):
    if isinstance(paths[0], list):
        new_paths = []
        for path in paths:
            if path not in new_paths:
                new_paths.append(path)
        return new_paths
    else:
        return paths

def check_sol(solutions, paths):
    if isinstance(solutions[0], list):
        solutions = remove_duplicate(solutions)
        new_solutions = []
        for sol in solutions:
            if paths_in_sol(sol, paths):
                new_solutions.append(sol)
        return new_solutions
        #return new_solutions[0] if len(new_solutions)==1 else new_solutions
    else:
        if paths_in_sol(solutions, paths):
            return solutions
        else:
            return []

#sol = [1,2,3,4,5,6,7,8,9,10,11]
#sol = [[1,2,3,4,5,6,7,8,9,10,11],[1,2,3,4,5,6,7,8,9,10,11,12]]
#x = [[1,2,3,4], [3,4,5,6,7,8], [1,2,3,4,5,6,7,8,9,10,11],[],[2,3,4],[4,5,6,7,9]]
#x = [[1,2,3,4], [3,4,5,6,7,8], [1,2,3,4,5,6,7,8,9,10,11],[],[2,3,4],[4,5,6,7]]
#sol = [1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 38, 39, 40]
#paths = [[-9, -8, -7, -6, -6, -5, -4, -3], [5, 6, 6, 7, 8, 9, 10, 11, 12, 13], [1, 2, 3, 4], [19, 20, 21, 22, 23, 25, 26, 27, 28], [15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26], [-16, -15, -14, -13, -12, -11, -10, -9], [9, 10], [18, 19, 20, 21, 22, 23, 25, 26], [9, 10, 11, 12, 13, 14, 15, 16, 17, 18], [-29, -28, -27, -26, -25, -23, -22, -21, -20], [10, 11, 12, 13, 14, 15, 16], [38, 39, 40], [23, 25, 26, 27], [-39, -38, -29, -28, -27, -26, -25, -23, -22], [-39, -38, -29, -28, -27, -26, -25, -23, -22, -21, -20], [9, 10, 11, 12, 13, 14, 15, 16, 17], [-9, -8, -7], [-23, -22, -21, -20, -19, -18, -17], [13, 14, 15, 16, 17, 18], [-40, -39, -38, -29, -28, -27, -26, -25], [-16, -15, -14, -13, -12, -11, -10], [14, 15, 16, 17, 18, 19, 20, 21], [6, 7, 8, 9], [8, 9, 10, 11, 12], [3, 4, 5, 6, 6, 7, 8, 9, 10, 11], [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9], [1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10], [3, 4, 5, 6, 6, 7, 8, 9, 10], [18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 38, 39], [-2, -1], [-9, -8, -7, -6, -6, -5, -4], [14, 15, 16, 17, 18, 19, 20], [4, 5, 6, 6, 7, 8], [-40, -39, -38, -29, -28], [40], [-10, -9, -8, -7, -6, -6, -5], [11, 12, 13, 14, 15], [-29, -28, -27, -26, -25, -23, -22, -21], [-13, -12, -11, -10, -9, -8, -7, -6, -6], [18, 19, 20, 21, 22], [-14, -13, -12], [19, 20, 21, 22, 23, 25], [1, 2, 3, 4, 5, 6], [21, 22, 23, 25, 26, 27, 28, 29], [29, 38, 39, 40], [27, 28, 29, 38, 39, 40], [25, 26, 27, 28, 29], [9, 10, 11], [7, 8, 9, 10, 11, 12, 13, 14, 15], [17, 18, 19, 20], [21, 22, 23, 25, 26], [7, 8, 9, 10, 11, 12, 13, 14, 15], [13, 14, 15, 16, 17, 18, 19], [20, 21, 22, 23, 25, 26, 27, 28, 29, 38, 39, 40], [-5, -4, -3, -2, -1], [15, 16, 17, 18], [21, 22, 23, 25], [1], [-15, -14, -13, -12, -11, -10, -9, -8, -7, -6], [-9, -8, -7, -6, -6, -5, -4, -3], [12, 13, 14, 15, 16, 17, 18, 19], [-26, -25, -23, -22, -21, -20, -19, -18, -17, -16, -15], [-7, -6, -6, -5], [-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6], [39, 40], [26, 27, 28, 29, 38, 39], [23, 25, 26, 27, 28, 29, 38], [17, 18, 19, 20, 21, 22], [6, 7, 8], [1, 2, 3, 4, 5, 6, 6, 7, 8], [21], [17, 18, 19, 20, 21, 22, 23], [2, 3, 4, 5], [15, 16, 17, 18, 19, 20, 21], [4, 5, 6, 6, 7], [-21, -20], [-39, -38, -29, -28, -27, -26, -25, -23], [39, 40], [28, 29], [7, 8, 9, 10, 11, 12, 13, 14, 15], [14, 15, 16, 17, 18, 19, 20], [-1], [-40], [3, 4, 5, 6, 6, 7, 8, 9, 10], [-28, -27, -26, -25, -23, -22], [15, 16, 17, 18, 19, 20], [22, 23, 25], [11, 12, 13, 14, 15, 16, 17], [40], [], [-38, -29, -28, -27, -26, -25, -23], [17, 18, 19, 20, 21, 22, 23], [6, 7, 8, 9, 10], [23, 25, 26, 27, 28, 29, 38], [6, 6, 7, 8, 9], [1, 2, 3, 4, 5, 6], [3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12], [29], [9, 10, 11, 12, 13, 14], [-6, -5, -4, -3, -2, -1], [14, 15, 16, 17, 18, 19, 20, 21, 22], [20, 21, 22, 23, 25, 26], [17, 18, 19, 20, 21, 22, 23, 25, 26, 27], [9, 10, 11, 12, 13, 14, 15, 16, 17], [-26, -25, -23, -22, -21, -20, -19, -18, -17], [8, 9, 10, 11, 12, 13, 14, 15], [22, 23, 25, 26, 27, 28, 29], [23, 25, 26, 27, 28, 29, 38, 39, 40], [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], [16, 17, 18, 19, 20, 21, 22], [-7, -6, -6, -5, -4, -3, -2, -1], [-13, -12, -11, -10, -9, -8], [21, 22, 23, 25, 26, 27, 28, 29, 38, 39, 40], [-16, -15, -14, -13, -12, -11, -10, -9], [39, 40], [-13, -12, -11, -10, -9, -8, -7, -6, -6, -5, -4, -3, -2, -1], [1, 2, 3, 4, 5, 6], [6, 6, 7, 8, 9, 10, 11], [-10, -9, -8, -7], [6, 6, 7, 8, 9, 10, 11, 12]]
#print(check_sol(sol, paths))
#print(mapping(sol, x[0]))
#print(paths_in_sol(sol, x))
#print(check_sol(sol, x))



def abs_sort_path(path):
    new_path = []
    #new_path = [abs(x) for x in path]
    for number in path:
        new_path.append(abs(number))
    new_path.sort()
    return new_path

def abs_sort_path_list(paths):
    new_paths = []
    for path in paths:
        new_paths.append(abs_sort_path(path))
    new_paths2 = []
    for path in new_paths:
        for number in path:
            new_paths2.append(number)
    new_paths2.sort()
    return new_paths2

def coverage_sol(solution):
    new_solution = abs_sort_path(solution)
    Dic_LU = {}
    for LU in new_solution:
        Dic_LU[LU]=0
    for LU in new_solution:
        Dic_LU[LU]=Dic_LU[LU]+1
    return Dic_LU

def coverage_list(paths):
    return coverage_sol(abs_sort_path_list(paths))

def coverage_list_ref(paths, reference):
    abs_paths = abs_sort_path_list(paths)
    Dic = {}
    for LU in reference:
        LU = abs(LU)
        if LU not in Dic:
            Dic[LU] = abs_paths.count(LU)
    return Dic

def coverage_list_ref_normalized(paths, reference):
    abs_paths = abs_sort_path_list(paths)
    abs_reference = abs_sort_path(reference)
    Dic = {}
    for LU in reference:
        LU = abs(LU)
        if LU not in Dic:
            Dic[LU] = abs_paths.count(LU) / abs_reference.count(LU)
    return Dic

# use LoxP_unit_count_Dict_list to count the coverage against a reference


#A=[1,2,2,11,3,5,6,7,8,8,8,9,10,11,3,-5,-22,-5,-9]
#print(coverage_sol(A))

#P = [[-6,-5,3,3,3,],[1,1,1,9,9,9],[-8,8,8,8,-9],[1,9,11,-30],[1,2,3,4,-6]]
#print(abs_sort_path_list(P))
#print(coverage_list(P))
#AAA = coverage_list_ref(P, [1,2,3,-5,6,7,8,11])
#print(AAA)

def plot_Dic(Dic):
    L_Dic = len(Dic)
    #plt.bar(*zip(*Dictionary.items()))
    plt.bar(range(L_Dic), Dic.values(), align='center')
    plt.xticks(range(L_Dic), list(Dic.keys()))
    plt.ylabel("Coverage")
    plt.xlabel("LoxP Unit")
    #plt.savefig('coverage.png')
    plt.show()
    return None


def read_length(reads):
    if reads == []:
        return reads
    if isinstance(reads[0], list):
        reads_L = [len(x) for x in reads]
        reads_L.sort()
        Dic = {}
        for R in reads_L:
            if R not in Dic:
                Dic[R] = reads_L.count(R)
        return Dic
    else:
        return reads


def divide_dictionary(Dic1, divisor):
    Dic_values = [x / divisor for x in Dic1.values()]
    Dic_keys = list(Dic1.keys())
    Dic_divided = {Dic_keys[i]: Dic_values[i] for i in range(len(Dic_keys))}
    return Dic_divided

def N50_reads(reads):
    if isinstance(reads[0], list):
        reads_L = [len(x) for x in reads]
        reads_L.sort()
        Dic_L = read_length(reads)
        Dic_L_keys = list(Dic_L.keys())
        Dic_L_values = list(Dic_L.values())

        sum_reads = 0
        for L in reads_L:
            sum_reads = sum_reads + L
        N50 = round(sum_reads / 2)
        #print("N50 =", N50)
        counter = 0
        cumulative = 0
        while cumulative < N50:
            #cumulative = cumulative + reads_L[counter]
            cumulative = cumulative + Dic_L_keys[counter] * Dic_L_values[counter]
            #print("cumulative =", cumulative)
            counter = counter + 1
        N50_reads = Dic_L_keys[counter - 1]
        L50_reads = counter - 1
        return N50_reads
    else:
        print("input is not a list of lists")
        return reads

def plot_read_length(read_length):
    plt.bar(range(len(read_length)), read_length.values(), align='center')
    plt.xticks(range(len(read_length)), read_length)
    plt.ylabel("Number of reads")
    plt.xlabel("Read length")
    plt.title("Read length distribution")
    #plt.text(10, 10, "Number of reads = 20")
    # plt.savefig('read_length.png')
    plt.show()
    plt.close()
    return None

def plot_read_length2(S_path, name="0"):
    # Process reads
    number_reads = len(S_path)
    R_L = read_length(S_path)
    #R_L_percentage = divide_dictionary(R_L, number_reads)
    # Statistics
    reads_L = [len(x) for x in S_path]
    reads_length_mean = round(statistics.mean(reads_L), 1)
    reads_length_median = statistics.median(reads_L)
    N50 = N50_reads(S_path)

    # Plot
    # Font size
    SMALL_SIZE = 15
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 24
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    Mx_value = max(R_L.values())
    plt.figure(figsize=(16, 9))
    plt.bar(R_L.keys(), R_L.values(), align='center')
    plt.xticks(range(max(R_L.keys())+1), range(max(R_L.keys())+1))
    plt.ylabel("Number of reads")
    plt.xlabel("Read length (LUs)")
    plt.title("Read length distribution")
    plt.vlines(N50, 0, Mx_value, colors="orange", linestyle="--", label="N50")
    plt.text(max(R_L.keys()) * 0.8, Mx_value * 0.9, "Reads length mean = " + str(reads_length_mean) + "\n" + "Reads length median = " + str(reads_length_median) + "\n" + "N50 = " + str(N50))
    plt.savefig("read_length_distribution_" + name + "_R_number_" + str(number_reads), dpi=300)
    plt.savefig("read_length_distribution_" + name + "_R_number_" + str(number_reads) + ".svg", format='svg', dpi=300)
    # plt.legend()
    #plt.show()
    plt.close()
    return None

#A1 = [[26, 29, 30, 31, -32, 33, 39, 41, 42, 43, 44, 1], [2, 3, 4, 37, 38, -26, -25, -24, -21, -20, -19, -18, -17, -16, -12, -6, -5, -11], [-11, 7, 8, 9, 10], [30, 31, -32, 33, 39, 41, 42, 43, 44, 44], [3, 4, 37, 38, -26, -25, -24, -21, -20, -19, -18, -17, -16, -12, -6, -5, -11, -10, -9], [-20, -19, -18, -17, -16, -12, -6, -5, -11, -10, -9, -8, -7], [30, 31, -32, 33, 39, 41, 42, 43, 44, 1, 5], [-33, 32, -31, -30, -29, -26, 28, 39, 40], [-8, -7, 11, -10, -9, -6, -5, -1, -44, -43, -42, -41], [-44, -43, -42, -41, -39, -33, 32, -31, -30, -29, -26, 28], [-1, -44, -43, -42, -41, -40, -39, -28, 26, 29, 30], [17, 18, 19, 20, 21, 24, 25, 26, -38, -37, -4, -3, -2, -1, -44, -43, -42, -41, -40, -39, -28], [32, -31, -30, -29, -26, 28, 39, 40, 41, 42, 43, 44, 44]]
#plot_read_length2(A1)

def plot_read_length_percentage(S_path, name="0"):
    # Process reads
    number_reads = len(S_path)
    R_L = read_length(S_path)
    R_L_percentage = divide_dictionary(R_L, number_reads)
    # Statistics
    reads_L = [len(x) for x in S_path]
    reads_length_mean = statistics.mean(reads_L)
    reads_length_median = statistics.median(reads_L)
    N50 = N50_reads(S_path)

    # Plot
    Mx_value = max(R_L_percentage.values())
    plt.figure(figsize=(20, 10))
    plt.bar(range(len(R_L_percentage)), R_L_percentage.values(), align='center')
    plt.xticks(range(len(R_L_percentage)), R_L_percentage)
    plt.ylabel("Number of reads in percentage")
    plt.xlabel("Read length")
    plt.title("Read length distribution")
    plt.vlines(N50, 0, Mx_value, colors="orange", linestyle="--", label="N50")
    plt.text(max(R_L_percentage.keys()) * 0.8, Mx_value * 0.9, "Number of reads = " + str(number_reads) + "\n" + "Reads length mean = " + str(round(reads_length_mean, 1)) + "\n" + "Reads length median = " + str(reads_length_median) + "\n" + "N50 = " + str(N50))
    plt.savefig("read_length_distribution_" + name + "_R_number_" + str(number_reads), dpi=200)
    # plt.legend()
    # plt.show()
    plt.close()
    return None

# Count the number of LUs in a list of reads.
def number_LU(reads):
    sum_LU = 0
    for read in reads:
        sum_LU = sum_LU + len(read)
    return sum_LU

# Plot statistics on the simulated chromosomes. E.g. Copy number (CN) of each LUs.
def count_LU_CN(Chr, Max_LU=0):
    # make all the LU positive
    new_chr = [abs(x) for x in Chr]
    # generate a reference dictionary
    LU_CN = {}
    if Max_LU < max(new_chr):
        Max_LU = max(new_chr)
    for i in range(1, Max_LU+1):
        LU_CN[i] = 0
    # count the LU
    for LU in new_chr:      # Note the LU should be already positive
        LU_CN[LU] += 1
    return LU_CN

def find_max_value(Chrs):
    MAX = 0
    for Chr in Chrs:
        new_chr = [abs(x) for x in Chr]
        MAX_temp = max(new_chr)
        if MAX_temp > MAX:
            MAX = MAX_temp
    return MAX

def count_LU_CN_multi(Chrs):
    Max_LU = find_max_value(Chrs)
    LU_CN_TOT = {}
    for Chr in Chrs:
        count_LU = count_LU_CN(Chr, Max_LU=Max_LU)
        for key, value in count_LU.items():
            if key not in LU_CN_TOT:
                LU_CN_TOT[key] = [value]
            else:
                LU_CN_TOT[key].append(value)
    return LU_CN_TOT

def calculate_LU_CN_percentage(Chrs, max_CN=5):
    num_chrs = len(Chrs)
    LU_CN = count_LU_CN_multi(Chrs)
    LU_CN_percentage = {}
    for key, value in LU_CN.items():
        # This changed all the CN bigger than max_CN into max_CN. So it can calculate the percentage of CN equal of bigger than max_CN
        CN_cleaned = []
        for CN in value:
            if CN <= max_CN:
                CN_cleaned.append(CN)
            else:
                CN_cleaned.append(max_CN)
        percentage = []
        for CN in range(max_CN+1):
            # This will store the percentage of copy number (CN) for each LU
            percentage.append(CN_cleaned.count(CN) / num_chrs)
        LU_CN_percentage[key] = percentage
    return LU_CN_percentage

def plot_LU_CN(Chrs, Plot="histogram", essential=[], CEN=[]):
    LU_CN_TOT = count_LU_CN_multi(Chrs)
    LU_CN_mean = {}
    LU_CN_SD = {}
    for key, value in LU_CN_TOT.items():
        LU_CN_mean[key] = statistics.mean(value)
        LU_CN_SD[key] = statistics.stdev(value)
    # Plot
    plt.figure(figsize=(20, 10))
    if Plot == "boxplot":
        plt.boxplot(LU_CN_TOT.values(), labels=LU_CN_TOT.keys(), showfliers=True, showmeans=True)
    elif Plot == "violinplot":
        plt.violinplot(LU_CN_TOT.values())
    else:
        plt.bar(LU_CN_mean.keys(), LU_CN_mean.values(), align='center', yerr=LU_CN_SD.values(), capsize=5)
    plt.xticks(range(1, len(LU_CN_TOT.values()) + 1))
    for esse in essential:
        plt.gca().get_xticklabels()[esse-1].set_color("red")
    if CEN != []:
        plt.gca().get_xticklabels()[CEN[0] - 1].set_color("gold")
    #plt.gca().get_xticklabels()[13].set_color("blue")
    #plt.gca().get_xticklabels()[31].set_color("blue")
    plt.ylabel("LU CN")
    plt.xlabel("LUs")
    plt.title("LoxP Unit Copy Number")
    #plt.savefig("LU_CN.png", dpi=300)
    #plt.savefig("LU_CN.svg", format='svg', dpi=300)
    plt.show()
    plt.close()
    return None

def plot_LU_CN_percentage(Chrs, max_CN=5, essential=[], CEN=[], filename="", SE=""):
    LU_CN_percentage = calculate_LU_CN_percentage(Chrs, max_CN=max_CN)
    #print("LU_CN_percentage =", LU_CN_percentage)
    LU_names = LU_CN_percentage.keys()
    CN_percentage_sorted = [[] for _ in range(max_CN+1)]
    for i in range(max_CN+1):
        for LU in LU_CN_percentage.values():
            CN_percentage_sorted[i].append(LU[i])
    Bottom = [0 for _ in range(len(LU_CN_percentage.values()))]
    #print("CN_percentage_sorted =", CN_percentage_sorted)
    Labels = [str(x) for x in range(max_CN+1)]
    Labels[-1] = "≥" + Labels[-1]
    # get discrete colormap
    # You can choose different gradients here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    Colours = cm.Reds(np.linspace(0, 1, max_CN+1))     # Colours = Blues, Reds, Greys, YlGn, YlOrBr, binary

    # Plot
    #plt.figure(figsize=(12, 6), dpi=300)
    plt.figure(figsize=(7.5, 2.75), dpi=300)    # figsize=(7.5, 3.5)
    # Font size
    SMALL_SIZE = 6
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels     # labelsize=4 for chromosomes longer than 100
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    for i in range(max_CN+1):
        plt.bar(LU_names, CN_percentage_sorted[i], bottom=Bottom, label=Labels[i], color=Colours[i])  # color=Colours[i]
        # This sums the two lists bottom and the last values CN_percentage_sorted[i]
        Bottom = [x + y for x, y in zip(Bottom, CN_percentage_sorted[i])]
    #plt.bar(LU_CN_mean.keys(), LU_CN_mean.values(), align='center', yerr=LU_CN_SD.values(), capsize=5)
    plt.xticks(range(1, len(LU_CN_percentage.values()) + 1))
    for esse in essential:
        plt.gca().get_xticklabels()[esse-1].set_color("red")
    if CEN != []:
        plt.gca().get_xticklabels()[CEN[0] - 1].set_color("gold")
    #plt.gca().get_xticklabels()[13].set_color("blue")
    #plt.gca().get_xticklabels()[31].set_color("blue")
    plt.ylabel("Percentage of LU CN")
    plt.xlabel("LUs")
    #plt.xticks(rotation=90)    # Use this for chromosomes longer than 100
    n_SE = ""
    if SE != "":
        n_SE = ". SE = " + str(SE)
    plt.title("Percentage of LU CN" + n_SE)
    plt.legend(loc=3)
    if filename != "":
        plt.savefig("SCRaMbLE_evolution_percentage_LU/percentage_LU_CN_" + filename + ".png", dpi=300, bbox_inches='tight')
        plt.savefig("SCRaMbLE_evolution_percentage_LU/percentage_LU_CN_" + filename + ".svg", format='svg', dpi=300, bbox_inches='tight')
    else:
        plt.savefig("SCRaMbLE_evolution_percentage_LU/percentage_LU_CN.png", dpi=300, bbox_inches='tight')
        plt.savefig("SCRaMbLE_evolution_percentage_LU/percentage_LU_CN.svg", format='svg', dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    return None

# SCRaMbLEs many synthetic chromosomes and plots how the LU CN change.
def SCRaMbLE_SIM_LU_CN(syn_chr, events=100, simulations=1000, essential=[], CEN=[], circular=False, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], max_CN=5):
    steps = 5
    snapshots = range(0, events+1, steps)
    snapshots_SCRaMbLEd = []
    SCRaMbLEd_chrs = [syn_chr[:] for _ in range(simulations)]
    for E in range(0, events+1, 1):
        print(E)
        if E in snapshots:
            snapshots_SCRaMbLEd.append(SCRaMbLEd_chrs[:])
        for s in range(simulations):
            # Perform SCRaMbLE on the synthetic chromosome
            SCRaMbLEd_chrs[s] = force_SCRaMLE_lin_cir(SCRaMbLEd_chrs[s], 1, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
    # Plot
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
    # Do the actual plotting
    for i in range(len(snapshots_SCRaMbLEd)):
        filename = lin_cir + "_chr_L" + str(len(syn_chr)) + "_" + random_seed + "_sim" + str(simulations) + "_P" + probability_str + "_SE" + str(snapshots[i])
        plot_LU_CN_percentage(snapshots_SCRaMbLEd[i], max_CN=max_CN, essential=essential, CEN=CEN, filename=filename, SE=str(snapshots[i]))
    return None

# test the code
if __name__ == "__main__":

    B = [[3, 3, 1, 1, 2, 3, 6, 4, 5, 6], [1, 1, 3, 5, 6, 6, 7], [3, 3, 5, 6, 1, 1], [5, 7, 3, 3, 1, 1, 7]]
    # print(count_LU_CN(B[0], Max_LU=0))
    # print(count_LU_CN_multi(B))
    # print(calculate_LU_CN_percentage(B))
    # plot_LU_CN(B, Plot="boxplot")
    #plot_LU_CN_percentage(B, max_CN=5)

    syn_chr = list(range(1, 45, 1))
    essential = [2, 7, 9, 10, 12, 19, 20, 24]  # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.

    SCRaMbLE_SIM_LU_CN(syn_chr, events=100, simulations=1000, essential=essential, CEN=[2], circular=True, mu=0, sigma=10, force=True, probability=[0, 2, 2, 1], max_CN=5)

    syn_chr = list(range(1, 101, 1))
    essential = [50]
    #SCRaMbLE_SIM_LU_CN(syn_chr, events=1000, simulations=500, essential=essential, CEN=[50], circular=True, mu=0,sigma=10, force=True, probability=[0, 2, 2, 2], max_CN=8)

    # I use the following website to create the giff: https://ezgif.com/maker

    """
    from SCRaMbLE_DNA_simulation import DNA_extraction
    S = list(range(1, 45, 1))
    number_reads = 1000000
    S_path = DNA_extraction(S, number_reads)
    #print("S_path =", S_path)
    N50 = N50_reads(S_path)
    print("N50 =", N50)
    R_L = read_length(S_path)
    R_L_values = [x/number_reads for x in R_L.values()]
    #R_L_keys = [x for x in R_L.keys()]
    R_L_keys = list(R_L.keys())
    R_L_percentage= {R_L_keys[i]: R_L_values[i] for i in range(len(R_L_keys))}
    print("R_L =", R_L)
    print("R_L_percentage =", R_L_percentage)
    plot_read_length(R_L)
    plot_read_length(R_L_percentage)
    plot_read_length_percentage(S_path)
    """

    """
    from SCRaMbLE_simulation_3 import SCRaMbLE4
    from SCRaMbLE_DNA_simulation import DNA_extraction_coverage
    from comparison_sol import half_pos_one
    from correction_MM import clean_paths

    segments = 44  #number of loxP segments
    SCRaMbLEd_chr = list(range(1, segments, 1))
    essential = [2,7,9,10,12,20]
    SCRaMbLEd_events = 20

    for simulations in range(20):
        S = SCRaMbLE4(SCRaMbLEd_chr, SCRaMbLEd_events,essential)
        S = half_pos_one(S)
        print("CHROMOSOME =",S)
        S_path_25 = DNA_extraction_coverage(S, 50)
        dic2 = coverage_list_ref_normalized(S_path_25, S)
        plot_Dic(dic2)
        cleaned = clean_paths(S_path_25, 5)
        print()
        print(S_path_25)
        print(cleaned)
        print()
        print("-------------------")
    """