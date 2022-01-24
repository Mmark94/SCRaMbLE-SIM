import random
from comparison_sol import invert
from SCRaMbLE_simulation_3_circular import smallest_SV
import numpy as np
import matplotlib.pyplot as plt

#segments = 100  #number of loxP segments

#SCRaMbLEd_chr = list(range(1, segments, 1))

#SCRaMbLEd_chr = [1, 2, 3, 4, 5, 6, -8, -7, 9, 10, 11, 12, -8, -7, 9, 10, 11, 12, 13, 14, 15, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 43, 44, 2, 3, 4, 5, 6, -8, -7, 9, 10, 11, 12, 13, 14, 15, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]

# 1 LoxP unit = 2.9Kb


def DNA_extraction_any_L(syn_chr, coverage):
    DNA_fragments = []
    T=0
    while T < coverage:
        pos1 = random.randrange(0, len(syn_chr))
        pos2 = random.randrange(0, len(syn_chr))
        if pos1 > pos2:
            temp = pos1
            pos1 = pos2
            pos2 = temp
        if pos1 == pos2:
            continue
        DNA_fragments.append(syn_chr[pos1:pos2])
        T += 1
    return DNA_fragments

#print(DNA_extraction_any_L(SCRaMbLEd_chr , 20))

# I should create a function for DNA extraction of circular chromosomes
def DNA_extraction(syn_chr, reads, mu=8, sigma=3, circular=False):
    DNA_fragments = []
    if len(syn_chr)<4:
        return [syn_chr]
    T=0
    while T < reads:
        pos1 = random.randrange(0, len(syn_chr))
        #pick a random number to decide the length of the fragment
        #R = random.choices(range(3, 10), [1, 2, 3, 4, 3, 2, 1], k=1)[0]
        R = int(random.gauss(mu, sigma))
        if R<1:
            R=1
        if random.random() > 0.5:
            pos2 = pos1 + R
        else:
            pos2 = pos1 - R
        if pos1 == pos2:
            continue
        # if the chromosome is circular
        if circular:
            if pos2 > len(syn_chr):
                pos2 = pos2 - len(syn_chr)
            if pos2 < 0:
                pos2 = len(syn_chr) + pos2
        else:
            if pos2 > len(syn_chr):
                pos2 = len(syn_chr)
            if pos2 < 0:
                pos2 = 0
        if pos1 > pos2:
            temp = pos1
            pos1 = pos2
            pos2 = temp
        # With thw two position it will create the fragment or read
        #if not(circular) and not(smallest_SV(pos1, pos2, syn_chr)):
        #    fragment = syn_chr[pos1:pos2]
        if circular:
            if smallest_SV(pos1, pos2, syn_chr):
                fragment = syn_chr[pos1:pos2]
            else:
                fragment = syn_chr[pos2:] + syn_chr[:pos1]
        else:
            fragment = syn_chr[pos1:pos2]
        # This will invert the fragment with probability 30%
        P = random.random()
        if P>0.3:
            DNA_fragments.append(fragment)
        else:
            DNA_fragments.append(invert(fragment))
        T += 1
    return DNA_fragments

#syn_chr = list(range(1, 45, 1))
#print(DNA_extraction(syn_chr , 20))
#A = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#reads = DNA_extraction(A, 8, mu=8, sigma=3, circular=True)
#print(reads)

def DNA_extraction_coverage(syn_chr, coverage, reads_len=8, sigma=3, circular=False):
    #The function firstly find out how many reads are necessary to have one fold coverage using the formula: NL/G=coverage   (N=number of reads, L=read length, G=genome size)
    if len(syn_chr) <= reads_len:
        reads_num = coverage
    else:
        reads_num = int(len(syn_chr) / reads_len)
    DNA_fragments = DNA_extraction(syn_chr, reads_num * coverage, reads_len, sigma, circular)
    return DNA_fragments


def coverage(paths):
    Dic = {}
    if isinstance(paths[0], int):
        for path in paths:
            if path not in Dic:
                print(path)
                Dic[path] = 1
            else:
                Dic[path] = Dic[path]+1
        return Dic
    else:
        for path2 in paths:
            for path in path2:
                if path not in Dic:
                    Dic[path] = 1
                else:
                    Dic[path] = Dic[path] + 1
        return Dic

#L=[1,2,3,4,3]
#L2=[[1,2,6,3,4,3],[3,3,2,1,5,6]]

#print(coverage(L))
#print(coverage(L2))

def coverage2(paths, solution):
    Dic = {}
    for elem in solution:
        elem = abs(elem)
        if elem not in Dic:
            Dic[elem] = 0
    if isinstance(paths[0], int):
        for path in paths:
            path = abs(path)
            if path not in Dic:
                print(path)
                Dic[path] = 1
            else:
                Dic[path] = Dic[path]+1
        return Dic
    else:
        for path2 in paths:
            for path in path2:
                path = abs(path)
                if path not in Dic:
                    Dic[path] = 1
                else:
                    Dic[path] = Dic[path] + 1
        return Dic

#solution=[1,2,3,4,5,6]

#print(coverage2(L2,solution))

#myDictionary = coverage2(L2,solution)
# look also in Mapping_coverage_MM.py
def plot_Dic(myDictionary):
    pos = np.arange(len(myDictionary.keys()))
    width = 1.0     # gives histogram aspect to the bar diagram

    ax = plt.axes()
    ax.set_xticks(pos + (width / 2))
    ax.set_xticklabels(myDictionary.keys())

    plt.bar(myDictionary.keys(), myDictionary.values(), width, color='g')
    plt.show()


def plot_Dic2(myDictionary):
    # Plot
    plt.figure(figsize=(20, 10))
    plt.bar(myDictionary.keys(), myDictionary.values(), align='center')
    #plt.xticks(range(max(myDictionary.keys()) +1), range(max(myDictionary.keys())+1))
    plt.ylabel("Number of LU")
    plt.xlabel("LU")
    plt.title("Coverage")
    #plt.text(max(R_L.keys()) * 0.8, Mx_value * 0.9,"Reads length mean =" + str(reads_length_mean) + "\n" + "Reads length median =" + str(reads_length_median) + "\n" + "N50 =" + str(N50))
    plt.savefig("Coverage.png", dpi=200)
    # plt.legend()
    plt.show()

# This function is from Mapping_coverage_MM.py
def divide_dictionary(Dic1, divisor):
    Dic_values = [x / divisor for x in Dic1.values()]
    Dic_keys = list(Dic1.keys())
    Dic_divided = {Dic_keys[i]: Dic_values[i] for i in range(len(Dic_keys))}
    return Dic_divided

# test the code
if __name__ == "__main__":
    syn_chr = list(range(1, 101, 1))
    coverage = 100000
    #DNA = DNA_extraction_coverage(syn_chr, coverage, reads_len=8, sigma=3, circular=True)
    DNA = DNA_extraction(syn_chr, coverage, mu=8, sigma=3, circular=False)
    cov = coverage2(DNA, syn_chr)
    cov_percentage = divide_dictionary(cov, coverage)
    print("cov =", cov)
    print("cov_percentage =", cov_percentage)
    #plot_Dic(cov)
    plot_Dic2(cov)
    plot_Dic2(cov_percentage)