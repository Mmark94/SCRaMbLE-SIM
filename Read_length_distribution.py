from SCRaMbLE_DNA_simulation import DNA_extraction_coverage
from Mapping_coverage_MM import read_length
from Mapping_coverage_MM import plot_read_length
from Mapping_coverage_MM import divide_dictionary
from Mapping_coverage_MM import N50_reads
import matplotlib.pyplot as plt



syn_chr = list(range(1, 100, 1))
number_reads = 50000     #coverage
reads_len = 30       # mu
sigma = 3

S_path = DNA_extraction_coverage(syn_chr, number_reads, reads_len, sigma)
R_L = read_length(S_path)
R_L_percentage= divide_dictionary(R_L, number_reads)
N50 = N50_reads(S_path)
print(R_L)
print(R_L_percentage)
print("N50 =", N50)

# Plot
#plot_read_length(R_L)
#plot_read_length(R_L_percentage)
Mx_value = max(R_L_percentage.values())
plt.figure(figsize=(20,10))
plt.bar(R_L_percentage.keys(), R_L_percentage.values(), align='center')
plt.xticks(range(max(R_L_percentage.keys())), range(max(R_L_percentage.keys())))
#plt.bar(range(len(R_L_percentage)), R_L_percentage.values(), align='center')
#plt.xticks(range(len(R_L_percentage)), R_L_percentage)
plt.ylabel("Number of reads in percentage")
plt.xlabel("Read length")
plt.title("Read length distribution")
plt.vlines(N50, 0, Mx_value, colors="orange", linestyle="--", label="N50")
plt.text(max(R_L_percentage.keys())*0.8, Mx_value*0.9, "reads length (mu) =" + str(reads_len) + "\n" + "reads length (sigma) =" + str(sigma) + "\n" + "N50 =" + str(N50))
#plt.savefig("read_length_distribution_R_mu"+str(reads_len)+"_R_sigma"+str(sigma)+"_sim_"+str(number_reads),dpi=200)
#plt.legend()
plt.show()