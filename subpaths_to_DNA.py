import pandas as pd
from comparison_sol import str_path
from comparison_sol import path_str
from comparison_sol import solution_to_bps_size
from Visualize_compare_paths import save_paths_in_file
import subprocess
import shlex
import os
import shutil
import signal

# Save and export the simulated subpaths into a file. Each subpath is a different row.
def subpath_to_file(subpaths, file_name="", solution=[]):
    if subpaths == []:
        return None
    if isinstance(subpaths[0], int):    # There is only one sequence in "subpaths"
        f = open("subpaths_SIM/subpaths_" + file_name + ".txt", "w+")
        f.write("single_read_" + path_str(subpaths) + "\t" + path_str(subpaths))
        f.close()
        return None
    else:
        f = open("subpaths_SIM/subpaths_" + file_name + ".txt", "w+")
        counter = 1
        for subpath in subpaths:
            if subpath == []:
                continue
            if counter < len(subpaths):
                # Each line is going to be: read_ID_subpath, tab and subpath
                # I record the solution of the paths in the first path
                if counter == 1 and solution != []:
                    f.write("read_" + str(counter) + "__" + path_str(subpath) + "__solution_" + path_str(solution) + "\t" + path_str(subpath) + "\n")
                else:
                    f.write("read_" + str(counter) + "__" + path_str(subpath) + "\t" + path_str(subpath) + "\n")
            else:
                # This is the last read, so do not add "\n"
                f.write("read_" + str(counter) + "__" + path_str(subpath) + "\t" + path_str(subpath))
            counter += 1
        f.close()
        return None

# Save and export the simulated subpaths into DNA sequences.
def subpath_to_DNA(file_name="", syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", Linux=True, debug=False, Paths=True):
    if Paths == True:
        file_name_s = "subpaths_" + file_name
    else:
        file_name_s = file_name
    path_file = "subpaths_SIM/" + file_name_s + ".txt"

    # Convert the subpath file into a DNA file file.txt => file.reseq.fasta using the script restructseq.v1.01.pl
    # perl restructseq.v1.01.pl IXR_BACnewseq.fa IXR_BACnewseq.loxpreg Test2.txt -spid Test2 -outdir ./
    perl_command = "perl ./restructseq.v1.01.pl " + syn_DNA + " " + syn_DNA_loxP + " " + path_file + " -spid " + file_name_s + " -outdir " + "subpaths_SIM/"
    if debug:
        print(perl_command)
    if Linux:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    else:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE)
    perl_script.communicate()   # I am not sure of the function of .communicate()
    if not debug:
        os.remove(path_file)
    return None

# convert a solution in a fasta format. For synIII use syn_DNA="yeast_chr03_9_02.fa", syn_DNA_loxP="yeast_chr03_9_02.loxpreg"
def solution_to_DNA(solution=[], paths_name="A_path", file_name="", syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", Linux=True, debug=False):
    save_paths_in_file(paths=solution, paths_name=paths_name, file_name=file_name, PATH="subpaths_SIM/")
    subpath_to_DNA(file_name=file_name, syn_DNA=syn_DNA, syn_DNA_loxP=syn_DNA_loxP, Linux=Linux, debug=debug, Paths=False)
    return None

# Use Canu to assemble the reads
def DNA_to_contig_Canu(file_name="", genome_size="100k", Linux=True, debug=False):
    # ~/canu-2.1.1/bin/canu -p "${!sample}"_canu -d "$path2"/"${!sample}"-Canu genomeSize=12.5m maxThreads=$NSLOTS mhapThreads=$NSLOTS maxMemory=$((4*NSLOTS)) useGrid=false -nanopore "$path1"/"${!sample}".fastq correctedErrorRate=0.13 minReadLength=5000 minOverlapLength=4000 stopOnLowCoverage=4 minInputCoverage=4
    # assembly_command = "canu -p " + file_name + "_canu -d " + file_name + "-Canu genomeSize=100k maxThreads=$NSLOTS mhapThreads=$NSLOTS maxMemory=$((4*NSLOTS)) useGrid=false -nanopore " + file_name + ".fastq correctedErrorRate=0.13 minReadLength=5000 minOverlapLength=4000 stopOnLowCoverage=4 minInputCoverage=4"
    assembly_command = "~/canu-2.1.1/bin/canu -p " + file_name + "_canu -d subpaths_SIM/" + file_name + "-Canu genomeSize=" + genome_size + " -nanopore subpaths_SIM/subpaths_" + file_name + ".reseq.fasta correctedErrorRate=0.13 minReadLength=5000 minOverlapLength=4000 stopOnLowCoverage=0.1 minInputCoverage=0.1"
    copy_command = "cp subpaths_SIM/" + file_name + "-Canu/" + file_name + "_canu.contigs.fasta subpaths_SIM/" + file_name + "_canu.contigs.fasta"
    if debug:
        print(assembly_command)
        print(copy_command)
    if Linux:
        assembly_process = subprocess.Popen(assembly_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        assembly_process.communicate()
        copy_process = subprocess.Popen(copy_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        copy_process.communicate()
    if not debug:
        shutil.rmtree("subpaths_SIM/" + file_name + "-Canu")
    return None

# Use Shasta to assemble the reads
def DNA_to_contig_Shasta(file_name="", Linux=True, debug=False):
    assembly_command = "~/shasta-Linux-0.8.0 --input subpaths_SIM/subpaths_" + file_name + ".reseq.fasta --config Nanopore-Oct2021 --threads 8 --Reads.minReadLength 10000 --assemblyDirectory subpaths_SIM/" + file_name + "-Shasta"
    copy_command = "cp subpaths_SIM/" + file_name + "-Shasta/Assembly.fasta subpaths_SIM/" + file_name + "_shasta.contigs.fasta"
    if debug:
        print(assembly_command)
        print(copy_command)
    if Linux:
        assembly_process = subprocess.Popen(assembly_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        assembly_process.communicate()
        copy_process = subprocess.Popen(copy_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        copy_process.communicate()
    if not debug:
        shutil.rmtree("subpaths_SIM/" + file_name + "-Shasta")
    return None

# Convert the assembled contig into a path
def contig_to_path(file_name="", syn_DNA_loxP_DNA="IXR_BACnewseq.loxpreg.fa", Linux=True, debug=False):
    # grapsubpath_test.sh
    # minimap2 -N 100 -t 8 -ax map-ont $f IXR_BACnewseq.loxpreg.fa > ${f%.*}.loxpreg.minimap.sam
    mapping_command = "~/minimap2/minimap2 -N 100 -t 8 -ax map-ont subpaths_SIM/" + file_name + "_canu.contigs.fasta " + syn_DNA_loxP_DNA + " > subpaths_SIM/" + file_name + ".loxpreg.minimap.sam"
    # perl grap_subpath.v1.02.pl -prefix ${f%.*} -lrds $f -luseq "$reference1" -aligns ${f%.*}.loxpreg.minimap.sam -gap 100 -mapqual 0 -coverage 0.5 -accuracy 0.2 -lrdscvg 0.5 -outdir .
    subpath_command = "perl grap_subpath.v1.02.pl -prefix " + file_name + " -lrds subpaths_SIM/" + file_name + "_canu.contigs.fasta -luseq " + syn_DNA_loxP_DNA + " -aligns subpaths_SIM/" + file_name + ".loxpreg.minimap.sam -gap 100 -mapqual 0 -coverage 0.5 -accuracy 0.2 -lrdscvg 0.5 -outdir subpaths_SIM/"
    if debug:
        print(mapping_command)
        print(subpath_command)
    if Linux:
        map_process = subprocess.Popen(mapping_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        map_process.communicate()
        subpath_calling_command = subprocess.Popen(subpath_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
        subpath_calling_command.communicate()
    if not debug:
        os.remove("subpaths_SIM/" + file_name + ".loxpreg.minimap.sam")
        os.remove("subpaths_SIM/" + file_name + ".subpath.info")
        os.remove("subpaths_SIM/" + file_name + ".subpath.log")
    return None

# This is similar to the function extract_subpaths() in get_sol8.py
# Read the assembled path
def import_path_contig(file_name=""):
    # the file name is usually ID_name + ".subpath.lst" or ID_name + ".subpath.merge.lst"
    # open the file with read_csv
    file = pd.read_csv("subpaths_SIM/" + file_name + ".subpath.lst", sep='\t', header=None)
    # extract all the subpath in one list
    LIST = []
    for read in file[2]:
        #print("read =", read)
        if read != read:    # This excludes the NaN values.
            continue
        LIST.append(str_path(read))
    return LIST

# From subpaths generate a file and convert it to fasta and assemble it with Canu and import the solution.
def SIM_subpath_DNA_Canu_solution(subpaths, file_name="", solution=[], syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", syn_DNA_loxP_DNA="IXR_BACnewseq.loxpreg.fa", genome_size="100k", Linux=True, debug=False, true_sol_size=False):
    if solution != []:
        solution_size = solution_to_bps_size(solution)
    # Save and export the simulated subpaths into DNA sequences.
    subpath_to_file(subpaths=subpaths, file_name=file_name, solution=solution)
    # Save and export the simulated subpaths into DNA sequences.
    subpath_to_DNA(file_name=file_name, syn_DNA=syn_DNA, syn_DNA_loxP=syn_DNA_loxP, Linux=Linux, debug=debug)
    # Use Canu to assemble the reads
    if true_sol_size:
        DNA_to_contig_Canu(file_name=file_name, genome_size=solution_size, Linux=Linux, debug=debug)
    else:
        DNA_to_contig_Canu(file_name=file_name, genome_size=genome_size, Linux=Linux, debug=debug)
    # Convert the assembled contig into a path
    contig_to_path(file_name=file_name, syn_DNA_loxP_DNA=syn_DNA_loxP_DNA, Linux=Linux, debug=debug)
    # Read the assembled path
    Solution = import_path_contig(file_name=file_name)
    return Solution


# From subpaths generate a file and convert it to fasta and assemble it with Canu and import the solution.
def SIM_subpath_DNA_Canu_solution_steps(subpaths=[], file_name="", solution=[], syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", syn_DNA_loxP_DNA="IXR_BACnewseq.loxpreg.fa", genome_size="100k", Linux=True, debug=False, true_sol_size=False, NO_subpath_to_DNA=False, NO_Canu=False, only_save_DNA=False):
    if not(NO_subpath_to_DNA):
        if solution != []:
            solution_size = solution_to_bps_size(solution)
        # Save and export the simulated subpaths into DNA sequences.
        subpath_to_file(subpaths=subpaths, file_name=file_name, solution=solution)
        # Save and export the simulated subpaths into DNA sequences.
        subpath_to_DNA(file_name=file_name, syn_DNA=syn_DNA, syn_DNA_loxP=syn_DNA_loxP, Linux=Linux, debug=debug)
        if only_save_DNA:
            return None
    if not(NO_Canu):
        # Use Canu to assemble the reads
        if true_sol_size:
            DNA_to_contig_Canu(file_name=file_name, genome_size=solution_size, Linux=Linux, debug=debug)
        else:
            DNA_to_contig_Canu(file_name=file_name, genome_size=genome_size, Linux=Linux, debug=debug)
    # Convert the assembled contig into a path
    contig_to_path(file_name=file_name, syn_DNA_loxP_DNA=syn_DNA_loxP_DNA, Linux=Linux, debug=debug)
    # Read the assembled path
    Solution = import_path_contig(file_name=file_name)
    return Solution

# test the code
if __name__ == "__main__":

    Paths = [[1,2,3], [2,3,4], [3,3,5,6], [3], [1,2,3]]
    Paths_test = [[1,2,3,3], [2,3,3,-5], [3,-5,-4,6], [-5,-4,6,7], [3,3], [3,-4,-5], [1,2,3], [2,3], [-4,6,7]]
    Paths_test = [[1, 2, 3, 3. -5], [2, 3, 3, -5. -4, 6], [3, -5, -4, 6, 7], [-5, -4, 6, 7], [3, 3], [3, -4, -5], [1, 2, 3], [2, 3], [-4, 6, 7]]
    Solution_test = [1,2,3,3,-5,-4,6,7] # 22810 bp

    #paths = [[1,2,3,4,5,6,7,8,9,10], [-2,-3,4,5,-4,-3,-2,-1,6,7,8,8,8,-7]]
    #solution_to_DNA(solution=paths[0], paths_name="AAA_path", file_name="TESTTTt", syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", Linux=False, debug=False)
    #Paths = [1,2,3,4,5,6,7]
    #subpath_to_file(Paths_test, file_name="S_test", solution=[1,2,3,3,-5,-4,6,7])
    #subpath_to_DNA(file_name="S_test", debug=True)
    #DNA_to_contig_Canu("S_test")
    #DNA_to_contig_Shasta(file_name="S_test", Linux=True, debug=True)
    #contig_to_path("S_test")
    #print(import_path_contig(file_name="S_test"))

    #print(SIM_subpath_DNA_Canu_solution(subpaths=Paths, file_name="subpaths_test"))
    #print(SIM_subpath_DNA_Canu_solution(Paths_test, file_name="S_test", solution=[1,2,3,3,-5,-4,6,7], syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", syn_DNA_loxP_DNA="IXR_BACnewseq.loxpreg.fa", Linux=True, debug=True))

    #################

    JS96 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    JS96_paths = [[-21, -20, -19], [-11, -10, -9, -8, -7, -6], [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [1, 2, 3, 4, 5, 6, 7, 8, 9], [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10], [1, 2], [-15, -14, -13, -12, -11], [35, 36, 37, 38, 39], [34, 35], [21, 22, 23, 24, 25, 26, 27, 28], [-15, -14, -13, -12, -11, -10, -9], [-23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11], [31, 32, 33, 34, 35], [-21, -20, -19, -18, -17], [1, 2, 3, 4], [6, 7, 8, 9, 10], [21, 22, 23, 24, 25, 26, 27, 28], [17, 18, 19], [9], [26, 27, 28, 29], [9, 10, 11, 12, 13, 14, 15, 16, 17, 18], [-31, -30, -29, -28, -27, -26, -25, -24, -23], [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33], [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44], [-27, -26, -25, -24, -23, -22], [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34], [11, 12, 13], [26, 27, 28, 29, 30, 31, 32], [18, 19, 20, 21, 22], [23, 24, 25, 26, 27, 28, 29], [1, 2, 3, 4, 5, 6], [-42, -41, -40, -39, -38], [23, 24, 25, 26, 27, 28, 29, 30, 31], [], [42, 43, 44], [5, 6, 7, 8, 9, 10, 11, 12, 13], [28, 29, 30, 31, 32, 33, 34, 35], [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11], [28, 29, 30, 31, 32, 33], [-35, -34, -33, -32, -31, -30, -29, -28], [27, 28, 29, 30, 31], [-22, -21, -20, -19, -18, -17, -16, -15], [], [22, 23, 24, 25, 26, 27, 28], [-33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22], [39, 40, 41, 42, 43, 44], [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], [1, 2, 3, 4, 5, 6], [31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42], [-28, -27, -26, -25, -24, -23, -22], [17, 18], [-11, -10, -9, -8, -7], [-19, -18, -17, -16, -15, -14, -13, -12], [-4, -3, -2, -1], [21, 22, 23, 24, 25], [8, 9, 10, 11], [39, 40, 41, 42, 43, 44], [1, 2, 3, 4, 5], [35, 36, 37, 38, 39, 40], [-39, -38, -37, -36, -35, -34, -33, -32, -31], [33, 34, 35, 36, 37, 38, 39], [27, 28, 29], [17, 18, 19, 20, 21, 22, 23, 24], [9, 10, 11, 12, 13, 14, 15, 16], [37, 38, 39, 40, 41, 42, 43], [9, 10, 11, 12, 13, 14, 15], [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [-11, -10, -9, -8, -7, -6, -5, -4, -3, -2], [1, 2, 3, 4, 5, 6], [10, 11, 12, 13, 14, 15], [7, 8], [13, 14, 15, 16, 17], [-23, -22, -21], [21, 22, 23, 24, 25, 26, 27, 28], [1, 2, 3, 4], [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [29, 30, 31, 32, 33, 34], [-43, -42, -41, -40, -39, -38, -37, -36, -35, -34, -33], [26, 27, 28, 29], [10, 11, 12, 13, 14, 15, 16, 17, 18], [33, 34, 35, 36, 37, 38, 39, 40, 41], [-40, -39, -38, -37], [1, 2, 3, 4, 5, 6, 7, 8, 9], [44], [3, 4, 5, 6, 7, 8, 9, 10, 11], [-25, -24, -23, -22, -21, -20, -19], [40, 41, 42, 43, 44], [6, 7, 8, 9, 10, 11, 12, 13, 14], [-13, -12, -11, -10, -9, -8, -7, -6, -5], [18, 19, 20, 21, 22, 23, 24, 25], [1, 2, 3, 4, 5], [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2], [-26, -25, -24, -23, -22, -21, -20], [6, 7, 8, 9, 10], [14, 15, 16, 17, 18, 19, 20, 21, 22], [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37], [-16, -15, -14, -13, -12, -11], [23, 24, 25, 26, 27, 28, 29], [-22, -21, -20], [-33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23], [-23, -22, -21, -20, -19, -18, -17], [8, 9, 10, 11, 12, 13], [34, 35, 36, 37, 38, 39, 40, 41], [1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [41, 42, 43, 44], [-33, -32, -31], [11, 12, 13, 14, 15], [9, 10, 11, 12, 13, 14, 15, 16, 17], [36, 37, 38, 39, 40, 41, 42, 43], [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], [-13, -12, -11, -10, -9, -8, -7, -6], [-30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18], [-30, -29, -28, -27, -26], [37, 38, 39, 40, 41, 42], [23, 24, 25, 26, 27], [7, 8, 9, 10, 11, 12, 13, 14], [1, 2, 3, 4, 5, 6], [-5, -4, -3, -2, -1], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30], [-42, -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29], [36, 37, 38, 39, 40, 41], [-31, -30, -29, -28, -27, -26, -25], [-30, -29, -28, -27, -26, -25, -24, -23], [28, 29, 30, 31, 32, 33, 34, 35, 36, 37], [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], [24, 25, 26, 27, 28, 29, 30], [27, 28, 29, 30], [20, 21, 22, 23, 24, 25], [24, 25, 26, 27, 28, 29, 30, 31, 32], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [-14], [-16, -15, -14, -13], [35, 36, 37], [10], [], [18, 19, 20, 21, 22, 23, 24, 25], [-39, -38, -37, -36, -35, -34, -33, -32, -31, -30], [-8, -7, -6, -5, -4, -3, -2], [36, 37, 38, 39, 40, 41], [4, 5, 6, 7, 8, 9], [-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5], [35, 36, 37, 38, 39, 40], [12, 13, 14], [1, 2, 3, 4], [24, 25, 26], [20, 21, 22], [1], [12, 13, 14, 15, 16, 17, 18]]

    #print(SIM_subpath_DNA_Canu_solution(JS96_paths, file_name="JS96_test", syn_DNA="IXR_BACnewseq.fa", syn_DNA_loxP="IXR_BACnewseq.loxpreg", syn_DNA_loxP_DNA="IXR_BACnewseq.loxpreg.fa", Linux=True, debug=True))

    #subpath_to_file(JS96_paths, file_name="JS96", solution=JS96_paths)
    #subpath_to_DNA(file_name="JS96", debug=True)
    #DNA_to_contig_Shasta(file_name="JS96", Linux=True, debug=True)
