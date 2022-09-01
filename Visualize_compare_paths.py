import pandas as pd
from comparison_sol import str_path
from comparison_sol import path_str
from comparison_sol import solution_to_bps_size
import subprocess
import shlex
import os
import shutil
import signal

def save_paths_in_file(paths, paths_name="A_path", file_name="A_test", PATH="./"):
    if paths == [] or paths == [[]] or paths == "[]" or paths == "[[]]":
        return None
    if isinstance(paths, str):
        paths = str_path(paths)
    if isinstance(paths_name, str):
        paths_name = [paths_name]
    file = open(PATH + file_name + ".txt", "w+")
    if isinstance(paths[0], list):
        for i in range(len(paths)):
            if paths_name == ["A_path"]:
                path_name = "path_" + str(i)
            else:
                path_name = paths_name[i]
            if isinstance(paths[i], str):
                sol_MM = path_str(paths[i])
            else:
                sol_MM = path_str(paths[i])
            file.write(path_name + "\t" + sol_MM + "\n")
    else:
        file.write(paths_name[0] + "\t" + path_str(paths) + "\n")
    file.close()
    return None

def arrow_plot_csv(file_name="", Linux=True, debug=False):
    perl_command = "perl ./bin/SCBINgrams.pl ./bin/linearGradient.rgb.tsv " + file_name + ".csv -prefix " + file_name + " -outdir arrow_plots"
    if debug:
        print(perl_command)
    if Linux:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    else:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE)
    perl_script.communicate()   # I am not sure of the function of .communicate()
    return None

def arrow_plot(paths, paths_name="A_path", file_name="", Linux=True, debug=False):
    # save paths in a txt file
    save_paths_in_file(paths, paths_name=paths_name, file_name="Arrows_" + file_name)
    # Generate the arrow plot
    perl_command = "perl ./bin/SCBINgrams.pl ./bin/linearGradient.rgb.tsv Arrows_" + file_name + ".txt -prefix " + file_name + " -outdir arrow_plots"
    if debug:
        print(perl_command)
    if Linux:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    else:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE)
    perl_script.communicate()   # I am not sure of the function of .communicate()
    return None

def dot_plot(path1, path2, paths_name=[""], Linux=True, debug=False):
    if paths_name == [""] or paths_name == "":
        paths_name = ["path1", "path2"]
    if paths_name[0] == paths_name[1]:
        paths_name[1] = paths_name[1] + "_2"
    # save path1 and path2 in a txt file
    save_paths_in_file(path1, paths_name=paths_name[0], file_name="A_" + paths_name[0])
    save_paths_in_file(path2, paths_name=paths_name[1], file_name="A_" + paths_name[1])
    # generate the dot plot
    perl_command = "perl ./bin/SCBdotplot_compr.pl A_" + paths_name[0] + ".txt A_" + paths_name[1] + ".txt -prefix dot_plot -outdir arrow_plots"
    if debug:
        print(perl_command)
    if Linux:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    else:
        perl_script = subprocess.Popen(perl_command, stdin=subprocess.PIPE)
    perl_script.communicate()   # I am not sure of the function of .communicate()
    if not debug:
        os.remove("A_" + paths_name[0] + ".txt")
        os.remove("A_" + paths_name[1] + ".txt")
    return None

# test the code
if __name__ == "__main__":

    paths = [[1,2,3,4,5,6,7,8,9,10], [-2,-3,4,5,-4,-3,-2,-1,6,7,8,8,8,-7]]
    save_paths_in_file(paths,  paths_name="A_path", file_name="A_test_A", PATH="subpaths_SIM/")
    #dot_plot(paths[0], paths[1], paths_name=["path_A_test", "path_B_test"], Linux=True, debug=True)

    #arrow_plot_csv(file_name="Syn9R_MM_sol_vs_Chantal-arrows_5", Linux=True, debug=True)
    #arrow_plot(paths, paths_name=["path_A_test", "path_B_test"], file_name="A_TEST", Linux=True, debug=True)

"""

# Dot plot of all syn9R strains
if __name__ == "__main__":
	compared_sol = pd.read_excel("Syn9R_MM_sol_vs_Chantal-16-scores_sorted.xlsx", engine="openpyxl")
	ID_MM = list(compared_sol["ID_MM"])
	ID_Chantal = list(compared_sol["ID_Chantal"])
	sol_MM = list(compared_sol["MM_seq"])
	sol_Chantal = list(compared_sol["Chantal_seq"])
	Percentage_score = list(compared_sol["Percentage_score"])

	JS094 = list(range(1, 45))

	for i in range(len(ID_MM)):
		ID_MM_new = ID_MM[i] + "_LR"
		ID_Chantal_new = ID_Chantal[i] + "_SR"
		if Percentage_score[i] != 1:
			print(i, Percentage_score[i])
			dot_plot(sol_MM[i], sol_Chantal[i], paths_name=[ID_MM_new, ID_Chantal_new], Linux=True, debug=False)

		dot_plot(JS094, sol_MM[i], paths_name=["JS094", ID_MM_new], Linux=True, debug=False)

"""
