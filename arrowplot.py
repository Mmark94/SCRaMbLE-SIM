from comparison_sol import str_path
#from Global_alignment_MM_functions import seq_to_alignement
from matplotlib import colors
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def reverse_path(path):
    rev_path = []
    for i in path:
        rev_path.insert(0, i)
    return rev_path

def arrowplot(path, max_LU=0, essential=[], CEN=[], names="", essential_col=True, filename="Arrowplot", REF=False):
    if isinstance(path, str):
        path = str_path(path)
    if isinstance(path[0], int):  # There is only one path in the list
        path = [path]
    # This will reverse the order of the paths. It is useful to plot the paths from the top to the bottom.
    path = reverse_path(path)

    if names == "":
        names = range(len(path))
    names = reverse_path(names)

    # get discrete colormap
    if max_LU == 0:
        #max_LU = max([max(x) for x in path])
        path_join = []
        for P in path:
            path_join = path_join + P
        path_abs = [abs(x) for x in path_join]
        max_LU = max(path_abs)
    # You can choose different gradients here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    colors_non_esse = cm.Blues(np.linspace(0, 1, max_LU - len(essential)))
    if CEN != []:
        colors_esse = cm.Reds(np.linspace(0, 1, len(essential)-1))      # -1 because there is one centromere in yellow
    else:
        colors_esse = cm.Reds(np.linspace(0, 1, len(essential)))
    # Choose the colour gradient. I like gist_rainbow and jet. jet might be a little be better for color blindness. Other options are: gist_rainbow, jet, hsv, viridis, plasma, magma, gnuplot2, rainbow
    colors_rainbow = cm.jet(np.linspace(0, 1, max_LU))
    ALL_colour = np.zeros(shape=(max_LU, 4))
    counter_esse = 0
    counter_non_esse = 0

    # Add the reference in the first line
    if REF:
        path.append(range(1, max_LU + 1))
        if names != "":
            names.append("ref")

    # Plot the arrows
    for LU in range(0, max_LU):
        if LU in CEN:
            ALL_colour[LU] = [1., 1., 0., 1.]   # this is yellow
        elif LU in essential:
            ALL_colour[LU] = colors_esse[counter_esse]
            counter_esse += 1
        else:
            ALL_colour[LU] = colors_non_esse[counter_non_esse]
            counter_non_esse += 1
    #print("ALL_colour =", ALL_colour)

    # Plot the arrow plot
    longest_path = max([len(x) for x in path])
    dpi = 250
    plt.figure(figsize=(longest_path, len(path)*1.8))
    #plt.figure(figsize=(8, len(path)*0.6), dpi=dpi)     # 8 inches
    for i in range(len(path)):  # y-axis
        #plt.scatter(range(len(path[i])), [i for _ in path[i]], color="black", label=">", s=50, zorder=3)
        for j in range(len(path[i])):  # x-axis
            if essential_col:
                colour = ALL_colour[abs(path[i][j]) - 1]
            else:
                colour = colors_rainbow[abs(path[i][j]) - 1]
            if [abs(path[i][j])] == CEN:    # If the triangle is the centromere, change the colour of the edge
                edgecolor = "red"
            else:
                edgecolor = "dimgrey"

            if path[i][j] >= 0:
                #plt.scatter(j, i, color=colour, marker=">", s=100, zorder=-1, edgecolor="dimgrey", label=abs(path[i][j]))
                # Use this if you want to center the triangles in each point (e.g. center at x=1)
                #triangle = plt.Polygon([[j-0.5,i+0.25], [j-0.5,i-0.25], [j+0.5,i]], facecolor=colour, edgecolor="dimgrey")
                # Use this if you want to center the triangles between two points (e.g. center at x=0.5)
                triangle = plt.Polygon([[j, i+0.25], [j, i-0.25], [j+1, i]], facecolor=colour, edgecolor=edgecolor, label=abs(path[i][j]))
                plt.gca().add_patch(triangle)
            else:
                #plt.scatter(j, i, color=colour, marker="<", s=100, zorder=-1, edgecolor="dimgrey", label=abs(path[i][j]))
                #triangle = plt.Polygon([[j+0.5,i+0.25], [j+0.5,i-0.25], [j-0.5,i]], facecolor=colour, edgecolor="dimgrey")
                triangle = plt.Polygon([[j+1, i+0.25], [j+1, i-0.25], [j, i]], facecolor=colour, edgecolor=edgecolor, label=abs(path[i][j]))
                plt.gca().add_patch(triangle)

    # Font size
    SMALL_SIZE = 35     #12
    MEDIUM_SIZE = 45    #15
    BIGGER_SIZE = 50    #18
    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)
    plt.rc('axes', labelsize=MEDIUM_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=BIGGER_SIZE)
    # Axis
    plt.yticks(ticks=range(len(path)), labels=names)
    plt.ylim([-0.5, len(path)-0.5])
    plt.xlim([-0.5, longest_path+0.5])
    plt.ylabel("LoxPsym Units (LUs)")
    # This will make sure that the labels in the legend don't get duplicated
    #handles, labels = plt.gca().get_legend_handles_labels()
    #by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    if essential_col:
        esse = "_essential"
    else:
        esse = ""
    plt.savefig("arrow_plots\\" + filename + esse + "_arrowplot.png", format="png", dpi=dpi, bbox_inches='tight')
    plt.savefig("arrow_plots\\" + filename + esse + "_arrowplot.svg", format="svg", dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.close()
    return None



if __name__ == '__main__':
    # These are for global alignment
    match_score = 2
    mismatch_score = -1
    gap_penalty = -2  # gap penalty (gamma)

    A = [1,2,3,9,4,5]
    B = [2,3,5,-1,7,4,5]
    arrowplot([A,B], essential_col=False, REF=False)
    #print(seq_to_alignement(A,B, match_score, mismatch_score, gap_penalty))
    JS731_SR = [1, 2, 3, 4, 5, 6, 7, -13, -12, -1, -44, -43, -42, -41, -40, 39, -38, -37, -36, 8, 9, 11, 14, 15, 16, 17, 19, 20, 21, -37, -36, 19, 20, 14, 15, 16, -24, -23, -22, 17, 18, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -17, -16, -15, -14, -11, -10, -9, -8, 36, 37, 38, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, 38, -39, 40, 41, 42, 43, 44]
    JS731_LR = [1, 2, 3, 4, 5, 6, 7, -13, -12, -1, -44, -43, -42, -41, -40, 39, -38, -37, -36, 8, 9, 11, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, -21, -20, -19, -17, -38, -37, -36, 8, 9, 10, 11, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, 38, -39, 40, 41, 42, 43, 44]

    # Syn9R SCRaMbLE
    synIXR = list(range(1, 45, 1))
    essential = [2, 7, 9, 10, 12, 20]
    CEN = [2]
    names = ["synIXR", "JS731_SR", "JS731_LR"]
    arrowplot([synIXR, JS731_SR, JS731_LR], max_LU=44, essential=essential, CEN=CEN, names=names, essential_col=True, filename="JS731", REF=False)
    #arrowplot([synIXR, JS731_SR, JS731_LR], max_LU=44, essential=essential, CEN=CEN, essential_col=False, filename="A_test")

    # Syn3 SCRaMbLE
    syn3 = list(range(1, 100, 1))
    essential_3 = [7, 9, 13, 20, 28, 36, 37, 48, 70, 73, 74, 84, 97]
    CEN3 = [37]
