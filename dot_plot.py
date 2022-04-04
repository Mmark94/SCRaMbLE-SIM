#from Global_alignment_MM_functions import make_matrix
from comparison_sol import str_path
import matplotlib.pyplot as plt
from matplotlib import colors
import argparse

# This function is from Global_alignment_MM_functions
def make_matrix(x, y):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*len(y) for _ in range(len(x))]

def reverse_path(path):
    rev_path = []
    for i in path:
        rev_path.insert(0, i)
    return rev_path

# Path1 is on the y-axis, path2 is on the x-axis
def dot_plot_two_paths(path1, path2, path1_name="A", path2_name="B", pixel=False):
    # If the input is a string, covert is in a list
    if isinstance(path1, str):
        path1 = str_path(path1)
    if isinstance(path2, str):
        path2 = str_path(path2)
    #path1 = reverse_path(path1)
    path1_abs = [abs(x) for x in path1]
    path2_abs = [abs(x) for x in path2]
    # Store the coordinates of the points
    points_x, points_y = [], []
    points_inv_x, points_inv_y = [], []

    unmapped_x, unmapped_y = [], []

    Matrix = make_matrix(path1, path2)
    for i in range(len(path1)):     # This is on the y-axis
        for j in range(len(path2)):     # This is on the x-axis
            if path1[i] == path2[j]:
                Matrix[i][j] = 1
                points_x.append(j)
                points_y.append(i)
            elif abs(path1[i]) == abs(path2[j]):
                Matrix[i][j] = 2
                points_inv_x.append(j)
                points_inv_y.append(i)
            # Check if there are LUs that are not present in the other path
            if abs(path2[j]) not in path1_abs:
                unmapped_x.append(j)
        if abs(path1[i]) not in path2_abs:
            unmapped_y.append(i)
    #for row in Matrix:
    #    print(row)

    # Plot the data on a grid
    # Here you can chose the colours to use in the image. Here there is a list of colours you can use: https://matplotlib.org/stable/gallery/color/named_colors.html
    colours_grid = ["white", "black", "blue"]
    Cmap = colors.ListedColormap(colours_grid)

    # Show the image
    # I want the dimension of the plot to be around 14
    scale = 14 / max(len(path1), len(path2))
    dpi = 120
    plt.figure(figsize=(scale*len(path2), scale*len(path1)), dpi=dpi)
    if pixel:
        plt.imshow(Matrix, cmap=Cmap, origin="lower")
    else:
        plt.scatter(points_x, points_y, color="black", label="+", s=50, zorder=3)
        plt.scatter(points_inv_x, points_inv_y, color="tab:blue", label="-", s=50, zorder=3)
        if unmapped_y != []:
            plt.scatter([0 for _ in range(len(unmapped_y))], unmapped_y, color="tab:red", label="deleted", s=100, marker="P", zorder=3)
        if unmapped_x != []:
            plt.scatter(unmapped_x, [0 for _ in range(len(unmapped_x))], color="tab:red", label="deleted", s=100, marker="P", zorder=3)
        # This will make sure that the labels in the legend don't get duplicated
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
    # plt.axis('off')
    #plt.rc('axes', axisbelow=True)
    #plt.rcParams['axes.axisbelow'] = True
    plt.xticks(range(len(path2)), path2)
    plt.yticks(range(len(path1)), path1)
    plt.grid(color="grey", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)   # , alpha=0.5, zorder=-1.0
    plt.xticks(rotation=90)
    plt.xlabel(path2_name)
    plt.ylabel(path1_name)
    plt.title(path2_name + " vs " + path1_name)
    filename = path1_name + "_vs_" + path2_name
    plt.savefig("arrow_plots\\" + filename + "_dot_plot.png", format="png", dpi=dpi)
    plt.savefig("arrow_plots\\" + filename + "_dot_plot.svg", format="svg", dpi=dpi)
    plt.show()
    plt.close()
    return None


if __name__ == '__main__':

    A = [1,2,3,9,4,5]
    B = [2,3,5,-1,7,4, 5]
    #dot_plot_two_paths(A, B, path1_name="A", path2_name="B", pixel=False)
    #dot_plot_two_paths(A, B, path1_name="A", path2_name="B", pixel=True)

    JS731_SR = [1, 2, 3, 4, 5, 6, 7, -13, -12, -1, -44, -43, -42, -41, -40, 39, -38, -37, -36, 8, 9, 11, 14, 15, 16, 17, 19, 20, 21, -37, -36, 19, 20, 14, 15, 16, -24, -23, -22, 17, 18, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -17, -16, -15, -14, -11, -10, -9, -8, 36, 37, 38, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, 38, -39, 40, 41, 42, 43, 44]
    JS731_LR = [1, 2, 3, 4, 5, 6, 7, -13, -12, -1, -44, -43, -42, -41, -40, 39, -38, -37, -36, 8, 9, 11, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, -21, -20, -19, -17, -38, -37, -36, 8, 9, 10, 11, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, -18, -17, 22, 23, 24, -16, -15, -14, -20, -19, 36, 37, 38, -39, 40, 41, 42, 43, 44]
    #dot_plot_two_paths(JS731_SR, JS731_LR, path1_name="JS731_SR", path2_name="JS731_LR", pixel=False)
    S6887 = [1, -4, -3, -2, 5, 6, 3, 2, -3, 4, 5, 6, 3, 4, 5, 7, 11, -13, -12, -8, -10, 9, 10, 11, 9, -14, -13, 14, -8, 12, 15]
    #dot_plot_two_paths(S6887, S6887, path1_name="S6887", path2_name="S6887")

    parser = argparse.ArgumentParser(description="Generate a dot plot between two paths or chromosomes")
    parser.add_argument("-path1", "--path1", type=str, metavar="", required=True, help="The first path (y-axis)")
    parser.add_argument("-path2", "--path2", type=str, metavar="", required=True, help="The second path (x-axis)")
    parser.add_argument("-name1", "--path1_name", type=str, metavar="", required=False, default="A", help="The name of the first path")
    parser.add_argument("-name2", "--path2_name", type=str, metavar="", required=False, default="B", help="The name of the second path")
    args = parser.parse_args()

    dot_plot_two_paths(path1=args.path1, path2=args.path2, path1_name=args.path1_name, path2_name=args.path2_name, pixel=False)
