import statistics
import numpy as np

def start_one(items: list, start: int):
    if start in items:
        return items[items.index(start):]+items[:items.index(start)]
    if -start in items:
        return items[items.index(-start):]+items[:items.index(-start)]
    if start not in items and -start not in items:
        if 2 not in items and -2 not in items:
            return items
        if 2 in items:
            return items[items.index(2):] + items[:items.index(2)]
        if -2 in items:
            return items[items.index(-2):] + items[:items.index(-2)]

#A=[-29,-41,-40,42,43,44,1,-4,-3,5,6,-11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,29,30,31,36,37,38,39,40,41,42,43,44,1,2,3,4,5]
#print(start_one(A,1))

# This function will change the start of a path. The new path will start from the element path[start].
def start_index(path: list, start: int):
    if start < len(path):
        return path[start:] + path[:start]
    else:
        print("The start is not in the path")
        print("path =", path, "start = ", start)
        return path
#A = [1,2,3]
#print(start_index(A, 0))
#print(start_index(A, 1))
#print(start_index(A, 2))
#print(start_index(A, 3))


def LoxP_unit_count(Path, unit):
    return Path.count(unit)+Path.count(-unit)

def LoxP_unit_count_list(Path, list_unit):
    LP_unit_count = 0
    for unit in list_unit:
        LP_unit_count = LP_unit_count + Path.count(unit) + Path.count(-unit)
    return LP_unit_count

def LoxP_unit_count_Dict(Path):
    Dic = {}
    for unit in Path:
        if abs(unit) not in Dic:
            Dic[abs(unit)] = Path.count(unit)+Path.count(-unit)
    return Dic

def LoxP_unit_count_Dict_list(Path, ORI):
    Dic = {}
    for unit in ORI:
        if abs(unit) not in Dic:
            Dic[abs(unit)] = Path.count(unit)+Path.count(-unit)
    return Dic

# Count how many different LU are in one chromosome
def count_different_LU(Chr):
    new_chr = [abs(x) for x in Chr]
    new_chr = sorted(new_chr)
    new_chr = sorted(new_chr, key=new_chr.count, reverse=True)
    # The dictionary keys cannot be repeated
    unique = list(dict.fromkeys(new_chr))
    return len(unique)

# Check that the LU number is decrescent in a chromosome. E.g. the number of 1 is less than 2, 2 < 3, ...
# This script is used to create minimal chromosomes. E.g. [1,1,2,3] is the same as [2,2,1,3]
def descescent_LU_num(Chr):
    # Count how many different LU are in one chromosome
    different_LU = count_different_LU(Chr)
    count_previous_LU = 10000
    for i in range(1, different_LU + 1):
        #print(i)
        count_previous_LU_temp = LoxP_unit_count(Chr, i)
        #print("count_previous_LU =", count_previous_LU, count_previous_LU_temp)
        # The previous LU number has to be bigger or equal than the current LU number
        if count_previous_LU < count_previous_LU_temp:
            return False
        count_previous_LU = count_previous_LU_temp
    return True
"""
A = [1,-1,2,3,3,-5]
B = [1,-1,2,3,-5]
C = [1,-1,1,2,3,3,4,-5]
D = [1,-1,1,2,-2,3,3,4,-5]
print(descescent_LU_num(A))
print(descescent_LU_num(B))
print(descescent_LU_num(C))
print(descescent_LU_num(D))
"""


def invert(x):
    inv_x=[]
    for element in x:
        if isinstance(element, int):
            inv_x.insert(0,-element)
        else:
            inv_x.insert(0,element)
    return inv_x

# This is not useful and not tested
def count_pos_neg(solution: list):
    pos = 0
    neg = 0
    for element in solution:
        if isinstance(element, int):
            if element >= 0:
                pos = pos + 1
            if element < 0:
                neg = neg + 1
    return pos, neg

def half_pos(x):
    inv_x=[]
    neg = 0
    pos = 0
    for element in x:
        if isinstance(element, int):
            if (element >= 0):
                pos = pos + 1
            if(element < 0):
                neg = neg + 1
    if neg > pos:
        for element in x:
            if isinstance(element, int):
                inv_x.insert(0, -element)
            else:
                inv_x.insert(0, element)
        return inv_x
    else:
        return x

def half_pos2(path):
    if path == []:
        return path
    if isinstance(path[0], list):
        new_paths=[]
        for p in path:
            new_paths.append(half_pos2(p))
        return new_paths
    else:
        neg = 0
        pos = 0
        for element in path:
            if isinstance(element, int):
                if element >= 0:
                    pos += 1
                elif element < 0:
                    neg += 1
        if neg > pos:
            return invert(path)
        elif neg == pos:
            # if there are the same number of LU positive and negative, the program will decide with the first LU, usually 1 or -1
            direction_start = path[0]
            direction_invert_start = invert(path)[0]
            if direction_start > 0 and direction_invert_start > 0:
                if direction_start <= direction_invert_start:
                    return path
                else:
                    return invert(path)
            elif direction_start > 0 and direction_invert_start < 0:
                return path
            elif direction_start < 0 and direction_invert_start > 0:
                return invert(path)
            elif direction_start < 0 and direction_invert_start < 0:
                if abs(direction_start) <= abs(direction_invert_start):
                    return path
                else:
                    return invert(path)
        else:
            return path


def half_pos_one(solution: list, circular=False):
    if solution == []:
        return solution
    if isinstance(solution[0], list):
        new_solution_list=[]
        for sol in solution:
            new_solution_list.append(half_pos_one(sol, circular=circular))
        return new_solution_list
    else:
        new_solution = solution[:]
        # make it start with one and then find the orientation and if it invert it, it makes it start with one again
        if circular:
            new_solution = start_one(new_solution, 1)
            new_solution = half_pos2(new_solution)
            new_solution = start_one(new_solution, 1)
        else:
            new_solution = half_pos2(new_solution)
        return new_solution

"""
solution_ov   = [-1,2,3,-4,-5,6]
SCRaMbLEd_chr = [-6,5,4,-3,-2,1]
#print(half_pos_one(solution_ov))
#print(half_pos_one(SCRaMbLEd_chr))
solution_ov   = [13, 14, -21, -20, -17, 18, -22, -21, -20, 17, -12, -11, -10, 7, 8, 9, -6, 2, 3, 4, 5, -1]
SCRaMbLEd_chr = [1, -5, -4, -3, -2, 6, -9, -8, -7, 10, 11, 12, -17, 20, 21, 22, -18, 17, 20, 21, -14, -13]
A=[1,2,3,4,5,6,-7,-8,-9,-10]
B=[-1,-2,-3,-4,-5,-6,7,8,9,10]
C=[3,4,5,6,7,1,3,4]
D=[-10,-9,-8,-4,1,-2]
E=[[1,2,3,4,5,6,-7,-8,-9,-10],
[-1,-2,-3,-4,-5,-6,7,8,9,10],
[3,4,5,6,7,1,3,4],
[-10,-9,-8,-4,1,-2]]
print(half_pos_one(A))
print(half_pos_one(B))
print(half_pos_one(C))
print(half_pos_one(D))
print(half_pos_one(E))
"""


# Convert the LUs string into integers. It also don't modify "*", the sign when the subpath calling program could not call the LU.
def list_int(sol):
    if sol == [""]:
        return []
    sol2 = []
    for i in sol:
        if i == "*" or i == "'*'":
            sol2.append("*")
            #sol2.append(0)     #if you want to substitute the "*" with zero (0)
        else:
            sol2.append(int(i))
    return sol2

#A=['-44', '-43', '-42', '-41', '-40', '-39', '7', '8', '9', '10', '11', '12', '16', '17', '18', '19', '20', '21', '22', '23', '-24', '25', '26', '27', '28', '29']
#print(list_int(A))

def str_path0(x: str):
    to_list = x.split(",")
    # list_int converts strings in integers
    path = list_int(to_list)
    return path

def str_path(x: str):
    if x == 1 or x == "1" or x == 0 or x == "0" or x == "":
        return []
    if x == "[]":
        return []
    # Sometimes when I retrieve reads from excel it cuts the last reads (limit of 32,767 characters in each cell)
    if x[0] == "[" and x[-1] != "]":
        print("Careful the reads file is incomplete!!!")
    x = x.replace(" ", "")
    if x.count("[[") > 0:
        L_of_L = []
        start = 2   # discard the first "[["
        for i in range(len(x)):
            if i == len(x)-1:  # if it is the last ] finish the function
                return L_of_L
            if x[i] == "]":
                #print(x[start:i])
                L_of_L.append(str_path0(x[start:i]))
                #start = i + 3   # discard "],["
                start = x.find("[", i) + 1  # discard "], ["
        #return L_of_L
    elif x.count("[") + x.count("]") == 2:
        return str_path0(x[1:-1])
    elif x.count("[") + x.count("]") > 2:
        print("check your Input, the characters [ and ] are more than two")
    elif x.count("[") == 0:
        return str_path0(x)

#A = "[[-43, -42], [39, 40, 41, 42, 43, 44], [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 1, 2, 3, 4], [30], [-3, -2, -1], [-37, -36, -35, -34, '*', -33, -32], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44], [-8, -7, -6, -5, -4], [1], [-37], [-3, -2, -1], [-39, -38], [7, 8, 9, 10, 11], [26, 27, 28, 29, 30], [-8, -7, -6], [20], [-30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14], [-37, -36, -35, -34, -33], [44, 1], [-30], [-24, -23, 23], [29], [33, 34, 35, 36, 37, 38], [1, '*', 3, 4, 5, 6, 7], [37, 38, 39, 40, 41, 42, 43], [24, 25], [-43, -42], [38, '*', 41, 42, 43, 44, 1, '*', 2, 3, 4, 5], [7], [-17, -17], [42, 43, 44], [34, 35, 36, 37, '*', 38, 39, 40, 41, 42, 43, 44], [4], [44], [-29, -28, -27, -26, -25, -24], [30], [-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4], [38, 39, 40, 41], [-44, -43, -42], [-14, -13, -12], [38, 39, 40, 41, 42, 43, 44], [16, 17, '*', 19, 20, 21, 22, 23, 24, 25], [-3, -2, -1], [-10], [-3, -2, -1], [-20, -19, -18, -17, -16, -15, -14, -13, -12, -10, -9, -8]]"
#print(str_path(A))
#A="1,2,3,4,*,5,6,-7,-8"
#B="[1,2,3,4,4,5,6,-7,-8]"
#C="[[1,2,3,4],[4,5,6,-7,-8]]"
#D = "[[1, 2, 3], [37, 38, 44, 44], [-1, -44, 36, 37, 38], [4, 5, -6], [10, 11, 8], []]"
#print(str_path(A))
#print(str_path(B))
#print(str_path(C))
#print(str_path(D))

def path_str(x: list):
    to_list = [str(i) for i in x]
    to_string = ",".join(to_list)
    return to_string
#B=[1, 2, 3, 4, '*', 5, 6, -7, -8]
#print(path_str(B))

def flatten_level_lists(list_of_lists):
    new_list = []
    if list_of_lists == []:
        new_list.append(list_of_lists)
    elif isinstance(list_of_lists[0], list):
        for i in list_of_lists:
            new_list = new_list + flatten_level_lists(i)
    else:
        new_list.append(list_of_lists)
    return new_list

#A = [[1,2,3],[2,3,4,5],[[4,5,6],[[2,3,4,5],[5,6,6,7,8]],[5,6,7,8]],[], [8,9,10]]
#B = one_level_lists(A)
#print(A)
#print(B)

def remove_duplicate(paths:list, list_of_list=False):
    if paths == []:
        return paths
    if isinstance(paths[0], list):
        new_paths = []
        for path in paths:
            #if path == invert(path):
            #    print("Palindromic path =", path)
            if path not in new_paths and invert(path) not in new_paths:
                new_paths.append(half_pos2(path))
        if list_of_list:        # This keeps the output type constant, always a list of list: [[1,2,3]]
            return new_paths
        else:
            return new_paths[0] if len(new_paths) == 1 else new_paths
    else:
        return paths

#print(remove_duplicate([[1,2,3],[-3,-2,-1]]))

# This function counts how many solutions are present in a list.
def count_solutions(solutions):
    solutions = remove_duplicate(solutions)
    if solutions == [] or solutions == [[]]:
        print("There are zero solutions.")
        return 0
    if isinstance(solutions[0], list):
        number_sol = len(solutions)
    else:
        number_sol = 1
    return number_sol


# Count and print how many multiple solutions are present in a simulation experiment. The input could either be the solutions (in which case the function will count how many solutions) or it could be directely the number of solutions.
def count_multiple_sol(solutions_LIST=[], solution_len_LIST=[], name="", only_solved=False, solved_L=[]):
    len_sol_ORI = max(len(solutions_LIST), len(solution_len_LIST))
    if only_solved and solved_L != []:     # Only use solved solutions
        solutions_LIST_only_solved = []
        solution_len_LIST_only_solved = []
        for i in range(0, max(len(solutions_LIST), len(solution_len_LIST))):
            if solved_L[i] == 1:
                if solutions_LIST != []:
                    solutions_LIST_only_solved.append(solutions_LIST[i])
                elif solution_len_LIST != []:
                    solution_len_LIST_only_solved.append(solution_len_LIST[i])
        # Replace the list of solutions with the one with only solved solutions
        solutions_LIST = solutions_LIST_only_solved[:]
        solution_len_LIST = solution_len_LIST_only_solved[:]
        len_only_solved = max(len(solutions_LIST_only_solved), len(solution_len_LIST_only_solved))
        percentage_only_solved = (len_only_solved / len_sol_ORI) * 100
        print("length solutions_LIST_only_solved", name, "=", len_only_solved, "percentage =", percentage_only_solved)
        print()
    if solutions_LIST == [] and solution_len_LIST == []:
        print("None of the solutions is correct!")
        return None
    if solution_len_LIST == []:
        if isinstance(solutions_LIST[0], str):
            solution_len_LIST = [count_solutions(str_path(sol)) for sol in solutions_LIST]
        else:
            solution_len_LIST = [count_solutions(sol) for sol in solutions_LIST]
    count_multi_sol = 0
    for sol in solution_len_LIST:
        if sol > 1:
            count_multi_sol += 1
    percentage_multiple_solution_ALL = count_multi_sol / len_sol_ORI * 100
    percentage_multiple_solution = count_multi_sol / len(solution_len_LIST) * 100
    multiple_solution_average = statistics.mean(solution_len_LIST)
    multiple_solution_std = statistics.stdev(solution_len_LIST)
    print("number of multiple solution solved in", name, "=", count_multi_sol)
    print(name + " percentage multiple solution ALL =", percentage_multiple_solution_ALL)
    print(name + " percentage multiple solution only on the solved =", percentage_multiple_solution)
    print("solution_len_LIST =", sorted(solution_len_LIST, reverse=True)[:11])
    print(name + " multiple solution average =", multiple_solution_average, "+-", multiple_solution_std)
    print(name + " multiple solution quantiles =", np.percentile(solution_len_LIST, [25, 50, 75]))
    if only_solved:
        print(name, "algorithm =", [percentage_only_solved, percentage_only_solved - percentage_multiple_solution_ALL, percentage_multiple_solution_ALL, 100-percentage_only_solved])
    print()
    return percentage_multiple_solution, multiple_solution_average


def comparison_sol(path1, path2):
    return path1 == path2

# This function will compare the two paths starting from their first element LU.
def compare_sol(path1: list, path2: list):
    path1_L = len(path1)
    path2_L = len(path2)
    path_L_min = min(path1_L, path2_L)
    path_L_max = max(path1_L, path2_L)
    if path1_L == 0 or path2_L == 0:
        return 0, 0
    same = 0
    i = 0
    while i < path_L_min:
        if path1[i] == path2[i]:
            same += 1
        i += 1
    return same, same/path_L_max

#A = [1,2,3,4,5,6,-7,-8,-9,-10]
#B = [1,2,-3,4,12,6,-7,-8,-9,-10]
#print(compare_sol(A, B))


# This function will compare the two paths starting from the first LU that they have in common
def compare_sol_no_start(path1: list, path2: list):
    path1_L = len(path1)
    path2_L = len(path2)
    path_L_min = min(path1_L, path2_L)
    path_L_max = max(path1_L, path2_L)
    # find an element (LU) in the path2 that match the first element in path1
    i = 0
    while i < path2_L:
        if path1[0] == path2[i]:
            # start_path2 = i
            break
        i += 1
    # modify the path2 to make it starts with the first element of path1. This could be dangerous as it changes the sequence of path2. And it is assumes that path2 is circular.
    new_path2 = path2[i:] + path2[:i]
    # count how many element are the same
    same = 0
    i = 0
    while i < path_L_min:
        if path1[i] == new_path2[i]:
            same += 1
        i += 1
    return same, same/path_L_max

#A = [4,5,6,-7,-8,-9,-10,1,2,-3]
#B = [1,2,3,4,12,6,-7,-8,-9,-10]
#print(compare_sol_no_start(A, B))

#C = [4,5,6,-7,-8,-9,-10,1,2,-3,-4]
#D = [1,2,3,4,12,6,-7,-8,-9,-10]
#print(compare_sol_no_start(C, D))

# This function will compare the two paths starting from every element (LU) of the path1 and output True if they are the same
def compare_sol_every_points_true(path1: list, path2: list):
    path1_L = len(path1)
    i = 0
    while i < path1_L:
        path1_new_start = start_index(path1, i)
        if path1_new_start == path2 or invert(path1_new_start) == path2:
            return True
        i += 1
    return False

# This function will compare the two paths starting from every element (LU) of the path1
def compare_sol_every_points(path1: list, path2: list):
    path1_L = len(path1)
    i = 0
    maximum_score = 0
    while i < path1_L:
        path1_new_start = start_index(path1, i)
        score = compare_sol(path1_new_start, path2)[0]
        if score > maximum_score:
            maximum_score = score
        i += 1
    return maximum_score

#A = [4,5,6,-7,-8,-9,-10,1,2,-3]
#B = [1,2,3,4,12,6,-7,-8,-9,-10]
#C = [8,8,9,10,11,12]
#D = [1,2,3,4,5,6,7,8,9,10,11,12]
#print(compare_sol_every_points(A, B))
#print(compare_sol_every_points(C, D))
#print(compare_sol_every_points(D, C))

# This function is the same as the one before "compare_sol_every_points" but it will output the best start
def find_start_path(path1: list, path2: list):
    path1_L = len(path1)
    i = 0
    maximum_score = 0
    best_start = 0
    while i < path1_L:
        path1_new_start = start_index(path1, i)
        score = compare_sol(path1_new_start, path2)[0]
        if score > maximum_score:
            maximum_score = score
            best_start = i
        i += 1
    return best_start

# This function is the same as the one before "compare_sol_every_points" but it will output the best start
def find_start_path_inv(path1: list, path2: list):
    path1_L = len(path1)
    i = 0
    maximum_score = 0
    best_start = 0
    while i < path1_L:
        path1_new_start = start_index(path1, i)
        score = compare_sol(path1_new_start, path2)[0]
        if score > maximum_score:
            maximum_score = score
            best_start = i
        i += 1
    # Try to find the start in the inverted path
    path1 = invert(path1)
    i = 0
    while i < path1_L:
        path1_new_start = start_index(path1, i)
        score = compare_sol(path1_new_start, path2)[0]
        if score > maximum_score:
            maximum_score = score
            best_start = -i     # For the inverted sequence I put the minus sign
        i += 1
    return best_start

# This try to find the best alignment between two paths. First it applies the function half_pos_one and then it try to find the start.
def align_two_paths(path1: list, path2: list, circular=True):
    new_path1 = half_pos_one(path1, circular=circular)
    new_path2 = half_pos_one(path2, circular=circular)
    start = find_start_path(new_path1, new_path2)
    new_path1_new_start = start_index(new_path1, start)
    return new_path1_new_start, new_path2

# This try to find the best alignment between two paths. First it applies the function half_pos_one and then it try to find the start.
def align_two_paths_inv(path1: list, path2: list, circular=True):
    new_path1 = half_pos_one(path1, circular=circular)
    new_path2 = half_pos_one(path2, circular=circular)
    start = find_start_path_inv(new_path1, new_path2)
    if start >= 0:
        new_path1_new_start = start_index(new_path1, start)
    else:
        new_path1 = invert(new_path1)       # Invert path1
        new_path1_new_start = start_index(new_path1, -start)        # Remember to change the sign of start
    return new_path1_new_start, new_path2

#A = [1,2,3,-4,-5,-6,-7,8,9,10]
#B = [-8,7,6,5,4,-3]
#Align = align_two_paths(A, B)
#Align_inv = align_two_paths_inv(A, B)
#print(Align)
#print(Align_inv)

# Find the best alignment and alignment score between a reference solution and a list of solutions.
def find_best_match_solutions(reference, solutions, circular=True):
    if solutions == [] or reference == []:
        return None
    if isinstance(solutions[0], int):
        solutions = [solutions]
    best_score = 0
    new_solutions = []
    counter = 1
    for sol in solutions:
        if reference == sol or reference == invert(sol):
            print("the solution number", counter, "is the reference!")
        new_ref, new_sol = align_two_paths_inv(reference, sol, circular=circular)
        new_solutions.append(new_sol)
        print("counter =", counter)
        #print("Ref =", new_ref)
        print("new_sol =", new_sol)
        score_MM, score_MM_L = compare_both(new_ref, new_sol)
        print(score_MM, score_MM_L)
        counter += 1
        if score_MM_L > best_score:
            best_score = score_MM_L
            best_sol = new_sol
    print("best_score =", best_score)
    print("Ref      =", reference)
    print("best_sol =", best_sol)
    return new_solutions

# Use the function compare_sol_every_points with the longest path as first.
def compare_longest(path1: list, path2: list):
    if len(path1) >= len(path2):
        return compare_sol_every_points(path1, path2)
    else:
        return compare_sol_every_points(path2, path1)

# Use the function compare_sol_every_points in both direction.
def compare_both(path1: list, path2: list):
    path_L_max = max(len(path1), len(path2))
    score1_2 = compare_sol_every_points(path1, path2)
    score2_1 = compare_sol_every_points(path2, path1)
    score_max = max(score1_2, score2_1)
    return score_max, score_max/path_L_max

def combine_score_len(path1: list, path2: list):
    length = max(len(path1), len(path2))
    score = compare_both(path1, path2)
    return score/length

#A = [1, -4, -3, 5, 6, -11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, -27, -26, -25, -24, 28, 3, -2, 4, 5, 6, 7, 8, 9, 10, -6, -5, -4, -14]
#B = [-2, 4, 5, 6, -10, -9, -8, -7, -6, -5, -4, -15, -14, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, -27, -26, -25, -24, 28, 3]
#C = [-2, 4, 5, 6, 7, 8, 9, 10, -6, -5, -4, -15, -14, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, -27, -26, -25, -24, 28, 3]
#print(compare_sol_every_points(A, B))
#print(compare_sol_every_points(A, C))

def check_sol_in_list_lin(sol, list_sol):
    if list_sol == []:
        return False
    if isinstance(list_sol[0], list):
        for i in list_sol:
            if i == sol or invert(i) == sol:
                return True
        # All the solution in the list do not match the real solution
        return False
    else:
        return sol == list_sol or sol == invert(list_sol)


# This compare two solutions or one solution to many. This is for circular chromosomes as it compare every point in a solution.
def check_sol_in_list_cir(sol, list_sol):
    if list_sol == []:
        return False
    if isinstance(list_sol[0], list):
        for i in list_sol:
            if compare_sol_every_points_true(i, sol):       # compare_sol_every_points_true also compare inverted sequences
                return True
        # All the solution in the list do not match the real solution
        return False
    else:
        return compare_sol_every_points_true(list_sol, sol)


# This compare two solutions or one solution to many
def check_sol_in_list_lin_circular(sol, list_sol, circular=False):
    if circular:
        return check_sol_in_list_cir(sol, list_sol)
    else:
        return check_sol_in_list_lin(sol, list_sol)

#A = [1,2,3]
#B = [-3,-2,-1]
#B = [[-3,-2,-1], []]
#C = [[-2,-1,-3], []]
#compared0 = check_sol_in_list_lin_circular(A, C, circular=False)
#compared1 = check_sol_in_list_lin_circular(A, C, circular=True)
#print("compared0 =", compared0)
#print("compared1 =", compared1)

# Calculate the difference in length between two lists.
def diff_len(path1: list, path2: list):
    diff_length = abs(len(path1) - len(path2))
    return diff_length



# Mapping

def mapping(solution, path):
    if solution == [] or path == []:
        return False
    L_path = len(path)
    c=0
    while c + L_path < len(solution) + 1:
        k_mer = solution[c:L_path+c]
        if path == k_mer or invert(path) == k_mer:
            return [c, L_path+c-1, True]
        c = c + 1
    return False

def mapping_every_points_true(solution: list, path: list):
    if solution == [] or path == []:
        return False
    T=0
    while T < len(solution):
        solution_new_start = start_index(solution, T)
        if mapping(solution_new_start, path):
            return True
        T = T + 1
    return False

# This function returns True if all the subpaths are compatible with the solution
def paths_in_sol(solution, paths, debug=False):
    if solution == [] or paths == [] or paths == [[]]:
        return False
    for path in paths:
        if path == []:
            continue
        if not(mapping_every_points_true(solution, path)) or not(mapping_every_points_true(solution, invert(path))):
            if debug:
                print("solution =", solution)
                print("read =", path)
            return False
    return True



# This is very slow
# By default it treats the solutions as circular.
def check_sol(solutions, paths, debug=False):
    if solutions == []:
        return solutions
    solutions = remove_duplicate(solutions)
    if isinstance(solutions[0], list):
        new_solutions = []
        for sol in solutions:
            if paths_in_sol(sol, paths, debug=debug):
                new_solutions.append(sol)
        return new_solutions
        #return new_solutions[0] if len(new_solutions)==1 else new_solutions
    else:
        if paths_in_sol(solutions, paths, debug=debug):
            return solutions
        else:
            return []

# By default it treats the solutions as circular.
# This function returns all the subpaths that are compatible with the solution
def extract_paths_from_sol(solution, paths):
    if solution == [] or solution == [[]] or paths == [] or paths == [[]] or isinstance(solution, int):
        return [paths, []]
    paths_compatible = []
    if isinstance(solution[0], list):
        solutions = remove_duplicate(solution)
        for sol in solutions:
            paths_compatible = paths_compatible + extract_paths_from_sol(sol, paths)[0]
    else:
        for path in paths:
            if path == []:
                continue
            # Check if the read can map in the circular solution
            if mapping_every_points_true(solution, path) or mapping_every_points_true(solution, invert(path)):
                paths_compatible.append(path)
    # Remove duplicates
    paths_compatible = remove_duplicate(paths_compatible, list_of_list=True)
    #print("paths_compatible =", paths_compatible)
    paths_not_compatible = []
    for path in paths:
        if path not in paths_compatible and invert(path) not in paths_compatible:
            paths_not_compatible.append(path)

    return [paths_compatible, paths_not_compatible]


# By default it treats the solutions as circular.
# This calculate the number of paths that do not fit in the solution (extract_paths_from_sol) and there are less non compatible paths than the threshold, it returns the solution. Otherwise, it excludes the solution.
def check_sol_threshold(solutions, paths, threshold=0):
    if solutions == [] or solutions == [[]] or paths == [] or paths == [[]] or isinstance(solutions, int):
        return solutions
    solutions = remove_duplicate(solutions)
    if isinstance(solutions[0], list):
        new_solutions = []
        for sol in solutions:
            paths_compatible, paths_not_compatible = extract_paths_from_sol(solution=sol, paths=paths)
            if len(paths_not_compatible) <= threshold:
                new_solutions.append(sol)
        return new_solutions
        #return new_solutions[0] if len(new_solutions)==1 else new_solutions
    else:
        paths_compatible, paths_not_compatible = extract_paths_from_sol(solution=solutions, paths=paths)
        if len(paths_not_compatible) <= threshold:
            return solutions
        else:
            return []

#sol = [[-2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37, 38, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, -26, -25, -24, -23, -22], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37, 38, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, -26, -25, -24, -23, -22], [-2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, -22, -21, -20], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, -22, -21, -20]]
#paths = [[7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13], [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28], [20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37, 38, 5, 6, 7, 8, 9, 10, 12, 13, 14], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 36], [13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 36, 37], [25, 26, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -10, -9, -8, -7, -6, -5, -38], [-2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [-3, -2, 3, 4, 5]]

#sol = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
#paths = [[1,2,3,5,6,7,8], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44], [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30], [8, 9, 10, 11, 12, 13, 14, 15, 16], [41, 42, 43, 44, 1], [7, 8, 9, 10, 11]]

#print(extract_paths_from_sol(sol, paths))
#print(sol)
#print(check_sol_threshold(sol, paths, threshold=0))
#print(check_sol_threshold(sol, paths, threshold=1))

#JS599 = [[1, 2, 3, 4, -7, 5, 35, 38, 39, 41, 43, 43, 44], [1, 2, 3, 4, -7, 5, 37, 38, 39, 41, 42, 44, 44], [1, 2, 3, 4, -7, 5, 37, 38, 39, 41, 42, 43, 44], [-33, -32, 29, 30, -22, -20, -19, -7, 8, 9, 10, 11, 12], [1, -12, -11, -10, -9, -8, 7, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44], [-12, -11, -10, -9, -8, 7, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44, 44], [2, 3, 4, 5, 36, 37, 38, 39, 41, 42, 43, 44, 44], [36, 37, 37, 39, 41, 42, 43, 44], [11, 8, 9, 10, 12, 19, 20, -31, -30, -29, 32, 33, 36], [1, 4, 5, -6, 44]]
#JS599_paths = [[27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44], [22, 23, 24, 25, 26, 27, 29, 30, 31, 32], [19, 20, 22, 23, 24, 25, 26, 27, 29, 30], [-11, -10, -9, -8, 7, 19, 20, 22, 23, 24], [9, 10, 12, 19, 20, -31, -30, -29], [20, -31, -30, -29, 32, 33, 36], [11, 8, 9, 10, 12, 19], [1, 2, 3, 4, -7, 5], [44, 1, 4, 5, -6], [25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43], [35, 38, 39, 41, 43, 43, 44, 1, 2, 3, 4], [30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44, 44], [24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35], [36, 37, 38, 39, 41, 42, 43, 44, 44], [37, 38, 39, 41, 42, 43, 44, 44, 2, 3, 4], [37, 38, 39, 41, 42, 44, 44, 1, 2, 3, 4], [37, 38, 39, 41, 42, 43, 44, 1, 2, 3, 4], [2, 3, 4, 5], [-12, -11, -10, -9, -8, 7, 19, 20, 22], [7, 19, 20, 22, -30, -29, 32, 33], [34, 35, 38, 39, 41, 43, 44, 1], [36, 37, 37, 39, 41, 42, 43], [37, 39, 41, 42, 43, 44]]
#JS599_clean = check_sol(JS599, JS599_paths)
#print(JS599)
#print(JS599_clean)

#JS599 = [[1, 2, 3, 4, -7, 5, 35, 38, 39, 41, 43, 43, 44], [1, 2, 3, 4, -7, 5, 37, 38, 39, 41, 42, 44, 44], [1, 2, 3, 4, -7, 5, 37, 38, 39, 41, 42, 43, 44], [-33, -32, 29, 30, -22, -20, -19, -7, 8, 9, 10, 11, 12], [1, -12, -11, -10, -9, -8, 7, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44], [-12, -11, -10, -9, -8, 7, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44, 44], [2, 3, 4, 5, 36, 37, 38, 39, 41, 42, 43, 44, 44], [36, 37, 37, 39, 41, 42, 43, 44], [11, 8, 9, 10, 12, 19, 20, -31, -30, -29, 32, 33, 36], [1, 4, 5, -6, 44]]
#JS599_paths = [[27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44], [22, 23, 24, 25, 26, 27, 29, 30, 31, 32], [19, 20, 22, 23, 24, 25, 26, 27, 29, 30], [-11, -10, -9, -8, 7, 19, 20, 22, 23, 24], [9, 10, 12, 19, 20, -31, -30, -29], [20, -31, -30, -29, 32, 33, 36], [11, 8, 9, 10, 12, 19], [1, 2, 3, 4, -7, 5], [44, 1, 4, 5, -6], [25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 38, 39, 41, 43], [35, 38, 39, 41, 43, 43, 44, 1, 2, 3, 4], [30, 31, 32, 33, 34, 35, 38, 39, 41, 43, 44, 44], [24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35], [36, 37, 38, 39, 41, 42, 43, 44, 44], [37, 38, 39, 41, 42, 43, 44, 44, 2, 3, 4], [37, 38, 39, 41, 42, 44, 44, 1, 2, 3, 4], [37, 38, 39, 41, 42, 43, 44, 1, 2, 3, 4], [2, 3, 4, 5], [-12, -11, -10, -9, -8, 7, 19, 20, 22], [7, 19, 20, 22, -30, -29, 32, 33], [34, 35, 38, 39, 41, 43, 44, 1], [36, 37, 37, 39, 41, 42, 43], [37, 39, 41, 42, 43, 44]]
#JS599_clean = extract_paths_from_sol(JS599[0], JS599_paths)
#print(JS599)
#print(JS599_clean[0])
#print(JS599_clean[1])

#JS603 = [[-1, -44, -43, -42, -41, -40, -39, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 20, 22, 23, -24, 25, 26, 16, 17, 18, 19, 20, 22, 23, -24, 25, 26, 27, -36, -37, 28, 29, 30, 31, 32, 33, -4, -3, -2, 28, 29, 30, 31, -41, -40, 7, 8, 9, 10, 11, -12, 16, 17, 18, 19, 20, 21, 22, 23, -24, 25, 26, 27, 28, 29, 30, 36, 37, 38, -6, -5, -4, -3, -2], [-1, -44, -43, -42, -41, -40, -39, 7, 8, 9, 10, 11, 12, 16, -44, -44, -43, -42, -41, -40, -39, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 20, 22, 23, -24, 25, 26, 27, -36, -37, 28, 29, 30, 31, 32, 33, -4, -3, -2, 28, 29, 30, 31, -41, -40, 7, 8, 9, 10, 11, -12, 16, 17, 18, 19, 20, 21, 22, 23, -24, 25, 26, 27, 28, 29, 30, 36, 37, 38, -6, -5, -4, -3, -2]]
#JS603_paths = [[17, 18, 19, 20, 21, 22, 23, -24, 25, 26, 27, 28, 29, 30, 36, 37, 38, -6, -5, -4, -3], [18, 19, 20, 22, 23, -24, 25, 26, 27, -36, -37, 28, 29, 30, 31, 32, 33, -4], [17, 18, 19, 20, 22, 23, -24, 25, 26, 27, -36, -37, 28, 29, 30], [16, 17, 18, 19, 20, 22, 23, -24, 25, 26, 27, -36, -37, 28], [32, 33, -4, -3, -2, 28, 29, 30, 31, -41, -40, 7], [37, 38, -6, -5, -4, -3, -2, -1, -44], [-44, -44, -43, -42, -41, -40, -39, 7], [-41, -40, 7, 8, 9, 10, 11], [-44, -43, -42, -41, -40, -39, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 20, 22, 23, -24, 25, 26], [20, 22, 23, -24, 25, 26, 27, -36, -37, 28, 29, 30, 31, 32, 33, -4, -3], [7, 8, 9, 10, 11, -12, 16, 17, 18, 19, 20, 21, 22, 23, -24, 25, 26, 27, 28, 29], [-4, -3, -2, -1, -44, -43, -42, -41, -40, -39, 7, 8, 9, 10, 11, 12, 16], [30, 31, 32, 33, -4, -3, -2, 28, 29]]
#JS603_cleaned = check_sol(JS603, JS603_paths)
#print(JS603)
#print(JS603_cleaned)
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


def NG50_calculator(Chr):
    if Chr == []:
        return 0, 0,  0, 0, 0, 0
    new_chr = [abs(x) for x in Chr]
    new_chr = sorted(new_chr)
    new_chr = sorted(new_chr, key=new_chr.count, reverse=True)
    unique = list(dict.fromkeys(new_chr))
    #print("new_chr =", new_chr)
    #print("unique =", unique)
    N50 = round(len(new_chr)/2)
    counter = 0
    cumulative = 0
    while cumulative < N50:
        cumulative = cumulative + new_chr.count(unique[counter])
        #print("cumulative =", cumulative)
        counter = counter + 1
    NG50 = new_chr.count(unique[counter - 1])
    LG50 = counter - 1
    # Easy duplication rate
    #duplication_rate = 1 - (len(unique) / len(new_chr))
    duplication_rate = len(new_chr) / len(unique)
    return len(new_chr), len(unique),  NG50, LG50, LG50 / len(unique), duplication_rate

#Chr = [1,2,2,3,4,6,13,1,5,9,16,4,7,4,4,5,6,-3,7,8,8,2,8,-11,8,9,10,11,11,12]
#Chr_a = NG50_calculator(Chr)
#print(Chr_a)

def essential_ratio_calculator(chromosome: int, essential=[]):
    if chromosome == []:
        return 0
    essential = [abs(x) for x in essential]
    chr_L = len(chromosome)
    num_essential = 0
    # Count the number of essential LUs
    for LU in chromosome:
        if abs(LU) in essential:
            num_essential += 1
    return num_essential / chr_L

# This function calculate the base pairs length of path/solution or subpaths.
def solution_to_bps_size(solution: list):
    if solution == [] or solution == [[]]:
        return 0
    syn9R_LU_size = {1: 6375, 2: 688, 3: 1305, 4: 5549, 5: 1628, 6: 849, 7: 5111, 8: 2546, 9: 1895, 10: 5185, 11: 1835, 12: 1644, 13: 157, 14: 1830, 15: 244, 16: 4041, 17: 4716, 18: 944, 19: 3259, 20: 4337, 21: 178, 22: 1964, 23: 134, 24: 1420, 25: 4155, 26: 1496, 27: 228, 28: 1008, 29: 2048, 30: 4475, 31: 176, 32: 1444, 33: 1081, 34: 1585, 35: 216, 36: 2428, 37: 4101, 38: 2924, 39: 954, 40: 997, 41: 1803, 42: 984, 43: 1005, 44: 9386}

    solution_Kb = 0
    if isinstance(solution[0], int):
        for LU in solution:
            solution_Kb = solution_Kb + syn9R_LU_size[abs(LU)]
        return [solution_Kb]
    else:
        solution_bps_list = []
        for sol in solution:
            if sol == []:
                continue
            sol_L = solution_to_bps_size(sol)
            solution_bps_list = solution_bps_list + sol_L
        return solution_bps_list

# This function output the start and end position of each LU in a path/solution or subpath.
def position_LU(solution: list):
    if solution == [] or solution == [[]]:
        return 0
    syn9R_LU_size = {1: 6375, 2: 688, 3: 1305, 4: 5549, 5: 1628, 6: 849, 7: 5111, 8: 2546, 9: 1895, 10: 5185, 11: 1835, 12: 1644, 13: 157, 14: 1830, 15: 244, 16: 4041, 17: 4716, 18: 944, 19: 3259, 20: 4337, 21: 178, 22: 1964, 23: 134, 24: 1420, 25: 4155, 26: 1496, 27: 228, 28: 1008, 29: 2048, 30: 4475, 31: 176, 32: 1444, 33: 1081, 34: 1585, 35: 216, 36: 2428, 37: 4101, 38: 2924, 39: 954, 40: 997, 41: 1803, 42: 984, 43: 1005, 44: 9386}

    chr_size = 0
    LU_pos_dic = {}
    if isinstance(solution[0], int):
        for LU in range(len(solution)):
            LU_pos = [chr_size]         # This keep record of the start position of the LU
            chr_size = chr_size + syn9R_LU_size[abs(solution[LU])]
            LU_pos.append(chr_size)     # This keep record of the end position of the LU
            LU_pos_dic[LU] = [solution[LU], LU_pos]     # Put the start and end position in a dictionary
            #print(solution[LU], LU_pos)
        return LU_pos_dic
    else:
        # If there are multiple solution, it put all the dictionaries in a list
        LU_pos_dic_list = []
        for sol in solution:
            sol_L = position_LU(sol)
            LU_pos_dic_list.append(sol_L)
        return LU_pos_dic_list

#A = [1,2,2,3,4,6,13,1,5,9,16,4,7,4,4,5,6,-3,7,8,8,2,8,-11,8,9,10,11,11,12]
#B = [[3,4,6,13,1,5,9,16,4,7,4],[8,8,2,8,-11,8,9,10,11], [8,8,2,8,1]]
#syn9R = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
#print(solution_to_bps_size(A))
#print(solution_to_bps_size(B))
#print(solution_to_bps_size(syn9R))
#print(position_LU(A))
#print(position_LU(B))
#print(position_LU(syn9R))

#--------------------
# These following scripts are not very useful as I did not use them

# Find the equation between two points
def equation2points(A, B):
    # find the slope m
    #m = (yA - yB)/(xA - xB)
    m = (A[1]-B[1]) / (A[0]-B[0])
    # find the equation
    #y − y1 = m(x − x1)
    c = -m * A[0] + B[0]
    #equation = "y = "+str(m)+"x + "+str(c)
    #print(equation)
    return m, c

def point_on_equation(A, B, y):
    m, c = equation2points(A,B)
    # find the x
    x = (y-c)/m
    return x

#A=[3,9]
#B=[5,6]
#m = equation2points(A,B)
#print("m =", m)
#pp = point_on_equation(A,B,10)
#print("pp =", pp)

#points = [10,9,8,7,6,5,4,3,2,1]

# This function try to find the half value of a set, but it removes something first. I want to remove the essential LUs and see how many different non-essential LUS remaing. This to calculate the non-essential half-life.
def half_plus(start, plus=10, count=10):
    H=start
    half_list=[H]
    for _ in range(count):
        #H = (H / 2) + plus
        H = ((H - plus)/ 2) + plus
        half_list.append(H)
    return half_list

#print(half_plus(100, 10, 10))
#print(half_plus(50, 10, 10))
#print(half_plus(20, 10, 10))

# this is from the function 1/x with some parameters. I am trying to find an equation and curve that fits our results.
def func(x, a, b, c, d=0):
    y = (c/((x+d)**a))+b
    return y
"""
a = 1.5
b = 10
c = 80000
d = 160
equ=[]
for i in range(1, 1011, 10):
    Y = func(i, a, b, c, d)
    equ.append(Y)
    #print(Y)
print(equ)
"""