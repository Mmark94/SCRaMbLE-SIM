from comparison_sol import str_path
from comparison_sol import count_solutions
from comparison_sol import remove_duplicate
#from overlap_lists import redundant
import pandas as pd
import subprocess
import shlex
import os
import time
import signal


# These functions are from short_reads_pipeline.py

def find_edges(LUs):
    sr_edges = []
    c = 0  # c = counter
    for i in range(len(LUs)):
        edge = str(c) + "_" + str(c + 1)
        c = c + 2
        sr_edges.append(edge)
    #print("sr_edges =", sr_edges)
    Dic_LUs_edges = dict(zip(LUs, sr_edges))
    #print("Dic_LUs_edges =", Dic_LUs_edges)
    return Dic_LUs_edges


# this invert (reverse) a node. For example, 2/1 becomes 1/2
def invert_node(node):
    pos_l = node.find("/")     # find the position of / in the node
    first = node[:pos_l]
    second = node[pos_l+1:]
    #print(first, second)
    inv_node = second + "/" + first
    return inv_node

# remove duplicates in the list of nodes
def reduce_nodes(nodes):
    new_nodes = []
    for i in nodes:
        if i not in new_nodes and invert_node(i) not in new_nodes:
            new_nodes.append(i)
    return new_nodes

# this make the smaller LU  first in the node. For example, 2/1 becomes 1/2
def sort_nodes(nodes):
    new_nodes = []
    for n in nodes:
        pos_l = n.find("/")  # find the position of / in the node
        first = n[:pos_l]
        second = n[pos_l + 1:]
        if int(first) <= int(second):
            new_nodes.append(n)
        else:
            new_nodes.append(invert_node(n))
    return new_nodes

# Find the Nodes or the LoxP junctions
def find_nodes_and_edges(path, debug=False):
    # Find the edges in the reference
    #path_abs = [abs(i) for i in path]
    MAX = abs(max(path, key=abs))   #this find the biggest number by absolute value
    ref = range(1, MAX+1)
    ref_edges = find_edges(ref)

    # Convert the path in the edges
    path_edges = []
    for LU in path:
        if LU >= 0:
            edge = ref_edges[LU]
        else:
            # if the LU is negative invert the edge
            pos_ = ref_edges[-LU].find("_")
            edge = ref_edges[-LU][pos_ + 1:] + "_" + ref_edges[-LU][:pos_]
        path_edges.append(edge)

    # Count the edges
    Dic_edges_counter = {}
    for LU in ref_edges:
        Dic_edges_counter[ref_edges[LU]] = 0

    for edge in path_edges:
        try:
            Dic_edges_counter[edge] = Dic_edges_counter[edge] + 1
        except:
            # it means that the LU is inverted (negative) and that the edge is inverted
            pos_ = edge.find("_")
            edge = edge[pos_ + 1:] + "_" + edge[:pos_]
            Dic_edges_counter[edge] = Dic_edges_counter[edge] + 1

    # Find the nodes
    path_nodes = []
    first = "0"
    for i in path_edges:
        # find position of "_"
        pos_ = i.find("_")
        second = i[:pos_]
        node = first + "/" + second
        path_nodes.append(node)
        first = i[pos_+1:]
    # remove first element
    path_nodes.pop(0)
    # remove duplicates in the list of nodes
    new_path_nodes = reduce_nodes(path_nodes)
    new_path_nodes = sort_nodes(new_path_nodes)     # I am not sure this is essential, maybe can be removed
    if debug == True:
        print("ref_edges =", ref_edges)
        print("path_edges =", path_edges)
        print("Dic_edges_counter =", Dic_edges_counter)
        print("path_nodes =", path_nodes)
    return Dic_edges_counter, new_path_nodes,

###############################################

# These are new functions

# Find novel loxP junctions.
def find_novel_loxP_junctions(path, ref=[], debug=False):
    if ref == []:
        # Find the reference
        #path_abs = [abs(i) for i in path]
        MAX = abs(max(path, key=abs))   #this find the biggest number by absolute value
        ref = range(1, MAX+1)

    # find loxP_junctions in the reference
    nodes_ref = find_nodes_and_edges(ref, debug=False)[1]

    # find loxP_junctions in the SCRaMbLEd chromosome
    nodes_SCRaMbLEd = find_nodes_and_edges(syn, debug=False)[1]

    # subtract loxP_junctions in the SCRaMbLEd chromosome to the loxP_junctions in the reference
    novel_nodes = []
    for N in nodes_SCRaMbLEd:
        if N not in nodes_ref:
            novel_nodes.append(N)

    if debug:
        print(find_nodes_and_edges(ref))
        print("nodes_ref       =", nodes_ref)
        print("nodes_SCRaMbLEd =", nodes_SCRaMbLEd)
        print("novel_nodes     =", nodes_SCRaMbLEd)
    return novel_nodes


# Classifying recombination events


# test the code
if __name__ == "__main__":
    ref = [1,2,3,4,5]
    syn = [1,2,5,-3,2]
    #print(find_nodes_and_edges(syn, debug=False))
    print(find_novel_loxP_junctions(syn, ref=[], debug=True))
