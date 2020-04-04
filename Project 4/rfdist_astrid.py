# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 12:41:31 2020

@author: qeast
"""

from Bio import Phylo

tree1 = Phylo.read('tree1.new', 'newick')
tree2 = Phylo.read('tree2.new', 'newick')
'''
Phylo.draw(tree1)
for clade in tree1.find_clades(order='level'):
    print(clade.name)


print(tree1)
print(tree2)
Phylo.draw(tree1)
Phylo.draw(tree2)

#tree1(rooted = True)
#print("TEST", tree1.get_terminals()[0])

#tree1.root_with_outgroup({'name': 'seq1'})  # Operates in-place
#print (tree1)
#Phylo.draw(tree1)

def find_internal_nodes(tree):
    int_nodes = []
    for clade in tree.find_clades(order='level'):
        #print(clade)
        #for child in clade:
      #      int_nodes.append(clade)
    #return parents

#find_internal_nodes(tree1)
'''


def find_clade_and_leaves_dict(tree):
    # Find internal clades in tree and the leaves in the subtrees for which they are root
    clades = tree.find_clades()    
    # Put these into dictionary, e.g. {Clade(branch_length=0.02441): ['seq10', 'seq7'], Clade(branch_length=0.06209): ['seq8', 'seq6'], ...}
    clade_and_leaves_dict = {}
    for clade in list(clades): # Depths first seach (important, else it doesn't work)
        if not clade.is_terminal(): # if nonterminal, it is internal; add as key to dict
            clade_and_leaves_dict[clade] = [] 
        if clade.is_terminal(): # if terminal, it is leaf; add as value to dict for all keys that are its parents
            parents = tree.get_path(clade)
            for parent in parents:
                if not parent.is_terminal():
                    clade_and_leaves_dict[parent].append(clade.name)
    return clade_and_leaves_dict


#Annotate internal nodes in tree1 with their DF-intervals 
def annotate_internal_nodes_by_DF_interval(tree_clade_and_leaves_dict, leaf_name_to_index_dict):
    for int_clade in tree_clade_and_leaves_dict.keys():
        leaf_names =  tree_clade_and_leaves_dict[int_clade]
        leaf_values = [leaf_name_to_index_dict[leaf_name] for leaf_name in leaf_names]
        if leaf_values != []:
            # Find min and max leaf values
            min_leaf_val = min(leaf_values) 
            max_leaf_val = max(leaf_values) 
            # Rename internal clade to form [min, max]
            #int_clade.name = "[" + str(min_leaf_val) + ", " + str(max_leaf_val) + "]"
            int_clade.name = str(min_leaf_val) + "." + str(max_leaf_val)
            # For testing, remove later:
            #print(leaf_names)
            #print(leaf_values)
            #print(int_clade.name)
            #print("************")


def days_algo(tree1, tree2):
    Phylo.draw(tree1)
    Phylo.draw(tree2)
    # Step 1
    leaf_to_root = tree1.get_terminals()[0]
    tree1.root_with_outgroup({'name': str(leaf_to_root)})
    tree2.root_with_outgroup({'name': str(leaf_to_root)})
    #Phylo.draw(tree1)
    #Phylo.draw(tree2)
    
    
    # Step 2 (and 3)
    depth_first_leaves = tree1.get_terminals()
    depth_first_leaves.remove(leaf_to_root) #remove root from list
    index_to_leaf_name_dict = {i+1: depth_first_leaves[i].name for i in range(0, len(depth_first_leaves))}
    leaf_name_to_index_dict = {depth_first_leaves[i].name : i+1 for i in range(0, len(depth_first_leaves))}
    #print(leaf_name_to_index_dict)
    #print(index_to_leaf_name_dict)
    
    
    # Step 4.1
    tree1_clade_and_leaves_dict = find_clade_and_leaves_dict(tree1)
    annotate_internal_nodes_by_DF_interval(tree1_clade_and_leaves_dict, leaf_name_to_index_dict)
    Phylo.draw(tree1)

    
    # Step 4.2
    tree2_clade_and_leaves_dict = find_clade_and_leaves_dict(tree2)
    annotate_internal_nodes_by_DF_interval(tree2_clade_and_leaves_dict, leaf_name_to_index_dict)
    # Remove non-interval nodes that are not an interval
    for int_clade in tree2_clade_and_leaves_dict.keys():
        if int_clade.name != None:
            min_max_list = [int(v) for v in int_clade.name.split('.')]
            min_val = min_max_list[0]
            max_val = min_max_list[1]
            if max_val - min_val + 1 != len(tree2_clade_and_leaves_dict[int_clade]):
                int_clade.name = None
    Phylo.draw(tree2)

    # Step 5
    # Count how many internal "interval nodes" from tree1 is also in tree2
    
    # Sort intervals in tree1
    tree1_clade_list = [float(c.name) for c in tree1_clade_and_leaves_dict.keys() if c.name]
    tree1_sorted_int_clades = sorted(tree1_clade_list)
    print(tree1_sorted_int_clades)
    
    # Sort intervals in tree2    
    tree2_clade_list = [float(c.name) for c in tree2_clade_and_leaves_dict.keys() if c.name]
    tree2_sorted_int_clades = sorted(tree2_clade_list)
    print(tree2_sorted_int_clades)

    # Compare the intervals in tree1 and tree2
    intersect = set(tree1_sorted_int_clades).intersection(tree2_sorted_int_clades)
    #print(intersect) 
    #print(len(intersect))  
    
    # Return Robinson-Foulds Distance
    # "number of intervals not found in both trees" = (intervals in tree1 - intersect) + (intervals in tree2 - intersect)
    return len(tree1_sorted_int_clades) + len(tree2_sorted_int_clades) - 2*len(intersect)
    





print(days_algo(tree1, tree2))