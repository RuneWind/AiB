import sys
from Bio import Phylo


'''
Returns a dictionary of the form {clade: list of names of clade's leaves} with all internal nodes in tree as keys in dictionary
e.g. {Clade(name = "A"): ["B", "C"], Clade(name = "D"): ["E", "F"], ...}
'''
def find_clade_and_leaves_dict(tree, root):
    # All clades in tree
    clades = tree.find_clades()    
    clade_and_leaves_dict = {}
    # Find internal clades in tree and the leaves in the subtrees for which they are root and put these into dict
    for clade in list(clades): # Depths first seach (important, else it doesn't work)
        if not clade.is_terminal(): # if nonterminal, it is internal; add as key to dict
            clade_and_leaves_dict[clade] = [] 
        if clade.is_terminal(): # if terminal, it is leaf; add as value to dict for all keys that are its parents
            if clade.name != root.name: # omit the root from dict keys
                parents = tree.get_path(clade)
                for parent in parents:
                    if not parent.is_terminal():
                        clade_and_leaves_dict[parent].append(clade.name)
    return clade_and_leaves_dict


'''
Annotates internal nodes in a tree with their DF-intervals
Input: Dictionary of form {clade: list of clade's leaves} 
and dictionary of form {"A": 1} where A is a leaf name and 1 is the corresponding index
'''
def annotate_internal_nodes_with_DF_intervals(tree_clade_and_leaves_dict, leaf_name_to_index_dict):
    for int_clade in tree_clade_and_leaves_dict.keys():
        leaf_names =  tree_clade_and_leaves_dict[int_clade]
        leaf_values = [leaf_name_to_index_dict[leaf_name] for leaf_name in leaf_names]
        if leaf_values != []:
            # Find min and max leaf values
            min_leaf_val = min(leaf_values) 
            max_leaf_val = max(leaf_values) 
            # Rename internal clade to form [min, max]
            int_clade.name = str(min_leaf_val) + "." + str(max_leaf_val)


'''
Returns the Robinson-Foulds distance (RF distance) by running Day's Algorithm
'''
def days_algo(tree1, tree2):
    # ***** Step 1 *****
    # ***** Root the 2 trees at the same leaf *****
    leaf_to_root = tree1.get_terminals()[0] # we root them at the "first" leaf in tree1
    tree1.root_with_outgroup({'name': str(leaf_to_root)})
    tree2.root_with_outgroup({'name': str(leaf_to_root)})
    
    
    # ***** Step 2 and 3 *****
    # ***** Make depth first numbering of leaves in tree1 and rename leaves in tree2 correspondingly *****
    depth_first_leaves = tree1.get_terminals()
    depth_first_leaves.remove(leaf_to_root) #remove root from list
    leaf_name_to_index_dict = {depth_first_leaves[i].name : i+1 for i in range(0, len(depth_first_leaves))}
    index_to_leaf_name_dict = {i+1: depth_first_leaves[i].name for i in range(0, len(depth_first_leaves))} # Delete? Maybe we don't need this

    
    # ***** Step 4.1 *****
    # ***** Annotate internal nodes in tree1 with their depth-first intervals *****
    tree1_clade_and_leaves_dict = find_clade_and_leaves_dict(tree1, leaf_to_root)
    annotate_internal_nodes_with_DF_intervals(tree1_clade_and_leaves_dict, leaf_name_to_index_dict)

    
    # ***** Step 4.2 ***** 
    # ***** Annotate internal nodes in tree2 with their depth-first intervals, IF the leaves of 
    # the subtree indeed IS an interval (we annotate all internal nodes and then remove non-intervals) *****
    tree2_clade_and_leaves_dict = find_clade_and_leaves_dict(tree2, leaf_to_root)
    annotate_internal_nodes_with_DF_intervals(tree2_clade_and_leaves_dict, leaf_name_to_index_dict)
    # Remove non-interval nodes
    for int_clade in tree2_clade_and_leaves_dict.keys():
        if int_clade.name != None:
            min_max_list = [int(v) for v in int_clade.name.split('.')]
            min_val = min_max_list[0]
            max_val = min_max_list[1]
            if max_val - min_val + 1 != len(tree2_clade_and_leaves_dict[int_clade]):
                int_clade.name = None

    
    # ***** Step 5 ***** 
    # ***** Count how many internal "interval nodes" the 2 trees share and calculate RF-distance *****
    # We don't look at the first "outer clade" and the clade that is "splitting" the root from the rest
    tree1_clades_to_compare = tree1.get_nonterminals()[2:]
    tree2_clades_to_compare = tree2.get_nonterminals()[2:]
    
    # Sort intervals in tree1
    tree1_clade_list = [float(c.name) for c in tree1_clades_to_compare if c.name]
    tree1_sorted_int_clades = sorted(tree1_clade_list)

    # Sort intervals in tree2    
    tree2_clade_list = [float(c.name) for c in tree2_clades_to_compare if c.name]
    tree2_sorted_int_clades = sorted(tree2_clade_list)

    # Compare the intervals in tree1 and tree2
    intersect = [c for c in tree1_sorted_int_clades if c in tree2_sorted_int_clades]
    #intersect = set(tree1_sorted_int_clades).intersection(tree2_sorted_int_clades)
        
    # Return Robinson-Foulds Distance
    # RF_cist = "number of intervals not found in both trees" = (intervals in tree1 - intersect) + (intervals in tree2 - intersect)
    rfdist = len(tree1_clades_to_compare) + len(tree2_clades_to_compare) - 2*len(intersect)
    return rfdist
    
    


'''
Code to run
'''

tree1 = Phylo.read(sys.argv[1], 'newick')
tree2 = Phylo.read(sys.argv[2], 'newick')

print("The RF-distance of ", sys.argv[1], " and ", sys.argv[2], " is ", days_algo(tree1, tree2))



'''
tree1 = Phylo.read('tree1.new', 'newick')
tree2 = Phylo.read('tree2.new', 'newick')
print("tree1:")
Phylo.draw(tree1)
print("tree2:")
Phylo.draw(tree2)
print("RF-dist of tree1 and tree2:", days_algo(tree1, tree2))

print("\n\n")

T1 = Phylo.read('T1.new', 'newick')
T2 = Phylo.read('T2.new', 'newick')
print("T1:")
Phylo.draw(T1)
print("T2:")
Phylo.draw(T2)
print("RF-dist of T1 and T2:", days_algo(T1, T2))
'''
