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
def days_algo(tree1, tree2):
    # Step 1
    leaf_to_root = tree1.get_terminals()[0]
    tree1.root_with_outgroup({'name': str(leaf_to_root)})
    tree2.root_with_outgroup({'name': str(leaf_to_root)})
    Phylo.draw(tree1)
    #Phylo.draw(tree2)
    
    for clade in tree1.find_clades(order='level'):
        print(clade)
        for c in clade.__iter__():
            print("sub: ", c)
    
    # Step 2 (and 3)
    depth_first_leaves = tree1.get_terminals()
    depth_first_leaves.remove(leaf_to_root) #remove root from list
    leaf_name_to_index_dict = {i+1: depth_first_leaves[i] for i in range(0, len(depth_first_leaves))}
    leaf_index_to_name_dict = {depth_first_leaves[i] : i+1 for i in range(0, len(depth_first_leaves))}
    #print(leaf_name_to_index_dict)
    #print(leaf_index_to_name_dict)
    
    # Step 4.1
    #subtrees = Phylo.subtrees(tree1, wait=False)
    #for subtree in subtrees:
     #   min_leaf = min([leaf_name_to_index_dict[name_of_leaf.name] for name_of_leaf in subtree.get_terminals()])
      #  print(min_leaf)
    
    # Name internal nodes
    
    
    
    
days_algo(tree1, tree2)