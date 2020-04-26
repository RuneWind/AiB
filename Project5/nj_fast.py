import sys
import numpy as np
from Bio import Phylo


'''
HELPER FUNCTION
Reads a file of this format:
5
A 0.00 0.23 0.16 0.20 0.17
B 0.23 0.00 0.23 0.17 0.24
C 0.16 0.23 0.00 0.20 0.11
D 0.20 0.17 0.20 0.00 0.21
E 0.17 0.24 0.11 0.21 0.00

representing a dissimilarity/distance matrix and returns: 

The distance matrix of form:
[[0.00 0.23 0.16 0.20 0.17]
[0.23 0.00 0.23 0.17 0.24]
[0.16 0.23 0.00 0.20 0.11]
[0.20 0.17 0.20 0.00 0.21]]

and 
   
A list of the letters / node names in the file:
['A', 'B', 'C', 'D', 'E']
'''
def parse_phylip_to_matrix_and_letters(filename):
    # Read file
    f= open(filename, "r")    
    f1 = f.readlines()
    f2 = list()
    for x in f1:
        f2.append(x.split())
    alph_size = int(f2[0][0])
    
    # Letters
    letters = list()
    for i in range(1, alph_size+1):
        letters.insert(i, f2[i][0])
        
    # Distance matrix
    dist_matrix = np.zeros((len(letters), len(letters)))  
    for i in range(len(letters)):
        inner_list = np.zeros(len(letters))
        for j in range(len(letters)):
            inner_list[j] = f2[i+1][j+1]
        dist_matrix[i] = inner_list
        
    return dist_matrix, letters


'''
Opens a file tree.newick and writes the input string to it
'''
def print_tree_to_file(tree):
    x = open("tree.newick", "w")
    x.write(tree)
    x.close()    
    
    
'''
Neighbor-joining algorithm
Input is a distance matrix for leaves in tree, letters are the leaf names
Returns a tree in newick format constructed from the input distance matrix
'''    
def nj(dist_matrix, letters):
    # We have a list of nodes "S", in each iteration of the algorithm
   # S_nodes are the actual names of our leaves and S is the corresponding indices in dist_matrix
    S_nodes = letters
    S = list(range(len(letters)))
    
    # While we have more than 3 nodes left in S
    while len(S) > 3:
        
        # Step 1 (a): Compute matrix N
        
        # Create list of row sums, i.e. r_i's, which we use to look up in
        row_sum_list = [1/(len(S)-2) * sum(dist_matrix[i][:]) for i in S]

        # Matrix N entry (i,j) is a number indicating "how cloase i and j are to each other and how far they are from the rest"        
        # N is symmetric, so we only fill out upper triangle (without diagonal, since we do not use it)
        # Fill out upper triangle in N matrix (rest is set to 0)
        N = [[dist_matrix[i][j] - (row_sum_list[i] + row_sum_list[j]) if j > i else 0 for j in S] for i in S] 
        
        # Step 1 (b): Find minimum entry in matrix N and its indices (i, j)
        N_upper = [N[i][j] for i in S for j in S if j > i] # only check upper triangle in N
        min_val = min(N_upper)
        min_i_j = np.where(N == min_val)
        listOfCordinates = list(zip(min_i_j[0], min_i_j[1]))
        min_coord = None
        for coord in listOfCordinates:
            if(coord[0] != coord[1]):
                min_coord = coord
                break
        i = int(min_coord[0])
        j = int(min_coord[1])
        
        
        # Step 2 and 3: Add new node k to the tree with weights
        
        # Calculate weights and find leaf names
        weight_ki = (1/2) * (dist_matrix[i][j] + row_sum_list[i] - row_sum_list[j])
        weight_kj = (1/2) * (dist_matrix[i][j] + row_sum_list[j] - row_sum_list[i]) 
        leaf_1 = S_nodes[i]
        leaf_2 = S_nodes[j]
        
        # Create Newick formatted node
        k = "(" + leaf_1 + ":" + str(round(weight_ki, 3)) + ", " + leaf_2 + ":" + str(round(weight_kj, 3)) + ")"
        
        
        # Step 4: Update dist_matrix (delete rows i and j and columns i and j and add new row and column for node k)
        
        # Create row for new node k to insert
        row = [(1/2) * (dist_matrix[i][m] + dist_matrix[j][m] - dist_matrix[i][j]) for m in S if m != i and m != j]
        
        # Remove row i and j
        dist_matrix = np.delete(dist_matrix, [i,j], 0)        
        # Remove column i and j
        dist_matrix = np.delete(dist_matrix, [i,j], 1)        

        # Insert row and column for new node k
        dist_matrix = np.row_stack((dist_matrix, row))
        dist_matrix = np.column_stack((dist_matrix, np.append(row, 0)))
        
        
        # Step 5: Delete i and j from S and add new node k to S
        S_nodes.remove(leaf_1)
        S_nodes.remove(leaf_2)
        S_nodes.append(k)
        S = list(range(len(S_nodes)))
        
    # Termination: We now have 3 leaves left: i, j and m
    # Add a new node v to the tree with weights
    i, j, m = S[0], S[1], S[2]
    weight_vi = (1/2) * (dist_matrix[i][j] + dist_matrix[i][m] - dist_matrix[j][m])        
    weight_vj = (1/2) * (dist_matrix[i][j] + dist_matrix[j][m] - dist_matrix[i][m])        
    weight_vm = (1/2) * (dist_matrix[i][m] + dist_matrix[j][m] - dist_matrix[i][j])       

    # Return Newick formatted tree
    v = "(" + S_nodes[i] + ":" + str(round(weight_vi, 3)) + ", " + S_nodes[j] + ":" + str(round(weight_vj, 3)) + ", " + S_nodes[m] + ":" + str(round(weight_vm, 3)) + ");" 
    return v


'''
Code to run
'''
# Read distance matrix
distance_matrix, letters = parse_phylip_to_matrix_and_letters(sys.argv[1])

# Construct tree and it print to file
tree = nj(distance_matrix, letters) 
print_tree_to_file(tree)

# Draw tree
#tree1 = Phylo.read("tree.newick", 'newick')
#Phylo.draw(tree1, branch_labels=lambda c: c.branch_length)


'''
Test code
'''
'''
# Read distance matrix
distance_matrix, letters = parse_phylip_to_matrix_and_letters("example_slide4.phy")

distance_matrix, letters = parse_phylip_to_matrix_and_letters("89_Adeno_E3_CR1.phy")

# Construct tree and it print to file
tree = nj(distance_matrix, letters) 
print_tree_to_file(tree)

# Draw tree
tree1 = Phylo.read("tree.newick", 'newick')
Phylo.draw(tree1, branch_labels=lambda c: c.branch_length)
'''













