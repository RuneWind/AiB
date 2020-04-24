import numpy as np
from Bio import Phylo


"""
HELPER FUNCTION
Reads a file of this format:
4 
A  10  2  5  2
C  2  10  2  5
G  5  2  10  2
T  2  5  2  10

representing a substitution matrix and returns a dictionary corresponding to the substitutionmatrix

Returning a dictionary of this format (if getAlphabet = False):
{"A": {"A": 10, "C": 2, "G": 5, "T": 2}, 
 "C": {"A": 2, "C": 10, "G": 2, "T": 5}, 
 "G": {"A": 5, "C": 2, "G": 10, "T": 2}, 
 "T": {"A": 2, "C": 5, "G": 2, "T": 10}}

If getAlphabet = True, we instead return a list of the alphabet letters:
['A', 'C', 'G', 'T']
"""

def parse_phylip(filename, getAlphabet = False):
    f= open(filename, "r")    
    f1 = f.readlines()
    f2 = list()
    for x in f1:
        f2.append(x.split())
    alph_size = int(f2[0][0])
    
    letters = list()
    for i in range(1, alph_size+1):
        letters.insert(i, f2[i][0])
    '''
    sub_dict = dict()
    for i in range(len(letters)):
        inner_dict = dict()
        for j in range(len(letters)):
            inner_dict[letters[j]] = f2[i+1][j+1]
        sub_dict[letters[i]] = inner_dict
    '''
    dist_matrix = np.zeros((len(letters), len(letters)))  
    for i in range(len(letters)):
        inner_list = np.zeros(len(letters))
        for j in range(len(letters)):
            inner_list[j] = f2[i+1][j+1]
        dist_matrix[i] = inner_list
    
    index_to_letter_dict = dict()
    for i in range(len(letters)):
        index_to_letter_dict[i] = letters[i]
    if(getAlphabet):
        return letters
    else:
        return dist_matrix, index_to_letter_dict
    
    
def get_d_and_rs(dist_matrix, S, i, j):
    d_ij = dist_matrix[i][j]
    r_i = 1/(len(S)-2) * sum([dist_matrix[i][m] for m in S])
    r_j = 1/(len(S)-2) * sum([dist_matrix[j][m] for m in S])
    return d_ij, r_i, r_j
    
    
def nj(dist_matrix, index_to_letter_dict):
    S_nodes = list(index_to_letter_dict.values())
    S = index_to_letter_dict.keys()
    while len(S) > 3:
        
        # Step 1 (a): Compute matrix N
        # Matrix N entry (i,j) is a number indicating "how cloase i and j are to each other and how far they are from the rest"
        N = np.zeros((len(S), len(S)))
        for i in range(len(S)):
            inner_list = dist_matrix[i]
            for j in range(len(inner_list)):
                d_ij, r_i, r_j = get_d_and_rs(dist_matrix, S, i, j)
                N[i][j] = d_ij - (r_i + r_j)
        
        # Step 1 (b): Find minimum entry in matrix N
        N_removed_diag = [N[i][j] for i in S for j in S if i != j]
        min_val = min(N_removed_diag)
        min_i_j = np.where(N == min_val)[0]
        i = int(min_i_j[0])
        j = int(min_i_j[1])
        
        # Step 2 and 3: Add new node k to the tree with weights
        d_ij, r_i, r_j = get_d_and_rs(dist_matrix, S, i, j)
        weight_ki = (1/2) * (d_ij + r_i - r_j)
        weight_kj = (1/2) * (d_ij + r_j - r_i) #d_ij - weight_ki equals (1/2) * (d_ij + r_j - r_i)
        leaf_1 = S_nodes[i]
        leaf_2 = S_nodes[j]
        k = "(" + leaf_1 + ":" + str(round(weight_ki, 3)) + ", " + leaf_2 + ":" + str(round(weight_kj, 3)) + ")"
        
        
        # Step 4: Update dist_matrix (delete rows i and j and columns i and j and add new row and column for node k)
        
        # Create row for new node k to insert
        row = []
        for m in range(len(dist_matrix)):
            if(m != i and m != j):
                d_km = (1/2) * (dist_matrix[i][m] + dist_matrix[j][m] - d_ij)
                row.append(d_km)        
        
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
        
    # We now have 3 leaves left: i, j and m
    # Add a new node v to the tree with weights
    i = S[0]
    j = S[1]
    m = S[2]
    weight_vi = (1/2) * (dist_matrix[i][j] + dist_matrix[i][m] - dist_matrix[j][m])        
    weight_vj = (1/2) * (dist_matrix[i][j] + dist_matrix[j][m] - dist_matrix[i][m])        
    weight_vm = (1/2) * (dist_matrix[i][m] + dist_matrix[j][m] - dist_matrix[i][j])       

    v = "(" + S_nodes[i] + ":" + str(round(weight_vi, 3)) + ", " + S_nodes[j] + ":" + str(round(weight_vj, 3)) + ", " + S_nodes[m] + ":" + str(round(weight_vm, 3)) + ")" 

    return v


def print_tree_to_file(tree):
    x = open("tree.newick", "w")
    x.write(tree)
    x.close()



'''
Code to run
'''        
# Read distance matrix
distance_matrix, index_to_letter_dict = parse_phylip("example_slide4.phy")
letters = parse_phylip("example_slide4.phy", True)

# Construct tree
tree = nj(distance_matrix, index_to_letter_dict) 
print_tree_to_file(tree)

# Draw tree
tree1 = Phylo.read("tree.newick", 'newick')
Phylo.draw(tree1, branch_labels=lambda c: c.branch_length)















