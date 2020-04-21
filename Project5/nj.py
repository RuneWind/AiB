import numpy as np


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
    
 
def nj(dist_matrix, index_to_letter_dict):
    S = index_to_letter_dict.values()
    while len(S) > 3:
        min_val = dist_matrix.min()
        return min_val
    
#print(parse_phylip("example_slide4.phy"))
distance_matrix, index_to_letter_dict = parse_phylip("example_slide4.phy")
letters = parse_phylip("example_slide4.phy", True)



print(nj(distance_matrix, index_to_letter_dict))















