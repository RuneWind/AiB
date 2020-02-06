# -*- coding: utf-8 -*-

import numpy as np

def read_fasta_file(filename):
    for record in SeqIO.parse(filename, "fasta"):  
        return record

#substitition matrix
#sub_mtrx = np.array([[10, 2, 5, 2], [2, 10, 2, 5], [5, 2, 10, 2], [2, 5, 2, 10]])
#print(sub_mtrx)

sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, 
            "C": {"A": 2, "C": 10, "G": 2, "T": 5}, 
            "G": {"A": 5, "C": 2, "G": 10, "T": 2}, 
            "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = -5

#Calculate cost of an optimal alignment for string str_A and str_B with substitution matrix sm and gap cost gp
def cost_of_optimal_alignment(sm, gp, str_A, str_B):
    #create array with None values
    T = np.full((len(str_A), len(str_B)), None)
    #iterate through rows
    for i in range(0, len(str_A)):
        #iterate through columns
        for j in range(0, len(str_B)):
            T[i,j] = calc_cost(i, j, T, str_A, str_B, sm, gp)
    return T[len(str_A) - 1, len(str_B) - 1]

#Calculate cost of one cell
def calc_cost(i, j, T, str_A, str_B, sm, gp):
    if(T[i,j] is None):
        cost = float("-inf")
        #get diagonal value
        if(i > 0 and j > 0):
            cost = calc_cost(i-1, j-1, T, str_A, str_B, sm, gp) + sm[str_A[i]][str_B[j]]
        #get above value
        if(i > 0 and j >= 0):
            cost = max(cost, calc_cost(i-1, j, T, str_A, str_B, sm, gp) + gp)
        #get left value
        if(i >= 0 and j > 0):
            cost = max(cost, calc_cost(i, j-1, T, str_A, str_B, sm, gp) + gp)
        #Left top corner
        if(i == 0 and j == 0):
            cost = max(cost, 0)
        return cost    
    return T[i,j]

#Find an optimal alignment based on an alignment matrix, T
#def backtrack(T, str_A, str_B, sm, gp):
#    last_cell = T[len(str_A) - 1, len(str_B) - 1]
#    upper = 


#Question 1
#String A and B
#string_A = "AATAAT"
#string_B = "AAGG"
#print(cost_of_optimal_alignment(sub_matrix, gap_cost, string_A, string_B))
#Optimal alignment: 10


#Question 2
fasta1 = read_fasta_file("seq1.fasta")
fasta2 = read_fasta_file("seq2.fasta")
print(cost_of_optimal_alignment(sub_matrix, gap_cost, fasta1.seq, fasta2.seq))
#Optimal alignment: 1336


 



