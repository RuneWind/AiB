# -*- coding: utf-8 -*-

import numpy as np
from Bio import SeqIO
import sys 

def read_fasta_file(filename):
    for record in SeqIO.parse(filename, "fasta"):  
        return record



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
    
    sub_matrix = dict()
    for i in range(len(letters)):
        inner_dict = dict()
        for j in range(len(letters)):
            inner_dict[letters[j]] = int(f2[i+1][j+1])
        sub_matrix[letters[i]] = inner_dict
    if(getAlphabet):
        return letters
    else:
        return sub_matrix
        
    
#substitition matrix
#sub_mtrx = np.array([[10, 2, 5, 2], [2, 10, 2, 5], [5, 2, 10, 2], [2, 5, 2, 10]])
#print(sub_mtrx)


'''
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, 
            "C": {"A": 2, "C": 10, "G": 2, "T": 5}, 
            "G": {"A": 5, "C": 2, "G": 10, "T": 2}, 
           "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = -5
'''




#Calculate cost of an optimal alignment for string str_A and str_B with substitution matrix sm and gap cost gc
def calculate_alignment_matrix(sm, gc, str_A, str_B):
    #create array with None values to fill with cost values
    T = np.full((len(str_A) + 1, len(str_B) + 1), None)
    #create array with None values to fill with # of predecessors values    
    P = np.full((len(str_A) + 1, len(str_B) + 1), None)

    #iterate through rows
    for i in range(0, len(str_A) + 1):
        #iterate through columns
        for j in range(0, len(str_B) + 1):
            T[i,j], P[i,j] = calc_cost(i, j, T, P, str_A, str_B, sm, gc)
    print(T[len(str_A),len(str_B)])
    return T, P

#Calculate cost of one cell
def calc_cost(i, j, T, P, str_A, str_B, sm, gc):
    if(T[i,j] is None):
        pred = 0
        diag_cost = above_cost = left_cost = zero_cost = float("inf")
        #get diagonal value
        if(i > 0 and j > 0):
            diag_cost = calc_cost(i-1, j-1, T, P, str_A, str_B, sm, gc) + sm[str_A[i-1]][str_B[j-1]]
        #get above value
        if(i > 0 and j >= 0):
            above_cost = calc_cost(i-1, j, T, P, str_A, str_B, sm, gc) + gc
        #get left value
        if(i >= 0 and j > 0):
            left_cost = calc_cost(i, j-1, T, P, str_A, str_B, sm, gc) + gc
        #Left top corner
        if(i == 0 and j == 0):
            zero_cost = 0
            pred = 1

        min_val = min(diag_cost, above_cost, left_cost, zero_cost)  
        #Set number of predecessors

        if(diag_cost == min_val):
            pred += P[i-1,j-1] 
        
        if(above_cost == min_val):
            pred += P[i-1,j] 
        if(left_cost == min_val):
            pred += P[i,j-1] 
        return min_val, pred  
    return T[i,j]




#Find an optimal alignment based on an alignment matrix, T
#last cell: T[len(str_A) - 1, len(str_B) - 1]
def backtrack(T, str_A, str_B, sm, gc, res_str_A, res_str_B, i, j):
    cell = T[i, j]
    #print("i: " + str(i) + " j:" + str(j))
    #diagonal cell - substitution
    if (i > 0 and j > 0 and cell == T[i-1, j-1] + sm[str_A[i-1]][str_B[j-1]]):
        backtrack(T, str_A, str_B, sm, gc, res_str_A + str_A[i-1] , res_str_B + str_B[j-1], i-1, j-1)
    #upper cell - insertion    
    elif (i > 0 and j >= 0 and cell == T[i-1, j] + gc):
        backtrack(T, str_A, str_B, sm, gc, res_str_A + str_A[i-1], res_str_B + "-",  i-1, j)
    #left cell - deletion
    elif (i >= 0 and j > 0 and cell == T[i, j-1] + gc):
        backtrack(T, str_A, str_B, sm, gc, res_str_A + "-", res_str_B + str_B[j-1], i, j-1)
    elif (i==0 and j==0):
        #print(res_str_A[::-1] + "\n" + res_str_B[::-1])
        x = open("alignment.fasta", "w")
        x.write(">seq1\n" + res_str_A[::-1] + "\n\n" + ">seq2\n" + res_str_B[::-1])
        x.close()

#Question 1 (optimal alignment of short sequences)
#String A and B
#string_A = "AATAAT"
#string_B = "AAGG"
#t, p = calculate_alignment_matrix(sub_matrix, gap_cost, string_A, string_B)
#print(t)
#print(t[len(string_A), len(string_B)])
#Optimal alignment cost is 20


#Question 2 (optimal alignment of fasta files)
#fasta1 = read_fasta_file("seq1.fasta")
#fasta2 = read_fasta_file("seq2.fasta")
#t2, p2 = calculate_alignment_matrix(sub_matrix, gap_cost, fasta1.seq, fasta2.seq)
#print(t2)
#print(t2[len(fasta1.seq), len(fasta1.seq)])
#Optimal alignment cost is 1346


# Question 3 (backtracking)
#backtrack(t, string_A, string_B, sub_matrix, gap_cost, "", "", len(string_A), len(string_B))
#backtrack(t2, fasta1.seq, fasta2.seq, sub_matrix, gap_cost, "", "", len(fasta1.seq), len(fasta2.seq))
#Optimal alignment is written to the file "alignment.fasta"


#Question 4 (number of optimal alignments)
#print(p)
#there is one optimal alignment
#print(p2)
#there are 256 optimal alignments


#Test case example (there should be 4 optimal alignments)
#string_1 = "TCCAGAGA"
#string_2 = "TCGAT"
#t3, p3 = calculate_alignment_matrix(sub_matrix, gap_cost, string_1, string_2)
#print(p3)
#we get 4 optimal alignments

sub_matrix = parse_phylip(sys.argv[1])
gap_cost = int(sys.argv[2])
str_A = read_fasta_file(sys.argv[3]).seq.upper()
str_B = read_fasta_file(sys.argv[4]).seq.upper()


t3, p3 = calculate_alignment_matrix(sub_matrix, gap_cost, str_A, str_B)
print(sys.argv)
if len(sys.argv)==6:
	if sys.argv[5]=="True":
		
		b = backtrack(t3, str_A, str_B, sub_matrix, gap_cost, "", "", len(str_A), len(str_B))
	


# Get the sub matrix
#print(parse_phylip("sub_matrix2.txt"))
# Get the alphabet - we use this to check if our sequences contain only the letters in the alphabet
#print(parse_phylip("sub_matrix2.txt", True))