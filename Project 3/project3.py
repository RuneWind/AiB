import numpy as np
import project2_linear as pa #pairwise alignment

# Example from slide 12 (Extend M3 with A)
#A = ["ACG-T", "ACGGT"]
#M = ["A--CGT", "ATTC-T", "CT-CGA", ""]
S = ["ACGT", "ATTCT", "CTCGA"]

#string_A = "AATAAT"
#string_B = "AAGG"
"""
sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, 
            "C": {"A": 2, "C": 10, "G": 2, "T": 5}, 
            "G": {"A": 5, "C": 2, "G": 10, "T": 2}, 
            "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = -5
"""

sub_matrix = {"A": {"A": 0, "C": 5, "G": 2, "T": 5}, 
            "C": {"A": 5, "C": 0, "G": 5, "T": 2}, 
            "G": {"A": 2, "C": 5, "G": 0, "T": 5}, 
            "T": {"A": 5, "C": 2, "G": 5, "T": 0}}

gap_cost = 5

# Example from week 4 exercise 4.5 (Extend M4 with A)
#A = ["-AC-GT", "G-TAGT"]
#M = ["A--CG-T", "ATTC--T", "CT-CG-A", "A--CGGT", ""]


def find_center_string(S):
    for pos_cen in S: # possible center
        sum_scores_list = []
        sum_scores = 0
        for s in S:
            sum_scores += pa.calculate_alignment_matrix(sub_matrix, gap_cost, pos_cen, s)[len(pos_cen), len(s)]
        sum_scores_list.append(sum_scores)
    index_of_min_score = np.argmin(sum_scores_list)
    return S[index_of_min_score]

def extend_M_with_A(M, A):
    M.append("")
    old_pos = 0
    pos = 0
    for j in range(0, len(A[0])):
        # Case: (Mis)match or deletion
        if A[0][j] != "-":
            pos = M[0][old_pos:].find(A[0][j]) + old_pos
            gaps = pos - old_pos - 1
            old_pos = pos
            # Case: (Mis)match
            if A[1][j] != "-":
                M[-1] = M[-1] + "-"*gaps + A[1][j]
            # Case: Deletion
            else:
               M[-1] = M[-1][:pos] + "-"
        # Case: Insertion
        else:
            if pos == 0:
                for k in range(0,len(M)-1):
                   M[k] = "-" + M[k]
                M[-1] = A[1][j] 
            if pos > 0:
                for k in range(0,len(M)-1):
                   M[k] = M[k][:pos+1] + "-" + M[k][pos+1:]
                M[-1] = M[-1][:pos+1] + A[1][j] 
            pos = pos + 1
            old_pos = pos
    return M

def multiple_align(S, center):
    M = []
    S.remove(center)
    for s in S:
        print("s: ", s)
        print("M ", M)
        A_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, center, s)
        # optimal alignment
        A = pa.backtrack(A_matrix, center, s, sub_matrix, gap_cost, "", "", len(center), len(s))
        print("A ", A)
        if(s != S[0]):
            M = extend_M_with_A(M, A)
        else:
            M = A
        print("M ny ", M)
    return M

def print_alignment_to_file(seq_list):
    #write alignment list to fasta file
    x = open("alignment.fasta", "w")
    for i in range(len(seq_list)):    
        x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
    x.close()

S = ["GTTCCGAAAGGCTAGCGCTAGGCGCC", "ATGGATTTATCTGCTCTTCG", "TGCATGCTGAAACTTCTCAACCA"]
center = find_center_string(S)
seq_list = multiple_align(S, center)
print_alignment_to_file(seq_list)

