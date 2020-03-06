import numpy as np
import project2_linear as pa #pairwise alignment

# Example from slide 12 (Extend M3 with A)
A = ["ACG-T", "ACGGT"]
M = ["A--CGT", "ATTC-T", "CT-CGA", ""]

string_A = "AATAAT"
string_B = "AAGG"

sub_matrix = {"A": {"A": 10, "C": 2, "G": 5, "T": 2}, 
            "C": {"A": 2, "C": 10, "G": 2, "T": 5}, 
            "G": {"A": 5, "C": 2, "G": 10, "T": 2}, 
            "T": {"A": 2, "C": 5, "G": 2, "T": 10}}
gap_cost = -5

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


print(extend_M_with_A(M, A))

print(find_center_string(["AG", "AA", "GG", "CC", "TT"]))