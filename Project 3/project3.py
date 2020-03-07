import numpy as np
import project2_linear as pa #pairwise alignment

# Example from slide 12 (Extend M3 with A)
#A = ["ACG-T", "ACGGT"]
#M = ["A--CGT", "ATTC-T", "CT-CGA", ""]
S = ["ACGT", "ATTCT", "CTCGA"]

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
    # List for holding new M-strings
    new_M = ["" for i in range(len(M))]
    # String holding the last M-string for the new M-list
    new_M_str = ""
    # Iterate through columns in A (opt alignment)
    for j in range(len(A[0])):
        # Upper and lower symbol in the j'th column
        upper_symbol = A[0][j]
        lower_symbol = A[1][j]
        # Case: insertion
        if upper_symbol == "-":
            # We want to add a column to M consisting of only "-"'s except 
            # the last row, which should contain lower_symbol:
            # Add "-" to new M-strings (except for last M-string)
            new_M = [old_str + "-" for old_str in new_M]
            # Add lower_symbol to last M-string
            new_M_str += lower_symbol
        # Case: (mis)match or deletion
        else:
            # Find first occurrance of upper_symbol in first M-string 
            # (we have cut off the part the string that we have already added to the new M-strings)
            pos = M[0].find(upper_symbol)
            # Add everything up until this symbol to all the new M-strings (except the last)
            prefixes_M = [m_str[:pos+1] for m_str in M]
            new_M = [new_M[i] + prefixes_M[i] for i in range(len(new_M))]
            # To the last string, add correspondingly many gaps (minus 1) and lower_symbol
            new_M_str += "-"*pos + lower_symbol
            # Update M to only contain suffix that we have not yet "transferred" to new_M yet
            M = [m_str[pos+1:] for m_str in M]
    # Add the last M-string to the updated M-strings in new_M
    new_M.append(new_M_str)
    # Return entire new M
    return new_M

def multiple_align(S, center):
    M = []
    S.remove(center)
    for s in S:
        A_matrix = pa.calculate_alignment_matrix(sub_matrix, gap_cost, center, s)
        # optimal alignment
        A = pa.backtrack(A_matrix, center, s, sub_matrix, gap_cost, "", "", len(center), len(s))
        if(s != S[0]):
            M = extend_M_with_A(M, A)
        else:
            M = A
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

