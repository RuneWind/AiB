import numpy as np


sub_m = {"A": {"A": 0, "C": 5, "G": 2, "T": 5}, 
            "C": {"A": 5, "C": 0, "G": 5, "T": 2}, 
            "G": {"A": 2, "C": 5, "G": 0, "T": 5}, 
            "T": {"A": 5, "C": 2, "G": 5, "T": 0}}
gc = 5

#Calculate cost of an optimal alignment for string str_A and str_B with substitution matrix sm and gap cost gc
def calculate_alignment_matrix(sub_m, gap_cost, strA, strB, strC):
    # Global vars 
    global str_A
    global str_B
    global str_C
    global sm
    global gc 
    global T 
    
    # Set global vars
    str_A = strA
    str_B = strB
    str_C = strC
    sm = sub_m
    gc = gap_cost
    
    T = np.full((len(str_A) + 1, len(str_B) + 1, len(str_C) + 1), None)
    #iterate through rows
    for i in range(0, len(str_A) + 1):
        #iterate through columns
        for j in range(0, len(str_B) + 1):
            for k in range(0, len(str_C) + 1):
                T[i, j, k] = calc_cost_nonrec(i, j, k)
    #print(T[len(str_A),len(str_B)])
    return T

def calc_cost_nonrec(i, j, k):
    if(T[i,j,k] is None):
        v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")
        #get diagonal value
        if(i > 0 and j > 0 and k > 0):
            v1 = T[i-1, j-1, k-1] + sm[str_A[i-1]][str_B[j-1]] + sm[str_A[i-1]][str_C[k-1]] + sm[str_B[j-1]][str_C[k-1]]
        #get above value
        if(i > 0 and j > 0 and k>=0):
            v2 = T[i-1, j-1, k] + sm[str_A[i-1]][str_B[j-1]] + 2*gc
        #get left value
        if(i > 0 and j >= 0 and k>0):
            v3 = T[i-1, j, k-1] + sm[str_A[i-1]][str_C[k-1]] + 2*gc
        if(i >= 0 and j > 0 and k>0):
            v4 = T[i, j-1, k-1] + sm[str_B[j-1]][str_C[k-1]] + 2*gc 
        if(i > 0 and j >= 0 and k >= 0):
            v5 = T[i-1, j, k] + 2*gc 
        if(i >= 0 and j > 0 and k >= 0):
            v6 = T[i, j-1, k] + 2*gc 
        if(i >= 0 and j >= 0 and k > 0):
            v7 = T[i, j,k-1] + 2*gc 
        #Left top corner
        if(i == 0 and j == 0 and k ==0):
            v0 = 0
        min_val =  min(v0, v1, v2, v3, v4, v5, v6, v7)
        return min_val  
    else:
        return T[i,j,k]

def backtrack_nonrec(T, str_A, str_B, str_C, i, j, k):
    print("Backtracking nonrecursively")
    res_str_A = ""    
    res_str_B = ""
    res_str_C = ""
    while(i >= 0 and j >= 0 and k>= 0):
        cell = T[i, j, k]
        #diagonal cell - substitution
        if (i > 0 and j > 0 and k > 0 and cell == T[i-1, j-1, k-1] + sm[str_A[i-1]][str_B[j-1]] + sm[str_A[i-1]][str_C[k-1]] + sm[str_B[j-1]][str_C[k-1]]):
            res_str_A += str_A[i-1]
            res_str_B += str_B[j-1]
            res_str_C += str_C[k-1]
            i -= 1
            j -= 1
            k -= 1
        #upper cell - insertion    
        elif (i > 0 and j > 0 and k>= 0 and cell == T[i-1, j-1, k] + sm[str_A[i-1]][str_B[j-1]] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += str_B[j-1]
            res_str_C += "-"
            i -= 1
            j -= 1
        
        elif (i > 0 and j >= 0 and k> 0 and cell == T[i-1, j, k-1] + sm[str_A[i-1]][str_C[k-1]] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += "-"
            res_str_C += str_C[k-1]
            i -= 1
            k -= 1
        elif (i >= 0 and j > 0 and k> 0 and cell == T[i, j-1, k-1] + sm[str_B[j-1]][str_C[k-1]] + 2*gc):
            res_str_A += "-"
            res_str_B += str_B[j-1]
            res_str_C += str_C[k-1]
            j -= 1
            k -= 1

        elif (i > 0 and j >= 0 and k>= 0 and cell == T[i-1, j, k] + 2*gc):
            res_str_A += str_A[i-1]
            res_str_B += "-"
            res_str_C += "-"
            i -= 1
        elif (i >= 0 and j > 0 and k>= 0 and cell == T[i, j-1, k] + 2*gc):
            res_str_A += "-"
            res_str_B += str_B[j-1]
            res_str_C += "-"
            j -= 1
        elif (i >= 0 and j >= 0 and k> 0 and cell == T[i, j, k-1] + 2*gc):
            res_str_A += "-"
            res_str_B += "-"
            res_str_C += str_C[k-1]
            k -= 1

        elif (i==0 and j==0 and k==0):
            #write resulting alignment to fasta file
            #x = open("alignment.fasta", "w")
            #x.write(">seq1\n" + res_str_A[::-1] + "\n\n" + ">seq2\n" + res_str_B[::-1])
            #x.close()
            return [res_str_A[::-1], res_str_B[::-1], res_str_C[::-1]]


'''
strA = "GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGATAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTCTCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTG"
strB = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAACGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA"
strC = "CGCTGGTGCAACTCGAAGACCTATCTCCTTCCCGGGGGGGCTTCTCCGGCATTTAGGCCTCGGCGTTTGGAAGTACGGAGGTTTTTCTCGGAAGAAAGTTCACTGGAAGTGGAAGAAATGGATTTATCTGCTGTTCGAATTCAAGAAGTACAAAATGTCCTTCATGCTATGCAGAAAATCTTGGAGTGTCCAATCTGTTT"
'''

strA = "GTTCCGAAAGGCTAGCGCTAGGCGCC"
strB = "ATGGATTTATCTGCTCTTCG"
strC = "TGCATGCTGAAACTTCTCAACCA"

D = calculate_alignment_matrix(sub_m, gc, strA, strB, strC)
print(D[len(strA), len(strB), len(strC)])

b = backtrack_nonrec(T, strA, strB, strC, len(strA), len(strB), len(strC))
print(b)
