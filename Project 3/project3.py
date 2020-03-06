# Example from slide 12 (Extend M3 with A)
A = ["ACG-T", "ACGGT"]
M = ["A--CGT", "ATTC-T", "CT-CGA", ""]

# Example from week 4 exercise 4.5 (Extend M4 with A)
#A = ["-AC-GT", "G-TAGT"]
#M = ["A--CG-T", "ATTC--T", "CT-CG-A", "A--CGGT", ""]



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