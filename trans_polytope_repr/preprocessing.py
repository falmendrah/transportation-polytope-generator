import numpy as np

#################################################
def prep_rep(M):
    """
    This function takes a standard representation of a polytope P and 
    returns the representation of the polytope Q as in step 3.1 of [De Loera and Onn, 2006]

    Input:
         - Array M=(A|b) where P = {y>=0 : Ay=b} 
    Output:
        - Array (C|d) where Q = {x>=0 : Cx=d} as in step 3.1 in paper
    """    
    A = M[:, :-1]
    b = M[:, -1]
    
    nrow, ncol = A.shape

    k = []

    for j in range(ncol):
        # The next variable is the j-th column but we replace 0s by 1s 
        col_ones = np.where(A[:, j] == 0, 1, abs(A[:, j]))
        bit_col = np.floor(np.log2(col_ones))
        k.append(int(np.max(bit_col)))

    nrow_new = nrow + np.sum(k)
    ncol_new = ncol + np.sum(k)
    C = np.zeros((nrow_new, ncol_new))

    # This defines the first sum(k) rows of the matrix C
    init_col = 0
    init_row = 0
    for i in range(len(k)):
        for j in range(k[i]):
                C[init_row,init_col] = 2
                C[init_row,init_col+1] = -1
                init_row +=1
                init_col +=1
        init_col +=1
        
    # This defines the last ncol rows of the matrix C    
    for i in range(nrow):
        init_col = 0
        for j in range(ncol):
            sign = np.sign(A[i,j])
            for s, t in enumerate(np.binary_repr(abs(A[i,j]))[::-1]):
                C[init_row, init_col+s] = sign*int(t)
            
            init_col +=(k[j]+1)
        init_row +=1
        
    d = np.concatenate((np.zeros(sum(k)),b))
    d = d[..., None]
    M = np.concatenate((C, d), axis=1)

    return(M.astype('int32'))
#########################################################