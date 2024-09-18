import numpy as np
from itertools import product

#######################################################
class slim_line_sum:
    """
    This class has the information of a slim line-sum 
    transportation polytope. 
    
    Attributes 
        U: 2-margin fixing entries i,j
        V: 2-margin fixing entries i,k
        W: 2-margin fixing entries j,k
    """
    def __init__(self, U, V, W):
        self.U = U
        self.V = V
        self.W = W

    def get_integer_points(self, relaxed_coord=None, all=False, n_relax=-1):
        """
        Input:
            - lower_bounds: array with the lower bounds of the integer
              points in my transportation polytope 
        """
        r, c = self.U.shape
        l = len(self.W[0])

        if(all):
            low_bounds = np.full((r,c,l), -1)
        else:
            low_bounds = np.zeros((r,c,l))

        if(relaxed_coord is not None):
            for coord in relaxed_coord:
                i,j,k = coord
                low_bounds[i,j,k] = n_relax

        up_bounds = np.zeros((r,c,l))
        for i in range(r):
            for j in range(c):
                for k in range(l):
                    up_bounds[i,j,k] = min(self.U[i,j], self.V[i,k], self.W[j,k])

        ranges = [range(int(a),int(b)+1) for a,b in zip(low_bounds.ravel(), up_bounds.ravel())]
        # This should be temporary
        sizes = [b-a+1 for a,b in zip(low_bounds.ravel(), up_bounds.ravel())]
        # print(sizes)

        integer_points = []

        # This is temporary to keep track of number of candidates
        counter = 0
        print('We are trying to do ', np.prod(sizes), ' iterations.')
        print('Checking integer points in T ...')

        for vector in product(*ranges):
            counter += 1
            candidate = np.array(vector).reshape(r,c,l)
            if self.verify_line_sums(candidate):
                integer_points.append(candidate)

        print(counter)
        return(integer_points)

    def verify_line_sums(self, x):
        r, c = self.U.shape
        l = len(self.W[0])

        flag = True

        for i in range(r):
            if(not flag):
                break
            for j in range(c):
                if np.sum(x[i,j,:]) != self.U[i,j]:
                    flag = False
                    break

        for i in range(r):
            if(not flag):
                break
            for k in range(l):
                if np.sum(x[i,:,k]) != self.V[i,k]:
                    flag = False
                    break

        for j in range(c):
            if(not flag):
                break
            for k in range(l):
                if np.sum(x[:,j,k]) != self.W[j,k]:
                    flag = False
                    break

        return(flag)
#######################################################


#######################################################
def as_slim_line_sum(P):
    """
    Input
        - P: A polytope in the form of plane-sum restricted entries
    Output
        - U: array of line-sums slicing through K-axis 
        - V: array of line-sums slicing through J-axis
        - W: array of line-sums slicing through I-axis
    """
    
    a_margin = P.u
    b_margin = P.v
    c_margin = P.w
    E = P.get_bounds()
    
    l = len(a_margin)
    m = len(b_margin)
    n = len(c_margin)
    
    r = l*m
    c = n+l+m
        
    U_bound = min(max(a_margin), max(b_margin))
    
    U = np.zeros((r,c))
    V = np.zeros((r,3))
    W = np.zeros((c,3))
    
    I = [(i,j) for i in range(l) for j in range(m)]
    J = [(0,t) for t in range(n)]+[(1,t) for t in range(l)]+[(2,t) for t in range(m)]
    K = [0,1,2]
    
    # Defining the rxc array U
    for row_idx, row_label in enumerate(I):
        i,j = row_label[0], row_label[1]
        
        for col_idx, col_label in enumerate(J[:n]):
            t = col_label[1]
            U[row_idx, col_idx] = E[i,j,t]
            
        for col_idx, col_label in enumerate(J[n:n+l], n):
            t = col_label[1]
            if t == i:
                U[row_idx, col_idx] = U_bound
            
        for col_idx, col_label in enumerate(J[n+l:], n+l):
            t = col_label[1]
            if t == j:
                U[row_idx, col_idx] = U_bound
                                
    # Defining the rx3 array V
    for row_idx, row_label in enumerate(I):
        i,j = row_label[0], row_label[1]
        for col_idx in range(3):
            if col_idx == 1:
                V[row_idx, col_idx] = sum(E[i, j, :])
            else:
                V[row_idx, col_idx] = U_bound
    
    # Defining the cx3 array W
    for row_idx, row_label in enumerate(J):
        i,j = row_label[0], row_label[1]
        
        if i == 0:
            W[row_idx, 0] = c_margin[j]
            W[row_idx, 1] = np.sum(E[:,:,j]) - c_margin[j] 
        
        if i == 1:
            W[row_idx, 0] = m*U_bound - a_margin[j]
            W[row_idx, 2] = a_margin[j]
        
        if i ==2:
            W[row_idx, 1] = a_margin[j] 
            W[row_idx, 2] = l*U_bound - b_margin[j]
            
    return(slim_line_sum(U, V, W))
#######################################################