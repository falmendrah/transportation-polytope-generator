import numpy as np
from itertools import product

#########################################################
class plane_sum_entry_forbidden:
    """
    This class has the information of a plane-sum entry-forbidden 
    transportation polytope. 
    
    Attributes 
        u: vector of [23] margin sums
        v: vector of [13] margin sums
        w: vector of [12] margin sums
        Enabled: List of enabled cells 
    """
    def __init__(self, u, v, w, Enabled):
        self.u = u
        self.v = v
        self.w = w
        self.Enabled = Enabled
    
    def get_bounds(self):
        """
        Given the 2-margin sums and the enabled cells, it returns a 
        3D array with the bounds on each entry. 
        """
        U = self.u[0]
        r = len(self.u)
        h = len(self.w)
        bounds = np.zeros((r,r,h))
        for a,b,c in self.Enabled:
            bounds[a,b,c] = U
        return bounds


    def get_integer_points(self):
        r, c, l = len(self.u), len(self.v), len(self.w)
        
        given_bounds = self.get_bounds()
        bounds = np.zeros((r,c,l))
        
        for i in range(r):
            for j in range(c):
                for k in range(l):
                    sum_bound = min(self.u[i], self.v[j], self.w[k])
                    bounds[i,j,k] = min(sum_bound, given_bounds[i,j,k])
        
        ranges = [range(int(j)+1) for j in bounds.ravel()]
        integer_points = []

        for index,vector in enumerate(product(*ranges)):
            candidate = np.array(vector).reshape(r,c,l)
            if self.verify_plane_sums(candidate):
                integer_points.append(candidate)

        return(integer_points)
    
    def verify_plane_sums(self, x):
        r, c, l = len(self.u), len(self.v), len(self.w)

        flag = True
        
        for i in range(r):
            if(not flag):
                break
            if np.sum(x[i,:,:]) != self.u[i]:
                flag = False
        
        for j in range(c):
            if(not flag):
                break
            if np.sum(x[:,j,:]) != self.v[j]:
                flag = False
                
        for k in range(l):
            if(not flag):
                break
            if np.sum(x[:,:,k]) != self.w[k]:
                flag = False
                
        return(flag)
#######################################################            
        

#######################################################                 
def as_plane_sum(M, U):
    """
    Input:
        - M: Array M=(A|b) representing a bounded polytope {x>=0 : Ax=b} 
        - U: Upper bound for the entries of the polytope represented by M
    Output:
        - Object of the class "plane_sum_entry_forbidden" encoding the
          information of M into a plane-sum entry-forbidden 
          transportation polytope
    """
    
    A = M[:, :-1]
    b = M[:, -1]
    b_arr = np.array(b).reshape(-1)
    nrow, ncol = A.shape
        
    r = []
    for j in range(ncol):
        pos_sum = sum([A[i,j] for i in range(nrow) if A[i,j]>0])
        neg_sum = sum([A[i,j] for i in range(nrow) if A[i,j]<0])
        r.append(max(pos_sum, abs(neg_sum))) 
    
    h = nrow+1
    
    # I might want to reshape u and v so they have shape 1xr
    u = np.full(sum(r), U)
    v = np.full(sum(r), U)
    
    w = []
    for k in range(nrow):
        w.append(b[k] + U*sum([abs(A[k,j]) for j in range(ncol) if A[k,j]<0])) 
    w.append(sum(r)*U -sum(w))
    w = np.array(w)

    # For every k in range(ncol), Y_pos[k] and Y_neg[k] will be the set 
    # of positive and negative (respectively) enabled triplets 
    # corresponding to the variable y_{k+1}. 
    # 'Enabled' is the set list of enabled triplets
    
    # We define the 1st 2 entries of the enabled triplets
    Enabled = []
    Y_pos = [[] for k in range(ncol)]
    Y_neg = [[] for k in range(ncol)]
    marker = 0
    for k in range(ncol):
        if r[k] > 1:
            for i in range(r[k]):
                Y_pos[k].append([marker, marker, None])
                if i == 0:
                    Y_neg[k].append([marker+r[k]-1, marker, None]) 
                else:
                    Y_neg[k].append([marker-1, marker, None])
                marker +=1
        else:
            Y_pos[k].append([marker, marker, None])
            Y_neg[k].append([marker, marker, None])
            marker +=1
            
        Enabled +=(Y_pos[k]+Y_neg[k])
        
    # We define the 3rd entry of the enabled triplets
    for i in range(nrow):
        for j in range(ncol):
            if A[i,j]>0:
                prop_enabled = Y_pos[j]
            else:
                prop_enabled = Y_neg[j]
            n = abs(A[i,j])
            
            counter = 0
            ind_counter = 0
            while counter < n:
                if prop_enabled[ind_counter][2] is None:
                    prop_enabled[ind_counter][2] = i
                    counter +=1
                ind_counter +=1
                
    for triplet in Enabled:
        if triplet[2] is None:
            triplet[2] = nrow
    
    return(plane_sum_entry_forbidden(u, v, w, Enabled))
#######################################################