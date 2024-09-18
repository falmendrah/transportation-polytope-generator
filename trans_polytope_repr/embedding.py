import numpy as np
from .plane_sum import as_plane_sum
from .slim_line_sum import as_slim_line_sum
from .preprocessing import prep_rep

#######################################################
def embed_in_plane_sum(M, U, y):
    """
    Recieves an integer point y inside a polytope P = {y>=0 : Ay=b}
    and returns its embedding in a plane-sum transportation 
    polytope Q. 

    Input
        - y: Integer point in the polytope P
        - M: Array M = (A|b) that represents the polytope P
        - U: Upper bound for the entries of the polytope P
    Output    
        -d: a dictionary with the following information
            x: integer point in plane-sum polytope Q
            real_coordinates: image of the injection sigma (we are 
                getting just ONE copy per coordinate y_k)
            projected_point: image of y under the coordinate-erasing
                projection
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
    
    P = as_plane_sum(M,U)
    E = P.Enabled
    x = np.zeros((np.sum(r),np.sum(r),h))  

    x_proj = []
    real_coord = []
    marker = 0
    for k in range(ncol):
        # This part of the code is used when for each point y in P
        # we consider exactly one real coordinate associated to 
        # each y_k
        # real_coord.append(E[marker])
        # x_proj.append(y[k]) 
        if r[k] > 1:
            for i in range(r[k]):
                a,b,c = E[marker]
                x[a,b,c] = y[k]
                real_coord.append(E[marker])
                x_proj.append(y[k])
                marker +=1
            for i in range(r[k]):
                a,b,c = E[marker]
                x[a,b,c] = U-y[k]
                marker +=1
        else:
            a,b,c = E[marker]
            x[a,b,c] = y[k]
            real_coord.append(E[marker])
            x_proj.append(y[k])
            a,b,c = E[marker+1]
            x[a,b,c] = U-y[k]
            marker +=2

    d = dict()
    d['point'] = x
    d['real_coordinates'] = real_coord
    d['projected_point'] = x_proj 

    return(d)
#######################################################


#######################################################
def embed_plane_sum_in_line_sum(P, y, y_real_coordinates=None):
    """
    Recieves an integer point y inside a plane-sum entry-forbidden 
    polytope P and returns the corresponding integer point x inside the
    slim line-sum transportation polytope T
    
    Input
        - y: Integer point in a plane-sum, entry-forbidden 
             transportation polytope P
        - P: Plane-sum, entry-forbidden transportation polytope
    Optional input
        - y_real_coordinates: list of coordinates in y that come
          from a polytope that was embedded in P   
    Output
        -d: a dictionary with the following information
            x: integer point in slim line-sum polytope T
            real_coordinates: image of the injection sigma (we are 
                getting just one copy per coordinate y_k)
            projected_point: image of y under the coordinate-erasing
                projection
    """
    l, m, n = y.shape
    
    U = P.u[0]
    e = P.get_bounds()
    
    r = l*m
    c = n+l+m
    
    I = [(i,j) for i in range(l) for j in range(m)]
    J = [(0,t) for t in range(n)]+[(1,t) for t in range(l)]+[(2,t) for t in range(m)]
    K = [0,1,2]
    
    x = np.zeros((r,c,3))
    real_coord = []
    x_proj = []

    # Define slice corresponding to y (encoding y)
    for row_idx, row_label in enumerate(I):
        i, j = row_label[0], row_label[1]
        for k in range(n):
            x[row_idx, k, 0] = y[i,j,k]
            if(y_real_coordinates is None):
                real_coord.append([row_idx, k, 0])
                x_proj.append(y[i,j,k])
            else:
                if([i,j,k] in y_real_coordinates):
                    real_coord.append([row_idx, k, 0])
                    x_proj.append(y[i,j,k])
    
    # Define entries of the form x_{I,(2,t),1} and x_{I,(2,t),3}
    for row_idx, row_label in enumerate(I):
        i, j = row_label[0], row_label[1]
        for t in range(l):
            if t == i:
                x[row_idx, t+n, 0] = U-np.sum(y[i,j,:])
                x[row_idx, t+n, 2] = np.sum(y[i,j,:])
                
    # Define entries of the form x_{I,(1,t),2}         
    for row_idx, row_label in enumerate(I):
        i, j = row_label[0], row_label[1]
        for k in range(n):
            x[row_idx, k, 1] = e[i,j,k]-y[i,j,k]
            
    # Define entries of the form x_{I,(3,t),2} and x_{I,(3,t),3}         
    for row_idx, row_label in enumerate(I):
        i, j = row_label[0], row_label[1]
        for t in range(m):
            if t == j:
                x[row_idx, t+n+l, 1] = np.sum(y[i,j,:])
                x[row_idx, t+n+l, 2] = U-np.sum(y[i,j,:])
    
    d = dict()
    d['point'] = x
    d['real_coordinates'] = real_coord
    d['projected_point'] = x_proj 

    return(d)
#######################################################

#######################################################
def slim_line_sum_representation(P, upper_bound):
    """
    Input
        - P: Array encoding a convex polytope P = {x>=0: Ax=b} in standard form
        - upper_bound: Upper bound on the entries of the integers points inside P 
    Output
        'slim_line_sum' object encoding the slim line sum representation of P
        as described by [De Loera and Onn, 2006]  
    """
    P_updated = np.array(P)
    if np.max(P_updated[:, :-1]) > 2:
        P_updated = prep_rep(P)
    
    P_plane_sum = as_plane_sum(P_updated, upper_bound)
    P_slim_line_sum = as_slim_line_sum(P_plane_sum)

    return P_slim_line_sum 
#######################################################


#######################################################
def slim_line_sum_representation(P, upper_bound):
    """
    Input
        - P: Array encoding a convex polytope P = {x>=0: Ax=b} in standard form
        - upper_bound: Upper bound on the entries of the integers points inside P 
    Output
        'slim_line_sum' object encoding the slim line sum representation of P
        as described by [De Loera and Onn, 2006]  
    """
    P_updated = np.array(P)
    if np.max(P_updated[:, :-1]) > 2:
        P_updated = prep_rep(P)
    
    P_plane_sum = as_plane_sum(P_updated, upper_bound)
    P_slim_line_sum = as_slim_line_sum(P_plane_sum)

    return P_slim_line_sum 
#######################################################


#######################################################
def embed_in_line_sum(y, M, upper_bound):
    """
    Recieves an integer point y inside a polytope P = {y>=0 : Ay=b}
    and returns its embedding in a slim line sum transportation 
    polytope T. 

    Input
        - y: Integer point in the polytope P
        - M: Array M = (A|b) that represents the polytope P
        - upper_bound: Upper bound for the entries of the polytope P
    Output
        -d: a dictionary with the following information
            x: integer point in slim line-sum polytope T
            real_coordinates: image of the injection sigma (we are 
                getting just one copy per coordinate y_k)
            projected_point: image of y under the coordinate-erasing
                projection
    """
    M_updated = np.array(M)
    if np.max(M_updated[:, :-1]) > 2:
        M_updated = prep_rep(M)
    
    M_plane_sum = as_plane_sum(M_updated, upper_bound)
    plane_sum_y = embed_in_plane_sum(M_updated, upper_bound, y)

    line_sum_y = embed_plane_sum_in_line_sum(M_plane_sum, plane_sum_y['point'])

    return line_sum_y 
#######################################################