#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Sneha Simkhada
#JACOBI INTERGRATION
#The system of equations was reordered to form a strictly diagonally dominant system, because we cannot have a zero at the diagonal.

import numpy as np

def jacobi(A,b,x):
    maxIt = 1000
    n = len(b)
    tol = 10**(-2)
    Dinv = np.zeros((n,n))
    L = np.zeros((n,n))
    U= np.zeros((n,n))
   
    
    #Assemble Dinv and c
    for i in range(n):
        for j in range(n):
            if i == j:
                Dinv[i,j] = 1/A[i,j]
            elif i < j:  #upper Traingle
                U[i,j] = -A[i,j]
            elif i>j:   #lower Traingle
                L[i,j] = -A[i,j]
    
    T = np.matmul(Dinv,L+U)
    C = np.matmul(Dinv,b)
    
    for i in range(maxIt):
        xnew = np.matmul(T,x) + C
        err = max(abs(xnew - x))
        if err <= tol:
            print('Solution is', xnew, 'after', i, 'iterations.')
            return xnew
        elif i == maxIt -1:
            print('System did not converge')   
        else:
            x =  xnew
        
    print(xnew)
    
A = np.array(([-1,0,0,np. sqrt(2)/2,1,0,0,0 ], 
              [0,-1,0,np. sqrt(2)/2,0,0,0,0 ], 
              [0,0,-1,0,0,0,1/2,0],
              [0,0,0,-np. sqrt(2)/2,0,-1,-1/2,0],
              [0,0,0,0,-1,0,0,1],
              [0,0,0,0,0,1,0,0],
             [0,0,0,-np. sqrt(2)/2,0,0,np. sqrt(3)/2,0],
             [0,0,0,0,0,0,-np. sqrt(3)/2,-1]))

b = np.array([0,0,0,0,0,10000,0,0])
n = len(b)
x = [1,1,1,1,1,1,1,1]
soln = jacobi(A,b,x)


# In[ ]:




