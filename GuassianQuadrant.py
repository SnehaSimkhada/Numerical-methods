#!/usr/bin/env python
# coding: utf-8

# In[8]:


#Sneha Simkhada
#Assignment 3
#Math 392

import numpy as np

def gauss_points(n):

    if n == 2:

        c = [1,1]

        r = [0.5773502692, -.5773502692]

        return (c,r)

    if n == 3:

        c = [.55555555, .888888888, .55555555]

        r = [.7745966692, 0, -.7745966692]

        return (c,r)

    if n == 4:
        
        c = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]
        
        r = [0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116]
        
        return (c,r)
    
    if n == 5:
        
        c = [0.2369268850,0.4786286705, 0.5688888889,0.4786286705,0.2369268850]
        
        r = [0.9061798459, 0.5384693101, 0.0000000000, -0.5384693101, -0.9061798459]
        
        return (c,r)

    else:

        print('Value not supported.')
        
#define function to be integrated on [a,b]
def f(x):
    return x**2


n = 5

#define a and b
a = 0
b = 1


#Get the Gauss points
[c,r] = gauss_points(n)

#initialize area
area = 0 

#be sure to transform to [-1,1]
def gauss_quad(f, a, b, n, r, c):
    x = np.zeros(n)
    for i in range(n):
        x[i] = (b+a)/2 + (b-a)/2 *r[i]
    return (b-a)/2 * (c[0]*f(x[0]) + c[1]*f(x[1]) + c[2]*f(x[2]) + c[3]*f(x[3])+c[4]*f(x[4]))

#compute area
area = gauss_quad(f,a,b,n,r,c)
print("Gaussian integral is ", area)


# In[ ]:





# In[ ]:




