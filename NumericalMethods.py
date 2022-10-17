#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# //example1
import math as m
a=-2
b=1.5  

for i in range (maxit):
    p = 1/2*(a+b)
    fa = 3*(a+1)*(a-0.5)*(a-1)
    fp = 3*(p+1)*(p-0.5)*(p-1)
    if fa*fp<0:
        b = p
    else:
        a = p
print(p);



# In[4]:


# //example2
import math as m
a=0
b=1 
maxit= 17
for i in range (maxit):
    p = 1/2*(a+b)
    fa = m.exp(a)- a**2 +3*a-2
    fp = m.exp(p)- p**2 +3*p-2
    if fa*fp<0:
        b = p
    else:
        a = p
print(p)


# In[6]:


# example 3 Newtons Method conversions 
import math as m 
def newtonsMethod(p):
    maxIt = 1000
    tol = 10**(-5)
    for i in range (1,maxIt):
        fp = 2*p**3+p**2-p+1
        fpprime = 6*p**2+2*p-1
        pnew = p- fp/fpprime
        if m.fabs(p-pnew)<tol:
            p = pnew
            print('Num iterations:', i)
            return p
        else:
            p = pnew
    
p0 = -1.2
root = newtonsMethod(p0)
print(root)


# In[3]:


# example 4 Newtons Method diversions
import math as m 
def newtonsMethod(p):
    maxIt = 1000
    tol = 10**(-5)
    for i in range (1,maxIt+1):
        fp = p**(1/3)
        fpprime = 1/(3*p**(2/3))
        pnew = p- fp/fpprime
        if m.fabs(p-pnew)<tol:
            p = pnew
            print('Num iterations:', i)
            return p
        elif i == maxIt:
            print('Method Diverges')
            return []
        p = pnew
    
p0 = -1.2
root = newtonsMethod(p0)
print(root)


# In[19]:


# example 5 Secant Method conversions
import math as m 

def f(x):
    return 2*x**3+x**2-x+1
    
def secantMethod(p0, p1):
    maxIt = 5
    tol = 10**(-5)
    for i in range (1,maxIt+1):
        fp0 = f(p0)
        fp1 = f(p1)
        pnew = p1 - fp1*(p1-p0)/(fp1-fp0)
        if m.fabs(p1-pnew)<tol:
            p = pnew
            print('Num iterations:', i)
            return p
        elif i == maxIt+1:
            print('Method Diverges')
            return []
        p0 = p1
        p1 = pnew
    
p0 = -1.2
p1 = -1.3
root = secantMethod(p0,p1)
print(root)


# In[9]:


# example 5 hornors method 

def horners(p0, P):
    #initializing Q to be an n-1 list of zeros
    n = len(P)
    Q = (n-1)*[0]
    Q[0] = P[0]
    for i in range( 1, n):
        if i == n-1:
            R = Q[i-1]*p0 + P[i]
        else:
            Q[i] = Q[i-1]*p0 + P[i]
    print('Q', Q, 'R', R)
    return R, Q

p0 = -2
P = [2, 0, -3, 3, -4]
[R,Q] = horners(p0, P)


# In[7]:


# hw Method conversions
import math as m 

def f(x):
    return x-2**-x
    
def steffensenMethod(p0):
    maxIt = 1000
    tol = 10**(-4)
    for i in range (1,maxIt+1):
        p1 = p0 + f(p0)
        p2 = p1 + f(p1)
        pnew = p2 - (pow((p2 - p1),2)/(p2 - (2*p1) + p0))
        if abs(pnew-p0)<tol:
            p = pnew
            print('Num iterations:', i)
            return p
        elif i == maxIt+1:
            print('Method Diverges')
            return []
        p0 = pnew
    
p0 = 0.5
root = steffensenMethod(p0)
print(root)


# In[8]:


#Neville's method

import numpy as np #numerical python

def neville(x, y, xr):
    #function for neville's method
    #initialize 
    n = len(x)
    Q = np.zeros([n,n])  #pass single argument of a list
    Q[:,0] = y #first col, all rows = y values
    for i in range(1,n):
        for j in range(1,i):
            Q[i,j] = ((xr -x[i-j])*Q[i,j-1]            -(xr -x[i])*Q[i-1, j-1])            / (x[i]- x[i-j])
    print(Q)
    return Q

xr = 8.4
x = [8.1, 8.3, 8.6, 8.7]
y = [16.9, 17.6, 18.5, 18.8]
Q = neville(x, y, xr)


# In[2]:


#Neville's method

import numpy as np #numerical python

def neville(x, y, xr):
    #function for neville's method
    #initialize 
    n = len(x)
    Q = np.zeros([n,n])  #pass single argument of a list
    Q[:,0] = y #first col, all rows = y values
    for i in range(1,n):
        for j in range(1,i):
            Q[i,j] = ((xr -x[i-j])*Q[i,j-1]            -(xr -x[i])*Q[i-1, j-1])            / (x[i]- x[i-j])
    print(Q)
    return Q

xr = 0
x = [-0.5,-0.25,0.25,0.5]
y = [1.9375,1.33203,0.800781,0.687500]
Q = neville(x, y, xr)


# In[3]:


#Central differenctial


import numpy as np
import matplotlib.pyplot as plt


#generates x value in [0,2pi]. Step size h
h = 0.1
x = np.arange(0,2*np.pi + h, h)
f = np.sin(x)
n = len(x)

#compute derivatives
#fwd diff
diffArr = np.diff(f)

#do backward diff at the end
backDiff = (f[n-1] - f[n-2])
fprime = np.append(diffArr,backDiff)
fprime = fprime/h


print('f',np.shape(f))
print('diffArr',np.shape(diffArr))
print(diffArr)

plt.plot(x,f)
plt.plot(x,fprime)
plt.grid()
plt.show()


# In[4]:


#trapeziodal rule
from scipy.integrate import trapz
import numpy as np
a =-1
b= 1
n =4
h = (b-a)/n
x = np.arange(a,b+h,h)
y = np.cos(x)**2
trapz(y,x)


# In[ ]:




