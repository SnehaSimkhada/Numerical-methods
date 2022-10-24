#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Sneha Simkhada
# Math 392
#Program 2


import numpy as np
import matplotlib.pyplot as plt


#forward difference for the first point
def forward(new_list,f,h):
    y=(f[1]-f[0])/h  #formula: y = (f(x0 + h) - f(x0))/h
    new_list.append(y)
    return new_list

#central difference for the all the center points point
def central(list_1,f,h,n):
    for i in range(1,n-1):
        cd = (f[i+1]-f[i-1])/(2*h)  #formula: y = (f(x0 + h) - f(x0-h))/2h
        list_1.append(cd)
    return list_1


#backward difference for the last point
def backward(f,list_2,h,n):
    bd=(f[-1]-f[-2])/(h)  #formula: y = (f(x0) - f(x0-h))/h
    list_2.append(bd)
    return list_2  
        
x = [805, 825, 845, 865, 885, 905, 925, 945, 965, 985]
f = [0.710, 0.763, 0.907, 1.336, 2.169, 1.598, 0.916, 0.672, 0.615, 0.606]    
h = x[1]-x[0]
n = len(x)
new_list=[]
list_1 = forward(new_list,f,h)
list_2= central(list_1,f,h,n)
list_3  = backward(f,list_2,h,n)
print(list_3)

plt.plot(x,f)
plt.plot(x,f, label = "Given function")
plt.plot(x,list_3)
plt.plot(x,list_3, label = "Derivative function")
plt.grid()
plt.legend()
plt.show()


# In[ ]:




