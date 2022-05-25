#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import math as m
from numpy.linalg import matrix_power
from numpy import linalg as la


# In[2]:


def xCC_4(L):
    
    X = []
    Y = []
    
    for x in np.arange(-0.9,0.9,0.05):
        
        t = 1/(m.sqrt(m.exp(2*x)+1))
        r = m.sqrt(1-t**2)
        q = 3
        l = L//q
        A0 = np.array([[1/t, r/t, 0, 0], [r/t, 1/t ,0 ,0],[0, 0, 1/t, r/t], [0, 0, r/t, 1/t]])
        B0 = np.array([[1/r, 0, 0, t/r], [0, 1/r, t/r, 0], [0, t/r, 1/r, 0], [t/r, 0, 0, 1/r]])
        C0 = B0@A0
        U = matrix_power(C0,50)
        #U = np.identity(4)
        Tl = matrix_power(C0,q)

        u = np.zeros((l,4), dtype=float)
        Q0,R0 = la.qr(U)


        for i in range(0,l):
            #print(i+1," Q0 \n", Q0)
            Q,R = la.qr(Tl@Q0)
            u[i] = np.diagonal(R)
            ##print(i+1," Q \n", Q)
            #print(i+1," R \n", R)
            #print(i+1," u \n", u)]

            for j in range(0,4):
                if u[i,j] < 0:
                    u[i,j] = -u[i,j]
                    Q[:,j] = -Q[:,j]

            Q0 = Q
            #print(i+1," Q0 \n", Q0)

        nu = np.log(u)/L
        gamma = nu.sum(axis=0)
        gamma_ = gamma.sum()
        #print(gamma,gamma_)
        
        #print("Valores analíticos: \n",np.log((1+t)*(1+r)/(r*t)),np.log((1-t)*(1+r)/(r*t)),np.log((1+t)*(1-r)/(r*t)),np.log((1-t)*(1-r)/(r*t)))
        X.append(x)
        Y.append(gamma[1])
    plt.plot(X,Y,label = str(L))
    plt.legend(title = "Valores de L")


# In[3]:


def LCC_4(x):
    
    X = []
    Y = []
    
    for L in range(3,90):
        
        t = 1/(m.sqrt(m.exp(2*x)+1))
        r = m.sqrt(1-t**2)
        q = 3
        l = L//q
        A0 = np.array([[1/t, r/t, 0, 0], [r/t, 1/t ,0 ,0],[0, 0, 1/t, r/t], [0, 0, r/t, 1/t]])
        B0 = np.array([[1/r, 0, 0, t/r], [0, 1/r, t/r, 0], [0, t/r, 1/r, 0], [t/r, 0, 0, 1/r]])
        C0 = B0@A0
        U = matrix_power(C0,50)
        #U = np.identity(4)
        Tl = matrix_power(C0,q)

        u = np.zeros((l,4), dtype=float)
        Q0,R0 = la.qr(U)


        for i in range(0,l):
            #print(i+1," Q0 \n", Q0)
            Q,R = la.qr(Tl@Q0)
            u[i] = np.diagonal(R)
            ##print(i+1," Q \n", Q)
            #print(i+1," R \n", R)
            #print(i+1," u \n", u)]

            for j in range(0,4):
                if u[i,j] < 0:
                    u[i,j] = -u[i,j]
                    Q[:,j] = -Q[:,j]

            Q0 = Q
            #print(i+1," Q0 \n", Q0)

        nu = np.log(u)/L
        gamma = nu.sum(axis=0)
        gamma_ = gamma.sum()
        #print(gamma,gamma_)
        
        #print("Valores analíticos: \n",np.log((1+t)*(1+r)/(r*t)),np.log((1-t)*(1+r)/(r*t)),np.log((1+t)*(1-r)/(r*t)),np.log((1-t)*(1-r)/(r*t)))
        X.append(L)
        Y.append(gamma[1])
    plt.plot(X,Y,label = str(x))
    plt.legend(title = "Valores de x")


# In[4]:


xCC_4(9)
xCC_4(15)
xCC_4(30)
xCC_4(60)
xCC_4(120)

plt.xlabel("Parâmetro de energia x")
plt.ylabel("Expoentes críticos?")
plt.show()

LCC_4(0.2)
LCC_4(0.3)
LCC_4(0.4)
LCC_4(0.5)
    
plt.xlabel("Comprimento da rede L")
plt.ylabel("Expoentes críticos?")
plt.show()


# In[37]:


L= 90
x= 0.3
q = 3
l = L//q
t = 1/(m.sqrt(m.exp(2*x)+1))
r = m.sqrt(1-t**2)
z1 = np.random.random(4) + np.random.random(4) * 1j
z2 = np.random.random(4) + np.random.random(4) * 1j
Vl = np.diag(z1)
Ul = np.diag(z2)
A0 = np.array([[1/t, r/t, 0, 0], [r/t, 1/t ,0 ,0],[0, 0, 1/t, r/t], [0, 0, r/t, 1/t]])
B0 = np.array([[1/r, 0, 0, t/r], [0, 1/r, t/r, 0], [0, t/r, 1/r, 0], [t/r, 0, 0, 1/r]])
C0 = B0@Vl@A0@Ul
U = matrix_power(C0,50)
#U = np.identity(4)
Tl = matrix_power(C0,q)
u = np.zeros((l,4), dtype=float)
Q0,R0 = la.qr(U)


for i in range(0,l):
    #print(i+1," Q0 \n", Q0)
    Q,R = la.qr(Tl@Q0)
    u[i] = np.diagonal(R)
    ##print(i+1," Q \n", Q)
    #print(i+1," R \n", R)
    #print(i+1," u \n", u)]

    for j in range(0,4):
        if u[i,j] < 0:
            u[i,j] = -u[i,j]
            Q[:,j] = -Q[:,j]

    Q0 = Q
    #print(i+1," Q0 \n", Q0)

nu = np.log(u)/L
gamma = nu.sum(axis=0)
gamma_ = gamma.sum()
print(gamma,gamma_)

print("Valores analíticos: \n",np.log((1+t)*(1+r)/(r*t)),np.log((-t+1)*(1+r)/(r*t)),np.log((1+t)*(-r+1)/(r*t)),np.log((-1+t)*(-1+r)/(r*t)))
#X.append(L)
#Y.append(gamma[0])
#plt.plot(X,Y,label = str(x))
#plt.legend(title = "Valores de x")


# In[ ]:





# In[ ]:





# In[ ]:




