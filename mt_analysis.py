#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 17:13:26 2020

@author: elemhunt
"""
'''
  MT model
  elastoplastic behaviour
  metal matrix reinforced with elastic spheres
'''
#%% Imports
import numpy as np
import math
import matplotlib.pyplot as plt
 
#%% Constants
k_i = 222.2     # k index in GPa
k_m = 62.5      # k matrix in GPa
mu_i = 166.67   # mu index in GPa
mu_m = 28.85    # mu matrix in GPa
f_i = 0.999      # volume fraction (between 0-1)
tol = 1E-6      # tolerance for the calculation
A =0.4          # 400MPa in GPa
N = 0.15 

#%% J and K Matrices

J = np.zeros([(6),(6)])
for i in range(3):
    for j in range(3):
        J[i,j] = (1./3)
        
K = np.zeros([(6),(6)])
for i in range(3):
    for j in range(3):
        K[i,j] = -(1./3)
for i in range(3):
    K[i,i] = (2./3)
for i in range(3,6):
    K[i,i] = 1

#%% Calculations
# Inital stress (Sigma_o) and Uniaxial Stress vector
Sigma_VEC0 =np.zeros([6])
Sigma_VEC0[0]=0.001
# Inclusion Matrix
M_i = (1/(3.*k_i))*J +(1/(2.*mu_i))* K
# Final solution storage lists
epsilon_L =[]
Sigma_L =[]
## MT calculation M_sec  
# Calculating 150 points
for n in range(150):
    M_M_sec = np.ones([6,6]) # inital M Matrix secant  value
    Sigma_EQ = 1.0
    Sigma_Old = 0.0
    while abs((Sigma_EQ - Sigma_Old)> tol): # performance of the secant method   
        # MT calculations parameters 
        k_P = (3*k_m + 4*mu_m)/3
        mu_P = ((5/3)*mu_m*(3*k_m + 4*mu_m)/(k_m + 2*mu_m))/2 
        MT_k = k_m +(f_i*(k_i - k_m))/(1.+(1-f_i)*(k_i - k_m)/k_P)
        MT_mu = mu_m +(f_i*(mu_i - mu_m))/(1.+(1-f_i)*((mu_i - mu_m)/mu_P))
        # M and B Matrix secant calculations via MT
        M_sec = (1/(3.*MT_k))*J +(1/(2*MT_mu))* K
        B_msec = (1./(1.-f_i))*(((M_M_sec-M_i)**-1)*(M_sec-M_i))
        Sigma_Old =Sigma_EQ
        # Matrix and Equivalent stress calculations
        Sigma_M = np.dot(B_msec,Sigma_VEC0)
        Sigma_M_Pr = np.dot(K,Sigma_M)    
        Sigma_EQ = math.sqrt((3/2)*np.dot( Sigma_M_Pr, Sigma_M_Pr))
        # Equivalent strain calculation
        Eps_EQ = ((Sigma_EQ)/A)**(1/N)
        # Mu Matrix secant and M Matrix secant calculations
        mu_M_sec = (Sigma_EQ)/((Sigma_EQ/mu_m)+(3*Eps_EQ))
        M_M_sec = (1/(3.*k_m))*J +(1/(2.*mu_M_sec))* K   
    #Epsilon_o calculation and storage
    Eps_o =np.dot(M_M_sec, Sigma_VEC0)
    epsilon_L.append(Eps_o[0])
    Sigma_L.append(Sigma_VEC0[0]) 
    #Sigma_o  increment
    Sigma_VEC0[0]=Sigma_VEC0[0]+0.1
    

    

    
      
#%% Plots
   
fig, ax = plt.subplots()
ax.plot(epsilon_L,Sigma_L,'+')

ax.grid(True)
ax.set_ylabel('SIGMA (GPa)')
ax.set_xlabel('EPSILON (GPa)')