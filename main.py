# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 16:18:44 2026

@author: stan
"""
from param import *
from solver import*
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 

#case 1 
#initialise the temperature
T_0_1=np.zeros(N1)
T_0_1[N1-1]=1
#compute Temeperature numerically , calls the solver 
Temp1=FCDS(T_0_1,pe_cell1,N1)
"--------------------"
#case 2
#initialise the temperature
T_0_2=np.zeros(N2)
T_0_2[N2-1]=1
#compute Temeperature numerically , calls the solver 
Temp2=FBDS1(T_0_2,pe_cell2,N2)
Temp2_2order=FBDS2(T_0_2,pe_cell2,N2)
"--------------------"
#case 3 
#initialise the temperature
T_0_3=np.zeros(N3)
T_0_3[N3-1]=1
#compute Temeperature numerically , calls the solver 
Temp3=FCDS(T_0_3,pe_cell3,N3)

#ploting results
x_dom=np.arange(N1)*delta_x1
x = np.linspace(0, 1, 100)

#case 1 
T_Analytical1=(np.exp(rho*u1*x/gamma)-1)/(np.exp(rho*u1*L/gamma)-1) #ANALYTICAL temperature for case 1 
plt.plot(x_dom, Temp1, "ro-", label="Numerical")
plt.plot(x, T_Analytical1, "b-", label="Analytical")
plt.xlabel("x")
plt.ylabel("T")
plt.title("Numerical vs Analytical Comparison for CASE 1 ")
plt.legend()
plt.grid()
plt.show()
#Table 1 
T_Analytical_interp1 = np.interp(x_dom, x, T_Analytical1)
# Absolute difference at each node
Difference = np.abs(Temp1 - T_Analytical_interp1)

# Maximum absolute error
print("Table 1 for case 1 : Central difference scheme for u=0.1 and 5 mesh")
print("--------------------------------------------")
Max_error = np.max(Difference)
table_1=pd.DataFrame({
    "Node": np.arange(1,N1+1),
    "Distance" :x_dom,
    "Numerical solution" :Temp1,
    "Analytical solution":T_Analytical_interp1,
    "Absolute error":Difference,
    
    })
table_1_display=table_1.round(6)
print(table_1_display.to_string(index=False))
table_1_display.to_csv("table_1_display.csv", index=False)
print(f"\nMaximum error = {Max_error:.6e}\n")
print("--------------------------------------------")    
#case 2 
#computing analytical temperature distribution for case 2 
T_Analytical2=(np.exp(rho*u2*x/gamma)-1)/(np.exp(rho*u2*L/gamma)-1) #ANALYTICAL temperature 
plt.plot(x_dom, Temp2, "ro-", label="Numerical 1st order ")
plt.plot(x_dom, Temp2_2order, "s--", label="Numerical 2nd order")
plt.plot(x, T_Analytical2, "b-", label="Analytical")
plt.xlabel("x")
plt.ylabel("T")
plt.title("Numerical 1st order and 2nd order vs Analytical Comparison for case 2 using 1st order backward scheme")
plt.legend()
plt.grid()
plt.show()
T_Analytical_interp2 = np.interp(x_dom, x, T_Analytical2)
# Absolute difference at each node
Difference = np.abs(Temp2 - T_Analytical_interp2)

# Maximum absolute error
print("Table 2 for case 2 : Backward difference 1st order scheme for u=2.5 and 5 mesh")
print("--------------------------------------------")
Max_error1 = np.max(Difference)
table_2_1st_order=pd.DataFrame({
    "Node": np.arange(1,N2+1),
    "Distance" :x_dom,
    "Numerical solution 01" :Temp2,
    "Analytical solution":T_Analytical_interp2,
    "Absolute error":Difference,
    
    })
table_2_1st_order_display=table_2_1st_order.round(6)
print(table_2_1st_order_display.to_string(index=False))
table_2_1st_order_display.to_csv("table_2_1st_order_display.csv", index=False)
print(f"\nMaximum error = {Max_error1:.6e}\n")
print("Table 3 for case 2 : Backward difference 2nd order scheme for u=2.5 and 5 mesh")
print("--------------------------------------------")
#table for second order backward scheme    
# Absolute difference at each node
Difference = np.abs(Temp2_2order - T_Analytical_interp2)

# Maximum absolute error
Max_error2 = np.max(Difference)
table_2_2nd_order=pd.DataFrame({
    "Node": np.arange(1,N2+1),
    "Distance" :x_dom,
    "Numerical solution O2" :Temp2_2order,
    "Analytical solution":T_Analytical_interp2,
    "Absolute error":Difference,
    
    })
table_2_2nd_order_display=table_2_2nd_order.round(6)
print(table_2_2nd_order_display.to_string(index=False))
table_2_2nd_order_display.to_csv("table_2_2nd_order_display.csv", index=False)
print(f"\nMaximum error = {Max_error2:.6e}\n")  
#case 3 : mesh refinement 
x_dom3=np.arange(N3)*delta_x3
T_Analytical3=(np.exp(rho*u3*x/gamma)-1)/(np.exp(rho*u3*L/gamma)-1) #ANALYTICAL temperature 
plt.plot(x_dom3, Temp3, "ro-", label="Numerical")
plt.plot(x, T_Analytical3, "b-", label="Analytical")
plt.xlabel("x")
plt.ylabel("T")
plt.title("Numerical vs Analytical Comparison for case 3 using Central differencing scheme after mesh refinement")
plt.legend()
plt.grid()
plt.show()
T_Analytical_interp3 = np.interp(x_dom3, x, T_Analytical3)
# Absolute difference at each node
Difference = np.abs(Temp3 - T_Analytical_interp3)
print("Table 4 for case 3 : Central differencing scheme after mesh refinement for u=2.5 and 20 mesh")  
print("--------------------------------------------")
# Maximum absolute error
Max_error3 = np.max(Difference)
table_3=pd.DataFrame({
    "Node": np.arange(1,N3+1),
    "Distance" :x_dom3,
    "Numerical solution" :Temp3,
    "Analytical solution":T_Analytical_interp3,
    "Absolute error":Difference,
    
    })
table_3_display=table_3.round(6)
print(table_3_display.to_string(index=False))
table_3_display.to_csv("table_3_display.csv", index=False)
print(f"\nMaximum error = {Max_error3:.6e}")
print("--------------------------------------------")  
    
import os

table_3_display.to_csv("table_3_display.csv", index=False)

print(os.path.abspath("table_3_display.csv"))





