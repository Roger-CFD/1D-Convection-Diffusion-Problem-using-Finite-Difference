
L=1.0
gamma=0.1


#case1 

u1=0.1
rho=1
n1=5
N1=n1+1
delta_x1=L/n1
pe_cell1=(u1*delta_x1)/gamma #cell peclet number
pe_L1=(u1*L)/gamma #domain peclet number  
print(f" --------\nINVESTIGATING THE PECLET NUMBER FOR CASE 1 \n -------- \nCell peclet number is {pe_cell1:.5g} \nDomain peclet number={pe_L1:.3g} \n")
if -2<=pe_cell1<=2 :
    print("FINITE Central difference scheme is stable")
else:
    print("FINITE Central difference will be unstable and produce ocillations, note that because of the iterative solving method will amplify the oscillations and will create overflow. \nConvection is dominant over diffusion for the whole domain since the peclet number is >>1 . \nBackward differencing should be used because it accounts for directionality of information transport \n")
#case 2 
u2=2.5
n2=5
N2=n2+1
delta_x2=L/n2
pe_cell2=(u2*delta_x2)/gamma #cell peclet number
pe_L2=(u2*L)/gamma #domain peclet number  
print(f" --------\nINVESTIGATING THE PECLET NUMBER FOR CASE 2 \n -------- \nCell peclet number is {pe_cell2:.5g} \nDomain peclet number={pe_L2:.3g} \n")
if -2<=pe_cell2<=2 :
    print("FINITE Central difference scheme is stable")
else:
    print("FINITE Central difference will be unstable and produce ocillations, note that because of the iterative solving method will amplify the oscillations and will create overflow. \nConvection is dominant over diffusion for the whole domain since the peclet number is >>1 . \nBackward differencing should be used because it accounts for directionality of information transport \n")
u3=2.5
n3=20
N3=n3+1
delta_x3=L/n3
pe_cell3=(u3*delta_x3)/gamma #cell peclet number
pe_L3=(u3*L)/gamma #domain peclet number  
print(f" --------\nINVESTIGATING THE PECLET NUMBER FOR CASE 3 \n -------- \nCell peclet number is {pe_cell3:.5g} \nDomain peclet number={pe_L3:.3g} \n")
if -2<=pe_cell3<=2 :
    print("FINITE Central difference scheme is stable")
else:
    print("FINITE Central difference will be unstable and produce ocillations, note that because of the iterative solving method will amplify the oscillations and will create overflow. \nConvection is dominant over diffusion for the whole domain since the peclet number is >>1 . \nBackward differencing should be used because it accounts for directionality of information transport \n")
