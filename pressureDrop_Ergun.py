# Purpose: Calculate the pressure drop across a packed bed

# Author:  Piotr T. Zaniewicz
# Date:    08/03/2023

# Imported libraries
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable

# =========================================================================================================== #
# - Description: - This script is used to calculate the pressure drop across a packed bed.                    #
#                  by using the Ergun equation.                                                               #
# =========================================================================================================== #
# Input data
dp = 0.010                                              # m         # diameter of catalyst particles        
d = 0.55                                                # m         # diameter of the bed                   

rho_1 = 34.60914526                                     # kg/m3     # density fluid                         # FROM ASPEN (avg. feed and product)
rho_2 = 36.38491739                                     # kg/m3     # density fluid                         # FROM ASPEN (avg. feed and product)
mu_1 = 0.0000291                                        # Pa.s      # viscosity fluid                       # FROM ASPEN (avg. feed and product)
mu_2 = 0.0000299                                        # Pa.s      # viscosity fluid                       # FROM ASPEN (avg. feed and product)

feed_1_flowrate = 4558.503765                           # L/min     # flowrate of feed 1
feed_1_flowrate_m3 = feed_1_flowrate/1000/60            # m3/s      # flowrate of feed 1
vs_1 = feed_1_flowrate_m3 / ((np.pi / 4) * d**2)        # m/s       # superficial velocity


feed_2_flowrate = 4414.515726                           # L/min     # flowrate of feed 2
feed_2_flowrate_m3 = feed_2_flowrate/1000/60            # m3/s      # flowrate of feed 2
vs_2 = feed_2_flowrate_m3 / ((np.pi / 4) * d**2)        # m/s       # superficial velocity
L_1 = 0                                                 # m         # length of the bed 1
L1 = 2.10                                               # m         # length of the bed 1
L_2 = 0                                                 # m         # length of the bed 2
L2 = 5.35                                               # m         # length of the bed 2
voidage = 0.40                                          # N/A       # void fraction 

# empty lists
lengthBed1 = []
lengthBed2 = []
pressureDropBed1 = []
pressureDropBed1.append(0)
pressureDropBed2 = []
pressureDropBed2.append(0)
dL_step = 0.01                                          # m         # length step size
Len_Bed1 = (dL_step, 2.10, dL_step)
Len_Bed2 = (dL_step, 5.35, dL_step)
dPdL_1__laminarPerLength = []                           # Pa/m      # Laminar flow term pressure drop per unit length in bed 1
dPdL_1__turbulentPerLength = []                         # Pa/m      # Turbulent flow term pressure drop per unit length in bed 1
dPdL_1__TotalPerLength = []                             # Pa/m      # Total pressure drop per unit length in bed 1
dPdL_1__laminarTotal = []                               # Pa        # Total pressure drop bed 1 from laminar flow term
dPdL_1__turbulentTotal = []                             # Pa        # Total pressure drop bed 1 from turbulent flow term
dPdL_1__Total = []                                      # Pa        # Total pressure drop bed 1
dPdL_2__laminarPerLength = []                           # Pa/m      # Laminar flow term pressure drop per unit length in bed 2
dPdL_2__turbulentPerLength = []                         # Pa/m      # Turbulent flow term pressure drop per unit length in bed 2
dPdL_2__TotalPerLength = []                             # Pa/m      # Total pressure drop per unit length in bed 2
dPdL_2__laminarTotal = []                               # Pa        # Total pressure drop bed 2 from laminar flow term
dPdL_2__turbulentTotal = []                             # Pa        # Total pressure drop bed 2 from turbulent flow term
dPdL_2__Total = []                                      # Pa        # Total pressure drop bed 2




# Calculations

# Create a loop to calculate the pressure drop across bed 1
while L1 >= L_1+dL_step:
    L_1 = L_1 + dL_step
    lengthBed1.append(L_1)

# Ergun equation - BED 1
    dPdL_1__laminarPerLength.append((150 * ( ((1 - voidage)**2) / (voidage**3) ) * ( (mu_1 * vs_1) / dp**2) ))            # Pa/m          # Laminar flow term pressure drop per unit length in bed 1
    dPdL_1__turbulentPerLength.append((1.75 * (rho_1 * vs_1**2 / dp) * (((1 - voidage)**2) / (voidage**3)) ))             # Pa/m          # Turbulent flow term pressure drop per unit length in bed 1
    dPdL_1__TotalPerLength.append(dPdL_1__laminarPerLength[-1] + dPdL_1__turbulentPerLength[-1])                                  # Pa/m          # Total pressure drop per unit length in bed 1

    dPdL_1__laminarTotal.append(dPdL_1__laminarPerLength[-1] * L_1)                                                           # Pa            # Total pressure drop bed 1 from laminar flow term
    dPdL_1__turbulentTotal.append(dPdL_1__turbulentPerLength[-1] * L_1)                                                       # Pa            # Total pressure drop bed 1 from turbulent flow term
    dPdL_1__Total.append((dPdL_1__laminarPerLength[-1] + dPdL_1__turbulentPerLength[-1]) * L_1)                                                                 # Pa            # Total pressure drop bed 1

    pressureDropBed1.append(dPdL_1__Total)

# Create a loop to calculate the pressure drop across bed 1
while L2 >= L_2+dL_step:
    L_2 = L_2 + dL_step
    lengthBed2.append(L_2)

# Ergun equation - BED 2
    dPdL_2__laminarPerLength.append((150 * ( ((1 - voidage)**2) / (voidage**3) ) * ( (mu_2 * vs_2) / dp**2) ))            # Pa/m          # Laminar flow term pressure drop per unit length in bed 2
    dPdL_2__turbulentPerLength.append((1.75 * (rho_2 * vs_2**2 / dp) * (((1 - voidage)**2) / (voidage**3)) ))             # Pa/m          # Turbulent flow term pressure drop per unit length in bed 2
    dPdL_2__TotalPerLength.append(dPdL_2__laminarPerLength[-1] + dPdL_2__turbulentPerLength[-1])                                  # Pa/m          # Total pressure drop per unit length in bed 2

    dPdL_2__laminarTotal.append(dPdL_2__laminarPerLength[-1] * L_2)                                                           # Pa            # Total pressure drop bed 2 from laminar flow term
    dPdL_2__turbulentTotal.append(dPdL_2__turbulentPerLength[-1] * L_2)                                                       # Pa            # Total pressure drop bed 2 from turbulent flow term
    dPdL_2__Total.append((dPdL_2__laminarPerLength[-1] + dPdL_2__turbulentPerLength[-1])*L_2)                                                               # Pa            # Total pressure drop bed 2
    
    pressureDropBed2.append(dPdL_2__Total)



# Print results
print('==================================================================================================================================')
print('==================================================================================================================================')

print('The (per m length) pressure drop in bed 1 corresponding to the:')
print('                                                                 laminar flow term       =  ',round(dPdL_1__laminarPerLength[-1], ),' Pa/m', '    or    ', round(dPdL_1__laminarPerLength[-1]/1e5, 3),'bar/m')
print('                                                                 turbulent flow term     =  ',round(dPdL_1__turbulentPerLength[-1], ),'Pa/m', '   or    ', round(dPdL_1__turbulentPerLength[-1]/1e5, 3),'bar/m')
print('                                                                 TOTAL                   =  ',round(dPdL_1__TotalPerLength[-1],),'Pa/m', '   or    ', round(dPdL_1__TotalPerLength[-1]/1e5, 3),'bar/m')

print('----------------------------------------------------------------------------------------------------------------------------------')
print('The (total) pressure drop in bed 1 corresponding to the :')
print('                                                                 laminar flow term       =  ',round(dPdL_1__laminarTotal[-1], ),' Pa', '     or    ', round(dPdL_1__laminarTotal[-1]/1e5, 3),'bar')
print('                                                                 turbulent flow term     =  ',round(dPdL_1__turbulentTotal[-1], ),'Pa', '     or    ', round(dPdL_1__turbulentTotal[-1]/1e5, 3),'bar')
print('                                                                 TOTAL                   =  ',round(dPdL_1__Total[-1], ),'Pa', '     or    ', round(dPdL_1__Total[-1]/1e5, 3),'bar')

if dPdL_1__laminarPerLength > dPdL_1__turbulentPerLength:
    dPdL_1 = dPdL_1__laminarPerLength
    print('_____ THEREFORE THE FLOW IS LAMINAR _____')
else:
    dPdL_1 = dPdL_1__turbulentPerLength
    print('_____ THEREFORE THE FLOW IS TURBULENT _____')

print('==================================================================================================================================')

print('The (per m length) pressure drop in bed 2 corresponding to the:')
print('                                                                 laminar flow term       =  ',round(dPdL_2__laminarPerLength[-1],),'  Pa/m', '   or    ', round(dPdL_2__laminarPerLength[-1]/1e5, 3),'bar/m')
print('                                                                 turbulent flow term     =  ',round(dPdL_2__turbulentPerLength[-1],),'Pa/m', '   or    ', round(dPdL_2__turbulentPerLength[-1]/1e5, 3),'bar/m')
print('                                                                 TOTAL                   =  ',round(dPdL_2__TotalPerLength[-1], ),'Pa/m', '   or    ', round(dPdL_2__TotalPerLength[-1]/1e5, 3),'bar/m')

print('----------------------------------------------------------------------------------------------------------------------------------')

print('The (total) pressure drop in bed 2 corresponding to the :')
print('                                                                 laminar flow term       =  ',round(dPdL_2__laminarTotal[-1], ),'  Pa', '    or    ', round(dPdL_2__laminarTotal[-1]/1e5, 3),'bar')
print('                                                                 turbulent flow term     =  ',round(dPdL_2__turbulentTotal[-1], ),'Pa', '    or    ', round(dPdL_2__turbulentTotal[-1]/1e5, 3),'bar')
print('                                                                 TOTAL                   =  ',round(dPdL_2__Total[-1], ),'Pa', '    or    ', round(dPdL_2__Total[-1]/1e5, 3),'bar')

if dPdL_2__laminarPerLength > dPdL_2__turbulentPerLength:
    dPdL_2 = dPdL_2__laminarPerLength
    print('_____ THEREFORE THE FLOW IS LAMINAR _____')
else:
    dPdL_2 = dPdL_2__turbulentPerLength
    print('_____ THEREFORE THE FLOW IS TURBULENT _____')
    
print('==================================================================================================================================')
print('==================================================================================================================================')


# Plot results
plt.figure(3)
plt.plot(dPdL_1__laminarTotal, 'b', label='Laminar flow term')
plt.plot(dPdL_1__turbulentTotal, 'r', label='Turbulent flow term')
plt.plot(dPdL_1__Total, 'k', label='Total pressure drop')
plt.xlabel('Bed length [m]')
plt.ylabel('Pressure drop [Pa]')
plt.title('Pressure drop in bed 1')
plt.ylim(0, )
plt.xlim(0, )
plt.legend()
plt.grid()

plt.figure(4)
plt.plot(dPdL_2__laminarTotal, 'b', label='Laminar flow term')
plt.plot(dPdL_2__turbulentTotal, 'r', label='Turbulent flow term')
plt.plot(dPdL_2__Total, 'k', label='Total pressure drop')
plt.xlabel('Bed length [m]')
plt.ylabel('Pressure drop [Pa]')
plt.title('Pressure drop in bed 2')
plt.ylim(0, )
plt.xlim(0, )
plt.legend()
plt.grid()



# plot table of printed results in console
table = PrettyTable()
table._set_double_border_style()
table.field_names = [" Pressure Drop Across Packed Bed " , "Bed 1", "Bed 2"]
table.add_row(["Laminar flow term (Pa)", round(dPdL_1__laminarTotal[-1], ), round(dPdL_2__laminarTotal[-1], )])
table.add_row(["Turbulent flow term (Pa)", round(dPdL_1__turbulentTotal[-1], ), round(dPdL_2__turbulentTotal[-1], )])
table.add_row(["Total (Pa)", round(dPdL_1__Total[-1], ), round(dPdL_2__Total[-1], )])

table2 = PrettyTable()
table2._set_double_border_style()
table2.field_names = [" Pressure Drop over Bed " , "Bed 1", "Bed 2"]
table2.add_row(["Laminar flow term (bar)", round(dPdL_1__laminarTotal[-1]/1e5, 3), round(dPdL_2__laminarTotal[-1]/1e5, 3)])
table2.add_row(["Turbulent flow term (bar)", round(dPdL_1__turbulentTotal[-1]/1e5, 3), round(dPdL_2__turbulentTotal[-1]/1e5, 3)])
table2.add_row(["Total (bar)", round(dPdL_1__Total[-1]/1e5, 3), round(dPdL_2__Total[-1]/1e5, 3)])

print(table,'\n-------------------------------------------------')
print(table2)

