# Plot Heat Capacities of components of the system

import numpy as np
import matplotlib.pyplot as plt


# composition
initialMoleFractionH2 = 0.714089
initialMoleFractionN2 = 0.238253
initialMoleFractionNH3 = 0.0213228
initialMoleFractionAr = 0.0262431

# temperature range
pressure = 225  #atm
T1 = int(input("Enter a starting temperature (501 - 900K): "))
T2 = int(input("Enter an ending temperature (501 - 900K): "))
FeedFlowRate = float(input("Enter the feed flow rate (kg/s): "))
CpH2 = []
CpN2 = []
CpAr = []
CpNH3 = []
Cp = []
temperature = []

yH2 = initialMoleFractionH2
yN2 = initialMoleFractionN2
yNH3 = initialMoleFractionNH3
yAr = initialMoleFractionAr

molecularWeightH2 = 2.016
molecularWeightN2 = 28.0134
molecularWeightNH3 = 17.0305
molecularWeightAr = 39.948

molecularWeight = molecularWeightH2 * initialMoleFractionH2 + molecularWeightN2 * initialMoleFractionN2 + molecularWeightNH3 * initialMoleFractionNH3 + molecularWeightAr * initialMoleFractionAr # g/mol

# Calculate the heat capacity of each component of the system
temp = T1

while temp <= T2 :

    HeatCapacityH2 = (
        #                 # specific_heat_hydrogen * yH2
        (
            4.184
            * (
                6.952
                - 4.576e-4 * temp
                + 9.563e-7 * temp**2
                - 2.079e-10 * temp**3
            )
        ))


    HeatCapacityN2 = (
        #        # + specific_heat_nitrogen * yN2
        + (
            4.184
            * (
                6.903
                - 3.753e-4 * temp
                + 1.93e-6 * temp**2
                - 6.861e-10 * temp**3
            )
        ))


    HeatCapacityAr = (
        #        # + specific_heat_argon * yAr
            + (4.9675
                ) * 4.184
        )


    HeatCapacityNH3 = (
        # + specific_heat_ammonia * yNH3
        + (4.184 *
            (
                6.5846 * 1
                - 6.1251e-3 * temp
                + 2.3663e-6 * temp**2
                - 1.5981e-9 * temp**3
                + (
                    96.1678
                    - 0.067571 * pressure
                    + (-0.2225 + 1.6847e-4 * pressure) * temp
                    + (1.289e-4 - 1.0095e-7 * pressure) * temp**2
                )
            )
        )
    ) 

    temperature.append(temp)
    CpH2.append(HeatCapacityH2)
    CpN2.append(HeatCapacityN2)
    CpAr.append(HeatCapacityAr)
    CpNH3.append(HeatCapacityNH3)
    Cp.append(HeatCapacityH2 *initialMoleFractionH2 + HeatCapacityN2 *initialMoleFractionN2 + HeatCapacityAr *initialMoleFractionAr + HeatCapacityNH3 *initialMoleFractionNH3)
    temp += 1


AverageHeatCapacity = sum(Cp)/len(Cp)
print ("The average heat capacity of the system across the user specified temperature range is:   ", AverageHeatCapacity, "kJ/mol*K") 
print("The heat transfer rate per mole is:   ", AverageHeatCapacity * (T2-T1), "kJ/kmol") # on a per mole basis
print("The heat transfer rate per unit mass is:    ", AverageHeatCapacity * (T2-T1) / molecularWeight, "kJ/kg") # on a per unit mass basis
print("Q, The heat transfer rate is:    ", AverageHeatCapacity * (T2-T1) * FeedFlowRate / molecularWeight, "kJ/S") # on a per unit mass basis
print("The Duty of the system is:    ", (AverageHeatCapacity * (T2-T1) * FeedFlowRate / molecularWeight) /1e3, "MW") # on a per unit mass basis


#Plot the heat capacities of the components of the system
fig1, ax1 = plt.subplots()
plt.plot(temperature, CpH2, label="H2")
plt.plot(temperature, CpN2, label="N2")
plt.plot(temperature, CpAr, label="Ar")
plt.plot(temperature, CpNH3, label="NH3")
plt.plot(temperature, Cp, label="Average Heat Capacity")
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Heat Capacity (J/mol*K)")
plt.ylim(0, )
plt.xlim(T1,T2)
plt.title("Heat Capacity of cold stream components")
plt.show()
