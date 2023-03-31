# =========================================================================================================== #
# - Author :     Piotr T. Zaniewicz                                                                           #
# - Date   :     09/03/2023                                                                                   #
# - Description: - This script is used to model the conversion of N2 to NH3 across a packed bed reactor.      #
#                - The script is used to generate the data for the ammonia synthesis reactor section of       #
#                  the design project report.                                                                 #
#                - The function main() is used to run the reactor simulation.                                 #
#                - The function conversion_pressureR1() is used to generate the data for the conversion       #
#                  of N2 to NH3 across both beds at varying pressures as defined in list var. pressurelist()  #
# =========================================================================================================== #
# =====================================   I N F O R M A T I O N   =========================================== #
# - This program calculates the conversion of N2 to NH3 across a packed bed reactor.                          #
# - The program uses the following assumptions:                                                               #
#   - The reactor is isothermal                                                                               #
#   - The reactor is isobaric                                                                                 #
# =========================================================================================================== #
# ==================================   I N P U T   V A R I A B L E S   ====================================== #
# ----------------------------------------- Initial Conditions ---------------------------------------------- #
upperTempLimit = 803.15               # Max K, cannot exceed (catalyst max temp)                              #
constantPressure = 225                # atm (assumed constatn as pressure drop across reactor is negligible   #
#                                     # and hence it's effects are ignored in this model)                     #
pressurelist = [150, 175, 200, 225, 250]                                                                      #
# ------------------------------------------- Reactor Parameters -------------------------------------------- #
StepSize = 0.001                                   # m                                                        #
bed1 = 2.10                                                                                                   #
bed2 = 5.35                                                                                                   #
chosenLength1 = int(bed1 / StepSize)               # m (chosen length of reactor bed 1)                       #
chosenLength2 = int(bed2 / StepSize)               # m (chosen length of reactor bed 2)                       #
MW_H2 = 2.016                                                                                                 #
MW_N2 = 28.0134                                                                                               #
MW_NH3 = 17.0305                                                                                              #
MW_Ar = 39.948                                                                                                #     
# =========================================================================================================== #
# =========================================================================================================== #
# ===================================   I M P O R T   L I B R A R I E S   =================================== #
import numpy as np
from reactorCalcs_1 import Fn2, ReactorUpdates
import matplotlib.pyplot as plt
import os
import pathlib
storagePath = os.path.join(pathlib.Path(__file__).parent.absolute(), "Figures")
import matplotlib.backends.backend_pdf
import csv
import dataclasses
# =========================================================================================================== #

bed1conversion = 0.14296               # conversion of bed 1 -
bed1feed = 1041.55

R2FN2 = 248.153 * (1 - bed1conversion)
R2F = bed1feed - 2 * Fn2 * bed1conversion

@dataclasses.dataclass
class ReactorConfig:
    StepSize: float
    incomingTemp: float
    pressure: float
    BedLengthcalc: float
    initialMoleFractionH2: float
    initialMoleFractionN2: float
    initialMoleFractionNH3: float
    initialMoleFractionAr: float
    baseLength: float
    upperTempLimit: float = 803.15               # Max K, cannot exceed (catalyst max temp)
    constantPressure: float = 225               # atm (assumed constatn as pressure drop across reactor is negligible

    def __post_init__(self):
        self.chosenLengthIndex = int(self.baseLength / self.StepSize)               # distance along reactor bed locator


R1Config = ReactorConfig(StepSize=0.001,
                         incomingTemp=673.15,
                         pressure=225,
                         BedLengthcalc=5,
                         initialMoleFractionH2=0.714089,
                         initialMoleFractionN2=0.238253,
                         initialMoleFractionNH3=0.0213228,
                         initialMoleFractionAr=0.0262431,
                         baseLength=2.10)

R2Config = ReactorConfig(StepSize=0.001,
                         incomingTemp=692,
                         pressure=225,
                         BedLengthcalc=10,
                         initialMoleFractionH2=0.656696,
                         initialMoleFractionN2=0.219138,
                         initialMoleFractionNH3=0.0959075,
                         initialMoleFractionAr=0.0281596,
                         baseLength=5.35)


# =========================================================================================================== #
# =============================================   C L A S S E S   =========================================== #
# ===================================   R E A C T O R   S E Q U E N C E   =================================== #
class Reactor(ReactorUpdates):

    def __init__(self, stepSize, incomingTemp, pressure, bedLength, R1InitialMoleFractionH2, R1InitialMoleFractionN2,
                 R1InitialMoleFractionNH3, R1InitialMoleFractionAr, Fn2, F):
        super().__init__(stepSize, incomingTemp, pressure, bedLength, R1InitialMoleFractionH2, R1InitialMoleFractionN2,
                         R1InitialMoleFractionNH3, R1InitialMoleFractionAr, Fn2, F)

    def run(self, iterations=1):
        for _ in range(iterations):
            self.updateEffFactor()
            self.updateFugacityN2()
            self.updateFugacityH2()
            self.updateFugacityNH3()
            self.updateHeatOfReaction()
            self.updateSpecificHeat()
            self.updateMoleFractionN2()
            self.updateMoleFractionH2()
            self.updateMoleFractionNH3()
            self.updateMoleFractionAr()
            self.updateActivationCoefficientN2()
            self.updateActivationCoefficientH2()
            self.updateActivationCoefficientNH3()
            self.updateReactionRateConstant()
            self.updateEquilibriumConstant()
            self.updateEquilibriumConversion()
            self.updateRateOfReactionNH3()
            self.updateConversionN2()
            self.updateTemp()
            self.updateStep()


#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #


# =====================   T E M P   P R E S S U R E   P R O G R A M M E   =====================#
def conversion_pressureR1():

    fig21, ax00 = plt.subplots(figsize=(10, 5))
    ax01 = ax00.twinx()
    for i in range(0, 5):
        pressure = pressurelist[i]

        R1 = Reactor(
            R1Config.StepSize,
            R1Config.incomingTemp,
            pressure,
            R1Config.BedLengthcalc,
            R1Config.initialMoleFractionH2,
            R1Config.initialMoleFractionN2,
            R1Config.initialMoleFractionNH3,
            R1Config.initialMoleFractionAr,
            248.153,
            1041.55,
        )
        R1.run(int(R1Config.BedLengthcalc / StepSize))

        ax00.plot(R1.steps, R1._temp, "--", color="black")
        ax01.plot(R1.steps, R1._conversionN2, label=str(pressurelist[i]) + "atm")
        ax00.plot(len([x for x in R1._temp if x < upperTempLimit]) * StepSize,
                  upperTempLimit,
                  "ro",
                  markersize=10,
                  markerfacecolor="None",
                  markeredgecolor="red",
                  markeredgewidth=3)

    plt.xlim(0, R1Config.BedLengthcalc)
    ax00.set_ylim(R1Config.incomingTemp,)
    ax00.set_xlim(0, R1Config.BedLengthcalc)
    ax00.set_xlabel("Length of Bed (m)")
    ax00.set_ylabel("Temperature (K)")
    ax01.set_ylabel("Conversion", color="black", rotation=270, labelpad=20)
    ax01.set_ylim(0,)
    plt.title("R-601 Temperature/Conversion Profile at varying pressures")
    plt.legend()
    #plt.savefig(os.path.join(storagePath, "TempProfilePressure__R-601.png"))

    fig22, ax02 = plt.subplots(figsize=(10, 5))
    ax03 = ax02.twinx()
    for i in range(0, 5):
        pressure = pressurelist[i]
        R2 = Reactor(
            R2Config.StepSize, 
            R2Config.incomingTemp, 
            pressure, 
            R2Config.BedLengthcalc,
            R2Config.initialMoleFractionH2,
            R2Config.initialMoleFractionN2, 
            R2Config.initialMoleFractionNH3,
            R2Config.initialMoleFractionAr, 
            R2FN2, 
            R2F
        )
        R2.run(int(R2Config.BedLengthcalc / StepSize))

        ax02.plot(R2.steps, R2._temp, "--", color="black")
        ax03.plot(R2.steps, R2._conversionN2, label=str(pressurelist[i]) + "atm")
        ax02.plot(len([x for x in R2._temp if x < upperTempLimit]) * StepSize,
                  upperTempLimit,
                  "ro",
                  markersize=10,
                  markerfacecolor="none",
                  markeredgecolor="r",
                  markeredgewidth=3)

    plt.xlim(0, R2Config.BedLengthcalc)
    ax02.set_ylim(R2Config.incomingTemp,)
    ax02.set_xlim(0, R2Config.BedLengthcalc)
    ax02.set_xlabel("Length of Bed (m)")
    ax02.set_ylabel("Temperature (K)")
    ax03.set_ylabel("Conversion", color="black", rotation=270, labelpad=20)
    ax03.set_ylim(0,)
    plt.title("R-602 Temperature/Conversion Profile at varying pressures")
    plt.legend()
    plt.savefig(os.path.join(storagePath, "TempProfilePressure__R-602.png"))

    pp = matplotlib.backends.backend_pdf.PdfPages(os.path.join(storagePath, "PressureVsConversion.pdf"))
    #group all figures together in a list
    figs = [fig21, fig22]
    for figs in figs:
        figs.set_size_inches(9.0, 5)
        figs.gca().grid(True, linestyle=':')
        pp.savefig(figs, bbox_inches="tight", dpi=500)
    pp.close()

    showFig = input("Show Temperature vs. Pressure figures? (y/n): ")
    if showFig == "y":
        plt.show()               # show the figure
    else:
        pass


#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #


# =====================   M A I N    P R O G R A M   =====================#
def main():
    R1 = Reactor(
        R1Config.StepSize,
        R1Config.incomingTemp,
        R1Config.constantPressure,
        R1Config.BedLengthcalc,
        R1Config.initialMoleFractionH2,
        R1Config.initialMoleFractionN2,
        R1Config.initialMoleFractionNH3,
        R1Config.initialMoleFractionAr,
        248.153,
        1041.55,
    )
    R1.run(int(R1Config.BedLengthcalc / R1Config.StepSize))

    R2 = Reactor(
        R2Config.StepSize,
        R2Config.incomingTemp,
        R2Config.constantPressure,
        R2Config.BedLengthcalc,
        R2Config.initialMoleFractionH2,
        R2Config.initialMoleFractionN2,
        R2Config.initialMoleFractionNH3,
        R2Config.initialMoleFractionAr,
        R2FN2,
        R2F
    )
    R2.run(int(R2Config.BedLengthcalc / R2Config.StepSize))

    # ----------------------- print various results to 3 d.p --------------------------------#
    print("R-601 Starting Temperature: ", round(R1.incomingTemp, 3), "K")
    print("R-601 Final Temperature: ", round(R1.temp, 3), "K")
    print("---------------------------------------------------------------------")
    print("R-601 Initial Mole Fraction H2: ", round(R1Config.initialMoleFractionH2, 3))
    print("R-601 Initial Mole Fraction N2: ", round(R1Config.initialMoleFractionN2, 3))
    print("R-601 Initial Mole Fraction NH3: ", round(R1Config.initialMoleFractionNH3, 3))
    print("R-601 Initial Mole Fraction Ar: ", round(R1Config.initialMoleFractionAr, 3))
    print("---------------------------------------------------------------------")
    print("R-601 Final Mole Fraction H2: ", round(R1.moleFractionH2, 5))
    print("R-601 Final Mole Fraction N2: ", round(R1.moleFractionN2, 5))
    print("R-601 Final Mole Fraction NH3: ", round(R1.moleFractionNH3, 5))
    print("R-601 Final Mole Fraction Ar: ", round(R1.moleFractionAr, 5))
    print("---------------------------------------------------------------------")
    print("R-601 EQUILIBRIUM CONSTANT: ", round(R1.equilibriumConstant, 3), "kmol^-1")
    print("R-601 REACTION RATE CONSTANT: ", round(R1.reactionRateConstant, 3), "kmol hr^-1")
    print("---------------------------------------------------------------------")
    print("R-601 FUGACITY H2: ", round(R1.fugacityH2, 3))
    print("R-601 FUGACITY N2: ", round(R1.fugacityN2, 3))
    print("R-601 FUGACITY NH3: ", round(R1.fugacityNH3, 3))
    print("---------------------------------------------------------------------")
    print("R-601 HEAT OF REACTION: ", round(R1.heatOfReaction), "kJ kmol^-1 NH3")
    print("R-601 CONVERSION N2: ", round(R1.conversionN2) * 100, "%")
    print("---------------------------------------------------------------------")
    print("R-601 NActivation Coefficient N2: ", round(R1.activationCoefficientN2, 3))
    print("R-601 Activation Coefficient H2: ", round(R1.activationCoefficientH2, 3))
    print("R-601 Activation Coefficient NH3: ", round(R1.activationCoefficientNH3, 3))
    print("---------------------------------------------------------------------")
    print("R-601 SPECIFIC HEAT: ", round(R1.specificHeat, 3), "kJ kmol^-1 K^-1")
    print("R-601 Max Conversion: ", round(R1._conversionN2[-1] * 100, 3), "%")
    print(
        "R-601 Average Molecular Weight: ",
        round(
            R1.moleFractionAr * MW_Ar + R1.moleFractionH2 * MW_H2 + R1.moleFractionN2 * MW_N2 +
            R1.moleFractionNH3 * MW_NH3, 3), "g mol^-1")
    print(
        "R-601 Average Iniial Molecular Weight: ",
        round(
            R1.initialMoleFractionAr * MW_Ar + R1.initialMoleFractionH2 * MW_H2 + R1.initialMoleFractionN2 * MW_N2 +
            R1.initialMoleFractionNH3 * MW_NH3, 3), "g mol^-1")
    print("---------------------------------------------------------------------")
    print("R-602 Starting Temperature: ", round(R2.incomingTemp, 3), "K")
    print("R-602 Final Temperature: ", round(R2.temp, 3), "K")
    print("---------------------------------------------------------------------")
    print("R-602 Initial Mole Fraction H2: ", round(R2.initialMoleFractionH2, 3))
    print("R-602 Initial Mole Fraction N2: ", round(R2.initialMoleFractionN2, 3))
    print("R-602 Initial Mole Fraction NH3: ", round(R2.initialMoleFractionNH3, 3))
    print("R-602 Initial Mole Fraction Ar: ", round(R2.initialMoleFractionAr, 3))
    print("---------------------------------------------------------------------")
    print("R-602 Final Mole Fraction H2: ", round(R2.moleFractionH2, 5))
    print("R-602 Final Mole Fraction N2: ", round(R2.moleFractionN2, 5))
    print("R-602 Final Mole Fraction NH3: ", round(R2.moleFractionNH3, 5))
    print("R-602 Final Mole Fraction Ar: ", round(R2.moleFractionAr, 5))
    print("---------------------------------------------------------------------")
    print("R-602 EQUILIBRIUM CONSTANT: ", round(R2.equilibriumConstant, 3), "kmol^-1")
    print("R-602 REACTION RATE CONSTANT: ", round(R2.reactionRateConstant, 3), "kmol hr^-1")
    print("---------------------------------------------------------------------")
    print("R-602 FUGACITY H2: ", round(R2.fugacityH2, 3))
    print("R-602 FUGACITY N2: ", round(R2.fugacityN2, 3))
    print("R-602 FUGACITY NH3: ", round(R2.fugacityNH3, 3))
    print("---------------------------------------------------------------------")
    print("R-602 HEAT OF REACTION: ", round(R2.heatOfReaction), "kJ kmol^-1 NH3")
    print("R-602 CONVERSION N2: ", round(R2.conversionN2) * 100, "%")
    print("---------------------------------------------------------------------")
    print("R-602 Activation Coefficient N2: ", round(R2.activationCoefficientN2, 3))
    print("R-602 Activation Coefficient H2: ", round(R2.activationCoefficientH2, 3))
    print("R-602 Activation Coefficient NH3: ", round(R2.activationCoefficientNH3, 3))
    print("---------------------------------------------------------------------")
    print("R-602 SPECIFIC HEAT: ", round(R2.specificHeat, 3), "kJ kmol^-1 K^-1")
    print("R-602 Max Conversion: ", round(R2._conversionN2[-1] * 100, 5), "%")
    print(
        "R-602 Average Molecular Weight: ",
        round(
            R2.moleFractionAr * MW_Ar + R2.moleFractionH2 * MW_H2 + R2.moleFractionN2 * MW_N2 +
            R2.moleFractionNH3 * MW_NH3, 3), "g mol^-1")
    print(
        "R-602 Average Iniial Molecular Weight: ",
        round(
            R2.initialMoleFractionAr * MW_Ar + R2.initialMoleFractionH2 * MW_H2 + R2.initialMoleFractionN2 * MW_N2 +
            R2.initialMoleFractionNH3 * MW_NH3, 3), "g mol^-1")
    print("---------------------------------------------------------------------")

    # save results to csv file
    with open(
            'AmmoniaReactormodel.csv',
            'w',
            newline='',
    ) as file:
        #save local variables to csv file
        writer = csv.writer(file)
        writer.writerow([
            "Temperature (K)", "conversionN2", "Mole Fraction N2", "Mole Fraction H2", "Mole Fraction NH3",
            "Mole Fraction Ar", "Activation Coefficient N2", "Activation Coefficient H2", "Activation Coefficient NH3",
            "Reaction Rate Constant", "Equilibrium Constant", "Rate of Reaction NH3", "Conversion N2",
            "Temperature (K)", "Steps", " ", "Temperature (K)", "conversionN2", "Mole Fraction N2", "Mole Fraction H2",
            "Mole Fraction NH3", "Mole Fraction Ar", "Activation Coefficient N2", "Activation Coefficient H2",
            "Activation Coefficient NH3", "Reaction Rate Constant", "Equilibrium Constant", "Rate of Reaction NH3",
            "Conversion N2", "Temperature (K)", "Steps"
        ])
        writer.writerows(
            zip(R1._temp, R1._conversionN2, R1._moleFractionN2, R1._moleFractionH2, R1._moleFractionNH3,
                R1._moleFractionAr, R1._activationCoefficientN2, R1._activationCoefficientH2,
                R1._activationCoefficientNH3, R1._reactionRateConstant, R1._equilibriumConstant, R1._rateOfReactionNH3,
                R1._conversionN2, R1._temp, R1._steps, " ", R2._temp, R2._conversionN2, R2._moleFractionN2,
                R2._moleFractionH2, R2._moleFractionNH3, R2._moleFractionAr, R2._activationCoefficientN2,
                R2._activationCoefficientH2, R2._activationCoefficientNH3, R2._reactionRateConstant,
                R2._equilibriumConstant, R2._rateOfReactionNH3, R2._conversionN2, R2._temp, R2._steps))

#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #

# =====================    P L O T   R E S U L T S  R - 6 0 1  =====================#
# --------------------- plot temperature and conversion ---------------------#
    temp_sequence = R1._temp
    conversion_sequence = R1._conversionN2
    # fig1, ax = plt.subplots(figsize=(8, 4))
    fig1, ax = plt.subplots()
    lns1 = ax.plot(temp_sequence, color="blue", label="Temperature (K)", linewidth=4, linestyle="dotted")
    ax.tick_params(axis="y", labelcolor="black")
    ax.set_ylim(R1._temp[0],)
    ax2 = ax.twinx()
    lns2 = ax2.plot(conversion_sequence, color="green", label="Conversion X", linewidth=2, linestyle="-")
    ax2.tick_params(axis="y", labelcolor="black")
    ax2.set_ylim(0,)
    ax.set_xlim(0, len(temp_sequence))
    ax.plot(chosenLength1, R1._temp[chosenLength1], "ro", markersize=7, label="Chosen Length")
    plt.xlabel("Length (m)")
    plt.ylabel("Temperature (K)", color="blue")
    ax2.set_xlim(0, len(conversion_sequence))
    plt.ylabel("Conversion", color="green", rotation=270, labelpad=15)
    LEGEND = lns1 + lns2
    labs = [l.get_label() for l in LEGEND]
    #relable x axis values
    ax.set_adjustable('datalim')
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    ax.set_xticklabels([0, 1, 2, 3, 4, 5])
    ax2.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    ax2.set_xticklabels([0, 1, 2, 3, 4, 5])
    ax.get_xgridlines()[0].set_visible(True)
    ax2.get_xgridlines()[0].set_visible(True)
    ax.get_xgridlines()[1].set_visible(True)
    ax2.get_xgridlines()[1].set_visible(True)
    ax.get_xgridlines()[2].set_visible(True)
    ax2.get_xgridlines()[2].set_visible(True)
    ax.get_xgridlines()[3].set_visible(True)
    ax2.get_xgridlines()[3].set_visible(True)
    ax.get_xgridlines()[4].set_visible(True)
    ax2.get_xgridlines()[4].set_visible(True)
    ax.get_xgridlines()[5].set_visible(True)
    ax2.get_xgridlines()[5].set_visible(True)
    ax.set_xlabel("Length (m)")
    ax2.set_xlabel("Length (m)")
    ax.set_ylabel("Temperature (K)", color="blue")
    ax2.set_ylabel("Conversion", color="green", rotation=270, labelpad=15)


    # --------------------- plot rate of reaction ---------------------#
    fig2, ax = plt.subplots()
    plt.plot(R1.steps[1:], R1._rateOfReactionNH3[1:], color="blue", label="Rate of reaction ")
    plt.xlabel("Length (m)")
    plt.ylabel(r"Rate of reaction (kmol$_{H_2}$ m$^{-3}$ hr$^{-1}$)")
    ax.plot(chosenLength1 * StepSize, R1._rateOfReactionNH3[chosenLength1], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot rate of formation NH3 ---------------------#

    fig10, ax = plt.subplots()
    plt.plot(R1.steps[1:], R1._rateOfReactionNH3[1:], color="blue", label="Rate of formation")
    plt.xlabel("Length (m)")
    # raise text to power
    plt.ylabel(r"Rate of formation (kmol$_{NH_3}$ m$^{-3}$ hr$^{-1}$)")
    ax.plot(chosenLength1 * StepSize, R1._rateOfReactionNH3[chosenLength1], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot mole fraction of each species ---------------------#
    fig4, ax = plt.subplots()
    plt.plot(R1.steps, R1._moleFractionH2, color="blue", label="yH2")
    plt.plot(R1.steps, R1._moleFractionN2, color="green", label="yN2")
    plt.plot(R1.steps, R1._moleFractionNH3, color="red", label="yNH3")
    plt.plot(R1.steps, R1._moleFractionAr, color="orange", label="yAr")
    ax.plot(chosenLength1 * StepSize, R1._moleFractionH2[chosenLength1], "ro", markersize=7)
    ax.plot(chosenLength1 * StepSize, R1._moleFractionN2[chosenLength1], "ro", markersize=7)
    ax.plot(chosenLength1 * StepSize, R1._moleFractionNH3[chosenLength1], "ro", markersize=7)
    ax.plot(chosenLength1 * StepSize, R1._moleFractionAr[chosenLength1], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("mole fraction")
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot equilibrium constant vs distance ---------------------#
    fig3, ax = plt.subplots()
    plt.plot(R1._steps[1:], R1._equilibriumConstant[1:], color="blue", label="Equilibrium constant")
    plt.xlabel("Length (m)")
    plt.ylabel("Equilibrium constant")
    ax.plot(chosenLength1 * StepSize, R1._equilibriumConstant[chosenLength1], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot reaction rate constant vs distance ---------------------#
    fig9, ax = plt.subplots()
    plt.plot(R1._steps[1:], R1._reactionRateConstant[1:], color="blue", label="Reaction Rate Constant")
    ax.plot(chosenLength1 * StepSize, R1._reactionRateConstant[chosenLength1], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("Reaction Rate Constant")
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot efficiency factor ---------------------#
    fig5, ax = plt.subplots()
    plt.plot(R1._steps[1:], R1._effFactor[1:], color="blue", label="Efficiency Factor")
    ax.plot(chosenLength1 * StepSize, R1._effFactor[chosenLength1], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("Efficiency Factor")
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)

    # --------------------- plot reaction rate constant vs T ---------------------#
    fig6, ax = plt.subplots()
    plt.plot(R1._temp[1:], R1._reactionRateConstant[1:], color="blue", label="Reaction Rate Constant")
    ax.plot(R1._temp[chosenLength1], R1._reactionRateConstant[chosenLength1], "ro", markersize=7)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Reaction Rate Constant")
    ax.set_ylim(0,)
    ax.set_xlim(temp_sequence[0], temp_sequence[-1])

    #  --------------------- plot equilibrium constant vs T ---------------------#
    fig8, ax = plt.subplots()
    plt.plot(R1._temp[1:], R1._equilibriumConstant[1:], color="blue", label="Equilibrium constant")
    ax.plot(R1._temp[chosenLength1], R1._equilibriumConstant[chosenLength1], "ro", markersize=7)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Equilibrium constant")
    ax.set_ylim(0,)
    ax.set_xlim(temp_sequence[0], temp_sequence[-1])

    #  plot activation coefficients ---------------------#
    fig7, ax = plt.subplots()
    plt.plot(
        R1.steps,
        R1._activationCoefficientN2,
        color="blue",
        label="Activation Coefficient N2",
    )
    plt.plot(
        R1.steps,
        R1._activationCoefficientH2,
        color="green",
        label="Activation Coefficient H2",
    )
    plt.plot(
        R1.steps,
        R1._activationCoefficientNH3,
        color="red",
        label="Activation Coefficient NH3",
    )
    plt.xlabel("Length (m)")
    plt.ylabel("Activation Coefficient")
    ax.set_ylim(0,)
    ax.set_xlim(0, R1.bedLength)
    ax.plot(chosenLength1 * StepSize, R1._activationCoefficientN2[chosenLength1], "ro", markersize=7)
    ax.plot(chosenLength1 * StepSize, R1._activationCoefficientH2[chosenLength1], "ro", markersize=7)
    ax.plot(chosenLength1 * StepSize, R1._activationCoefficientNH3[chosenLength1], "ro", markersize=7)

    # --------------------- Save PDF ---------------------#
    pp = matplotlib.backends.backend_pdf.PdfPages(os.path.join(storagePath, "R-601_ONLY_ALL.pdf"))
    #group all figures together in a list
    figs = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10]
    for figs in figs:
        figs.set_size_inches(9.0, 5)
        figs.gca().grid(True, linestyle=':')
        figs.legend(loc='upper center', bbox_to_anchor=(0.5, 0.98), shadow=True, ncol=4)
        pp.savefig(figs, bbox_inches="tight", dpi=300)
    pp.close()

    # #  =========================================================================================================================================================================== #
    # #  =========================================================================================================================================================================== #

    # =====================    P L O T   R E S U L T S  R - 6 0 2  =====================#
    # --------------------- plot temperature and conversion ---------------------#
    temp_sequence2 = R2._temp
    conversion_sequence2 = R2._conversionN2
    fig11, ax = plt.subplots()
    lns1 = ax.plot(temp_sequence2, color="blue", label="Temperature (K)", linewidth=4, linestyle="dotted")
    ax.tick_params(axis="y", labelcolor="black")
    ax.set_ylim(R2._temp[0],)
    ax2 = ax.twinx()
    lns2 = ax2.plot(conversion_sequence2, color="green", label="Conversion X", linewidth=2, linestyle="-")
    ax2.tick_params(axis="y", labelcolor="black")
    ax2.set_ylim(0,)
    ax.set_xlim(0, len(temp_sequence2))
    ax.plot(chosenLength2, R2._temp[chosenLength2], "ro", markersize=7, label="Chosen Length")
    plt.xlabel("Length (m)")
    plt.ylabel("Temperature (K)", color= "blue")
    ax2.set_xlim(0, len(conversion_sequence))
    plt.ylabel("Conversion", color="green",rotation=270, labelpad=15)
    LEGEND = lns1 + lns2
    labs = [l.get_label() for l in LEGEND]
    #relable x axis values
    ax.set_adjustable('datalim')
    ax.set_xticklabels([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
    ax2.set_xticklabels([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax2.set_xticks([0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
    ax.get_xgridlines()[0].set_visible(True)
    ax2.get_xgridlines()[0].set_visible(True)
    ax.get_xgridlines()[1].set_visible(True)
    ax2.get_xgridlines()[1].set_visible(True)
    ax.get_xgridlines()[2].set_visible(True)
    ax2.get_xgridlines()[2].set_visible(True)
    ax.get_xgridlines()[3].set_visible(True)
    ax2.get_xgridlines()[3].set_visible(True)
    ax.get_xgridlines()[4].set_visible(True)
    ax2.get_xgridlines()[4].set_visible(True)
    ax.get_xgridlines()[5].set_visible(True)
    ax2.get_xgridlines()[5].set_visible(True)
    ax.get_xgridlines()[6].set_visible(True)
    ax2.get_xgridlines()[6].set_visible(True)
    ax.get_xgridlines()[7].set_visible(True)
    ax2.get_xgridlines()[7].set_visible(True)
    ax.get_xgridlines()[8].set_visible(True)
    ax2.get_xgridlines()[8].set_visible(True)
    ax.get_xgridlines()[9].set_visible(True)
    ax2.get_xgridlines()[9].set_visible(True)
    ax.get_xgridlines()[10].set_visible(True)
    ax2.get_xgridlines()[10].set_visible(True)
    ax.set_xlabel("Length (m)")
    ax2.set_xlabel("Length (m)")
    ax.set_ylabel("Temperature (K)", color="blue")
    ax2.set_ylabel("Conversion", color="green", rotation=270, labelpad=15)

    # --------------------- plot rate of reaction ---------------------#
    fig12, ax = plt.subplots()
    plt.plot(R2.steps[1:], R2._rateOfReactionNH3[1:], color="blue", label="Rate of reaction ")
    plt.xlabel("Length (m)")
    plt.ylabel(r"Rate of reaction (kmol$_{H_2}$ m$^{-3}$ hr$^{-1}$)")
    ax.plot(chosenLength2 * StepSize, R2._rateOfReactionNH3[chosenLength2], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot rate of formation NH3 ---------------------#

    fig13, ax = plt.subplots()
    plt.plot(R2.steps[1:], R2._rateOfReactionNH3[1:], color="blue", label="Rate of formation")
    plt.xlabel("Length (m)")
    # raise text to power
    plt.ylabel(r"Rate of formation (kmol$_{NH_3}$ m$^{-3}$ hr$^{-1}$)")
    ax.plot(chosenLength2 * StepSize, R2._rateOfReactionNH3[chosenLength2], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot mole fraction of each species ---------------------#
    fig14, ax = plt.subplots()
    plt.plot(R2.steps, R2._moleFractionH2, color="blue", label="yH2")
    plt.plot(R2.steps, R2._moleFractionN2, color="green", label="yN2")
    plt.plot(R2.steps, R2._moleFractionNH3, color="red", label="yNH3")
    plt.plot(R2.steps, R2._moleFractionAr, color="orange", label="yAr")
    ax.plot(chosenLength2 * StepSize, R2._moleFractionH2[chosenLength2], "ro", markersize=7)
    ax.plot(chosenLength2 * StepSize, R2._moleFractionN2[chosenLength2], "ro", markersize=7)
    ax.plot(chosenLength2 * StepSize, R2._moleFractionNH3[chosenLength2], "ro", markersize=7)
    ax.plot(chosenLength2 * StepSize, R2._moleFractionAr[chosenLength2], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("mole fraction")
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot equilibrium constant vs distance ---------------------#
    fig15, ax = plt.subplots()
    plt.plot(R2._steps[1:], R2._equilibriumConstant[1:], color="blue", label="Equilibrium constant")
    plt.xlabel("Length (m)")
    plt.ylabel("Equilibrium constant")
    ax.plot(chosenLength2 * StepSize, R2._equilibriumConstant[chosenLength2], "ro", markersize=7)
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot reaction rate constant vs distance ---------------------#
    fig16, ax = plt.subplots()
    plt.plot(R2._steps[1:], R2._reactionRateConstant[1:], color="blue", label="Reaction Rate Constant")
    ax.plot(chosenLength2 * StepSize, R2._reactionRateConstant[chosenLength2], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("Reaction Rate Constant")
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot efficiency factor ---------------------#
    fig17, ax = plt.subplots()
    plt.plot(R2._steps[1:], R2._effFactor[1:], color="blue", label="Efficiency Factor")
    ax.plot(chosenLength2 * StepSize, R2._effFactor[chosenLength2], "ro", markersize=7)
    plt.xlabel("Length (m)")
    plt.ylabel("Efficiency Factor")
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)

    # --------------------- plot reaction rate constant vs T ---------------------#
    fig18, ax = plt.subplots()
    plt.plot(R2._temp[1:], R2._reactionRateConstant[1:], color="blue", label="Reaction Rate Constant")
    ax.plot(R2._temp[chosenLength2], R2._reactionRateConstant[chosenLength2], "ro", markersize=7)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Reaction Rate Constant")
    ax.set_ylim(0,)
    ax.set_xlim(temp_sequence2[0], temp_sequence2[-1])

    #  --------------------- plot equilibrium constant vs T ---------------------#
    fig19, ax = plt.subplots()
    plt.plot(R2._temp[1:], R2._equilibriumConstant[1:], color="blue", label="Equilibrium constant")
    ax.plot(R2._temp[chosenLength2], R2._equilibriumConstant[chosenLength2], "ro", markersize=7)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Equilibrium constant")
    ax.set_ylim(0,)
    ax.set_xlim(temp_sequence2[0], temp_sequence2[-1])

    #  plot activation coefficients ---------------------#
    fig20, ax = plt.subplots()
    plt.plot(
        R2.steps,
        R2._activationCoefficientN2,
        color="blue",
        label="Activation Coefficient N2",
    )
    plt.plot(
        R2.steps,
        R2._activationCoefficientH2,
        color="green",
        label="Activation Coefficient H2",
    )
    plt.plot(
        R2.steps,
        R2._activationCoefficientNH3,
        color="red",
        label="Activation Coefficient NH3",
    )
    plt.xlabel("Length (m)")
    plt.ylabel("Activation Coefficient")
    ax.set_ylim(0,)
    ax.set_xlim(0, R2.bedLength)
    ax.plot(chosenLength2 * R2Config.StepSize, R2._activationCoefficientN2[chosenLength2], "ro", markersize=7)
    ax.plot(chosenLength2 * R2Config.StepSize, R2._activationCoefficientH2[chosenLength2], "ro", markersize=7)
    ax.plot(chosenLength2 * R2Config.StepSize, R2._activationCoefficientNH3[chosenLength2], "ro", markersize=7)

    # --------------------- Save PDF ---------------------#
    pp = matplotlib.backends.backend_pdf.PdfPages(os.path.join(storagePath, "R-602_ONLY_ALL.pdf"))
    #group all figures together in a list
    figs = [fig11, fig12, fig13, fig14, fig15, fig16, fig17, fig18, fig19, fig20]
    for figs in figs:
        figs.set_size_inches(9.0, 5)
        figs.gca().grid(True, linestyle=':')
        figs.legend(loc='upper center', bbox_to_anchor=(0.5, 0.98), shadow=True, ncol=4)
        pp.savefig(figs, bbox_inches="tight", dpi=300)
    pp.close()

    # --------------------- Show Plots? --------------------#
    showFig = input("Show individual figures? (y/n): ")
    if showFig == "y":
        plt.show()               # show the figure
    else:
        matplotlib.pyplot.close("all")               # close all figures

#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #
#  =========================================================================================================================================================================== #

#  C O M B I N E D   P L O T S

# --------------------- plot temperature and conversion ---------------------#
    fig23, ax4 = plt.subplots()
    ax4.plot(list(np.arange(0, R1Config.chosenLengthIndex, R1Config.StepSize))[:R1Config.chosenLengthIndex],
             R1._temp[:R1Config.chosenLengthIndex],
             color="blue",
             label="Temperature (K)",
             linewidth=4,
             linestyle="dotted")
    ax4.plot(list(R1Config.baseLength +
                  np.arange(0, R2Config.chosenLengthIndex, R2Config.StepSize))[:R2Config.chosenLengthIndex],
             R2._temp[:R2Config.chosenLengthIndex],
             color="blue",
             linewidth=4,
             linestyle="dotted")
    ax4.tick_params(axis="y", labelcolor="blue")
    ax4.set_ylim(R1._temp[0],)
    plt.ylabel("Temperature (K)", color="blue")
    ax5 = ax4.twinx()
    ax5.plot(list(np.arange(0, R1Config.chosenLengthIndex, R1Config.StepSize))[:R1Config.chosenLengthIndex],
             R1._conversionN2[:R1Config.chosenLengthIndex],
             color="green",
             label="Conversion X",
             linewidth=2,
             linestyle="-")
    ax5.plot(list(R1Config.baseLength +
                  np.arange(0, R2Config.chosenLengthIndex, R2Config.StepSize))[:R2Config.chosenLengthIndex],
             R2._conversionN2[:R2Config.chosenLengthIndex],
             color="green",
             linewidth=2,
             linestyle="-")
    ax5.tick_params(axis="y", labelcolor="black")
    ax5.set_ylim(0,)
    ax5.axvline(x=bed1, color='black', linestyle='--')
    plt.xlabel("Length (m)")
    ax5.set_xlim(0, R1Config.baseLength + R2Config.baseLength)
    plt.ylabel("Conversion", color="green", rotation=270, labelpad=15)

    #relable x axis values
    # ax.set_adjustable('datalim')
    # ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    # ax.set_xticklabels([0, 1, 2, 3, 4, 5])
    # ax2.set_xticks([0, 1000, 2000, 3000, 4000, 5000])
    # ax2.set_xticklabels([0, 1, 2, 3, 4, 5])
    # ax.get_xgridlines()[0].set_visible(True)
    # ax2.get_xgridlines()[0].set_visible(True)
    # ax.get_xgridlines()[1].set_visible(True)
    # ax2.get_xgridlines()[1].set_visible(True)
    # ax.get_xgridlines()[2].set_visible(True)
    # ax2.get_xgridlines()[2].set_visible(True)
    # ax.get_xgridlines()[3].set_visible(True)
    # ax2.get_xgridlines()[3].set_visible(True)
    # ax.get_xgridlines()[4].set_visible(True)
    # ax2.get_xgridlines()[4].set_visible(True)
    # ax.get_xgridlines()[5].set_visible(True)
    # ax2.get_xgridlines()[5].set_visible(True)
    # ax.set_xlabel("Length (m)")
    # ax2.set_xlabel("Length (m)")
    # ax.set_ylabel("Temperature (K)")
    # ax2.set_ylabel("Conversion")

    # --------------------- plot rate of reaction ---------------------#
    fig24, ax6 = plt.subplots()
    plt.plot(R1.steps[1:R1Config.chosenLengthIndex],
             R1._rateOfReactionNH3[1:R1Config.chosenLengthIndex],
             color="blue",
             label="Rate of reaction ")
    plt.plot([item + R1.steps[R1Config.chosenLengthIndex] for item in R2.steps[1:R2Config.chosenLengthIndex]],
             R2._rateOfReactionNH3[1:R2Config.chosenLengthIndex],
             color="blue")
    plt.xlabel("Length (m)")
    plt.ylabel(r"Rate of reaction (kmol$_{H_2}$ m$^{-3}$ hr$^{-1}$)")
    ax6.axvline(x=bed1, color='black', linestyle='--')

    ax6.set_ylim(0,)
    ax6.set_xlim(0, R1Config.baseLength + R2Config.baseLength)

    # --------------------- plot mole fraction of each species ---------------------#
    fig25, ax7 = plt.subplots()
    plt.plot(R1.steps[1:R1Config.chosenLengthIndex],
             R1._moleFractionH2[1:R1Config.chosenLengthIndex],
             color="blue",
             label="yH2")
    plt.plot(R1.steps[1:R1Config.chosenLengthIndex],
             R1._moleFractionN2[1:R1Config.chosenLengthIndex],
             color="green",
             label="yN2")
    plt.plot(R1.steps[1:R1Config.chosenLengthIndex],
             R1._moleFractionNH3[1:R1Config.chosenLengthIndex],
             color="red",
             label="yNH3")
    plt.plot(R1.steps[1:R1Config.chosenLengthIndex], R2._moleFractionAr[1:R1Config.chosenLengthIndex], color="orange", label="yAr")
    ax7.axvline(x=bed1, color='black', linestyle='--')
    plt.plot([item + R1.steps[R1Config.chosenLengthIndex] for item in R2.steps[1:R2Config.chosenLengthIndex]],
             R2._moleFractionH2[1:R2Config.chosenLengthIndex],
             color="blue")
    plt.plot([item + R1.steps[R1Config.chosenLengthIndex] for item in R2.steps[1:R2Config.chosenLengthIndex]],
             R2._moleFractionN2[1:R2Config.chosenLengthIndex],
             color="green")
    plt.plot([item + R1.steps[R1Config.chosenLengthIndex] for item in R2.steps[1:R2Config.chosenLengthIndex]],
             R2._moleFractionNH3[1:R2Config.chosenLengthIndex],
             color="red")
    plt.plot([item + R1.steps[R1Config.chosenLengthIndex] for item in R2.steps[1:R2Config.chosenLengthIndex]],
             R2._moleFractionAr[1:R2Config.chosenLengthIndex],
             color="orange")
    plt.xlim(0, R1Config.baseLength + R2Config.baseLength)
    plt.xlabel("Length (m)")
    plt.ylabel("mole fraction")
    ax7.set_ylim(0,)

    figs = [fig23, fig24, fig25]
    for figs in figs:
        figs.set_size_inches(9.0, 5)
        figs.gca().grid(True, linestyle=':')
        #add title left of dashed line
        figs.text(0.21, 0.83, "R-601", fontsize=15, transform=figs.transFigure, fontweight='bold')
        #add title right of dashed line
        figs.text(0.6, 0.83, "R-602", fontsize=15, transform=figs.transFigure, fontweight='bold')
        figs.legend(loc='upper center', bbox_to_anchor=(0.5, 0.98), shadow=True, ncol=4)


# --------------------- Save PDF ---------------------#
    pp = matplotlib.backends.backend_pdf.PdfPages(os.path.join(storagePath, "COMBINED_CONSECUTIVE_PLOTS.pdf"))
    #group all figures together in a list
    figs = [fig23, fig24, fig25]
    for figs in figs:
        figs.set_size_inches(9.0, 5)
        figs.gca().grid(True, linestyle=':')
        figs.legend(loc='upper center', bbox_to_anchor=(0.5, 0.98), shadow=True, ncol=4)
        pp.savefig(figs, bbox_inches="tight", dpi=300)
    pp.close()

    # --------------------- Show Plots? --------------------#
    showFig = input("Show Combined figures? (y/n): ")
    if showFig == "y":
        plt.show()               # show the figures
    else:
        pass

    # =====================   R U N   P R O G R A M    =====================#

if __name__ == "__main__":
    main()
    conversion_pressureR1()
# =====================   E N D   O F   P R O G R A M    =====================#
