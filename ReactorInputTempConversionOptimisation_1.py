# =========================================================================================================== #
# - Author :     Piotr T. Zaniewicz                                                                           #
# - Date   :     02/03/2023                                                                                   #
#                                                                                                             #                    
# - Description: - This script is used to optimise the inlet temperature of a reactor to maximise the         #
#                  conversion of N2 to NH3.                                                                   #
#                - The script is used to generate the data for the ammonia synthesis reactor section of       #
#                  the design project report.                                                                 #
# =========================================================================================================== #
#                                                                                                             #
# --------------------------------------   I N S T R U C T I O N S   ---------------------------------------- #
#                                                                                                             #
# - Input temp should not exceed 755K as this will cause the equilibrium constant to be negative              #
#   and the results are not valid for negative equilibrium constants                                          #
#                                                                                                             #
# - Value for BedLengthcalc is the bed length of the reactor in the design project report                     #
#                                                                                                             #
# - You should adjust thus value to match the bed length of the reactor you are using                         #
# which was determined in the ReactorConversionAcrossBed_1.py script                                          #
#                                                                                                             #
# - Running the script will calculate the optimal inlet temperature for the reactor                           #
#   to maximise the conversion of N2 to NH3.                                                                  #
#                                                                                                             #
# - The script will also plot the final conversion at each inlet temperature in the range specified           #
#   in the PlotTempRange variable.                                                                            # 
#                                                                                                             #
# - A second plot will be generated showing the final mole composition of the reactor at each inlet           #
#   temperature in the range specified in the PlotTempRange variable.                                         #
#                                                                                                             #
# - The script will also print the:     - optimal inlet temperature                                           #       
#                                       - resulting final conversion                                          #
#                                       - & resultant end mole composition                                    #
#   to the console.                                                                                           #
#                                                                                                             #
# - The value of 3.15m is the bed length of bed 1 in the design project report                                #
# =========================================================================================================== #
# =========================================================================================================== #
# ==================================   I N P U T   V A R I A B L E S   ====================================== #
# ----------------------------------------- Initial Conditions ---------------------------------------------- #
upperTempLimit = 803.15                 # Max K, cannot exceed (catalyst max temp)
constantPressure = 225                  # atm

# ---------------------------------------- Initial Mole Fractions ------------------------------------------- #
initialMoleFractionH2 = 0.714089
initialMoleFractionN2 = 0.238253
initialMoleFractionNH3 = 0.0213228
initialMoleFractionAr = 0.0262431

# -------------------------------------------- Temp Range Plot -----------------------------------------------#
PlotTempRange = [550, 740]              # [start, end] K (MUST BE INTEGER VALUE)
BedLengthcalc = 2.1                   # m
StepSize_dL = 0.001                     # m

# =========================================================================================================== #
# =========================================================================================================== #
# ===================================   I M P O R T   L I B R A R I E S   =================================== #
import numpy as np
from reactorCalcs_1 import ReactorUpdates
import matplotlib.pyplot as plt
import os
import pathlib
storagePath = os.path.join(pathlib.Path(__file__).parent.absolute(), "Figures")
import matplotlib.backends.backend_pdf
# =========================================================================================================== #

# =============================================   C L A S S E S   =========================================== #
# ===================================   R E A C T O R   S E Q U E N C E   =================================== #
class Reactor(ReactorUpdates):
    def __init__(
        self,
        stepSize,
        incomingTemp,
        pressure,
        bedLength,
        initialMoleFractionH2,
        initialMoleFractionN2,
        initialMoleFractionNH3,
        initialMoleFractionAr,
        Fn2,
        F
    ):
        super().__init__(
            stepSize,
            incomingTemp,
            pressure,
            bedLength,
            initialMoleFractionH2,
            initialMoleFractionN2,
            initialMoleFractionNH3,
            initialMoleFractionAr,
            248.153,
            1041.55,
        )

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


# =========================================   M A I N   P R O G R A M   ======================================#
def main():
    temperatureloop = []
    tempnow = PlotTempRange[0]
    temperatureloop.append(tempnow)

    outTempLimit = []
    outTempLimit.append(upperTempLimit)

    reactorconversioniterative = []
    reactorconversioniterative.append(0)

    reactortempfinaliterative = []
    reactortempfinaliterative.append(PlotTempRange[0])

    EquilibriumConstantiterative = []
    EquilibriumConstantiterative.append(0)

    reactionrateconstantiterative = []
    reactionrateconstantiterative.append(0)

    iterationsToMAX = 0

    yH2 = []
    yN2 = []
    yNH3 = []
    yAr = []
    yN2.append(initialMoleFractionN2)
    yH2.append(initialMoleFractionH2)
    yNH3.append(initialMoleFractionNH3)
    yAr.append(initialMoleFractionAr)
    CatMaxTemp = []
    CatMaxTemp.append(upperTempLimit)
    maxConversionRate = 0
    iterationsToMAX = 0 
    equilibriumConversionIterative = []



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    while tempnow < PlotTempRange[1]:

        R1 = Reactor(
            StepSize_dL,
            tempnow,
            constantPressure,
            BedLengthcalc,
            initialMoleFractionH2,
            initialMoleFractionN2,
            initialMoleFractionNH3,
            initialMoleFractionAr,
            248.153,
            1041.55,
        )
        R1.run(int(BedLengthcalc / StepSize_dL))

        reactorconversioniterative.append(R1.conversionN2)

        temperatureloop.append(tempnow)
        tempnow = tempnow + 1
        reactortempfinaliterative.append(R1.temp)
        outTempLimit.append(upperTempLimit)

        yN2.append(R1.moleFractionN2)
        yH2.append(R1.moleFractionH2)
        yNH3.append(R1.moleFractionNH3)
        yAr.append(R1.moleFractionAr)

        EquilibriumConstantiterative.append(R1.equilibriumConstant)
        reactionrateconstantiterative.append(R1.reactionRateConstant)

        equilibriumConversionIterative.append(R1.equilibriumConversion)



        CatMaxTemp.append(upperTempLimit)

        if R1.temp <= upperTempLimit:
            if reactorconversioniterative[-1] > maxConversionRate:
                maxConversionRate = reactorconversioniterative[-1]
                iterationsToMAX = tempnow - PlotTempRange[0] - 1
        else:
            pass

        Max_R1_Conversion = max(reactorconversioniterative)
        Optimal_Inlet_Temperature = temperatureloop[
            reactorconversioniterative.index(max(reactorconversioniterative))
        ]
        Resultant_Outlet_Temperature = reactortempfinaliterative[
            reactorconversioniterative.index(max(reactorconversioniterative))
        ]

        Max_R1_Conversion_Below_TempLim = reactorconversioniterative[iterationsToMAX]
        Optimal_Inlet_Temperature_Below_TempLim = temperatureloop[
            temperatureloop.index(PlotTempRange[0] + iterationsToMAX)
        ]
        Resultant_Outlet_Temperature_Below_TempLim = reactortempfinaliterative[
            temperatureloop.index(Optimal_Inlet_Temperature_Below_TempLim)
        ]


# --------------------- Print Results ---------------------#
    print(
        "\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    print(
        "Max R-601 Conversion:",
        reactorconversioniterative[
            reactorconversioniterative.index(max(reactorconversioniterative))
        ],
    )
    print(
        "Optimal Inlet Temperature: ",
        temperatureloop[
            reactorconversioniterative.index(max(reactorconversioniterative))
        ],
        "K",
    )
    print(
        "Resultant Outlet Temperature: ",
        reactortempfinaliterative[
            reactorconversioniterative.index(max(reactorconversioniterative))
        ],
        "K",
    )
    print(
        "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    print(
        "Max R-601 Conversion: below T_out Limit ",
        reactorconversioniterative[
            temperatureloop.index(PlotTempRange[0] + iterationsToMAX)
        ],
    )
    print(
        "Optimal Inlet Temperature: ",
        temperatureloop[temperatureloop.index(PlotTempRange[0] -1 + iterationsToMAX)],
        "K",
    )
    print(
        "Resultant Outlet Temperature: ",
        reactortempfinaliterative[
            temperatureloop.index(PlotTempRange[0] -1 + iterationsToMAX)
        ],
        "K",
    )
    print(
        "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )

    print(
        "Resultant Mole Fractions: ",
        "    N2: ",
        yN2[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        "    H2: ",
        yH2[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        "    NH3: ",
        yNH3[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        "   Ar: ",
        yAr[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
    )
    print(
        "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
    )

    # =====================    P L O T   R E S U L T S    =====================#
    # ---------------------- plot Tin, Tout & Conversion ----------------------#
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot()
    lns0 = ax.plot(
        temperatureloop,
        reactortempfinaliterative,
        color="blue",
        label="Final Temperature",
    )
    lns1 = ax.plot(
        iterationsToMAX + PlotTempRange[0],
        Resultant_Outlet_Temperature_Below_TempLim,
        marker="o",
        color="none",
        label=" Optimum Reaction Conditions",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    lns2 = ax.plot(
        temperatureloop,
        temperatureloop,
        color="black",
        linestyle="--",
        label="Tin = Tout",
    )
    lns3 = ax.plot(
        iterationsToMAX + PlotTempRange[0],
        Optimal_Inlet_Temperature_Below_TempLim,
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    plt.xlabel(
        "Temperature In (K)", color="black", rotation=0, fontsize=11, fontweight="light"
    )
    plt.ylabel(
        "Temperature Out (K)",
        color="blue",
        rotation=90,
        labelpad=15,
        fontsize=11,
        fontweight="light",
    )
    lns6 = ax.plot(
        temperatureloop,
        CatMaxTemp,
        color="red",
        linestyle="dotted",
        label="Catalyst Max Temperature",
    )
    plt.ylim(
        PlotTempRange[0],
    )
    plt.xlim(PlotTempRange[0], PlotTempRange[1])
    ax1 = ax.twinx()
    lns4 = ax1.plot(
        temperatureloop, reactorconversioniterative, color="green", label="Conversion"
    )
    lns7 = ax1.plot(
        temperatureloop[1:], np.array(equilibriumConversionIterative)-np.array(initialMoleFractionNH3)*np.array(initialMoleFractionAr), color="orange", label="Equilibrium Conversion"
    )


    lns5 = ax1.plot(
        iterationsToMAX + PlotTempRange[0],
        Max_R1_Conversion_Below_TempLim,
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    plt.ylabel(
        "Conversion",
        color="green",
        rotation=270,
        labelpad=15,
        fontsize=11,
        fontweight="light",
    )
    LEGEND = lns0 + lns1 + lns2 + lns3 + lns4 + lns5 + lns6 +lns7
    labs = [l.get_label() for l in LEGEND]
    ax.legend(LEGEND, labs, loc="upper left", bbox_to_anchor=(0.0, 0.78), fontsize=7)
    plt.ylim(0, 1)
    ax1.set_yticks(np.arange(0, 1.01, 0.1))
    plt.title("R-601 Temperature(in/out) vs Conversion")


    # --------------------- Plot Mole Fractions --------------------#
    fig1, ax = plt.subplots(figsize=(8, 4))
    plt.plot(temperatureloop, yN2, color="blue", label="N2")
    plt.plot(temperatureloop, yH2, color="red", label="H2")
    plt.plot(temperatureloop, yNH3, color="green", label="NH3")
    plt.plot(temperatureloop, yAr, color="black", label="Ar")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Mole Fraction")
    plt.xlim(PlotTempRange[0], PlotTempRange[1])
    plt.ylim(0, 1)
    plt.title("R-601 Mole Fractions vs Temperature")
    plt.plot(
        iterationsToMAX + PlotTempRange[0],
        yN2[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    plt.plot(
        iterationsToMAX + PlotTempRange[0],
        yH2[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    plt.plot(
        iterationsToMAX + PlotTempRange[0],
        yNH3[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
    )
    plt.plot(
        iterationsToMAX + PlotTempRange[0],
        yAr[temperatureloop.index(PlotTempRange[0] + iterationsToMAX)],
        label="Optimum Reaction Conditions",
        marker="o",
        markersize=7,
        markeredgewidth=2,
        markeredgecolor="red",
        markerfacecolor="None",
        color="none",
    )
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97),fancybox=True, shadow=False, ncol=5)


    # --------------------- Plot Equilibrium Constant --------------------#
    fig2, ax3 = plt.subplots(figsize=(8, 4))
    plt.plot((range(PlotTempRange[0],PlotTempRange[1])), reactionrateconstantiterative[1:], color="blue", label="Reaction Rate Constant")
    plt.plot(iterationsToMAX + PlotTempRange[0], reactionrateconstantiterative[iterationsToMAX], color="None", marker="o", markersize=7, markeredgewidth=2, markeredgecolor="red", markerfacecolor="None", label="Optimum Reaction Conditions")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Reaction Rate Constant")
    plt.xlim(PlotTempRange[0], PlotTempRange[1])
    plt.ylim(0, )
    # ax3.twinx()      # create a second y-axis
    # ax3.plot((range(PlotTempRange[0],PlotTempRange[1])), EquilibriumConstantiterative[1:], color="red", label="Equilibrium Constant")
    # plt.ylabel("Reaction Rate Constant")
    plt.legend()
    plt.title("R-601 Equilibrium Constant vs Temperature")

# --------------------- Save PDF ---------------------#

    pp = matplotlib.backends.backend_pdf.PdfPages(
            os.path.join(storagePath, "REACTOR_1__TEMP_VS_CONVERSION_OPTIMISATION.pdf")
    )
    #group all figures together in a list
    figs = [fig, fig1, fig2]
    for figs in figs:
        figs.set_size_inches(10.0, 5)
        figs.gca().grid(True, linestyle=':')
        # figs.legend(loc='upper center', bbox_to_anchor=(0.5, 0.9), shadow=True, ncol=4)
        pp.savefig(figs, bbox_inches="tight", dpi=300)
    pp.close()


# --------------------- Show Plots? --------------------#

    showFig = input("Show figures? (y/n): ")
    if showFig == "y":
        plt.show()  # show the figure
    else:
        quit()
        



# =====================   R U N   P R O G R A M    =====================#
    plt.show()

if __name__ == "__main__":
    main()
# =====================   E N D   O F   P R O G R A M    =====================#
