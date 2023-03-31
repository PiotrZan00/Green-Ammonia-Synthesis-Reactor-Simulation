# =========================================================================================================== #
# - Author :     Piotr T. Zaniewicz                                                                           #
# - Date   :     02/03/2023                                                                                   #
#                                                                                                             #
# - Description: This script contains the calculations for the reactor bed 1                                  #
# =========================================================================================================== #
#                                                                                                             #
# ===================================   I M P O R T   L I B R A R I E S   =================================== #
import numpy as np
from reactorUtils import ReactorBase
# =========================================================================================================== #
# ------------------------------------------ C O N S T A N T S ---------------------------------------------- #
ko = 8.85e14        # Arrhenius constant              -                                                       #
E = 1.7056e5        # activation energy               - Activation energy for ammonia synthesis with          # 
#                                                       catalyst(J/mol)                                       #
R = 8.314           # Universal Gas Constant:         - R = 8.314 J/mol-K                                     #
alpha = 0.5         # Temkin parameter:               - can range from: 0.5 - 0.75                            #
#                                                       (0.5 is most common and is used in this calculation)  #
# ==================================   I N P U T   V A R I A B L E S   ====================================== #
diameter_internal = 0.55 # internal diameter of packed bed - m                                                #
A = np.pi * (diameter_internal / 2) ** 2            # cross-sectional area of packed bed    - m^2                        #
print(A)
F = 1041.55        # total feed molar flowrate              - kmol/hr                                          #
Fn2 = 248.153         # total feed flowrate for nitrogen    - kmol/hr                                         #
stepSize = 0.001    # increment size for each step           - 1/iterations                                   #
#                                                                                                             #
# =========================================================================================================== #
# ===================================   I M P O R T   L I B R A R I E S   =================================== #
import numpy as np
from reactorUtils import ReactorBase
# =========================================================================================================== #

# =============================================   C L A S S E S   =========================================== #

class FugacityCalcs(ReactorBase):
    def calcFugacityN2(self):
        fugacity = (
            0.93431737
            + (0.2028538e-3) * self.temp
            + (0.295896e-3) * self.pressure
            - (0.270727e-6) * (self.temp**2)
            + (0.4775207e-6) * self.pressure**2
        )
        return fugacity

    def calcFugacityH2(self):
        fugacity = np.exp(
            np.exp(-3.8402 * (self.temp**0.125) + 0.541) * self.pressure
            - np.exp(-0.1263 * (self.temp**0.5) - 15.98) * self.pressure**2
            + 300 * np.exp(-0.011901 * self.temp - 5.941) * np.exp(-self.pressure / 300)
        )
        return fugacity

    def calcFugacityNH3(self):
        fugacity = (
            0.1438996
            + 0.002028538 * self.temp
            - (0.4487672e-3) * self.pressure
            - (0.1142945e-5) * (self.temp**2)
            + (0.2761216e-6) * self.pressure**2
        )
        return fugacity


class MoleFractionCalcs(ReactorBase):
    def calcMoleFractionN2(self):
        moleFraction = ((self.initialMoleFractionN2 * self.F) - (self.Fn2 * self.conversionN2)) / (
            self.F - (2 * self.conversionN2 * self.Fn2)
        )
        return moleFraction

    def calcMoleFractionH2(self):
        moleFraction = (
            ((self.initialMoleFractionH2 * self.F)
            - (3  * self.conversionN2 * self.initialMoleFractionN2 * self.F))
        ) / ( self.F - (2 * self.conversionN2 * self.Fn2))
        return moleFraction

    def calcMoleFractionNH3(self):
        moleFraction = (
            (self.initialMoleFractionNH3 * self.F)
            + (2 * self.conversionN2 * self.Fn2)
        ) / ( self.F - (2 * self.conversionN2 * self.Fn2))
        return moleFraction

    def calcMoleFractionAr(self):
        moleFraction = (
            (self.F * self.initialMoleFractionAr)
            + (0)
        ) / (self.F - (2 * self.conversionN2 * self.Fn2))
        return moleFraction


class ActivationCoefficientCalcs(ReactorBase):
    def calcActivationCoefficientH2(self):
        return self.fugacityH2 * self.moleFractionH2 * self.pressure

    def calcActivationCoefficientN2(self):
        return self.fugacityN2 * self.moleFractionN2 * self.pressure

    def calcActivationCoefficientNH3(self):
        return self.fugacityNH3 * self.moleFractionNH3 * self.pressure

    def calcActivationCoefficientAr(self):
        pass


class ReactorCalcs(FugacityCalcs, MoleFractionCalcs, ActivationCoefficientCalcs):
    def calcEffFactor(self):
        effFactorCoeff = [
            -8.2125534,
            0.03774149,
            6.190112,
            -5.354571e-5,
            -20.86963,
            2.379142e-8,
            27.88403,
        ]

        effFactor = (
            effFactorCoeff[0]
            + effFactorCoeff[1] * self.temp
            + effFactorCoeff[2] * self.conversionN2
            + effFactorCoeff[3] * (self.temp**2)
            + effFactorCoeff[4] * (self.conversionN2**2)
            + effFactorCoeff[5] * (self.temp**3)
            + effFactorCoeff[6] * (self.conversionN2**3)
        )
        return effFactor

    def calcHeatOfReaction(self):   # J/mol_NH3
        heatOfReaction = 4.184 * (
            -self.pressure
            * (0.54526 + (340.609 / self.temp) + (459.734 * 10**6) / (self.temp**3))
            - 5.34685 * self.temp
            - 0.0002525 * (self.temp**2)
            + 0.00000169167 * (self.temp**3)
            - 9157.09
        )
        return heatOfReaction


    def calcSpecificHeat(self):         
        #                     UNITS: kcal/kmol/K
        specificHeat = (
            #                 # specific_heat_hydrogen * yH2
            (
                4.184
                * (
                    6.952
                    - 4.576e-4 * self.temp
                    + 9.563e-7 * self.temp**2
                    - 2.079e-10 * self.temp**3
                )
                * self.moleFractionH2
            )
            #        # + specific_heat_nitrogen * yN2
            + (
                4.184
                * (
                    6.903
                    - 3.753e-4 * self.temp
                    + 1.93e-6 * self.temp**2
                    - 6.861e-10 * self.temp**3
                )
                * self.moleFractionN2
            )
            #        # + specific_heat_argon * yAr
                + (4.9675
                   * self.moleFractionAr) * 4.184
            # + specific_heat_ammonia * yNH3
            + (4.184 *
                (
                    6.5846 * 1
                    - 6.1251e-3 * self.temp
                    + 2.3663e-6 * self.temp**2
                    - 1.5981e-9 * self.temp**3
                    + (
                        96.1678
                        - 0.067571 * self.pressure
                        + (-0.2225 + 1.6847e-4 * self.pressure) * self.temp
                        + (1.289e-4 - 1.0095e-7 * self.pressure) * self.temp**2
                    )
                )
            )
            * self.moleFractionNH3
        ) 
        return specificHeat


    def calcReactionRateConstant(self):
        return ko * np.exp(-E / (R * (self.temp)))
        # UNITS: J/mol/K

    def calcEquilibriumConstant(self):
        return 10 ** (
            -2.691122 * (np.log10(self.temp))
            - (5.519265e-5) * self.temp
            + (1.848863e-7) * (self.temp**2)
            + (2001.6 / self.temp)
            + 2.689
        )

    def calcRateOfReactionNH3(self):
        return (
            2
            * self.reactionRateConstant
            * (
                self.equilibriumConstant**2
                * self.activationCoefficientN2
                * (
                    (
                        self.activationCoefficientH2**3
                        / self.activationCoefficientNH3**2
                    )
                    ** alpha
                )
                - (
                    (
                        self.activationCoefficientNH3**2
                        / self.activationCoefficientH2**3
                    )
                    ** (1 - alpha)
                )
            )
        )

    def calcChangeInTempAcrossBed(self):
        return (
            self.effFactor * (-self.heatOfReaction) * A * self.rateOfReactionNH3
        ) / (self.F * self.specificHeat)

    def calcChangeOfConversionAcrossBed(self):
        return (self.effFactor * self.rateOfReactionNH3 * A) / (self.Fn2 * 2)

    def calcNewConversion(self):
        return self.conversionN2 + (stepSize * self.calcChangeOfConversionAcrossBed())

    def calcNewTemp(self):
        newTemp = self.temp + (stepSize * self.calcChangeInTempAcrossBed())
        return newTemp
    
    def calcNewEquilibriumConversion(self):
        return (self.equilibriumConstant*100)/(self.equilibriumConstant*100 + 1)


class ReactorUpdates(ReactorCalcs):
    def updateEffFactor(self):
        self.effFactor = self.calcEffFactor()

    def updateFugacityN2(self):
        self.fugacityN2 = self.calcFugacityN2()

    def updateFugacityH2(self):
        self.fugacityH2 = self.calcFugacityH2()

    def updateFugacityNH3(self):
        self.fugacityNH3 = self.calcFugacityNH3()

    def updateHeatOfReaction(self):
        self.heatOfReaction = self.calcHeatOfReaction()

    def updateSpecificHeat(self):
        self.specificHeat = self.calcSpecificHeat()

    def updateMoleFractionN2(self):
        self.moleFractionN2 = self.calcMoleFractionN2()

    def updateMoleFractionH2(self):
        self.moleFractionH2 = self.calcMoleFractionH2()

    def updateMoleFractionNH3(self):
        self.moleFractionNH3 = self.calcMoleFractionNH3()

    def updateMoleFractionAr(self):
        self.moleFractionAr = self.calcMoleFractionAr()

    def updateActivationCoefficientN2(self):
        self.activationCoefficientN2 = self.calcActivationCoefficientN2()

    def updateActivationCoefficientH2(self):
        self.activationCoefficientH2 = self.calcActivationCoefficientH2()

    def updateActivationCoefficientNH3(self):
        self.activationCoefficientNH3 = self.calcActivationCoefficientNH3()

    def updateReactionRateConstant(self):
        self.reactionRateConstant = self.calcReactionRateConstant()

    def updateEquilibriumConstant(self):
        self.equilibriumConstant = self.calcEquilibriumConstant()

    def updateRateOfReactionNH3(self):
        self.rateOfReactionNH3 = self.calcRateOfReactionNH3()

    def updateConversionN2(self):
        self.conversionN2 = self.calcNewConversion()

    def updateTemp(self):
        self.temp = self.calcNewTemp()

    def updateStep(self):
        self.steps = self.step + self._stepSize

    def updateEquilibriumConversion(self):
        self.equilibriumConversion = self.calcNewEquilibriumConversion()
