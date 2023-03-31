# =========================================================================================================== #
# - Author :     Piotr T. Zaniewicz                                                                           #
# - Date   :     02/03/2023                                                                                   #
#                                                                                                             #
# - Description: - Required definitions for simulation .py script files
# =========================================================================================================== #

class ReactorBase:

    def __init__(self, stepSize, incomingTemp, pressure, bedLength, initialMoleFractionH2, initialMoleFractionN2,
                 initialMoleFractionNH3, initialMoleFractionAr, Fn2, F):
        #self._conversionN2 = [0]
        self._stepSize = stepSize
        self.incomingTemp = incomingTemp
        self._temp = [incomingTemp]
        self._steps = [0]
        self.pressure = pressure
        self.bedLength = bedLength
        self._effFactor = [0]
        self._heatOfReaction = [0]
        self._specificHeat = [0]

        self._moleFractionH2 = [initialMoleFractionH2]
        self._moleFractionN2 = [initialMoleFractionN2]
        self._moleFractionNH3 = [initialMoleFractionNH3]
        self._moleFractionAr = [initialMoleFractionAr]

        self._conversionN2 = [0]
        self._reactionRateConstant = [0]
        self._equilibriumConstant = [0]
        self._equilibriumConversion = [0]
        self._rateOfReactionNH3 = [0]

        self._activationCoefficientH2 = [0]
        self._activationCoefficientN2 = [0]
        self._activationCoefficientNH3 = [0]
        # self._activationCoefficientAr = [0]

        self._fugacityH2 = [self.calcFugacityH2()]
        self._fugacityN2 = [self.calcFugacityN2()]
        self._fugacityNH3 = [self.calcFugacityNH3()]
        # self._fugacityAr = [self.calcFugacityAr()]
        self.Fn2 = Fn2
        self.F = F

    @property
    def temp(self):
        return self._temp[-1]

    @temp.setter
    def temp(self, value):
        self._temp.append(value)

    @property
    def effFactor(self):
        return self._effFactor[-1]

    @effFactor.setter
    def effFactor(self, value):
        self._effFactor.append(value)

    @property
    def heatOfReaction(self):
        return self._heatOfReaction[-1]

    @heatOfReaction.setter
    def heatOfReaction(self, value):
        self._heatOfReaction.append(value)

    @property
    def specificHeat(self):
        return self._specificHeat[-1]

    @specificHeat.setter
    def specificHeat(self, value):
        self._specificHeat.append(value)

    @property
    def initialMoleFractionN2(self):
        return self._moleFractionN2[0]

    @property
    def initialMoleFractionH2(self):
        return self._moleFractionH2[0]

    @property
    def initialMoleFractionNH3(self):
        return self._moleFractionNH3[0]

    @property
    def initialMoleFractionAr(self):
        return self._moleFractionAr[0]

    @property
    def moleFractionH2(self):
        return self._moleFractionH2[-1]

    @moleFractionH2.setter
    def moleFractionH2(self, value):
        self._moleFractionH2.append(value)

    @property
    def moleFractionN2(self):
        return self._moleFractionN2[-1]

    @moleFractionN2.setter
    def moleFractionN2(self, value):
        self._moleFractionN2.append(value)

    @property
    def moleFractionNH3(self):
        return self._moleFractionNH3[-1]

    @moleFractionNH3.setter
    def moleFractionNH3(self, value):
        self._moleFractionNH3.append(value)

    @property
    def moleFractionAr(self):
        return self._moleFractionAr[-1]

    @moleFractionAr.setter
    def moleFractionAr(self, value):
        self._moleFractionAr.append(value)

    @property
    def conversionN2(self):
        return self._conversionN2[-1]

    @conversionN2.setter
    def conversionN2(self, value):
        self._conversionN2.append(value)

    @property
    def reactionRateConstant(self):
        return self._reactionRateConstant[-1]

    @reactionRateConstant.setter
    def reactionRateConstant(self, value):
        self._reactionRateConstant.append(value)

    @property
    def equilibriumConstant(self):
        return self._equilibriumConstant[-1]

    @equilibriumConstant.setter
    def equilibriumConstant(self, value):
        self._equilibriumConstant.append(value)

    @property
    def rateOfReactionNH3(self):
        return self._rateOfReactionNH3[-1]

    @rateOfReactionNH3.setter
    def rateOfReactionNH3(self, value):
        self._rateOfReactionNH3.append(value)

    @property
    def activationCoefficientH2(self):
        return self._activationCoefficientH2[-1]

    @activationCoefficientH2.setter
    def activationCoefficientH2(self, value):
        self._activationCoefficientH2.append(value)

    @property
    def activationCoefficientN2(self):
        return self._activationCoefficientN2[-1]

    @activationCoefficientN2.setter
    def activationCoefficientN2(self, value):
        self._activationCoefficientN2.append(value)

    @property
    def activationCoefficientNH3(self):
        return self._activationCoefficientNH3[-1]

    @activationCoefficientNH3.setter
    def activationCoefficientNH3(self, value):
        self._activationCoefficientNH3.append(value)

    @property
    def activationCoefficientAr(self):
        return self._activationCoefficientAr[-1]

    @activationCoefficientAr.setter
    def activationCoefficientAr(self, value):
        self._activationCoefficientAr.append(value)

    @property
    def fugacityH2(self):
        return self._fugacityH2[-1]

    @fugacityH2.setter
    def fugacityH2(self, value):
        self._fugacityH2.append(value)

    @property
    def fugacityN2(self):
        return self._fugacityN2[-1]

    @fugacityN2.setter
    def fugacityN2(self, value):
        self._fugacityN2.append(value)

    @property
    def fugacityNH3(self):
        return self._fugacityNH3[-1]

    @fugacityNH3.setter
    def fugacityNH3(self, value):
        self._fugacityNH3.append(value)

    @property
    def fugacityAr(self):
        return self._fugacityAr[-1]

    @fugacityAr.setter
    def fugacityAr(self, value):
        self._fugacityAr.append(value)

    @property
    def step(self):
        return self._steps[-1]

    @property
    def steps(self):
        return self._steps

    @steps.setter
    def steps(self, value):
        self._steps.append(value)

    @property
    def equilibriumConversion(self):
        return self._equilibriumConversion[-1]

    @equilibriumConversion.setter
    def equilibriumConversion(self, value):
        self._equilibriumConversion.append(value)
