# -*- coding: utf-8 -*-
from casadi import MX, vertcat, horzcat, sum1, Function, jacobian, inv

'''Module for symbolic variables creation'''

class SymbolicVariables:
    
    def __init__(self, data, Constraints=None, Model=None):
        """
        :param Data: problem data. For SDRPE problem, data is a list with four data sets (DR:data[0], PE:data[1], Param:data[2], xParam:data[3])
        :param Constrains: problem constraints
        :param Model: model which parameters will be estimated
        """
        if Constraints is not None:
            self.constraints = Constraints

        if Model is not None:
            self.model = Model

        # Data Reconciliation data
        self.symbols_DR = data[0].symbols
        self.NX = data[0].NE
        self.xObs = data[0].obs
        self.uxObs = data[0].u

        # y parameter estimation data
        self.NY = data[1].NE
        self.yObs = data[1].obs
        self.uyObs = data[1].u

        # Parameters data
        self.symbols_Param = data[2].symbols
        self.NP = data[2].NV

        # x parameter estimation data
        self.symbols_xPE = data[3].symbols

    def _SDRPEproblem(self,Method):

        """
        :param method: SDREP solution method; must be decoupled or coupled
        """

        self.xe = []; self.ux = []; self.xr = []; self.param = []; self.ye = []; self.uy = []; xaux = []

        self.method = Method

        # creation of the symbolic parameters
        for i in range(self.NP):
            self.param = vertcat(self.param, MX.sym(self.symbols_Param[i]))

        # creation of the symbolic variables for DR
        for i in range(self.NX):
            self.xe = vertcat(self.xe, MX.sym('xe' + str(i)))
            self.ux = vertcat(self.ux, MX.sym('ux' + str(i)))
            self.xr = vertcat(self.xr, MX.sym(self.symbols_DR[i] + 'r'))

        # creation of the symbolic variables for PE
        for i in range(self.NY):
            self.ye   = vertcat(self.ye, MX.sym('ye' + str(i)))
            self.uy   = vertcat(self.uy, MX.sym('uy' + str(i)))
            xaux = vertcat(xaux, MX.sym(self.symbols_xPE[i] + 'r'))

        # construction of the variable that contains the reconciled variables that are used in the model

        self.xePE = []
        self.positions = [] # positions in reconciled variables(xr) where are "x of parameter estimation"
        i = 0
        for k in range(self.NX):
            if i == xaux.size()[0]:
                i = i - 1
            if xaux[i].str() == self.xr[k].str():
                self.xePE = vertcat(self.xePE, self.xr[k])
                self.positions.append(k)
                i = i + 1

        # symbolic: model and constraints
        self.SymModel       = self.model(self.xr, self.param)
        self.SymConstraints = self.constraints(self.xr, self.param)

        # executable: model and constraint
        self.ExecModel       = Function('Model', [self.xr, self.param], [self.SymModel])
        self.ExecConstraints = Function('Constraints', [self.xr, self.param], [self.SymConstraints])

        # optimization problem definition
        if Method == 'coupled':
            self.SymFO      = sum1((self.xe - self.xr) ** 2 / (self.ux ** 2)) + sum1((self.ye - self.SymModel) ** 2 / (self.uy ** 2))
            self.ExecFO     = Function('Objective_Function', [self.xr, self.xe, self.ux, self.ye, self.uy], [self.SymFO])
            self.OptProblem = {'x': vertcat(self.xr, self.param),'p': vertcat(self.xe, self.ux, self.ye, self.uy),'f': self.SymFO,'g': self.SymConstraints}
            self.OptData    = vertcat(self.xObs, self.uxObs, self.yObs, self.uyObs)

        if Method == 'decoupled':
            self.SymFO_DR = sum1((self.xe - self.xr) ** 2 / (self.ux ** 2))  # objective function for DR
            self.ExecFO_DR = Function('Objective_Function_DR', [self.xr, self.xe, self.ux], [self.SymFO_DR])

            self.SymFO_PE = sum1((self.ye - self.SymModel) ** 2 / (self.uy ** 2))  # objective function for PE
            self.ExecFO_PE = Function('Objective_Function', [self.xr, self.ye, self.uy, self.param], [self.SymFO_PE])

            self.OptProblem_DR = {'x': vertcat(self.xr),'p': vertcat(self.xe, self.ux, self.param),'f': self.SymFO_DR,'g': self.SymConstraints}
            self.OptProblem_PE = {'x': vertcat(self.param),'p': vertcat(self.xr, self.ye, self.uy),'f': self.SymFO_PE}