# -*- coding: utf-8 -*-

"""Calculation engine to perform simultaneous data reconciliation and parameter estimation"""

from pandas import read_excel, DataFrame
from casadi import vertcat, qpsol, nlpsol
from numpy import array, mean, linspace, log10
import matplotlib.pyplot as plt
from Variables import grandezas
from Symbolic import SymbolicVariables
from scipy.stats import norm

class Rec:

    def __init__(self, constraints=None, model=None):
        """
           :param constrains: problem constraints
           :param Model: model which parameters will be estimated
        """

        self.constraints = constraints
        self.model = model

    def dataRead(self, dataType, file, sheet, **kwargs):
        """
           :param dataType: can be DR (data reconciliaiton) or PE (parameter estimation)
           :param file: excel file with the data set
           :param sheet: sheet name in the file
        """

        if not isinstance(dataType, str):
            raise TypeError('dataType should be a string')

        if dataType != 'DR' and dataType != 'PE':
            raise NameError('dataType should be DR or PE')

        # Kwargs to the pandas
        skiprows = kwargs.get('skiprows');
        usecols = kwargs.get('usecols')

        if kwargs.get('skipfooter'):
            skipfooter = kwargs.get('skipfooter')
        else:
            skipfooter = 0

        # data importation
        data = read_excel(file, sheet_name=sheet, skiprows=skiprows, skipfooter=skipfooter, usecols=usecols)

        # data frame creation
        if dataType == 'DR':

            # data frame creation for Data Reconciliation
            self.__dfDR = DataFrame(data, columns=['symbols', 'x', 'ux'])
            # variable creation
            self.xDR = grandezas(self.__dfDR['x'].values, self.__dfDR['ux'].values, self.__dfDR['symbols'].values)

        elif dataType == 'PE':
            # data frame creation for parameter estimation
            self.__dfPE = DataFrame(data, columns=['symbols_x', 'x', 'ux', 'symbols_y', 'y', 'uy', 'symbols_p'])
            # variables creation
            self.xPE = grandezas(self.__dfPE['x'].values, self.__dfPE['ux'].values, self.__dfPE['symbols_x'].values)
            self.yPE = grandezas(self.__dfPE['y'].values, self.__dfPE['uy'].values, self.__dfPE['symbols_y'].values)
            self.parameters = grandezas(None, None, self.__dfPE['symbols_p'].dropna().values)


    def simultaneousDataReconciliationParameterEstimation(self, initial_estimative, lbg=None, ubg=None,
                                                          algoritmo='ipopt', method='coupled',
                                                          first=None):
        """
           :param initial_estimative: list with the initial estimative for reconciled values and parameters
           :param lbg: list with the lower bound for constraints
           :param ubg: list with the upper bound for constraints
           :param algoritmo: optimization algorithm. Available: ipot; sqpmethod; qpoases; bonmin; qrqp
           :param method: the SDRPE solution method: Available: coupled; decoupled
           :param first: for decoupled method it's necessary to inform which problem will be solved first. "DR" or "PE".
        """

        # initial_estimative: first the estimates for the reconciled. Then the estimates for the parameters.

        if self.model is None:
            raise ValueError('To execute SDRPE it is necessary to inform the model of the depended quantities')
        if self.constraints is None:
            raise ValueError('To execute Data Reconciliation it is necessary to inform the problem constraints')
        if method != 'coupled' and method != 'decoupled':
            raise ValueError('just two methods are available: (i) decoupled and (ii) coupled.')

        # Initialization of symbolic variables class
        self.__SymbolicProblem = SymbolicVariables(data=[self.xDR, self.yPE, self.parameters, self.xPE],
                                                   Constraints=self.constraints, Model=self.model)
        # Construction of the symbolic problem
        self.__SymbolicProblem._SDRPEproblem(Method=method)

        # ----------------------------
        #    coupled problem
        # ----------------------------

        if method == 'coupled':

            if algoritmo == 'qpoases' or algoritmo == 'qrqp':
                S = qpsol('S', algoritmo, self.__SymbolicProblem.OptProblem)
            elif algoritmo == 'ipopt' or algoritmo == 'bonmin' or algoritmo == 'sqpmethod':
                S = nlpsol('S', algoritmo, self.__SymbolicProblem.OptProblem)

                # solving the problem
            Opt = S(x0=initial_estimative, p=self.__SymbolicProblem.OptData, lbg=lbg,
                    ubg=ubg)

            # setting residuals and reconciled values
            self.xDR._SETreconciled(Opt['x'].elements()[0:self.xDR.NE])

            # setting estimated parameters
            self.parameters._SETparameters(Opt['x'].elements()[self.xDR.NE:])

            # setting model estimatives
            self.yPE._SETestimated(
                self.__SymbolicProblem.ExecModel(self.xDR.obs, self.parameters.estimative).elements())

            self.FobjOtim = float(Opt['f'])  # objective function value at optimal point
            self.imbalances = array(Opt['g'].elements())  # imbalances


        # ------------------------
        #    Decoupled problem
        # ------------------------

        elif method == 'decoupled':

            if first is None:
                raise ValueError('To use decoupled method it is necessary to inform which problem will be solved first')

            if algoritmo == 'qpoases' or algoritmo == 'qrqp':
                S_DR = qpsol('S', algoritmo, self.__SymbolicProblem.OptProblem_DR)
                S_PE = qpsol('S', algoritmo, self.__SymbolicProblem.OptProblem_PE)
            elif algoritmo == 'ipopt' or algoritmo == 'bonmin' or algoritmo == 'sqpmethod':
                S_DR = nlpsol('S', algoritmo, self.__SymbolicProblem.OptProblem_DR)
                S_PE = nlpsol('S', algoritmo, self.__SymbolicProblem.OptProblem_PE)

            # initial estimatives
            x0_DR = initial_estimative[0:self.xDR.NE] # intial estimative for the reconciled values
            x0_PE = initial_estimative[self.xDR.NE:]  # intial estimative for the parameters

            fobj_DR = 10e30
            fobj_PE = 10e30  # initial values for the objective functions
            tol_DR = 10
            tol_PE = 10  # initial values for the tolerance
            optimized = False  # the problem are still not optimized
            i = 0  # iterations number

            # -----------------------------------
            # Decoupled problem solve algorithm
            # -----------------------------------

            if first == 'PE':

                while not optimized:

                    if tol_DR >= 0.0001 or tol_PE >= 0.0001:

                        '''Solving parameter estimation problem'''
                        Opt_PE = S_PE(x0=x0_PE, p=vertcat(x0_DR, self.yPE.obs, self.yPE.u))
                        tol_PE = abs(fobj_PE - float(Opt_PE['f']))  # tolerance calculation
                        fobj_PE = float((Opt_PE['f']))  # new objective function value
                        x0_PE = Opt_PE['x'].elements()  # new estimative for parameters

                        '''Solving data reconciliation problem'''
                        Opt_DR = S_DR(x0=x0_DR, p=vertcat(self.xDR.obs, self.xDR.u, x0_PE), lbg=lbg, ubg=ubg)
                        tol_DR = abs(fobj_DR - float(Opt_DR['f']))  # toleracne calculation
                        fobj_DR = float((Opt_DR['f']))  # new objective function value
                        x0_DR = Opt_DR['x'].elements()  # new estimative for reconciled values

                    else:
                        optimized = True

                    i = i + 1
                    if i > 500:
                        optimized = True

            elif first == 'DR':

                while not optimized:

                    if tol_DR >= 0.0001 or tol_PE >= 0.0001:

                        '''data reconciliation problem'''
                        Opt_DR = S_DR(x0=x0_DR, p=vertcat(self.xDR.obs, self.xDR.u, x0_PE), lbg=lbg, ubg=ubg) # solve the problem
                        tol_DR = abs(fobj_DR - float(Opt_DR['f']))  # tolerance calculation
                        fobj_DR = float((Opt_DR['f']))  # new objective function value
                        x0_DR = Opt_DR['x'].elements()  # new estimative for reconciled values

                        '''parameter estimation problem'''
                        Opt_PE = S_PE(x0=x0_PE, p=vertcat(x0_DR, self.yPE.obs, self.yPE.u)) # solve the problem
                        tol_PE = abs(fobj_PE - float(Opt_PE['f']))  # tolerance calculation
                        fobj_PE = float((Opt_PE['f']))  # new objective function value
                        x0_PE = Opt_PE['x'].elements()  # new estimative for parameters

                    else:
                        optimized = True

                    i = i + 1
                    if i > 500:
                        optimized = True

            # setting residuals and reconciled values
            self.xDR._SETreconciled(Opt_DR['x'].elements())

            # setting estimated parameters
            self.parameters._SETparameters(Opt_PE['x'].elements())

            # setting model estimatives
            self.yPE._SETestimated(
                self.__SymbolicProblem.ExecModel(self.xDR.obs, self.parameters.estimative).elements())

            self.FobjOtim_DR = float(Opt_DR['f'])  # DR objective function value at optimal point
            self.FobjOtim_PE = float(Opt_PE['f'])  # PE objective function value at optimal point
            self.imbalances = array(Opt_DR['g'].elements())  # imbalances


    def report(self, nameReport):

        with open(nameReport + '.txt', 'w') as arquivo:
            arquivo.write(('{:#^' + str(max([107, 36])) + '}' + "\n").format(nameReport))


            arquivo.write('\n')
            arquivo.write('    Data reconciliation:\n')
            arquivo.write('\n')
            arquivo.write(
                '      Quantity   |' + '      X obs     |' + '      U exp     |' + '      X rec     |' + '   Deviations   |' + '   Deviations(%)  |')
            arquivo.write('\n')
            arquivo.write('    {:-^118}'.format(''))

            for i in range(self.xDR.NE):
                arquivo.write('\n')

                arquivo.write('       ' + '{0:^5}'.format(self.xDR.symbols[i]) + '     |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.xDR.obs[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.xDR.u[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.xDR.reconciled[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:+4.5f}'.format(self.xDR.deviationsABS[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:+4.5f}'.format(self.xDR.deviationsRELATIVE[i])) + '%    |')

            # write constraints

            arquivo.write('\n')
            arquivo.write('\n')
            arquivo.write('    Constraints:\n')
            arquivo.write('{0:^101}'.format('     Imbalance') + '{: ^8}'.format('|') + '  Value' + '     |')
            arquivo.write('\n')
            arquivo.write('       {:-^115}'.format(''))

            for i in range(len(self.imbalances)):
                arquivo.write('\n')
                arquivo.write('          ' + '{0:94}'.format(str(self.__SymbolicProblem.SymConstraints[i])) + '|')
                arquivo.write('  {0:+4.8f}'.format(self.imbalances[i]) + '   |')


            arquivo.write('\n')
            arquivo.write('\n')
            arquivo.write('    Parameter Estimation:\n')
            arquivo.write('\n')
            arquivo.write(
                '      Quantity   |' + '      y obs     |' + '      U obs     |' + '      y est     |' + '   Deviations   |' + '   Deviations(%)  |' )
            arquivo.write('\n')
            arquivo.write('    {:-^118}'.format(''))

            for i in range(self.yPE.NE):
                arquivo.write('\n')

                arquivo.write('       ' + '{0:^5}'.format(self.yPE.symbols[i]) + '     |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.yPE.obs[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.yPE.u[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:4.5f}'.format(self.yPE.estimated[i])) + '   |')
                arquivo.write('   ' + '{0:^10}'.format('{0:+4.5f}'.format(self.yPE.deviationsABS[i])) + '   |')
                arquivo.write(
                    '   ' + '{0:^10}'.format('{0:+4.5f}'.format(self.yPE.deviationsRELATIVE[i])) + '%    |')


            # writing parameters
            arquivo.write('\n')
            arquivo.write('\n')
            arquivo.write('    Parameters value:\n')
            arquivo.write('\n')
            arquivo.write('       Parameter  |' + '      Value     |')
            arquivo.write('\n')
            arquivo.write('    {:-^51}'.format(''))

            for i in range(self.parameters.NV):
                arquivo.write('\n')

                arquivo.write('        ' + '{0:^5}'.format(self.parameters.symbols[i]) + '     |')
                arquivo.write('   ' + '{0:^10}'.format('{:+10.3e}'.format(self.parameters.estimative[i])) + '   |')


            arquivo.write('\n')
            arquivo.write('\n')
            arquivo.write('\n')
            if self.__SymbolicProblem.method == 'coupled':
                arquivo.write('    * Objective function value at optimal point = ')
                arquivo.write(str(self.FobjOtim))
            else:
                arquivo.write('    * Objective DR function value at optimal point = ')
                arquivo.write(str(self.FobjOtim_DR))
                arquivo.write('\n')
                arquivo.write('    * Objective PE function value at optimal point = ')
                arquivo.write(str(self.FobjOtim_PE))
            arquivo.close()

    def charts(self):

        aux = []
        for i in range(self.xDR.NE):  # create a column of samples
            aux.append(i + 1)
        self.__dfDR['samples'] = aux

        '''Scatter Plot'''

        # absolute deviation
        for i, variables in enumerate(self.xDR.symbols):
            plt.scatter(i, self.xDR.deviationsABS[i], s=30, c='DarkBlue')
            if self.xDR.NE < 20:
                plt.text(i + 0.2, self.xDR.deviationsABS[i], variables)
        plt.grid()
        plt.axhline(mean(self.xDR.deviationsABS), label='residues mean', color='red', zorder=6)
        plt.legend(loc='best')
        plt.xlabel('Samples')
        plt.ylabel('Residues')
        plt.savefig('Absolute residues')
        plt.close()

        # relative deviation
        for i, variables in enumerate(self.xDR.symbols):
            plt.scatter(i, self.xDR.deviationsRELATIVE[i], s=30, c='green')
            if self.xDR.NE < 20:
                plt.text(i + 0.2, self.xDR.deviationsRELATIVE[i] + 0.2, variables)
        plt.grid()
        plt.axhline(mean(self.xDR.deviationsRELATIVE), label='residues mean', color='red', zorder=6)
        plt.legend(loc='best')
        plt.xlabel('Samples', fontsize=15)
        plt.ylabel('Residues', fontsize=15)
        a = min(self.xDR.deviationsRELATIVE)
        b = max(self.xDR.deviationsRELATIVE)
        plt.ylim(a - 0.1 * abs(a), b + 0.1 * abs(b))
        plt.savefig('Relative residues')
        plt.close()

        # Observed by reconciled
        for i, variables in enumerate(self.xDR.symbols):
            plt.scatter(self.xDR.obs[i], self.xDR.reconciled[i], s=30, c='DarkBlue')
        plt.grid()
        xaux = linspace(min(self.xDR.obs) * 0.8, max(self.xDR.reconciled) * 1.2, self.xDR.NE)
        plt.plot(xaux, xaux, label='observed = reconciled', color='red', zorder=6)
        plt.legend(loc='best')
        plt.xlabel('Observed data')
        plt.ylabel('Reconciled data')
        plt.savefig('Observed by reconciled')
        plt.close()

        # deviations histogram
        plt.hist(self.xDR.deviationsRELATIVE, bins=int(1 + 3.3 * log10(self.xDR.NE)))  # Sturges rule
        # plot pdf
        xmin, xmax = plt.xlim()
        x = linspace(xmin, xmax, 1000)
        mu, std = norm.fit(self.xDR.deviationsRELATIVE)  # mean and std deviation
        p = norm.pdf(x, mu, std)*100  # se for relativo tem que multuplicar por 100
        plt.plot(x, p, 'k', linewidth=2, c='red', label='normal PDF')

        plt.legend(loc='best')
        plt.xlabel('Residues', fontsize=15)
        plt.ylabel('Frequency', fontsize=15)
        plt.savefig('residues histogram')
        plt.close()

        # deviation autocorrelation
        plt.acorr(self.xDR.deviationsABS, normed=True, maxlags=None, color='blue')
        plt.axhline(1.96 / (len(self.xDR.deviationsABS) ** 0.5), ls='--', color='red')  # intervalo de confiança
        plt.axhline(-1.96 / (len(self.xDR.deviationsABS) ** 0.5), ls='--', color='red')  # intervalo de confiança
        plt.grid()
        plt.xlabel('Lag', fontsize=15)
        plt.ylabel('Autocorrelation of the residues', fontsize=15)
        plt.xlim(0, self.xDR.NE)
        plt.savefig('residues autocorrelation')
        plt.close()

        # comparing the observations with the estimates
        plt.scatter(self.xPE.obs, self.yPE.estimated, s=30, label='Estimated', c='DarkBlue')
        plt.scatter(self.xPE.obs, self.yPE.obs, s=30, label='Observed', c='Red')
        plt.grid()
        plt.legend(loc='best')
        plt.xlabel('Temperature / (°C)', fontsize=15)
        plt.ylabel('Specific heat / (kJ / (kg °C))', fontsize=15)
        plt.savefig('Temperature x Specific heat')
        plt.close()



