# -*- coding: utf-8 -*-
from numpy import array, reshape, ones, diag

class grandezas:

    def __init__(self,observedValues, uncertaintys, symbols):

        if observedValues is not None:
            self.obs = array(observedValues).transpose()
            self.NE = len(observedValues)
        if uncertaintys is not None:
            self.u = array(uncertaintys).transpose()
            self.variance = diag(uncertaintys**2)
        if symbols is not None:
            self.symbols = symbols
            self.NV = len(symbols)

    def _SETreconciled(self, reconciled):
        self.reconciled = array(reconciled).transpose()
        self.deviationsABS = self.obs-self.reconciled
        self.deviationsRELATIVE = 100*(self.obs-self.reconciled)/self.obs

    def _SETparameters(self, estimative):
        self.estimative = array(estimative).transpose()

    def _SETestimated(self, estimated):
        self.estimated = array(estimated).transpose()
        self.deviationsABS = self.obs - self.estimated
        self.deviationsRELATIVE = 100 * (self.obs - self.estimated) / self.obs


