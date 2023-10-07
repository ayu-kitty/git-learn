# File: nca
# Author: jiang tao
# Time: 2023/7/6 8:44
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.linear_model import LinearRegression


class PharmacokineticsException(Exception):
    pass


class RawDataError(PharmacokineticsException):
    pass


class NonCompart():
    """Non-compartmental Analysis for Pharmacokinetic Data"""
    def __init__(self, x, y, dose=0, r2Adj=0.3):
        """
        :param x: times[ndarray], h
        :param y: concentration[ndarray], ng/mL
        :param dose: give amount, mg
        :param r2Adj: Minimum adjusted R-square value to determine terminal slope automatically
        """
        self.x = x
        self.y = y
        self.dose = dose
        self.r2Adj = r2Adj

    def getCTmax(self):
        Cmax = self.y[self.y.argmax()]
        Tmax = self.x[self.y.argmax()]
        return Cmax, Tmax

    def calcAUCObs(self):
        """AUC using all the given points, including trailing zero concentrations"""
        return trapezoid(self.y, self.x)

    def calcAUCInf(self, ke, Tlast):
        """AUC from 0 to Time last, to infinity"""
        idx = int(np.argwhere(self.x == Tlast)[0][0])
        Clast = self.y[idx]
        AUClast = trapezoid(self.y[:idx+1], self.x[:idx+1])
        AUCinf = AUClast + Clast/ke
        return AUCinf

    def findBestSlope(self):
        """lambda_z negative of the best-fit terminal slope"""
        cmax_idx = self.y.argmax()
        Cmax = self.y[cmax_idx]
        x1 = self.x[cmax_idx:]
        y1 = self.y[cmax_idx:]
        # removes zeros
        nonzero_idx = np.flatnonzero(y1)
        times, cons = x1[nonzero_idx], y1[nonzero_idx]
        candidates = []
        n_max = times.size
        if n_max < 3:
            return
        for i in range(3, n_max):
            estimator = LinearRegression()
            to, Co = times[-i:][-1], cons[-i:][-1]
            X, y = times[-i:].reshape(-1, 1), np.log(cons[-i:] / Cmax).reshape(-1, 1)
            estimator.fit(X, y)
            coef = estimator.coef_[0][0]
            intercept = estimator.intercept_[0]
            r2 = estimator.score(X, y)
            r2Adj = 1 - (1 - r2) * (i - 1) / (i - 2)
            if coef < 0:
                candidates.append((coef, intercept, r2Adj, to, Co, estimator, i))
        candidates.sort(key=lambda x: x[2], reverse=True)
        if len(candidates) == 0:
            return
        r2AdjMax, r2AdjMax_n_points = candidates[0][2], candidates[0][-1]
        if r2AdjMax <= self.r2Adj:
            return
        if len(candidates) >= 2:
            r2AdjMax1, r2AdjMax1_n_points = candidates[1][2], candidates[1][-1]
            if abs(r2AdjMax-r2AdjMax1) <= 0.0001 and r2AdjMax1_n_points > r2AdjMax_n_points:
                return candidates[1]
        return candidates[0]

    def calcMRTObs(self, Tlast):
        """mean residence time (MRT) to TLST, for extravascular administration"""
        idx = int(np.argwhere(self.x == Tlast)[0][0])
        AUClast = trapezoid(self.y[:idx+1], self.x[:idx+1])
        ct_data = self.x[:idx+1] * self.y[:idx+1]
        AUMCObs = trapezoid(y=ct_data, x=self.x[:idx+1])
        MRTObs = AUMCObs / AUClast
        return MRTObs

    def calcMRTInf(self):
        """mean residence time (MRT) infinity using CLST, for extravascular administration"""
        # AUMCinf = AUMCObs  + Co/(ke*ke) + to*Co/ke
        ...

    def calcOtherParams(self, ke, AUCinf):
        """T1/2z(hr), Vz/F(L/kg), CLz/F(L/hr/kg)"""
        halfLife = np.log(2)/ke
        VzF = self.dose * 1e7 / (ke * AUCinf)
        CLzF = self.dose * 1e7 / AUCinf
        return halfLife, VzF, CLzF

    def NCA(self):
        Cmax, Tmax = self.getCTmax()
        AUCObs = self.calcAUCObs()
        params = self.findBestSlope()
        if params is None:
            return np.zeros(9)
        coef, intercept, r2Adj, to, Co, estimator, n_points = params
        ke, Tlast = -coef, to
        AUCInf = self.calcAUCInf(ke, Tlast)
        MRTObs = self.calcMRTObs(Tlast)
        halfLife, VzF, CLzF = self.calcOtherParams(ke, AUCInf)
        return Cmax, Tmax, ke, AUCObs, AUCInf, VzF, CLzF, MRTObs, halfLife

if __name__ == '__main__':
    ...
    # data = pd.read_excel("winolin.xlsx", sheet_name=5)
    # x = data["T"].to_numpy()
    # y = data["C"].to_numpy()
    # data = pd.read_excel(r"D:\luming_project\药动学资料-LK\药动学资料-LK\Winnonlin分析数据\Winnonlin.xls", sheet_name="CA", header=0,
    #                      skiprows=[1, ])
    # grps = data.subject.drop_duplicates()
    # data = data.set_index(['subject'])
    # data1 = []
    # import time
    # start = time.time()
    # for g in grps:
    #     _ = data.loc[g, :].reset_index(drop=True)
    #     x = _["time"].to_numpy()
    #     y = _["concentration"].to_numpy()
    #     nca = NonCompart(x, y)
    #     row = nca.NCA()
    #     data1.append(row)
    # end = time.time()
    # print(end-start, "seconds")
    # df = pd.DataFrame(data=data1, columns=["Cmax", "Tmax", "ke", "AUCObs", "AUCInf", "VzF", "CLzF", "MRTObs", "halfLife"], index=grps)
    # df.to_excel("pharmCARes.xlsx")

