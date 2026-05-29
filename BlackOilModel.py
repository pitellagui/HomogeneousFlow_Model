import numpy as np
import unitsconverter
from barril.units import Scalar

class GasPhase:
    """
    Classe para calcular compressibilidade, viscosidade e massa específica na fase gás.

    Parâmetros:
    ----------
    zcorrelation : str
        Informa a correlação a ser usada para calcular o fator Z - suporta somente Papay (1985)
    μcorrelation : str
        Informa a correlação a ser usada para calcular a viscosidade do gás - suporta Lee et al (1966) e Dempsey (1965)
    P : float
        Pressao absoluta
    Ppc : float
        Pressao pseudocritica
    Ppr : float
        Pressao pseudoreduzida
    T : float
        Temperatura
    Tpc : float
        Temperatura pseudocritica
    Tpr : float
        Temperatura pseudoreduzida
    """
    
    def __init__(self, zcorrelation, μcorrelation, dg=None,P=None, Ppc=None, Ppr=None, T=None, Tpc=None, Tpr=None):
        self.dg = dg
        self.P, self.Ppc, self.Ppr, self.T, self.Tpc, self.Tpr = P, Ppc, Ppr, T, Tpc, Tpr
        self.zcorrelation, self.μcorrelation = zcorrelation, μcorrelation
        self.Cpr, self.Cg, self.Z, self.ρ = None, None, None, None
        self.initialproperties()
        self.zfactcorrelselector()
        self.ρ_g()
        self.μfactcorrelselector()
        self.compressisoterm()

    # initial properties #

    def initialproperties(self):
        if self.dg is not None and self.Ppc is None and self.Tpc is None:
            if self.dg < 0.75:
                self.Ppc = 677 + 15*self.dg - 37.5*self.dg**2
                self.Tpc = 168 + 325*self.dg - 12.5*self.dg**2
            else:
                self.Ppc = 706 - 51.7*self.dg - 11.1*self.dg**2
                self.Tpc = 187 + 330*self.dg - 71.5*self.dg**2
        if self.Ppr is None and self.Tpr is None:
            self.Ppr = self.P/self.Ppc
            self.Tpr = self.T/self.Tpc
        if self.dg is None:
            pass
        self.Mg = self.dg * 28.96 # lb / lbmol
        self.Psc = 14.7 # psia
        self.Tsc = 60 # F
    def Mg(self):
        self.Mg = self.dg * 28.96 # lb / lbmol
    def ρ_g(self):
        R = 10.73 # psia·ft³/ (lb·mol·°R)
        self.ρ = (self.P * self.Mg)/(self.Z * R * self.T)
        
    def Bg(self):
        return self.Psc/(self.Tsc) * self.Z * unitsconverter.Temperature(self.T, 'R', 'F')/self.P
        
    def zfactcorrelselector(self):
        if self.zcorrelation == "Papay":
            self.papay()
        else:
            raise Exception('Correlação não encontrada')

    def μfactcorrelselector(self):
        if self.μcorrelation == "Lee":
            self.lee()
        elif self.μcorrelation == "Dempsey":
            self.dempsey()
        else:
            raise Exception('Correlação não encontrada')
            
    # Z factor correlations #

    def papay(self, Ppr=None):
        if Ppr is None:
            self.Z = 1 - (3.53*self.Ppr)/(10**(0.9813*self.Tpr)) + (0.274*self.Ppr**2)/(10**(0.8157*self.Tpr))
        else:
            return 1 - (3.53*Ppr)/(10**(0.9813*self.Tpr)) + (0.274*Ppr**2)/(10**(0.8157*self.Tpr))

    # μ correlations #
    
    def lee(self):
        x_v = 3.448 + (986.4/self.T) + (0.01009*self.Mg)
        y_v = 2.4 - (0.2*x_v)
        k_v = ((9.379+0.0160*self.Mg) * (self.T**1.5)) / (209.2 + (19.26*self.Mg) + self.T)

        factor = x_v * ((self.ρ/62.4)**y_v)
        self.μ = 1e-4 * k_v * np.exp(factor)
        
    def dempsey(self):
        self.μΔgN2, self.μΔgCO2, self.μΔgH2S = 0, 0, 0

        μgncorrigida = (1.709e-5 - 2.062e-6*self.dg)*unitsconverter.Temperature(self.T, 'R', 'F') + 8.188e-3-6.15e-3*np.log10(self.dg)
        
        self.μg = μgncorrigida + self.μΔgN2 + self.μΔgCO2 + self.μΔgH2S
        
        a0, a1, a2, a3, a4, a5 = -2.4621, 2.9705, -2.8626e-1, 8.0542e-3, 2.8086, -3.4980
        a6, a7, a8, a9, a10, a11 = 3.6037e-1, -1.0443e-2, -7.9339e-1, 1.3964, 1.4914e-1, 4.4102e-3
        a12, a13, a14, a15 = 8.3939e-2, -1.8641e-1, 2.0336e-2, -6.0958e-4

        lnexpression = a0 + a1*self.Ppr + a2*self.Ppr**2 + a3*self.Ppr**3 + self.Tpr*(a4 + a5*self.Ppr + a6*self.Ppr**2
                    + a7*self.Ppr**3) + self.Tpr**2*(a8 + a9*self.Ppr + a10*self.Ppr**2 +a11*self.Ppr**3) \
                    + self.Tpr**3*(a12 + a13*self.Ppr + a14*self.Ppr**2 + a15*self.Ppr**3)
                    
        self.μ = np.exp(lnexpression)*self.μg/self.Tpr

    def compressisoterm(self):
        Ppr = self.Ppr
        dx = 0.5 * 10 ** (-12)
        self.Cpr = 1/self.Ppr - 1/self.Z * (self.papay(Ppr+dx)-self.papay(Ppr))/dx
        self.Cg = self.Cpr/self.Ppc

    def output(self):
        print(f'Massa específica do gás: {self.ρ:.4f} lbm/ft³')
        print(f'Compressibilidade do gás: {self.Cg:.4e} psi⁻¹')
        print(f'Viscosidade do gás: { self.μ:.4f} cp')
        return self.ρ, self.Cg, self.μ

class OilPhase:
    """
    Classe para calcular compressibilidade, viscosidades e massa específica na fase óleo.

    Parâmetros:
    ----------
    soluratiocorrelation : str
        Informa a correlação a ser usada para calcular a razão de solubilidade - suporta somente Standing (1947)
    oilvolcorrelation : str
         Informa a correlação a ser usada para calcular a compressibilidade isotérmica do óleo e o fator volume-formação do óleo - suporta somente Standing (1947)
    oilviscosityselector : str
        Informa a correlação a ser usada para calcular as viscosidades do óleo morto e saturado - suporta somente Brill e Beggs (1974)
    do : float
        densidade relativa do óleo
    dg : float
        densidade relativa do gás
    P : float
        Pressao absoluta
    Pb : float
        Pressao de bolha
    T : float
        Temperatura
    units : list
        Lista com unidades de [pressão, temperatura], respectivamente
    """
    def __init__(self, soluratiocorrelation, oilvolcorrelation, oilviscosityselector, do=None, dg=None, P=None, Pb=None, T=None):
        
        self.do, self.dg, self.P, self.Pb, self.T = do, dg, P, Pb, T
        self.API, self.Rs = None, None
        self.soluratiocorrelation, self.oilvolcorrelation, self.oilviscosityselector = soluratiocorrelation, oilvolcorrelation, oilviscosityselector
        
        if Pb is None:
            self.Pb = self.standingbubblepress()
    
            
        self.initialproperties()
        self.soluratio()
        self.oilvolforcorrelationselector()
        self.ρ_o()
        self.oilvisccorrelationselector()
        self.oilclassifier()
        
    
    def initialproperties(self):
        self.API = 141.5/self.do - 131.5
        
    def standingbubblepress(self):
        self.Pb = 18.2 * ((self.Rs)/self.dg) * (10**0.00091 * self.T - 0.0125 * self.API) - 1.4

    
    # Solubility ratio #
    
    def soluratio(self):
        if self.P > self.Pb:
            self.Rs = self.dg * (((self.Pb/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83) # self.Rs = self.Rsb
            
        else:
            self.solurationcorrelationselector()
            
            
    def solurationcorrelationselector(self):
        if self.soluratiocorrelation == 'Standing':
            self.standingsoluration()
    
    def standingsoluration(self):
        self.Rs = self.dg * (((self.P/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83)
    
    # Oil volume formation factor #

    def oilvolforcorrelationselector(self):
        if self.oilvolcorrelation == "Standing":
            self.standingdRsdP()
            self.standingdBodP()
            # condição....
            if self.P > self.Pb:
                self.standingCo()
                self.standingBo()
            else:
                self.standingBo()
                self.standingCo()

    def standingdRsdP(self):
        self.dRsdP = self.dg * (((1/18.2) + 1.4) * 10**(0.0125*self.API - 0.00091*self.T))**(1/0.83)
    
    def standingBo(self):
        if self.P > self.Pb:
            self.Bo = self.Bob * np.exp(-self.Co * (self.P - self.Pb))
            
        else:
            self.Bo = 0.9759 + 0.00012*(self.Rs * (self.dg/self.do)**(0.5) + 1.25*self.T)**(1.2)
            
    def standingdBodP(self):
        self.dBodP = 0.00012*(self.dRsdP * (1.2)*(self.dg/self.do)**(0.5) + 1.25*self.T)**(0.2)
    
    def standingCo(self):
        if self.P >= self.Pb:
            self.ρ_ob()
            num = self.ρob + 0.004347 * (self.P - self.Pb) - 79.1
            dem = 0.0007141 * (self.P - self.Pb) - 12.938
            self.Co = 1e-6 * np.exp(num/dem)
        else:
            self.Co = -1/self.Bo * self.dBodP + GasPhase("Papay", "Lee", dg=self.dg, P=self.P, T=unitsconverter.Temperature(self.T, 'F', 'R')).Bg()/self.Bo * self.dRsdP
    # Oil density #

    def ρ_ob(self):
        raise NotImplementedError('Função não implementada!')
    
    def ρ_o(self):
        if self.P > self.Pb:
            self.ρo = self.ρob * np.exp(self.Co*(self.P-self.Pb))

        else:
            self.ρo = (62.4 * self.do + 0.0136 * self.Rs * self.dg)/self.Bo

    # Oil viscosity #
    
    def oilvisccorrelationselector(self):
        if self.oilviscosityselector == "Beggs":
           self.μbeggsod() 
           self.μbeggsob()
    
    def μbeggsod(self):
        self.μod = (10 ** ((10 ** (3.0324 - 0.02023 * self.API)) * self.T ** (-1.163))) - 1
    
    def μbeggsob(self):
        coeff = 10.715*(self.Rs+100)**(-0.515)
        exp = 5.44*(self.Rs+150)**(-0.338)
        self.μob = coeff*self.μod**exp

    def oilclassifier(self):
        if self.API < 40:
            type = "Black-Oil"
        elif 40 <= self.API <= 45:
            if self.Bo >= 2:
                type = "Óleo volátil"
            else:
                type = "Necessário determinar RGO para caracterizar o fluído"
        elif 45 < self.API <= 55:
            type = "Óleo volátil"

        else:
            type = "Modelo não aplicável - Óleo não é Black-Oil ou óleo volátil"
            raise Exception(type)
        print('* Tipo de fluido: ', type, '*')
    
    def output(self):
        print(f'Massa específica do óleo: {self.ρo:.2f} lbm/ft³')
        print(f'Compressibilidade do óleo: {self.Co:.4e} psi⁻¹')
        print(f'Viscosidade do óleo morto: {self.μod:.3f} cp')
        print(f'Viscosidade do óleo saturado: {self.μob:.3f} cp')
        return self.ρo, self.Co, self.μod, self.μob
    
def WaterPhase(P, T, S):
    S_value = Scalar('dimensionless', S, '%')
    salinity = Scalar('concentration', S_value.GetValue() * 10000, 'mg/L')

    rho_w_value = 62.368 + 0.438603 * S_value.GetValue('%') + 1.60074e-3 * S**2
    rho_w = Scalar('density', rho_w_value, 'lbm/ft3')

    A = 8.15839 + (-6.12265e-2) * T.GetValue('degF') + (1.91663e-4) * T.GetValue('degF')**2 + (-2.1654e-7) * T.GetValue('degF')**3
    B = 1.01021e-2 + (-7.44241e-5) * T.GetValue('degF') + (3.0553e-7) * T.GetValue('degF')**2 + (-2.94883e-10) * T.GetValue('degF')**3
    C = ((-9.02505) + (0.130237) * T.GetValue('degF') + (-8.53425e-4) * T.GetValue('degF')**2 + (2.34122e-6) * T.GetValue('degF')**3 + (-2.37049e-9) * T.GetValue('degF')**4) * 1e-7

    Rsw_value = A + B * P.GetValue('psia') + C * P.GetValue('psia')**2
    Rsw = Scalar('standard volume per standard volume', Rsw_value, 'scf/stb')

    if (1000 <= P.GetValue('psia') <= 20000) and (200 <= T.GetValue('degF') <= 270) and (0 <= salinity.GetValue('mg/L') <= 200):
        Cw_value = 1 / (7.033 * P.GetValue('psia') + 0.5415 * S_value.GetValue() - 537.0 * T.GetValue('degF') + 403300)
    else:
        Cw_value = 0

    Cw = Scalar('compressibility', Cw_value, '1/psi')

    VwT = -1.0001e-2 + 1.33391e-4 * T.GetValue('degF') + 5.50654e-7 * T.GetValue('degF')**2
    Vwp = -1.95301e-9 * P.GetValue('psia') * T.GetValue('degF') - 1.72834e-13 * P.GetValue('psia')**2 * T.GetValue('degF') - 3.58922e-7 * P.GetValue('psia') - 2.25341e-10 * P.GetValue('psia')**2

    Bw_value = (1 + VwT) * (1 + Vwp)
    Bw = Scalar('volume per standard volume', Bw_value, 'bbl/stb')

    A = 109.527 + (-8.40564) * S + (0.313314) * S**2 + (8.72213e-3) * S**3
    B = (-1.12166) + (2.63951e-2) * S + (-6.79461e-4) * S**2 + (-5.47119e-5) * S**3 + (-1.55586e-6) * S**4

    mi_w1_value = A * (T.GetValue('degF')**B) * (0.9994 + 4.0295e-5 * P.GetValue('psia') + 3.1062e-9 * P.GetValue('psia')**2)
    mi_w1 = Scalar('dynamic viscosity', mi_w1_value, 'cP')

    mi_w_value = mi_w1.GetValue('cP') * (0.9994 + 4.0295e-5 * P.GetValue('psia') + 3.1062e-9 * P.GetValue('psia')**2)
    mi_w = Scalar('dynamic viscosity', mi_w_value, 'cP')

    mi_w_ratio_value = mi_w.GetValue('cP') / mi_w1.GetValue('cP')
    mi_w_ratio = Scalar('dimensionless', mi_w_ratio_value, 'unitless')

    return rho_w, Rsw, Cw, Bw, mi_w1, mi_w, mi_w_ratio