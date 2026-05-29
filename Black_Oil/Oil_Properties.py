import numpy as np
from barril.units import Array, Scalar

def Calculate_oil_properties(P, T, do, rho_o, Pb, Rs,dg,Bg):
    rho_water = Scalar('density', 999.1, 'kg/m3')
    if do is None and rho_o is not None:
        do_value = rho_o.GetValue()/rho_water.GetValue()
        do = Scalar('dimensionless', do_value, 'unitless')

    API_value = 141.5/do.GetValue() - 131.5
    API = Scalar('dimensionless', API_value, 'unitless')

    rho_ob = Scalar('density', do.GetValue() *
                    rho_water.GetValue('kg/m3'), 'kg/m3')

    if Pb is None and Rs is not None:
        Pb_values = []
        for i in range(len(P)):
            a = 0.00091 * T.GetValues('degF')[i] - 0.0125 * API.GetValue()
            Calculate_Pb = 18.2 * \
                (((Rs.GetValues('scf/stb')
                 [i] / dg.GetValue()) ** 0.83) * (10 ** a) - 1.4)
        Pb_values.append(Calculate_Pb)
        Pb = Array('pressure', Pb_values, 'psia')

    Rs_values = []
    Rsb_values = []
    Co_values = []
    Bo_values = []
    Bob_values = []
    Bt_values = []
    rho_o_values = []
    mi_od_values = []
    mi_ob_values = []
    mi_o_values = []

    for i in range(len(P)):
        Calculate_Rsb = dg.GetValue() * (((Pb.GetValues('psia')[i] / 18.2) + 1.4) * 10**(
            0.0125 * API.GetValue() - 0.00091 * T.GetValues('degF')[i]))**(1/0.83)
        if P[i] <= Pb[i]:
            Calculate_Rs = dg.GetValue() * (((P.GetValues('psia')[i] / 18.2) + 1.4) * 10**(
                0.0125 * API.GetValue() - 0.00091 * T.GetValues('degF')[i]))**(1/0.83)
        else:
            Calculate_Rs = Calculate_Rsb
        Rs_values.append(Calculate_Rs)
        Rsb_values.append(Calculate_Rsb)
    Rs = Array('standard volume per standard volume',
               np.array(Rs_values), 'scf/stb')
    Rsb = Array('standard volume per standard volume',
                np.array(Rsb_values), 'scf/stb')

    for i in range(len(P)):
        if P[i] < Pb[i]:
            Bo_valuess = 0.9759 + 0.00012 * ((Rs.GetValues('scf/stb')[i] * (
                dg.GetValue()/do.GetValue())**0.5 + 1.25 * T.GetValues('degF')[i])**1.2)
            dRs_dP = dg.GetValue() * (1/0.83) * (1/18.2) * np.power(10**(0.0125*API.GetValue() - 0.00091 *
                                                                         T.GetValues('degF')[i]), 1/0.83) * np.power((P.GetValues('psia')[i]/18.2 + 1.4), (1/0.83 - 1))
            dBo_dP = 0.00012 * 1.2 * (dg.GetValue()/do.GetValue())**0.5 * ((Rs.GetValues('scf/stb')[
                i]*(dg.GetValue()/do.GetValue())**0.5 + 1.25*T.GetValues('degF')[i])**0.2) * dRs_dP
            Calculate_Co = (1/Bo_valuess)*dBo_dP + \
                (Bg.GetValues('bbl/stb')[i]/Bo_valuess)*dRs_dP
        else:
            Calculate_Co = 1e-6*np.exp((rho_ob.GetValue('lbm/ft3') + 0.004347*(P.GetValues('psia')[i] - Pb.GetValues(
                'psia')[i]) - 79.1)/(0.0007141*(P.GetValues('psia')[i] - Pb.GetValues('psia')[i]) - 12.938))
        Co_values.append(Calculate_Co)
    Co = Array('compressibility', np.array(Co_values), '1/psi')

    for i in range(len(P)):
        Calculate_Bob = 0.9759 + 0.00012 * ((Rsb.GetValues('scf/stb')[i] * (
            dg.GetValue()/do.GetValue())**0.5 + 1.25 * T.GetValues('degF')[i])**1.2)
        Bob_values.append(Calculate_Bob)
        Bob = Array('standard volume per standard volume',
                    Bob_values, 'stb/stb')
        if P[i] > Pb[i]:
            Calculate_Bo = Bob.GetValues('stb/stb')[i] * np.exp(-Co.GetValues(
                '1/psi')[i] * (P.GetValues('psia')[i] - Pb.GetValues('psia')[i]))
        else:
            Calculate_Bo = 0.9759 + 0.00012 * ((Rs.GetValues('scf/stb')[i] * (
                dg.GetValue()/do.GetValue())**0.5 + 1.25 * T.GetValues('degF')[i])**1.2)
        Bo_values.append(Calculate_Bo)
    Bo = Array('volume per standard volume', np.array(Bo_values), 'bbl/stb')

    for i in range(len(P)):
        if P[i] >= Pb[i]:
            Calculate_Bt = Bo.GetValues('bbl/stb')[i]
        else:
            Calculate_Bt = Bo.GetValues('bbl/stb')[i] + (Rsb.GetValues(
                'scf/stb')[i] - Rs.GetValues('scf/stb')[i]) * Bg.GetValues('bbl/stb')[i]
        Bt_values.append(Calculate_Bt)
    Bt = Array('volume per standard volume', np.array(Bt_values), 'bbl/stb')

    for i in range(len(P)):
        if P[i] <= Pb[i]:
            Calculate_rho_o = (62.4 * do.GetValue() + 0.0136 * Rs.GetValues(
                'scf/stb')[i] * dg.GetValue()) / Bo.GetValues('bbl/stb')[i]
        else:
            Calculate_rho_o = rho_ob.GetValue('lbm/ft3') * np.exp(Co.GetValues(
                '1/psi')[i] * (P.GetValues('psia')[i] - Pb.GetValues('psia')[i]))
        rho_o_values.append(Calculate_rho_o)
    rho_o = Array('density', np.array(rho_o_values), 'lbm/ft3')

    for i in range(len(P)):
        A = 10 ** (0.43 + (8.33 / API.GetValue()))
        Calculate_mi_od = 0.32 + \
            (1.8e7 / (API.GetValue() ** 4.53)) * \
            ((360 / (T.GetValues('degR')[i] - 260)) ** A)
        mi_od_values.append(Calculate_mi_od)
    mi_od = Array('dynamic viscosity', mi_od_values, 'cP')

    for i in range(len(P)):
        a = 10 ** (-7.4e-4 * Rsb.GetValues('scf/stb')
                   [i] + 2.2e-7 * Rsb.GetValues('scf/stb')[i]**2)
        b = 0.68 / (10 ** (8.62e-5 * Rsb.GetValues('scf/stb')[i])) + 0.25 / (10 ** (
            1.1e-3 * Rsb.GetValues('scf/stb')[i])) + 0.062 / (10 ** (3.74e-3 * Rsb.GetValues('scf/stb')[i]))
        Calculate_mi_ob = a * (mi_od.GetValues('cP')[i] ** b)
        mi_ob_values.append(Calculate_mi_ob)
    mi_ob = Array('dynamic viscosity', mi_ob_values, 'cP')

    for i in range(len(P)):
        if P[i] <= Pb[i]:
            Calculate_mi_o = mi_ob[i]
        else:
            Calculate_mi_o = mi_ob[i] + (0.001 * (P.GetValues('psi')[i] - Pb.GetValues('psi')[i])) * (
                0.024 * mi_ob.GetValues('cP')[i]**1.6 + 0.038 * mi_ob.GetValues('cP')[i]**0.56)
        mi_o_values.append(Calculate_mi_o)
    mi_o = Array('dynamic viscosity', mi_o_values, 'cP')

    return do, API, Pb, Rs, Co, Bo, Bt, rho_o, mi_od, mi_ob, mi_o