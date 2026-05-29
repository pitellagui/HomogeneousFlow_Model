import numpy as np
from barril.units import Scalar


def Calculate_gas_properties(P, T, dg, rho_g):
    if dg is None and rho_g is not None:
        rho_ar = Scalar('density', 1.225, 'kg/m3')
        rho_g = Scalar('density', rho_g, 'kg/m3')
        dg = Scalar('dimensionless', rho_g.GetValue('kg/m3') / rho_ar.GetValue('kg/m3'), 'unitless')

    if dg is None:
        raise ValueError("É necessário informar dg ou rho_g.")

    Mar = 28.96
    Mg_value = Mar * dg.GetValue()
    Mg = Scalar('mass per mol', Mg_value, 'lb/lbmol')

    if float(dg.GetValue()) < 0.75:
        P_pc_val = 677 + 15.0 * dg.GetValue() - 37.5 * dg.GetValue()**2
        T_pc_val = 168 + 325 * dg.GetValue() - 12.5 * dg.GetValue()**2
    else:
        P_pc_val = 706 - 51.7 * dg.GetValue() - 11.1 * dg.GetValue()**2
        T_pc_val = 187 + 330 * dg.GetValue() - 71.5 * dg.GetValue()**2

    P_pc_value = Scalar('pressure', P_pc_val, 'psia')
    T_pc_value = Scalar('temperature', T_pc_val, 'degR')

    P_pc = Scalar('pressure', P_pc_value.GetValue('Pa'), 'Pa')
    T_pc = Scalar('temperature', T_pc_value.GetValue('K'), 'K')

    P_pr_value = P.GetValue('Pa') / P_pc_value.GetValue('Pa')
    P_pr = Scalar('dimensionless', P_pr_value, 'unitless')

    T_pr_value = T.GetValue('K') / T_pc_value.GetValue('K')
    T_pr = Scalar('dimensionless', T_pr_value, 'unitless')

    R = Scalar('molar gas constant', 8.314, 'J/mol.K')

    Z_value = 1 - (3.53 * P_pr.GetValue()) / (10**(0.9813 * T_pr.GetValue())) + (0.274 * P_pr.GetValue()**2) / (10**(0.8157 * T_pr.GetValue()))
    Z = Scalar('dimensionless', Z_value, 'unitless')

    rho_g_value = (P.GetValue('Pa') * Mg.GetValue('kg/mol')) / (Z.GetValue() * R.GetValue('J/mol.K') * T.GetValue('K'))
    rho_g = Scalar('density', rho_g_value, 'kg/m3')

    dZdP_pr = -(3.53 / (10**(0.9813 * T_pr.GetValue()))) + (2 * 0.274 * P_pr.GetValue() / (10**(0.8157 * T_pr.GetValue())))

    C_pr_value = (1 / P_pr.GetValue()) - ((1 / Z.GetValue()) * dZdP_pr)
    C_pr = Scalar('dimensionless', C_pr_value, 'unitless')

    Cg_value = C_pr.GetValue() / P_pc.GetValue('Pa')
    Cg = Scalar('compressibility', Cg_value, '1/Pa')

    Psc = Scalar('pressure', 14.7, 'psia')
    Tsc = Scalar('temperature', 60, 'degF')

    Bg_value = (Psc.GetValue('psia') / Tsc.GetValue('degF')) * (Z.GetValue() * T.GetValue('degF') / P.GetValue('psia'))
    Bg = Scalar('volume per standard volume', Bg_value, 'm3/sm3')

    Eg_value = 1 / Bg.GetValue('m3/sm3')
    Eg = Scalar('standard volume per volume', Eg_value, 'sm3/m3')

    xv = 3.448 + 986.4 / T.GetValue('degR') + 0.01009 * Mg.GetValue('lb/lbmol')
    yv = 2.4 - 0.2 * xv
    kv = ((9.379 + 0.0160 * Mg.GetValue('lb/lbmol')) * (T.GetValue('degR')**1.5)) / (209.2 + 19.26 * Mg.GetValue('lb/lbmol') + T.GetValue('degR'))

    mi_g_value = 1e-4 * kv * np.exp(xv * (rho_g.GetValue('lbm/ft3') / 62.4)**yv)
    mi_g = Scalar('dynamic viscosity', mi_g_value, 'cP')

    return Mg, dg, P_pc, T_pc, P_pr, T_pr, rho_g, Z, Cg, Bg, Eg, mi_g