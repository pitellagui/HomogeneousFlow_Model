import numpy as np
from barril.units import Array, Scalar

def Calculate_gas_properties(P, T, dg, rho_g):
    if dg is None and rho_g is not None:
        rho_ar = Scalar('density', 1.225, 'kg/m3')
        rho_g = Scalar('density', rho_g, 'kg/m3')
        dg = Scalar('dimensionless', rho_g.GetValue(
            'kg/m3')/rho_ar.GetValue('kg/m3'), 'unitless')

    Mar = 28.96
    Mg_value = Mar * dg.GetValue()
    Mg = Scalar('mass per mol', Mg_value, 'lb/lbmol')

    if float(dg.GetValue()) < 0.75:
        P_pc_val = 677 + 15.0 * dg.GetValue() - 37.5 * dg.GetValue()**2
        T_pc_val = 168 + 325 * dg.GetValue() - 12.5 * dg.GetValue()**2
    elif float(dg.GetValue()) >= 0.75:
        P_pc_val = 706 - 51.7 * dg.GetValue() - 11.1 * dg.GetValue()**2
        T_pc_val = 187 + 330 * dg.GetValue() - 71.5 * dg.GetValue()**2

    P_pc_value = Scalar('pressure', P_pc_val, 'psia')
    T_pc_value = Scalar('temperature', T_pc_val, 'degR')
    P_pc = Scalar('pressure', P_pc_value.GetValue('Pa'), 'Pa')
    T_pc = Scalar('temperature', T_pc_value.GetValue('K'), 'K')
    P_pr = Array('dimensionless', P.GetValues('Pa') /
                 P_pc_value.GetValue('Pa'), 'unitless')
    T_pr = Array('dimensionless', T.GetValues('K') /
                 T_pc_value.GetValue('K'), 'unitless')

    R = Scalar('molar gas constant', 8.314, 'J/mol.K')
    Z_values = 1 - (3.53*P_pr.GetValues())/(10**(0.9813*T_pr.GetValues())) + \
        (0.274*P_pr.GetValues()**2)/(10**(0.8157*T_pr.GetValues()))
    Z = Array('dimensionless', Z_values, 'unitless')
    rho_g_values = (P.GetValues('Pa') * Mg.GetValue('kg/mol')) / \
        (Z * R.GetValue('J/mol.K') * T.GetValues('K'))
    rho_g = Array('density', rho_g_values, 'kg/m3')

    dZdP_pr = - (3.53 / (10**(0.9813 * T_pr.GetValues()))) + (2 *
                                                              0.274 * P_pr.GetValues() / (10**(0.8157 * T_pr.GetValues())))
    C_pr_values = (1 / P_pr.GetValues()) - ((1 / Z) * dZdP_pr)
    C_pr = Array('dimensionless', C_pr_values, 'unitless')

    Cg_values = C_pr.GetValues()/P_pc.GetValue('Pa')
    Cg = Array('compressibility', Cg_values, '1/Pa')

    Psc = Scalar('pressure', 14.7, 'psia')
    Tsc = Scalar('temperature', 60, 'degF')

    Bg_values = (Psc.GetValue('psia')/Tsc.GetValue('degF')) * \
        (Z.GetValues()*T.GetValues('degF')/P.GetValues('psia'))
    Bg = Array('volume per standard volume', Bg_values, 'm3/sm3')

    Eg_values = 1/Bg.GetValues('m3/sm3')
    Eg = Array('standard volume per volume', Eg_values, 'sm3/m3')

    xv = 3.448 + 986.4 / T.GetValues('degR') + \
        0.01009 * Mg.GetValue('lb/lbmol')
    yv = 2.4 - 0.2 * xv
    kv = ((9.379 + 0.0160 * Mg.GetValue('lb/lbmol')) * (T.GetValues('degR') **
          1.5)) / (209.2 + 19.26 * Mg.GetValue('lb/lbmol') + T.GetValues('degR'))
    mi_g_values = 1e-4 * kv * \
        np.exp(xv * (rho_g.GetValues('lbm/ft3') / 62.4) ** yv)
    mi_g = Array('dynamic viscosity', mi_g_values, 'cP')

    return Mg, dg, P_pc, T_pc, P_pr, T_pr, rho_g, Z, Cg, Bg, Eg, mi_g