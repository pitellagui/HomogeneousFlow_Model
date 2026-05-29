import numpy as np
from barril.units import Array, Scalar

def Calculate_water_properties(P, T, S):
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