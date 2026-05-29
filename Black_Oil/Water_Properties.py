import numpy as np
from barril.units import Array, Scalar

def Calculate_water_properties(P, T, S):
    S_value = Scalar('dimensionless', S, '%')
    salinity = Scalar('concentration', S_value.GetValue() * 10000, 'mg/L')

    rho_w_value = 62.368 + 0.438603 * S_value.GetValue('%') + 1.60074e-3 * S**2
    rho_w = Scalar('density', rho_w_value, 'lbm/ft3')

    A = 8.15839 + (-6.12265e-2)*T.GetValues('degF') + (1.91663e-4) * \
        T.GetValues('degF')**2 + (-2.1654e-7)*T.GetValues('degF')**3
    B = 1.01021e-2 + (-7.44241e-5)*T.GetValues('degF') + (3.0553e-7) * \
        T.GetValues('degF')**2 + (-2.94883e-10)*T.GetValues('degF')**3
    C = ((-9.02505) + (0.130237)*T.GetValues('degF') + (-8.53425e-4)*T.GetValues('degF') **
         2 + (2.34122e-6)*T.GetValues('degF')**3 + (-2.37049e-9)*T.GetValues('degF')**4) * 1e-7
    Rsw_values = A + B*P.GetValues('psia') + C*P.GetValues('psia')**2
    Rsw = Array('standard volume per standard volume', Rsw_values, 'scf/stb')

    Cw_values = []
    for i in range(len(P)):
        if (1000 <= P.GetValues('psia')[i] <= 20000) and (200 <= T.GetValues('degF')[i] <= 270) and (0 <= salinity.GetValue('mg/L') <= 200):
            Cw_values.append(1 / (7.033 * P.GetValues('psia')[
                             i] + 0.5415 * S_value.GetValue() - 537.0 * T.GetValues('degF')[i] + 403300))
        else:
            Cw_values.append(0)
    Cw = Array('compressibility', np.array(
        [v for v in Cw_values if v is not None], np.float64), '1/psi')

    VwT = -1.0001e-2 + 1.33391e-4 * \
        T.GetValues('degF') + 5.50654e-7 * T.GetValues('degF')**2
    Vwp = -1.95301e-9 * P.GetValues('psia') * T.GetValues('degF') - 1.72834e-13 * P.GetValues(
        'psia')**2 * T.GetValues('degF') - 3.58922e-7 * P.GetValues('psia') - 2.25341e-10 * P.GetValues('psia')**2
    Bw_values = (1 + VwT) * (1 + Vwp)
    Bw = Array('volume per standard volume', Bw_values, 'bbl/stb')

    A = 109.527 + (-8.40564)*S + (0.313314)*S**2 + (8.72213e-3)*S**3
    B = (-1.12166) + (2.63951e-2)*S + (-6.79461e-4) * \
        S**2 + (-5.47119e-5)*S**3 + (-1.55586e-6)*S**4
    mi_w1_values = A * (T.GetValues('degF') ** B) * (0.9994 + 4.0295e-5 *
                                                     P.GetValues('psia') + 3.1062e-9 * P.GetValues('psia')**2)
    mi_w1 = Array('dynamic viscosity', mi_w1_values, 'cP')

    mi_w_values = mi_w1 * (0.9994 + 4.0295e-5 * P +
                           3.1062e-9 * P.GetValues('psia')**2)
    mi_w = Array('dynamic viscosity', mi_w_values, 'cP')

    mi_w_ratio_values = mi_w.GetValues()/mi_w1.GetValues()
    mi_w_ratio = Array('dimensionless', mi_w_ratio_values, 'unitless')

    return rho_w, Rsw, Cw, Bw, mi_w1, mi_w, mi_w_ratio