import numpy as np
from barril.units import Scalar


def Calculate_oil_properties(P, T, do, rho_o, Pb, Rs, dg, Bg):
    rho_water = Scalar('density', 999.1, 'kg/m3')

    if do is None and rho_o is not None:
        do_value = rho_o.GetValue('kg/m3') / rho_water.GetValue('kg/m3')
        do = Scalar('dimensionless', do_value, 'unitless')

    if do is None:
        raise ValueError("É necessário informar do ou rho_o.")

    API_value = 141.5 / do.GetValue() - 131.5
    API = Scalar('dimensionless', API_value, 'unitless')

    rho_ob = Scalar('density', do.GetValue() * rho_water.GetValue('kg/m3'), 'kg/m3')

    if Pb is None and Rs is not None:
        a = 0.00091 * T.GetValue('degF') - 0.0125 * API.GetValue()
        Pb_value = 18.2 * (((Rs.GetValue('scf/stb') / dg.GetValue())**0.83) * (10**a) - 1.4)
        Pb = Scalar('pressure', Pb_value, 'psia')

    if Pb is None:
        raise ValueError("É necessário informar Pb ou Rs para calcular as propriedades do óleo.")

    Calculate_Rsb = dg.GetValue() * (((Pb.GetValue('psia') / 18.2) + 1.4) * 10**(0.0125 * API.GetValue() - 0.00091 * T.GetValue('degF')))**(1 / 0.83)
    Rsb = Scalar('standard volume per standard volume', Calculate_Rsb, 'scf/stb')

    if P.GetValue('psia') <= Pb.GetValue('psia'):
        Calculate_Rs = dg.GetValue() * (((P.GetValue('psia') / 18.2) + 1.4) * 10**(0.0125 * API.GetValue() - 0.00091 * T.GetValue('degF')))**(1 / 0.83)
    else:
        Calculate_Rs = Rsb.GetValue('scf/stb')

    Rs = Scalar('standard volume per standard volume', Calculate_Rs, 'scf/stb')

    if P.GetValue('psia') < Pb.GetValue('psia'):
        Bo_valuess = 0.9759 + 0.00012 * ((Rs.GetValue('scf/stb') * (dg.GetValue() / do.GetValue())**0.5 + 1.25 * T.GetValue('degF'))**1.2)

        dRs_dP = dg.GetValue() * (1 / 0.83) * (1 / 18.2) * np.power(10**(0.0125 * API.GetValue() - 0.00091 * T.GetValue('degF')), 1 / 0.83) * np.power((P.GetValue('psia') / 18.2 + 1.4), (1 / 0.83 - 1))

        dBo_dP = 0.00012 * 1.2 * (dg.GetValue() / do.GetValue())**0.5 * ((Rs.GetValue('scf/stb') * (dg.GetValue() / do.GetValue())**0.5 + 1.25 * T.GetValue('degF'))**0.2) * dRs_dP

        Calculate_Co = (1 / Bo_valuess) * dBo_dP + (Bg.GetValue('bbl/stb') / Bo_valuess) * dRs_dP
    else:
        Calculate_Co = 1e-6 * np.exp((rho_ob.GetValue('lbm/ft3') + 0.004347 * (P.GetValue('psia') - Pb.GetValue('psia')) - 79.1) / (0.0007141 * (P.GetValue('psia') - Pb.GetValue('psia')) - 12.938))

    Co = Scalar('compressibility', Calculate_Co, '1/psi')

    Calculate_Bob = 0.9759 + 0.00012 * ((Rsb.GetValue('scf/stb') * (dg.GetValue() / do.GetValue())**0.5 + 1.25 * T.GetValue('degF'))**1.2)
    Bob = Scalar('standard volume per standard volume', Calculate_Bob, 'stb/stb')

    if P.GetValue('psia') > Pb.GetValue('psia'):
        Calculate_Bo = Bob.GetValue('stb/stb') * np.exp(-Co.GetValue('1/psi') * (P.GetValue('psia') - Pb.GetValue('psia')))
    else:
        Calculate_Bo = 0.9759 + 0.00012 * ((Rs.GetValue('scf/stb') * (dg.GetValue() / do.GetValue())**0.5 + 1.25 * T.GetValue('degF'))**1.2)

    Bo = Scalar('volume per standard volume', Calculate_Bo, 'bbl/stb')

    if P.GetValue('psia') >= Pb.GetValue('psia'):
        Calculate_Bt = Bo.GetValue('bbl/stb')
    else:
        Calculate_Bt = Bo.GetValue('bbl/stb') + (Rsb.GetValue('scf/stb') - Rs.GetValue('scf/stb')) * Bg.GetValue('bbl/stb')

    Bt = Scalar('volume per standard volume', Calculate_Bt, 'bbl/stb')

    if P.GetValue('psia') <= Pb.GetValue('psia'):
        Calculate_rho_o = (62.4 * do.GetValue() + 0.0136 * Rs.GetValue('scf/stb') * dg.GetValue()) / Bo.GetValue('bbl/stb')
    else:
        Calculate_rho_o = rho_ob.GetValue('lbm/ft3') * np.exp(Co.GetValue('1/psi') * (P.GetValue('psia') - Pb.GetValue('psia')))

    rho_o = Scalar('density', Calculate_rho_o, 'lbm/ft3')

    A = 10**(0.43 + (8.33 / API.GetValue()))

    Calculate_mi_od = 0.32 + (1.8e7 / (API.GetValue()**4.53)) * ((360 / (T.GetValue('degR') - 260))**A)
    mi_od = Scalar('dynamic viscosity', Calculate_mi_od, 'cP')

    a = 10**(-7.4e-4 * Rsb.GetValue('scf/stb') + 2.2e-7 * Rsb.GetValue('scf/stb')**2)

    b = 0.68 / (10**(8.62e-5 * Rsb.GetValue('scf/stb'))) + 0.25 / (10**(1.1e-3 * Rsb.GetValue('scf/stb'))) + 0.062 / (10**(3.74e-3 * Rsb.GetValue('scf/stb')))

    Calculate_mi_ob = a * (mi_od.GetValue('cP')**b)
    mi_ob = Scalar('dynamic viscosity', Calculate_mi_ob, 'cP')

    if P.GetValue('psia') <= Pb.GetValue('psia'):
        Calculate_mi_o = mi_ob.GetValue('cP')
    else:
        Calculate_mi_o = mi_ob.GetValue('cP') + (0.001 * (P.GetValue('psia') - Pb.GetValue('psia'))) * (0.024 * mi_ob.GetValue('cP')**1.6 + 0.038 * mi_ob.GetValue('cP')**0.56)

    mi_o = Scalar('dynamic viscosity', Calculate_mi_o, 'cP')

    return do, API, Pb, Rs, Co, Bo, Bt, rho_o, mi_od, mi_ob, mi_o