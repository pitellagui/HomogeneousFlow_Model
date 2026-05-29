import numpy as np
from barril.units import Array, Scalar

def Calculate_homogeneousFlow_properties(Bg, Bsw, Rs, Rsw, Vo_sc, Bo, Vw_sc, Bw, rho_o, rho_w, mu_o, mu_w, Vg, rho_G, mu_G, D_h, epsilon):

    GOR = Vg_sc / Vo_sc
    if Bsw == 0:
        Vg_sc = (GOR - Rs)*Vo_sc*Bg
    else:
        Vg_sc = (Vg_sc - Vo_sc * Rs - Vw_sc * Rsw) * Bg

    Vl = Vo_sc * Bo + Vw_sc * Bw

    A_p = np.pi * D_h.GetValue("m")**2 / 4

    Vt = Vl.GetValue("m3/s") + Vg.GetValue("m3/s")

    v_m = Vt / A_p

    lambda_L = Vl.GetValue("m3/s") / Vt

    lambda_G = 1 - lambda_L
    
    rho_L = rho_o * (Vo_sc * Bo)/Vl + rho_w * (Vw_sc * Bw)/Vl

    rho_m = lambda_L * rho_L.GetValue("kg/m3") + lambda_G * rho_G.GetValue("kg/m3")

    mu_L = mu_o * (Vo_sc * Bo)/Vl + mu_w * (Vw_sc * Bw)/Vl

    mu_m = lambda_L * mu_L.GetValue("Pa.s") + lambda_G * mu_G.GetValue("Pa.s")

    m_L = rho_L.GetValue("kg/m3") * Vl.GetValue("m3/s")

    m_G = rho_G.GetValue("kg/m3") * Vg.GetValue("m3/s")

    m_M = m_L + m_G

    x = m_G / m_M

    Re_m = rho_m * v_m * D_h.GetValue("m") / mu_m

    f_TP = 0.0055 * (1 + (20000 * epsilon.GetValue("m") / D_h.GetValue("m") + 1000000 / Re_m)**(1 / 3))

    return {
        "Q_L": Vl,
        "Q_G": Vg,
        "Q_total": Scalar(Vt, "m3/s"),
        "A_p": Scalar(A_p, "m2"),
        "v_m": Scalar(v_m, "m/s"),
        "lambda_L": lambda_L,
        "lambda_G": lambda_G,
        "rho_m": Scalar(rho_m, "kg/m3"),
        "mu_m": Scalar(mu_m, "Pa.s"),
        "m_dot_L": Scalar(m_L, "kg/s"),
        "m_dot_G": Scalar(m_G, "kg/s"),
        "m_dot_m": Scalar(m_M, "kg/s"),
        "x": x,
        "Re_m": Re_m,
        "f_TP": f_TP
    }

def Calculate_dpdl_t(rho_m, f_TP, m_dot_m, x, v_m, D_h, theta, A_p, M, R, T, rho_G):
    g = Scalar(9.81, "m/s2")
    dpdl_f = -f_TP * rho_m.GetValue("kg/m3") * v_m.GetValue("m/s")**2 / (2 * D_h.GetValue("m"))

    dpdl_g = -rho_m.GetValue("kg/m3") * g.GetValue("m/s2") * np.sin(theta.GetValue("rad"))

    aceleration_cte = -((m_dot_m.GetValue("kg/s") / A_p.GetValue("m2"))**2) * x * M.GetValue("kg/mol") / (
        R.GetValue("J/mol.K") * T.GetValue("K") * rho_G.GetValue("kg/m3")**2 )

    dpdl_t = (dpdl_f + dpdl_g) / (1 - aceleration_cte)

    dpdl_a = aceleration_cte * dpdl_t

    return {
        "dpdl_friccao": Scalar(dpdl_f, "Pa/m"),
        "dpdl_gravidade": Scalar(dpdl_g, "Pa/m"),
        "valor_aceleracao": aceleration_cte,
        "dpdl_aceleracao": Scalar(dpdl_a, "Pa/m"),
        "dpdl_total": Scalar(dpdl_t, "Pa/m")
    }