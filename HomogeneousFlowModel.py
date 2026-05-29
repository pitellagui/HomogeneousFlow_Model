import numpy as np
from barril.units import Scalar


def Calculate_HomogeneousFlowProperties(
    Bg,
    Rs,
    Rsw,
    Vo_sc,
    Bo,
    Vw_sc,
    Bw,
    rho_o,
    rho_w,
    mu_o,
    mu_w,
    Vg_sc,
    rho_G,
    mu_G,
    D_h,
    epsilon
):
    if (Vo_sc.GetValue("m3/s") + Vw_sc.GetValue("m3/s")) == 0:
        raise ValueError("A soma de Vo_sc e Vw_sc não pode ser zero.")

    Bsw_value = Vw_sc.GetValue("m3/s") / (Vo_sc.GetValue("m3/s") + Vw_sc.GetValue("m3/s"))

    GOR_value = Vg_sc.GetValue("m3/s") / Vo_sc.GetValue("m3/s")

    if Bsw_value == 0:
        Vg_value = (GOR_value - Rs.GetValue("sm3/sm3")) * Vo_sc.GetValue("m3/s") * Bg.GetValue("m3/sm3")
    else:
        Vg_value = (Vg_sc.GetValue("m3/s") - Vo_sc.GetValue("m3/s") * Rs.GetValue("sm3/sm3") - Vw_sc.GetValue("m3/s") * Rsw.GetValue("sm3/sm3")) * Bg.GetValue("m3/sm3")

    Vl_value = Vo_sc.GetValue("m3/s") * Bo.GetValue("bbl/stb") + Vw_sc.GetValue("m3/s") * Bw.GetValue("bbl/stb")

    A_p_value = np.pi * D_h.GetValue("m")**2 / 4

    Vt_value = Vl_value + Vg_value

    v_m_value = Vt_value / A_p_value

    lambda_L_value = Vl_value / Vt_value

    lambda_G_value = 1 - lambda_L_value

    rho_L_value = rho_o.GetValue("kg/m3") * (Vo_sc.GetValue("m3/s") * Bo.GetValue("bbl/stb")) / Vl_value + rho_w.GetValue("kg/m3") * (Vw_sc.GetValue("m3/s") * Bw.GetValue("bbl/stb")) / Vl_value

    rho_m_value = lambda_L_value * rho_L_value + lambda_G_value * rho_G.GetValue("kg/m3")

    mu_L_value = mu_o.GetValue("Pa.s") * (Vo_sc.GetValue("m3/s") * Bo.GetValue("bbl/stb")) / Vl_value + mu_w.GetValue("Pa.s") * (Vw_sc.GetValue("m3/s") * Bw.GetValue("bbl/stb")) / Vl_value

    mu_m_value = lambda_L_value * mu_L_value + lambda_G_value * mu_G.GetValue("Pa.s")

    m_L_value = rho_L_value * Vl_value

    m_G_value = rho_G.GetValue("kg/m3") * Vg_value

    m_M_value = m_L_value + m_G_value

    x_value = m_G_value / m_M_value

    Re_m_value = rho_m_value * v_m_value * D_h.GetValue("m") / mu_m_value

    f_TP_value = 0.0055 * (1 + (20000 * epsilon.GetValue("m") / D_h.GetValue("m") + 1000000 / Re_m_value)**(1 / 3))

    Bsw = Scalar("dimensionless", Bsw_value, "unitless")
    GOR = Scalar("dimensionless", GOR_value, "unitless")
    Vg = Scalar("volume flow rate", Vg_value, "m3/s")
    Vl = Scalar("volume flow rate", Vl_value, "m3/s")
    Vt = Scalar("volume flow rate", Vt_value, "m3/s")
    A_p = Scalar("area", A_p_value, "m2")
    v_m = Scalar("velocity", v_m_value, "m/s")
    lambda_L = Scalar("dimensionless", lambda_L_value, "unitless")
    lambda_G = Scalar("dimensionless", lambda_G_value, "unitless")
    rho_L = Scalar("density", rho_L_value, "kg/m3")
    rho_m = Scalar("density", rho_m_value, "kg/m3")
    mu_L = Scalar("dynamic viscosity", mu_L_value, "Pa.s")
    mu_m = Scalar("dynamic viscosity", mu_m_value, "Pa.s")
    m_L = Scalar("mass flow rate", m_L_value, "kg/s")
    m_G = Scalar("mass flow rate", m_G_value, "kg/s")
    m_M = Scalar("mass flow rate", m_M_value, "kg/s")
    x = Scalar("dimensionless", x_value, "unitless")
    Re_m = Scalar("dimensionless", Re_m_value, "unitless")
    f_TP = Scalar("dimensionless", f_TP_value, "unitless")

    return {
        "Bsw": Bsw,
        "GOR": GOR,
        "Q_G_calc": Vg,
        "Q_L": Vl,
        "Q_G": Vg,
        "Q_total": Vt,
        "A_p": A_p,
        "v_m": v_m,
        "lambda_L": lambda_L,
        "lambda_G": lambda_G,
        "rho_L": rho_L,
        "rho_m": rho_m,
        "mu_L": mu_L,
        "mu_m": mu_m,
        "m_dot_L": m_L,
        "m_dot_G": m_G,
        "m_dot_m": m_M,
        "x": x,
        "Re_m": Re_m,
        "f_TP": f_TP
    }

def Calculate_dpdl_t(rho_m, f_TP, m_dot_m, x, v_m, D_h, theta, A_p, M, T, rho_G):
    g_value = 9.81
    R = Scalar("molar gas constant", 8.314, "J/mol.K")

    dpdl_f_value = -f_TP.GetValue() * rho_m.GetValue("kg/m3") * v_m.GetValue("m/s")**2 / (2 * D_h.GetValue("m"))

    dpdl_g_value = -rho_m.GetValue("kg/m3") * g_value * np.sin(theta.GetValue("rad"))

    aceleration_cte_value = -((m_dot_m.GetValue("kg/s") / A_p.GetValue("m2"))**2) * x.GetValue() * M.GetValue("kg/mol") / (R.GetValue("J/mol.K") * T.GetValue("K") * rho_G.GetValue("kg/m3")**2)

    dpdl_t_value = (dpdl_f_value + dpdl_g_value) / (1 - aceleration_cte_value)

    dpdl_a_value = aceleration_cte_value * dpdl_t_value

    dpdl_f = Scalar("pressure per length", dpdl_f_value, "Pa/m")
    dpdl_g = Scalar("pressure per length", dpdl_g_value, "Pa/m")
    aceleration_cte = Scalar("dimensionless", aceleration_cte_value, "unitless")
    dpdl_a = Scalar("pressure per length", dpdl_a_value, "Pa/m")
    dpdl_t = Scalar("pressure per length", dpdl_t_value, "Pa/m")

    return {
        "dpdl_friccao": dpdl_f,
        "dpdl_gravidade": dpdl_g,
        "valor_aceleracao": aceleration_cte,
        "dpdl_aceleracao": dpdl_a,
        "dpdl_total": dpdl_t
    }
