from Black_Oil.Water_Properties import Calculate_water_properties
from Black_Oil.Gas_Properties import Calculate_gas_properties
from Black_Oil.Oil_Properties import Calculate_oil_properties
from Homogeneous_Flow_Model import Calculate_homogeneousFlow_properties, Calculate_dpdl_t
from barril.units import Scalar


P = Scalar("pressure", 1200, "psia")
T = Scalar("temperature", 100, "degF")
Pb = Scalar("pressure", 2792, "psia")
dg = Scalar("dimensionless", 0.72, "unitless")
do = Scalar("dimensionless", 0.81, "unitless")
S = 0.01
rho_g_input = None
rho_o_input = None
Rs_input = None
D_h = Scalar("length", 0.05, "m")
epsilon = Scalar("length", 0.000045, "m")
theta = Scalar("plane angle", 0, "rad")
Vo_sc = Scalar("volume flow rate", 0.01, "m3/s")
Vw_sc = Scalar("volume flow rate", 0.0, "m3/s")
Vg_sc = Scalar("volume flow rate", 1.0, "m3/s")

Mg, dg, P_pc, T_pc, P_pr, T_pr, rho_g, Z, Cg, Bg, Eg, mu_g = Calculate_gas_properties(P, T, dg, rho_g_input)

do, API, Pb, Rs, Co, Bo, Bt, rho_o, mu_od, mu_ob, mu_o = Calculate_oil_properties(P, T, do, rho_o_input, Pb, Rs_input, dg, Bg)

rho_w, Rsw, Cw, Bw, mu_w1, mu_w, mu_w_ratio = Calculate_water_properties(P, T, S)

flow_properties = Calculate_homogeneousFlow_properties(
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
    rho_g,
    mu_g,
    D_h,
    epsilon
)

dpdl_results = Calculate_dpdl_t(
    flow_properties["rho_m"],
    flow_properties["f_TP"],
    flow_properties["m_dot_m"],
    flow_properties["x"],
    flow_properties["v_m"],
    D_h,
    theta,
    flow_properties["A_p"],
    Mg,
    T,
    rho_g
)

print("=== GAS ===")
print("Mg =", Mg)
print("dg =", dg)
print("P_pc =", P_pc)
print("T_pc =", T_pc)
print("P_pr =", P_pr)
print("T_pr =", T_pr)
print("rho_g =", rho_g)
print("Z =", Z)
print("Cg =", Cg)
print("Bg =", Bg)
print("Eg =", Eg)
print("mi_g =", mu_g)


print("\n=== WATER ===")
print("rho_w =", rho_w)
print("Rsw =", Rsw)
print("Cw =", Cw)
print("Bw =", Bw)
print("mi_w1 =", mu_w1)
print("mi_w =", mu_w)
print("mi_w_ratio =", mu_w_ratio)


print("\n=== OIL ===")
print("do =", do)
print("API =", API)
print("Pb =", Pb)
print("Rs =", Rs)
print("Co =", Co)
print("Bo =", Bo)
print("Bt =", Bt)
print("rho_o =", rho_o)
print("mi_od =", mu_od)
print("mi_ob =", mu_ob)
print("mi_o =", mu_o)


print("\n=== HOMOGENEOUS FLOW ===")
print("GOR =", flow_properties["GOR"])
print("Q_G_calc =", flow_properties["Q_G_calc"])
print("Q_L =", flow_properties["Q_L"])
print("Q_G =", flow_properties["Q_G"])
print("Q_total =", flow_properties["Q_total"])
print("A_p =", flow_properties["A_p"])
print("v_m =", flow_properties["v_m"])
print("lambda_L =", flow_properties["lambda_L"])
print("lambda_G =", flow_properties["lambda_G"])
print("rho_L =", flow_properties["rho_L"])
print("rho_m =", flow_properties["rho_m"])
print("mu_L =", flow_properties["mu_L"])
print("mu_m =", flow_properties["mu_m"])
print("m_dot_L =", flow_properties["m_dot_L"])
print("m_dot_G =", flow_properties["m_dot_G"])
print("m_dot_m =", flow_properties["m_dot_m"])
print("x =", flow_properties["x"])
print("Re_m =", flow_properties["Re_m"])
print("f_TP =", flow_properties["f_TP"])


print("\n=== PRESSURE GRADIENT ===")
print("dpdl_friccao =", dpdl_results["dpdl_friccao"])
print("dpdl_gravidade =", dpdl_results["dpdl_gravidade"])
print("valor_aceleracao =", dpdl_results["valor_aceleracao"])
print("dpdl_aceleracao =", dpdl_results["dpdl_aceleracao"])
print("dpdl_total =", dpdl_results["dpdl_total"])