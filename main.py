from BlackOilModel import GasPhase, OilPhase, WaterPhase
import unitsconverter
from HomogeneousFlowModel import Calculate_HomogeneousFlowProperties, Calculate_dpdl_t
from barril.units import Scalar

# Dados iniciais
densidaderelativagas = 0.84
densidaderelativaoleo = 0.86
Pressaodebolha = 5000 # psia
pressaolinha = 3626 # psia
temperaturalinha = 122 # °F

gas = GasPhase(zcorrelation="Papay", μcorrelation='Lee', dg=densidaderelativagas, P=pressaolinha, T=unitsconverter.Temperature(temperaturalinha, 'F','R'))
oil = OilPhase("Standing", "Standing", "Beggs",densidaderelativaoleo, densidaderelativagas, pressaolinha, Pressaodebolha, temperaturalinha)
gas.output()
oil.output()

Bg = Scalar(gas.Bg(), "m3/sm3")
rho_G = Scalar(gas.ρ, "lbm/ft3")
mu_G = Scalar(gas.μ, "cP")
Mg = Scalar(gas.Mg, "lb/lbmol")
Rs = Scalar(oil.Rs, "sm3/sm3")
Bo = Scalar(oil.Bo, "bbl/stb")
rho_o = Scalar(oil.ρo, "lbm/ft3")
mu_o = Scalar(oil.μob, "cP")

P = Scalar(pressaolinha, "psia")
T = Scalar(temperaturalinha, "degF")
dg = Scalar(densidaderelativagas, "unitless") 
do = Scalar(densidaderelativaoleo, "unitless")

rho_w, Rsw, Cw, Bw, mi_w1, mi_w, mi_w_ratio = WaterPhase(P, T, 0)

D_h = Scalar("length", 0.05, "m")
epsilon = Scalar("length", 0.000045, "m")
theta = Scalar("plane angle", 0, "rad")
Vo_sc = Scalar("volume flow rate", 0.01, "m3/s")
Vw_sc = Scalar("volume flow rate", 0.0, "m3/s")
Vg_sc = Scalar("volume flow rate", 1.0, "m3/s")

homogeneous_output = Calculate_HomogeneousFlowProperties(Bg, Rs, Rsw, Vo_sc, Bo, Vw_sc, Bw, rho_o, rho_w, mu_o, mi_w, Vg_sc, rho_G, mu_G, D_h, epsilon)
Bsw, GOR, Q_G_calc, Q_L, Q_G, Q_total, A_p, v_m, lambda_L, lambda_G, rho_L, rho_m, mu_L, mu_m, m_dot_L, m_dot_G, m_dot_m, x, Re_m, f_TP = homogeneous_output["Bsw"], homogeneous_output["GOR"], homogeneous_output["Q_G_calc"], homogeneous_output["Q_L"], homogeneous_output["Q_G"], homogeneous_output["Q_total"], homogeneous_output["A_p"], homogeneous_output["v_m"], homogeneous_output["lambda_L"], homogeneous_output["lambda_G"], homogeneous_output["rho_L"], homogeneous_output["rho_m"], homogeneous_output["mu_L"], homogeneous_output["mu_m"], homogeneous_output["m_dot_L"], homogeneous_output["m_dot_G"], homogeneous_output["m_dot_m"], homogeneous_output["x"], homogeneous_output["Re_m"], homogeneous_output["f_TP"]

dpdl_output = Calculate_dpdl_t(rho_m, f_TP, m_dot_m, x, v_m, D_h, theta, A_p, Mg, T, rho_G)
dpdl_friccao, dpdl_gravidade, valor_aceleracao, dpdl_aceleracao, dpdl_total = dpdl_output["dpdl_friccao"], dpdl_output["dpdl_gravidade"], dpdl_output["valor_aceleracao"], dpdl_output["dpdl_aceleracao"], dpdl_output["dpdl_total"]
print(f"Perda de carga por fricção: {dpdl_friccao.GetValue('Pa/m'):.2f} Pa/m"
      f"\nPerda de carga por gravidade: {dpdl_gravidade.GetValue('Pa/m'):.2f} Pa/m"
      f"\nValor da aceleração: {valor_aceleracao.GetValue()}"
      f"\nPerda de carga por aceleração: {dpdl_aceleracao.GetValue('Pa/m'):.2f} Pa/m"
      f"\nPerda de carga total: {dpdl_total.GetValue('Pa/m'):.2f} Pa/m")
