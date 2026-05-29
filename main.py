from Black_Oil.Water_Properties import Calculate_water_properties
from Black_Oil.Gas_Properties import Calculate_gas_properties
from Black_Oil.Oil_Properties import Calculate_oil_properties
from Homogeneous_Flow_Model import Calculate_homogeneousFlow_properties, Calculate_dpdl_t

import numpy as np 
from barril.units import Array, Scalar

P = Array('pressure', np.array([   0, 1200, 2400, 3600, 4800, 6000],np.float64), "psia")
T = Array('temperature', np.array([ 70, 100, 130, 160, 190, 220],np.float64), "degF")

Pb = Array('pressure', np.array([2800, 2792, 2784, 2776, 2768, 2760], np.float64), "psia")

dg = Scalar('dimensionless', 0.72, "unitless")
do = Scalar('dimensionless', 0.81, "unitless")
S = 0.01

if len(Pb) < len(P):
  Pb_value = np.pad(Pb, (0, len(P) - len(Pb)), mode='edge')
  Pb = Array('pressure', Pb_value, "psia")
rho_o = None
Rs = None
Mg, dg, P_pc, T_pc, P_pr, T_pr, rho_g, Z, Cg, Bg, Eg, mi_g = Calculate_gas_properties(P,T,dg, rho_g=None)
do, API, Pb, Rs, Co, Bo, Bt, rho_o, mi_od, mi_ob, mi_o = Calculate_oil_properties(P, T, do, rho_o, Pb, Rs,dg,Bg)
rho_w, Rsw, Cw, Bw, mi_w1, mi_w, mi_w_ratio = Calculate_water_properties(P, T, S)
#Calculate_homogeneousFlow_properties
#Calculate_dpdl_t

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
print("mi_g =", mi_g)

print("\n=== WATER ===")
print("rho_w =", rho_w)
print("Rsw =", Rsw)
print("Cw =", Cw)
print("Bw =", Bw)
print("mi_w1 =", mi_w1)
print("mi_w =", mi_w)
print("mi_w_ratio =", mi_w_ratio)

print("\n=== OIL ===")
print("do =", do)
print("API =", API)
print("Pb =", Pb)
print("Rs =", Rs)
print("Co =", Co)
print("Bo =", Bo)
print("Bt =", Bt)
print("rho_o =", rho_o)
print("mi_od =", mi_od)
print("mi_ob =", mi_ob)
print("mi_o =", mi_o)