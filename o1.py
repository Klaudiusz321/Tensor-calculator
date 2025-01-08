import sympy as sp



a, chi, u, v, w = sp.symbols('a chi u v w', real=True, positive=True)
c = sp.Symbol('c', real=True, positive=True)
rho_0 = sp.Symbol('rho_0', real=True, positive=True)
kappa = sp.Symbol('kappa', real=True, positive=True)

G = sp.Symbol('G', real=True, positive=True)


U = sp.Function('U')(chi)
V = sp.Function('V')(chi)
W = sp.Function('W')(chi)


P_chi = sp.Function('P')(chi)
R_chi = sp.Function('R')(chi)

def e_tau(u, a):
    return sp.Matrix([sp.exp(-u) / a, 0, 0, 0])

def e_chi(u, v, a):
    return sp.Matrix([0, sp.exp(u - v) / a, 0, 0])

def e_theta(a, chi):
    return sp.Matrix([0, 0, 1/(a*chi), 0])

def e_zeta(w, a):
    return sp.Matrix([0, 0, 0, sp.exp(w)/a])



etau   = e_tau(u, a)
echi   = e_chi(u, v, a)
etheta = e_theta(a, chi)
ezeta  = e_zeta(w, a)

print(f"[{e_chi}, {e_tau}, {e_theta}, {e_zeta}]")

Vu = sp.Matrix([c, 0, 0, 0])
scalar_product_Vu_Vu = Vu.dot(Vu)
print("\nIloczyn skalarny Vu ze sobą =", scalar_product_Vu_Vu.simplify())

g = sp.Matrix([
    [-1,  0,  0,  0],
    [ 0,  1,  0,  0],
    [ 0,  0,  1,  0],
    [ 0,  0,  0,  1]
])

Vd = g * Vu
print("\nKowariantne składowe Vd =", Vd)

scalar_product_Vu_Vd = Vu.dot(Vd)
print("\nIloczyn skalarny Vu i Vd =", scalar_product_Vu_Vd.simplify())

H_ud = g + (1/c**2)*Vd*Vd.T
print("\nH_ud =")
sp.pprint(H_ud)


T_ud = (rho_0 * P_chi / c**2)*H_ud  +  (rho_0 * P_chi / c**4)*(Vd*Vu.T)
T_ud_simplified = T_ud.expand()

print("\nTensor energii-pędu T_{ud}:")
sp.pprint(T_ud_simplified)

Vd_Tud_Vu = Vd.T * T_ud * Vu
print("\nIloczyn Vd^T * T_ud * Vu =")
sp.pprint(Vd_Tud_Vu[0])


rho_0_eff = rho_0 * c**2

G_eff = (kappa * c**4) / (8*sp.pi * a**2 * rho_0_eff)
G_ud = g - (8*sp.pi * G_eff / c**4)*T_ud  

Gud_simplified = G_ud.expand()
EQud = sp.MatrixSymbol('EQud', 4, 4)
print("\n\"G_ud\" = g - (8 pi G_eff / c^4)*T_ud =")
sp.pprint(G_ud)

G_up = G_ud * g.inv()
print("\n\"G_up\" = G_ud * g.inv() =")
sp.pprint(G_up)

EQtau = sp.MatrixSymbol('EQtau', 4, 4)
EQchi = sp.MatrixSymbol('EQchi', 4, 4)
EQtheta = sp.MatrixSymbol('EQtheta', 4, 4)
EQzeta = sp.MatrixSymbol('EQzeta', 4, 4)

