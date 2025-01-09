import sympy as sp



a, chi, u, v, w, = sp.symbols('a chi u v w', real=True, positive=True)
c = sp.Symbol('c', real=True, positive=True)
rho_0 = sp.Symbol('rho_0', real=True, positive=True)
kappa = sp.Symbol('kappa', real=True, positive=True)
P = sp.Symbol('P', real=True, positive=True)
G = sp.Symbol('G', real=True, positive=True)
chi = sp.Symbol('chi', real=True, positive=True)

U = sp.Function('U')(chi)
V = sp.Function('V')(chi)
W = sp.Function('W')(chi)


P_chi = sp.Function('P')(chi)
R_chi = sp.Function('R')(chi)

def e_tau(U, a, chi):
    return sp.Matrix([(1/a)*(sp.exp(-U)), 0, 0, 0])

def e_chi(U, V, a):
    return sp.Matrix([0, (1/a)*(sp.exp(U-V)) / a, 0, 0])

def e_phi(a, chi):
    return sp.Matrix([0, 0, 1/(a*chi), 0])

def e_zeta(W, a):
    return sp.Matrix([0, 0, 0, (1/a)*sp.exp(W)])


etau   = e_tau(u, a, chi)
echi   = e_chi(u, v, a)
ephi = e_phi(a, chi)
ezeta  = e_zeta(w, a)

print(f"[{e_chi}, {e_tau}, {e_phi}, {e_zeta}]")



g1 = sp.diag(-1, 1, 1, 1)

g = sp.Matrix([
    [-1,  0,  0,  0],
    [ 0,  1,  0,  0],
    [ 0,  0,  1,  0],
    [ 0,  0,  0,  1]
])

ig = g.inv()



Vu = c * sp.Matrix([(1/a)*(sp.exp(-U)), 0, 0, 0])
Vd = g * Vu


print(Vu)
print(Vd)

Vd[0] = sp.simplify(-a * Vd[0] * a)

expr = (c * sp.exp(-U)) / a

simplified_expr = sp.simplify(-a * expr * a)


print("Statystyczna predkosc cieczy", simplified_expr)



H_ud = ig * ((1 / c**2) * Vd * Vd.T + g)


H_ud = H_ud.applyfunc(lambda x: sp.simplify(x))
H_ud = H_ud.subs(a**2 * sp.exp(-2 * U), 1)
print("\nH_ud =")
sp.pprint(H_ud)
##do tad jest ok


T_ud = (rho_0 / c**2)* P_chi *H_ud  + ig *(1/ c**2)*(Vd*Vu.T)*P_chi
T_ud = T_ud.applyfunc(lambda x: sp.simplify(x))
T_ud = T_ud.subs(sp.exp(-2 * U), rho_0)

T_ud_simplified = sp.simplify(T_ud)


print("\nTensor energii-pędu T_{ud}:")
sp.pprint(T_ud_simplified)



Eud = (c**2 / 3) * (T_ud_simplified.trace() + (Vd.T * T_ud_simplified * Vd)[0, 0])


Eud = Eud.subs(a**2 * sp.exp(-2 * U), 1)
Eud = Eud.subs(sp.exp(-U) * sp.exp(U), 1)

Eud = Eud.subs((Vu.T * Vd)[0, 0], -c**2)

Eud = sp.simplify(Eud*3/(-c**4 - c**2 +3))

print("\nEud =")
sp.pprint(Eud)




rho_0 = rho_0*c**2



EQud = (Eud * sp.eye(4) - (8 * sp.pi * G / c**4) * T_ud)

EQud = EQud.subs(G, kappa * c**4 / (rho_0 * a**2 * 8 * sp.pi))
EQud = EQud.subs(rho_0, rho_0 * c**2)


EQud_simplified = sp.simplify(EQud)

print("EQud (po uproszczeniu):")
sp.pprint(EQud_simplified)





EQtau =  etau.T * g * EQud * etau
EQchi = echi.T*g*H_ud*EQud*echi
EQphi = ephi.T*g*H_ud*EQud*ephi
EQzeta = ezeta.T*g*H_ud*EQud*ezeta





sp.pprint(EQtau)
sp.pprint(EQchi)        
sp.pprint(EQphi)
sp.pprint(EQzeta)

# # Rozwiązywanie równania EQtau = 0 względem R(x)
# EQtau_simplified = sp.simplify(EQtau[0, 0])
# SOL_tau = sp.solve(EQtau_simplified, R_chi)
# print("Rozwiązanie dla R(x):")
# sp.pprint(SOL_tau)

# # Rozwiązywanie równania EQchi = 0 względem P(x)
# EQchi_simplified = sp.simplify(EQchi[0, 0])
# SOL_chi = sp.solve(EQchi_simplified, P_chi)
# print("Rozwiązanie dla P(x):")
# sp.pprint(SOL_chi)

# SOL_tau = sp.solve(EQtau_simplified, R_chi)
# SOL_chi = sp.solve(EQchi_simplified, P_chi)
# SOL_phi_zeta = sp.solve([EQphi[0, 0], EQzeta[0, 0]], [U.diff(chi, 2), W.diff(chi, 2)])