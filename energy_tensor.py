import sympy as sp

a, chi, u, v, w = sp.symbols('a chi u v w', real=True, positive=True)
c = sp.symbols('c', real=True, positive=True)
rho_0 = sp.symbols('rho_0', real=True, positive=True)
P = sp.symbols('P', real=True, positive=True)
G = sp.Symbol('G', real=True, positive=True)
kappa = sp.Symbol('kappa', real=True, positive=True)
U = sp.Function('U')(chi)
V = sp.Function('V')(chi)
W = sp.Function('W')(chi)


def e_tau(u, a):
    return sp.Matrix([sp.exp(-u) / a, 0, 0, 0])

def e_x(u, v, a):
    return sp.Matrix([0, sp.exp(u - v) / a, 0, 0])

def e_phi(a, chi):
    return sp.Matrix([0, 0, 1 / (a * chi), 0])

def e_xi(w, a):
    return sp.Matrix([0, 0, 0, sp.exp(w) / a])

etau = e_tau(u, a)
ex = e_x(u, v, a)
ephi = e_phi(a, chi)
exi = e_xi(w, a)

print("Baza ortonormalna:")
print("e_tau:", etau)
print("e_x:", ex)
print("e_phi:", ephi)
print("e_xi:", exi)


Vu = sp.Matrix([c,0,0,0])

scalar_product_Vu_Vu = Vu.dot(Vu)
print("\nIloczyn skalarny Vu ze sobą:", scalar_product_Vu_Vu.simplify())


g = sp.Matrix([[-1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]])

matrix_trace = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

Vd = g * Vu
print("\nKowariantne składowe Vd:", Vd)


scalar_product_Vu_Vd = Vu.dot(Vd)
print("\nIloczyn skalarny Vu i Vd:", scalar_product_Vu_Vd.simplify())

H_ud =  g +(1 / c**2) * Vd * Vd.T

print("\nH_ud:",H_ud)

T_ud = (rho_0 * P / c**2) * H_ud + (rho_0 * P / c**4) * Vd * Vu.T


T_ud_simplified = T_ud.expand()
print("\nTensor energii-pędu T^{ud}:")
sp.pprint(T_ud_simplified)

Vd_Tud_Vu = Vd.T * T_ud * Vu
print("\nIloczyn Vd * T^{ud} * Vu:")
sp.pprint(Vd_Tud_Vu[0])

rho_0 = rho_0 * c**2
G = (kappa * c**4) / (8 * sp.pi * a**2 * rho_0)

EQ = g - (8 * sp.pi * G / c**4) * T_ud

print(EQ)



