import sympy as sp


a, chi, u, v, w = sp.symbols('a chi u v w', real=True, positive=True)
c = sp.symbols('c', real=True, positive=True)


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


Vu = c * etau

scalar_product_Vu_Vu = Vu.dot(Vu)
print("\nIloczyn skalarny Vu ze sobą:", scalar_product_Vu_Vu.simplify())


g = sp.Matrix([[-1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]])


Vd = g * Vu
print("\nKowariantne składowe Vd:", Vd)


scalar_product_Vu_Vd = Vu.dot(Vd)
print("\nIloczyn skalarny Vu i Vd:", scalar_product_Vu_Vd.simplify())
