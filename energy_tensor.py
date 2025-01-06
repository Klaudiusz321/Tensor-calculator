import sympy as sp

a, c = sp.symbols('a c')
tau, chi, theta, zeta = sp.symbols('tau chi theta zeta')
U = sp.Function('U')(chi)
V = sp.Function('V')(chi)
W = sp.Function('W')(chi)
chi_var = sp.symbols('chi')

e_tau = sp.Matrix([a * sp.exp(-U), 0, 0, 0])
e_chi = sp.Matrix([0, a * sp.exp(U) - a * V, 0, 0])
e_theta = sp.Matrix([0, 0, 1 / (a * chi_var), 0])
e_zeta = sp.Matrix([0, 0, 0, a * sp.exp(W)])

g = sp.diag(-1, 1, 1, 1)

Vu = c * e_tau
Vu_dot_Vu = sp.simplify(Vu.dot(g * Vu))
print("Vu · Vu =", Vu_dot_Vu)

Vd = g * Vu
Vu_dot_Vd = sp.simplify(Vu.dot(Vd))
print("Vu · Vd =", Vu_dot_Vd)

G_00 = (-chi_var*sp.exp(2*U)*sp.diff(U, chi)**2*sp.diff(W, chi) + chi_var*sp.exp(2*U)*sp.diff(V, chi)*sp.diff(W, chi) + chi_var*sp.exp(2*U)*sp.diff(W, chi)**2 - chi_var*sp.exp(2*U)*sp.diff(W, chi, chi) + sp.exp(2*U)*sp.diff(U, chi) - sp.exp(2*U)*sp.diff(V, chi) - sp.exp(2*U)*sp.diff(W, chi)) * sp.exp(-2*V)/(a**2 * chi_var)

G_11 = (-chi_var*sp.exp(2*U)*sp.diff(U, chi)*sp.diff(W, chi) + sp.exp(2*U)*sp.diff(U, chi) - sp.exp(2*U)*sp.diff(W, chi)) * sp.exp(-2*V)/(a**2 * chi_var)

G_22 = (2*sp.exp(2*U)*sp.diff(U, chi)**2 - sp.exp(2*U)*sp.diff(U, chi)*sp.diff(V, chi) - 2*sp.exp(2*U)*sp.diff(U, chi)*sp.diff(W, chi) + sp.exp(2*U)*sp.diff(U, chi, chi) + sp.exp(2*U)*sp.diff(V, chi)*sp.diff(W, chi) + sp.exp(2*U)*sp.diff(W, chi)**2 - sp.exp(2*U)*sp.diff(W, chi, chi)) * sp.exp(-2*V)/a**2

G_33 = (2*chi_var*sp.exp(2*U)*sp.diff(U, chi)**2 - chi_var*sp.exp(2*U)*sp.diff(U, chi)*sp.diff(V, chi) + chi_var*sp.exp(2*U)*sp.diff(U, chi, chi) + 2*sp.exp(2*U)*sp.diff(U, chi) - sp.exp(2*U)*sp.diff(V, chi)) * sp.exp(-2*V)/(a**2 * chi_var)

print("Niezerowe składowe tensora Einsteina (G^i_j):")
print("G^0_0 =", sp.simplify(G_00))
print("G^1_1 =", sp.simplify(G_11))
print("G^2_2 =", sp.simplify(G_22))
print("G^3_3 =", sp.simplify(G_33))
