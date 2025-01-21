import sympy as sp

# Definicja zmiennych
rho = sp.Symbol('rho', positive=True)  # Promień
P = sp.Function('P')(rho)             # Ciśnienie jako funkcja promienia
Phi = sp.Function('Phi')(rho)         # Potencjał jako funkcja promienia

# Równanie różniczkowe ciśnienia
dP_drho = sp.Eq(sp.diff(P, rho), - (rho + P) * sp.diff(Phi, rho))

# Przykładowa zależność Phi (ułatwia rozwiązanie)
Phi_expr = sp.log(rho)  # Przykładowa funkcja Phi
dP_drho = dP_drho.subs(sp.diff(Phi, rho), sp.diff(Phi_expr, rho))

# Rozwiązanie symboliczne
solution = sp.dsolve(dP_drho)
print("Rozwiązanie symboliczne:")
print(solution)
