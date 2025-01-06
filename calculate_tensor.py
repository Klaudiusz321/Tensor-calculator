import sympy as sp

def generate_index_riemann(n):
    index = []
    for a in range(n):
        for b in range(a, n):
            for c in range(n):
                for d in range(c, n):
                    if (a * n + b) <= (c * n + d):
                        index.append((a, b, c, d))
    return index

def generate_index_ricci(n):
    index = []
    for i in range(n):
        for j in range(i, n):
            index.append((i, j))
    return index

def generate_index_christoffel(n):
    index = []
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                index.append((a, b, c))
    return index

def lower_indices(Riemann, g, n):
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd

def write_einstein_components(G_upper, G_lower, n):
    print("Niezerowe składowe tensora Einsteina (G^i_j):")
    for i in range(n):
        for j in range(n):
            val = custom_simplify(G_upper[i, j])
            if val != 0:
                print(f"G^{{{i}}}_{{{j}}} = {val}")
    print("")
    
    print("Niezerowe składowe tensora Einsteina (G_ij):")
    for i in range(n):
        for j in range(n):
            val = custom_simplify(G_lower[i, j])
            if val != 0:
                print(f"G_{{{i}{j}}} = {val}")
    print("")

def write_metric_components(g, n):
    print("Niezerowe składowe metryki (g_{ij}):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(g[i, j])
            if val != 0:
                print(f"g_{i}{j} = {val}")
    print("")

def write_christoffel_symbols(Gamma, n):
    print("Niezerowe symbole Christoffela (Γ^a_{bc}):")
    ch_index = generate_index_christoffel(n)
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        if custom_simplify(val) != 0:
            print(f"Γ^{a}_{{{b}{c}}} = {val}")
    print("")

def write_full_riemann_components(R_abcd, n):
    print("Niezerowe komponenty tensora Riemanna (R_{abcd}):")
    riemann_index = generate_index_riemann(n)
    for (a, b, c, d) in riemann_index:
        val = R_abcd[a][b][c][d]
        if val != 0:
            print(f"R_{a}{b}{c}{d} = {val}")
    print("")

def write_ricci_components(Ricci, n):
    print("Niezerowe komponenty tensora Ricciego (R_{ij}):")
    ricci_index = generate_index_ricci(n)
    for (i, j) in ricci_index:
        val = Ricci[i, j]
        if val != 0:
            print(f"R_{{{i}{j}}} = {val}")
    print("")

def custom_simplify(expr):
    from sympy import simplify, factor, expand, trigsimp, cancel, ratsimp
    
    expr_simpl = expand(expr)
    expr_simpl = trigsimp(expr_simpl)
    expr_simpl = factor(expr_simpl)
    expr_simpl = simplify(expr_simpl)
    expr_simpl = cancel(expr_simpl)
    expr_simpl = ratsimp(expr_simpl)
    
    return expr_simpl

def wczytaj_metryke(filename):
    symbol_assumptions = {
        'a':    dict(real=True, positive=True),
        'tau':  dict(real=True),
        'psi':  dict(real=True),
        'theta':dict(real=True),
        'phi':  dict(real=True),
    }

    def create_symbol(sym_name):
        if sym_name in symbol_assumptions:
            return sp.Symbol(sym_name, **symbol_assumptions[sym_name])
        else:
            return sp.Symbol(sym_name)

    wspolrzedne = []
    parametry = []
    metryka = {}

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.split('#')[0].strip()
                if not line:
                    continue

                if ';' in line:
                    wsp_, prm_ = line.split(';')
                    wsp_strs = [sym.strip() for sym in wsp_.split(',') if sym.strip()]
                    wspolrzedne = [create_symbol(s) for s in wsp_strs]

                    prm_ = prm_.strip()
                    if prm_:
                        par_strs = [sym.strip() for sym in prm_.split(',') if sym.strip()]
                        parametry = [create_symbol(s) for s in par_strs]
                else:
                    dat = line.split(maxsplit=2)
                    if len(dat) == 3:
                        try:
                            i, j, expr = int(dat[0]), int(dat[1]), dat[2]
                            # Tworzymy słownik symboli
                            symbols_dict = {str(sym): sym for sym in wspolrzedne + parametry}
                            metryka[(i, j)] = sp.sympify(expr, locals=symbols_dict)
                        except ValueError:
                            print(f"Error: Incorrect data in line: {line}")
    except FileNotFoundError:
        print(f"Błąd: File has not been found: {filename}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return wspolrzedne, parametry, metryka

def oblicz_tensory(wspolrzedne, metryka):
    n = len(wspolrzedne)

    g = sp.Matrix(n, n, lambda i, j: metryka.get((i, j), metryka.get((j, i), 0)))
    g_inv = g.inv()

    Gamma = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                Gamma_sum = 0
                for lam in range(n):
                    partial_mu  = sp.diff(g[nu, lam], wspolrzedne[mu])
                    partial_nu  = sp.diff(g[mu, lam], wspolrzedne[nu])
                    partial_lam = sp.diff(g[mu, nu], wspolrzedne[lam])
                    Gamma_sum += g_inv[sigma, lam] * (partial_mu + partial_nu - partial_lam)
                Gamma[sigma][mu][nu] = custom_simplify(sp.Rational(1, 2) * Gamma_sum)

    Riemann = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for rho in range(n):
        for sigma in range(n):
            for mu in range(n):
                for nu in range(n):
                    term1 = sp.diff(Gamma[rho][nu][sigma], wspolrzedne[mu])
                    term2 = sp.diff(Gamma[rho][mu][sigma], wspolrzedne[nu])
                    sum_term = 0
                    for lam in range(n):
                        sum_term += (Gamma[rho][mu][lam] * Gamma[lam][nu][sigma]
                                     - Gamma[rho][nu][lam] * Gamma[lam][mu][sigma])
                    Riemann[rho][sigma][mu][nu] = custom_simplify(term1 - term2 + sum_term)

    R_abcd = lower_indices(Riemann, g, n)

    Ricci = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            Ricci[mu, nu] = custom_simplify(sum(Riemann[rho][mu][rho][nu] for rho in range(n)))
            Ricci[mu, nu] = custom_simplify(Ricci[mu, nu])

    Scalar_Curvature = custom_simplify(sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n)))
    Scalar_Curvature = custom_simplify(Scalar_Curvature)

    return g, Gamma, R_abcd, Ricci, Scalar_Curvature

def compute_einstein_tensor(Ricci, Scalar_Curvature, g, g_inv, n):
    G_lower = sp.zeros(n, n)  # G_{ij}
    G_upper = sp.zeros(n, n)  # G^i_j

    # Oblicz G_{ij} = R_{ij} - 1/2 g_{ij} R
    for mu in range(n):
        for nu in range(n):
            G_lower[mu, nu] = custom_simplify(Ricci[mu, nu] - sp.Rational(1, 2) * g[mu, nu] * Scalar_Curvature)

    # Oblicz G^i_j = g^{ik} G_{kj}
    for mu in range(n):
        for nu in range(n):
            sum_term = 0
            for alpha in range(n):
                sum_term += g_inv[mu, alpha] * G_lower[alpha, nu]
            G_upper[mu, nu] = custom_simplify(sum_term)

    return G_upper, G_lower

def wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature, G_upper, G_lower, n):
    write_metric_components(g, n)
    write_christoffel_symbols(Gamma, n)
    write_full_riemann_components(R_abcd, n)
    write_ricci_components(Ricci, n)
    write_einstein_components(G_upper, G_lower, n)

    print("Skalarowa krzywizna R:")
    sp.pprint(Scalar_Curvature)
    print("")

if __name__ == "__main__":

    filename = r"C:\Users\sorak\Desktop\metric.txt"

    wspolrzedne, parametry, metryka = wczytaj_metryke(filename)
    print("Coordinates:", wspolrzedne)
    print("Parameters:", parametry)
    print("Metric (non zero elements):", metryka)
    print("")

    if wspolrzedne and metryka:
        g, Gamma, R_abcd, Ricci, Scalar_Curvature = oblicz_tensory(wspolrzedne, metryka)
        
        # Oblicz odwrotną metrykę
        try:
            g_inv = g.inv()
        except Exception as e:
            print("Błąd przy obliczaniu odwrotnej metryki:", e)
            exit(1)
        
        # Oblicz tensor Einsteina
        G_upper, G_lower = compute_einstein_tensor(Ricci, Scalar_Curvature, g, g_inv, len(wspolrzedne))
        
        # Wyświetl wszystkie tensory
        wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature, G_upper, G_lower, len(wspolrzedne))
