import sympy as sp

def generate_index_riemann(n):
    # Generuje indeksy (a,b,c,d) uwzględniając symetrie Riemanna
    index = []
    for a in range(n):
        for b in range(a, n):
            for c in range(n):
                for d in range(c, n):
                    if (a * n + b) <= (c * n + d):
                        index.append((a, b, c, d))
    return index

def generate_index_ricci(n):
    # Generuje indeksy (i,j) dla Ricciego
    index = []
    for i in range(n):
        for j in range(i, n):
            index.append((i, j))
    return index

def generate_index_christoffel(n):
    # Generuje indeksy (a,b,c) dla symboli Christoffela Γ^a_{bc}, b ≤ c
    index = []
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                index.append((a, b, c))
    return index

def lower_indices(Riemann, g, n):
    # Obniżenie indeksów tensora Riemanna: R_{abcd} = g_{aρ} R^ρ_{bcd}
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd

def write_metric_components(g, n):
    print("non-zero components of the metric tensor:")
    for i in range(n):
        for j in range(i, n):
            val = g[i,j]
            if val != 0:
                print(f"g_{i}{j}:    {val}")
    print("")

def write_christoffel_symbols(Gamma, n):
    print("non-zero components of Christoffel symbols:")
    ch_index = generate_index_christoffel(n)
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        if val != 0:
            print(f"Γ^{a}_{b}{c}:    {val}")
    print("")

def write_full_riemann_components(R_abcd, n):
    print("non-zero components of the Riemann tensor:")
    riemann_index = generate_index_riemann(n)
    for (a, b, c, d) in riemann_index:
        val = R_abcd[a][b][c][d]
        if val != 0:
            print(f"R_{a}{b}{c}{d}:    {val}")
    print("")

def write_ricci_components(Ricci, n):
    print("non-zero components of the Ricci tensor:")
    ricci_index = generate_index_ricci(n)
    for (i, j) in ricci_index:
        val = Ricci[i, j]
        if val != 0:
            print(f"R_{i}{j}:    {val}")
    print("")

def custom_simplify(expr):
    # Można stopniowo usuwać kolejne kroki upraszczania, jeśli pojawiają się ostrzeżenia lub problemy.
    expr = sp.simplify(expr, rational=True)
    return expr

def wczytaj_metryke(filename):
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
                    wspolrzedne = [sp.Symbol(sym.strip()) for sym in wsp_.split(',')]
                    prm_ = prm_.strip()
                    if prm_:
                        parametry = [sp.Symbol(sym.strip()) for sym in prm_.split(',')]
                else:
                    dat = line.split(maxsplit=2)
                    if len(dat) == 3:
                        try:
                            i, j, expr = int(dat[0]), int(dat[1]), dat[2]
                            metryka[(i, j)] = sp.sympify(expr)
                        except ValueError:
                            print(f"Błąd: Nieprawidłowe dane w linii: {line}")
    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku: {filename}")
    except Exception as e:
        print(f"Wystąpił nieoczekiwany błąd: {e}")

    return wspolrzedne, parametry, metryka

def oblicz_tensory(wspolrzedne, metryka):
    n = len(wspolrzedne)
    g = sp.Matrix(n, n, lambda i, j: metryka.get((i, j), metryka.get((j, i), 0)))
    g_inv = g.inv()

    # Obliczenie symboli Christoffela
    Gamma = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                Gamma_sum = 0
                for lam in range(n):
                    partial_mu = sp.diff(g[nu, lam], wspolrzedne[mu])
                    partial_nu = sp.diff(g[mu, lam], wspolrzedne[nu])
                    partial_lam = sp.diff(g[mu, nu], wspolrzedne[lam])
                    Gamma_sum += g_inv[sigma, lam] * (partial_mu + partial_nu - partial_lam)
                Gamma[sigma][mu][nu] = custom_simplify(sp.Rational(1, 2) * Gamma_sum)

    # Tensor Riemanna R^ρ_{σμν}
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

    # Obniżenie indeksów w tensorze Riemanna: R_{abcd}
    R_abcd = lower_indices(Riemann, g, n)

    # Tensor Ricciego: R_ij = R^ρ_{iρj}
    Ricci = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            Ricci[mu, nu] = sum(Riemann[rho][mu][rho][nu] for rho in range(n))
            Ricci[mu, nu] = custom_simplify(Ricci[mu, nu])

    # Skalar krzywizny: R = g^{μν}R_{μν}
    Scalar_Curvature = sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n))
    Scalar_Curvature = custom_simplify(Scalar_Curvature)

    return g, Gamma, R_abcd, Ricci, Scalar_Curvature

def wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature, n):
    write_metric_components(g, n)
    write_christoffel_symbols(Gamma, n)
    write_full_riemann_components(R_abcd, n)
    write_ricci_components(Ricci, n)
    print("scalar curvature:")
    sp.pprint(Scalar_Curvature)
    print("")

# Przykład użycia - proszę dostosować ścieżkę do pliku metric.txt
filename = r"C:\Users\sorak\Desktop\metric.txt"
wspolrzedne, parametry, metryka = wczytaj_metryke(filename)
print("Współrzędne:", wspolrzedne)
print("Parametry:", parametry)
print("Metryka:", metryka)
print("")

if wspolrzedne and metryka:
    g, Gamma, R_abcd, Ricci, Scalar_Curvature = oblicz_tensory(wspolrzedne, metryka)
    wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature, len(wspolrzedne))