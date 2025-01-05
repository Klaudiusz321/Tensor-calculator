import sympy as sp
from sympy import latex


def wyswietl_tensory(gdd, Gammaudd, Riemdddd, Ricci, Scalar_Curvature, Gud, n, T):
    # -- już masz: wyświetlanie T, G, skalar krzywizny, itd. --
    # (fragmenty pomijam)

    # ---------------------------------------------
    # 1. WYŚWIETLANIE SYMBOLI CHRISTOFFELA Γ^σ_{μν}
    # ---------------------------------------------
    print("\nSymbole Christoffela (Gamma^sigma_{mu nu}):")
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                val = Gammaudd[sigma, mu, nu]
                # Jeśli chcesz zobaczyć też 0, usuń warunek:
                if val != 0:
                    gamma_str = format_tensor_text("Γ", [sigma], [mu, nu])
                    print(f"{gamma_str} = {val}")
    print("")

    # -----------------------------------------------------
    # 2. WYŚWIETLANIE TENORA RIEMANNA (R^rho_{sigma mu nu})
    # -----------------------------------------------------
    print("Tensor Riemanna (R^rho_{sigma mu nu}):")
    for rho in range(n):
        for sigma in range(n):
            for mu in range(n):
                for nu in range(n):
                    val = Riemdddd[rho, sigma, mu, nu]
                    if val != 0:
                        r_str = format_tensor_text("R", [rho], [sigma, mu, nu])
                        print(f"{r_str} = {val}")
    print("")

    # ------------------------------------------------
    # 3. WYŚWIETLANIE TENORA RICCI (R_{mu nu})
    # ------------------------------------------------
    print("Tensor Ricciego (R_{mu nu}):")
    for mu in range(n):
        for nu in range(n):
            val = Ricci[mu, nu]
            if val != 0:
                ricci_str = format_tensor_text("R", [], [mu, nu])
                print(f"{ricci_str} = {val}")
    print("")



# Definicje symboli
# Stałe fizyczne
G, c, m, n0, beta = sp.symbols('G c m n0 beta', real=True, positive=True)

# Bezrozmiarowa zmienna radialna
chi = sp.Symbol('chi', real=True)

# Funkcja opisująca bezwymiarową gęstość liczby cząstek
nu = sp.Function('nu')(chi)

# Funkcje metryczne
psi = sp.Function('psi')(chi)
gamma_func = sp.Function('gamma')(chi)
mu_func = sp.Function('mu')(chi)
nu_func = sp.Function('nu_func')(chi)  # Unikaj konfliktu nazw z bezwymiarową gęstością

# Gęstość energii i ciśnienie
rho = m * c**2 * (nu)**(beta**2 + 1)
p = beta**2 * rho

# Definicja czteroprędności
A = sp.Function('A')(chi)

# Funkcja do formatowania tensorów w formacie tekstowym
def format_tensor_text(name, upper_indices, lower_indices):
    tensor_str = name
    if upper_indices:
        tensor_str += "^" + "".join(str(i) for i in upper_indices)
    if lower_indices:
        tensor_str += "_" + "".join(str(j) for j in lower_indices)
    return tensor_str

# Funkcja do formatowania tensorów w formacie LaTeX
def format_tensor_latex(name, upper_indices, lower_indices):
    tensor_str = name
    if upper_indices:
        tensor_str += "^{" + " ".join(str(i) for i in upper_indices) + "}"
    if lower_indices:
        if upper_indices:
            tensor_str += "_{" + "\\phantom{" + str(upper_indices[0]) + "}" + "".join(str(j) for j in lower_indices) + "}"
        else:
            tensor_str += "_{" + "".join(str(j) for j in lower_indices) + "}"
    return tensor_str

# Funkcja do generowania indeksów Riemanna (dla n wymiarów)
def generate_index_riemann(n):
    index = []
    for a in range(n):
        for b in range(a, n):
            for c in range(n):
                for d in range(c, n):
                    if (a * n + b) <= (c * n + d):
                        index.append((a, b, c, d))
    return index

# Funkcja do generowania indeksów Ricciego
def generate_index_ricci(n):
    index = []
    for i in range(n):
        for j in range(i, n):
            index.append((i, j))
    return index

# Funkcja do generowania indeksów symboli Christoffela
def generate_index_christoffel(n):
    index = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                index.append((a, b, c))
    return index
def write_einstein_components_latex(G, n):
    print(r"Non-zero Einstein Tensor Components ($G_{ij}$):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(G[i, j])  # Simplify each component
            if val != 0:
                expr_latex = latex(val)  # Convert to LaTeX
                print(rf"G_{{{i}{j}}} &= {expr_latex} \\\\")
    print(r"\end{align*}")

# Funkcja do obniżania indeksów Riemanna (nieużywana w tym programie)
def lower_indices(Riemann, g, n):
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd

def write_einstein_components(G, n):
    print("Non zero einstein tensor components (G_ij):")
    for i in range(n):
        for j in range(i,n):
            val = G[i,j]
            if val !=0:
                print(f"G{{{i}{j}}} = {val}")
    print("")

def write_metric_components(g, n):
    print("non zero metric tensor (g_{ij}):")
    for i in range(n):
        for j in range(i, n):
            val = g[i, j]
            if val != 0:
                print(f"g_{i}{j} = {val}")
    print("")

def write_christoffel_symbols(Gamma, n):
    print("Non zero Christoffela (Γ^a_{bc}):")
    ch_index = generate_index_christoffel(n)
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        if val != 0:
            print(f"Γ^{a}_{{{b}{c}}} = {val}")
    print("")

def write_full_riemann_components(R_abcd, n):
    print("Non zero Riemann (R_{abcd}):")
    riemann_index = generate_index_riemann(n)
    for (a, b, c, d) in riemann_index:
        val = R_abcd[a][b][c][d]
        if val != 0:
            print(f"R_{a}{b}{c}{d} = {val}")
    print("")

def write_ricci_components(Ricci, n):
    print("Non zero Ricci (R_{ij}):")
    ricci_index = generate_index_ricci(n)
    for (i, j) in ricci_index:
        val = Ricci[i, j]
        if val != 0:
            print(f"R_{{{i}{j}}} = {val}")
    print("")




def custom_simplify(expr):
    from sympy import simplify, factor, expand, expand_trig, expand_log, trigsimp

    expr_simpl = expand(expr)
    expr_simpl = expand_trig(expr_simpl)
    expr_simpl = expand_log(expr_simpl)


    expr_simpl = trigsimp(expr_simpl)
    expr_simpl = factor(expr_simpl)
    expr_simpl = simplify(expr_simpl)

    return expr_simpl


def wczytaj_metryke(filename):

    symbol_assumptions = {
        'a':    dict(real=True, positive=True),
        'psi':  dict(real=True),
        'gamma':dict(real=True),
        'mu':   dict(real=True),
        'nu_func': dict(real=True),
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

    # Tworzenie macierzy metryki
    gdd = sp.Matrix(n, n, lambda i, j: metryka.get((i, j), metryka.get((j, i), 0)))
    g_inv = gdd.inv()

    # Obliczanie symboli Christoffela
    Gammaudd = sp.MutableDenseNDimArray.zeros(n, n, n)
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                Gamma_sum = 0
                for lam in range(n):
                    partial_mu  = sp.diff(gdd[nu, lam], wspolrzedne[mu])
                    partial_nu  = sp.diff(gdd[mu, lam], wspolrzedne[nu])
                    partial_lam = sp.diff(gdd[mu, nu], wspolrzedne[lam])
                    Gamma_sum += g_inv[sigma, lam] * (partial_mu + partial_nu - partial_lam)
                Gammaudd[sigma, mu, nu] = (sp.Rational(1, 2) * Gamma_sum).simplify()

    # Obliczanie tensora Riemanna
    Riemdddd = sp.MutableDenseNDimArray.zeros(n, n, n, n)
    for rho in range(n):
        for sigma in range(n):
            for mu in range(n):
                for nu in range(n):
                    term1 = sp.diff(Gammaudd[rho, nu, sigma], wspolrzedne[mu])
                    term2 = sp.diff(Gammaudd[rho, mu, sigma], wspolrzedne[nu])
                    sum_term = 0
                    for lam in range(n):
                        sum_term += (Gammaudd[rho, mu, lam] * Gammaudd[lam, nu, sigma]
                                     - Gammaudd[rho, nu, lam] * Gammaudd[lam, mu, sigma])
                    Riemdddd[rho, sigma, mu, nu] = (term1 - term2 + sum_term).simplify()

    # Obliczanie tensora Ricciego
    Ricci = sp.MutableDenseNDimArray.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            Ricci[mu, nu] = sum(Riemann[rho][mu][rho][nu] for rho in range(n))
            Ricci[mu, nu] = custom_simplify(Ricci[mu, nu])


    Scalar_Curvature = sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n))
    Scalar_Curvature = custom_simplify(Scalar_Curvature)

    return g, Gamma, R_abcd, Ricci, Scalar_Curvature


# def write_full_riemann_components_latex(R_abcd, n):
#     riemann_index = generate_index_riemann(n)
#     print(r"\begin{align*}")
#     for (a, b, c, d) in riemann_index:
#         val = R_abcd[a][b][c][d]
#         if val != 0:
#             expr_latex = latex(val)
#             print(rf"R_{{{a}{b}{c}{d}}} &= {expr_latex} \\\\")
#     print(r"\end{align*}")

def compute_einstein_tensor(Ricci, Scalar_Curvature, g, n):
    G = sp.zeros(n,n)
    for mu in range(n):
        for nu in range(n):
            G[mu, nu] = Ricci[mu, nu] - sp.Rational(1, 2)*g[mu, nu]* Scalar_Curvature
            G[mu, nu] = custom_simplify(G[mu, nu])
    return G

# def write_ricci_components_latex(Ricci, n):
#     ricci_index = generate_index_ricci(n)
#     print(r"\begin{align*}")
#     for (i, j) in ricci_index:
#         val = Ricci[i, j]
#         if val != 0:
#             expr_latex = latex(val)
#             print(rf"R_{{{i}{j}}} &= {expr_latex} \\\\")
#     print(r"\end{align*}")

def write_christoffel_symbols_latex(Gamma, n):
    ch_index = generate_index_christoffel(n)
    print(r"\begin{align*}")
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        if val != 0:
            expr_latex = latex(val)
            print(rf"\Gamma^{{{a}}}_{{{b}{c}}} &= {expr_latex} \\\\")
    print(r"\end{align*}")

def write_scalar_curvature_latex(Scalar_Curvature):
    scalar_latex = latex(Scalar_Curvature)
    print(r"\begin{equation*}")
    print(rf"R = {scalar_latex}")
    print(r"\end{equation*}")




def wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature,G, n):
    write_metric_components(g, n)
    write_christoffel_symbols(Gamma, n)
    write_full_riemann_components(R_abcd, n)
    write_ricci_components(Ricci, n)
    write_einstein_components(G, n)

    print("Skalar krzywizny R:")
    sp.pprint(Scalar_Curvature)
    print("")

    # print("\n--- LaTeX Code for Christoffel Symbols ---")
    # write_christoffel_symbols_latex(Gamma, n)

    # print("\n--- LaTeX Code for Riemann Tensor ---")
    # write_full_riemann_components_latex(R_abcd, n)
    #
    # print("\n--- LaTeX Code for Ricci Tensor ---")
    # write_ricci_components_latex(Ricci, n)

    # print("\n--- LaTeX Code for Scalar Curvature ---")
    # write_scalar_curvature_latex(Scalar_Curvature)
    # print("")





if __name__ == "__main__":

    filename = r"C:\Users\sorak\Desktop\metric.txt"

    # Wczytanie metryki
    wspolrzedne, parametry, metryka = wczytaj_metryke(filename)
    print("Coordinates:", wspolrzedne)
    print("Parameters:", parametry)
    print("Metric (non zero elements):", metryka)
    print("")

    if wspolrzedne and metryka:
        gdd, Gammaudd, Riemdddd, Ricci, Scalar_Curvature, Gud, T = oblicz_tensory(wspolrzedne, metryka)

        # Wyświetlenie tensorów energii-pędu i tensora Einsteina
        wyswietl_tensory(gdd, Gammaudd, Riemdddd, Ricci, Scalar_Curvature, Gud, len(wspolrzedne), T)

# Uruchomienie programu
if __name__ == "__main__":
    main()
