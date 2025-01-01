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

# Funkcja do obniżania indeksów Riemanna (nieużywana w tym programie)
def lower_indices(Riemann, g, n):
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd

# Funkcja do tworzenia tensorów energii-pędu
def energy_momentum_tensor(gdd, A, rho, p, n):
    T = sp.MutableDenseNDimArray.zeros(n, n)
    # Zakładamy, że czteroprędność ma tylko pierwszy składnik niezerowy
    T[0, 0] = rho
    for i in range(1, n):
        T[i, i] = p
    return T

# Funkcja do obliczania tensora Einsteina
def compute_einstein_tensor(Ricci, Scalar_Curvature, g_inv, n):
    Gud = sp.MutableDenseNDimArray.zeros(n, n)
    for i in range(n):
        for j in range(n):
            Gud[i, j] = Ricci[i, j] - (sp.Rational(1, 2) * g_inv[i, j] * Scalar_Curvature)
            Gud[i, j] = Gud[i, j].simplify()
    return Gud

# Funkcja do wczytywania metryki z pliku
def wczytaj_metryke(filename):
    # Implementacja funkcji wczytującej metrykę z pliku
    # Zakładamy, że plik zawiera linie w formacie: i j expr
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
                            metryka[(i, j)] = sp.sympify(expr, locals={"Symbol": create_symbol})
                        except ValueError:
                            print(f"Error: Incorrect data in line: {line}")
    except FileNotFoundError:
        print(f"Błąd: File has not been found: {filename}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return wspolrzedne, parametry, metryka

# Funkcja do obliczania tensorów
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
            Ricci[mu, nu] = sum(Riemdddd[rho, mu, rho, nu] for rho in range(n)).simplify()

    # Obliczanie skalarnej krzywizny
    Scalar_Curvature = sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n)).simplify()

    # Definicja czteroprędności
    # Zakładamy, że A = exp(-psi(chi))
    A_expr = sp.exp(-psi)

    # Obliczenie tensora energii-pędu
    T = energy_momentum_tensor(gdd, A_expr, rho, p, n)

    # Obliczenie tensora Einsteina
    Gud = compute_einstein_tensor(Ricci, Scalar_Curvature, g_inv, n)

    return gdd, Gammaudd, Riemdddd, Ricci, Scalar_Curvature, Gud, T


# Funkcja do wczytywania metryki z pliku (ta sama co wcześniej)
# [Została już zdefiniowana powyżej]

# Funkcja do wyświetlania tensorów
def wyswietl_tensory(gdd, Gammaudd, Riemdddd, Ricci, Scalar_Curvature, Gud, n, T):
    # Tekstowe Wyjścia
    print("Tensor energii-pędu (T^i_j):")
    for i in range(n):
        for j in range(n):
            val = T[i, j]
            if val != 0:
                tensor_str = format_tensor_text("T", [i], [j])
                print(f"{tensor_str} = {val}")
    print("")

    print("Tensor Einsteina (G^i_j):")
    for i in range(n):
        for j in range(n):
            val = Gud[i, j]
            if val != 0:
                tensor_str = format_tensor_text("G", [i], [j])
                print(f"{tensor_str} = {val}")
    print("")

    print("Skalar krzywizny R:")
    sp.pprint(Scalar_Curvature)
    print("")

    # Możesz dodać inne wyświetlenia tensorów według potrzeb

# Funkcja główna
def main():
    # Definicje stałych fizycznych (przykładowe wartości)
    G_val, c_val, m_val, n0_val, beta_val = 6.67430e-11, 3e8, 1.0, 1.0, 0.1  # Jednostki SI

    # Ustalanie parametru długości a
    a_val = (c_val**2) / (G_val * m_val * n0_val)

    # Ścieżka do pliku z metryką
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
