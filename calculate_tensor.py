import sympy as sp
from sympy import latex


def generate_index_riemann(n):
    """Generate independent indices for Riemann tensor components."""
    index = []
    for a in range(n):
        for b in range(a, n):
            for c in range(n):
                for d in range(c, n):
                    if (a*n + b) <= (c*n + d):
                        index.append((a, b, c, d))
    return index


def generate_index_ricci(n):
    """Generate independent indices for Ricci tensor components."""
    index = []
    for i in range(n):
        for j in range(i, n):
            index.append((i, j))
    return index


def generate_index_christoffel(n):
    """Generate independent indices for Christoffel symbols."""
    index = []
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                index.append((a, b, c))
    return index


def lower_indices(Riemann, g, n):
    """Lower the first index of Riemann tensor using metric tensor."""
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd


def write_einstein_components(G, n):
    """Print non-zero Einstein tensor components."""
    print("Non-zero Einstein tensor components (G_{ij}):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(G[i, j])
            if val != 0:
                print(f"G_{{{i}{j}}} = {val}")
    print("")


def write_metric_components(g, n):
    """Print non-zero metric tensor components."""
    print("Non-zero metric tensor components (g_{ij}):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(g[i, j])
            if val != 0:
                print(f"g_{{{i}{j}}} = {val}")
    print("")


def write_christoffel_symbols(Gamma, n):
    """Print non-zero Christoffel symbols."""
    print("Non-zero Christoffel symbols (Γ^a_{bc}):")
    ch_index = generate_index_christoffel(n)
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        if custom_simplify(val) != 0:
            print(f"Γ^{a}_{{{b}{c}}} = {val}")
    print("")


def write_full_riemann_components(R_abcd, n):
    """Print non-zero Riemann tensor components."""
    print("Non-zero Riemann tensor components (R_{abcd}):")
    riemann_index = generate_index_riemann(n)
    for (a, b, c, d) in riemann_index:
        val = R_abcd[a][b][c][d]
        if val != 0:
            print(f"R_{{{a}{b}{c}{d}}} = {val}")
    print("")


def write_ricci_components(Ricci, n):
    """Print non-zero Ricci tensor components."""
    print("Non-zero Ricci tensor components (R_{ij}):")
    ricci_index = generate_index_ricci(n)
    for (i, j) in ricci_index:
        val = Ricci[i, j]
        if val != 0:
            print(f"R_{{{i}{j}}} = {val}")
    print("")


def custom_simplify(expr):
    """Apply multi-step symbolic simplification to expression."""
    from sympy import simplify, factor, expand, trigsimp, cancel, ratsimp
    
    expr_simpl = expand(expr)
    expr_simpl = trigsimp(expr_simpl)
    expr_simpl = factor(expr_simpl)
    expr_simpl = simplify(expr_simpl)
    expr_simpl = cancel(expr_simpl)
    expr_simpl = ratsimp(expr_simpl)
    
    return expr_simpl


def load_metric(filename):
    """
    Load metric tensor from a text file.
    
    File format:
    - First line: coordinates ; parameters (e.g., "t, r, theta, phi ; M")
    - Following lines: i j expression (e.g., "0 0 -1")
    - Comments start with #
    
    Returns:
        tuple: (coordinates, parameters, metric_dict)
    """
    symbol_assumptions = {
        'a':    dict(real=True, positive=True),
        'tau':  dict(real=True),
        'psi':  dict(real=True),
        'theta': dict(real=True),
        'phi':  dict(real=True),
    }

    def create_symbol(sym_name):
        if sym_name in symbol_assumptions:
            return sp.Symbol(sym_name, **symbol_assumptions[sym_name])
        else:
            return sp.Symbol(sym_name)

    coordinates = []
    parameters = []
    metric = {}

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.split('#')[0].strip()
                if not line:
                    continue

                if ';' in line:
                    coord_, param_ = line.split(';')
                    coord_strs = [sym.strip() for sym in coord_.split(',') if sym.strip()]
                    coordinates = [create_symbol(s) for s in coord_strs]

                    param_ = param_.strip()
                    if param_:
                        param_strs = [sym.strip() for sym in param_.split(',') if sym.strip()]
                        parameters = [create_symbol(s) for s in param_strs]
                else:
                    dat = line.split(maxsplit=2)
                    if len(dat) == 3:
                        try:
                            i, j, expr = int(dat[0]), int(dat[1]), dat[2]
                            # Create symbol dictionary
                            symbols_dict = {str(sym): sym for sym in coordinates + parameters}
                            metric[(i, j)] = sp.sympify(expr, locals=symbols_dict)
                        except ValueError:
                            print(f"Error: Incorrect data in line: {line}")
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return coordinates, parameters, metric


def compute_tensors(coordinates, metric):
    """
    Compute all geometric tensors from the metric.
    
    Args:
        coordinates: List of coordinate symbols
        metric: Dictionary of metric components
        
    Returns:
        tuple: (g, Gamma, R_abcd, Ricci, Scalar_Curvature)
    """
    n = len(coordinates)

    # Build metric matrix
    g = sp.Matrix(n, n, lambda i, j: metric.get((i, j), metric.get((j, i), 0)))
    g_inv = g.inv()

    # Compute Christoffel symbols
    Gamma = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for sigma in range(n):
        for mu in range(n):
            for nu in range(n):
                Gamma_sum = 0
                for lam in range(n):
                    partial_mu = sp.diff(g[nu, lam], coordinates[mu])
                    partial_nu = sp.diff(g[mu, lam], coordinates[nu])
                    partial_lam = sp.diff(g[mu, nu], coordinates[lam])
                    Gamma_sum += g_inv[sigma, lam] * (partial_mu + partial_nu - partial_lam)
                Gamma[sigma][mu][nu] = custom_simplify(sp.Rational(1, 2) * Gamma_sum)

    # Compute Riemann tensor
    Riemann = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for rho in range(n):
        for sigma in range(n):
            for mu in range(n):
                for nu in range(n):
                    term1 = sp.diff(Gamma[rho][nu][sigma], coordinates[mu])
                    term2 = sp.diff(Gamma[rho][mu][sigma], coordinates[nu])
                    sum_term = 0
                    for lam in range(n):
                        sum_term += (Gamma[rho][mu][lam] * Gamma[lam][nu][sigma]
                                     - Gamma[rho][nu][lam] * Gamma[lam][mu][sigma])
                    Riemann[rho][sigma][mu][nu] = custom_simplify(term1 - term2 + sum_term)

    # Lower indices on Riemann tensor
    R_abcd = lower_indices(Riemann, g, n)

    # Compute Ricci tensor
    Ricci = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            Ricci[mu, nu] = custom_simplify(sum(Riemann[rho][mu][rho][nu] for rho in range(n)))
            Ricci[mu, nu] = custom_simplify(Ricci[mu, nu])

    # Compute scalar curvature
    Scalar_Curvature = custom_simplify(sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n)))
    Scalar_Curvature = custom_simplify(Scalar_Curvature)

    return g, Gamma, R_abcd, Ricci, Scalar_Curvature


def compute_einstein_tensor(Ricci, Scalar_Curvature, g, n):
    """Compute Einstein tensor from Ricci tensor and scalar curvature."""
    G = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            G[mu, nu] = custom_simplify(Ricci[mu, nu] - sp.Rational(1, 2) * g[mu, nu] * Scalar_Curvature)
            G[mu, nu] = custom_simplify(G[mu, nu])
    return G


def display_tensors(g, Gamma, R_abcd, Ricci, Scalar_Curvature, G, n):
    """Display all computed tensor components."""
    write_metric_components(g, n)
    write_christoffel_symbols(Gamma, n)
    write_full_riemann_components(R_abcd, n)
    write_ricci_components(Ricci, n)
    write_einstein_components(G, n)

    print("Scalar curvature R:")
    sp.pprint(Scalar_Curvature)
    print("")


if __name__ == "__main__":
    # Path to metric file
    filename = "metric.txt"

    # Load metric from file
    coordinates, parameters, metric = load_metric(filename)
    print("Coordinates:", coordinates)
    print("Parameters:", parameters)
    print("Metric (non-zero elements):", metric)
    print("")

    # Compute and display all tensors
    if coordinates and metric:
        g, Gamma, R_abcd, Ricci, Scalar_Curvature = compute_tensors(coordinates, metric)
        G = compute_einstein_tensor(Ricci, Scalar_Curvature, g, len(coordinates))

        display_tensors(g, Gamma, R_abcd, Ricci, Scalar_Curvature, G, len(coordinates))
