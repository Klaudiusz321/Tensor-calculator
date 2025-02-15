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
def write_einstein_components_latex(G, n):
    print(r"Non-zero Einstein Tensor Components ($G_{ij}$):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(G[i, j])  # Simplify each component
            if val != 0:
                expr_latex = latex(val)  # Convert to LaTeX
                print(rf"G_{{{i}{j}}} &= {expr_latex} \\\\")
    print(r"\end{align*}")

def lower_indices(Riemann, g, n):
    R_abcd = [[[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    R_abcd[a][b][c][d] = sum(g[a, i] * Riemann[i][b][c][d] for i in range(n))
    return R_abcd

<<<<<<< HEAD
def write_scalar_curvatre(scalar_curvature, n):
   
    
    latex_str = sp.latex(scalar_curvature)
    
    
    cleaned_latex = process_latex(latex_str)
    
    print("Scalar curvature R:")
    print(f"$R = {cleaned_latex}$")
    print("Curvature sclar R:")
    sp.pprint(scalar_curvature)
    print("")

def write_einstein_components(G_upper, G_lower, n):
    print("Non-zero Einstein tensor components (textual format and LaTeX):")
    
    
    for i in range(n):
        for j in range(n):
            val = custom_simplify(G_upper[i, j])
            if val != 0:
                
                print(f"G^({i})_({j}) = {val}")
               
                latex_str = sp.latex(val)
                print(f"G^{{{i}}}_{{{j}}} = {latex_str}")
    
    print("\nNon-zero Einstein tensor components (G_ij):")
    for i in range(n):
        for j in range(n):
            val = custom_simplify(G_lower[i, j])
            if val != 0:
               
                print(f"G_({i}{j}) = {val}")
                
                latex_str = sp.latex(val)
                print(f"G_{{{i}{j}}} = {latex_str}")


def write_metric_components(g, n):
    print("Metric tensor components (textual format and LaTeX):")
    for i in range(n):
        for j in range(i, n):
            val = custom_simplify(g[i, j])
            if val != 0:
                # Tekst
                print(f"g_({i}{j}) = {val}")
                # LaTeX
                latex_str = sp.latex(val)
                print(f"g_{{{i}{j}}} = {latex_str}")
                # Dodanie pustej linii
                print()



def write_christoffel_symbols(Gamma, n):
    print("Non-zero Christoffel symbols (textual format and LaTeX):")
    ch_index = generate_index_christoffel(n)
    for (a, b, c) in ch_index:
        val = Gamma[a][b][c]
        simplified_val = custom_simplify(val)
        if simplified_val != 0:
            # Tekstowy format
            print(f"Γ^({a})_({b}{c}) = {simplified_val}")
            
            # Generowanie LaTeX
            latex_str = sp.latex(simplified_val)
            print(f"\\Gamma^{{{a}}}_{{{b}{c}}} = {latex_str}")
            
            # Dodanie pustej linii
            print()


def write_full_riemann_components(R_abcd, n):
    print("Non-zero components of the Riemann tensor (textual format and LaTeX):")
=======
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
>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
    riemann_index = generate_index_riemann(n)
    for (a, b, c, d) in riemann_index:
        val = R_abcd[a][b][c][d]
        if val != 0:
<<<<<<< HEAD
            # Tekstowy format
            print(f"R_({a}{b}{c}{d}) = {val}")
            
            # Generowanie LaTeX
            latex_val = sp.latex(val)
            print(f"R_{{{a}{b}{c}{d}}} = {latex_val}")
            
            # Dodanie pustej linii
            print()



def write_ricci_components(Ricci, n):
    print("Non-zero components of the Ricci tensor (textual format and LaTeX):")
=======
            print(f"R_{a}{b}{c}{d} = {val}")
    print("")

def write_ricci_components(Ricci, n):
    print("Non zero Ricci (R_{ij}):")
>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
    ricci_index = generate_index_ricci(n)
    for (i, j) in ricci_index:
        val = Ricci[i, j]
        if val != 0:
<<<<<<< HEAD
            # Tekstowy format
            print(f"R_({i}{j}) = {val}")
            
            # Generowanie LaTeX
            latex_val = sp.latex(val)
            print(f"R_{{{i}{j}}} = {latex_val}")
            
            # Dodanie pustej linii
            print()
=======
            print(f"R_{{{i}{j}}} = {val}")
    print("")

>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de



def custom_simplify(expr):
<<<<<<< HEAD
    from sympy import simplify, factor, expand, trigsimp, cancel, ratsimp
    
    expr_simpl = expand(expr)
    expr_simpl = trigsimp(expr_simpl)
    expr_simpl = factor(expr_simpl)
    expr_simpl = simplify(expr_simpl)
    expr_simpl = cancel(expr_simpl)
    expr_simpl = ratsimp(expr_simpl)
    
    return expr_simpl
def process_latex(latex_str):
   
    def remove_function_argument(latex):
        result = ""
        i = 0
        while i < len(latex):
            if latex[i].isalpha():
                start = i
                while i < len(latex) and latex[i].isalpha():
                    i += 1
                func_name = latex[start:i]
                
                if i < len(latex) and latex[i] == '(':
                    i += 1  
                    arg_start = i
                    paren_count = 1
                    while i < len(latex) and paren_count > 0:
                        if latex[i] == '(':
                            paren_count += 1
                        elif latex[i] == ')':
                            paren_count -= 1
                        i += 1
                    arg = latex[arg_start:i-1].strip()
                    if arg in ['\\chi', 'chi']:
                       
                        result += func_name
                    else:
                        
                        result += f"{func_name}({arg})"
                else:
                   
                    result += func_name
            else:
                
                result += latex[i]
                i += 1
        return result

 
    def replace_derivative(latex):
        search_str = "\\frac{d}{d \\chi}"
        while search_str in latex:
            index = latex.find(search_str)
            
            after = latex[index + len(search_str):]
           
            var_end = index + len(search_str)
            while var_end < len(latex) and (latex[var_end].isalpha() or latex[var_end] == '\\'):
                var_end += 1
            var = latex[index + len(search_str):var_end].strip()
            
            if var.startswith('\\'):
                var_name = ""
                j = 0
                while j < len(var) and (var[j].isalpha() or var[j] == '\\'):
                    var_name += var[j]
                    j += 1
                var_replaced = f"{var_name}'"
            else:
                var_replaced = f"{var}'"
          
            latex = latex[:index] + var_replaced + latex[var_end:]
        return latex

    
    latex_str = remove_function_argument(latex_str)
   
    latex_str = replace_derivative(latex_str)

    return latex_str

def wczytaj_metryke(filename):
=======
    from sympy import simplify, factor, expand, expand_trig, expand_log, trigsimp

    expr_simpl = expand(expr)
    expr_simpl = expand_trig(expr_simpl)
    expr_simpl = expand_log(expr_simpl)


    expr_simpl = trigsimp(expr_simpl)
    expr_simpl = factor(expr_simpl)
    expr_simpl = simplify(expr_simpl)

    return expr_simpl


def wczytaj_metryke(filename):

>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
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
<<<<<<< HEAD
                           
=======
                            # Tworzymy słownik symboli
>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
                            symbols_dict = {str(sym): sym for sym in wspolrzedne + parametry}
                            metryka[(i, j)] = sp.sympify(expr, locals=symbols_dict)
                        except ValueError:
                            print(f"Error: Incorrect data in line: {line}")
    except FileNotFoundError:
        print(f"Błąd: File has not been found: {filename}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return wspolrzedne, parametry, metryka

<<<<<<< HEAD
=======


>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
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
<<<<<<< HEAD
            Ricci[mu, nu] = custom_simplify(sum(Riemann[rho][mu][rho][nu] for rho in range(n)))
            Ricci[mu, nu] = custom_simplify(Ricci[mu, nu])

    Scalar_Curvature = custom_simplify(sum(g_inv[mu, nu] * Ricci[mu, nu] for mu in range(n) for nu in range(n)))
    Scalar_Curvature = custom_simplify(Scalar_Curvature)

    return g, Gamma, R_abcd, Ricci, Scalar_Curvature

def compute_einstein_tensor(Ricci, Scalar_Curvature, g, g_inv, n):
    G_lower = sp.zeros(n, n)  
    G_upper = sp.zeros(n, n) 

    
    for mu in range(n):
        for nu in range(n):
            G_lower[mu, nu] = custom_simplify(Ricci[mu, nu] - sp.Rational(1, 2) * g[mu, nu] * Scalar_Curvature)


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
    write_scalar_curvatre(Scalar_Curvature, n)
    # print("Curvature sclar R:")
    # sp.pprint(Scalar_Curvature)
    print("")

=======
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
>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de

if __name__ == "__main__":

<<<<<<< HEAD
=======
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

>>>>>>> 4e0811047bac8be511ea2be8f04b9cf3687c08de
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
            print("error:", e)
            exit(1)
        
        # Oblicz tensor Einsteina
        G_upper, G_lower = compute_einstein_tensor(Ricci, Scalar_Curvature, g, g_inv, len(wspolrzedne))
        
        # Wyświetl wszystkie tensory
        wyswietl_tensory(g, Gamma, R_abcd, Ricci, Scalar_Curvature, G_upper, G_lower, len(wspolrzedne))
