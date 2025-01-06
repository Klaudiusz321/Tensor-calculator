import sympy as sp
from calculate_tensor import compute_einstein_tensor


def define_basis_vectors(parametry):
    a, chi = sp.symbols('a chi', positive=True)

    e_tau = sp.Matrix([sp.exp(-parametry[0])/a, 0, 0, 0])
    e_x = sp.Matrix([0, sp.exp(parametry[0]-parametry[1])/a, 0, 0])
    e_phi = sp.Matrix([0, 0, 1/(a * chi), 0])
    e_xi = sp.Matrix([0, 0, 0, sp.exp(parametry[2])/a])

    return [e_tau, e_x, e_phi, e_xi]

def compute_velocity_vector(c, basis_vector):
    return c * basis_vector


def check_velocity_norm(Vu, g):
    return (Vu.T * g * Vu)[0]

if __name__ == "__main__":

    with open("einstein_tensor.txt", "r") as file:
        einstein_tensor = sp.sympify(file.read())


    parametry = [sp.Symbol('u'), sp.Symbol('v'), sp.Symbol('w')]
    basis_vectors = define_basis_vectors(parametry)

 
    c = sp.Symbol('c', positive=True)
    Vu = compute_velocity_vector(c, basis_vectors[0])

  
    g = sp.diag(-1, 1, 1, 1)


    norm_Vu = check_velocity_norm(Vu, g)
    print("Norma Vu:", norm_Vu)

    T = sp.Matrix([
        [einstein_tensor[0, 0] * Vu[0], einstein_tensor[0, 1] * Vu[0]],
        [einstein_tensor[1, 0] * Vu[1], einstein_tensor[1, 1] * Vu[1]]
    ])
    print("Tensor energii-pędu:", T)
