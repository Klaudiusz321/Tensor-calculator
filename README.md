# Tensor Calculator for General Relativity

## 📖 Project Description

This program calculates key geometric objects used in **General Relativity (GR)**. Based on a given spacetime metric, the program automatically computes:

- **Metric tensor** $g_{\mu\nu}$
- **Christoffel symbols** $\Gamma^\sigma_{\mu\nu}$
- **Riemann tensor** $R_{\rho\sigma\mu\nu}$
- **Ricci tensor** $R_{\mu\nu}$
- **Scalar curvature (Ricci scalar)** $R$
- **Einstein tensor** $G_{\mu\nu}$

These quantities are fundamental in describing spacetime curvature and form the basis of Einstein's field equations.

---

## 🛠️ Requirements

- **Python 3.x**
- **SymPy** library (symbolic computation)

### Installing required packages

```bash
pip install -r requirements.txt
```

or manually:

```bash
pip install sympy
```

---

## 📂 Project Structure

```
Tensor-calculator/
├── calculate_tensor.py   # Main program file
├── energy_tensor.py      # Helper module (extension)
├── metric.txt            # Example metric file (Schwarzschild)
└── README.md             # Documentation
```

---

## 🚀 Getting Started

### Step 1: Prepare the metric file

The repository includes an example `metric.txt` file with the Schwarzschild metric (black hole). You can use it directly or modify it.

#### Metric file format:

```
# Comments start with #

# First line (without comment): coordinates ; parameters
t, r, theta, phi ; M

# Metric tensor elements: index_i index_j expression
0 0 -(1 - 2*M/r)
1 1 1/(1 - 2*M/r)
2 2 r**2
3 3 r**2 * sin(theta)**2
```

#### Format rules:

| Element | Description | Example |
|---------|-------------|---------|
| **Coordinates** | Variable names separated by commas | `t, r, theta, phi` |
| **Separator** | Semicolon `;` separates coordinates from parameters | `;` |
| **Parameters** | Constants used in the metric (optional) | `M, a, Lambda` |
| **Metric elements** | `i j expression` - indices start from 0 | `0 0 -1` |
| **Comments** | Lines starting with `#` | `# this is a comment` |

#### Available mathematical functions:
- `sin(x)`, `cos(x)`, `tan(x)` - trigonometric functions
- `exp(x)` - exponential function
- `sqrt(x)` - square root
- `**` - exponentiation (e.g., `r**2`)
- `log(x)` - natural logarithm

### Step 2: Set the metric file path

In `calculate_tensor.py`, change the metric file path if needed:

```python
filename = "metric.txt"
```

### Step 3: Run the program

```bash
python calculate_tensor.py
```

---

## 📊 Example Output

The program will display:
1. Loaded coordinates and parameters
2. Non-zero metric tensor components
3. Non-zero Christoffel symbols
4. Non-zero Riemann tensor components
5. Non-zero Ricci tensor components
6. Non-zero Einstein tensor components
7. Scalar curvature

---

## 🔬 How It Works - Algorithm Details

### 1. Loading the metric (`load_metric`)

The function parses the text file and:
- Creates SymPy symbols for coordinates (e.g., $\tau, \psi, \theta, \phi$)
- Creates symbols for parameters (e.g., $a$ - scale factor)
- Builds a dictionary of metric elements

### 2. Computing tensors (`compute_tensors`)

#### a) Metric tensor $g_{\mu\nu}$
A symmetric matrix describing the geometry of spacetime.

#### b) Christoffel symbols (connection coefficients)
Computed using the formula:

$$\Gamma^\sigma_{\mu\nu} = \frac{1}{2} g^{\sigma\lambda} \left( \partial_\mu g_{\nu\lambda} + \partial_\nu g_{\mu\lambda} - \partial_\lambda g_{\mu\nu} \right)$$

Christoffel symbols describe how vectors change during parallel transport.

#### c) Riemann tensor
A measure of spacetime curvature:

$$R^\rho{}_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho_{\nu\sigma} - \partial_\nu \Gamma^\rho_{\mu\sigma} + \Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}$$

#### d) Ricci tensor
Contraction of the Riemann tensor:

$$R_{\mu\nu} = R^\rho{}_{\mu\rho\nu}$$

#### e) Scalar curvature
Contraction of the Ricci tensor:

$$R = g^{\mu\nu} R_{\mu\nu}$$

#### f) Einstein tensor
Appears in Einstein's field equations:

$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} g_{\mu\nu} R$$

### 3. Expression simplification (`custom_simplify`)

The program applies multi-step symbolic simplification:
- `expand` - expression expansion
- `trigsimp` - trigonometric simplification
- `factor` - factorization
- `simplify` - general simplification
- `cancel` - fraction reduction
- `ratsimp` - rational expression simplification

---

## 📐 Metric Examples

### Minkowski metric (flat spacetime)
```
t, x, y, z ;
0 0 -1
1 1 1
2 2 1
3 3 1
```

### Schwarzschild metric (black hole)
```
t, r, theta, phi ; M
0 0 -(1 - 2*M/r)
1 1 1/(1 - 2*M/r)
2 2 r**2
3 3 r**2 * sin(theta)**2
```

### FLRW metric (cosmology)
```
t, r, theta, phi ; a
0 0 -1
1 1 a**2
2 2 a**2 * r**2
3 3 a**2 * r**2 * sin(theta)**2
```

---

## ⚠️ Notes

- Symbolic computations may be time-consuming for complex metrics
- Ensure the metric is properly defined (symmetric, non-singular)
- Indexing starts from 0

---

## 👨‍💻 Author

Project created as part of a Python course at university.

---

## 📚 References

1. Carroll, S. M. - "Spacetime and Geometry: An Introduction to General Relativity"
2. Misner, Thorne, Wheeler - "Gravitation"
3. SymPy Documentation: https://docs.sympy.org/
