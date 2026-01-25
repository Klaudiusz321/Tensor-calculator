# Tensor Calculator - Kalkulator Tensorów dla Ogólnej Teorii Względności

## 📖 Opis projektu

Program służy do obliczania kluczowych obiektów geometrycznych używanych w **Ogólnej Teorii Względności (OTW)**. Na podstawie zadanej metryki czasoprzestrzeni, program automatycznie oblicza:

- **Tensor metryczny** $g_{\mu\nu}$
- **Symbole Christoffela** $\Gamma^\sigma_{\mu\nu}$
- **Tensor Riemanna** $R_{\rho\sigma\mu\nu}$
- **Tensor Ricciego** $R_{\mu\nu}$
- **Skalar krzywizny (skalar Ricciego)** $R$
- **Tensor Einsteina** $G_{\mu\nu}$

Te wielkości są fundamentalne w opisie zakrzywienia czasoprzestrzeni i stanowią podstawę równań pola Einsteina.

---

## 🛠️ Wymagania

- **Python 3.x**
- Biblioteka **SymPy** (obliczenia symboliczne)

### Instalacja wymaganych pakietów

```bash
pip install -r requirements.txt
```

lub ręcznie:

```bash
pip install sympy
```

---

## 📂 Struktura projektu

```
Tensor-calculator/
├── calculate_tensor.py   # Główny plik programu
├── energy_tensor.py      # Moduł pomocniczy (rozszerzenie)
├── metric.txt            # Przykładowy plik z metryką Schwarzschilda
└── README.md             # Dokumentacja
```

---

## 🚀 Instrukcja uruchomienia

### Krok 1: Przygotuj plik z metryką

W repozytorium znajduje się przykładowy plik `metric.txt` z metryką Schwarzschilda (czarna dziura). Możesz go użyć bezpośrednio lub zmodyfikować.

#### Format pliku metryki:

```
# Komentarze zaczynają się od #

# Pierwsza linia (bez komentarza): współrzędne ; parametry
t, r, theta, phi ; M

# Elementy tensora metrycznego: indeks_i indeks_j wyrażenie
0 0 -(1 - 2*M/r)
1 1 1/(1 - 2*M/r)
2 2 r**2
3 3 r**2 * sin(theta)**2
```

#### Zasady formatu:

| Element | Opis | Przykład |
|---------|------|----------|
| **Współrzędne** | Nazwy zmiennych oddzielone przecinkami | `t, r, theta, phi` |
| **Separator** | Średnik `;` oddziela współrzędne od parametrów | `;` |
| **Parametry** | Stałe używane w metryce (opcjonalne) | `M, a, Lambda` |
| **Elementy metryki** | `i j wyrażenie` - indeksy od 0 | `0 0 -1` |
| **Komentarze** | Linie zaczynające się od `#` | `# to jest komentarz` |

#### Dostępne funkcje matematyczne:
- `sin(x)`, `cos(x)`, `tan(x)` - funkcje trygonometryczne
- `exp(x)` - funkcja wykładnicza
- `sqrt(x)` - pierwiastek kwadratowy
- `**` - potęgowanie (np. `r**2`)
- `log(x)` - logarytm naturalny

### Krok 2: Ustaw ścieżkę do pliku metryki

W pliku `calculate_tensor.py` zmień ścieżkę do pliku metryki:

```python
filename = r"C:\Users\sorak\Desktop\Tensor-calculator\metric.txt"
```

Na ścieżkę do Twojego pliku z metryką.

### Krok 3: Uruchom program

```bash
python calculate_tensor.py
```

---

## 📊 Przykładowe wyjście

Program wyświetli:
1. Wczytane współrzędne i parametry
2. Niezerowe komponenty tensora metrycznego
3. Niezerowe symbole Christoffela
4. Niezerowe komponenty tensora Riemanna
5. Niezerowe komponenty tensora Ricciego
6. Niezerowe komponenty tensora Einsteina
7. Skalar krzywizny

---

## 🔬 Jak działa program - szczegółowy opis algorytmu

### 1. Wczytywanie metryki (`wczytaj_metryke`)

Funkcja parsuje plik tekstowy i:
- Tworzy symbole SymPy dla współrzędnych (np. $\tau, \psi, \theta, \phi$)
- Tworzy symbole dla parametrów (np. $a$ - promień wszechświata)
- Buduje słownik elementów metryki

### 2. Obliczanie tensorów (`oblicz_tensory`)

#### a) Tensor metryczny $g_{\mu\nu}$
Macierz symetryczna opisująca geometrię czasoprzestrzeni.

#### b) Symbole Christoffela (symbole połączenia)
Obliczane ze wzoru:

$$\Gamma^\sigma_{\mu\nu} = \frac{1}{2} g^{\sigma\lambda} \left( \partial_\mu g_{\nu\lambda} + \partial_\nu g_{\mu\lambda} - \partial_\lambda g_{\mu\nu} \right)$$

Symbole Christoffela opisują jak wektory zmieniają się przy transporcie równoległym.

#### c) Tensor Riemanna
Miara krzywizny czasoprzestrzeni:

$$R^\rho{}_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho_{\nu\sigma} - \partial_\nu \Gamma^\rho_{\mu\sigma} + \Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}$$

#### d) Tensor Ricciego
Kontrakcja tensora Riemanna:

$$R_{\mu\nu} = R^\rho{}_{\mu\rho\nu}$$

#### e) Skalar krzywizny
Kontrakcja tensora Ricciego:

$$R = g^{\mu\nu} R_{\mu\nu}$$

#### f) Tensor Einsteina
Występuje w równaniach pola Einsteina:

$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} g_{\mu\nu} R$$

### 3. Upraszczanie wyrażeń (`custom_simplify`)

Program stosuje wieloetapowe upraszczanie symboliczne:
- `expand` - rozwinięcie wyrażeń
- `trigsimp` - uproszczenia trygonometryczne
- `factor` - faktoryzacja
- `simplify` - ogólne uproszczenie
- `cancel` - redukcja ułamków
- `ratsimp` - uproszczenie wyrażeń wymiernych

---

## 📐 Przykłady metryk

### Metryka Minkowskiego (płaska czasoprzestrzeń)
```
t, x, y, z ;
0 0 -1
1 1 1
2 2 1
3 3 1
```

### Metryka Schwarzschilda (czarna dziura)
```
t, r, theta, phi ; M
0 0 -(1 - 2*M/r)
1 1 1/(1 - 2*M/r)
2 2 r**2
3 3 r**2 * sin(theta)**2
```

### Metryka FLRW (kosmologia)
```
t, r, theta, phi ; a
0 0 -1
1 1 a**2
2 2 a**2 * r**2
3 3 a**2 * r**2 * sin(theta)**2
```

---

## ⚠️ Uwagi

- Obliczenia symboliczne mogą być czasochłonne dla złożonych metryk
- Upewnij się, że metryka jest poprawnie zdefiniowana (symetryczna, nieosobliwa)
- Indeksowanie zaczyna się od 0

---

## 👨‍💻 Autor

Projekt wykonany w ramach kursu Python na studiach.

---

## 📚 Bibliografia

1. Carroll, S. M. - "Spacetime and Geometry: An Introduction to General Relativity"
2. Misner, Thorne, Wheeler - "Gravitation"
3. Dokumentacja SymPy: https://docs.sympy.org/
