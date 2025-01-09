import numpy as np
import matplotlib.pyplot as plt

# Zakres temperatur
temperatures = np.linspace(-50, 150, 250)  # od -50°C do 150°C

# Parametry rezystorów
R0 = 0.1 # Rezystancja początkowa w omach

# Nikiel - współczynnik temperaturowy α dla niklu
alpha_nikiel = 0.006  # Współczynnik temperaturowy [1/°C]
R_nikiel = R0 * (1 + alpha_nikiel * temperatures)

# NTC - współczynnik B dla termistora NTC
B_ntc = 3950  # Stała B [K]
T0 = 298.15  # Temperatura odniesienia 25°C w kelwinach
R_ntc = R0 * np.exp(B_ntc * (1 / (temperatures + 273.15) - 1 / T0))

# PT - współczynnik temperaturowy α i β dla platyny (Pt100)
alpha_pt = 0.00385
beta_pt = -5.8e-7
R_pt = R0 * (1 + alpha_pt * temperatures + beta_pt * temperatures**2)

# Konstantan - praktycznie stała rezystancja
R_konstantan = np.full_like(temperatures, R0)

# Tworzenie wykresu
plt.figure(figsize=(12, 8))
plt.plot(temperatures, R_nikiel, label="Nikiel", color="blue")
plt.plot(temperatures, R_ntc, label="NTC", color="red")
plt.plot(temperatures, R_pt, label="PT", color="green")
plt.plot(temperatures, R_konstantan, label="Konstantan", color="orange")

# Opisy osi i tytuł
plt.title("Zależność rezystancji od temperatury dla różnych rezystorów")
plt.xlabel("Temperatura [°C]")
plt.ylabel("Rezystancja [Ω]")
plt.legend()
plt.grid(True)
plt.show()
