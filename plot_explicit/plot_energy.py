import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("solvation_energies.csv")
energy_8A3 = df["8A3"]
t = np.linspace(1,50,1)

plt.plot(energy_8A3, t)
plt.show()
