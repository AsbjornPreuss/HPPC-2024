import numpy as np
import matplotlib.pyplot as plt
x, y, z, t = np.loadtxt('sir_out.txt', skiprows=5, unpack=True)
plt.plot(t, x, label="Infected")
plt.plot(t, y, label="Recovered")
plt.plot(t, z, label="Susceptible")
plt.xlabel("Time [days]")
plt.ylabel("Number of people [#]")
plt.legend()
plt.tight_layout()
stepsize = input("Stepsize of plot [1]: ")
plt.savefig(f"../Report/Images/SIR_plot_dt{stepsize}.png", dpi=400)
plt.show()

