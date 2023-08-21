# Code to simulate the two-capacitor circuit
# Gnd --)|--o--\/\/\/--o--|(---Gnd
#       C1  V1    R    V2 C2
# One of the capacitance (C1) is varied with a user defined function
# Steady state is V1 = V2
# Intuition is that (Q1=C1*V1) + (Q2=C2*V2) is conserved

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd

# Defining the US-active duration
tbegin = 1
tdur = 2.75

def c1(t):
    # Defining the sinusoidal capacitance for US-active duration
    if t > tbegin and t < tbegin + tdur:
        return 1 + 0.5 * np.sin( 2 * np.pi * t )
    else:
        return 1


def dc1(t):
    # Defining the dc1/dt function for c1
    if t > tbegin and t < tbegin + tdur:
        return np.pi * np.cos( 2 * np.pi * t )
    else:
        return 0

# Rest of the circuit elements
c2 = 1
R = 1 # This value gives best match to NEURON's output for Ra = 1e6 Ohm.cm. why is that?

# Defining the derivative function
def dX(t, X):
    # derivative function for v1 and v2
    v1, v2 = X

    dv1 = (1/c1(t)) * ( -(dc1(t) + 1/R) * v1 + v2/R )
    dv2 = (v1 - v2) / (c2*R)

    return [dv1, dv2]


# Initialization

v1_0 = -65
v2_0 = -65

tspan = [0, 5]

X0 = [v1_0, v2_0]

# Solving the ODE to obtain v1 and v2
sol = solve_ivp( dX, tspan, X0, method = 'RK45', dense_output = True)

# Displaying the results
dt = 1/40 # 40 sample per Capacitance cycle 
tsim = np.arange(0, 5, dt )
v1, v2 = sol.sol(tsim)

# Loading the NEURON output for comparison
v1_nrn = pd.read_csv("v1_nrn.dat", sep="\t", names= ['t', 'v'])
v2_nrn = pd.read_csv("v2_nrn.dat", sep="\t", names= ['t', 'v'])

plt.figure("V1")
plt.plot(tsim, v1, linewidth = 2, label = "Python")
plt.plot(v1_nrn.t, v1_nrn.v, linewidth = 2, label = 'NEURON')
plt.xlabel("Time")
plt.ylabel("Potential")
plt.ylim([-150, 0])
plt.legend()
plt.savefig("v1_nrn_vs_py.png")


plt.figure("V2")
plt.plot(tsim, v2, linewidth = 2, label = "Python")
plt.plot(v2_nrn.t, v2_nrn.v, linewidth = 2, label = 'NEURON')
plt.xlabel("Time")
plt.ylabel("Potential")
plt.ylim([-150, 0])
plt.legend()
plt.savefig("v2_nrn_vs_py.png")


# Plotting the Total charge
c1_py = np.array( [c1(t) for t in tsim] )
Q = v1*c1_py + v2 # v1c1 + v2c2

c1_nrn =  np.array( [c1(t) for t in v1_nrn.t] )
Q_nrn = v1_nrn.v * c1_nrn + v2_nrn.v

plt.figure("Q")
plt.plot(tsim, Q, linewidth = 2, label="Python")
plt.plot(v1_nrn.t, Q_nrn, linewidth = 2, label = "NEURON")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Charge")
plt.savefig("Q_nrn_vs_py.png")

#plt.show()

