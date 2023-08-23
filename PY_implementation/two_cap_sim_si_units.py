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
w = 1e3 # Frequency in Hz
C0 = np.pi * 1e-12 # Base capacitance in F

def c1(t):
    # Defining the sinusoidal capacitance for US-active duration. 
    # t in seconds
    if t > tbegin and t < tbegin + tdur:
        return C0 + 0.5 * np.sin( 2 * np.pi * w * t )
    else:
        return C0


def dc1(t):
    # Defining the dc1/dt function for c1
    # t in ms again
    if t > tbegin and t < tbegin + tdur:
        return np.pi * w * np.cos( 2 * np.pi * w * t )
    else:
        return 0

# Rest of the circuit elements
c2 = C0
R = 4e9/np.pi # Resistance in Ohm 

# Defining the derivative function
def dX(t, X):
    # derivative function for v1 and v2
    v1, v2 = X

    dv1 = (1/c1(t)) * ( -(dc1(t) + 1/R) * v1 + v2/R )
    dv2 = (v1 - v2) / (c2*R)

    return [dv1, dv2]


# Initialization

# Voltages in mV
v1_0 = -65e-3
v2_0 = -65e-3 

tspan = [0, 5e-3] # time in seconds

X0 = [v1_0, v2_0]

# Solving the ODE to obtain v1 and v2
sol = solve_ivp( dX, tspan, X0, method = 'RK45', dense_output = True, rtol = 1e-9 )

# Displaying the results
dt = 1/(w*40) # 40 sample per Capacitance cycle 
tsim = np.arange(0, 5e-3, dt ) # tsim in seconds
v1, v2 = sol.sol(tsim) # voltages in V

plt.figure("SI unit solution")
plt.plot(tsim, v1, label = "V1")
plt.plot(tsim, v2, label = "V2")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.legend()
plt.savefig("soln_si_units.png")
# Calculating the Total charge in nC/cm2
c1_py = np.array( [c1(t) for t in tsim] )
Q = v1*c1_py + v2*c2 # v1c1 + v2c2; will get the result in Coulombs
area = 2 * np.pi * 10e-4 * 10e-4 # area in cm2

# converting voltage to mV
v1 *= 1e3
v2 *= 1e3 

# converting tsim to ms 
tsim *= 1e3

# converting Charge units
Q = Q * 1e9 / area # Charge converted into nC/cm2

# Loading the NEURON output for comparison
v1_nrn = pd.read_csv("v1_nrn.dat", sep="\t", names= ['t', 'v'])
v2_nrn = pd.read_csv("v2_nrn.dat", sep="\t", names= ['t', 'v'])

def c1_nrnfn(t):
    # Defining the sinusoidal capacitance in NEURON units
    # t in ms
    if t > tbegin and t < tbegin + tdur:
        return 1 + 0.5 * np.sin( 2 * np.pi * t )
    else:
        return 1

c1_nrn =  np.array( [c1_nrnfn(t) for t in v1_nrn.t] )
Q_nrn = v1_nrn.v * c1_nrn + v2_nrn.v

plt.figure("V1")
plt.plot(tsim, v1, linewidth = 2, label = "Python")
plt.plot(v1_nrn.t, v1_nrn.v, linewidth = 2, label = 'NEURON')
plt.xlabel("Time")
plt.ylabel("V1")
plt.ylim([-150, 0])
plt.legend()
plt.savefig("v1_nrn_vs_py.png")


plt.figure("V2")
plt.plot(tsim, v2, linewidth = 2, label = "Python")
plt.plot(v2_nrn.t, v2_nrn.v, linewidth = 2, label = 'NEURON')
plt.xlabel("Time")
plt.ylabel("V2")
plt.ylim([-150, 0])
plt.legend()
plt.savefig("v2_nrn_vs_py.png")
plt.figure("Q")
plt.plot(tsim, Q, linewidth = 2, label="Python")
plt.plot(v1_nrn.t, Q_nrn, linewidth = 2, label = "NEURON")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Charge")
plt.savefig("Q_nrn_vs_py.png")

plt.show()

