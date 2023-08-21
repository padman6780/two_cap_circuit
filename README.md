# two_cap_circuit
Simulation of two-capacitor circuit using NEURON, Python, and MATLAB

Study conducted by [Michael Hines](https://github.com/nrnhines)

Code to simulate the two-capacitor circuit

```
 Gnd --)|--o--\/\/\/--o--|(---Gnd
       C1  V1    R    V2 C2
```
One of the capacitance (C1) is varied with a user defined function
Steady state is V1 = V2
Intuition is that (Q1=C1*V1) + (Q2=C2*V2) is conserved

The neuron implementation is available [here](https://github.com/nrnhines/dcmdt.git)
