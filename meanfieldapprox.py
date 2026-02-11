# Graphically solving the mean-field approximation of the Ising model 
# from Huang's Introduction to Statistical Physics.

import matplotlib.pyplot as plt
import numpy as np


# Self-consistent equation m = tanh(\beta (h + kappa \epsilon m))
k = 4 # \kappa number of nearest neighbours
E = 0.4 # \epsilon energy of spin > 0 for ferromagnet
beta = 1 # thermal constant

m = np.linspace(-1.2,1.2)

# If h, the external magnetic field is non-zero
h = 0.4
lhs = m
rhs = np.tanh(beta*(h+k*E*m))

plt.figure(1)
plt.title(r'Non-zero external field $h\neq 0$')
plt.plot(m, lhs, label=r"LHS $m$")
plt.plot(m, rhs, label=r"RHS $\beta(h+\kappa\epsilon m)$")
plt.grid(True)
plt.legend()
plt.xlabel('m')
plt.savefig("./isingnonzerofield.jpg")

# If h = 0, slope greater than one
h=0
rhs = np.tanh(beta*(h+k*E*m))

plt.figure(2)
plt.title(r'No external field $h = 0$ Case B')
plt.plot(m, lhs, label=r"LHS $m$")
plt.plot(m, rhs, label=r"RHS $\beta(h+\kappa\epsilon m)$")
plt.grid(True)
plt.legend()
plt.xlabel('m')
plt.savefig("./isingzerofield_largeslope.jpg")

# If h = 0, slope greater than one
h=0
E = 0.2
rhs = np.tanh(beta*(h+k*E*m))

plt.figure(3)
plt.title(r'No external field $h = 0$ Case B')
plt.plot(m, lhs, label=r"LHS $m$")
plt.plot(m, rhs, label=r"RHS $\beta(h+\kappa\epsilon m)$")
plt.grid(True)
plt.legend()
plt.xlabel('m')
plt.savefig("./isingzerofield_smallslope.jpg")



