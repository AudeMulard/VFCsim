#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:48:31 2017

@author: aude
"""

from fipy import *

nx = 50
dx = 1.
mesh = Grid1D(nx=nx, dx=dx)

phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=0.)

D = 1.

valueLeft = 1
valueRight = 0

phi.constrain(valueRight, mesh.facesRight)
phi.constrain(valueLeft, mesh.facesLeft)

eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D)

timeStepDuration = 0.9 * dx**2 / (2 * D)
steps = 100

phiAnalytical = CellVariable(name="analytical value", mesh=mesh)

if __name__ == '__main__':
    viewer = Viewer(vars=(phi, phiAnalytical), datamin=0., datamax=1.)
    viewer.plot()

x = mesh.cellCenters[0]
t = timeStepDuration * steps

try:
    from scipy.special import erf
    phiAnalytical.setValue(1 - erf(x / (2 * numerix.sqrt(D * t))))
except ImportError:
    print("The Scipy library is not available to test the solution to the transient diffusion equation")

for step in range(steps):
    eqX.solve(var=phi,
              dt=timeStepDuration)
    if __name__ == '__main__':
        viewer.plot()

print(phi.allclose(phiAnalytical, atol = 7e-4))

if __name__ == '__main__':
    raw_input("Explicit transient diffusion. Press <return> to proceed")