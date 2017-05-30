#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 09:03:34 2017

@author: aude
"""

from fipy import *

if __name__ == "__main__":
    nx = 1000
else:
    nx = 20
dx=0.25
mesh = Grid1D(nx=nx, dx=0.25)
phi = CellVariable(name=r"$\phi$", mesh=mesh)


x=mesh.cellCenters[0]
phi.setValue(0.)
phi.setValue(1., where=x > nx*dx/2)

if __name__ == "__main__":
    viewer = Viewer(vars=(phi,), datamin=-1., datamax=1.5)
    
PHI = phi.arithmeticFaceValue
D = a = epsilon = 1.
eq = (TransientTerm() == DiffusionTerm(coeff=D * a**2 * (1 - 6 * PHI * (1 - PHI))) - DiffusionTerm(coeff=(D, epsilon**2)))

dexp = -5
elapsed = 0.
duration = 1.
    
while elapsed < duration:
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(phi, dt=dt)
    if __name__ == "__main__":
        viewer.plot()
    elif (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3) and elapsed > 10.:
        break

print (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3)
