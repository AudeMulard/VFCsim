#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:54:31 2017

@author: aude
"""

from fipy import *

if __name__ == '__main__':
    nx = ny = 1000
else:
    nx = ny = 20

mesh = Grid2D(nx=nx, ny=ny, dx=0.25, dy=0.25)
phi = CellVariable(name=r"$\phi$", mesh=mesh)

phi.setValue(GaussianNoiseVariable(mesh=mesh,
                                   mean=0.5,
                                   variance=0.01))

if __name__ == '__main__':
    veiwer = Viewer(vars=(phi,), datamin=0., datamax=1.)

PHI = phi.arithmeticFaceValue
D = a = epsilon = 1.
eq = (TransientTerm() == DiffusionTerm(coeff=D * a**2 * (1 - 6 * PHI * (1 - PHI))) - DiffusionTerm(coeff=(D, epsilon**2)))

dexp = -5
elapsed = 0.
if __name__ == '__main__':
    duration = 1000.
else:
    duration = 1000.

while elapsed < duration:
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(phi, dt=dt)
    if __name__ == '__main__':
        veiwer.plot()
    elif (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3):
        break

print(max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3)
