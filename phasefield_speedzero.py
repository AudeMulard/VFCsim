#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:48:05 2017

@author: am2548
"""

from fipy import *
#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.02 * W #width of controle volume
nx = L / dx #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
Mobility = 1. #ratio of the two viscosities
viscosity1 = viscosity2 * Mobility
permeability1 = permeability2 = 1.
beta1 = - viscosity1 / permeability1
beta2 = - viscosity2 / permeability2

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh)

#New values
beta = Variable(name='beta')
beta.setValue = beta1 * phi + beta2 * (1-phi)

#Parameters
Cahn_number = 0.001
epsilon = Cahn_number * W
M = Mobility * epsilon**2
l = 1.

#Cahn-Hilliard equationimport numpy
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (3 * PHI**2 - 2 * PHI + 1/2)
eq = (TransientTerm() == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(0.5)
	phi.setValue(1, where=x > L/2)

initialize(phi)
viewer = Viewer(vars=(phi,), datamin=-1., datamax=2.)

for i in range(500):
    eq.solve(var=phi, dt = 1e-6)
    if __name__ == '__main__':
        viewer.plot()
