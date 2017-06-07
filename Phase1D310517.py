#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 17:48:36 2017

@author: aude
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
dx = 0.25 #width of controle volume
nx = 1000 #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
Mobility = 0.1 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon =1. #code starts going crazy below epsilon=0.1
l = 1. #this is lambda from Hamouda's paper
M = Mobility * epsilon**2 #M in Hamouda's paper
viscosity1 = viscosity2 * Mobility
permeability1 = permeability2 = 1.
beta1 = viscosity1 / permeability1
beta2 = viscosity2 / permeability2

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh)

#New values
beta = CellVariable(mesh=mesh, name='beta', value = beta2 * phi + beta1 * (1.-phi))

#Parameters
#Cahn_number = 0.001
#epsilon = Cahn_number * W

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1
eq = (TransientTerm()  == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))
#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(0.)
	phi.setValue(1., where=x > nx*dx*2/3)

initialize(phi)

viewer = Viewer(vars=(phi,), datamin=-0.5, datamax=1.5)

dexp = 1.
elapsed = 0.
duration = 1500.
while elapsed < duration:
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt)
    if __name__ == '__main__':
        viewer.plot()


print((max(phi.globalValue) > 1.00001) and (min(phi.globalValue) < -0.00001))
raw_input("pause")


