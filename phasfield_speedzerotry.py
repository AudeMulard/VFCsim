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
dx = 0.25 #width of controle volume
nx = 1000 #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)
#mesh = Grid2D(dx=dx, dy=dx, nx=nx, ny=nx)#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
Mobility = 1. #ratio of the two viscosities
epsilon =1. #code starts going crazy below epsilon=0.1
l = 1. #this is lambda from Hamouda's paper
M = Mobility * epsilon**2 #M in Hamouda's paper
viscosity1 = viscosity2 * Mobility
permeability1 = permeability2 = 1.
beta1 = - viscosity1 / permeability1
beta2 = - viscosity2 / permeability2

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
#x = mesh.cellCenters[0]
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld=1)
#phi = CellVariable(name=r'$\phi$', mesh=mesh, value =GaussianNoiseVariable(mesh=mesh, mean=0.5, variance=0.01))
#New values
beta = Variable(name='beta', value = beta1 * phi + beta2 * (1.-phi))

#Parameters
#Cahn_number = 0.001
#epsilon = Cahn_number * W




#Cahn-Hilliard equationimport numpy
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1
#D = 10.
#a = 1.
#coeff1=D * a**2 *(1-6*PHI *(1-PHI))
#eq = (TransientTerm() == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(D, epsilon**2)))
U = FaceVariable(mesh=mesh, rank = 1)
eq = (TransientTerm() + ConvectionTerm(coeff=U) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(0.)
	phi.setValue(1., where=x > nx*dx/2)

initialize(phi)

viewer = Viewer(vars=(phi,), datamin=0., datamax=20.)


timeStep = 10.
for i in range(20):
    phi.updateOld()
    res = 1e+10
    while res > 1e-7:
        res = eq.sweep(var=phi, dt=timeStep)

if __name__ == '__main__':
    viewer.plot()

x = mesh.cellCenters[0]    
U.setValue((3.,))
displacement = 1.
velocity1 = 3.
timeStep = 0.1 * dx / velocity1
elapsed = 0.

phi.updateOld()
res = 1e+10


print(U)

while elapsed < displacement/velocity1:
    phi.updateOld()
    res = 1e+10
    while res > 1e-4:
        res = eq.sweep(var=phi, dt=timeStep)
    elapsed +=timeStep
    if __name__ == '__main__':
        viewer.plot()
        print(U)
