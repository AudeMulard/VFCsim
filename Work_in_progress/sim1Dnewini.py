#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 5 09:30:59 2017

@author: aude

New smooth initial conditions
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
Mobility = 1. #ratio of the two viscosities; M_c in Hamouda's paper
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
phi = CellVariable(name=r'$\phi$', mesh=mesh)

#New values
beta = Variable(name='beta', value = beta1 * phi + beta2 * (1.-phi))

#Parameters
#Cahn_number = 0.001
#epsilon = Cahn_number * W

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (3.* PHI*(PHI-1.) + 0.5)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1
eq = (TransientTerm()  == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))
#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
#	phi.setValue(0.)
#	phi.setValue(1., where=x > nx*dx/2)
#    phi.setValue(1-0.5*(1- numerix.tanh((x - nx*dx/2)/(2*numerix.sqrt(M)))))
    phi.setValue(0.5*(1+ numerix.tanh((x-nx*dx/2)/(numerix.sqrt(2)*epsilon))))

initialize(phi)

viewer = Viewer(vars=(phi,), datamin=0., datamax=1.5)

dexp = 1.
elapsed = 0.
duration = 100.
while elapsed < duration:
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt)
    if __name__ == '__main__':
        viewer.plot()
    elif (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3) and elapsed > 10.:
        break


