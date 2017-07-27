#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:21:51 2017

@author: aude

Correction du code avec le bon coefficient a
"""



from fipy import *
import random
import csv
import os, sys
import numpy



U = 0.8
Mobility = 0.2 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 0.5 #code starts going crazy below epsilon=0.1
l = 0.1 #this is lambda from Hamouda's paper
duration = 150000. #stabilisation phase
sweeps = 41 #stabilisation vitesse


#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.25 #width of controle volume
nx = 150 #number of controle volume
dy = 1.
ny = 60
mesh = Grid2D(dx=dx, nx=nx, dy=dy, ny=ny)
startpoint=0.5*nx*dx

parameters=csv.writer(open('parameters_'+sys.argv[1]+'.csv','w'), delimiter=' ', quotechar='|')
parameters.writerow(['U']+['Mobility']+['epsilon']+['l']+['dx']+['nx']+['dy']+['ny'])
parameters.writerow([U]+[Mobility]+[epsilon]+[l]+[dx]+[nx]+[dy]+[ny])

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
M = Mobility * epsilon**2 #M in Hamouda's paper
viscosity1 = viscosity2 * Mobility
permeability1 = permeability2 = 1.
beta1 = viscosity1 / permeability1
beta2 = viscosity2 / permeability2

#Variable of the fluids
pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X Velocity')
yVelocity = CellVariable(mesh=mesh, name='Y Velocity')
velocity = FaceVariable(mesh=mesh, rank=1)

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld=1.)

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1

eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))




#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------
phi.faceGrad.constrain([0], mesh.facesRight)
#Phase
x = mesh.cellCenters[0]
y = mesh.cellCenters[1]

def initialize(phi):
#    phi.setValue(0.)
#    phi.setValue(1., where=(x > 0.2*nx*dx +numerix.sin(3*y)))
#    phi.setValue(1-0.5*(1-numerix.tanh((x-nx*dx/2)/(2*numerix.sqrt(M*2*epsilon**2/l)))))
    for i in range(ny):
        a = numpy.random.normal(startpoint, 0.005*nx*dx)
#        a = 0.1*nx*dx + 0.15*(numerix.sin(0.6*numerix.pi/2*(i+3)*dy)+numerix.sin(4*numerix.pi/2*i*dy)+numerix.sin(2*numerix.pi/2*i*dy+numerix.pi/2))
#        phi.setValue(1-0.5*(1-numerix.tanh((x-a*nx*dx)/(2*numerix.sqrt(M*2*epsilon**2/l)))), where=(y<(i+1)*dy) & (y>(i*dy)))
        phi.setValue(0.5*(1+numerix.tanh((x-a)/(2*epsilon))), where=(y<(i+1)*dy) & (y>(i*dy)))



initialize(phi)


beta = CellVariable(mesh=mesh, name=r'$\beta$', value = beta2 * phi + beta1 * (1.-phi))
#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])
yVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[1])

coeff = 1./ beta.arithmeticFaceValue
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable


#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi), datamin=0., datamax=1.)
viewer2 = Viewer(vars = (xVelocity), datamin=0.5, datamax=1.)
viewer3 = Viewer(vars = (yVelocity), datamin=0., datamax=0.2)
viewer4 = Viewer(vars = (pressure), datamin=0., datamax=nx*dx*U)



#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------

#Phase

dexp = 1.
elapsed = 0.

while elapsed < duration:
    phi.updateOld()
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt, solver=LinearGMRESSolver())
    if __name__ == '__main__':
        viewer.plot(filename='phi%d_' % elapsed +sys.argv[1]+'.png')


raw_input("pause")
