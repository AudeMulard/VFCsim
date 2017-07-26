#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:21:51 2017

@author: aude

Augmentation de la vitesse pour diminuer la tension de surface
"""



from fipy import *
import random
import csv
import os, sys
import numpy


U = 0.
Mobility = 1 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 0.35 #code starts going crazy below epsilon=0.1
l = 0.02 #this is lambda from Hamouda's paper
duration = 0. #stabilisation phase
sweeps = 41 #stabilisation vitesse
startpoint=0.1
k=0.05

#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.15 #width of controle volume
nx = 150 #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

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

velocity = FaceVariable(mesh=mesh, rank=1)

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld=1.)

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1

eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))




#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------
#phi.faceGrad.constrain([0], mesh.facesRight | mesh.facesLeft)
#Phase
x = mesh.cellCenters[0]

def initialize(phi):
    phi.setValue(0.5*(1+numerix.tanh((x-nx*dx/2)/(2*epsilon))))



initialize(phi)

beta = CellVariable(mesh=mesh, name=r'$\beta$', value = beta2 * phi + beta1 * (1.-phi))
#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])

coeff = 1./ beta.arithmeticFaceValue
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence + k*(1.-phi)

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable


#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi), datamin=-1., datamax=2.)
viewer2 = Viewer(vars = (xVelocity), datamin=0., datamax=4.)

viewer4 = Viewer(vars = (pressure), datamin=0., datamax=80.)



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
        viewer.plot()

 



#Pressure and velocity

#xVelocity.constrain(U, mesh.facesLeft)

pressureCorrection.constrain(0., mesh.facesRight)

pressureRelaxation = 0.8
velocityRelaxation = 0.5


for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred value
#    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
#    xmat = xVelocityEq.matrix
    ##update the ap coefficient from the matrix diagonal
#    ap[:] = xmat.takeDiagonal()
    ##update the face velocities based on starred values with the Rhi-Chow correction
    #cell pressure gradient
    presgrad = pressure.grad
    #face pressure gradient
    facepresgrad = _FaceGradVariable(pressure)
    #
    velocity[0] = xVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
#    velocity[0, mesh.facesLeft.value] = U
#    velocity[0, mesh.facesRight.value] = U
    ##solve the pressure correction equation
    pressureCorrectionEq.cacheRHSvector()
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    ## update the pressure using the corrected value
    pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
    ## update the velocity using the corrected pressure
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / beta)
#    xVelocity[0]=U
#    xVelocity[nx-1]=U
    viewer2.plot()
    viewer4.plot()


displacement = 90.
timeStep = 0.8 * dx #less than one space step per time step
elapsed = 0.

while elapsed < displacement/k:
    phi.updateOld()
    res = 1e+10
    while res > 1e-6:
        res = eq.sweep(var=phi, dt=timeStep, solver=LinearGMRESSolver())
    beta.setValue(beta2 * phi + beta1 * (1.-phi))
#    raw_input("pause")
    for sweep in range(sweeps):
        ##Solve the Stokes equations to get starred value
#        xVelocityEq.cacheMatrix()
        xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
#        xmat = xVelocityEq.matrix
        ##update the ap coefficient from the matrix diagonal
#        ap[:] = xmat.takeDiagonal()
        ##update the face velocities based on starred values with the Rhi-Chow correction
        #cell pressure gradient
        presgrad = pressure.grad
        #face pressure gradient
        facepresgrad = _FaceGradVariable(pressure)
        #
        velocity[0] = xVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
#        velocity[0, mesh.facesLeft.value] = U
#        velocity[0, mesh.facesRight.value] = U
        ##solve the pressure correction equation
        pressureCorrectionEq.cacheRHSvector()
        pres = pressureCorrectionEq.sweep(var=pressureCorrection)
        rhs = pressureCorrectionEq.RHSvector
        ## update the pressure using the corrected value
        pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
        ## update the velocity using the corrected pressure
        xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / beta)
#        xVelocity[0]=U
#        xVelocity[nx-1]=U
    elapsed +=timeStep
    viewer.plot(filename='phi%d_' % elapsed +'.png')
    viewer2.plot(filename='XVelocity%d_' % elapsed +'.png')
    viewer4.plot(filename='pressure%d_' % elapsed +'.png')
    TSVViewer(vars=(phi, xVelocity, pressure,beta)).plot(filename='essaidonne%d_'% elapsed +'.tsv')
    print(elapsed)


raw_input("pause")
