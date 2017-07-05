#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:56:58 2017

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
nx = 100 #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
Mobility = 0.1 #ratio of the two viscosities
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
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld=1)
#New values
beta = CellVariable(mesh=mesh, name='beta', value = beta2 * phi + beta1 * (1-phi))
#beta.setValue = beta1 * phi + beta2 * (1-phi)

#Parameters
#Cahn_number = 0.001
epsilon = 1.
M = Mobility * epsilon**2
l = 1.
fluxRight=1.
#phi.constrain(1., mesh.facesRight)
#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------

xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])


ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
X = mesh.cellCenters[0]
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

#Phase
x = mesh.cellCenters[0]
def initialize(phi):
#    phi.setValue(GaussianNoiseVariable(mesh=mesh, mean=0.5, variance=0.01), where=(x > nx*dx/2-epsilon/2) | (x < nx*dx/2+epsilon/2))
#    phi.setValue(1., where=x > nx*dx/2)
#    phi.setValue(0., where=x < nx*dx/2)
    phi.setValue(1-0.5*(1- numerix.tanh((x - nx*dx/2)/(2*numerix.sqrt(M*2*epsilon**2/l)))))
    
initialize(phi)
phi.faceGrad.constrain([0], mesh.facesRight)

#Velocity and pressure
Q = 1. #rate of injection
#U = Q / (b*W)
U = 0.8 #if more, it gets unstable, I should change the time step
#xVelocity.constrain(U, where=mesh.facesRight | mesh.facesLeft)
X = mesh.faceCenters
#pressureCorrection.constrain(0., mesh.facesLeft)

#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi,), datamin=-1., datamax=2.)
viewer2 = Viewer(vars = (xVelocity,), datamin=-1., datamax=3.)


#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------

#Phase
dexp = 1.
elapsed = 0.
duration=100.
while elapsed < duration:
    phi.updateOld()
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt, solver=LinearGMRESSolver())
    if __name__ == '__main__':
        viewer.plot()


