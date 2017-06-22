#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 10:55:13 2017

@author: aude
"""

from fipy import *


U = 0.8
Mobility = 0.2 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon =1. #code starts going crazy below epsilon=0.1
l = 0.01 #this is lambda from Hamouda's paper
alpha=0.001
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

#New values


#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1

eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))



#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

#Phase
x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(0.)
	phi.setValue(1., where=x > nx*dx/2)

initialize(phi)

phi.faceGrad.constrain([0], mesh.facesRight)


beta = CellVariable(mesh=mesh, name='beta', value = beta2 * phi + beta1 * (1.-phi))

#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0] - ImplicitSourceTerm(alpha/U))

ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence -1+phi
#pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence
#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue
#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi), datamin=-1., datamax=2.)
viewer2 = Viewer(vars = (xVelocity), datamin=-1., datamax=3.)
viewer3 = Viewer(vars=(pressure), datamin=0., datamax=250.)


#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------


dexp = 1.
elapsed = 0.
duration = 1500.
while elapsed < duration:
    phi.updateOld()
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt)
    if __name__ == '__main__':
        viewer.plot()

 


#Pressure and velocity
pressureRelaxation = 0.8
velocityRelaxation = 0.5
xVelocity.constrain(U, mesh.facesLeft)
pressureCorrection.constrain(0., mesh.facesRight)

sweeps = 41
for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred value
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix
    ##update the ap coefficient from the matrix diagonal
    ap[:] = xmat.takeDiagonal()
    #
    ##update the face velocities based on starred values with the Rhi-Chow correction
    #cell pressure gradient
    presgrad = pressure.grad
    #face pressure gradient
    facepresgrad = _FaceGradVariable(pressure)
    #
    velocity[0] = xVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
    velocity[0, mesh.facesLeft.value] = U
    velocity[0, mesh.facesRight.value] = U
    ##solve the pressure correction equation
    pressureCorrectionEq.cacheRHSvector()
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    ## update the pressure using the corrected value
    pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
    ## update the velocity using the corrected pressure
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)
    xVelocity[0]=U
    if sweep%10 == 0:
        viewer2.plot()
        viewer3.plot()
#        print 'sweep:',sweep,', x residual:',xres, ', p residual:',pres, ', continuity:',max(abs(rhs))



displacement = 50.
timeStep = 0.8 * dx / U #less than one space step per time step
elapsed = 0.

while elapsed < displacement/U:
    phi.updateOld()
    res = 1e+10
    while res > 1e-6:
        res = eq.sweep(var=phi, dt=timeStep)
    beta.setValue(beta2 * phi + beta1 * (1.-phi))
    for sweep in range(sweeps):
        ##Solve the Stokes equations to get starred value
        xVelocityEq.cacheMatrix()
        xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
        xmat = xVelocityEq.matrix
        ##update the ap coefficient from the matrix diagonal
        ap[:] = xmat.takeDiagonal()
        #
        ##update the face velocities based on starred values with the Rhi-Chow correction
        #cell pressure gradient
        presgrad = pressure.grad
        #face pressure gradient
        facepresgrad = _FaceGradVariable(pressure)
        #
        velocity[0] = xVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
        velocity[0, mesh.facesLeft.value] = U
        velocity[0, mesh.facesRight.value] = U
        ##solve the pressure correction equation
        pressureCorrectionEq.cacheRHSvector()
        pres = pressureCorrectionEq.sweep(var=pressureCorrection)
        rhs = pressureCorrectionEq.RHSvector
        ## update the pressure using the corrected value
        pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
        ## update the velocity using the corrected pressure
        xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)
        xVelocity[0]=U
    elapsed +=timeStep
    viewer.plot()
    viewer2.plot()
    viewer3.plot()


raw_input("pause")