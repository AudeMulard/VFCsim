#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:21:51 2017

@author: aude
"""



from fipy import *
import random
from math import sqrt


U = 0.8
Mobility = 0.2 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 0.6 #code starts going crazy below epsilon=0.1
l = 0.3 #this is lambda from Hamouda's paper
duration = 1500. #stabilisation phase
sweeps = 100 #stabilisation vitesse
alpha=0.1

#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.15 #width of controle volume
nx = 300 #number of controle volume
dy = 0.5
ny = 100
mesh = Grid2D(dx=dx, nx=nx, dy=dy, ny=ny)

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
xVelocity = CellVariable(mesh=mesh, name='X Velocity', value=U)
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
    phi.setValue(0.)
    for i in range(50):
        a = random.gauss(0.2, 0.01)
        phi.setValue(1., where=(x > nx*dx * a ) & (y<2*(i+1)*dy) & (y>2*(i*dy)))


initialize(phi)


beta = CellVariable(mesh=mesh, name=r'$\beta$', value = beta2 * phi + beta1 * (1.-phi))
#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0] - ImplicitSourceTerm(alpha/(numerix.sqrt(xVelocity*xVelocity+yVelocity*yVelocity))))
yVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[1] - ImplicitSourceTerm(alpha/(numerix.sqrt(xVelocity*xVelocity+yVelocity*yVelocity))))


coeff = 1./ beta.arithmeticFaceValue
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable


#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi), datamin=0., datamax=1.)
viewer2 = Viewer(vars = (xVelocity), datamin=0., datamax=1.)
viewer3 = Viewer(vars = (yVelocity), datamin=0., datamax=1.)
viewer4 = Viewer(vars = (pressure), datamin=0., datamax=250.)



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
    eq.solve(var=phi, dt = dt)
    if __name__ == '__main__':
        viewer.plot()

 

#Pressure and velocity

xVelocity.constrain(U, mesh.facesLeft)
yVelocity.constrain(0, mesh.facesTop | mesh.facesBottom)
pressureCorrection.constrain(0., mesh.facesRight)

pressureRelaxation = 0.8
velocityRelaxation = 0.5

for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred value
#    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
#    xmat = xVelocityEq.matrix
    yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation)
    ##update the ap coefficient from the matrix diagonal
#    ap[:] = xmat.takeDiagonal()
    ##update the face velocities based on starred values with the Rhi-Chow correction
    #cell pressure gradient
    presgrad = pressure.grad
    #face pressure gradient
    facepresgrad = _FaceGradVariable(pressure)
    #
    velocity[0] = xVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
    velocity[1] = yVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[1].arithmeticFaceValue-facepresgrad[1])
    velocity[0, mesh.facesLeft.value] = U
#    velocity[0, mesh.facesRight.value] = U
    ##solve the pressure correction equation
    pressureCorrectionEq.cacheRHSvector()
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    ## update the pressure using the corrected value
    pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
    ## update the velocity using the corrected pressure
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / beta)
    yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / beta)
    xVelocity[0]=U
#    xVelocity[nx-1]=U




displacement = 100.
timeStep = 0.6 * dx / U #less than one space step per time step
elapsed = 0.

while elapsed < displacement/U:
    phi.updateOld()
    res = 1e+10
    while res > 1e-6:
        res = eq.sweep(var=phi, dt=timeStep)
    beta.setValue(beta2 * phi + beta1 * (1.-phi))
    for sweep in range(sweeps):
        ##Solve the Stokes equations to get starred value
#        xVelocityEq.cacheMatrix()
        xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
#        xmat = xVelocityEq.matrix
        yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation)
        ##update the ap coefficient from the matrix diagonal
#        ap[:] = xmat.takeDiagonal()
        ##update the face velocities based on starred values with the Rhi-Chow correction
        #cell pressure gradient
        presgrad = pressure.grad
        #face pressure gradient
        facepresgrad = _FaceGradVariable(pressure)
        #
        velocity[0] = xVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
        velocity[1] = yVelocity.arithmeticFaceValue + 1. / beta.arithmeticFaceValue * (presgrad[1].arithmeticFaceValue-facepresgrad[1])
        velocity[0, mesh.facesLeft.value] = U
#        velocity[0, mesh.facesRight.value] = U
        ##solve the pressure correction equation
        pressureCorrectionEq.cacheRHSvector()
        pres = pressureCorrectionEq.sweep(var=pressureCorrection)
        rhs = pressureCorrectionEq.RHSvector
        ## update the pressure using the corrected value
        pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
        ## update the velocity using the corrected pressure
        xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / beta)
        yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / beta)
        xVelocity[0]=U
#        xVelocity[nx-1]=U
    elapsed +=timeStep
    viewer.plot(filename="phi%d.png" % elapsed)
    viewer2.plot(filename="XVelocity%d.png" % elapsed)
    viewer4.plot(filename="pressure%d.png" % elapsed)
    viewer3.plot(filename="YVelocity%d.png" % elapsed)
    TSVViewer(vars=(phi, xVelocity, yVelocity, pressure,beta)).plot(filename="essaidonne%d.tsv" % elapsed)
    print(elapsed)


