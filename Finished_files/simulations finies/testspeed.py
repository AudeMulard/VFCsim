#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:48:38 2017
Simulation 1D
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
beta = CellVariable(mesh=mesh, name='beta', value = beta2 * phi + beta1 * (1.-phi))

#Parameters
#Cahn_number = 0.001
#epsilon = Cahn_number * W

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
## blows up when mobility is between 2.2 and 2.3 and 0.7 and 0.8, while l and epsilon=1
eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])

ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
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
	phi.setValue(0.)
	phi.setValue(1., where=x > nx*dx/2)

initialize(phi)

phi.faceGrad.constrain([0], mesh.facesRight)



#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi,), datamin=-1., datamax=2.)
viewer2 = Viewer(vars = (xVelocity), datamin=-1., datamax=3.)


#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------
"""
#Phase
timeStep = 10.
for i in range(50):
    phi.updateOld()
    res = 1e+10
    while res > 1e-7:
        res = eq.sweep(var=phi, dt=timeStep)
    if __name__ == '__main__':
        viewer.plot()
"""

dexp = 1.
elapsed = 0.
duration = 1500.
while elapsed < duration:
    phi.updateOld()
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    res = eq.sweep(var=phi, dt = dt)
    if __name__ == '__main__':
        viewer.plot()
#    print(res)
 

viewer.plot(filename="myImage")
#Pressure and velocity
#Velocity and pressure
Q = 1. #rate of injection
#U = Q / (b*W)
U = 1. #if more, it gets unstable, I should change the time step
X = mesh.faceCenters
#pressureCorrection.constrain(0., mesh.facesLeft)
pressureRelaxation = 0.8
velocityRelaxation = 0.5
xVelocity.constrain(U, mesh.facesLeft)
pressureCorrection.constrain(0., mesh.facesRight)
"""
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
    #velocity[..., mesh.exteriorFaces.value]=0.
    #velocity[0].constrain(U, mesh.facesRight | mesh.facesLeft)
    velocity[0, mesh.facesLeft.value] = U
    velocity[0, mesh.facesRight.value] = U
    #
    ##solve the pressure correction equation
    pressureCorrectionEq.cacheRHSvector()
    ## left bottom point must remain at pressure 0, so no correction
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    #
    ## update the pressure using the corrected value
    pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
    ## update the velocity using the corrected pressure
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)
    xVelocity[0]=U
#    xVelocity[nx-1]=U
    if sweep%10 == 0:
        viewer2.plot()
#        print 'sweep:',sweep,', x residual:',xres, ', p residual:',pres, ', continuity:',max(abs(rhs))

displacement = 125.
timeStep = 0.1* dx / U
elapsed = 0.

while elapsed < displacement/U:
    phi.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(var=phi, dt=timeStep)
    elapsed +=timeStep
    viewer.plot()
    viewer2.plot()
"""
raw_input("pause")
