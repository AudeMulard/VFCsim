#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 08:56:05 2017

@author: aude
"""


from fipy import *
import random

#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.05 #width of controle volume
dy = 0.10
nx = 200
ny = 50 #number of controle volume
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

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
yVelocity = CellVariable(mesh=mesh, name='Y Velocity')
velocity = FaceVariable(mesh=mesh, rank=1)

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld=1)
#New values
beta = CellVariable(mesh=mesh, name='beta', value = beta2 * phi + beta1 * (1-phi), hasOld=1)

#Parameters
#Cahn_number = 0.001
epsilon = 0.2
M = Mobility * epsilon**2
l = 1.

phi.faceGrad.constrain([0.], mesh.facesRight)
#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (6.* PHI*(PHI-1.) + 1.)
eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) - DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------

xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad.dot([1.,0.]))
yVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad.dot([0.,1.]))

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
x,y = mesh.cellCenters
def initialize(phi):
    phi.setValue(0.)
    for i in range(ny):
        a = random.gauss(0.5, 0.1)
        phi.setValue(1., where=(x > nx*dx * a ) & (y<(i+1)*dy) & (y>(i*dy)))
    
    
#    a = random.gauss(0.5, 0.01)
#    phi.setValue(a, where=(x > nx*dx/2-5*epsilon) & (x < nx*dx/2+5*epsilon))




initialize(phi)


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
viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)
#viewer2 = Viewer(vars = (xVelocity, yVelocity), datamin=-1., datamax=3.)

#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------
raw_input("pause")


#Phase

timeStep = 10.
for i in range(1):
    phi.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(var=phi, dt=timeStep)
    if __name__ == '__main__':
        viewer.plot()       
    raw_input("pause")


"""    
def updatephi():
    phi.updateOld()
    beta.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(var=phi, dt=timeStep)
    if __name__ == '__main__':
        viewer.plot()

for i in range(50):
    updatephi()


if __name__ == '__main__':
    viewer.plot()

raw_input("pause")
"""
#Pressure and velocity
pressureRelaxation = 0.8
velocityRelaxation = 0.5
X, Y = mesh.faceCenters
xVelocity.constrain(U, mesh.facesLeft)
#xVelocity.constrain(U, mesh.facesRight)
pressureCorrection.constrain(0., mesh.facesRight)
sweeps = 10



for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred values
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix
    yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation)
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
    velocity[1] = yVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[1].arithmeticFaceValue-facepresgrad[1])
    #velocity[..., mesh.exteriorFaces.value]=0.
    #velocity[0].constrain(U, mesh.facesRight | mesh.facesLeft)
    #velocity[0, mesh.facesLeft.value] = U
    #velocity[0, mesh.facesRight.value] = U
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
    yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / ap * mesh.cellVolumes)
    xVelocity[0]=U


viewer2 = Viewer(vars = (xVelocity, yVelocity), datamin=-1., datamax=3.)
viewer2.plot()

raw_input("pause")

x = mesh.cellCenters[0]    

displacement = 3.
#velocity1 = 1.
timeStep = .1 * dx / U
elapsed = 0.
while elapsed < displacement/U:
    phi.updateOld()
    beta.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(var=phi, dt=timeStep)
    if __name__ == '__main__':
        viewer.plot()       
    if elapsed==0%140:
        raw_input("pause")
    elapsed +=timeStep
    for sweep in range(sweeps):
        ##Solve the Stokes equations to get starred values
        xVelocityEq.cacheMatrix()
        xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
        xmat = xVelocityEq.matrix
        yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation)
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
        velocity[1] = yVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[1].arithmeticFaceValue-facepresgrad[1])
        #velocity[..., mesh.exteriorFaces.value]=0.
        #velocity[0].constrain(U, mesh.facesRight | mesh.facesLeft)
        #velocity[0, mesh.facesLeft.value] = U
        #velocity[0, mesh.facesRight.value] = U
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
        yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / ap * mesh.cellVolumes)
        xVelocity[0]=U

#    if sweep%10 == 0:
#        viewer2.plot()



viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)
viewer.plot()    
viewer2 = Viewer(vars = (xVelocity, yVelocity), datamin=-1., datamax=3.)
viewer2.plot()

raw_input("pause")

