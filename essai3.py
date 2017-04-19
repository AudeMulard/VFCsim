# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
dx = 0.02 * W #width of controle volume
nx = L / dx #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity2 = 1.
Mobility = 0.75 #ratio of the two viscosities
viscosity1 = viscosity2 * Mobility
permeability1 = permeability2 = 1.
beta1 = - viscosity1 / permeability1
beta2 = - viscosity2 / permeability2


#Variable of the fluids
pressure = CellVariable(mesh=mesh, name='pressure')
xVelocity = CellVariable(mesh=mesh, name='X Velocity')
velocity = FaceVariable(mesh=mesh, rank=1)

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld = 1)

#New values
beta = Variable(name='beta')
beta.setValue = beta1 * phi + beta2 * (1-phi)

#Parameters
Cahn_number = 0.001
epsilon = Cahn_number * W
M = Mobility * epsilon**2
l = 1.

#Cahn-Hilliard equationimport numpy
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (3 * PHI**2 - 3 * PHI + 1/2)
eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) + DiffusionTerm(coeff=(M, l)))
#eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) + DiffusionTerm(coeff=(M, l)))
#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(1.)
	phi.setValue(0., where=x > L/2)

initialize(phi)

print(phi)
for i in range(10):
    res = eq.sweep(var=phi, dt=1e-6)
    print(phi)

viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)
"""
#Boundary conditions
Q = 1. #rate of injection
Uinf = Q / (b*W)
xVelocity.constrain(Uinf, where=mesh.exteriorFaces)
pressure.constrain(0., where=mesh.facesLeft)

#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
#viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)
#viewer2 = Viewer(vars = (phi, pressure, xVelocity))


#--------------------------------import numpy---------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------

#Equations
pressureCorrection = CellVariable(mesh=mesh)
pressureCorrection.constrain(0., mesh.facesLeft)

xVelocityEq = ImplicitSourceTerm(beta) - pressure.grad[0] #Darcy's law

ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence
#velocity and pressure field:
pressureRelaxation = 0.8
velocityRelaxation = 0.5

#Rhie-Chow correction
from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------

sweeps=300

xVelocityEq.cacheMatrix()
xVelocityEq.solve(var=xVelocity)
xmat = xVelocityEq.matrix
ap[:] = -xmat.takeDiagonal()
##Rhie-chow correction
presgrad = pressure.grad
facepresgrad = _FaceGradVariable(pressure)
velocity[0] = xVelocity.arithmeticFaceValue \
     + contrvolume / ap.arithmeticFaceValue * \
     (presgrad[0].arithmeticFaceValue-facepresgrad[0])

pressureCorrectionEq.cacheRHSvector()
## left bottom point must remain at pressure 0, so no correction
pressureCorrectionEq.solve(var=pressureCorrection)
rhs = pressureCorrectionEq.RHSvector
print(rhs)

## update the pressure using the corrected value
pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
## update the velocity using the corrected pressure
xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)

#viewer2.plot()
if __name__ == '__main__':
    print max(abs(rhs))


dexp = -5
print(phi)
timeStep = 1e-6
phi.updateOld()
res = 1e+10
eq.solve(phi, dt=timeStep)
elapsed = 0.
duration = 1000.
for i in range(10):
    elapsed += timeStep
    res = eq.sweep(var=phi, dt=1e-6)
    viewer.plot()
    print(phi)
"""
