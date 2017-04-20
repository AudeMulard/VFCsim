#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:51:00 2017

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
dx = 0.02 * W #width of controle volume
nx = L / dx #number of controle volume
mesh = Grid1D(dx=dx, nx=nx)

#-----------------------------------------------------------------------
#---------------------Description of the fluids-------------------------
#-----------------------------------------------------------------------

#Parameters of the fluids
viscosity = 1.
permeability = 1.
beta = - viscosity / permeability


#Variable of the fluids
pressure = CellVariable(mesh=mesh, name='pressure')
xVelocity = CellVariable(mesh=mesh, name='X Velocity')
velocity = FaceVariable(mesh=mesh, rank=1)

#-----------------------------------------------------------------------
#-------------------------Boundary Conditions---------------------------
#-----------------------------------------------------------------------

#Boundary conditions
Q = 1. #rate of injection
Uinf = Q / (b*W)
xVelocity.constrain(Uinf, where=mesh.facesLeft | mesh.facesRight)
pressure.constrain(0., where=mesh.facesLeft)

#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (pressure, xVelocity, velocity), xmin=0., xmax=1., ymin=0., ymax=1., colorbar=True)
#viewer2 = Viewer(vars = (phi, pressure, xVelocity))


#-----------------------------------------------------------------------
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

sweeps = 300

xVelocityEq.cacheMatrix()
xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
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
pres = pressureCorrectionEq.sweep(var=pressureCorrection)
rhs = pressureCorrectionEq.RHSvector

## update the pressure using the corrected value
pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
## update the velocity using the corrected pressure
xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)

if __name__ == '__main__':
    print('x residual:' , xres , ', p residual:', pres, ', continuity:', max(abs(rhs)))
    viewer.plot()
