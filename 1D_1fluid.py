#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:27:06 2017

@author: am2548
"""

from fipy import *

#Parameters
#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.25 #width of controle volume
nx = 10 #number of controle volume
#L = 1.0
#N = 50
#dL = L / N
viscosity = 1.
permeability = 1.
U = 1.
beta = permeability/viscosity
pressureRelaxation = 0.8 #some numerical parameters
velocityRelaxation = 0.5


#mesh
#mesh = Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)
mesh = Grid1D(nx=nx, dx=dx)

#Variables
pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X Velocity')

velocity = FaceVariable(mesh=mesh, rank=1)

xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])


ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
X = mesh.cellCenters[0]
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

#no-slip boundary conditions
xVelocity.constrain(U, mesh.facesRight | mesh.facesLeft)
X = mesh.faceCenters
pressureCorrection.constrain(0., mesh.facesLeft)

#Viewers
if __name__ == '__main__':
    viewer = Viewer(vars=(pressure, xVelocity))

sweeps = 80
for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred values
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
    if __name__ == '__main__':
        if sweep%10 == 0:
            print(xVelocity)
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
    velocity[0].constrain(U, mesh.facesRight | mesh.facesLeft)
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
    if __name__ == '__main__':
            if sweep%10 == 0:
                print 'sweep:',sweep,', x residual:',xres, ', p residual:', pres, ', continuity:', max(abs(rhs))
                viewer.plot()            
