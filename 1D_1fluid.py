#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:27:06 2017

@author: am2548
"""

from fipy import *

#Parameters
L = 1.0
N = 50
dL = L / N
viscosity = 1
permeability = 1
U = 1.
beta = permeability/viscosity
pressureRelaxation = 0.8
velocityRelaxation = 0.5

if __name__ == "__main__":
    sweeps = 300
else:
    sweeps = 5

#mesh
mesh = Grid1D(nx=N, dx=dL)

#Variables
pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X Velocity')

velocity = FaceVariable(mesh=mesh, rank=1)

xVelocityEq = ImplicitSourceTerm(beta) + pressure.grad[0.]


ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
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
    viewer = Viewer(vars=(pressure, xVelocity, velocity), xmin=0., xmax=1., ymin=0., colorbar=True)

#iterations
for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred values
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity,
                             underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix
    ##update the ap coefficient from the matrix diagonal
    ap[:] = -xmat.takeDiagonal()
    #
    ##update the face velocities based on starred values with the Rhi-Chow correction
    #cell pressure gradient
    presgrad = pressure.grad
    #face pressure gradient
    facepresgrad = _FaceGradVariable(pressure)
    #
    velocity[0] = xVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
    velocity[0, mesh.facesLeft.value] = U
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
    #
    if __name__ == '__main__':
        if sweep%10 == 0:
            print 'sweep:',sweep,', x residual:',xres, ', p residual:', pres, ', continuity:', max(abs(rhs)), xVelocity
            viewer.plot()

#test values in the last cell
print numerix.allclose(pressure.globalValue[...,-1], 162.790867927)
print numerix.allclose(xVelocity.globalValue[...,-1], 0.265072740929)
