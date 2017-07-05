#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 09:06:19 2017

@author: aude
"""

from fipy import *

L = 1.0
N = 50
dL = L/N
viscosity = 1
U = 1.
pressureRelaxation = 0.8
velocityRelaxation=0.5
sweeps=300

mesh=Grid2D(nx=N, ny=N, dx=dL, dy=dL)

pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X velocity')
yVelocity = CellVariable(mesh=mesh, name='Y velocity')

velocity = FaceVariable(mesh=mesh, rank=1)

xVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.grad.dot([1.,0.])
yVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.grad.dot([0., 1.])

ap = CellVariable(mesh=mesh, value=1.)
coeff= 1./ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

from fipy.variables.faceGradVariable import _FaceGradVariable
volume=CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

xVelocity.constrain(0., mesh.facesRight | mesh.facesLeft | mesh.facesBottom)
xVelocity.constrain(U, mesh.facesTop)
yVelocity.constrain(0., mesh.exteriorFaces)
X, Y = mesh.faceCenters
pressureCorrection.constrain(0., mesh.facesLeft & (Y < dL))

viewer = Viewer(vars=(pressure, xVelocity, yVelocity, velocity), xmin=0., xmax=1., ymin=0., ymax=1., colorbar=True)

for sweep in range(sweeps):
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix
    yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velocityRelaxation)
    ap[:] = -xmat.takeDiagonal()
    presgrad = pressure.grad
    facepresgrad = _FaceGradVariable(pressure)
    velocity[0] = xVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue - facepresgrad[0])
    velocity[1] = yVelocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[1].arithmeticFaceValue - facepresgrad[1])
#    velocity[..., mesh.exteriorFaces.value]=0.
#    velocity[0, mesh.facesTop.value]=U
    pressureCorrectionEq.cacheRHSvector()
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    pressure.setValue(pressure+pressureRelaxation*pressureCorrection)
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / ap * mesh.cellVolumes)
    yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / ap * mesh.cellVolumes)
    if __name__ == '__main__':
        if sweep%10 ==0:
            print('sweep:',sweep,', x residual:',xres, ', y residual',yres, ', p residual:',pres, ', continuity:',max(abs(rhs)))
            viewer.plot()
