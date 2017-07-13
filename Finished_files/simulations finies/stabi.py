#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:21:51 2017

@author: aude
"""



from fipy import *
import random
import math

U = 0.8
Mobility = 0.2 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 1. #code starts going crazy below epsilon=0.1
l = 0.01 #this is lambda from Hamouda's paper
#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
dx = 0.25 #width of controle volume
nx = 150 #number of controle volume
dy = 1.
ny = 60
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
xVelocity = CellVariable(mesh=mesh, name='X Velocity')
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
    for i in range(30):
        a = random.gauss(0.2, 0.01)
        phi.setValue(1., where=(x > nx*dx * a ) & (y<2*(i+1)*dy) & (y>2*(i*dy)))


initialize(phi)


beta = CellVariable(mesh=mesh, name=r'$\beta$', value = beta2 * phi + beta1 * (1.-phi))
#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])
yVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[1])

ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue


#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)



#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------

#Phase

dexp = 1.
elapsed = 0.
duration = 150000.
while elapsed < duration:
    phi.updateOld()
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt)
    viewer.plot(filename="phase%d .png" % elapsed)
    if math.floor(elapsed)%10==0:
        TSVViewer(vars=(phi)).plot(filename="essaidonne150000%d.tsv" % elapsed)

raw_input("pause")
