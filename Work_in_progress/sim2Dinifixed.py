#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:21:51 2017

@author: aude

Correction du code avec le bon coefficient a
"""



from fipy import *
import random
import csv
import os, sys
import numpy



U = 0.8
Mobility = float(sys.argv[2]) #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 0.5 #code starts going crazy below epsilon=0.1
l = 0.1 #this is lambda from Hamouda's paper
duration = 0. #stabilisation phase
sweeps = 41 #stabilisation vitesse


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
dy = 0.5
ny = 120
mesh = Grid2D(dx=dx, nx=nx, dy=dy, ny=ny)
startpoint=0.1*nx*dx

parameters=csv.writer(open('parameters_'+sys.argv[1]+'.csv','w'), delimiter=' ', quotechar='|')
parameters.writerow(['U']+['Mobility']+['epsilon']+['l']+['dx']+['nx']+['dy']+['ny'])
parameters.writerow([U]+[Mobility]+[epsilon]+[l]+[dx]+[nx]+[dy]+[ny])

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
    a=[3.95549130019566, 3.7376223025404847, 3.781426716516598, 3.752031637116009, 3.8265329627205777, 3.7453348696609656, 3.876211474472181, 3.6540959635199766, 3.902913270765488, 3.828518166716552, 3.8624706215113895, 3.666325590059528, 3.7411330168887256, 3.7992747548764005, 3.805905761314782, 3.780207462350296, 3.6783079524287414, 3.709191395241718, 3.8417931071966303, 3.691299078907945, 3.866427397160816, 3.608713934134947, 3.7238804025000936, 3.7713944548794647, 3.6724960473085617, 3.774127062991821, 3.901775570080976, 3.8649996153271844, 3.7730868571629115, 3.8085580315703655, 3.7236826403727474, 3.848306192792678, 3.7049484646901645, 3.6711520946230536, 3.737741504063741, 3.7634200930098944, 3.628756535273059, 3.699370024784554, 3.8467052047875896, 3.7692561244591585, 3.8022586956288835, 3.6575031613575826, 3.5407583783741488, 3.8196313903547425, 3.5880087248604116, 3.9161125290293186, 3.611592382037223, 3.6873675476428085, 3.5635597940458137, 3.680940003766946, 3.8578012480278043, 3.8988745097464745, 3.8640648566874014, 3.713205096516372, 3.8780636033919134, 3.576765748011475, 3.822955691090964, 3.8197973976888786, 3.8752626032573576, 3.7186550811060073, 3.6796599166562114, 3.697542259955018, 3.824815629234484, 3.9024390756067384, 3.4673932926246795, 3.754207164140537, 3.8093312821166156, 3.819086304690003, 3.6373397569591353, 3.4741643255261376, 3.7695838334347522, 3.6988994855043202, 3.7760772039588195, 3.802456331557675, 3.80949970857435, 3.9686487070316576, 3.7607049343967915, 3.7278587641925656, 3.7977219276710175, 3.730603212760668, 3.7146073633852574, 3.5661892437649314, 3.7917343004980544, 3.7579008904701636, 3.62609893180622, 3.793818852355885, 3.770115782474133, 3.9030838642229546, 3.871271665574576, 3.629421958288578, 3.722892728493727, 3.683700329464345, 3.8525545366375633, 3.714062117102295, 3.776223741820282, 3.7265868218494576, 3.8889982839293533, 3.637744686666094, 3.862814239355127, 3.5499089676820983, 3.5603759616683517, 3.7244587578743045, 3.6576768959871093, 3.6654514561358886, 3.823250675252249, 3.83218859968082, 3.650678463714495, 3.8268178650111953, 3.755033146898565, 3.6835788338952824, 3.72345895417868, 3.70778462632419, 3.701452274041286, 3.5967228354791794, 3.713728644415371, 3.7429404200834386, 3.81077459333936, 3.709396993400353, 3.689159536881579, 3.687686367269336]
    for i in range(ny):
        phi.setValue(0.5*(1+numerix.tanh((x-a[i])/(2*epsilon))), where=(y<(i+1)*dy) & (y>i*dy))


initialize(phi)


beta = CellVariable(mesh=mesh, name=r'$\beta$', value = beta2 * phi + beta1 * (1.-phi))
#-----------------------------------------------------------------------
#-------------------------Velocity and pressure-------------------------
#-----------------------------------------------------------------------


xVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[0])
yVelocityEq = (ImplicitSourceTerm(coeff=beta) + pressure.grad[1])

coeff = 1./ beta.arithmeticFaceValue
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

#Remove oscillations
from fipy.variables.faceGradVariable import _FaceGradVariable


#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi), datamin=0., datamax=1.)
viewer2 = Viewer(vars = (xVelocity), datamin=0.5, datamax=1.)
viewer3 = Viewer(vars = (yVelocity), datamin=0., datamax=0.2)
viewer4 = Viewer(vars = (pressure), datamin=0., datamax=nx*dx*U)



#-----------------------------------------------------------------------
#---------------------------Initialization------------------------------
#-----------------------------------------------------------------------

#Phase

dexp = 1.
elapsed = 0.

while elapsed < duration:
    phi.updateOld()
    dt = min(30, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(var=phi, dt = dt, solver=LinearGMRESSolver())
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
    velocity[0, mesh.facesRight.value] = U
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
    xVelocity[nx-1]=U




displacement = 35.
timeStep = float(sys.argv[3]) * dx / U #less than one space step per time step
elapsed = 0.
 
while elapsed < displacement/U:
    phi.updateOld()
    res = 1e+10
    while res > 1e-6:
        res = eq.sweep(var=phi, dt=timeStep, solver=LinearGMRESSolver())
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
        velocity[0, mesh.facesRight.value] = U
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
        xVelocity[nx-1]=U
    elapsed +=timeStep
    viewer.plot(filename='phi%d_' % elapsed +sys.argv[1]+'.png')
    viewer2.plot(filename='XVelocity%d_' % elapsed +sys.argv[1]+'.png')
    viewer4.plot(filename='pressure%d_' % elapsed +sys.argv[1]+'.png')
    viewer3.plot(filename='YVelocity%d_' % elapsed +sys.argv[1]+'.png')
    TSVViewer(vars=(phi, xVelocity, yVelocity, pressure,beta)).plot(filename='essaidonne%d_'% elapsed +sys.argv[1]+'.tsv')
    print(elapsed)



