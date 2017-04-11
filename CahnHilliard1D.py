from fipy import *


#-----------------------------------------------------------------------
#------------------------Geometry and mesh------------------------------
#-----------------------------------------------------------------------

#Space
L = 1. #length
W = 1. #width: characteristic length
b = 1. #gap

#Mesh
nx = 400 #number of controle volume
dx = L/nx #width of controle volume
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
velocity = CellVariable(mesh=mesh, name='velocity')

#-----------------------------------------------------------------------
#------------------------Phase-field model------------------------------
#-----------------------------------------------------------------------

#Order Parameter
phi = CellVariable(name=r'$\phi$', mesh=mesh, hasOld = 1)

#Parameters
Cahn_number = 0.001
epsilon = Cahn_number * W
M = Mobility * epsilon**2
l = 1.

#New values
beta = CellVariable(mesh=mesh, name='beta')
beta.setValue = beta1 * phi + beta2 * (1-phi)

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (3 * PHI**2 - 3 * PHI + 1/2)
eq = (TransientTerm() + ConvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) + DiffusionTerm(coeff=(M, l)))

#-----------------------------------------------------------------------
#---------------Initialization and Boundary Conditions------------------
#-----------------------------------------------------------------------

x = mesh.cellCenters[0]
def initialize(phi):
	phi.setValue(1.)
	phi.setValue(0., where=x > L/2)

initialize(phi)

#Boundary conditions
Q = 1. #rate of injection
Uinf = Q / (b*W)
velocity.constrain(Uinf, where=mesh.facesRight | mesh.facesLeft)

#-----------------------------------------------------------------------
#-------------------------------Viewers---------------------------------
#-----------------------------------------------------------------------

#Viewer
viewer = Viewer(vars = (phi,), datamin=0., datamax=1.)
viewer2 = Viewer(vars = (phi, pressure, velocity))


#-----------------------------------------------------------------------
#----------------Computation of velocity and pressure-------------------
#-----------------------------------------------------------------------

sweeps = 300
pressureRelaxation = 0.8
velocityRelaxation = 0.5
pressureCorrection = CellVariable(mesh=mesh)
Velocity = FaceVariable(mesh=mesh)
from fipy.variables.faceGradVariable import _FaceGradVariable
velocityEq = (ImplicitSourceTerm(coeff=beta) == pressure.grad[0.,])
ap = CellVariable(mesh=mesh, value=1.)

from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - Velocity.divergence
X = mesh.faceCenters[0]
pressureCorrection.constrain(0., mesh.facesLeft)

#-----------------------------------------------------------------------
#------------Phase field formation: Initial equilibrium-----------------

timeStep = 1e-6
for i in range(10):
	phi.updateOld()
	res = 1e+10
	while res > 1e-5:
		eq.sweep(var=phi, dt=timeStep)

#-----------------------------------------------------------------------
#----------------------------Dynamics-----------------------------------
#-----------------------------------------------------------------------


#Resolution de l'equation de phase field
##voir examples.phase.simple

timeStep = .1 * dx / velocity
elapsed = 0
while elapsed < displacement/velocity:
	phi.updateOld()
	res = 1e+10
	while res > 1e-5:
		eq.sweep(var=phi, dt=timeStep)
	elapsed += timeStep

#Resolution des eq de mouvements couplees sur staggered grid avec SIMPLE

#Variables
pressure = CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X Velocity')

xVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.grad.dot([1.,])


ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue * mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - xVelocity.grad[1.,]



#no-slip boundary conditions
xVelocity.constrain(0., mesh.facesRight | mesh.facesLeft | mesh.facesBottom)
xVelocity.constrain(U, mesh.facesTop)
yVelocity.constrain(0., mesh.exteriorFaces)
X = mesh.faceCenters[0]
pressureCorrection.constrain(0., mesh.facesLeft)

#Viewers
if __name__ == '__main__':
    viewer = Viewer(vars=(pressure, xVelocity, yVelocity), xmin=0., xmax=1., ymin=0., colorbar=True)

#iterations
for sweep in range(sweeps):
    ##Solve the Stokes equations to get starred values
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity,
                             underRelaxation=velocityRelaxation)
    xmat = xVelocityEq.matrix
    yres = yVelocityEq.sweep(var=yVelocity,
                             underRelaxation=velocityRelaxation)
    ##update the ap coefficient from the matrix diagonal
    ap[:] = -xmat.takeDiagonal()
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
            print 'sweep:',sweep,', x residual:',xres, ', y residual:',yres, ', p residual:', pres, ', continuity:', max(abs(rhs))
            viewer.plot()


#Algorithme global
for i in range(300):
	solve phi
	solve u,v
	viewer.plot()

