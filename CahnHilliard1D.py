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
Mobility = 1.5 #ratio of the two viscosities
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
phi = CellVariable(name=r'$\phi$', mesh=mesh)

#Parameters
Cahn_number = 0.001
epsilon = Cahn_number * W
M = Mobility * epsilon**2
l = 1.

#New values
beta = beta1 * phi + beta2 * (1-phi)

#Cahn-Hilliard equation
PHI = phi.arithmeticFaceValue #result more accurate by non-linear interpolation
coeff1 = Mobility * l * (3 * PHI**2 - 3 * PHI + 1/2)
eq = (TransientTerm() + AdvectionTerm(velocity) == DiffusionTerm(coeff=coeff1) + DiffusionTerm(coeff=(M, l)))

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
velocityEq = (ImplicitSourceTerm(coeff=beta) == pressure.grad)
ap = CellVariable(mesh=mesh, value=1.)

from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

coeff = 1./ ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence
X = mesh.faceCenters
pressureCorrection.constrain(0., mesh.facesLeft)

def solveStokes(phi):
	for sweep in range(sweeps):
		velocityEq.cacheMatrix()
		xres = velocityEq.sweep(var=velocity, underRelaxation=velocityRelaxation)
		xmat = velocityEq.matrix
		ap[:] = - xmat.takeDiagonal()
		presgrad = pressure.grad
		facepresgrad = _FaceGradVariable(pressure)
		Velocity=velocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad.arithmeticFaceValue-facepresgrad)
		pressureCorrectionEq.cacheRHSvector()
		pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
		velocity.setValue(velocity - pressureCorrection.grad / ap * mesh.cellVolumes)
		viewer2.plot()
		

#-----------------------------------------------------------------------
#--------------------------Time Resolution------------------------------
#-----------------------------------------------------------------------

dt = numerix.exp(-5)
time = 0.

while time < 1000:
	time +=dt
	eq.solve(phi, dt=dt)
	solveStokes(phi)
	viewer.plot()


