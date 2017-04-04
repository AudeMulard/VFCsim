from fipy import *

#------------------------------------------------------
#--------------------Initialization--------------------
#------------------------------------------------------

#Mesh
nx = 1000
dx=0.25
mesh = Grid2D(nx=nx, ny=ny, dx=dx, dy=dy)

#Variables
phi = CellVariable(name=r"$\phi$", mesh=mesh)
pressure = CellVariable(name='pressure', mesh=mesh)
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X Velocity')
yVelocity = CellVariable(mesh=mesh, name='Y Velocity')
Velocity = CellVariable(mesh=mesh)

#Parameters
beta = -1
pressureRelaxation = 0.8
velocityRelaxation = 0.5
sweeps = 300
M = 1.5
lambda = 1.
dexp = -5
elapsed = 0
duration=1000

#Boundary conditions
velocity.constrain(0., mesh.facesRight | mesh.facesLeft)
X = mesh.faceCenters

#InitialValues
phi.setValue(GaussianNoiseVariable(mesh=mesh, mean=0.5,variance=0.01)

#Viewers
if __name__ == '__main__':
	viewer1 = Viewer(vars=(pressure, velocity, Velocity), xmin=0., xmax=nx*dx)
	viewer2 = Viewer(vars=(phi,), datamin=0., datamax=1.)
	
#------------------------------------------------------
#--------------Pressure-Velocity coupling--------------
#------------------------------------------------------



#Darcy's law
velocityEq = ( ImplicitSourceTerm(coeff=beta) == pressure.grad.dot([1.,0.]) )

#Correction with continuity
ap = CellVariable(mesh=mesh, value=1.)
coeff = 1./ ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume=volume.arithmeticFaceValue

for sweep in range(sweeps):
	## solve the Stokes equations to get starred values
	velocityEq.cacheMatrix()
	xres = velocityEq.sweep(var=velocity, underRelaxation=velocityRelaxation)
	xmat = velocityEq.matrix
	##update the ap coefficient from the matrix diagonal
	ap[:] = -xmat.takeDiagonal()
	##update face velocities based on starred values with Rhie-Chow correction
	# cell pressure gradient
	presgrad = pressure.grad
	#face pressure gradient
	facepresgrad = _FaceGradVariable(pressure)
	#
	Velocity = velocity.arithmeticFaceValue + contrvolume / ap.arithmeticFaceValue * (presgrad[0].arithmeticFaceValue-facepresgrad[0])
	##solve pressure correction equation
	pressureCorrectionEq.cacheRHSvector()
	pres = pressureCorrectionEq.sweep(var=pressureCorrection)
	rhs = pressureCorrectionEq.RHSvector
	##update the pressure using the corrected value
	pressure.setValue(pressure + pressureRelaxation * pressureCorrection)
	##update the velocity using the corrected pressure
	velocity.setValue(velocity - pressureCorrection.grad[0]/ap * mesh.cellVolumes)
	

#------------------------------------------------------
#----------------Cahn-Hilliard evolution---------------
#------------------------------------------------------

PHI = phi.arithmeticFaceValue
coeff = M * lambda * (3 * PHI**2 - 3 * PHI +1/2)
eq = (TransientTerm() + ConvectionTerm(Velocity) == DiffusionTerm(coeff=coeff) + DiffusionTerm(coeff=(M, lambda)))

while elapsed < duration:
	dt = min(100, numerix.exp(dexp))
	elapsed += dt
	dexp += 0.01
	eq.solve(phi, dt=dt)
	if __name__ == '__main__':
		viewer.plot()
	elif (max(phi.globalValue) > 0.7) and (min(phi.globalValue) < 0.3) and elapsed > 10.:
		break


