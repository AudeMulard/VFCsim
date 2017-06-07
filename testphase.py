from fipy import *

nx = 400
dx = 5e-6
L = nx * dx
mesh = Grid1D(dx=dx, nx=nx)
phase = CellVariable(name="phase", mesh=mesh, hasOld=1)
x = mesh.cellCenters[0]
phase.setValue(1.)
phase.setValue(0., where=x > L/2)

Lv = 2350
Tm = 1728.
T = Variable(value=Tm)
enthalpy = Lv * (T-Tm) / Tm

delta = 1.5 *dx
sigma = 3.7e-5
beta = 0.33
kappa = 6 * sigma * delta
W = 6 * sigma/delta
Mphi = Tm * beta / (6. * Lv * delta)

if __name__ == '__main__':
    displacement = L * 0.5
else:
    displacement = L * 0.025

viewer = Viewer(vars=(phase,))
viewer.plot()


mPhi = -((1 - 2 * phase) * W + 30 * phase * (1 - phase) * enthalpy)
dmPhidPhi = 2 * W - 30 * (1 - 2 * phase) * enthalpy
S1 = dmPhidPhi * phase * (1 - phase) + mPhi * (1 - 2 * phase)
S0 = mPhi * phase * (1 - phase) - S1 * phase

eq = TransientTerm(coeff=1/Mphi) == DiffusionTerm(coeff=kappa) + S0 + ImplicitSourceTerm(coeff = S1)

timeStep = 1e-6
for i in range(10):
    phase.updateOld()
    res = 1e+10
    while res > 1e-5:
        res = eq.sweep(var=phase, dt=timeStep)

viewer.plot()

T.setValue(T() - 1)

velocity = beta * abs(Tm - T())
timeStep = .1 * dx / velocity
elapsed = 0

while elapsed < displacement / velocity:
    phase.updateOld()
    res = 1e+10
    while res > 1e-5:
        res = eq.sweep(var=phase, dt=timeStep)
    elapsed += timeStep
    viewer.plot()

raw_input("pause")
