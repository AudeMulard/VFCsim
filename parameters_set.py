import math

U = 0.8
Mobility = 0.1 #ratio of the two viscosities; M_c in Hamouda's paper
epsilon = 0.25 #code starts going crazy below epsilon=0.1
l = 0.1 #this is lambda from Hamouda's paper
duration = 1. #stabilisation phase
sweeps = 41 #stabilisation vitesse
startpoint=0.1
w=10

#Mesh
dx = 0.15 #width of controle volume
nx = 500 #number of controle volume
dy = 1.
ny = 60


#print(1-((2*math.pi/w)**2*2*math.sqrt(2)*l/(U*epsilon)))
print(2*math.pi*math.sqrt(2*math.sqrt(2)*l/((1-Mobility)*U*epsilon)))
