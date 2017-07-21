#!/bin/bash -x

date
mkdir -p Results_simulations/sim_"${2}"
cd Results_simulations/sim_"${2}"
mkdir -p phi pressure xVelocity yVelocity data
<<<<<<< HEAD
cd /home/aude/
#python VFCsim/Work_in_progress/"${1}".py "${2}"
=======
cd 
python VFCsim/Work_in_progress/"${1}".py "${2}"
>>>>>>> 650a562d65a9ec2a39b9019c9d33802bd220a1cf
mv phi*_"${2}".png VFCsim/Results_simulations/sim_"${2}"/phi
mv pressure*_"${2}".png VFCsim/Results_simulations/sim_"${2}"/pressure
mv XVelocity*_"${2}".png VFCsim/Results_simulations/sim_"${2}"/xVelocity
mv YVelocity*_"${2}".png VFCsim/Results_simulations/sim_"${2}"/yVelocity
mv essaidonne*_"${2}".tsv VFCsim/Results_simulations/sim_"${2}"/data
mv parameters_"${2}".csv VFCsim/Results_simulations/sim_"${2}"/data
for i in {1..9} 
do mv VFCsim/Results_simulations/sim_"${2}"/phi/phi"${i}"_"${2}".png VFCsim/Results_simulations/sim_"${2}"/phi/phi00"${i}"_"${2}".png
done
for i in {10..99} 
do mv VFCsim/Results_simulations/sim_"${2}"/phi/phi"${i}"_"${2}".png VFCsim/Results_simulations/sim_"${2}"/phi/phi0"${i}"_"${2}".png
done
cd VFCsim/Results_simulations/sim_"${2}"/phi/
mencoder mf:// -mf w=800:h=600:fps=04:type=png -ovc lavc -oac copy -o sim_"${2}".avi
cd
mv VFCsim/Results_simulations/sim_"${2}"/phi/sim_"${2}".avi VFCsim/Videos/sim_"${2}".avi
date
