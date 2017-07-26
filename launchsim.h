#!/bin/bash -x

date
mkdir -p Results_simulations/sim_"${2}"
cd Results_simulations/sim_"${2}"
mkdir -p phi pressure xVelocity yVelocity data

<<<<<<< HEAD
cd
python VFCsim/Work_in_progress/"${1}".py "${2}"
=======
cd /home/aude/
python VFCsim/Work_in_progress/"${1}".py "${2}" "${3}"
>>>>>>> e51b9a71474474c7a3a823196cf983194b291cff

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
cd /home/aude/
mv VFCsim/Results_simulations/sim_"${2}"/phi/sim_"${2}".avi VFCsim/Videos/sim_"${2}".avi
date
cd VFCsim/
