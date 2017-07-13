#!/bin/bash -x


mkdir -p Results_simulations/sim_"${2}"
cd Results_simulations/sim_"${2}"
mkdir -p phi pressure xVelocity yVelocity data
cd
python VFCsim/Work_in_progress/"${1}".py
mv phi*.png VFCsim/sim_"${2}"/phi
mv pressure*.png VFCsim/Results_simulations/sim_"${2}"/pressure
mv XVelocity*.png VFCsim/Results_simulations/sim_"${2}"/xVelocity
mv YVelocity*.png VFCsim/Results_simulations/sim_"${2}"/yVelocity
mv essaidonne*.tsv VFCsim/Results_simulations/sim_"${2}"/data
cd

