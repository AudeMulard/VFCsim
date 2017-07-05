#!/bin/bash -x


mkdir -p sim_"${2}"
cd sim_"${2}"
mkdir -p phi pressure xVelocity yVelocity data
cd
python VFCsim/"${1}".py
mv phi*.png VFCsim/sim_"${2}"/phi
mv pressure*.png VFCsim/sim_"${2}"/pressure
mv XVelocity*.png VFCsim/sim_"${2}"/xVelocity
mv YVelocity*.png VFCsim/sim_"${2}"/yVelocity
mv essaidonne*.tsv VFCsim/sim_"${2}"/data
cd

