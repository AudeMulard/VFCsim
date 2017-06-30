source activate MYFIPYENV
cd VFCsim/
mkdir sim50
cd sim50
mkdir phi pressure xVelocity yVelocity data
cd
cd /home/aude
python VFCsim/sim2D300617.py
cd VFCsim/sim50
mv phi*.png phi/phi*.png
mv pressure*.png pressure/pressure*.png
mv xVelocity*.png xVelocity/xVelocity*.png
mv yVelocity*.png yVelocity/yVelocity*.png
