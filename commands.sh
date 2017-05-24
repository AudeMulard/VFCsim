date > chrono1
cat chrono1 >> chrono2
python cahnhillarparex.py
date > chrono1
cat chrono1 >> chrono2
python cahnhillarparex.py --pysparse
date > chrono1
cat chrono1 >> chrono2
mpirun -np 1 python cahnhillarparex.py --trilinos
date > chrono1
cat chrono1 >> chrono2
mpirun -np 2 python cahnhillarparex.py --trilinos
date > chrono1
cat chrono1 >> chrono2
