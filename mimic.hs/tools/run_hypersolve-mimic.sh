rm -rf adapt$1
mkdir adapt$1
mkdir adapt$1/inputs
curr=$1
prev=$(($1 - 1))
echo $prev

python3 hypersolve_computehessian.py -m starship_v3.su2 -s restart_a$prev.snap -f temperature

mv out.hdf5 adapt$1/inputs/solution.h5 

cd adapt$1/inputs

ln -s ../../starship_v3.su2 mesh.su2

cd ../

cp -r ../run_hypersolve-mimic.sh
cp -r ../hypersolve_computehessian.py
cp -r ../metric.xml .
cp -r ../hs-starship.json .

sed -i "s/restart_a${prev}\.snap/restart_a${curr}.snap/g" hs.json

mpiexec -np 10 /nobackupp19/dekelsch/mimic-hls/utilities/readSU2/build/./readSU2

inf interpolate --source ../starship_v3.su2 --target mimic_mesh.su2 --snap ../restart_a$prev.snap

mpiexec -np 1000 hs run hs-starship.json
