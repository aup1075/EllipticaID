module unload gcc ; module unload aocc ; module unload openmpi ;module load cpu/0.15.4 ; module load gcc/9.2.0 ;module load openmpi/3.1.6 ;module load hdf5/1.10.6 ;module load boost ;module load openblas ; module load fftw ; module load sdsc ; module load gsl/2.5

cd arrangements
ln -s ../repos/EllipticaID EllipticaID
cd ExternalLibraries
ln -s ../../repos/EllipticaID/ExternalLibraries-SuiteSparse SuiteSparse
ln -s ../../repos/EllipticaID/ExternalLibraries-Elliptica Elliptica
cd ../../
cp repos/EllipticaID/utils/expanse.ini simfactory/mdb/machines/
cp repos/EllipticaID/utils/expanse.run simfactory/mdb/runscripts/
cp repos/EllipticaID/utils/expanse.sub simfactory/mdb/submitscripts/
cp repos/EllipticaID/utils/expanse.cfg simfactory/mdb/optionlists/
