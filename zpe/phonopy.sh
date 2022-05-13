#!/bin/sh

for dir_name in `ls -d */`
do
    echo $dir_name
    cd $dir_name
    if [ -d phonopy/ ]; then
        echo 'phonopy dir is already created'
    else
        logend=`tail -n 1 log | awk '{print $1}'`
        if [ $logend == 'reached' ]; then
            echo 'relaxation is finished'
            mkdir phonopy/
            cd phonopy
            cat >INCAR<<!
  PREC = Accurate
IBRION = -1
 ENCUT = 500
 EDIFF = 1.0e-08
ISMEAR = 1; SIGMA = 0.2
 IALGO = 38
 LREAL = .FALSE.
 LWAVE = .FALSE.
LCHARG = .FALSE.
  NPAR = 16
KSPACING = 0.314
!
            cp ../CONTCAR POSCAR
            phonopy -d --dim 2 2 2 --tolerance=2e-1
            for fname in `ls POSCAR-*`
            do
                disp="${fname/POSCAR/disp}"
                id="${fname/POSCAR-/}"
                mkdir $disp
                mv $fname $disp/POSCAR
                cd $disp
                cp ../INCAR .
                cp ../../POTCAR .
                cp ../../script.sh .
                sed -i "s/#SBATCH -J r/#SBATCH -J ${id}p/g" script.sh
                sed -i "s/#SBATCH -n 8/#SBATCH -n 16/g" script.sh
                sbatch script.sh
                cd .. 
            done
            cd ..
        else 
            echo 'relaxation is not yet finished, skipping this folder'
        fi
    fi
    cd ..
done

    