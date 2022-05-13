import os
import re


def get_system(number):
    system_1 = {
        "partition": "lenovo",
        "modules": "module load vasp/6.1.1",
        "N": 1,
        "n": 8,
        "mpirun": "$(which mpirun)",
        "npar": 4,
    }

    system_2 = {
        "partition": "cpu",
        "modules": "module load ",
        "N": 1,
        "n": 8,
        "mpirun": "$(which mpirun)",
        "npar": 4,
    }

    if number == 1:
        return system_1
    if number == 2:
        return system_2


def compose_potcar(structure, specific_path, dir_name):
    formula = structure.formula
    split = re.split(r'(\d+)', formula)
    potcar_paths = list()
    for specie in split:
        specie = specie.replace(' ', '')
        if specie.isalpha():
            potcar = os.path.join(specific_path, f'POTCAR_{specie}')
            assert os.path.isfile(potcar)
            potcar_paths.append(potcar)
    with open(os.path.join(dir_name, 'POTCAR'), 'w') as outfile:
        for potcar in potcar_paths:
            with open(potcar, 'r') as infile:
                for line in infile:
                    outfile.write(line)


def get_slurm_script(properties, system, path):
    script = f"""#!/bin/sh
#SBATCH -o out -e err
#SBATCH -p {system['partition']}
#SBATCH -J r{properties['id']}
#SBATCH -N {system['N']}
#SBATCH -n {system['n']}
#SBATCH -t 06:00:00

{system['modules']}

{system['mpirun']} $(which vasp_std) > log

""".format()

    with open(os.path.join(path, 'script.sh'), 'w') as f:
        f.write(script)


def write_incar(option, properites, pressure, path):
    relaxation_incar = f"""SYSTEM = EA{properites['id']} {properites['stability']} {properites['formula']} {properites['space_group']}
PREC = Normal
ENCUT = 600
EDIFF = 1e-8
IBRION = 2
ISIF = 3
NSW = 666
ISMEAR = 1 ; SIGMA = 0.2
POTIM = 0.250
ISTART = 0
LCHARG = FALSE
LWAVE = FALSE
KSPACING = 0.314
PSTRESS = {pressure}
NPAR = 4
""".format()

    if option == 'relaxation':
        incar = relaxation_incar
    elif option == 'phonons:':
        pass

    with open(os.path.join(path, 'INCAR'), 'w') as f:
        f.write(incar)


def launch_all(zpe_path):
    script = """#!/bin/sh

for dir_name in `ls -d */`
do
    echo $dir_name
    cd $dir_name
    sbatch script.sh
    cd ..
done
"""

    with open(os.path.join(zpe_path, 'launch_all.sh'), 'w') as f:
        f.write(script)


def phonopy(zpe_path):
    script = r"""#!/bin/sh

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

    """

    with open(os.path.join(zpe_path, 'phonopy.sh'), 'w') as f:
        f.write(script)


def gather(zpe_path):
    script = r"""#!/bin/sh

for dir_name in `ls -d */`
do
    cd $dir_name
    stat=$dir_name": "
    if [ -d phonopy/ ]; then 
        cd phonopy/
        phonstat=""
        if [ ! -f total_dos.pdf ]; then
            for disp in `ls -d disp-*/`
            do
                cd $disp
                if [ -f log ]; then
                    err=`cat err`
                    log=`tail -n 1 log | awk '{print $2}'`
                    if [ "$log" == "F=" ]; then
                        problems=`grep 'very serious problems' log`
                        if [ "$problems" == "" ]; then
                            phonstat=$phonstat""
                            else
                            phonstat=$phonstat"problems in $disp "
                            #sed -i "s/NPAR = 4/NPAR = 8/g" INCAR
                            #sbatch script.sh
                        fi
                    else
                        if [ "$err" != "" ]; then
                            phonstat=$phonstat" error in $disp"
                        else
                            phonstat=$phonstat"$disp still working "
                        fi
                    fi
                else
                    phonstat=$phonstat"$disp not yet started "
                fi
                cd ..
            done
        fi
        cd ..
    else
        stat=$stat"no phonopy dir created yet"
    fi
    if [ "$phonstat" == "" ]; then
        phonstat="all done"
        cd phonopy
        if [ ! -f FORCE_SETS ]; then
            disp_s=`ls -d disp-*/ | head -n 1`
            disp_f=`ls -d disp-*/ | tail -n 1`
            disp_s="${disp_s/disp-/}"
            disp_f="${disp_f/disp-/}"
            disp_s="${disp_s////}"
            disp_f="${disp_f////}"
            disp_range=disp-\{${disp_s}..${disp_f}\}
            cat>phon.sh<<!
#!/bin/sh
phonopy -f $disp_range/vasprun.xml
!
            chmod 777 phon.sh
            ./phon.sh
            if [ ! -f mesh.conf ]; then
                atom_name=`sed '6q;d' POSCAR`
                cat>mesh.conf <<!
ATOM_NAME = $atom_name
DIM = 2 2 2
MP = 10 10 10
!
            fi
            if [ ! -f total_dos.pdf ]; then
                phonopy -s -p mesh.conf --tolerance=2e-1
            fi
            if [ ! -f thermal_properties.yaml ]; then
                phonopy -t -p mesh.conf --tolerance=2e-1 > thermal.log
            fi
        fi
        cd ..
    fi
    printf "$stat $phonstat \n"
    cd ..
done

"""

    with open(os.path.join(zpe_path, 'gather.sh'), 'w') as f:
        f.write(script)


def clear(zpe_path):
    script = """#!/bin/sh

for dir_name in `ls -d */`
do
    cd $dir_name
    if [ -d phonopy/ ]; then
        cd phonopy/
    if [ -f total_dos.pdf ]; then
        if [ -f thermal_properties.yaml ]; then
            echo deleting disps in $dir_name
            rm -rf disp-*/
        fi
        else 
            echo work is not finished in $dir_name
        fi
        cd ..
    fi
    cd ..
done
    """

    with open(os.path.join(zpe_path, 'clear.sh'), 'w') as f:
        f.write(script)
