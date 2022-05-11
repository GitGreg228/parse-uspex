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
    script = f"""#!/bin/sh

    """.format()