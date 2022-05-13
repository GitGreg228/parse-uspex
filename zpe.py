import argparse
import os
import json

from tqdm import tqdm
from utils import listdirs, split_poscars, parse_ech, load_ech, collect_zpe, get_convex_hulls, boolean_string
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import IStructure
from inputs import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    parser.add_argument('--ths', nargs='+', default=['0.01'], help='Thresholds on metastable structures')
    parser.add_argument('--tol', type=float, default=0.2, help='Tolerance')
    parser.add_argument('--tol_min', type=float, default=0.01, help='Minimum tolerance')
    parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step')
    parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum tolerance')
    parser.add_argument('--system', type=int, default=1, help='Slurm config in inputs.py')
    parser.add_argument('--press', type=int, default=2000, help='pressure in kBar for subsequent VASP calculations')
    parser.add_argument('--temp', type=int, default=0, help='Temperature on which Convex Hull is calculated')
    parser.add_argument('--plot', type=boolean_string, default=False, help='Make 3d plots (some GUI required)')
    args = parser.parse_args()
    system = get_system(args.system)

    full_dirs = list()
    for dir in listdirs(args.path):
        if 'extended_convex_hull' in os.listdir(dir) and 'extended_convex_hull_POSCARS' in os.listdir(dir):
            full_dirs.append(dir)
    output = ', '.join(full_dirs)
    print(f'Will work in {output}')
    for dir in full_dirs:
        print(f'Working in {dir} ...')
        specific_path = os.path.join(dir, '../Specific/')
        poscars = split_poscars(dir)
        zpe_path = os.path.join(dir, 'ZPE')
        if not os.path.isdir(zpe_path):
            os.mkdir(zpe_path)
            launch_all(zpe_path)
            phonopy(zpe_path)
            gather(zpe_path)
            clear(zpe_path)
            structures = parse_ech(dir, args.ths, poscars, args.tol_min, args.tol_step, args.tol_max, dump_dir=zpe_path)
            for comp in structures.keys():
                # comp_path = os.path.join(zpe_path, comp)
                # if not os.path.isdir(comp_path):
                #     os.mkdir(comp_path)
                for cat in structures[comp].keys():
                    if cat and comp:
                        print(f'Working with the {comp} {cat} structures')
                        for structure in structures[comp][cat]:
                            symmetry = structure['symmetry']
                            space_group = int(symmetry[list(symmetry)[-1]][0])
                            formula = structure['formula']
                            if cat == 'stable':
                                stability = cat
                            else:
                                stability = structure["fitness"]
                            properties = {
                                'formula': formula,
                                'space_group': space_group,
                                'id': structure["id"],
                                'stability': stability,
                                          }
                            name = f'EA{structure["id"]}_{stability}_{formula}_{space_group}'
                            dir_name = os.path.join(zpe_path, name)
                            if not os.path.isdir(dir_name):
                                os.mkdir(dir_name)
                            tmp_structure = IStructure.from_dict(structure['structure'])
                            analyzer = SpacegroupAnalyzer(tmp_structure, symprec=args.tol)
                            primitive = analyzer.get_primitive_standard_structure()
                            poscar = os.path.join(dir_name, 'POSCAR')
                            primitive.to(filename=poscar, fmt='POSCAR')
                            with open(poscar, 'r') as f:
                                lines = f.readlines()
                            with open(poscar, 'w') as f:
                                lines[0] = f'EA{structure["id"]} {stability} {formula} {space_group}\n'
                                f.writelines(lines)
                            compose_potcar(tmp_structure, specific_path, dir_name)
                            get_slurm_script(properties, system, dir_name)
                            write_incar('relaxation', properties, args.press, dir_name)
        else:
            if not os.path.isfile(os.path.join(zpe_path, 'structures.json')):
                _ = parse_ech(dir, args.ths, poscars, args.tol_min, args.tol_step, args.tol_max, dump_dir=zpe_path)
            print(f'Working with ZPE calculations in {zpe_path}')
            structures = load_ech(zpe_path)
            zpe_structures = collect_zpe(zpe_path, structures)
            zpe_structures = get_convex_hulls(zpe_structures, args.temp, args.plot)
            with open(os.path.join(zpe_path, 'zpe_structures.json'), 'w', encoding='utf-8') as f:
                json.dump(zpe_structures, f, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    main()
