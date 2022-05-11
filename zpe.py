import argparse
import os

from tqdm import tqdm
from utils import listdirs, analyze_symmetry, split_poscars, parse_ech
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import IStructure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    parser.add_argument('--ths', nargs='+', default=['0.01'], help='Thresholds on metastable structures')
    parser.add_argument('--tol', type=float, default=0.2, help='Tolerance')
    parser.add_argument('--tol_min', type=float, default=0.01, help='Minimum tolerance')
    parser.add_argument('--tol_step', type=float, default=0.01, help='Tolerance step')
    parser.add_argument('--tol_max', type=float, default=0.5, help='Maximum tolerance')
    args = parser.parse_args()

    full_dirs = list()
    for dir in listdirs(args.path):
        if 'extended_convex_hull' in os.listdir(dir) and 'extended_convex_hull_POSCARS' in os.listdir(dir):
            full_dirs.append(dir)
    output = ', '.join(full_dirs)
    print(f'Will work in {output}')
    for dir in full_dirs:
        print(f'Working in {dir} ...')
        poscars = split_poscars(dir)
        zpe_path = os.path.join(dir, 'ZPE')
        if not os.path.isdir(zpe_path):
            os.mkdir(zpe_path)
        structures = parse_ech(dir, args.ths, poscars, args.tol_min, args.tol_step, args.tol_max, dump_dir=zpe_path)
        for comp in structures.keys():
            comp_path = os.path.join(zpe_path, comp)
            if not os.path.isdir(comp_path):
                os.mkdir(comp_path)
            for cat in structures[comp].keys():
                if cat and comp:
                    print(f'Working with the {comp} {cat} structures')
                    for structure in structures[comp][cat]:
                        symmetry = structure['symmetry']
                        space_group = symmetry[list(symmetry)[-1]][0]
                        formula = structure['formula']
                        if cat == 'stable':
                            name = f'{cat}_{formula}_{space_group}'
                        else:
                            name = f'{structure["fitness"]}_{formula}_{space_group}'
                        dir_name = os.path.join(comp_path, name)
                        if not os.path.isdir(dir_name):
                            os.mkdir(dir_name)
                        tmp_structure = IStructure.from_dict(structure['structure'])
                        analyze_symmetry(tmp_structure, args.tol_min, args.tol_step, args.tol_max, save_dir=dir_name)


if __name__ == '__main__':
    main()
