import argparse
import os

from tqdm import tqdm
from utils import listdirs, analyze_symmetry, split_poscars, parse_ech, reduce_structures, super_reduce_structures
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
    parser.add_argument('--red', type=str, choices=['none', 'reduce', 'super'], default='reduce', help='How to reduce structures')
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
        poscars_path = os.path.join(dir, 'POSCARS')
        if not os.path.isdir(poscars_path):
            os.mkdir(poscars_path)
        structures = parse_ech(dir, args.ths, poscars, args.tol_min, args.tol_step, args.tol_max, dump_dir=poscars_path, reduce=False)
        if args.red == 'reduce':
            structures = reduce_structures(structures)
        if args.red == 'super':
            structures = super_reduce_structures(structures)
        for comp, value in structures.items():
            comp_path = os.path.join(poscars_path, comp)
            if not os.path.isdir(comp_path):
                os.mkdir(comp_path)
            for cat in structures[comp].keys():
                if cat == 'stable':
                    cat_path = comp_path
                else:
                    cat_path = os.path.join(comp_path, cat)
                if not os.path.isdir(cat_path):
                    os.mkdir(cat_path)
                poscars = os.path.join(poscars_path, f'{comp}_{cat}_POSCARS')
                if cat and comp:
                    with open(poscars, 'w'):
                        pass
                    print(f'Working with the {comp} {cat} structures')
                    for structure in tqdm(structures[comp][cat]):
                        tmp_structure = IStructure.from_dict(structure['structure'])
                        analyzer = SpacegroupAnalyzer(tmp_structure, symprec=args.tol)
                        num = str(analyzer.get_space_group_number())
                        formula = structure['formula']
                        fitness = structure['fitness']
                        if fitness == 0:
                            fname = os.path.join(cat_path, f'{formula}_{num}')
                        else:
                            fit = "{:0.4f}".format(structure['fitness'])
                            fname = os.path.join(cat_path, f'{fit}_{formula}_{num}')
                        if not os.path.isdir(fname):
                            os.mkdir(fname)
                        analyze_symmetry(tmp_structure, args.tol_min, args.tol_step, args.tol_max, save_dir=fname)
                        with open(os.path.join(fname, 'POSCAR'), 'w') as f:
                            lines = structure['poscar']
                            f.writelines(lines)
                        if fitness == 0:
                            lines[0] = f'EA {comp} stable {formula} ({num})\n'
                        else:
                            lines[0] = f'EA {comp} {str(fitness)} {formula} ({num})\n'
                        with open(poscars, 'a') as f:
                            f.writelines(lines)


if __name__ == '__main__':
    main()
