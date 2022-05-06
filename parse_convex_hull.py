import argparse
import os
import json

from tqdm import tqdm
from utils import listdirs, Structure, analyze_symmetry, split_poscars
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import IStructure


def parse_ech(dirname, ths, poscars, tol_min, tol_step, tol_max, dump_dir=''):
    structures = dict()
    for category in ['single', 'binary', 'ternary']:
        structures[category] = {'stable': list()}
    ths_f = list()
    for th in ths:
        for key, value in structures.items():
            structures[key][str(th)] = list()
        ths_f.append(float(th))
    ths_f.sort()
    max_th = max(ths_f)
    with open(os.path.join(dirname, 'extended_convex_hull'), 'r') as f:
        lines = f.readlines()
    new_structures = dict()
    for i, line in tqdm(enumerate(lines)):
        if '[' in line:
            S = Structure(line)
            if not S.comp_cat in new_structures:
                new_structures[S.comp_cat] = dict()
            stab_cat = S.stability(ths_f)
            if not stab_cat in new_structures[S.comp_cat]:
                new_structures[S.comp_cat][stab_cat] = list()
            _ = S.struc(poscars)
            _ = S.symm(tol_min, tol_step, tol_max)
            if S.fit > max_th:
                break
            else:
                new_structures[S.comp_cat][stab_cat].append(S.as_dict())
    if dump_dir:
        with open(os.path.join(dump_dir, 'structures.json'), 'w', encoding='utf-8') as f:
            json.dump(new_structures, f, ensure_ascii=False, indent=4)
    return new_structures


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
        poscars_path = os.path.join(dir, 'POSCARS')
        if not os.path.isdir(poscars_path):
            os.mkdir(poscars_path)
        structures = parse_ech(dir, args.ths, poscars, args.tol_min, args.tol_step, args.tol_max, dump_dir=poscars_path)
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
                        # analyzer.get_primitive_standard_structure().to(fmt='poscar', filename=os.path.join(cat_path, fname, 'POSCAR'))
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