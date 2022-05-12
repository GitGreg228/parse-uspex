import re
import os
import json
import yaml
import numpy as np

from functools import reduce
from tqdm import tqdm
from math import gcd
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.cif import CifWriter
from pymatgen.core.structure import IStructure
from pymatgen.analysis.structure_matcher import StructureMatcher
from scipy.spatial import ConvexHull


def parse_formula(structure):
    formula = structure.formula
    split = formula.split()
    if len(split) == 1:
        return formula
    else:
        split = re.split(r'(\d+)', formula)
        numeric = list()
        alpha = list()
        for word in split:
            if word.isnumeric():
                numeric.append(int(word))
            else:
                alpha.append(word)
        GCD = reduce(gcd, numeric)
        result = str()
        for i in range(len(alpha)):
            result = result + alpha[i]
            if not i == len(numeric):
                if not numeric[i] // GCD == 1:
                    result = result + str(numeric[i] // GCD)
        return result.replace(' ', '')


def formula_from_comp(comp, atoms):
    assert len(comp) == len(atoms)
    formula = str()
    for i, a in enumerate(atoms):
        if comp[i] > 0:
            formula = formula + a
            if comp[i] > 1:
                formula = formula + r'$_{' + str(comp[i]) + r'}$'
    return formula


def analyze_symmetry(structure, tol_min, tol_step, tol_max, save_dir=''):
    prev_num = 0
    tols = dict()
    for tol in np.arange(tol_min, tol_max, tol_step):
        analyzer = SpacegroupAnalyzer(structure, symprec=tol)
        num = analyzer.get_space_group_number()
        if num > prev_num:
            symbol = analyzer.get_space_group_symbol()
            tols["{:0.2f}".format(tol)] = (num, symbol) #str(num) + ' (' + symbol + ')'
            if save_dir:
                prim = analyzer.find_primitive()
                symm = analyzer.get_symmetrized_structure()
                try:
                    path = os.path.join(save_dir, f'{parse_formula(prim)}_{str(num)}.cif')
                    CifWriter(symm, symprec=tol).write_file(path)
                    path = os.path.join(save_dir, f'{parse_formula(prim)}_{str(num)}.vasp')
                    Poscar(prim).write_file(path)
                except TypeError:
                    pass
        prev_num = num
    if save_dir:
        with open(os.path.join(save_dir, 'symm.json'), 'w', encoding='utf-8') as f:
            json.dump(tols, f, ensure_ascii=False, indent=4)
            f.close()
    return tols


def listdirs(rootdir):
    files = list()
    for file in os.listdir(rootdir):
        d = os.path.join(rootdir, file)
        if os.path.isdir(d):
            files = files + listdirs(d)
            if 'results' in d:
                files.append(d)
    return files


def get_comp(lst):
    k = 0
    l = len(lst)
    categories = ['single', 'binary', 'ternary', 'quaternary', 'complex']
    cat = categories[:l]
    for i in lst:
        if i == 0:
            k = k + 1
    if k == 0:
        return cat[-1]
    elif k == 1:
        return cat[-2]
    elif k == 2:
        return cat[-3]
    elif k == 3:
        return cat[-4]


class Structure(object):
    id = int()
    comp = tuple()
    comp_r = tuple()
    comp_cat = str()
    e = float()
    v = float()
    fit = float()
    coords = list()
    stab_cat = str()
    structure = dict()
    formula = str()
    symm = dict()
    subcat = dict()
    X = tuple()
    Y = float()
    poscar = list()

    def __init__(self, line):
        split = line.split()
        self.id = int(split[0])
        self.comp = tuple([int(c) for c in split[split.index('[')+1:split.index(']')]])
        self.comp_r = tuple([c//reduce(gcd, self.comp) for c in self.comp])
        self.e = float(split[split.index(']')+1])
        self.v = float(split[split.index(']')+2])
        self.fit = float(split[split.index(']')+3])
        self.coords = [float(c) for c in split[split.index(']')+5:]]
        self.comp_cat = get_comp(self.comp)
        dim = split.index(']') - split.index('[') - 1
        if dim == 2:
            self.X = (float(split[-2]))
        if dim == 3:
            self.X = (float(split[-3]), float(split[-2]))
        self.Y = (float(split[-1]))

    def stability(self, ths_f):
        if self.fit == 0:
            self.stab_cat = 'stable'
        else:
            for th in ths_f:
                if self.fit <= th:
                    self.stab_cat = str(th)
                    break
        return self.stab_cat

    def struc(self, poscars):
        self.poscar = poscars[str(self.id)]
        with open('tmp_POSCAR', 'w') as f:
            f.writelines(self.poscar)
        self.structure = IStructure.from_file('tmp_POSCAR')
        self.formula = parse_formula(self.structure)
        os.remove('tmp_POSCAR')
        return self.structure

    def symm(self, tol_min, tol_step, tol_max):
        self.symm = analyze_symmetry(self.structure, tol_min, tol_step, tol_max)
        return self.symm

    def subcat(self):
        """
        deprecated
        :return:
        """
        split = re.split(r'(\d+)', self.structure.formula)
        alpha = list()
        for word in split:
            if not word.isnumeric():
                if word.split():
                    alpha.append(word.split()[0])
        self.subcat = '-'.join(alpha)
        return self.subcat

    def as_dict(self):
        struc_dict = {
            'id': self.id,
            'composition': self.comp,
            'composition reduced': self.comp_r,
            'composition category': self.comp_cat,
            'energy, eV/atom': self.e,
            'volume, A^3/atom': self.v,
            'fitness': self.fit,
            'stability category': self.stab_cat,
            'structure': self.structure.as_dict(),
            'formula': self.formula,
            'symmetry': self.symm,
            'convex hull x': self.X,
            'convex hull y': self.Y,
            'poscar': self.poscar
        }
        return struc_dict


def split_poscars(dirname):
    with open(os.path.join(dirname, 'extended_convex_hull_POSCARS'), 'r') as f:
        lines = f.readlines()
    idxs = list()
    eas = list()
    for i, line in enumerate(lines):
        if 'EA' in line:
            idxs.append(i)
            eas.append(line.split()[0])
    idxs.append(len(lines))
    poscars = dict()
    for i in range(len(idxs) - 1):
        poscars[eas[i].replace('EA', '')] = lines[idxs[i]:idxs[i+1]]
    return poscars


def reduce_structures(structures):
    matcher = StructureMatcher()
    for comp in structures.keys():
        for cat in structures[comp].keys():
            if cat != 'stable':
                for structure in structures[comp][cat]:
                    for stable_structure in structures[comp]['stable']:
                        if structure['formula'] == stable_structure['formula'] or comp == 'single':
                            formula = structure['formula']
                            fitness = structure['fitness']
                            s1 = IStructure.from_dict(structure['structure'])
                            s2 = IStructure.from_dict(stable_structure['structure'])
                            if matcher.fit(s1, s2):
                                print(f'{comp} {formula} with fitness {fitness} is the same as stable one')
                                structures[comp][cat].remove(structure)
                for i, structure in enumerate(structures[comp][cat]):
                    for j, other_structure in enumerate(structures[comp][cat]):
                        if i < j:
                            if structure['formula'] == other_structure['formula'] or comp == 'single':
                                f1 = structure['fitness']
                                f2 = other_structure['fitness']
                                if f2 > f1:
                                    formula = structure['formula']
                                    s1 = IStructure.from_dict(structure['structure'])
                                    s2 = IStructure.from_dict(other_structure['structure'])
                                    if matcher.fit(s1, s2):
                                        print(f'{comp} {formula} with fitness {f2} is the same as another one at {f1}')
                                        structures[comp][cat].remove(other_structure)
    return structures


def parse_ech(dirname, ths, poscars, tol_min, tol_step, tol_max, dump_dir='', reduce=True):
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
    if reduce:
        return reduce_structures(new_structures)
    else:
        return new_structures


def load_ech(zpe_path):
    with open(os.path.join(zpe_path, 'structures.json'), 'r', encoding='utf-8') as f:
        structures = json.load(f)
    return structures


def get_structure(structures, id):
    desired_structure = dict()
    for comp in structures.keys():
        for cat in structures[comp].keys():
            for structure in structures[comp][cat]:
                if structure['id'] == id:
                    desired_structure = structure
                    break
            if desired_structure:
                break
        if desired_structure:
            break
    return desired_structure


def get_zpe(thermal):
    with open(thermal, 'r') as f:
        therm = yaml.safe_load(f)
    natom = therm['natom']
    zpe = therm['zero_point_energy'] * 0.01036410
    return zpe/natom


def get_ref_energies(system, zpe_structures):
    ref_energies = list()
    for specie in system:
        energies = list()
        for structure in zpe_structures:
            comp_cat = structure['composition category']
            formula = structure['formula']
            if comp_cat == 'single':
                if specie in formula:
                    energies.append(structure['new energy, eV/atom'])
        ref_energies.append(min(energies))
    return ref_energies


def get_new_y(structure, ref_energies):
    comp = structure['composition reduced']
    e = structure['new energy, eV/atom']
    natoms = sum(comp)
    y = natoms * e
    for i in range(len(comp)):
        y = y - comp[i] * ref_energies[i]
    y = y / natoms
    structure.update({'ZPE convex hull y': y})


def collect_zpe(zpe_path, structures):
    zpe_ids = list()
    dirs = list()
    for fname in os.listdir(zpe_path):
        dir_ = os.path.join(zpe_path, fname)
        if os.path.isdir(dir_):
            if fname.startswith('EA'):
                id = int(fname.split('_')[0].replace('EA', ''))
                zpe_ids.append(id)
                dirs.append(dir_)
    zpe_structures = list()
    zpe_ids_true = list()
    zpe_ids_false = list()
    for i, id in enumerate(zpe_ids):
        thermal = os.path.join(dirs[i], 'phonopy', 'thermal_properties.yaml')
        structure = get_structure(structures, id)
        comp_cat = structure['composition category']
        symmetry = structure['symmetry']
        space_group = symmetry[list(symmetry)[-1]][1]
        formula = structure['formula']
        if structure['stability category'] == 'stable':
            stab_cat = 'stable'
        else:
            stab_cat = structure['fitness']
        if os.path.isfile(thermal):
            zpe = get_zpe(thermal)
            structure.update({'zpe': zpe, 'new energy, eV/atom': structure['energy, eV/atom'] + zpe})
            zpe_structures.append(structure)
            zpe_ids_true.append(f"{space_group}-{formula} (EA{str(id)}, {stab_cat})")
            # print(f"{comp_cat} {space_group}-{formula} (EA{str(id)}, {stab_cat}) has ZPE = {round(zpe, 4)} eV/atom")
        else:
            zpe_ids_false.append(f"{space_group}-{formula} (EA{str(id)}, {stab_cat})")
    print(f'\nWill work with {len(zpe_ids)} structures: {", ".join(zpe_ids_true)}\n')
    print(f'Will NOT work with {len(zpe_ids)} structures: {", ".join(zpe_ids_false)}\n')
    system_cat = 1
    for structure in zpe_structures:
        system_cat = len(structure['composition'])
        break
    if system_cat == 1:
        assert False
    system = list()
    for structure in zpe_structures:
        if system_cat == 3:
            if structure['composition category'] == 'ternary':
                formula = IStructure.from_dict(structure['structure']).formula
                split = re.split(r'(\d+)', formula)
                for specie in split:
                    specie = specie.replace(' ', '')
                    if specie.isalpha():
                        system.append(specie)
                print(f'The program understands this system as ternary {"-".join(system)}')
                break
        if system_cat == 2:
            if structure['composition category'] == 'binary':
                formula = IStructure.from_dict(structure['structure']).formula
                split = re.split(r'(\d+)', formula)
                for specie in split:
                    specie = specie.replace(' ', '')
                    if specie.isalpha():
                        system.append(specie)
                print(f'The program understands this system as binary {"-".join(system)}')
                break
    ref_energies = get_ref_energies(system, zpe_structures)
    for structure in zpe_structures:
        get_new_y(structure, ref_energies)
    return zpe_structures


def get_zpe_ids(zpe_structures):
    zpe_ids = list()
    for structure in zpe_structures:
        zpe_ids.append(structure['id'])
    return zpe_ids


def area(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)


def inside(x1, y1, x2, y2, x3, y3, x, y):
    A = area(x1, y1, x2, y2, x3, y3)
    A1 = area(x, y, x2, y2, x3, y3)
    A2 = area(x1, y1, x, y, x3, y3)
    A3 = area(x1, y1, x2, y2, x, y)
    if A == A1 + A2 + A3:
        return True
    else:
        return False


def get_simplex(point, points):
    points = np.asarray(points)
    x, y = point
    _hull = ConvexHull(points)
    result = None
    for triangle in points[_hull.simplices]:
        if triangle[:,-1].all() <= 0:
            # print(triangle[:,-1])
            x1 = triangle[0, 0]
            x2 = triangle[1, 0]
            x3 = triangle[2, 0]
            y1 = triangle[0, 1]
            y2 = triangle[1, 1]
            y3 = triangle[2, 1]
            if inside(x1, y1, x2, y2, x3, y3, x, y):
                result = triangle
                break
    return result


def get_convex_hulls(zpe_structures):
    old_coords = list()
    for structure in zpe_structures:
        X = structure['convex hull x']
        Y = structure['convex hull y']
        coords = X + [Y]
        old_coords.append(coords)
    hull = ConvexHull(old_coords)
    old_stable = list()
    for i in hull.vertices:
        structure = zpe_structures[i]
        y = structure['convex hull y']
        if y <= 0:
            symmetry = structure['symmetry']
            space_group = symmetry[list(symmetry)[-1]][1]
            formula = structure['formula']
            id = structure['id']
            old_stable.append(f"{space_group}-{formula} (EA{str(id)})")
    print(f'\nOld stable structures are {", ".join(old_stable)}')
    new_coords = list()
    for structure in zpe_structures:
        X = structure['convex hull x']
        Y = structure['ZPE convex hull y']
        coords = X + [Y]
        new_coords.append(coords)
    hull = ConvexHull(new_coords)
    new_stable = list()
    for i in hull.vertices:
        structure = zpe_structures[i]
        y = structure['ZPE convex hull y']
        if y <= 0:
            symmetry = structure['symmetry']
            space_group = symmetry[list(symmetry)[-1]][1]
            formula = structure['formula']
            id = structure['id']
            new_stable.append(f"{space_group}-{formula} (EA{str(id)})")
    print(f'New stable structures are {", ".join(new_stable)}\n')
    """
    for i in range(len(zpe_structures)):
        if i not in hull.vertices:
            structure = zpe_structures[i]
            formula = structure['formula']
            x = structure['convex hull x']
            print(formula, get_simplex(x, new_coords))
    """