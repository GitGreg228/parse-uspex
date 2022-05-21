import re
import os
import json
import yaml
import numpy as np
import pandas as pd

from functools import reduce
from tqdm import tqdm
from math import gcd
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.cif import CifWriter
from pymatgen.core.structure import IStructure
from pymatgen.analysis.structure_matcher import StructureMatcher
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'


def parse_formula(structure):
    formula = structure.formula
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


def get_comp_cat(lst):
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
        self.comp_cat = get_comp_cat(self.comp)
        dim = split.index(']') - split.index('[') - 1
        if dim == 2:
            self.X = (float(split[-2]))
        if dim == 3:
            self.X = (float(split[-3]), float(split[-2]))
        self.Y = (float(split[-1]))

    def stability(self, ths_f):
        if self.fit <= 0:
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

    def symm(self, tol_min=0.01, tol_step=0.01, tol_max=0.5, save_dir=''):
        self.symm = analyze_symmetry(self.structure, tol_min, tol_step, tol_max, save_dir=save_dir)
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
    """
    Reduces structures using pymatgen matcher to find same ones
    :param structures:
    :return:
    """
    matcher = StructureMatcher()
    new_structures = dict()
    for comp in structures.keys():
        new_structures[comp] = dict()
        for cat in structures[comp].keys():
            new_structures[comp][cat] = list()
    comp_dict = dict()
    for comp in structures.keys():
        for cat in structures[comp].keys():
            for structure in structures[comp][cat]:
                comp_r = '-'.join([str(c) for c in structure['composition reduced']])
                if comp_r not in comp_dict.keys():
                    comp_dict[comp_r] = list()
                comp_dict[comp_r].append(structure)
    for comp_r in comp_dict.keys():
        structure_list = comp_dict[comp_r]
        fit_list = [structure['fitness'] for structure in structure_list]
        fittest_idx = fit_list.index(min(fit_list))
        duplicates = list()
        for i, structure in enumerate(structure_list):
            for j, other_structure in enumerate(structure_list):
                f1 = structure['fitness']
                f2 = other_structure['fitness']
                if f2 > f1:
                    formula = structure['formula']
                    s1 = IStructure.from_dict(structure['structure'])
                    s2 = IStructure.from_dict(other_structure['structure'])
                    if matcher.fit(s1, s2):
                        print(f'{formula} with fitness {f2} is the same as another one at {f1}')
                        if i not in duplicates:
                            duplicates.append(i)
                        if j not in duplicates:
                            duplicates.append(j)
        new_structure_list = [structure_list[fittest_idx]]
        for i, structure in enumerate(structure_list):
            if i not in duplicates:
                new_structure_list.append(structure)
        comp_dict[comp_r] = new_structure_list
    for comp_r in comp_dict.keys():
        for structure in comp_dict[comp_r]:
            comp = structure['composition category']
            cat = structure['stability category']
            new_structures[comp][cat].append(structure)
    return new_structures


def super_reduce_structures(structures):
    """
    Reduces structures only by composition
    :param structures:
    :return:
    """
    comp_list = list()
    new_structures = dict()
    for comp in structures.keys():
        new_structures[comp] = dict()
        for cat in structures[comp].keys():
            new_structures[comp][cat] = list()
    for comp in structures.keys():
        for structure in structures[comp]['stable']:
            comp_r = '-'.join([str(c) for c in structure['composition reduced']])
            comp_list.append(comp_r)
            new_structures[comp]['stable'].append(structure)
    for comp in structures.keys():
        for cat in structures[comp].keys():
            if cat != 'stable':
                for structure in structures[comp][cat]:
                    comp_r = '-'.join([str(c) for c in structure['composition reduced']])
                    if not comp_r in comp_list:
                        comp_list.append(comp_r)
                        new_structures[comp][cat].append(structure)
    return new_structures


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


def get_comp(structure, system):
    comp = list()
    for element in system:
        comp.append([specie.symbol for specie in structure.species].count(element))
    GCD = reduce(gcd, comp)
    comp_r = [c//GCD for c in comp]
    return comp, comp_r


def get_x(comp):
    frac_0 = comp[0] / sum(comp)
    frac_1 = comp[1] / sum(comp)
    frac_2 = comp[2] / sum(comp)
    vec_0 = np.asarray([np.sqrt(3), 0]) * frac_0
    vec_1 = np.asarray([- np.sqrt(3) / 2, 1.5]) * frac_1
    vec_2 = np.asarray([- np.sqrt(3) / 2, 1.5]) * frac_2
    return vec_0


class VaspDir(object):
    id = int()
    # struc = Structure

    def __init__(self, dir_path, zpe_path):
        self.id = int(os.path.basename(dir_path).split('_')[0].replace('EA', ''))
        with open(os.path.join(zpe_path, '..', 'extended_convex_hull')) as f:
            lines = f.readlines()
        for line in lines:
            try:
                if self.id == int(line.split()[0]):
                    structure = Structure(line)
                    break
            except ValueError:
                pass
        poscars = {str(self.id): Poscar(IStructure.from_file(os.path.join(dir_path, 'CONTCAR'))).__str__()}
        # print(structure.comp_r, get_x(structure.comp_r), structure.X)
        structure.struc(poscars)
        structure.symm(save_dir=dir_path)
        self.structure = structure.as_dict()


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
    k = 0.01036410  # 1 kJ/mol = 0.01036410 eV
    zpe = therm['zero_point_energy'] * k
    thermal_properties = therm['thermal_properties']
    thermal_dict = dict()
    for each in thermal_properties:
        thermal_dict[str(int(each['temperature']))] = each['free_energy'] * k / natom
    return zpe/natom, thermal_dict


def get_system(structures):
    system = list()
    for structure in structures:
        species = set()
        species = [specie.symbol for specie in structure.species if not (specie in species or species.add(specie))]
        if len(species) > len(system):
            system = species
    comp_cat = get_comp_cat(system)
    print(f'The program understands this system as {comp_cat} {"-".join(system)}')
    return system


def get_ref_energies(system, zpe_structures):
    old_ref_energies = list()
    ref_energies = list()
    ref_gibbs = list()
    for specie in system:
        old_energies = list()
        energies = list()
        gibbs = list()
        for structure in zpe_structures:
            comp_cat = structure['composition category']
            formula = structure['formula']
            if comp_cat == 'single':
                if specie in formula:
                    old_energies.append(structure['energy, eV/atom'])
                    energies.append(structure['new energy, eV/atom'])
                    gibbs.append(structure['G(T)'])
        old_ref_energies.append((min(old_energies)))
        ref_energies.append(min(energies))
        ref_gibbs.append(gibbs[energies.index(min(energies))])  # need some fixes
    return old_ref_energies, ref_energies, ref_gibbs


def get_new_y(structure, ref_energies, old_ref_energies, ref_gibbs):
    comp = structure['composition reduced']
    e = structure['new energy, eV/atom']
    e_old = structure['energy, eV/atom']
    natoms = sum(comp)
    y = natoms * e
    y_old = natoms * e_old
    for i in range(len(comp)):
        y = y - comp[i] * ref_energies[i]
        y_old = y_old - comp[i] * old_ref_energies[i]
    y = y / natoms
    y_old = y_old / natoms
    structure.update({'convex hull y': round(y_old, 4)})
    structure.update({'ZPE convex hull y': round(y, 4)})
    for T in structure['G(T)'].keys():
        e = structure['G(T)'][T]
        y = natoms * e
        for i in range(len(comp)):
            y = y - comp[i] * ref_gibbs[i][T]
        y = y / natoms
        structure.update({f'T = {T} K convex hull y': round(y, 4)})
        # print(structure['formula'], y)


def collect_zpe(zpe_path):
    structures = list()
    structures_ = list()
    zpe_ids = list()
    dirs = list()
    zpe_structures = list()
    for fname in os.listdir(zpe_path):
        dir_ = os.path.join(zpe_path, fname)
        if os.path.isdir(dir_):
            if fname.startswith('EA'):
                dirs.append(dir_)
                structures_.append(IStructure.from_file(os.path.join(dir_, 'CONTCAR')))
    system = get_system(structures_)
    for structure in structures_:
        print(get_comp(structure, system))
    print('Processing VASP and phonopy calculations')
    for i in tqdm(range(len(dirs))):
        vaspdir = VaspDir(dirs[i], zpe_path)
        structures.append(vaspdir.structure)
        zpe_ids.append(vaspdir.id)
    zpe_ids_true = list()
    zpe_ids_false = list()
    for i, id in enumerate(zpe_ids):
        thermal = os.path.join(dirs[i], 'phonopy', 'thermal_properties.yaml')
        structure = structures[i]
        comp_cat = structure['composition category']
        symmetry = structure['symmetry']
        space_group = symmetry[list(symmetry)[-1]][1]
        formula = structure['formula']
        if os.path.isfile(thermal):
            energy = float(Oszicar(os.path.join(dirs[i], 'OSZICAR')).final_energy)
            natoms = IStructure.from_file(os.path.join(dirs[i], 'POSCAR')).num_sites
            zpe, thermal_dict = get_zpe(thermal)
            gibbs = {t: energy/natoms + thermal_dict[t] for t in thermal_dict.keys()}
            structure.update({'energy, eV/atom': energy/natoms, 'ZPE, eV/atom': zpe,
                              'new energy, eV/atom': energy/natoms + zpe, 'F(T)': thermal_dict, 'G(T)': gibbs})
            zpe_ids_true.append(f"{space_group}-{formula} (EA{str(id)})")
            zpe_structures.append(structure)
        else:
            zpe_ids_false.append(f"{space_group}-{formula} (EA{str(id)})")
    print(f'\nWill work with {len(zpe_ids)} structures: {", ".join(zpe_ids_true)}\n')
    print(f'Will NOT work with {len(zpe_ids_false)} structures: {", ".join(zpe_ids_false)}\n')
    old_ref_energies, ref_energies, ref_gibbs = get_ref_energies(system, zpe_structures)
    for structure in zpe_structures:
        get_new_y(structure, ref_energies, old_ref_energies, ref_gibbs)
    return zpe_structures, system


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
    if round(A, 6) == round(A1 + A2 + A3, 6):
        return True
    else:
        return False


def get_simplex(point, points):
    points = np.asarray(points)
    x, y = point
    _hull = ConvexHull(points)
    result = None
    for triangle in points[_hull.simplices]:
        if np.all(triangle[:, -1] <= 0.0) and np.any(np.round(triangle[:, -1], 3) < 0.0):
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


def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, structure, epsilon=1e-6):
    ndotu = planeNormal.dot(rayDirection)
    if abs(ndotu) < epsilon:
        raise RuntimeError(f"{structure['formula']} - no intersection or line is within plane")
    w = rayPoint - planePoint
    si = -planeNormal.dot(w) / ndotu
    Psi = w + si * rayDirection + planePoint
    return Psi


def distance_to_simplex(point, triangle, structure):
    p0, p1, p2 = triangle
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    ux, uy, uz = u = [x1 - x0, y1 - y0, z1 - z0]  # first vector
    vx, vy, vz = v = [x2 - x0, y2 - y0, z2 - z0]  # sec vector
    u_cross_v = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]  # cross product
    planeNormal = np.array(u_cross_v)
    planePoint = np.array([x0, y0, z0])  # Any point on the plane
    rayDirection = np.array([0, 0, 1])
    rayPoint = np.array([point[0], point[1], 0])  # Any point along the ray
    Psi = LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, structure)
    return Psi


class ExtendedConvexHull(object):
    t = float()
    # system = list()
    zpe_structures = list()
    old_coords = list()
    new_coords = list()
    old_stable = list()
    new_stable = list()
    old_hull = ConvexHull
    new_hull = ConvexHull

    def __init__(self, zpe_structures, system, t=0):
        self.system = system # get_system(zpe_structures)
        self.t = int(t)
        self.zpe_structures = zpe_structures
        old_coords = list()
        for structure in zpe_structures:
            X = structure['convex hull x']
            Y = structure['convex hull y']
            coords = list(X) + [Y]
            old_coords.append(coords)
        self.old_coords = np.asarray(old_coords)
        new_coords = list()
        for structure in zpe_structures:
            X = structure['convex hull x']
            if t == 0:
                Y = structure['ZPE convex hull y']
            else:
                Y = structure[f'T = {str(int(t))} K convex hull y']
            coords = list(X) + [Y]
            new_coords.append(coords)
        self.new_coords = np.asarray(new_coords)

    def get_stable(self):
        self.old_stable = list()
        self.new_stable = list()
        hull = ConvexHull(self.old_coords)
        self.old_hull = hull
        for i in hull.vertices:
            structure = self.zpe_structures[i]
            y = structure['convex hull y']
            if y <= 0:
                symmetry = structure['symmetry']
                space_group = symmetry[list(symmetry)[-1]][1]
                formula = structure['formula']
                id = structure['id']
                self.old_stable.append(f"{space_group}-{formula} (EA{str(id)})")
        print(f'\nOld stable structures are {", ".join(self.old_stable)}')
        hull = ConvexHull(self.new_coords)
        self.new_hull = hull
        for i in hull.vertices:
            structure = self.zpe_structures[i]
            if self.t == 0:
                y = structure['ZPE convex hull y']
            else:
                y = structure[f'T = {str(int(self.t))} K convex hull y']
            if y <= 0:
                symmetry = structure['symmetry']
                space_group = symmetry[list(symmetry)[-1]][1]
                formula = structure['formula']
                id = structure['id']
                self.new_stable.append(f"{space_group}-{formula} (EA{str(id)})")
        if self.t:
            print(f'At T = {self.t} K, stable structures are {", ".join(self.new_stable)}\n')
        else:
            print(f'New stable structures are {", ".join(self.new_stable)}\n')

    def get_new_fitness(self):
        for structure in self.zpe_structures:
            x = structure['convex hull x']
            y = structure['convex hull y']
            if self.t == 0:
                y = structure['ZPE convex hull y']
            else:
                y = structure[f'T = {str(int(self.t))} K convex hull y']
            triangle = get_simplex(x, self.new_coords)
            dist = distance_to_simplex(x, triangle, structure)
            new_fitness = y - dist[-1]
            # print(structure['formula'], triangle)
            if new_fitness > 0:
                if round(new_fitness, 4) > 0:
                    new_fitness = round(new_fitness, 4)
                else:
                    new_fitness = '%.2E' % new_fitness
            if self.t == 0:
                structure.update({'ZPE fitness': new_fitness})
            else:
                structure.update({f'T = {str(int(self.t))} K fitness': new_fitness})
                # ructure['formula'], round(y - dist[-1], 4))
        return self.zpe_structures

    def plot(self, ax):
        """
        points = self.old_coords
        x_old = points[:, 0]
        y_old = points[:, 1]
        z_old = points[:, 2]
        ax.scatter(x_old, y_old, z_old, c=z_old, marker='^')
        """
        x_stable, x_unstable = list(), list()
        y_stable, y_unstable = list(), list()
        z_stable, z_unstable = list(), list()
        for i, point in enumerate(self.new_coords):
            x_new, y_new, z_new = point
            structure = self.zpe_structures[i]
            if i in self.new_hull.vertices:
                x_stable.append(x_new)
                y_stable.append(y_new)
                z_stable.append(z_new)
            else:
                x_unstable.append(x_new)
                y_unstable.append(y_new)
                z_unstable.append(z_new)
        ax.scatter(x_stable, y_stable, z_stable, c=z_stable)
        ax.scatter(x_unstable, y_unstable, z_unstable, c='gray', s=0.7)
        xlim = [0.5, 0.5]
        ylim = [0, 0]
        zlim = [0, 0.1]
        single_points = list()
        for i, name in enumerate(self.zpe_structures):
            structure = self.zpe_structures[i]
            if structure['composition category'] == 'single':
                if i in self.new_hull.vertices and self.new_coords[i][-1] == 0:
                    single_points.append(self.new_coords[i])
            else:
                x_new, y_new, z_new = self.new_coords[i]
                if x_new < xlim[0]:
                    xlim[0] = x_new
                if x_new > xlim[1]:
                    xlim[1] = x_new
                if y_new < ylim[0]:
                    ylim[0] = y_new
                if y_new > ylim[1]:
                    ylim[1] = y_new
                if z_new < zlim[0]:
                    zlim[0] = z_new
                if z_new > zlim[1]:
                    zlim[1] = z_new
        for i, name in enumerate(self.zpe_structures):
            x_new, y_new, z_new = self.new_coords[i]
            structure = self.zpe_structures[i]
            formula = formula_from_comp(structure['composition reduced'], self.system)
            # formula = structure['formula']
            if i in self.new_hull.vertices:
                ha = 'left'
                va = 'bottom'
                ax.text(x_new, y_new, z_new, f"{formula}", fontsize=10, color='k', ha=ha, va=va)
            else:
                if self.t == 0:
                    ax.text(x_new, y_new, z_new,
                            f"{formula} ({structure['ZPE fitness']})", fontsize=7, color='gray')
                else:
                    pass
                    ax.text(x_new, y_new, z_new,
                            f"{formula} ({structure[f'T = {str(int(self.t))} K fitness']})", fontsize=7,
                            color='gray')
        single_points = np.asarray(single_points)
        main_x = np.append(single_points[:, 0], single_points[0, 0])
        main_y = np.append(single_points[:, 1], single_points[0, 1])
        main_z = np.append(single_points[:, 2], single_points[0, 2])
        ax.plot(main_x, main_y, main_z, 'k')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)
        points = self.new_coords
        for triangle in points[self.new_hull.simplices]:
            if np.all(triangle[:, -1] <= 0):
                x1 = triangle[:, 0]
                y1 = triangle[:, 1]
                z1 = triangle[:, 2]
                verts = [list(zip(x1, y1, z1))]
                srf = Poly3DCollection(verts, alpha=.15, edgecolor='gray')
                plt.gca().add_collection3d(srf)


def get_convex_hulls(zpe_structures, system, temp, plot):
    if '0' not in temp:
        temp = ['0'] + temp
    for _t in temp:
        ech = ExtendedConvexHull(zpe_structures, system, t=_t)
        ech.get_stable()
        zpe_structures = ech.get_new_fitness()
        if plot:
            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot(111, projection='3d')
            ech.plot(ax)
            ax.set_axis_off()
            plt.tight_layout()
            plt.show()
    return zpe_structures


def save_zpe(zpe_structures, zpe_path, temp):
    save_dict = dict()
    save_dict.update(
        {'EA': [structure['id'] for structure in zpe_structures]})
    save_dict.update({'Space Group': [structure['symmetry'][list(structure['symmetry'])[-1]][-1] for structure in zpe_structures]})
    save_dict.update({'Formula': [structure['formula'] for structure in zpe_structures]})
    save_dict.update({'X1': [structure['convex hull x'][0] for structure in zpe_structures]})
    save_dict.update({'X2': [structure['convex hull x'][1] for structure in zpe_structures]})
    save_dict.update({'E': [round(structure['energy, eV/atom'], 4) for structure in zpe_structures]})
    save_dict.update({'ZPE': [round(structure['ZPE, eV/atom'], 4) for structure in zpe_structures]})
    save_dict.update({'E+ZPE': [round(structure['new energy, eV/atom'], 4) for structure in zpe_structures]})
    for _t in temp:
        if int(_t) > 0:
            save_dict.update({f'F at T = {str(_t)}': [round(structure['F(T)'][str(_t)], 4) for structure in zpe_structures]})
            save_dict.update({f'E+F at T = {str(_t)}': [round(structure['G(T)'][str(_t)], 4) for structure in zpe_structures]})
    save_dict.update({'Y': [round(structure['convex hull y'], 4) for structure in zpe_structures]})
    save_dict.update({'ZPE Y': [round(structure['ZPE convex hull y'], 4) for structure in zpe_structures]})
    for _t in temp:
        if int(_t) > 0:
            save_dict.update({f'Y at T = {str(_t)}': [round(structure[f'T = {str(_t)} K convex hull y'], 4) for
                                                      structure in zpe_structures]})
    save_dict.update({'Old fitness': [structure['fitness'] for structure in zpe_structures]})
    save_dict.update({'ZPE fitness': [round(float(structure['ZPE fitness']), 10) for structure in zpe_structures]})
    for _t in temp:
        if int(_t) > 0:
            save_dict.update({f'Fitness at T = {str(_t)}': [round(float(structure[f'T = {str(_t)} K fitness']), 10) for structure in zpe_structures]})
    name = os.path.join(zpe_path, f'convex_hull.csv')
    df = pd.DataFrame.from_dict(save_dict)
    df.to_csv(name, index=False, header=True, sep='\t', float_format='%.4f')
    print(df)


def plot_ch(df, fit_label, ):
    pass
