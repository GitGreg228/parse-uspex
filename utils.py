import re
import os
import json
import numpy as np

from functools import reduce
from math import gcd
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.cif import CifWriter
from pymatgen.core.structure import IStructure


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
