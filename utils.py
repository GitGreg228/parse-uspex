import re
import os
import json
import yaml
import numpy as np
import pandas as pd

from copy import deepcopy
from functools import reduce
from itertools import combinations
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
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from adjustText import adjust_text


symtable={1:r'P1',2:r'P-1',3:r'P2',4:r'P2_{1}',5:r'C2',6:r'Pm',7:r'Pc',8:r'Cm',9:r'Cc', 10:r'P2/m',
         11:r'P2_{1}/m',12:r'C2/m',13:r'P2/c',14:r'P2_{1}/c',15:r'C2/c',16:r'P222',17:r'P222_{1}',18:r'P2_{1}2_{1}2',19:r'P2_{1}2_{1}2_{1}',20:r'C222_{1}',
         21:r'C222', 22:r'F222',23:r'I222',24:r'I2_{1}2_{1}2_{1}',25:r'Pmm2',26:r'Pmc2_{1}',27:r'Pcc2',28:r'Pma2',29:r'Pca2_{1}', 30:r'Pnc2',
         31:r'Pmn2_{1}',32:r'Pba2',33:r'Pna2_{1}',34:r'Pnn2',35:r'Cmm2',36:r'Cmc2_{1}',37:r'Ccc2',38:r'Amm2',39:r'Aem2',40:r'Ama2',
         41:r'Aea2',42:r'Fmm2',43:r'Fdd2',44:r'Imm2',45:r'Iba2',46:r'Ima2',47:r'Pmmm',48:r'Pnnn',49:r'Pccm',50:r'Pban',
         51:r'Pmma',52:r'Pnna',53:r'Pmna',54:r'Pcca',55:r'Pbam',56:r'Pccn',57:r'Pbcm',58:r'Pnnm',59:r'Pmmn',60:r'Pbcn',
         61:r'Pbca',62:r'Pnma',63:r'Cmcm',64:r'Cmce',65:r'Cmmm',66:r'Cccm',67:r'Cmme',68:r'Ccce',69:r'Fmmm',70:r'Fddd',
         71:r'Immm',72:r'Ibam',73:r'Ibca',74:r'Imma',75:r'P4',76:r'P4_{1}',77:r'P4_{2}',78:r'P4_{3}',79:r'I4',80:r'I4_{1}',
         81:r'P-4',82:r'I-4',83:r'P4/m',84:r'P4_{2}/m',85:r'P4/n',86:r'P4_{2}/n',87:r'I4/m',88:r'I4_{1}/a',89:r'P422',90:r'P42_{1}2',
         91:r'P4_{1}22', 92:r'P4_{1}2_{1}2',93:r'P4_{2}22',94:r'P4_{2}2_{1}2',95:r'P4_{3}22',96:r'P4_{3}2_{1}2',97:r'I422',98:r'I4_{1}22',99:r'P4mm',100:r'P4bm',
         101:r'P4_{2}cm',102:r'P4_{2}nm',103:r'P4cc',104:r'P4nc',105:r'P4_{2}mc',106:r'P4_{2}bc',107:r'I4mm',108:r'I4cm',109:r'I4_{1}md',110:r'I4_{1}cd',
         111:r'P-42m',112:r'P-42c',113:r'P-42_{1}m',114:r'P-42_{1}c',115:r'P-4m2',116:r'P-4c2',117:r'P-4b2',118:r'P-4n2',119:r'I-4m2',120:r'I-4c2',
         121:r'I-42m',122:r'I-42d',123:r'P4/mmm',124:r'P4/mcc',125:r'P4/nbm',126:r'P4/nnc',127:r'P4/mbm',128:r'P4/mnc',129:r'P4/nmm',130:r'P4/ncc',
         131:r'P4_{2}/mmc',132:r'P4_{2}/mcm',133:r'P4_{2}/nbc',134:r'P4_{2}/nnm',135:r'P4_{2}/mbc',136:r'P4_{2}/mnm',137:r'P4_{2}/nmc',138:r'P4_{2}/ncm',139:r'I4/mmm',140:r'I4/mcm',
         141:r'I4_{1}/amd',142:r'I4_{1}/acd',143:r'P3',144:r'P31',145:r'P3_{2}',146:r'R3',147:r'P-3',148:r'R-3',149:r'P312',150:r'P321',
         151:r'P3112',152:r'P3_{1}21',153:r'P3_{2}12',154:r'P3_{2}21',155:r'R32',156:r'P3m1',157:r'P31m', 158:r'P3c1',159:r'P31c',160:r'R3m',
         161:r'R3c',162:r'P-31m',163:r'P-31c',164:r'P-3m1',165:r'P-3c1',166:r'R-3m',167:r'R-3c',168:r'P6',169:r'P6_{1}',170:r'P6_{5}',
         171:r'P6_{2}',172:r'P6_{4}',173:r'P6_{3}',174:r'P-6',175:r'P6/m',176:r'P6_{3}/m',177:r'P622',178:r'P6_{1}22',179:r'P6_{5}22',180:r'P6_{2}22',
         181:r'P6_{4}22',182:r'P6_{3}22',183:r'P6mm',184:r'P6cc',185:r'P6_{3}cm',186:r'P6_{3}mc',187:r'P-6m2',188:r'P-6c2',189:r'P-62m',190:r'P-62c',
         191:r'P6/mmm',192:r'P6/mcc',193:r'P6_{3}/mcm',194:r'P6_{3}/mmc',195:r'P23',196:r'F23',197:r'I23',198:r'P2_{1}3',199:r'I2_{1}3',200:r'Pm-3',
         201:r'Pn-3',202:r'Fm-3',203:r'Fd-3',204:r'Im-3',205:r'Pa-3',206:r'Ia-3',207:r'P432',208:r'P4_{2}32',209:r'F432',210:r'F4_{1}32',
         211:r'I432',212:r'P4_{3}32',213:r'P4_{1}32',214:r'I4_{1}32',215:r'P-43m',216:r'F-43m',217:r'I-43m',218:r'P-43n',219:r'F-43c',220:r'I-43d',
         221:r'Pm-3m',222:r'Pn-3n',223:r'Pm-3n',224:r'Pn-3m',225:r'Fm-3m',226:r'Fm-3c',227:r'Fd-3m',228:r'Fd-3c',229:r'Im-3m',230:r'Ia-3d'}


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


def get_X(comp_r):
    if len(comp_r) == 2:
        return np.round(comp_r[1] / (comp_r[0] + comp_r[1]), 6)
    elif len(comp_r) == 3:
        if not comp_r[0] and not comp_r[1]:
            n_AB = 0
            n_ABC = 0
        else:
            n_AB = comp_r[1] / (comp_r[0] + comp_r[1])
            n_ABC = 1 - comp_r[2] / (comp_r[0] + comp_r[1] + comp_r[2])
        x1, y1 = n_AB * np.sqrt(3) * np.cos(np.pi / 3) - np.sqrt(3) / 2, n_AB * np.sqrt(3) * np.sin(np.pi / 3) - 0.5
        x2, y2 = np.sqrt(3) / 2, -0.5
        a = n_AB * np.sqrt(3)
        b = np.sqrt(np.square(x2 - x1) + np.square(y2 - y1))
        c = np.sqrt(3)
        alpha = np.arccos((np.square(c) + np.square(b) - np.square(a))/(2*b*c))
        X1, X2 = np.sqrt(3) - b * n_ABC * np.cos(alpha) - np.sqrt(3) / 2, b * n_ABC * np.sin(alpha) - 0.5
        return np.round(X1, 6), np.round(X2, 6)


class VaspDir(object):
    id = int()
    system = list()
    # struc = Structure

    def __init__(self, dir_path, system, tol_min=0.01, tol_step=0.01, tol_max=0.5):
        self.id = int(os.path.basename(dir_path).split('_')[0].replace('EA', ''))
        self.system = system
        self.structure = dict()
        pmg_struc = IStructure.from_file(os.path.join(dir_path, 'CONTCAR'))
        self.structure['structure'] = pmg_struc.as_dict()
        self.structure['composition'], self.structure['composition reduced'] = get_comp(pmg_struc, system)
        self.structure['convex hull x'] = get_X(self.structure['composition reduced'])
        self.structure['composition category'] = get_comp_cat(self.structure['composition reduced'])
        self.structure['symmetry'] = analyze_symmetry(pmg_struc, tol_min, tol_step, tol_max, save_dir=dir_path)
        self.structure['formula'] = parse_formula(pmg_struc)
        self.structure['id'] = int(os.path.basename(dir_path).split('_')[0].replace('EA', ''))


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


def get_system(structures, p=True):
    system = list()
    for structure in structures:
        if isinstance(structure, dict):
            _structure = IStructure.from_dict(structure['structure'])
        else:
            _structure = structure
        species = [specie.symbol for specie in _structure.species]
        for specie in species:
            if specie not in system:
                system.append(specie)
    if 'H' in system:
        system.append(system.pop(system.index('H')))
    comp_cat = get_comp_cat(system)
    if p:
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
    print('Processing VASP and phonopy calculations')
    for i in tqdm(range(len(dirs))):
        vaspdir = VaspDir(dirs[i], system)
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


def inside(x1, y1, x2, y2, x3, y3, x, y, n):
    A = area(x1, y1, x2, y2, x3, y3)
    A1 = area(x, y, x2, y2, x3, y3)
    A2 = area(x1, y1, x, y, x3, y3)
    A3 = area(x1, y1, x2, y2, x, y)
    if round(A, n) == round(A1 + A2 + A3, n):
        return True
    else:
        return False


def get_simplex(point, points, formula):
    points = np.asarray(points)
    x, y = point
    _hull = ConvexHull(points)
    result = np.zeros((3, 3))
    n = 6
    f = False
    while result.any() == 0 and n > 0:
        for triangle in points[_hull.simplices]:
            if round(get_normal(triangle)[-1], 3) != 0:
                if np.all(triangle[:, -1] <= 0.0) and np.any(np.round(triangle[:, -1], 3) < 0.0):
                    x1 = triangle[0, 0]
                    x2 = triangle[1, 0]
                    x3 = triangle[2, 0]
                    y1 = triangle[0, 1]
                    y2 = triangle[1, 1]
                    y3 = triangle[2, 1]
                    if inside(x1, y1, x2, y2, x3, y3, x, y, n):
                        f = True
                        result = triangle
                        break
        if f:
            break
        n = n - 1
        print(f'Problems with {formula}, decreasing n')
    return result


def get_simplex_2d(point, points):
    points = points[ConvexHull(points).vertices]
    x, y = points[:, 0], points[:, 1]
    x = x[y <= 0]
    y = y[y <= 0]
    y = y[np.argsort(x)]
    x = x[np.argsort(x)]
    for i in range(np.size(x) - 1):
        if x[i] <= point <= x[i + 1]:
            result = np.array([[x[i], y[i]], [x[i+1], y[i+1]]])
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


def get_normal(triangle):
    p0, p1, p2 = triangle
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    ux, uy, uz = u = [x1 - x0, y1 - y0, z1 - z0]  # first vector
    vx, vy, vz = v = [x2 - x0, y2 - y0, z2 - z0]  # sec vector
    u_cross_v = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]  # cross product
    planeNormal = np.array(u_cross_v)
    return planeNormal


def distance_to_simplex(point, triangle, structure):
    planeNormal = get_normal(triangle)
    planePoint = np.array(triangle[0])  # Any point on the plane
    rayDirection = np.array([0, 0, 1])
    rayPoint = np.array([point[0], point[1], 0])  # Any point along the ray
    Psi = LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, structure)
    return Psi


def distance_to_simplex_2d(x, y, segment):
    x1, y1, x2, y2 = segment[0, 0], segment[0, 1], segment[1, 0], segment[1, 1]
    a = (y2 - y1) / (x2 - x1)
    b = y1 - a * x1
    return a * x + b


def get_stable_(hull, zpe_structures):
    stable = list()
    stable_idx = list()
    for i in hull.vertices:
        if i not in stable_idx:
            stable_idx.append(i)
            structure = zpe_structures[i]
            y = structure['convex hull y']
            if y <= 0:
                symmetry = structure['symmetry']
                space_group = symmetry[list(symmetry)[-1]][1]
                formula = structure['formula']
                id = structure['id']
                stable.append(f"{space_group}-{formula} (EA{str(id)})")
    return stable


def get_simplex_comp(simplex, structures):
    systems = list()
    for idx in simplex:
        formula = structures[idx]['formula']
        split = re.split(r'(\d+)', formula)
        system = str()
        for word in split:
            if word.isalpha():
                system = system + word
        if get_comp_cat(structures[idx]['composition reduced']) != 'single':
            systems.append(system)
    return set(systems)


def get_triangle(zpe_structures, system, ax):
    x1 = [structure['convex hull x'][0] for structure in zpe_structures]
    x2 = [structure['convex hull x'][1] for structure in zpe_structures]
    minx1, maxx1, minx2, maxx2 = 0.5, -0.5, 1, -0.5
    for i in range(len(x1)):
        comp = get_comp_cat(zpe_structures[i]['composition reduced'])
        if not comp == 'single':
            if x1[i] > maxx1:
                maxx1 = x1[i]
            if x1[i] < minx1:
                minx1 = x1[i]
            if x2[i] > maxx2:
                maxx2 = x2[i]
            if x2[i] < minx2:
                minx2 = x2[i]
    ax.plot((minx1 - 0.05, 0.866), (-0.5, -0.5), 'k', zorder=2)
    ax.plot((minx1 - 0.05, minx1 - 0.1), (-0.5, -0.5), 'k:', zorder=2)
    r = 0.866 - minx1
    ax.plot((0.866 - (r + 0.05) * 0.5, 0.866), ((r + 0.05) * 0.866 - 0.5, -0.5), 'k', zorder=2)
    ax.plot((0.866 - (r + 0.05) * 0.5, 0.866 - (r + 0.1) * 0.5), ((r + 0.05) * 0.866 - 0.5, (r + 0.1) * 0.866 - 0.5), 'k:', zorder=2)
    ax.plot((0.866 - (r + 0.1) * 0.5, minx1 - 0.1), ((r + 0.1) * 0.866 - 0.5, -0.5), 'k', zorder=2)
    ax.scatter((0.866 - (r + 0.1) * 0.5, minx1 - 0.1), ((r + 0.1) * 0.866 - 0.5, -0.5), c='k', marker='D', zorder=1)
    ax.annotate(system[0], xy=(minx1 - 0.1, -0.5), xytext=(minx1 - 0.1 - 0.015, -0.5 - 0.015), va='top', ha='right', fontweight='bold')
    ax.annotate(system[2], xy=(0.866, -0.5), xytext=(0.866+0.015, -0.5-0.015), va='top', ha='left', fontweight='bold')
    ax.annotate(system[1], xy=(0.866 - (r + 0.1) * 0.5, (r + 0.1) * 0.866 - 0.5),  xytext=(0.866 - (r + 0.1) * 0.5, (r + 0.1) * 0.866 - 0.5 + 0.015), va='bottom', ha='center', fontweight='bold')
    minx1, maxx1, minx2, maxx2 = minx1 - 0.15, maxx1 + 0.1, minx2 - 0.1, maxx2 + 0.15
    ax.xlim(minx1, maxx1)
    ax.ylim(minx2, maxx2)
    return (minx1, maxx1), (minx2, maxx2)


def get_phase_diagram(points, hull, xlim, ylim, ax):
    pairs = list()
    for simplex in hull.simplices:
        point = points[simplex]
        point0, point1, point2 = point[0][:2], point[1][:2], point[2][:2]
        _points = list()
        for point in (point0, point1, point2):
            if xlim[0] < point[0] < xlim[1] and ylim[0] < point[1] < ylim[1]:
                _points.append(point)
        for pair in combinations(_points, 2):
            ax.plot((pair[0][0], pair[1][0]), (pair[0][1], pair[1][1]), linestyle=':', color='gray', zorder=0)


class ExtendedConvexHull(object):
    t = float()
    # system = list()
    zpe_structures = list()
    old_coords = list()
    new_coords = list()
    old_stable = list()
    new_stable = list()
    dim = int()
    old_hull = ConvexHull
    new_hull = ConvexHull

    def __init__(self, zpe_structures, system, t=0, antiseeds=list()):
        self.system = system # get_system(zpe_structures)
        self.t = int(t)
        self.zpe_structures = zpe_structures
        old_coords = list()
        self.dim = len(system)
        for structure in zpe_structures:
            if structure['id'] not in antiseeds:
                X = structure['convex hull x']
                Y = structure['convex hull y']
                if Y <= 0:
                    if isinstance(X, float):
                        coords = [X, Y]
                        self.dim = 2
                    else:
                        coords = list(X) + [Y]
                        self.dim = len(X) + 1
                    old_coords.append(coords)
        self.old_coords = np.asarray(old_coords)
        self.old_hull = ConvexHull(self.old_coords)
        new_coords = list()
        for structure in zpe_structures:
            if structure['id'] not in antiseeds:
                X = structure['convex hull x']
                if t == 0:
                    Y = structure['ZPE convex hull y']
                else:
                    Y = structure[f'T = {str(int(t))} K convex hull y']
                if Y <= 0:
                    if isinstance(X, float):
                        coords = [X, Y]
                    else:
                        coords = list(X) + [Y]
                    new_coords.append(coords)
        self.new_coords = np.asarray(new_coords)
        self.new_hull = ConvexHull(self.new_coords)

    def get_new_fitness(self, sections=list()):
        section_id = list()
        section_old_fitness = list()
        section_fitness = list()
        for section in sections:
            if len(section) > 2:
                for structure in section:
                    if structure['id'] not in section_id:
                        section_id.append(structure['id'])
                        section_old_fitness.append(structure['fitness'])
                        if self.t == 0:
                            section_fitness.append(structure['ZPE fitness'])
                        else:
                            section_fitness.append(structure[f'T = {str(int(self.t))} K fitness'])
        self.old_stable = list()
        for structure in self.zpe_structures:
            if structure['id'] in section_id:
                idx = section_id.index(structure['id'])
                structure['fitness'] = section_old_fitness[idx]
            else:
                x = structure['convex hull x']
                y = structure['convex hull y']
                if self.dim == 3:
                    triangle = get_simplex(x, self.old_coords, structure['formula'])
                    if triangle.any() == 0:
                        print(f'No triangle for {structure["formula"]}')
                    dist = distance_to_simplex(x, triangle, structure)
                    new_fitness = y - dist[-1]
                elif self.dim == 2:
                    segment = get_simplex_2d(x, self.old_coords)
                    dist = distance_to_simplex_2d(x, y, segment)
                    new_fitness = y - dist
                structure.update({'fitness': np.round(new_fitness, 4)})
            if structure['fitness'] <= 0:
                symmetry = structure['symmetry']
                space_group = symmetry[list(symmetry)[-1]][1]
                formula = structure['formula']
                id = structure['id']
                self.old_stable.append(f"{space_group}-{formula} (EA{str(id)})")
        self.new_stable = list()
        for structure in self.zpe_structures:
            if structure['id'] in section_id:
                idx = section_id.index(structure['id'])
                new_fitness = section_fitness[idx]
                if self.t == 0:
                    structure['ZPE fitness'] = new_fitness
                else:
                    structure[f'T = {str(int(self.t))} K fitness'] = new_fitness
            else:
                x = structure['convex hull x']
                if self.t == 0:
                    y = structure['ZPE convex hull y']
                else:
                    y = structure[f'T = {str(int(self.t))} K convex hull y']
                if self.dim == 3:
                    triangle = get_simplex(x, self.new_coords, structure['formula'])
                    if triangle.any() == 0:
                        print(f'No triangle for {structure["formula"]}')
                    dist = distance_to_simplex(x, triangle, structure)
                    new_fitness = y - dist[-1]
                elif self.dim == 2:
                    segment = get_simplex_2d(x, self.new_coords)
                    dist = distance_to_simplex_2d(x, y, segment)
                    new_fitness = y - dist
                if new_fitness > 0:
                    if round(new_fitness, 4) > 0:
                        new_fitness = round(new_fitness, 4)
                    else:
                        new_fitness = '%.2E' % new_fitness
                if self.t == 0:
                    structure.update({'ZPE fitness': new_fitness})
                else:
                    structure.update({f'T = {str(int(self.t))} K fitness': new_fitness})
            if float(new_fitness) <= 0:
                symmetry = structure['symmetry']
                space_group = symmetry[list(symmetry)[-1]][1]
                formula = structure['formula']
                id = structure['id']
                self.new_stable.append(f"{space_group}-{formula} (EA{str(id)})")
        print(f'\nOld stable structures are {", ".join(self.old_stable)}')
        if self.t:
            print(f'At T = {self.t} K, stable structures are {", ".join(self.new_stable)}\n')
        else:
            print(f'New stable structures are {", ".join(self.new_stable)}\n')
        return self.zpe_structures

    def plot(self, zpe_path=str(), plot=False, th=0.03, press=None):
        """
        points = self.old_coords
        x_old = points[:, 0]
        y_old = points[:, 1]
        z_old = points[:, 2]
        ax.scatter(x_old, y_old, z_old, c=z_old, marker='^')
        """
        if self.dim == 2:
            plt.figure(figsize=(16, 9))
            x_stable, x_unstable = list(), list()
            y_stable, y_unstable = list(), list()
            xlim = [1, 0]
            ylim = [0.05, 0]
            x_labels = list()
            x_ticks = list()
            points = list()
            labels = list()
            texts = list()
            for i, structure in enumerate(self.zpe_structures):
                y = structure['ZPE convex hull y'] if self.t == 0 else structure[
                    f'T = {str(self.t)} K convex hull y']
                fit = round(float(structure['ZPE fitness']), 8) if self.t == 0 else round(
                    float(structure[f'T = {str(int(self.t))} K fitness']), 8)
                x_new, y_new = [structure['convex hull x'], y]
                symm = structure['symmetry'][list(structure['symmetry'].keys())[-1]][-1]
                comp_cat = get_comp_cat(structure['composition reduced'])
                if '-' in symm:
                    split = symm.split('-')
                    symm = split[0] + '\overline{' + split[1][0] + '}' + split[1][1:]
                label = r'$' + symm + r'$'
                points.append([x_new, y_new])
                if fit <= 0:
                    x_stable.append(x_new)
                    y_stable.append(y_new)
                elif 0 < fit <= th:
                    x_unstable.append(x_new)
                    y_unstable.append(y_new)
                    label = label + f'\n({np.round(fit, 3)})'
                labels.append(label)
                if y < 0 and comp_cat != 'single' and fit <= th:
                    if x_new < xlim[0]:
                        xlim[0] = x_new
                    if x_new > xlim[1]:
                        xlim[1] = x_new
                    if y_new < ylim[0]:
                        ylim[0] = y_new
                    if y_new > ylim[1]:
                        ylim[1] = y_new
            xlim[0] = 0.99 * xlim[0]
            xlim[1] = 1.01 * xlim[1]
            for i, structure in enumerate(self.zpe_structures):
                x_new = structure['convex hull x']
                y_new = structure['ZPE convex hull y'] if self.t == 0 else structure[
                    f'T = {str(self.t)} K convex hull y']
                fit = round(float(structure['ZPE fitness']), 8) if self.t == 0 else round(
                    float(structure[f'T = {str(int(self.t))} K fitness']), 8)
                comp_r = structure['composition reduced']
                if 0.95 * xlim[0] <= x_new <= xlim[1]:
                    if fit <= 0:
                        texts.append(plt.text(x_new, y_new, labels[i], color='k', fontsize=20, ha='left', va='top'))
                    elif 0 < fit <= th:
                        texts.append(
                            plt.text(x_new, y_new, labels[i], color='darkslategrey', fontsize=15, ha='right', va='bottom'))
                    if x_new not in x_ticks:
                        if y_new <= 0 and fit <= th:
                            # if comp_r[0] < 3:
                            x_ticks.append(x_new)
                            x_labels.append(formula_from_comp(comp_r, self.system))
            points = np.asarray(points)
            hull = ConvexHull(points)
            for simplex in hull.simplices:
                plt.plot(points[simplex, 0], points[simplex, 1], linestyle='--', color='darkgreen', lw=3)
            plt.plot(points[simplex, 0], points[simplex, 1], linestyle='--', color='darkgreen', lw=3,
                     label='Convex hull')
            plt.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'ro', markersize=10,
                     label='At convex hull')
            plt.scatter(x_unstable, y_unstable, s=100, color='darkorange', label='Not at convex hull + fitness')
            plt.legend(loc='upper left', fontsize=20)
            plt.xticks(x_ticks, labels=x_labels, rotation=90, fontsize=20)
            plt.yticks(fontsize=20)
            plt.xlim(xlim[0], xlim[1])
            plt.ylabel(r'$E$, eV/atom', fontsize=20)
            adjust_text(texts, arrowprops=dict(arrowstyle='-'))
            plt.grid(axis='x')
            if zpe_path:
                if self.t == 0:
                    plt.title(f"{'-'.join(self.system)} convex hull (ZPE) at P = {press} GPa", fontsize=20)
                    fname = os.path.join(zpe_path, f"{'-'.join(self.system)}_convex_hull_ZPE_P={press}GPa.pdf")
                    plt.savefig(fname, bbox_inches='tight')
                    plt.savefig(fname.replace('.pdf', '.png'), dpi=300)
                else:
                    plt.title(f"{'-'.join(self.system)} convex hull (T = {self.t} K, P = {press} GPa)", fontsize=20)
                    fname = os.path.join(zpe_path, f"{'-'.join(self.system)}_convex_hull_T={self.t}K_P={press}GPa.pdf")
                    plt.savefig(fname, bbox_inches='tight')
                    plt.savefig(fname.replace('.pdf', '.png'), dpi=300)
        elif self.dim == 3:
            x_stable, x_unstable = list(), list()
            y_stable, y_unstable = list(), list()
            z_stable, z_unstable = list(), list()
            xlim = [0.5, 0.5]
            ylim = [0, 0]
            zlim = [0, 0.1]
            single_points = list()
            fitness_lst = list()
            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot(111, projection='3d')
            for i, structure in enumerate(self.zpe_structures):
                if self.t == 0:
                    y = structure['ZPE convex hull y']
                    fit = round(float(structure['ZPE fitness']), 8)
                else:
                    y = structure[f'T = {str(self.t)} K convex hull y']
                    fit = round(float(structure[f'T = {str(int(self.t))} K fitness']), 8)
                x_new, y_new, z_new = list(structure['convex hull x']) + [y]
                if fit <= 0:
                    x_stable.append(x_new)
                    y_stable.append(y_new)
                    z_stable.append(z_new)
                elif 0 < fit <= th:
                    x_unstable.append(x_new)
                    y_unstable.append(y_new)
                    z_unstable.append(z_new)
                    fitness_lst.append(fit)
                if y == 0:
                    single_points.append([x_new, y_new, z_new])
                else:
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
                formula = formula_from_comp(structure['composition reduced'], self.system)
                if fit <= 0:
                    ax.text(x_new, y_new, z_new, f"{formula}", fontsize=10, color='k')
                else:
                    ax.text(x_new, y_new, z_new, f"{formula} ({fit})", fontsize=7, color='gray')
            ax.scatter(x_stable, y_stable, z_stable, c=z_stable)
            ax.scatter(x_unstable, y_unstable, z_unstable, c='gray', s=0.7)
            single_points = np.asarray(single_points)
            main_x = np.append(single_points[:, 0], single_points[0, 0])
            main_y = np.append(single_points[:, 1], single_points[0, 1])
            main_z = np.append(single_points[:, 2], single_points[0, 2])
            ax.plot(main_x, main_y, main_z, 'k')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_zlim(zlim)
            points = self.new_coords
            for simplex in self.new_hull.simplices:
                triangle = points[simplex]
                if np.all(triangle[:, -1] <= 0):
                    x1 = triangle[:, 0]
                    y1 = triangle[:, 1]
                    z1 = triangle[:, 2]
                    verts = [list(zip(x1, y1, z1))]
                    srf = Poly3DCollection(verts, alpha=.15, edgecolor='gray')
                    plt.gca().add_collection3d(srf)
            ax.set_axis_off()
            if plot:
                plt.show()
            # 2d projection
            plt.figure(figsize=(12, 12))
            font = {'weight': 'normal', 'size': 25}
            matplotlib.rc('font', **font)
            xlim, ylim = get_triangle(self.zpe_structures, self.system, plt)
            plt.scatter(x_stable, y_stable, c='k', marker='D', zorder=4)
            im = plt.scatter(x_unstable, y_unstable, c=fitness_lst, cmap='plasma', vmin=0, vmax=th, zorder=3)
            get_phase_diagram(self.new_coords, self.new_hull, xlim, ylim, plt)
            texts = list()
            appr_labels = list()
            for structure in self.zpe_structures:
                if get_comp_cat(structure['composition reduced']) != 'single':
                    fit = round(float(structure['ZPE fitness']), 8) if self.t == 0 \
                        else round(float(structure[f'T = {str(int(self.t))} K fitness']), 8)
                    label = formula_from_comp(structure['composition reduced'], self.system)
                    if fit <= 0:
                        appr_labels.append(label)
                        texts.append(plt.text(structure['convex hull x'][0], structure['convex hull x'][1],
                                              label, fontsize=20, color='k', zorder=2))
            for structure in self.zpe_structures:
                if get_comp_cat(structure['composition reduced']) != 'single':
                    fit = round(float(structure['ZPE fitness']), 8) if self.t == 0 \
                        else round(float(structure[f'T = {str(int(self.t))} K fitness']), 8)
                    label = formula_from_comp(structure['composition reduced'], self.system)
                    if 0 < fit <= th:
                        if label not in appr_labels:
                            appr_labels.append(label)
                            label = label + f'\n({round(fit, 3)})'
                            texts.append(plt.text(structure['convex hull x'][0], structure['convex hull x'][1],
                                                  label, fontsize=10, color='darkslategrey', zorder=3))
            adjust_text(texts, arrowprops=dict(arrowstyle='-'))
            plt.axis('off')
            # cbar = plt.colorbar(im, orientation='horizontal', ticks=[0, 0.005, 0.01, 0.015, 0.02])
            # cbar.ax.set_xticklabels([0, 0.005, 0.01, 0.015, 0.02], fontsize=25)
            # cbar.ax.set_title('$E_{\mathrm{hull}}$, eV/atom')
            if self.t == 0:
                plt.title(f'ZPE convex hull at P = {press} GPa', fontweight='bold')
                fname = os.path.join(zpe_path, f"{'-'.join(self.system)}_convex_hull_ZPE_P={press}GPa.pdf")
            else:
                plt.title(f'T = {self.t} K, P = {press} GPa convex hull', fontweight='bold')
                fname = os.path.join(zpe_path, f"{'-'.join(self.system)}_convex_hull_T={self.t}K_P={press}GPa.pdf")
            plt.tight_layout()
            plt.savefig(fname)
            plt.savefig(fname.replace('.pdf', '.png'), dpi=300)


def get_binary_systems(zpe_structures):
    sections = [list(), list(), list()]
    systems = [list(), list(), list()]
    structures = deepcopy(zpe_structures)
    for structure in structures:
        if len(structure['composition reduced']) == 3:
            if structure['composition reduced'][0] == 0:
                _structure = deepcopy(structure)
                del _structure['composition reduced'][0]
                _structure['convex hull x'] = get_X(_structure['composition reduced'])
                sections[0].append(_structure)
        if len(structure['composition reduced']) == 3:
            if structure['composition reduced'][1] == 0:
                _structure = deepcopy(structure)
                del _structure['composition reduced'][1]
                _structure['convex hull x'] = get_X(_structure['composition reduced'])
                sections[1].append(_structure)
                # print(1, _structure['formula'], _structure['composition reduced'])
        if len(structure['composition reduced']) == 3:
            if structure['composition reduced'][2] == 0:
                _structure = deepcopy(structure)
                del _structure['composition reduced'][2]
                _structure['convex hull x'] = get_X(_structure['composition reduced'])
                sections[2].append(_structure)
    for i, section in enumerate(sections):
        systems[i] = get_system(section)
    return sections, systems


def get_convex_hulls(zpe_structures, system, zpe_path, temp, plot, press):
    press = int(press / 10)
    antiseeds = dict()
    sections = list()
    for _t in temp:
        antiseeds[_t] = list()
    if len(system) == 3:
        sections, systems = get_binary_systems(zpe_structures)
        binary_echs = list()
        for i, _system in enumerate(systems):
            if len(sections[i]) > 2:
                for _t in temp:
                    _ech = ExtendedConvexHull(sections[i], _system, t=_t)
                    # _ech.get_stable()
                    sections[i] = _ech.get_new_fitness()
                    _ech.plot(zpe_path=zpe_path, press=press)
                    binary_echs.append(_ech)
                save_dict = save_zpe(sections[i], zpe_path, temp, _system, press)
                for _t in temp:
                    if int(_t) == 0:
                        fit = save_dict['ZPE fitness']
                    else:
                        fit = save_dict[f'Fitness at T = {str(_t)}']
                    for i, _fit in enumerate(fit):
                        if _fit > 0:
                            antiseeds[f'{str(_t)}'].append(int(save_dict['EA'][i]))
    if '0' not in temp:
        temp = ['0'] + temp
    for _t in temp:
        ech = ExtendedConvexHull(zpe_structures, system, t=_t, antiseeds=antiseeds[_t])
        zpe_structures = ech.get_new_fitness(sections)
        ech.plot(zpe_path=zpe_path, plot=plot, press=press)
    return zpe_structures


def save_zpe(zpe_structures, zpe_path, temp, system, press):
    press = int(press)
    save_dict = dict()
    save_dict.update(
        {'EA': [structure['id'] for structure in zpe_structures]})
    save_dict.update({'Space Group': [structure['symmetry'][list(structure['symmetry'])[-1]][-1] for structure in zpe_structures]})
    save_dict.update({'Formula': [structure['formula'] for structure in zpe_structures]})
    if isinstance(zpe_structures[0]['convex hull x'], float):
        save_dict.update({'X': [structure['convex hull x'] for structure in zpe_structures]})
    else:
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
    name = os.path.join(zpe_path, f'{"-".join(system)}_convex_hull_{press}GPa.csv')
    df = pd.DataFrame.from_dict(save_dict)
    df.to_csv(name, index=False, float_format='%.8f')
    print(df)
    save_dict.update({'system': system})
    save_dict.update({'T': temp})
    save_dict.update({'structures': zpe_structures})
    return save_dict


