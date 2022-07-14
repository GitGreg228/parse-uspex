import re
import os
import json
import yaml
import numpy as np
import pandas as pd

import argparse
from functools import reduce
from tqdm import tqdm
from math import gcd
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.cif import CifWriter
from pymatgen.core.structure import IStructure, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from adjustText import adjust_text

from utils import get_comp_cat, formula_from_comp, listdirs, get_comp


def get_triangle(save_dict, compositions):
    x1, x2 = save_dict['X1'], save_dict['X2']
    minx1, maxx1, minx2, maxx2 = 0.5, -0.5, 1, -0.5
    for i in range(len(x1)):
        comp = compositions[i]
        if get_comp_cat(comp) != 'single':
            if x1[i] > maxx1:
                maxx1 = x1[i]
            if x1[i] < minx1:
                minx1 = x1[i]
            if x2[i] > maxx2:
                maxx2 = x2[i]
            if x2[i] < minx1:
                minx2 = x2[i]
    minx1, maxx1, minx2, maxx2 = minx1 - 0.05, maxx1 + 0.1, minx2 - 0.2, maxx2 + 0.05
    return minx1, maxx1, minx2, maxx2


def getprop(save_dict, _t):
    x1, x2 = save_dict['X1'], save_dict['X2']
    system = save_dict['system']
    structures = [IStructure.from_dict(structure['structure']) for structure in save_dict['structures']]
    compositions = [get_comp(structure, system)[1] for structure in structures]
    formulas = [formula_from_comp(composition, system) for composition in compositions]
    if int(_t) == 0:
        fit = save_dict['ZPE fitness']
    else:
        fit = save_dict[f'Fitness at T = {str(_t)}']


def plot_ch_2d(save_dict):
    temp = save_dict['T']
    x1, x2 = save_dict['X1'], save_dict['X2']
    system = save_dict['system']
    structures = [IStructure.from_dict(structure['structure']) for structure in save_dict['structures']]
    compositions = [get_comp(structure, system)[1] for structure in structures]
    formulas = [formula_from_comp(composition, system) for composition in compositions]
    minx1, maxx1, minx2, maxx2 = get_triangle(save_dict, compositions)
    for i, _t in enumerate(temp):
        plt.figure(figsize=(9, 9))
        if int(temp[i]) == 0:
            y = save_dict['ZPE Y']
            fit = save_dict['ZPE fitness']
        else:
            y = save_dict[f'Y at T = {str(_t)}']
            fit = save_dict[f'Fitness at T = {str(_t)}']
        stable, unstable, unstable_fit = list(), list(), list()
        for i, f in enumerate(fit):
            if f == 0:
                stable.append([x1[i], x2[i]])
            else:
                unstable.append([x1[i], x2[i]])
                unstable_fit.append(fit[i])
        stable = np.asarray(stable).transpose()
        unstable = np.asarray(unstable).transpose()
        plt.scatter(stable[0], stable[1], c='k')
        plt.scatter(unstable[0], unstable[1], c=unstable_fit, marker='D')
        plt.axis('off')
        plt.xlim(minx1, maxx1)
        plt.ylim(minx2, maxx2)
        fs = [10 if f > 0 else 15 for f in fit]
        texts = [plt.text(x1[i], x2[i], formulas[i], ha='center', va='center', fontsize=fs[i]) for i in range(len(x1))]
        adjust_text(texts)
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    args = parser.parse_args()

    full_dirs = list()
    for dir in listdirs(args.path):
        if 'extended_convex_hull' in os.listdir(dir) and 'extended_convex_hull_POSCARS' in os.listdir(dir):
            full_dirs.append(dir)
    output = ', '.join(full_dirs)
    print(f'Will work in {output}')
    for dir in full_dirs:
        print(f'Working in {dir} ...')
        zpe_path = os.path.join(dir, 'ZPE')
        with open(os.path.join(zpe_path, 'save_dict.json'), 'r') as f:
            save_dict = json.load(f)
        plot_ch_2d(save_dict)


if __name__ == '__main__':
    main()
