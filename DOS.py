#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

def main():
    fig, ax = plt.subplots(1, 1, dpi=300, figsize=(3,6))
    
    folder = os.path.join(os.path.dirname(__file__), 'Rhodium_CO', 'gas phase CO')
    e, s, px, py, pz = read_data(folder)
    plot_dos(ax, e, s, px, py, pz, 'CO DOS')

    plt.tight_layout()

def plot_dos(ax, e, s, px, py, pz, title=None):
    colors = ['#2a9d8f', '#e9c46a', '#f4a261', '#e76f51']
    labels = ['s', 'px', 'py', 'pz']
    base = np.zeros_like(s)
    top = np.zeros_like(s)
    for l,c,label in zip([s, px, py, pz], colors, labels):
        top = np.add(top, l)
        ax.fill_betweenx(e, base, top, color=c, label=label, alpha=0.8)
        ax.plot(top, e, linewidth=0.5, color=c)
        base = np.add(base, l)

    ax.grid(linestyle='--', alpha=0.5, zorder=-1)
    ax.set_ylim(-30,15)
    ax.set_xlim(0, 12)
    ax.legend(loc='lower right')
    ax.set_xlabel('Density of States')
    ax.set_ylabel('Energy E - E$_{f}$ [eV]')
    
    if title:
        ax.set_title(title)

def read_data(folder):
    data_c = np.loadtxt(os.path.join(folder, 'DOSCAR.lobster'),
                        skiprows = 908, max_rows=901)

    data_o = np.loadtxt(os.path.join(folder, 'DOSCAR.lobster'),
                        skiprows = 1810, max_rows=901)

    e = data_c[:,0]
    s = data_c[:,1] + data_o[:,1]
    py = data_c[:,2] + data_o[:,2]
    pz = data_c[:,3] + data_o[:,3]
    px = data_c[:,4] + data_o[:,4]
    
    return e,s,px,py,pz

if __name__ == '__main__':
    main()