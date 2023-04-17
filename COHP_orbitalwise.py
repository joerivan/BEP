#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import numpy as np
import matplotlib.pyplot as plt
import os

#
# Plot the lm-decomposed COHP for the CO molecule in the gas phase
#
# COHPCAR is built up as follows:
# list of interactions (including average (Average), total (e.g. No.1:O2->C1) and orbitalwise (e.g. No.1:O2[2s]->C1[2s]) interactions)
# col0    col1           col2           col3           col4          ... col(2*nints+1)   col(2*nints+2)  ... col(2*nints+2*nints)
# energy  pcohp_int1_up  icohp_int1_up  pcohp_int2_up  icohp_int2_up ... pcohp_int1_down  icohp_int1_down ... icohp_nints_down
#
# both spin polarized COHPCARs and non-spin polarized COHPCARs can be read into the script

# script parameters
colors = ['#2a9d8f', # s-color
          '#e9c46a', '#f4a261', '#e76f51', #p-colors
          '#264653'] # total color

def main():
    filename = os.path.join(os.path.dirname(__file__), 'Rhodium_CO', 'gas phase CO', 'COHPCAR.lobster')
    data, types = read_data(filename)
    metadata = np.loadtxt(filename, skiprows=1, max_rows=1)
    fermi = metadata[5]
    spinpol = int(metadata[1])-1
    nrints = int(metadata[0])
    
    # process the data and produce a result
    energies = data[:,0]
    energies_dft_zero = energies+fermi
    pools = [np.zeros_like(energies) for j in range(0,4)]
    orbitalpools = ['ss', 'sp', 'ps', 'pp']
    for i,datatype in enumerate(types):
        j = 1+(2*(i+1))
        k = 1+(2*(i+1))+(2*(nrints))
        if datatype['type'] == 'orbitalwise':
            intlabel = datatype['orbital1'][1] + datatype['orbital2'][1]
            poolid = orbitalpools.index(intlabel)
            if spinpol:
                pools[poolid] = pools[poolid] + data[:,j] + data[:,k]
            else:
                pools[poolid] = pools[poolid] + data[:,j]
    if spinpol:
        cohp = data[:,3] + data[:,3+2*nrints]
        icohp = data[:,4] + data[:,4+2*nrints]
    else:
        cohp = data[:,3]
        icohp = data[:,4]
    index = np.where(energies == 0.0)
    icohp_fermi = icohp[index[0][0]]
    
    # plot results
    fig, ax = plt.subplots(1, 4, dpi=300, sharey=True, figsize=(8,4))
    ax[0].fill(pools[0], energies_dft_zero, alpha=0.8, color=colors[0])
    ax[0].plot(pools[0], energies_dft_zero, color=colors[0])
    ax[1].fill(pools[1] + pools[2], energies_dft_zero, alpha=0.8, color=colors[1])
    ax[1].plot(pools[1] + pools[2], energies_dft_zero, color=colors[1])
    ax[2].fill(pools[3], energies_dft_zero, alpha=0.8, color=colors[3])
    ax[2].plot(pools[3], energies_dft_zero, color=colors[3])
    ax[3].fill(cohp, energies_dft_zero, alpha=0.8, color=colors[4])
    ax[3].plot(cohp, energies_dft_zero, color=colors[4])
    ax[3].plot(icohp, energies_dft_zero, linestyle='--', linewidth=1.0, alpha=0.8, color=colors[4])
    
    xlim1 = 25
    for i in range(0,4):    
        ax[i].grid(linestyle='--')
        ax[i].set_xlabel('COHP')
        ax[i].set_xlim(-xlim1,xlim1)
        ax[i].set_ylim(-35,5)
        ax[i].plot(np.linspace(-xlim1, xlim1), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='E$_\mathrm{F}$')
    
    ax[0].set_ylabel('E$_\mathrm{DFT}$ [eV]')
    ax[0].legend(loc='upper left')
    ax[0].set_title('(s-s)')
    ax[1].set_title('(s-p)')
    ax[2].set_title('(p-p)')
    ax[3].set_title('total')
    xlim2 = 40
    ax[3].set_xlim(-xlim2,xlim2)
    ax[3].plot(np.linspace(-xlim2, xlim2), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='E$_\mathrm{F}$')
    ax[3].text(-35, fermi+0.25, format(icohp_fermi, ".2f"), color=colors[4])
    plt.tight_layout()

def read_data(filename): 
    f = open(filename)
    f.readline()                          # skip first line
    nrints = int(f.readline().split()[0]) # read nr of interactions
    f.readline()                          # skip another line

    print('Exploring %i interactions' % nrints)

    # loop over interactions and collect the types
    types = []
    for i in range(1, nrints):
        line = f.readline()
        m = re.match(r'No.([0-9]+):([A-Za-z]{1,2})([0-9]+)\[([0-9][spdf]_?[a-z^2]*)\]->([A-Za-z]{1,2})([0-9]+)\[([0-9][spdf]_?[a-z-^2]*)\].*', line)
        
        if m:
            types.append({
                'type' : 'orbitalwise',
                'interaction_id' : int(m.group(1)),
                'element1' : m.group(2),
                'atomid1' : int(m.group(3)),
                'orbital1' : m.group(4),
                'element2' : m.group(5),
                'atomid2' : int(m.group(6)),
                'orbital2' : m.group(7),
            })
        else:
            m = re.match(r'No.([0-9]+):([A-Za-z]{1,2})([0-9]+)->([A-Za-z]{1,2})([0-9]+).*', line)
            if m:
                types.append({
                    'type' : 'total',
                    'interaction_id' : int(m.group(1)),
                    'element1' : m.group(2),
                    'atomid1' : int(m.group(3)),
                    'element2' : m.group(4),
                    'atomid2' : int(m.group(5)),
                })
            else:
                raise Exception('Cannot parse line: %s' % line)
    f.close()

    # read all the data
    data = np.loadtxt(filename, skiprows=nrints+2)
    
    return data, types

if __name__ == '__main__':
    main()
