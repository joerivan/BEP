#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import numpy as np
import matplotlib.pyplot as plt
import os

#
# Plot the COHP for CO on a surface and compare it with
# CO in the gas phase
#
# note that in this script is assumed that the C-O interaction is the first interaction in the list of interactions (check lobsterin)
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
    fig, ax = plt.subplots(1, 1, dpi=300, figsize=(6,6))
    
    # CO gas phase total
    filename = os.path.join(os.path.dirname(__file__), 'Rhodium_CO', 'gas phase CO', 'COHPCAR.lobster')
    data, types = read_data(filename)
    metadata = np.loadtxt(filename, skiprows=1, max_rows=1)
    fermi = metadata[5]
    spinpol = int(metadata[1])-1
    nrints = int(metadata[0])

    energies = data[:,0]
    energies_dft_zero = energies+fermi
    if spinpol:
        cohp = data[:,3] + data[:,3+2*nrints]
        icohp = data[:,4] + data[:,4+2*nrints]
    else:
        cohp = data[:,3]
        icohp = data[:,4]
    index = np.where(energies == 0.0)
    icohp_fermi = icohp[index[0][0]]

    xlim = 40 
    ax.fill(cohp, energies_dft_zero, alpha=0.8, color=colors[4])
    ax.plot(cohp, energies_dft_zero, color=colors[4])
    ax.plot(icohp, energies_dft_zero, linestyle='--', linewidth=1.0, alpha=0.8, color=colors[4])
    ax.set_ylim(-35,5)
    ax.set_xlim(-xlim,xlim)
    ax.plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='E$_\mathrm{F}$')
    ax.text(-35, fermi+0.75, format(icohp_fermi, ".2f"), color=colors[4])
    ax.legend(loc='lower right')
    
    # # CO gas phase orbitalwise
    # pools = process_gas_data(data, types, fermi, spinpol, nrints)
    # xlim = 25
    # ax[1,0].plot(pools[0], energies_dft_zero, label='s-s', linewidth=1.0, color=colors[0])#'#555') '--', 
    # ax[1,0].plot(pools[1] + pools[2], energies_dft_zero, label='s-p', linewidth=1.0, color=colors[1])#'#AAA') '-.', 
    # ax[1,0].plot(pools[3], energies_dft_zero, label='p-p', linewidth=1.0, color=colors[3])#
    # # use only for verification:
    # # ax[1,0].plot(pools[0]+pools[1]+pools[2]+pools[3], energies_dft_zero, label='total', linewidth=1.0, color='black')
    # ax[1,0].plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='_nolegend_')
    # ax[1,0].set_ylim(-35,5)
    # ax[1,0].set_xlim(-xlim,xlim)
    # ax[1,0].legend(loc='lower right')

    # # Co(1121) / CO total
    # filename = os.path.join(os.path.dirname(__file__), 'data', 'co1121-co', 'COHPCAR.lobster')
    # data, types = read_data(filename)
    # metadata = np.loadtxt(filename, skiprows=1, max_rows=1)
    # fermi = metadata[5]
    # spinpol = int(metadata[1])-1
    # nrints = int(metadata[0])
    
    # energies = data[:,0]
    # energies_dft_zero = energies+fermi
    # if spinpol:
    #     cohp = data[:,3] + data[:,3+2*nrints]
    #     icohp = data[:,4] + data[:,4+2*nrints]
    # else:
    #     cohp = data[:,3]
    #     icohp = data[:,4]
    # index = np.where(energies == 0.0)
    # icohp_fermi = icohp[index[0][0]]

    # xlim = 28
    # ax[0,1].fill(cohp, energies_dft_zero, alpha=0.8, color=colors[4])
    # ax[0,1].plot(cohp, energies_dft_zero, color=colors[4])
    # ax[0,1].plot(icohp, energies_dft_zero, linestyle='--', linewidth=1.0, alpha=0.8, color=colors[4])
    # ax[0,1].set_ylim(-30,14.5)
    # ax[0,1].set_xlim(-xlim,xlim)
    # ax[0,1].plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='_nolegend_')
    # ax[0,1].text(-25, fermi+0.75, format(icohp_fermi, ".2f"), color=colors[4])
    
    # # Co(1121) / CO orbitalwise
    # pools = process_surface_data(data, types, 51, 49, fermi, spinpol, nrints)
    # xlim = 17
    # ax[1,1].plot(pools[0], energies_dft_zero, label='s-s', linewidth=1.0, color=colors[0])
    # ax[1,1].plot(pools[1] + pools[2], energies_dft_zero, label='s-p', linewidth=1.0, color=colors[1])
    # ax[1,1].plot(pools[3], energies_dft_zero, label='p-p', linewidth=1.0, color=colors[3])
    # # use only for verification:
    # # ax[1,1].plot(pools[0]+pools[1]+pools[2]+pools[3], energies_dft_zero, label='total', linewidth=1.0, color='black')
    # ax[1,1].set_ylim(-30,14.5)
    # ax[1,1].set_xlim(-xlim,xlim)
    # ax[1,1].plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='_nolegend_')
    # #ax[1,1].legend(loc='lower right')
    
    # # store value for next plot
    # co1121_pp = pools[3]
    # co1121_energies_dft_zero = energies_dft_zero
    
    # # Co(0001) / CO total
    # filename = os.path.join(os.path.dirname(__file__), 'data', 'co0001-co', 'COHPCAR.lobster')
    # data, types = read_data(filename)
    # metadata = np.loadtxt(filename, skiprows=1, max_rows=1)
    # fermi = metadata[5]
    # spinpol = int(metadata[1])-1
    # nrints = int(metadata[0])
    
    # energies = data[:,0]
    # energies_dft_zero = energies+fermi
    # if spinpol:
    #     cohp = data[:,3] + data[:,3+2*nrints]
    #     icohp = data[:,4] + data[:,4+2*nrints]
    # else:
    #     cohp = data[:,3]
    #     icohp = data[:,4]
    # index = np.where(energies == 0.0)
    # icohp_fermi = icohp[index[0][0]]

    # xlim = 42
    # ax[0,2].fill(cohp, energies_dft_zero, alpha=0.8, color=colors[4])
    # ax[0,2].plot(cohp, energies_dft_zero, color=colors[4])
    # ax[0,2].plot(icohp, energies_dft_zero, linestyle='--', linewidth=1.0, alpha=0.8, color=colors[4])  
    # ax[0,2].set_ylim(-30,15)
    # ax[0,2].set_xlim(-xlim,xlim)
    # ax[0,2].plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='_nolegend_')
    # ax[0,2].text(-25, fermi+0.75, format(icohp_fermi, ".2f"), color=colors[4])
    
    # pools = process_surface_data(data, types, 48, 46, fermi, spinpol, nrints)
    # xlim = 20
    # # ax[1,2].plot(pools[0], energies_dft_zero, label='_nolegend_', linewidth=1.0, color=colors[0])
    # # ax[1,2].plot(pools[1] + pools[2], energies_dft_zero, label='_nolegend_', linewidth=1.0, color=colors[1])
    # ax[1,2].plot(pools[3], energies_dft_zero, label='p-p$_{\mathrm{Co(0001)}}$', linewidth=1.0, color=colors[3])
    # ax[1,2].plot(co1121_pp, co1121_energies_dft_zero, label='p-p$_{\mathrm{Co(11\overline{2}1)}}$', linewidth=1.0, color='green')
    # # use only for verification:
    # # ax[1,2].plot(pools[0]+pools[1]+pools[2]+pools[3], energies_dft_zero, label='total', linewidth=1.0, color='black')
    # ax[1,2].set_ylim(-30,14.5)
    # ax[1,2].set_xlim(-xlim,xlim)
    # ax[1,2].plot(np.linspace(-xlim, xlim), np.linspace(fermi, fermi), color='darkorange', linewidth=0.5, label='_nolegend_')
    # ax[1,2].legend(loc='upper left')
    
    # for i in range(0,3):    
    #     for j in range(0,2):
    ax.grid(linestyle='--')
    ax.set_xlabel('COHP')
        
    # ax[0,0].set_ylabel('E$_\mathrm{DFT}$ [eV]')
    # ax[1,0].set_ylabel('E$_\mathrm{DFT}$ [eV]')
    # ax[0,0].set_title('C-O total')
    # ax[0,1].set_title('C-O / Co(11$\overline{2}$1) total')
    # ax[0,2].set_title('C-O / Co(0001) total')
    # ax[1,0].set_title('orbitalwise')
    # ax[1,1].set_title('orbitalwise')
    # ax[1,2].set_title('orbitalwise')
    # plt.tight_layout()
    
    # #filename2 = os.path.join(os.path.dirname(__file__), 'img', 'cohp_lobster_plot_surface.png')
    # #plt.savefig(filename2)

def process_gas_data(data, types, fermi, spinpol, nrints):
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
    return pools

def process_surface_data(data, types, at1, at2, fermi, spinpol, nrints):
    # process the data and produce a result
    energies = data[:,0]
    energies_dft_zero = energies+fermi
    orbitalpools = ['ss', 'sp', 'ps', 'pp', 'sd', 'ds', 'pd', 'dp', 'dd']
    pools = [np.zeros_like(energies) for j in range(0,len(orbitalpools))]
    total = np.zeros_like(energies)
    count = 0
    for i,datatype in enumerate(types):
        j = 1+(2*(i+1))
        k = 1+(2*(i+1))+(2*(nrints))
        if datatype['type'] == 'orbitalwise' and datatype['atomid1'] == at1 and datatype['atomid2'] == at2:
            intlabel = datatype['orbital1'][1] + datatype['orbital2'][1]
            poolid = orbitalpools.index(intlabel)
            if spinpol:
                pools[poolid] = pools[poolid] + data[:,j] + data[:,k]
            else:
                pools[poolid] = pools[poolid] + data[:,j]
            count = count + 1
            total = total + pools[poolid]
    return pools

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
        m = re.match(r'No.([0-9]+):([A-Za-z]{1,2})([0-9]+)\[([0-9][spdf]_?[a-z^2-]*)\]->([A-Za-z]{1,2})([0-9]+)\[([0-9][spdf]_?[a-z-^2]*)\].*', line)
        
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

