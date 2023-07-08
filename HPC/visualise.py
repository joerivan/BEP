import re
import os
import numpy as np
import matplotlib.pyplot as plt
import ase.io
from PIL import Image

def main():
    
    # Setting paths to files
    poscar = os.path.join(os.path.dirname(__file__),'data/POSCAR')
    params = os.path.join(os.path.dirname(__file__),'data/param.txt')
    
    # Settings for linear translation, 0 is equilibrium position
    param = np.loadtxt(params, max_rows=2)
    max_distance = param[0]
    steps = int(param[1])

    # Setting distance Rh to C in base structure
    base_dist = 1.389121355
      

    # Writing the files
    cnt = 0
    for d in np.linspace(0,max_distance,steps,2):
        cnt += 1

        tpath = os.path.join(os.path.dirname(__file__), 'output', '%i' % cnt)
        
        # Making the plots
        fig = plt.figure(dpi=600, figsize=[7,7])
        ax1 = plt.subplot(122)
        ax2 = plt.subplot(221)
        ax3 = plt.subplot(223)
        
        c_index = atom_index(poscar, 'C') 
        o_index = atom_index(poscar, 'O') 
        
        # Plot DOS
        e,sigma,pi = read_data_DOS(tpath, c_index, o_index)
        plot_DOS(ax1, e, sigma, pi)
        
        # Find and label peaks DOS
        # when peaks get to small remove the label from the labels list
        # and add the labels manually. See example in function plot_peaks
        
        idos_list = []
        peak_points_list = []
        
        thresholds = [0.9,1.2]
        
        count = 0
        for m in [sigma, pi]:
            idos_list.append(integrate_dos(m, e))
            peak_points_list.append(find_peaks(thresholds[count],m, e))
            count = count + 1
        
        labels = [['$3\sigma$','$4\sigma$','$5\sigma$','$6\sigma$','','','','','',''],
                  ['$1\pi$','$2\pi$','','','','','','']]
        
        dos = [sigma, pi]
        peak_x_list = []
        for n in range(2):
            peak_x, max_index = find_peaks_x(dos[n],idos_list[n],peak_points_list[n])
            if n == 1:
                for o in range(len(peak_x)):
                    peak_x[o] = peak_x[o] + dos[0][max_index[o]] 
            peak_x_list.append(peak_x)
            
        print(cnt)
        print(peak_x_list)

        
        for p in range(2):
            plot_peaks(idos_list[p], ax1, peak_points_list[p], e, 12, labels[p], peak_x_list[p])
     
        
        # Plot COHPs
        data, types, metadata = read_data_COHP(tpath)
        plot_COHP_total(ax2, data, types, metadata)
        plot_COHP_orbital(ax3, data, types, metadata)
        
        fig.tight_layout()

        # Saving just the plots
        plt.savefig(os.path.join(tpath, '%i.png' %cnt))
        
        # Finding distances from CONTCAR and adding to plot
        distance = np.loadtxt(os.path.join(tpath, 'param.txt')) + base_dist
        co_dist = find_bondlength(tpath)
        
        # Making images of the distances
        distance = format(distance, '.2f')
        co_dist = format(co_dist, '.2f')
        latex_image(r'|\vec{r}_{Rh-C}|=',distance,'Rh-C')   
        latex_image(r'|\vec{r}_{C-O}|=',co_dist,'C-O') 
        
        # Adding distances to image
        plot_image = Image.open(os.path.join(tpath, '%i.png' %cnt))
        img = add_distances(plot_image)
        
        # Saving image and closing off
        img.save(os.path.join(os.path.dirname(__file__),'output','images',
                              '%i.png' %cnt))
        plt.close('all')
        
        
def atom_index(poscar,atom_name):
    """
    Finds the index of an atom based on the poscar
    
    poscar : POSCAR file
    atom_name : string with a atom to find
    """

    atoms = np.loadtxt(poscar,dtype=str, skiprows=5, max_rows=1)
    nr_atoms = np.loadtxt(poscar, skiprows=6, max_rows=1)
    
    i_count = 0
    
    for i in atoms:
        if i == atom_name:
            ind_atom = i_count
        else:
            i_count += 1
    
    index_atom = 0
    
    for j in range(0,ind_atom):
        index_atom = index_atom + nr_atoms[j]
        
    
    return int(index_atom)

def read_data_DOS(folder,c_index,o_index):
    """
    Returns the extracted data of DOS
    
    folder : folder which contains a DOSCAR.lobster file
    """
    skip_c = 908 + c_index*902
    skip_o = 908 + o_index*902
    data_c = np.loadtxt(os.path.join(folder, 'DOSCAR.lobster'),
                        skiprows=skip_c, max_rows=901)
    data_o = np.loadtxt(os.path.join(folder, 'DOSCAR.lobster'),
                        skiprows=skip_o, max_rows=901)
    
    e = data_c[:,0]
    s = data_c[:,1] + data_o[:,1]
    py = data_c[:,2] + data_o[:,2]
    pz = data_c[:,3] + data_o[:,3]
    px = data_c[:,4] + data_o[:,4]
    
    sigma = s + pz
    pi = px + py
    
    return e,sigma,pi


def read_data_COHP(folder):
    """
    Returns extracted data of COHP
    
    folder : folder containing COHP.lobster file
    """
    filename = os.path.join(folder,'COHPCAR.lobster')
    
    f = open(filename)
    f.readline()                          # skip first line
    nrints = int(f.readline().split()[0]) # read nr of interactions
    f.readline()                          # skip another line

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
    
    metadata = np.loadtxt(filename, skiprows=1, max_rows=1)

    # read all the data
    data = np.loadtxt(filename, skiprows=nrints+2)
    
    return data, types, metadata

def plot_DOS(ax,e,sigma,pi, title=None, 
             colors=['#fe6100','#648fff']):
    """
    Plots the DOS
    """
    labels = ['$\sigma$','$\pi$']
    bottom = np.zeros_like(sigma)
    top = np.zeros_like(sigma)
    
    for l,c,label in zip([sigma,pi], colors, labels):
        top = np.add(top,l)
        ax.fill_betweenx(e, bottom, top, color=c, label=label, alpha=0.8)
        ax.plot(top, e, linewidth=0.9, color=c)
        bottom = np.add(bottom,l)
    
    ax.grid(linestyle='--', alpha=0.5, zorder=-1)
    ax.set_ylim(-30,15)
    ax.set_xlim(0,12)
    ax.legend(loc='lower right')
    ax.set_xlabel('pDOS (-)')
    ax.set_ylabel('Energy E - E$_{f}$ [eV]')
    
    if title:
        ax.set_title(title)

def find_peaks(threshold, dos, energies_dft_zero):
    """
    finds the peaks, at what y (energy) they are
    """
    # find peak areas
    dos = abs(dos) 

    # find start and end points of peaks
    dos_bin = [1 if x>=threshold else 0 for x in dos]
    peak_points = []
    for i in range(len(dos_bin)-1):
        if dos_bin[i] == 0 and dos_bin[i+1] == 1:
            peak_points.append(i+1)
        if dos_bin[i] == 1 and dos_bin[i+1] == 0:
            peak_points.append(i)
    return peak_points

def find_peaks_x(dos,idos,peak_points):
    """
    find at what value of x (pDOS) the peaks are
    """
    count = 0
    peak_x = []
    max_index = []
    for i in range(int(len(peak_points)/2.0)):
        index1 = peak_points[count]
        index2 = peak_points[count+1]
        delta_peak = idos[index2]-idos[index1]
        if abs(delta_peak) > 0.00001:
            label_loc = max(dos[index1:index2])
            peak_x.append(label_loc+0.1)
            max_index.append(np.where(dos == label_loc)[0][0])

        count = count + 2
    return peak_x, max_index
    
def plot_peaks(idos, ax, peak_points, energies_dft_zero, xlim, labels, peak_x):
    """
    plots the peak labels
    """
    lim_plot = 11
    count = 0
    for i in range(int(len(peak_points)/2.0)):
        index1 = peak_points[count]
        index2 = peak_points[count+1]
        delta_peak = idos[index2]-idos[index1]
        
        if int(len(peak_points)/2.0) > len(peak_x):
            peak_x.append(0)
        
        if abs(delta_peak) > 0.00001:
            if peak_x[i] < lim_plot:
              height = energies_dft_zero[index1]+abs((energies_dft_zero[index1]-energies_dft_zero[index2])/2.0)-0.5
              ax.text(peak_x[i], height, labels[i], color='black', size=11)
              print(height)
            elif peak_x[i] >= lim_plot:
              peak_x[i] = lim_plot
              height = energies_dft_zero[index1]+abs((energies_dft_zero[index1]-energies_dft_zero[index2])/2.0)-0.5
              ax.text(peak_x[i], height, labels[i], color='black', size=11)
              print(height)
        count = count + 2

        # Add label for 6sigma manually
        #ax.text(1.8, 5.5, '$6\sigma$', color='black', size=11)
        
    
def integrate_dos(dos, energies_dft_zero):
    """
    integrated the dos
    """
    idos = np.zeros_like(energies_dft_zero)
    dx = (abs(energies_dft_zero[0]) + abs(energies_dft_zero[-1])) / float(len(energies_dft_zero))
    for i in range(len(dos)):
        idos[i] = idos[i-1] + dos[i] * dx  
    return idos

def plot_COHP_total(ax, data, types, metadata, title=None,
              colors=['#785EF0','#000000']):
    """
    Plots total COHP
    """    
    spinpol = int(metadata[1])-1
    nrints = int(metadata[0])
    
    energies = data[:,0]
    energies_dft_zero = energies
    if spinpol:
        cohp = data[:,3] + data[:,3+2*nrints]
        icohp = data[:,4] + data[:,4+2*nrints]
    else:
        cohp = data[:,3]
        icohp = data[:,4]

    xlim = 40 
    ax.fill(cohp, energies_dft_zero, alpha=0.8, color=colors[0], label='COHP')
    ax.plot(cohp, energies_dft_zero, color=colors[0], linewidth=0.7)
    ax.plot(icohp, energies_dft_zero, linestyle='--', linewidth=1.0, alpha=0.8,
            color=colors[1], label='iCOHP')
    ax.set_ylim(-30,10)
    ax.set_xlim(-xlim,xlim)
    ax.legend(loc='lower right')
    
    ax.grid(linestyle='--', alpha=0.5, zorder=-1)
    ax.set_xlabel('COHP (-)')
    ax.set_ylabel('Energy E - E$_{f}$ [eV]')
    
    if title:
        ax.set_title(title)
        
def plot_COHP_orbital(ax, data, types, metadata, title=None,
              colors=['#fe6100','#648fff']):
    """
    Plots COHP orbitals wise
    """    
    spinpol = int(metadata[1])-1
    nrints = int(metadata[0])
    
    energies = data[:,0]
    energies_dft_zero = energies
    pools = [np.zeros_like(energies) for j in range(0,16)]
    orbitalpools = ['ss', 'sp_z', 'p_zs', 'p_zp_z', 'p_xp_x', 'p_xp_y', 'p_yp_x', 'p_yp_y',
                    'sp_x', 'p_xs', 'sp_y', 'p_ys', 'p_xp_z', 'p_zp_x', 'p_yp_z', 'p_zp_y']
    for i,datatype in enumerate(types):
        j = 1+(2*(i+1))
        k = 1+(2*(i+1))+(2*(nrints))
        if datatype['type'] == 'orbitalwise':
            intlabel = datatype['orbital1'][1:] + datatype['orbital2'][1:]
            poolid = orbitalpools.index(intlabel)
            if spinpol:
                pools[poolid] = pools[poolid] + data[:,j] + data[:,k]
            else:
                pools[poolid] = pools[poolid] + data[:,j]
    
    xlim = 40 
    
    sigmas = pools[0] + pools[1] + pools[2] + pools[3]
    pis = pools[4] + pools[5] + pools[6] + pools[7]
    
    labels = ['$\sigma$','$\pi$']
    
    ax.fill(sigmas, energies_dft_zero, alpha=0.6, color=colors[0],
            label=labels[0], zorder=4)
    ax.plot(sigmas, energies_dft_zero, color=colors[0], zorder=2,
            linewidth=0.7)
    ax.fill(pis, energies_dft_zero, alpha=0.6, color=colors[1],
            label=labels[1], zorder=3)
    ax.plot(pis, energies_dft_zero, color=colors[1], zorder=1,
            linewidth=0.7)

    ax.set_ylim(-30,10)
    ax.set_xlim(-xlim,xlim)
    ax.legend(loc='lower right')
    
    ax.grid(linestyle='--', alpha=0.5, zorder=-1)
    ax.set_xlabel('COHP (-)')
    ax.set_ylabel('Energy E - E$_{f}$ [eV]')
    
    if title:
        ax.set_title(title)

def find_bondlength(path):
    """
    Finds the bondlength between C and O for a CONTCAR in a given path
    """
    contcar = os.path.join(path,"CONTCAR")
    
    c_index = atom_index(contcar, 'C') 
    o_index = atom_index(contcar, 'O') 
    
    struc = ase.io.read(contcar)
    bond_dist = struc.positions[o_index][2]-struc.positions[c_index][2]
    
    return bond_dist
    
def latex_image(tex, value, name):
    """ 
    Generates a latex image with matplotlib and saves it 
    """
    plt.figure(figsize=(8,8))
    plt.axis('off')
    plt.text(0.05, 0.35, f'${tex}$ {value}', size=250)
    
    main_dir = os.path.dirname(__file__)
    temp_dir = os.path.join(main_dir, "temp")
    
    if os.path.isdir(temp_dir) is False:
        os.mkdir(temp_dir)

    plt.savefig(os.path.join(temp_dir, 'tex_%s.png' %name), bbox_inches = 'tight')
    plt.close()
  
def add_distances(img):
    """
    Adds a white pace next to plot where the distance from surface and bond
    length are shown
    """
    size = img.size
    
    base_img = Image.new(mode="RGB", size = (size[0]+2100,size[1]),
                          color=(255,255,255))
    
    base_img.paste(img, (0,0), mask=img)
    
    height_text = 420
    
    Rh_C = Image.open(os.path.join(os.path.dirname(__file__), 'temp/tex_Rh-C.png')).convert("RGBA")
    width, height = Rh_C.size
    ratio = width/height
    new_height = height_text
    new_width = int(ratio*new_height) 
    Rh_C = Rh_C.resize((new_width,new_height))
    
    C_O = Image.open(os.path.join(os.path.dirname(__file__), 'temp/tex_C-O.png')).convert("RGBA")
    width, height = C_O.size
    ratio = width/height
    new_height = height_text
    new_width = int(ratio*new_height) 
    C_O = C_O.resize((new_width,new_height))

    base_img.paste(Rh_C, (4200,900), mask=Rh_C)
    base_img.paste(C_O, (4200,1400), mask=C_O)
    
    return base_img
            
  
if __name__ == '__main__':
    main()
