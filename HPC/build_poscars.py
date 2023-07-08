import ase.io
import os
import numpy as np
from ase.constraints import FixAtoms, FixedLine, FixedPlane

def main():
    
    # Setting paths to files
    poscar = os.path.join(os.path.dirname(__file__),'data/POSCAR')
    params = os.path.join(os.path.dirname(__file__),'data/param.txt')
    
    # Settings for linear translation, 0 is equilibrium position
    param = np.loadtxt(params, max_rows=2)
    max_distance = param[0]	# Take first value from param file for the max distance
    steps = int(param[1])	# Take second value from param file for numder of steps
    
    # Writing the files
    cnt = 0
    for d in np.linspace(0, max_distance, steps):
        cnt += 1
        
        CO_on_Rh = change_distance(poscar, d)
        
        tpath = os.path.join(os.path.dirname(__file__), 'output', '%i' % cnt)
        if not os.path.exists(tpath):
            os.mkdir(tpath)
        ase.io.vasp.write_vasp(os.path.join(tpath, 'POSCAR'), CO_on_Rh,
                               direct=True)
        
        
        param_path = os.path.join(os.path.dirname(__file__),'output/%i/param.txt' % cnt)
        f = open(param_path, 'w')
        f.write('%f\n' % d)
        f.close()
    
    
def change_distance(poscar,dist):
    """
    Changes the distance of a CO molecule to a surface only in z direction

    poscar : POSCAR of CO on a surface
    dist : increase in distance
    """
    
    ref_struc = ase.io.read(poscar)

    index_C = atom_index(poscar, 'C')
    index_O = atom_index(poscar, 'O')

    new_struc = ref_struc.copy()
    
    new_struc.positions[index_C][2] = new_struc.positions[index_C][2] + dist
    new_struc.positions[index_O][2] = new_struc.positions[index_O][2] + dist

    constraint_c = FixAtoms(
        [atom.index for atom in new_struc if atom.symbol == 'C']
        )
    constraint_o = FixedLine(
        [atom.index for atom in new_struc if atom.symbol == 'O'],
        direction=[0,0,1]
        )
    constraint_rh = FixedPlane(
        [atom.index for atom in new_struc if atom.symbol == 'Rh'],
        direction=[0,0,1]
        )
    
    new_struc.set_constraint([constraint_c,constraint_o,constraint_rh])
    
    return new_struc


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
    
    
    
if __name__ == '__main__':
    main()
