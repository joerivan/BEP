import bpy
import numpy as np
import os
import time

atom_radii = {
    'H': 0.2,
    'C': 0.5,
    'O': 0.5,
    'Rh': 1.2
}

atom_colors = {
    'H': 'ffffff', 
    'C': '101010',
    'O': 'DD0000',
    'Rh': '0A7D8C'
}

settings = {
        'resolution': 512,
        'camera_location': (-10,0,2.5),
        'camera_rotation': (85/90*np.pi/2,0,-np.pi/2),
        'camera_scale' : 7.5
    }

def main():
    # set the scene
    set_environment(settings)
    
    # clear any remaining objects
    prune_scene()
    
    # read molecule file and load it
    root = "rootfolder"
    
    for i in range(0,1):
        mol = read_contcar(os.path.join(root, 'positions/%i/CONTCAR' %(i+1)))
        create_atoms(mol)
        create_bonds(mol)   
        render_scene(os.path.join(root, 'images/%i.png' % (i+1)))
        prune_scene()
        
        
        
def prune_scene():
    for obj in bpy.data.objects:
        if obj.name.startswith('atom') or \
           obj.name.startswith('bond') or \
           obj.name.startswith('isosurface'):
            bpy.ops.object.select_all(action='DESELECT')
            obj.select_set(True)
            bpy.ops.object.delete(use_global=False)

    # finally delete all materials
    for material in bpy.data.materials:
        material.user_clear()
        bpy.data.materials.remove(material)

def create_atoms(mol):
    """
    Create atoms
    """
    for i,at in enumerate(mol):
        scale = atom_radii[at[0]]
        bpy.ops.surface.primitive_nurbs_surface_sphere_add(
            radius=scale, 
            enter_editmode=False, 
            align='WORLD', 
            location=at[1])
        obj = bpy.context.view_layer.objects.active
        obj.name = "atom-%s-%03i" % (at[0],i)
        bpy.ops.object.shade_smooth()
        
        # set a material
        mat = create_material(at[0], atom_colors[at[0]])
        print(mat)
        obj.data.materials.append(mat)

def create_bonds(mol):
    """
    Create bonds between atoms
    """
    # set default orientation of bonds (fixed!)
    z = np.array([0,0,1])
    
    # add new bonds material if it does not yet exist
    matbond = create_material('bond', '555555')
    
    for i,at1 in enumerate(mol):
        r1 = np.array(at1[1])
        for j,at2 in enumerate(mol[i+1:]):
            r2 = np.array(at2[1])
            dist = np.linalg.norm(r2 - r1)
            
            # only create a bond if the distance is less than 1.5 A
            if dist < 2.3:
                axis = np.cross(z,r2-r1)
                axis /= np.linalg.norm(axis)
                angle = np.arccos(np.dot(r2-r1,z)/dist)
                
                bpy.ops.surface.primitive_nurbs_surface_cylinder_add(
                    enter_editmode=False, 
                    align='WORLD',
                    location=tuple((r1 + r2) / 2)
                )
                
                obj = bpy.context.view_layer.objects.active
                obj.scale = (0.2, 0.2, dist/2)
                obj.rotation_mode = 'AXIS_ANGLE'
                obj.rotation_axis_angle = (angle, axis[0], axis[1], axis[2])
                
                obj.name = "bond-%s-%03i-%s-%03i" % (at1[0],i,at2[0],j)
                bpy.ops.object.shade_smooth()
                obj.data.materials.append(matbond)

def set_environment(settings):
    """
    Specify canvas size, remove default objects, reset positions of
    camera and light, define film and set background
    """
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.device = 'GPU'
    bpy.context.scene.render.resolution_x = settings['resolution']
    bpy.context.scene.render.resolution_y = settings['resolution']
    #bpy.context.scene.render.tile_x = settings['resolution']
    #bpy.context.scene.render.tile_y = settings['resolution']
    
    # remove cube
    if 'Cube' in bpy.data.objects:
        o = bpy.data.objects['Cube']
        bpy.data.objects.remove(o, do_unlink=True)
    
    # set camera into default position
    bpy.data.objects['Camera'].location = tuple(settings['camera_location'])
    bpy.data.objects['Camera'].rotation_euler = tuple(settings['camera_rotation'])
    bpy.data.objects['Camera'].data.clip_end = 1000
    #bpy.data.objects['Camera'].data.type = 'ORTHO'
    bpy.data.objects['Camera'].data.ortho_scale = settings['camera_scale']
    
    # set lights
    if 'Light' not in bpy.data.objects:
        bpy.ops.object.light_add(type='POINT', radius=1, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.context.object.name = 'Light'
        bpy.data.objects['Light'].data.type = 'AREA'
        bpy.data.objects['Light'].data.energy = 1e4
        bpy.data.objects['Light'].location = (-10,10,10)
        bpy.data.objects['Light'].rotation_euler = tuple(np.radians([55, 0, 225]))
        bpy.data.objects['Light'].data.shape = 'DISK'
        bpy.data.objects['Light'].data.size = 10

    # set film
    bpy.context.scene.render.film_transparent = True
    
    # set background
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[1].default_value = 1

def create_material(name, color, alpha=1.0):
    """
    Build a new material
    """
    # early exit if material already exists
    if name in bpy.data.materials:
        return bpy.data.materials[name]
    
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    
    # set base color
    mat.node_tree.nodes["Principled BSDF"].inputs[0].default_value = hex2rgbtuple(color)
    
    # subsurface modifier
    mat.node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0.2
    
    # set subsurface radii
    mat.node_tree.nodes["Principled BSDF"].inputs[2].default_value = (0.3,0.3,0.3)

    # set subsurface color
    mat.node_tree.nodes["Principled BSDF"].inputs[3].default_value = hex2rgbtuple('000000')

    # metallic
    if name=='Rh':
        mat.node_tree.nodes["Principled BSDF"].inputs[6].default_value = 0.3
    else:
        mat.node_tree.nodes["Principled BSDF"].inputs[6].default_value = 0

    # roughness
    if name=='Rh':
        mat.node_tree.nodes["Principled BSDF"].inputs[9].default_value = 0.0
    else:
        mat.node_tree.nodes["Principled BSDF"].inputs[9].default_value = 0.0
    
    # alpha
    mat.node_tree.nodes["Principled BSDF"].inputs[21].default_value = alpha

    return mat

def render_scene(outputfile, samples=512):
    """
    Render the scene
    """
    bpy.context.scene.cycles.samples = samples
    
    print('Start render')
    start = time.time()
    bpy.data.scenes['Scene'].render.filepath = outputfile
    bpy.ops.render.render(write_still=True)
    end = time.time()
    print('Finished rendering frame in %.1f seconds' % (end - start))

def read_contcar(filename):
    f = open(filename)
    
    # skip first line, name of system
    f.readline()
    
    # extract scaling factor
    scaling_factor = float(f.readline())
    
    # find the unitcell matrix
    matrix = np.zeros((3,3))
    for i in range(0,3):
        row = f.readline().split()
        for j in range(0,3):
            matrix[i,j] = float(row[j])
    matrix = scaling_factor * matrix
    
    # find the atoms and thei occurances
    atoms = f.readline().split()
    atoms_occurance = f.readline().split()
    
    # skip selective dynamics and direct coordinates line
    f.readline()
    f.readline()
    
    mol = []
    for k in range(0,len(atoms)):
        for l in range(0,int(atoms_occurance[k])):
            # get fractional coordinates
            c1c2c3 = f.readline().split()
            c1 = float(c1c2c3[0])
            c2 = float(c1c2c3[1])
            c3 = float(c1c2c3[2])
            c1c2c3 = [c1,c2,c3]
            
            xyz = matrix.transpose() @ c1c2c3
            s_xyz = shift_struc(xyz)
            
            if s_xyz[2] > -6:
                append_mol(mol,atoms[k],s_xyz) 
            
            sz = 2
            if atoms[k] == "Rh":
                for x in np.arange(-sz,sz+1):
                    for y in np.arange(-sz,sz+1):
                        if x == 0 and y == 0:
                            print(' ')
                        else:
                            dup_c1c2c3 = c1c2c3.copy()
                            dup_c1c2c3[0] = dup_c1c2c3[0] + x
                            dup_c1c2c3[1] = dup_c1c2c3[1] + y
                            xyz = np.matmul(matrix.transpose(),dup_c1c2c3)
                            s_xyz = shift_struc(xyz)
                
                            if s_xyz[2] > -6:
                                append_mol(mol,atoms[k],s_xyz) 
                            
    f.close()
    
    return mol

def append_mol(mol,atom,xyz):
    
    mol.append(
                [atom,
                (
                    xyz[0],
                    xyz[1],
                    xyz[2]
                )]
            )

def shift_struc(xyz):
    """
    Places the structure on a predefined location
    uses the zero position of the carbon adsorbed
    """
    x_set = 1.3519881656286596
    y_set = 0.7805707313668949
    z_set = 29.652349912212305
    
    xyz[0] = xyz[0]-x_set
    xyz[1] = xyz[1]-y_set
    xyz[2] = xyz[2]-z_set
    
    return xyz


def hex2rgbtuple(hexcode):
    """
    Convert 6-digit color hexcode to a tuple of floats
    """
    hexcode += "FF"
    hextuple = tuple([int(hexcode[i:i+2], 16)/255.0 for i in [0,2,4,6]])
    
    return tuple([color_srgb_to_scene_linear(c) for c in hextuple])

def color_srgb_to_scene_linear(c):
    """
    Convert RGB to sRGB
    """
    if c < 0.04045:
        return 0.0 if c < 0.0 else c * (1.0 / 12.92)
    else:
        return ((c + 0.055) * (1.0 / 1.055)) ** 2.4

if __name__ == '__main__':
    main()