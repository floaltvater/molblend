# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

### TODO
# - need to clean up file reading (e.g., that everything returns the same thing
#   in the same format).

if "bpy" in locals():
    import imp
    imp.reload(mb_utils)
    imp.reload(mb_io_files)
else:
    from molblend import mb_utils
    from molblend import mb_io_files

import time

import bpy
import bmesh

from molblend.elements_default import ELEMENTS as ELEMENTS_DEFAULT
from molblend.mb_helper import debug_print
from molblend.mb_io_files import *

#--- Read file functions -----------------------------------------------------#
def is_inside_of_planes(planes, l0, flip=False):
    #raise IOError("mask not correctly implemented")
    for n, p0 in planes:
        vec = p0 - l0
        if (vec).dot(n) < 0:
            # point lies outside of plane
            return flip
    else:
        return not flip # if flip, point must be inside to not be selected

#--- Unit cell functions -----------------------------------------------------#

def draw_unit_cell(molecule, all_unit_vectors, draw_style='ARROWS'):
    # TODO implement different drawing styles
    all_obs = []
    
    # if unit_vectors is a list of unit vectors reassign them to 
    # all_unit_vectors
    unit_vectors = all_unit_vectors[0]

    me = bpy.data.meshes.new("unit_cell_{}".format(molecule.index))
    coords = (
        (0,0,0), #0, O
        unit_vectors[0], #1, x
        unit_vectors[0] + unit_vectors[1], #2, xy
        unit_vectors[1], #3, y
        unit_vectors[2], #4, z
        unit_vectors[0] + unit_vectors[2], #5, xz
        unit_vectors[1] + unit_vectors[2], #6, yz
        unit_vectors[0] + unit_vectors[1] + unit_vectors[2], #7, xyz
        )
    faces = (
        (0, 3, 2, 1),
        (0, 4, 6, 3),
        (0, 1, 5, 4),
        (7, 6, 4, 5),
        (7, 5, 1, 2),
        (7, 2, 3, 6),
        )
    me.from_pydata(coords, [], faces)
    uc_cube = bpy.data.objects.new("unit_cell_{}".format(molecule.index), me)
    bpy.context.scene.objects.link(uc_cube)
    uc_cube.draw_type = 'WIRE'
    
    vg = []
    vg.append(uc_cube.vertex_groups.new('a'))
    vg[-1].add([1], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('b'))
    vg[-1].add([3], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('c'))
    vg[-1].add([4], 1, 'REPLACE')
    
    all_obs.append(uc_cube)
    
    if 'ARROWS' in draw_style:
        radius = 0.1
        # get material
        material = bpy.data.materials.get('axes')
        if not material:
            material = mb_utils.new_material('axes', color=(1,0,0))
        
        # add sphere at origin
        bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0), size=radius, 
                                             segments=8, ring_count=8)
        ob = bpy.context.object
        ob.name = "unit_cell_origin"
        ob.parent = uc_cube
        ob.parent_type = 'VERTEX'
        ob.data.materials.append(None)
        ob.material_slots[0].material = material
        all_obs.append(ob)
        
        arrow_mesh = mb_utils.get_arrow_data(material=material, radius=radius)
        
        for di, vec in enumerate(unit_vectors):
            ob = bpy.data.objects.new(('a', 'b', 'c')[di], arrow_mesh)
            ob.parent = uc_cube
            ob.parent_type = 'VERTEX'
            bpy.context.scene.objects.link(ob)
            bpy.context.scene.objects.active = ob
            
            c = ob.constraints.new('STRETCH_TO')
            c.name = "mb.stretch"
            c.rest_length = 1.0
            c.volume = 'NO_VOLUME'
            c.target = uc_cube
            c.subtarget = vg[di].name
            all_obs.append(ob)
    
    if all_unit_vectors and len(all_unit_vectors) > 1:
        print("animating uc")
        n_frames = len(all_unit_vectors)
        #bpy.context.scene.frame_end = n_frames - 1
        # add one shape key per frame
        for i, unit_vectors in enumerate(all_unit_vectors):
            shape_key = uc_cube.shape_key_add("{}.{}".format(uc_cube.name, i), 
                                              from_mix=False)
            #shape_key.frame = i+1
            #shape_key.value = 1.
            coords = (
                (0,0,0), #0, O
                unit_vectors[0], #1, x
                unit_vectors[0] + unit_vectors[1], #2, xy
                unit_vectors[1], #3, y
                unit_vectors[2], #4, z
                unit_vectors[0] + unit_vectors[2], #5, xz
                unit_vectors[1] + unit_vectors[2], #6, yz
                unit_vectors[0] + unit_vectors[1] + unit_vectors[2], #7, xyz
                )
            for vert, coord in zip(uc_cube.data.vertices, coords):
                shape_key.data[vert.index].co = coord
        uc_cube.data.shape_keys.use_relative = False
        
        anim_data = uc_cube.data.shape_keys.animation_data_create()
        action = bpy.data.actions.new(name="uc_{}".format(molecule.index))
        anim_data.action = action
        
        fcu = action.fcurves.new(data_path="eval_time")
        #for nf in range(n_frames):
        for nf, kb in enumerate(uc_cube.data.shape_keys.key_blocks):
            fcu.keyframe_points.add(1)
            fcu.keyframe_points[-1].co = nf + 1, kb.frame
            fcu.keyframe_points[-1].interpolation = 'LINEAR'
                    
        #uc_cube.data.shape_keys.eval_time = len(all_unit_vectors) + 1
    
    return all_obs


#--- Other helper functions --------------------------------------------------#



def get_world_coordinates(ob):
    return ob.matrix_world.decompose()[0]


#--- Main functions ----------------------------------------------------------#

def import_molecule(report,
                    filepath,
                    modepath,
                    n_q,
                    qvec,
                    molecule,
                    refine_atoms,
                    refine_bonds,
                    bond_type,
                    scale_distances,
                    bond_guess,
                    use_center,
                    mask_planes,
                    mask_flip,
                    draw_uc,
                    supercell,
                    ):
    
    try:
        start = time.time()
        
        all_obs = []
        all_obs.append(molecule.objects.parent.object)
        
        structure = mb_io_files.MB_Structure.from_file(
            filepath,
            unit_fac=scale_distances,
            modefilepath=modepath,
            n_q=n_q,
            qvec=qvec,
            )
        
        # some sanity checks
        if not structure.all_atoms:
            debug_print("ERROR: No atoms found in {}. ".format(filepath) +
                "Please check file format and/or MolBlend code.",
                level=1)
            report({'ERROR'}, "No atoms found in {}.".format(filepath) +
                "Please check file format and/or MolBlend code.")
            return False
        
        if structure.axes and not len(structure.axes) == structure.nframes:
            raise IOError(("Number of unit vectors ({}) and frames ({})"
                            " does not match").format(len(structure.axes), 
                                                      len(structure.nframes)))
        
        if draw_uc and structure.axes:
            # read unit cell and create cube
            unit_cell_obs = draw_unit_cell(molecule, structure.axes)
            all_obs.extend(unit_cell_obs)
        elif draw_uc and not structure.axes:
            debug_print("WARNING: No unit cell vectors read. Currently "
                        "only supports Quantum Espresso input and output "
                        "and Abinit output.",
                        level=1)
        
        if sum(supercell) > 3:
            structure.create_supercell(supercell)
        
        if bond_guess:
            structure.guess_bonds(tol=0.2)
            debug_print("guess", level=4)
            debug_print(time.time() - start, level=5)
        
        # read normal modes
        
        #if modepath:
            #debug_print("Reading modes file {}".format(modepath), level=4)
            #if filepath.rsplit('.')[-1] == 'guo':
                #frequencies, modes = read_guo_modes(modepath)
            #elif modepath.rsplit('.')[-1] == 'yaml':
                #frequencies, modes = read_phonopy_mode(modepath, all_frames)
            #else:
                ##frequencies, modes = read_modes(modepath, n_q)
                #frequencies, modes = read_manual_modes(modepath)
            #if not modes:
                #debug_print("ERROR: Couldn't read any normal modes in file "
                            #"{}".format(modepath),
                            #level=1)
                #report({'ERROR'}, "Couldn't read any normal modes in file "
                                #"{}".format(modepath))
                #return False
            #else:
                #molecule.max_mode = len(frequencies)
                #modes_col = molecule.modes_col
                #m = modes_col.add()
                #m.index = 0
                #m.name = "mode_0"
                #for i, freq in enumerate(frequencies):
                    #m = modes_col.add()
                    #m.index = i+1
                    #m.name = "mode_{}".format(i+1)
                    #m.freq = freq
        #else:
            #modes = None

        ## replace the index i in data with the corresponding mode
        #if modepath and modes:
            #for all_atoms in all_frames:
                #for index, data in all_atoms.items():
                    ## data = [element, atom_name, location, i]
                    #try:
                        #data[-1] = [mode[data[-1]] for mode in modes]
                    #except (IndexError, ValueError) as e:
                        #debug_print("Error: Modes couldn't be matched with "
                                    #"atoms. ({})".format(e),
                                    #level=1)
                        #report({'ERROR'}, "Modes couldn't be matched with "
                                    #"atoms. ({})".format(e))
                        #modes = None
                        #return False
                #if not modes:
                    #break
        
        debug_print("read", level=4)
        debug_print(time.time() - start, level=5)
        
        # Create parent empty at center of mass
        center_of_mass = structure.get_center_of_mass()
        molecule.objects.parent.object.location = center_of_mass
        debug_print("center", level=4)
        debug_print(time.time() - start, level=5)
        
        # add all atoms to scene
        atom_obs = {}
        error = set()
        debug_print("", level=4)
        for index, atom in sorted(structure.all_atoms.items()):
            debug_print("\ratom {}".format(index), level=4, end='')
            new_atom = mb_utils.add_atom(bpy.context, 
                                         atom["coords"][0]-center_of_mass, 
                                         atom["element"],
                                         atom["name"],
                                         molecule)
            all_obs.append(new_atom)
            # adjust index to be the same as in imported file
            if new_atom.mb.index < index:
                new_atom.mb.index = index
                molecule.atom_index = index + 1
            elif new_atom.mb.index > index:
                error.add("WARNING: Indeces will not be the same as imported.")
            
            atom_obs[index] = new_atom
            
            #if modes and frequencies and len(modes) == len(frequencies):
            ## bpy.data.actions.new is very slow. Only make one action per atom
            ## store mode_vecs in atom object and add drivers to fcurves that
            ## point to this list
                
                #m = new_atom.mb.modes.add()
                #m.name = 0
                #m.index = 0
                #m.freq = 0.0
                #m.vec = (0, 0, 0)
                
                #for i, mode_vec in enumerate(data[3]):
                    #m = new_atom.mb.modes.add()
                    #m.name = i+1
                    #m.index = i+1
                    #m.freq = frequencies[i]
                    #m.vec = mode_vec
                
                #anim_data = new_atom.animation_data_create()
                #atom_id = '{}.{}'.format(new_atom.mb.get_molecule().index, 
                                         #new_atom.mb.index)
                #action = bpy.data.actions.new(name="mode_{}".format(atom_id))
                #anim_data.action = action
                ##for dim in range(3):
                    ##fcu = action.fcurves.new(data_path="location", index=dim)
                    ##fcu.keyframe_points.add(3)
                    ##for p in range(3):
                        ##loc = new_atom.location[dim]
                        ##fcu.keyframe_points[p].co = 1.0 + 10*p, loc
                        ##fcu.keyframe_points[p].interpolation = 'BEZIER'
                ## make new group
                #ag = action.groups.new("Location")
                #for dim in range(3):
                    #fcu = action.fcurves.new(data_path="location", index=dim)
                    #fcu.group = ag
                    #fcu.keyframe_points.add(3)
                    #for p in range(3):
                        #loc = new_atom.location[dim]
                        #fcu.keyframe_points[p].co = 1.0 + 10*p, loc
                        #fcu.keyframe_points[p].interpolation = 'BEZIER'
            #elif len(all_frames) > 1:
            if structure.nframes > 1:
                anim_data = new_atom.animation_data_create()
                atom_id = '{}.{}'.format(new_atom.mb.get_molecule().index, 
                                         new_atom.mb.index)
                action = bpy.data.actions.new(name="frames_{}".format(atom_id))
                anim_data.action = action
                ag = action.groups.new("Location")
                
                for dim in range(3):
                    fcu = action.fcurves.new(data_path="location", index=dim)
                    fcu.group = ag
                    for nf in range(structure.nframes):
                        coord = structure.all_atoms[index]["coords"][nf]
                        loc = (coord-center_of_mass)[dim]
                        fcu.keyframe_points.add(1)
                        fcu.keyframe_points[-1].co = nf + 1, loc
                        fcu.keyframe_points[-1].interpolation = 'LINEAR'
        
        debug_print("", level=4)

        debug_print("atoms", level=4)
        debug_print(time.time() - start, level=5)
        
        # add bonds to scene
        debug_print("", level=4)
        for index1, other in structure.bonds.items():
            for index2 in other:
                debug_print("\rbond {}-{}".format(index1, index2), level=4, 
                            end='')
                new_bond = mb_utils.add_bond(bpy.context, atom_obs[index1], 
                                             atom_obs[index2], 
                                             bond_type=bond_type)
                all_obs.append(new_bond)
        debug_print("", level=4)
        
        debug_print("bonds", level=4)
        debug_print(time.time() - start, level=5)
        
        if error:
            debug_print('\n'.join(error), level=1)
        
        if use_center:
            molecule.objects.parent.object.location -= center_of_mass
        
        for ob in unit_cell_obs[-3:]:
            mb_utils.check_ob_dimensions(ob)
        
        # select all objects and make parent active
        bpy.ops.object.select_all(action="DESELECT")
        bpy.context.scene.objects.active = molecule.objects.parent.object
        for ob in all_obs:
            ob.select = True
    except:
        # if something bad happend, delete all objects and re-raise
        for ob in all_obs:
            bpy.context.scene.objects.unlink(ob)
        raise
        
    return True


