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

# set path to mcubes. Needs to be of the same python version as Blender.
# If installed with pip you can get the path with "pip show mcubes".
# pymcubes path ###############################################################
mcubes_path = r"/usr/local/lib/python3.5/dist-packages"
###############################################################################

if "bpy" in locals():
    import importlib
    importlib.reload(mb_utils)
    importlib.reload(mb_io_files)
else:
    from molblend import mb_utils
    from molblend import mb_io_files

import math
import logging
import numpy as np

import bpy
import bmesh
from mathutils import Vector, Matrix

from molblend.elements_default import ELEMENTS as ELEMENTS_DEFAULT

logger = logging.getLogger(__name__)
A_per_Bohr = 0.529177249

# Try to load mcubes module for iso surface creation
try:
    import mcubes
except ImportError:
    import os
    if os.path.exists(mcubes_path):
        for module, fn in bpy.path.module_names(mcubes_path):
            if module == "mcubes":
                import sys
                if not mcubes_path in sys.path:
                    sys.path.append(mcubes_path)
                try:
                    import mcubes
                except ImportError:
                    logger.warning("mcubes could not be imported. Is it installed for python 3.5?")
    else:
        logger.error("Path to mcubes library doesn't exists. If you want to use isosurfaces, please install pymcubes and adjust the mcubes_path in mb_import_export.py.")
found_mcubes = ("mcubes" in locals())

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

def import_cube_iso(context,
                    report,
                    filepath,
                    iso_val='VOLFRAC',
                    vol_frac=0.7,
                    absolute=100,
                    origin_to_com=False,
                    ):
    """
    Format specification from 
    http://h5cube-spec.readthedocs.io/en/latest/cubeformat.html
    """
    if not found_mcubes:
        report({'ERROR'}, "mcubes was not found. Please install and/or set correct"
                          " path to mcubes in {}".format(__file__))
        return False
    bpy.ops.object.select_all(action="DESELECT")
    
    with open(filepath, "r") as fin:
        next(fin)
        next(fin)
        ls = next(fin).split()
        nat = int(ls[0])
        dset_ids_present = (nat < 0)
        nat = abs(nat)
        origin = Vector(list(map(float, ls[1:4])))
        nval = int(ls[-1]) if len(ls) == 5 else 1
        if nval != 1 and dset_ids_present:
            report({'ERROR'}, "{}".filepath +
                   "NVAL != 1 and NAT < 0 is not compatible.")
            return False
        nvoxel = np.zeros(3, dtype=int)
        voxel_vec = np.zeros((3,3))
        for i in range(3):
            n, x, y, z = next(fin).split()
            nvoxel[i] = int(n)
            voxel_vec[i,:] = list(map(float, (x, y, z)))
        
        if (nvoxel<0).all():
            unit = 1
        elif (nvoxel<0).any():
            msg = (
                "{} ".format(filepath) +
                "seems to contain mixed units (+/- mixed in lines 4-6). "
                "Please make sure either all units are in Bohr (+) or "
                "Angstrom (-)."
                )
            report({'ERROR'}, msg)
            return False
        else:
            unit = A_per_Bohr
        voxel_vec *= unit
        origin *= unit
        
        # skip atom info
        for n in range(nat):
            next(fin)
        if dset_ids_present:
            orbitals = []
            ls = list(map(int, next(fin).split()))
            m = ls[0]
            orbitals.extend(ls[1:])
            while len(orbital) < m:
                orbitals.extend(list(map(int, next(fin).split())))
        else:
            m = 1
        
        all_data = np.zeros(np.product(nvoxel)*m*nval)
        pos = 0
        for line in fin:
            ls = line.split()
            all_data[pos:(pos+len(ls))] = list(map(float, ls))
            pos += len(ls)
    
    all_data = all_data.reshape(list(nvoxel)+[m*nval])
    all_data = np.rollaxis(all_data, -1, 0)
    
    red = (.8, 0, 0)
    blue = (0, 0, .8)
    n_gt_1 = (len(all_data) > 1)
    for n, data in enumerate(all_data):
        plusminus = (data.min() < 0 and data.max() > 0)
        for fac in (1, -1): # for positive and negative valued wavefunction
            part = data[data*fac > 0]*fac
            if len(part) > 0:
                if iso_val == 'VOLFRAC':
                    flat = np.sort(part.flatten())[::-1]
                    cs = np.cumsum(flat)
                    # find first index to be larger than the percent volume we want
                    # to enclose
                    cut = cs[-1]*vol_frac
                    idx = np.argmax(cs > cut)
                    # linearly interpolate (this should be good enough for most 
                    # purposes
                    rat = (cs[idx]-cut) / (cs[idx]-cs[idx-1])
                    iso = rat * (flat[idx-1]-flat[idx]) + flat[idx]
                elif iso_val == 'ABSOLUTE':
                    iso = absolute
                # Use marching cubes to obtain the surface mesh of these ellipsoids
                verts, faces = mcubes.marching_cubes(data*fac, iso)
                
                # displace by have a voxel to account for the voxel volume
                verts += np.array((0.5, 0.5, 0.5))
                # convert to cartesian coordinates
                verts = verts.dot(voxel_vec)
                
                verts = verts.tolist()
                faces = faces.astype(int).tolist()
                base = bpy.path.display_name_from_filepath(filepath)
                orb = ("_{}".format(n) if n_gt_1 else "")
                ext = ("_p" if fac == 1 else "_n") if plusminus else ""
                name = "{}{}{}".format(base, orb, ext)
                me = bpy.data.meshes.new(name)
                me.from_pydata(verts, [], faces)

                ob = bpy.data.objects.new(name, me)
                context.scene.objects.link(ob)
                context.scene.objects.active = ob
                ob.select = True
                ob.location = origin
                
                bpy.ops.object.mode_set(mode='EDIT')
                bpy.ops.mesh.remove_doubles(threshold=0.06)
                bpy.ops.mesh.normals_make_consistent()
                bpy.ops.mesh.delete_loose()
                bpy.ops.object.mode_set(mode='OBJECT')
                
                bpy.ops.object.shade_smooth()
                
                mod = ob.modifiers.new("Subsurf", 'SUBSURF')
                mod.show_viewport = False
                mod.show_render = False
                
                mod = ob.modifiers.new("Remesh", 'REMESH')
                mod.octree_depth = 8
                mod.scale = 0.99
                mod.use_smooth_shade = True
                mod.use_remove_disconnected = False
                mod.mode = 'SMOOTH'
                mod.show_viewport = False
                mod.show_render = False
                
                if len(ob.material_slots) < 1:
                    ob.data.materials.append(None)
                # get or create element material per molecule
                material = bpy.data.materials.get(name)
                if not material:
                    color = (red if fac == 1 else blue)
                    # get material color from elements list, and Default if not an element
                    material = mb_utils.new_material(name, color=color)
                # finally, assign material to first slot.
                ob.material_slots[0].link = 'DATA'
                ob.material_slots[0].material = material
    
    if origin_to_com:
        bpy.ops.object.origin_set(type="ORIGIN_GEOMETRY")
                

def import_modes(context,
                 report,
                 modefilepath,
                 file_format,
                 molecule):
    
    logger.info("Reading modes file {}".format(modefilepath))
    try:
        qpts = mb_io_files.MB_Modes.from_file(modefilepath, 
                                              file_format)
    except:
        raise
    
    if not qpts or not qpts[0].modes:
        msg = "No modes found in {}\n".format(bpy.path.basename(modefilepath))
        msg += "Did you chose the correct file format?"
        report({'ERROR'}, msg)
        return False
    # TODO Check for correct number of evecs
    nat = len(molecule.objects.atoms)
    for qmode in qpts:
        if qmode.qvecs_format and not molecule["unit_cells"]:
            msg = "Can't convert qvecs to crystal coordinates because no unit"
            msg += " cell information is present"
            logger.error(msg)
            report({'ERROR'}, msg)
        for mode in qmode.modes:
            if nat % len(mode.evecs) != 0:
                msg = "number of displacement vectors "
                msg += "{}".format(len(mode.evecs))
                msg += " is different than number of atoms {}".format(nat)
                msg += " in active molecule."
                logger.error(msg)
                report({'ERROR'}, msg)
                return False
    
    mb_utils.clear_modes(molecule)
    molecule.active_mode = 0
    
    # This is only used for Quantum ESPRESSO, to calculate the crystal unit 
    # cell. The cartesian unit cell is not written in the mode output, so the
    # conversion can't be done when reading the data.
    uc = Matrix(molecule["unit_cells"][0])*1.889725989
    k_uc = Matrix([uc[(dim+1)%3].cross(uc[(dim+2)%3]) for dim in range(3)])
    fac = 2 * math.pi / uc[0].dot(uc[1].cross(uc[2]))
    k_uc = k_uc * fac
    inv_k_uc = k_uc.inverted()
    
    for nq, qmode in enumerate(qpts):
        qm = molecule.qpts.add()
        qm.iqpt = qmode.iqpt
        if qmode.qvecs_format == "QE":
            qm.qvec = (Vector(qmode.qvec) * 2 * math.pi / uc[0][0]) * inv_k_uc
        else:
            qm.qvec = qmode.qvec
        
        mode_name_fmt = 'mode_{}.{}.{{}}.{{}}'.format(molecule.name, nq)
        
        # add mode 0 as the equilibrium position
        m = qm.modes.add()
        m.freq = "equilibrium"
        for atom in molecule.objects.atoms:
            d = m.evecs.add()
            d.real = Vector((0,0,0))
            d.imag = Vector((0,0,0))
        
        for nm, mode in enumerate(qmode.modes):
            m = qm.modes.add()
            m.freq = mode.freq
            for disp in mode.evecs:
                d = m.evecs.add()
                d.real = disp.real
                d.imag = disp.imag
    
    for atom in molecule.objects.atoms:
        mb_utils.create_mode_action(context, atom, molecule)
    
    if len(qpts) > 1 and file_format == "PHONOPY":
        report({'WARNING'}, "Please check if q!=0 modes have been imported"
               " properly. If not, please notify Flo to fix it.")
        
    
def import_molecule(context,
                    report,
                    filepath,
                    molecule,
                    refine_atoms,
                    refine_bonds,
                    bond_material,
                    bond_type,
                    auto_unit,
                    scale_distances,
                    bond_guess,
                    put_origin,
                    parent_center,
                    mask_planes,
                    mask_flip,
                    draw_uc,
                    supercell,
                    ):
    
    try:
        
        all_obs = []
        all_obs.append(molecule.objects.parent)
        
        structure = mb_io_files.MB_Structure.from_file(
            filepath,
            auto_unit=auto_unit,
            unit_fac=scale_distances,
            )
        
        # some sanity checks
        if not structure.all_atoms:
            msg = "No atoms found in {}. ".format(filepath)
            msg += "Please check file format and/or MolBlend code."
            logger.error(msg)
            report({'ERROR'}, msg)
            return False
        
        logger.debug("Found {} frames in {}".format(structure.nframes,
                                                    filepath))

        molecule["unit_cells"] = structure.axes
        
        if draw_uc and molecule["unit_cells"]:
            # read unit cell and create cube
            unit_cell_obs = mb_utils.draw_unit_cell(molecule, context)
            all_obs.extend(unit_cell_obs)
        elif draw_uc and not molecule["unit_cells"]:
            msg = "No unit cell vectors read."
            logger.warning(msg)
            report({'WARNING'}, msg)
        
        if sum(supercell) > 3:
            structure.create_supercell(supercell)
        
        if mask_planes:
            structure.apply_mask(mask_planes, mask_flip)
        
        if bond_guess:
            structure.guess_bonds(tol=0.2)
        
        if parent_center:
            center_of_mass = structure.get_center_of_mass()
        else:
            center_of_mass = Vector(structure.origin)
        
        molecule.objects.parent.location = center_of_mass
        
        # add all atoms to scene
        atom_obs = {}
        warning = set()
        for index, atom in sorted(structure.all_atoms.items()):
            new_atom = mb_utils.add_atom(context, 
                                         atom["coords"][0]-center_of_mass, 
                                         atom["element"],
                                         atom["name"],
                                         atom["id"],
                                         molecule)
            new_atom.mb.supercell = atom.get("supercell", (0,0,0))
            all_obs.append(new_atom)
            new_atom.mb.index = index
            atom_obs[index] = new_atom
            
            if structure.nframes > 1:
                anim_data = new_atom.animation_data_create()
                atom_id = '{}.{}'.format(new_atom.mb.get_molecule().name, 
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
        molecule.atom_index = index
        # add bonds to scene
        for index1, other in structure.bonds.items():
            for index2 in other:
                try:
                    new_bond = mb_utils.add_bond(context, atom_obs[index1],
                                                atom_obs[index2],
                                                bond_type=bond_type)
                except KeyError as err:
                    if not err.args: 
                        err.args=('',)
                    msg = "Atom index {} was found in CONECT"
                    msg += " but not in the list of atoms"
                    msg = msg.format(err.args[0])
                    err.args = (msg, )
                    report({'ERROR'}, msg)
                    raise
                all_obs.append(new_bond)
        molecule.bond_material = bond_material
        
        if warning:
            logger.warning('\n'.join(warning))
        
        if put_origin:
            molecule.objects.parent.location -= center_of_mass
        
        # select all objects and make parent active
        bpy.ops.object.select_all(action="DESELECT")
        context.scene.objects.active = molecule.objects.parent
        #for ob in all_obs:
            #ob.select = True
        
        return True

    except:
        # if something bad happend, delete all objects and re-raise
        logger.debug("Trying to delete all newly imported objects.")
        report({'ERROR'}, "Something went wrong in import. Check console.")
        for ob in all_obs:
            context.scene.objects.unlink(ob)
            bpy.data.objects.remove(ob)
        logger.exception('')
        return False
        


