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

import math
import logging

import bpy
import bmesh
from mathutils import Vector, Matrix

from molblend.elements_default import ELEMENTS as ELEMENTS_DEFAULT

logger = logging.getLogger(__name__)

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

def import_modes(context,
                 report,
                 modefilepath,
                 file_format,
                 molecule):
    
    logger.info("Reading modes file {}".format(modefilepath))
    try:
        qpts = mb_io_files.modes_from_file(modefilepath, 
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
    
    uc = Matrix(molecule["unit_cells"][0])*1.889725989
    #print(uc)
    k_uc = Matrix([uc[(dim+1)%3].cross(uc[(dim+2)%3]) for dim in range(3)])
    fac = 2 * math.pi / uc[0].dot(uc[1].cross(uc[2]))
    k_uc = k_uc * fac
    
    inv_k_uc = k_uc.inverted()
    for nq, qmode in enumerate(qpts):
        qm = molecule.qpts.add()
        qm.nqpt = qmode.nqpt
        if qmode.qvecs_format == "QE":
            qm.qvec = (Vector(qmode.qvec) * 2 * math.pi / uc[0][0]) * inv_k_uc
        else:
            qm.qvec = qmode.qvec
        
        mode_name_fmt = 'mode_{}.{}.{{}}.{{}}'.format(molecule.index, nq)
        
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

def import_molecule(context,
                    report,
                    filepath,
                    molecule,
                    refine_atoms,
                    refine_bonds,
                    bond_material,
                    bond_type,
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
            unit_fac=scale_distances,
            )
        
        # some sanity checks
        if not structure.all_atoms:
            msg = "No atoms found in {}. ".format(filepath)
            msg += "Please check file format and/or MolBlend code."
            logger.error(msg)
            report({'ERROR'}, msg)
            return False
        
        if structure.axes and not len(structure.axes) == structure.nframes:
            raise IOError(("Number of unit vectors ({}) and frames ({})"
                            " does not match").format(len(structure.axes), 
                                                      len(structure.nframes)))

        molecule["unit_cells"] = structure.axes
        
        if draw_uc and molecule["unit_cells"]:
            # read unit cell and create cube
            unit_cell_obs = mb_utils.draw_unit_cell(molecule)
            all_obs.extend(unit_cell_obs)
            for ob in unit_cell_obs[-3:]:
                mb_utils.check_ob_dimensions(ob)
        elif draw_uc and not molecule["unit_cells"]:
            msg = "No unit cell vectors read."
            logger.warning(msg)
            self.report({'WARNING'}, msg)
        
        if sum(supercell) > 3:
            structure.create_supercell(supercell)
        
        if mask_planes:
            structure.apply_mask(mask_planes, mask_flip)
        
        if bond_guess:
            structure.guess_bonds(tol=0.2)
        
        if parent_center:
            center_of_mass = structure.get_center_of_mass()
        else:
            center_of_mass = Vector((0,0,0))
        
        molecule.objects.parent.location = center_of_mass
        
        # add all atoms to scene
        atom_obs = {}
        error = set()
        for index, (old_index, atom) in enumerate(sorted(structure.all_atoms.items())):
            new_atom = mb_utils.add_atom(context, 
                                         atom["coords"][0]-center_of_mass, 
                                         atom["element"],
                                         atom["name"],
                                         molecule)
            new_atom.mb.supercell = atom.get("supercell", (0,0,0))
            all_obs.append(new_atom)
            new_atom.mb.index = index
            molecule.atom_index = index + 1
            if new_atom.mb.index != old_index:
                error.add("WARNING: Indeces will not be the same as imported.")
            
            atom_obs[old_index] = new_atom
            
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
        
        # add bonds to scene
        for index1, other in structure.bonds.items():
            for index2 in other:
                new_bond = mb_utils.add_bond(context, atom_obs[index1],
                                             atom_obs[index2],
                                             bond_type=bond_type)
                all_obs.append(new_bond)
        molecule.bond_material = bond_material
        
        if error:
            logger.error('\n'.join(error))
        
        if put_origin:
            molecule.objects.parent.location -= center_of_mass
        
        # select all objects and make parent active
        bpy.ops.object.select_all(action="DESELECT")
        context.scene.objects.active = molecule.objects.parent
        for ob in all_obs:
            ob.select = True
    except:
        # if something bad happend, delete all objects and re-raise
        for ob in all_obs:
            context.scene.objects.unlink(ob)
            bpy.data.objects.remove(ob)
        raise
        
    return True


