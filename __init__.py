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


#  Author: Florian Altvater (florian.altvater@gmail.com
#
#  Start of project              : 2014-03-28

#  Acknowledgements 
#  ================
#  Fluid Designer scripts
#  pbd/xyz-import addons by Clemens Barth

bl_info = {
    "name": "MolBlend",
    "description": "Work with chemical structures",
    "author": "Florian Brown-Altvater",
    "version": (0,1),
    "blender": (2,7),
    "location": "View3D > Tool Shelf > MolBlend",
    "warning": "under heavy development",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Add Mesh"
}

if "bpy" in locals():
    import imp
    imp.reload(mb_datastructure)
    imp.reload(mb_operators)
else:
    from molblend import mb_datastructure
    from molblend import mb_operators

import bpy
from bpy.types import Panel
#from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.props import (StringProperty,
                       BoolProperty,
                       EnumProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       PointerProperty,
                       CollectionProperty)

# TODO BEFORE PUBLISHING
# - test all file formats, with multiple frames etc.

#-- Known Problems/Missing Features ------------------------------------------#
# I/O:
# - Reading input files is not the most stable. Mainly because some input files
#   don't follow strict formatting (for example XYZ files).
# - PDB files only support < 100k atoms.
# - Ignores double bonds in PDB files
# - Bond guessing can get slow for big structures (or large supercells)
# - QE input files are read in Angstrom, and scale_distances can screw that up
#
# Structures:
# - No double/triple bond drawing.
# - For bond guessing to work (reasonably), the atomic coordinates must be in
#   Angstrom

# -----------------------------------------------------------------------------
#                               GUI
### FIXES
# TODO When name of Molecule changes, change names of parent, dipoles etc.
# TODO Duplicating whole molecules doesn't duplicate the MB_Molecule, so atom colors, drawing styles etc are still linked.
# TODO Initialize on startup/load (like in manuelbastioni).
# TODO clean up deleting objects. There should be ONE function to call when an object needs to be deleted
#      that takes care of removing it from the molecule, and from Blender itself.
# TODO add unit cell objects to Molecule
# TODO clean up import structure (put Blender methods in mb_util etc. like draw_unit_cell
# TODO Override Delete operator to delete objects from collections when deleted
# TODO Make PEP 8 compatible
# TODO Go over comments
# TODO Make debug_print more consistent (and add more for higher level verbosity)
# TODO Clean up import_export
# TODO Let user select multiple files to import
# TODO Duplication introduces lots of errors because all the properties are the same
# TODO copy pasting from other files might make biiig problems!
# TODO introduce better Error handling: when anything crashes, clean up afterwards! (Like on import)
# TODO When MolBlend is not running and the user deletes an atom/bond, it doesn't take care of deleting it from the molecules collections. So add a function that checks at each start of MolBlend if a molecule's objects still exist, and if not, delete the name from the collection. Also needs to update all bonds and bonded_atoms lists. => Very expensive!!!
# TODO when switching draw styles, check that all molecule parts are on the same layer!
# TODO Cycles material nodesocketcolor is not updated when diffuse_color is changed (driver issue). Call update_all_meshes before rendering!
# TODO Geometries when adding atom

### FEATURES
# TODO add abinit input file
# TODO switch to ID properties. Easier to link to objects that way.
# TODO implement masking objects again
# TODO import q!=0 modes with phase
# TODO display bond angles/dihedrals?
# TODO add option to work in different unit in Blender (so that bond guessing works as well)
# TODO add GUI option to adjust bond tolerance and update bondlist/create new bonds automatically
# TODO add GUI option to ignore hydrogen atoms (or custom set of atoms)
# TODO add unit cell transforms (shift by lattice parameter, add unit cell, remove unit cell, modify super cell size)
# TODO frame import
# TODO export operator
# TODO "charge" objects
# TODO add frequencies to molecule, and think about plotting spectrum in blender
# TODO Think about modes after atoms are moved
# TODO Operator to combine molecule into mesh: apply and delete all drivers!
# TODO Allow to turn off hover for add atom operator. Get's really slow with lots of objects.
# TODO array modifier for unit cells. Have bonds go into next unit cell. Add apply operator as well
# TODO BoolProp to draw unit cell
# TODO clean up code and add comments
# TODO add more analytic tools for phonons like showing frequencies/spectrum etc.
# TODO allow more file formats for vibrational modes
# TODO add frames operator (also for import)
# TODO allow to use NURBS or META for atoms and bonds (http://johnnygizmo.blogspot.nl/2014/06/striping-curve-in-blendercycles.html for uv shader)
# TODO add operator that can "unify" and separate molecules
# TODO add alt+select operator to select all connected atoms
# TODO rethink groups and molecules
# TODO replace refine with global variable, that could even update current meshes after confirmation
# TODO add bond length tool set
# TODO streamline deleting atoms, so that unlinking from scene and deleting from molecule is one single function
# TODO change delete shortcut (x and entf) from keymap to include deleting atoms from molecules, splitting molecules etc.
# TODO display list of atoms of a molecule and allow changing the order (will affect how it is exported)

### MAYBES
# TODO make bond guessing more efficient
# TODO each new molecule is drawn with default values. Get values from last active object? Or from last changed values?
# TODO Haven't included Scene management (like delete globally etc.)
# TODO maybe add driver to vertex group that scales their radius, so that one bond can transition between two "molecules"
# TODO maybe add shift modifier to command add_atom to connect to active atom as an alternative to hovering.
# TODO include actual bond lengths in list (to allow ctrl to use that length)
# TODO maybe allow atom scale/radius atomistic property, not molecular
# TODO maybe implement change of units within Blender
# TODO maybe create a different mesh for refine value (maybe), allowing different molecule to have different refines


class MolBlendPanel():
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "MolBlend"

class MB_PT_initialize(MolBlendPanel, Panel):
    bl_label = "MolBlend"
    
    @classmethod
    def poll(cls, context):
        return not context.scene.mb.is_initialized
    
    def draw(self, context):
        layout = self.layout
        layout.operator("mb.initialize")

class MB_PT_import(MolBlendPanel, Panel):
    bl_label = "Import"
    
    @classmethod
    def poll(cls, context):
        return context.scene.mb.is_initialized
    
    def draw(self, context):
        initialized = len(context.scene.mb.elements) > 0
        
        mb = context.scene.mb.globals
        layout = self.layout
        layout.operator("mb.import_molecules")


class MB_PT_tools(MolBlendPanel, Panel):
    bl_label = "Tools"
    
    @classmethod
    def poll(cls, context):
        return context.scene.mb.is_initialized
    
    def draw(self, context):
        layout = self.layout
       
        modal = mb_operators.MB_OT_modal_add.is_running()
        label = "ESC to stop" if modal == True else "Draw atoms"
        layout.operator("mb.modal_add", text=label)
        #--- tools -----------------------------------------------------------#
        split = layout.split(0.4)
        col = split.column()
        col.label("Add")
        col.label("Geometry")
        col = split.column()
        col.prop(context.scene.mb.globals, "element_to_add", text="")
        col.prop(context.scene.mb.globals, "geometry_to_add", text="")
        layout.operator("mb.select_bonded")
        layout.operator("mb.select_molecule")
        layout.operator("mb.combine_molecules")
        layout.prop(context.scene.mb.globals, "show_bond_lengths")

        layout.separator()
        layout.operator("mb.make_static")
        

class MB_PT_atom(MolBlendPanel, Panel):
    bl_label = "Atom properties"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized
                and hasattr(active_ob, 'mb') and active_ob.mb.type == 'ATOM')
   
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        active_ob.mb.draw_properties(context, layout, active_ob)


class MB_PT_molecule_properties(MolBlendPanel, Panel):
    bl_label = "Molecule properties"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized 
                and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule())
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        mol = active_ob.mb.get_molecule()
        mol.draw_properties(layout)


class MB_PT_molecule_dipole(MolBlendPanel, Panel):
    bl_label = "Molecule dipole"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized 
                and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule())
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        mol = active_ob.mb.get_molecule()
        mol.draw_dipole_props(layout)


class MB_PT_molecule_unit_cell(MolBlendPanel, Panel):
    bl_label = "Molecule unit cell"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized 
                and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule())
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        mol = active_ob.mb.get_molecule()
        mol.draw_unit_cell_props(layout)


class MB_PT_vibration_properties(MolBlendPanel, Panel):
    bl_label = "Molecule vibrations"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized
                and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule())
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        mol = active_ob.mb.get_molecule()
        mol.draw_vibrations(layout)


class MB_PT_molecule_draw_styles(MolBlendPanel, Panel):
    bl_label = "Molecule draw styles"
    
    @classmethod
    def poll(cls, context):
        active_ob = context.object or context.scene.mb.modal_last_active
        return (context.scene.mb.is_initialized
                and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule())
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object or context.scene.mb.modal_last_active
        mol = active_ob.mb.get_molecule()
        mol.draw_styles(layout)


class MB_PT_view(MolBlendPanel, Panel):
    bl_label = "View"
    
    def draw(self, context):
        layout = self.layout
        #row.label("View")
        row = layout.row()
        row.label(icon='AXIS_TOP')
        row.operator("view3d.viewnumpad", text="Top").type = 'TOP'
        row.operator("view3d.viewnumpad", text="Bottom").type = 'BOTTOM'
        row = layout.row()
        row.label(icon='AXIS_FRONT')
        row.operator("view3d.viewnumpad", text="Front").type = 'FRONT'
        row.operator("view3d.viewnumpad", text="Back").type = 'BACK'
        row = layout.row()
        row.label(icon='AXIS_SIDE')
        row.operator("view3d.viewnumpad", text="Right").type = 'RIGHT'
        row.operator("view3d.viewnumpad", text="Left").type = 'LEFT'
        row = layout.row()
        row.label(icon='MANIPUL')
        row.operator("screen.region_quadview", text="Quadview")


def register():
    bpy.utils.register_module(__name__)
    mb_datastructure.register()


def unregister():
    bpy.utils.unregister_module(__name__)
    mb_datastructure.unregister()


if __name__ == "__main__":
    register()
