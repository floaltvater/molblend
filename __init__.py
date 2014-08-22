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
#

bl_info = {
    "name": "MolBlend",
    "description": "Work with chemical structures",
    "author": "Florian Altvater",
    "version": (0,1),
    "blender": (2,7),
    "location": "View3D > Tool Shelf > MolBlend",
    "warning": "under development",
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

#import os
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

# -----------------------------------------------------------------------------
#                               GUI
### FIXES
# TODO make import popup smaller (or remove completely)
# TODO normals on dipole arrow head (and probably cell vectors) are pointing inwards!
# TODO check what happens when molecule is duplicated! Maybe fix driver issues and stuff
# TODO introduce better Error handling: when anything crashes, clean up afterwards! (Like on import)
# TODO When MolBlend is not running and the user deletes an atom/bond, it doesn't take care of deleting it from the molecules collections. So add a function that checks at each start of MolBlend if a molecule's objects still exist, and if not, delete the name from the collection. Also needs to update all bonds and bonded_atoms lists. => Very expensive!!!
# TODO when switching draw styles, check that all molecule parts are on the same layer!
# TODO Cycles material nodesocketcolor is not updated when diffuse_color is changed (driver issue). Call update_all_meshes before rendering!
# TODO clean up read functions for QE files (f.seek not pretty)
# TODO copy pasting from other files might make biiig problems!

### FEATURES
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
    
class MB_PT_tools(MolBlendPanel, Panel):
    bl_label = "Tools"
    
    def draw(self, context):
        layout = self.layout
        if not context.window_manager.mb.is_running:
            ### Start screen
            layout.operator("mb.start")
        else:
            ### Stop screen
            layout.operator("mb.stop")
       
        #--- tools ------------------------------------------------------------#
        #layout.label("Main tools")
        row = layout.row()
        row.active = context.window_manager.mb.is_running
        row.prop(context.window_manager.mb.globals, "active_tool", expand=True)
        #col = row.column()
        #col.operator("view3d.select_border")
        row = layout.row()
        col = row.column()
        col.label("Add")
        col.label("Geometry")
        col = row.column()
        col.prop(context.window_manager.mb.globals, "element_to_add", text="")
        col.prop(context.window_manager.mb.globals, "geometry_to_add", text="")
        
        #row = layout.row()
        #row.operator("mb.group_selected", text="Group")
        #row = layout.row()
        #row.prop(context.window_manager.mb.globals, "group_selected_extend")
        #row = layout.row()
        #row.operator("mb.hover", text="Hover")

class MB_PT_atom(MolBlendPanel, Panel):
    bl_label = "Atom properties"
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object
        if active_ob and hasattr(active_ob, 'mb') and active_ob.mb.type == 'ATOM':
            active_ob.mb.draw_properties(context, layout, active_ob)

class MB_PT_molecule_properties(MolBlendPanel, Panel):
    bl_label = "Molecule properties"
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object
        if active_ob and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule():
            mol = active_ob.mb.get_molecule()
            mol.draw_properties(layout)

class MB_PT_molecule_draw_styles(MolBlendPanel, Panel):
    bl_label = "Molecule draw styles"
    
    def draw(self, context):
        layout = self.layout
        active_ob = context.object
        if active_ob and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule():
            mol = active_ob.mb.get_molecule()
            mol.draw_styles(layout)

class MB_PT_import(MolBlendPanel, Panel):
    bl_label = "Import"
    
    def draw(self, context):
        initialized = len(context.scene.mb.elements) > 0
        
        mb = context.scene.mb.globals
        layout = self.layout
        #row = layout.row()
        #row.label("Import")
        row = layout.row()
        row.prop(mb.import_props, "filepath")
        row = layout.row()
        row.prop(mb.import_props, "modes")
        row = layout.row()
        row.active = mb.import_props.modes
        row.prop(mb.import_props, "modes_path")
        row = layout.row()
        row.active = initialized
        row.operator("mb.import_molecule", text="Import")

class MB_PT_export(MolBlendPanel, Panel):
    bl_label = "Export"
    
    def draw(self, context):
        initialized = len(context.scene.mb.elements) > 0
        
        mb = context.scene.mb.globals
        layout = self.layout
        
        #row = layout.row()
        #row.label("Export")
        row = layout.row()
        row.prop(mb.export_props, "filepath")
        row = layout.row()
        row.prop(mb.export_props, "file_type")
        row = layout.row()
        row.prop(mb.export_props, "selection_only")
        row = layout.row()
        row.prop(mb.export_props, "length_unit")
        row = layout.row()
        row.active = (mb.export_props.length_unit == 'OTHER')
        row.prop(mb.export_props, "length_unit_other")
        row = layout.row()
        row.active = initialized
        row.operator("mb.export_molecule", text="Export")

#class MB_PT_global(MolBlendPanel, Panel):
    #bl_label = "Global Settings"
    
    #def draw(self, context):
        #layout = self.layout
        #layout.label("Global settings")
        
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
    #mb_operators.register()
    mb_datastructure.register()
    #bpy.types.INFO_MT_file_import.append(menu_func)
    #bpy.types.INFO_MT_file_export.append(menu_func_export) 
    
def unregister():
    bpy.utils.unregister_module(__name__)
    #mb_operators.unregister()
    mb_datastructure.unregister()
    #bpy.types.INFO_MT_file_import.remove(menu_func)
    #bpy.types.INFO_MT_file_export.remove(menu_func_export)

if __name__ == "__main__":
    register()

'''
structure:

add an atom:
    if alone:
        make new molecule
    if attached:
        add to molecule

remove an atom:
    if it breaks a molecule in two:
        break a bond
    delete from molecule

make a bond between two atoms:
    merge atoms into one molecule

break a bond:
    separate molecules


'''