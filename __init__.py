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
    "blender": (2,6),
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
#                                                                           GUI
# TODO When Blender is closed while operator is running, then the MolBlend keyconfig will be active on reload, even if addon is not running
# TODO add frames operator (also for import)
# TODO allow to use NURBS or META for atoms and bonds
# TODO add operator that can "unify" and separate molecules
# TODO add alt+select operator to select all connected atoms
# TODO rethink groups and molecules
# TODO replace refine with global variable, that could even update current meshes after confirmation
# TODO fix use of refine value
# TODO when switching draw styles, check that all molecule parts are on the same layer!
# TODO maybe create a different mesh for refine value (maybe), allowing different molecule to have different refines
# TODO need to make a different mesh for each molecule, element, refine and bond_radius value (maybe)
# TODO streamline deleting atoms, so that unlinking from scene and deleting from molecule is one single function
# TODO change delete shortcut (x and entf) from keymap to include deleting atoms from molecules, splitting molecules etc.
# TODO maybe allow atom scale/radius atomistic property, not molecular
# TODO make atom radii a global setting. Add covalent/vdW switch to atom prop, and add custom scale prop per atom to drivers
# TODO currently one atom can be part of many molecules => add driver that checks for all bond radii
# TODO maybe add driver to vertex group that scales their radius, so that one bond can transition between two "molecules"
# TODO maybe add shift modifier to command add_atom to connect to active atom as an alternative to hovering.
# TODO Cycles material is not updated when diffuse_color is changed (driver issue). Call update_all_meshes before rendering!
# TODO include actual bond lengths list
# TODO find out how to make import for big files faster (like only linking objects to scene after all have been created, or creating objects differently)
# TODO add bond length tool set
# TODO maybe implement change of units within Blender

class MB_PT_tools(Panel):
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_label = "MolBlend - tools"
    
    def draw(self, context):
        layout = self.layout
        if not context.window_manager.mb.is_running:
            ### Start screen
            layout.operator("mb.start")
        else:
            ### Stop screen
            layout.operator("mb.stop")
       
        #--- tools ------------------------------------------------------------#
        box = layout.box()
        box.active = context.window_manager.mb.is_running
        box.label("Main tools")
        box.prop(context.window_manager.mb.globals, "active_tool", expand=True)
        row = box.row()
        col = row.column()
        col.operator("view3d.select_border")
        row = box.row()
        col = row.column()
        col.label("Add")
        col.label("Geometry")
        col = row.column()
        col.prop(context.window_manager.mb.globals, "element_to_add", text="")
        col.prop(context.window_manager.mb.globals, "geometry_to_add", text="")
        
        #row = box.row()
        #row.operator("mb.group_selected", text="Group")
        #row = box.row()
        #row.prop(context.window_manager.mb.globals, "group_selected_extend")
        #row = box.row()
        #row.operator("mb.hover", text="Hover")
        
        #--- atom props -------------------------------------------------------#
        box = layout.box()
        active_ob = context.object
        if context.window_manager.mb.is_running and active_ob and active_ob.mb.type == 'ATOM':
            box.active = True
        else:
            box.active = False
        box.label("Atom properties")
        
        row = box.row()
        col = row.column()
        col.label("Element")
        col.label("Atom radius")
        col.label("Atom color")
        
        col = row.column()
        if active_ob and active_ob.mb.type == 'ATOM':
            col.prop(active_ob.mb, "element", text="")
            #if context.window_manager.mb.is_running:
            col.prop(context.scene.mb.elements[active_ob.mb.element], "covalent", text="")
            col.prop(active_ob.material_slots[0].material, "diffuse_color", text="")

        #--- molecule props ---------------------------------------------------#
        box = layout.box()
        active_ob = context.object
        if (context.window_manager.mb.is_running 
            and active_ob and active_ob.mb.type in ('ATOM', 'BOND')):
            box.active = True
        else:
            box.active = False
        box.label("Molecule properties")
        
        row = box.row()
        col = row.column()
        col.label("Atom scale")
        col.label("Bond radius")
        col.label("Radius type")
        col.label("Draw style")
        col.label("Bond material")
        col.label("Bond color")
        col = row.column()
        if active_ob and active_ob.mb.type in ('ATOM', 'BOND'):
            mol = active_ob.mb.get_molecule()
            col.prop(mol.atom_scales[mol.draw_style], "val", text="")
            col.prop(mol, "bond_radius", text="")
            col.prop(mol, "radius_type", text="")
            col.prop(mol, "draw_style", text="")
            col.prop(mol, "bond_material", text="")
            col.prop(mol, "bond_color", text="")
        
        row = box.row()
        row.operator("mb.center_mol_parent")
        
        #--- global settings --------------------------------------------------#
        box = layout.box()
        box.active = context.window_manager.mb.is_running
        box.label("Global settings")
        
        #--- view setting -----------------------------------------------------#        
        box = layout.box()
        box.active = context.window_manager.mb.is_running
        row = box.row()
        row.label("View")
        row = box.row()
        row.label(icon='AXIS_TOP')
        row.operator("view3d.viewnumpad", text="Top").type = 'TOP'
        row.operator("view3d.viewnumpad", text="Bottom").type = 'BOTTOM'
        row = box.row()
        row.label(icon='AXIS_FRONT')
        row.operator("view3d.viewnumpad", text="Front").type = 'FRONT'
        row.operator("view3d.viewnumpad", text="Back").type = 'BACK'
        row = box.row()
        row.label(icon='AXIS_SIDE')
        row.operator("view3d.viewnumpad", text="Right").type = 'RIGHT'
        row.operator("view3d.viewnumpad", text="Left").type = 'LEFT'
        row = box.row()
        row.label(icon='MANIPUL')
        row.operator("screen.region_quadview", text="Quadview")
        
        #--- import-export-----------------------------------------------------#
        box = layout.box()
        box.active = context.window_manager.mb.is_running
        row = box.row()
        row.label("Import - Export")
        row = box.row()
        row.operator("mb.import_molecule", text="Import")
        row = box.row()
        row.operator("mb.export_molecule", text="Export")
    
#class MB_PT_settings(Panel):
    #bl_space_type = "VIEW_3D"
    #bl_region_type = "TOOLS"
    #bl_label = "MolBlend - settings"
    
    #def draw(self, context):
        #layout = self.layout
        #box = layout.box()
    
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