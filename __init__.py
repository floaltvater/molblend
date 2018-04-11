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
    "blender": (2,7,9),
    "location": "View3D > Tool Shelf > MolBlend",
    "warning": "under heavy development",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Add Mesh"
}

if "bpy" in locals():
    import importlib
    importlib.reload(mb_datastructure)
    importlib.reload(mb_operators)
else:
    from molblend import mb_datastructure
    from molblend import mb_operators

import logging
import os

import bpy
from bpy.types import Panel
from bpy.props import (StringProperty,
                       BoolProperty,
                       EnumProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       PointerProperty,
                       CollectionProperty)

log_path = os.path.join(bpy.context.user_preferences.filepaths.temporary_directory, "molblend_log.txt")
log_is_writeable = True

try:
    test_writing = open(log_path, 'w')
    test_writing.close()
except:
    print("WARNING: Writing permission error for {0}".format(log_path))
    print("The log will be redirected to the console (here)")
    log_is_writeable = False

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_formatter = logging.Formatter('%(asctime)s:%(name)-12s:%(levelname)-8s:%(filename)-20s: %(message)s')
if logger.getEffectiveLevel() == logging.DEBUG:
    term_formatter = logging.Formatter('%(asctime)s:%(levelname)-8s:%(filename)-20s: %(message)s')
else:
    term_formatter = logging.Formatter('%(name)s: %(levelname)-8s: %(message)s')

chandler = logging.StreamHandler()
chandler.setLevel(logging.WARNING)
chandler.setFormatter(term_formatter)
logger.addHandler(chandler)

if log_is_writeable:
    fhandler = logging.FileHandler(log_path, mode ='w')
    fhandler.setLevel(logging.INFO)
    fhandler.setFormatter(file_formatter)
    logger.addHandler(fhandler)

logger.info("Writing log to file {}".format(log_path))


class MolBlendPanel():
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "MolBlend"

class MolBlendPropsPanel():
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Molecule"

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
        layout.operator("mb.import_cube_iso")


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
        layout.operator("mb.apply_scale")

class MB_PT_molecule_properties(MolBlendPropsPanel, Panel):
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


class MB_PT_molecule_draw_styles(MolBlendPropsPanel, Panel):
    bl_label = "Draw styles"
    
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


class MB_PT_atom(MolBlendPropsPanel, Panel):
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


class MB_PT_molecule_dipole(MolBlendPropsPanel, Panel):
    bl_label = "Dipole"
    
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


class MB_PT_molecule_unit_cell(MolBlendPropsPanel, Panel):
    bl_label = "Unit cell"
    
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


class MB_PT_vibration_properties(MolBlendPropsPanel, Panel):
    bl_label = "Vibrations"
    
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
