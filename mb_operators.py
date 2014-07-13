# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENCE BLOCK *****


if "bpy" in locals():
    import imp
    imp.reload(mb_utils)
    imp.reload(mb_import_export)
else:
    from molblend import mb_utils
    from molblend import mb_import_export

import bpy
from bpy.types import (Operator,
                       PropertyGroup,
                       Menu)
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       IntVectorProperty,
                       BoolVectorProperty,
                       PointerProperty,
                       CollectionProperty,
                       EnumProperty)
from bpy_extras.io_utils import ImportHelper, ExportHelper
import os
from mathutils import Vector
from .helper import debug_print

#from operator import attrgetter
#import time
#from . import seqgen
#import os
#import blf

class MB_OT_start(Operator):
    bl_idname = 'mb.start'
    bl_label = 'Start MolBlend'
    bl_options = {'REGISTER'}
    
    def draw(self, context):
        row = self.layout.row()
        row.label("Python scripts auto execute needs "+
                  "to be enabled in order for this "+
                  "script to run.")
        row = self.layout.row()
        row.prop(context.user_preferences.system, "use_scripts_auto_execute")
    
    def invoke(self, context, event):
        # check if python scripts can be executed. Needed for drivers
        if not context.user_preferences.system.use_scripts_auto_execute:
            return context.window_manager.invoke_props_dialog(self)
        else:
            return self.execute(context)
    
    def execute(self, context):
        if not context.user_preferences.system.use_scripts_auto_execute:
            self.report({'ERROR'}, "Python scripts auto execute not enabled")
            return {'CANCELLED'}
        
        debug_print('Start MolBlend', 1)
        wm = context.window_manager
        wm.mb.is_running = True
        
        # store last active keyconfig name
        wm.mb.last_keyconfig = wm.keyconfigs.active.name
        
        keyconfig = mb_utils.create_new_keyconfig(context)
        wm.keyconfigs.active = keyconfig
        
        # initialize elements
        mb_utils.initialize_elements(context)
        # initialize atom scales
        default_scales = {'BALLS': 1.0, 'BAS': 0.3, 'STICKS': 0.001}
        if not len(context.scene.mb.globals.atom_scales):
            for style in ('BALLS', 'BAS', 'STICKS'):
                atom_scale = context.scene.mb.globals.atom_scales.add()
                atom_scale.name = style
                atom_scale.val = default_scales[style]
        
        # don't show parent lines
        context.space_data.show_relationship_lines = False
        
        bpy.ops.mb.modal('INVOKE_DEFAULT')
        return {'FINISHED'}
        
class MB_OT_modal(Operator):
    bl_idname = 'mb.modal'
    bl_label = 'Modal MolBlend'
    bl_options = {'REGISTER'}
    
    
    def __del__(self):
        pass
    
    def modal(self, context, event):
        #print("modal")
        wm = context.window_manager
        if not wm.mb.is_running:
            debug_print('Stopping modal background', 3)
            return {'FINISHED'}
        else:
            if wm.mb.globals.active_tool == 'ADD_ATOM':
                if bpy.ops.object.select_all.poll():
                    bpy.ops.object.select_all(action='DESELECT')
                    hover_ob = mb_utils.return_cursor_object(context, event, mb_type='ATOM')
                    if hover_ob is not None:
                        hover_ob.select = True
                    context.scene.objects.active = hover_ob
            return {'PASS_THROUGH'}
    
    def invoke(self, context, event):
        debug_print('Start modal background', 3)
        if context.area.type == 'VIEW_3D':
            context.window_manager.modal_handler_add(self)

            ## Add the region OpenGL drawing callback
            ## draw in view space with 'POST_VIEW' and 'PRE_VIEW'
            #self._handle = context.region.callback_add(add_dna.draw_callback, (self, context), 'POST_PIXEL')
        return {'RUNNING_MODAL'}

class MB_OT_stop(Operator):
    bl_idname = 'mb.stop'
    bl_label = 'Stop MolBlend'
    bl_options = {'REGISTER'}
    
    def invoke(self, context, event):
        debug_print('Stop MolBlend', 2)
        wm = context.window_manager
        
        wm.mb.is_running = False
        wm.keyconfigs.active = wm.keyconfigs[wm.mb.last_keyconfig]
        
        # turn on parent lines again (TODO check if it was enabled before starting MolBlend)
        context.space_data.show_relationship_lines = True
        
        #kmitems = wm.keyconfigs.active.keymaps['3D View'].keymap_items
        #for i in wm.mb.added_kmitems:
            #kmitems.remove(kmitems.from_id(i))
        #wm.mb.added_kmitems.clear()
        #for i in wm.mb.deactivated_kmitems:
            #kmi = kmitems.from_id(i)
            #kmi.active = True
        #wm.mb.deactivated_kmitems.clear()
        return {'FINISHED'}

class MB_OT_left_click(Operator):
    bl_idname = 'mb.left_click'
    bl_label = 'Left Mouse Click Operator'
    bl_options = {'REGISTER'}
    
    shift = BoolProperty(default=False)
    ctrl = BoolProperty(default=False)
    alt = BoolProperty(default=False)
    #oskey = BoolProperty(default=False)
    action = StringProperty(name='Left mouse action',
                                  default='NONE')
    
    def invoke(self, context, event):
        #context.window_manager.modal_handler_add(self)
        self.properties.action = context.window_manager.mb.globals.active_tool
        self.execute(context)
        #return {'RUNNING_MODAL'}
        return {'FINISHED'}
    
    def execute(self, context):
        action = self.properties.action
        # add atom operator is started as soon as enum item is selected
        if action == 'ADD_ATOM':
            bpy.ops.mb.add_atom('INVOKE_DEFAULT',
                                shift=self.properties.shift,
                                ctrl=self.properties.ctrl,
                                alt=self.properties.alt)
        if action == 'MOVE':
            bpy.ops.transform.translate('INVOKE_DEFAULT')
        elif action == 'ROTATE':
            bpy.ops.transform.rotate('INVOKE_DEFAULT')
        elif action == 'SELECT':
            bpy.ops.mb.select('INVOKE_DEFAULT',
                              shift=self.properties.shift,
                              ctrl=self.properties.ctrl)
        elif action == 'BORDER_SELECT':
            pass
            #bpy.ops.mb.border_select('INVOKE_DEFAULT', extend=self.properties.shift)
        #return {'RUNNING_MODAL'}
        return {'FINISHED'}

class MB_OT_right_click(Operator):
    bl_idname = 'mb.right_click'
    bl_label = 'Right Mouse Click Operator'
    bl_options = {'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        return {'FINISHED'}

    def invoke(self, context, event):
        bpy.ops.wm.call_menu(name="MB_MT_right_click_menu")
        return {'FINISHED'}

class MB_OT_set_tool(Operator):
    '''
    Simple operator that allows to use keyboard shortcuts to set active tool
    '''
    bl_idname = "mb.set_tool"
    bl_label = "Set tool"
    bl_options = {'REGISTER'}
    
    add_atom = BoolProperty(default=False)
    select = BoolProperty(default=False)
    
    def execute(self, context):
        if self.add_atom and not self.select:
            context.window_manager.mb.globals.active_tool = 'ADD_ATOM'
        elif self.select and not self.add_atom:
            context.window_manager.mb.globals.active_tool = 'SELECT'
        else:
            debug_print('WARNING: mb.set_tool properties are set incorrectly.', 1)
        return {'FINISHED'}
    
class MB_OT_add_atom(Operator):
    '''
    Adds an atom at the current mouse position
    '''
    bl_idname = "mb.add_atom"
    bl_label = "Add atom"
    bl_options = {'UNDO', 'REGISTER'}
    
    element = StringProperty(name="Element", default="C")
    #coord_2d = FloatVectorProperty(name="2D position of new atom", 
                                   #size=2, default=(0.0,0.0))
    coord_3d = FloatVectorProperty(name="3D position of new atom", 
                                   size=3, default=(0.0,0.0,0.0), subtype='XYZ')
    depth_location = FloatVectorProperty(size=3, default=(0.0,0.0,0.0), subtype='XYZ')
    
    new_bond_name = StringProperty()
    new_atom_name = StringProperty()
    
    shift = BoolProperty(default=False)
    ctrl = BoolProperty(default=False)
    alt = BoolProperty(default=False)
    
    geometry = EnumProperty(name="Geometry", default='SINGLE', items=mb_utils.enums.geometries,
                      description="Geometry the new bond should be in relative to existing bonds. Press CTRL to activate.")
    
    def mb_atom_objects(self, context):
        items = [(" ", " ", "no bond")]
        items.extend([(ob.name, ob.name, "") for ob in context.scene.objects if ob.mb.type == 'ATOM' and not ob.name == self.new_atom_name])
        return items
    
    first_atom_name = EnumProperty(items=mb_atom_objects,
                                   description="Name of atom to bond the new atom to")
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, "element")
        col = layout.column()
        col.prop(self, "coord_3d", text="Location")
        
        col = layout.column()
        col.prop(self, "first_atom_name", text="Bond to")
        
    def modal(self, context, event):
        #if event.type in ['MOUSEMOVE', 'INBETWEEN_MOUSEMOVE', 'WHEELUPMOUSE',
                #'WHEELDOWNMOUSE', 'WHEELINMOUSE', 'WHEELOUTMOUSE']:
                #pass
        mouse_2d = event.mouse_region_x, event.mouse_region_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(context, mouse_2d, depth=self.depth_location)
        if event.type == 'MOUSEMOVE':
            new_atom = context.scene.objects.get(self.new_atom_name)
            context.scene.objects.active = new_atom
            new_atom.select = True
            #if not new_atom:
                #return {'RUNNING_MODAL'}
            new_bond = context.scene.objects.get(self.new_bond_name)
            first_atom = context.scene.objects.get(self.first_atom_name)
            
            if not event.ctrl and not event.alt:
                hover_ob = mb_utils.return_cursor_object(context, event, exclude=[new_atom], mb_type='ATOM')
                if hover_ob:
                    new_atom.draw_bounds_type = 'SPHERE'
                    new_atom.draw_type = 'BOUNDS'
                    if new_bond:
                        new_bond.constraints["stretch"].target = hover_ob
                else:
                    new_atom.draw_type = 'SOLID'
                    if new_bond:
                        new_bond.constraints["stretch"].target = new_atom
            else:
                new_atom.draw_type = 'SOLID'
                if new_bond:
                    new_bond.constraints["stretch"].target = new_atom
                if event.ctrl and new_bond:
                    # lock into certain angles
                    self.coord_3d = mb_utils.get_fixed_geometry(context, first_atom, new_atom, self.coord_3d, self.geometry)
                if event.alt and new_bond:
                    # constrain length
                    length = 1.0
                    self.coord_3d = mb_utils.get_fixed_length(context, first_atom, self.coord_3d, length)
            
            new_atom.location = self.coord_3d
            # sometimes, when bond is exactly along axis, the dimension goes to zero due to the stretch constraint
            # check for this occurence and fix it
            if new_bond:
                mb_utils.check_bond_dimension(new_bond)
        
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            # check if new atom is above already existing atom
            new_atom = context.object
            first_atom = context.scene.objects.get(self.first_atom_name)
            hover_ob = mb_utils.return_cursor_object(context, event, exclude=[new_atom], mb_type='ATOM')
            if hover_ob:
                mol = new_atom.mb.get_molecule()
                mol.remove_object(new_atom)
                mol.atom_index -= 1
                context.scene.objects.unlink(new_atom)
                bpy.data.objects.remove(new_atom)
                new_bond = context.scene.objects.get(self.new_bond_name)
                if new_bond and hover_ob:
                    first_atom.mb.remove_bond(new_bond)
                    context.scene.objects.unlink(new_bond)
                    bpy.data.objects.remove(new_bond)
                    mb_utils.add_bond(context, first_atom, hover_ob)
                    
            return {'FINISHED'}
        
        return {'RUNNING_MODAL'}
    
    def invoke(self, context, event):
        self.element = context.window_manager.mb.globals.element_to_add
        self.geometry = context.window_manager.mb.globals.geometry_to_add
        hover_ob = mb_utils.return_cursor_object(context, event, mb_type='ATOM')
        
        if hover_ob:
            self.first_atom_name = hover_ob.name
            self.depth_location = hover_ob.location
        else:
            self.first_atom_name = " "
            self.depth_location = context.scene.cursor_location.copy()
        mouse_2d = event.mouse_region_x, event.mouse_region_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(context, mouse_2d, depth=self.depth_location)
        
        ret_exe = self.execute(context)
        
        if 'FINISHED' in ret_exe:
            if context.area.type == 'VIEW_3D':
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                return {'FINISHED'}
        else:
            return ret_exe
    
    def execute(self, context):
        first_atom = context.scene.objects.get(self.first_atom_name)
        # first_atom is only a non-atom if it was manually selected by user.
        ## so throw a warning and stop operator
        if first_atom and first_atom.mb.type != 'ATOM':
            self.first_atom_name = " "
            first_atom = None
        
        if self.first_atom_name.strip() and not first_atom:
            debug_print("Object \"{}\" not found.".format(self.first_atom_name), 1)
            return {'CANCELLED'}
        
        if first_atom:
            molecule = first_atom.mb.get_molecule()
        else:
            molecule = context.scene.mb.new_molecule()
        
        # create a new atom object with the molecule's properties
        new_atom = mb_utils.add_atom(context, self.coord_3d, self.element, self.element, molecule)
        self.new_atom_name = new_atom.name
        
        # add a bond if atom is added to existing molecule
        if first_atom:
            new_bond = mb_utils.add_bond(context, first_atom, new_atom)
            self.new_bond_name = new_bond.name
        
        context.scene.objects.active = new_atom
        new_atom.select = True
        return {'FINISHED'}

class MB_OT_delete(Operator):
    '''
    Delete objects
    '''
    bl_idname = "mb.delete"
    bl_label = "Delete"
    bl_options = {'UNDO', 'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        if context.selected_objects:
            return True
        else:
            False
    
    def execute(self, context):
        for ob in context.selected_objects:
            # If an atom gets deleted, delete all its bonds as well
            if ob.mb.type == 'ATOM' and len(ob.mb.bonds) > 0:
                for bond in ob.mb.bonds:
                    bond_ob = bond.get_object()
                    bond_ob.select = True
                    # make sure that the other bonded atom forgets the bond
                    for a in bond_ob.mb.bonded_atoms:
                        if not a.name == ob.name:
                            a.get_object().mb.remove_bond(bond.get_object())
            # If a bond gets deleted, delete reference in bonded_atoms
            elif ob.mb.type == 'BOND':
                for a in ob.mb.bonded_atoms:
                    a.get_object().mb.remove_bond(ob)
        # now collect all objects per molecule
        molecules = {}
        for ob in context.selected_objects:
            try:
                molecules[ob.mb.get_molecule()].append(ob)
            except KeyError:
                molecules[ob.mb.get_molecule()] = [ob]
        for mol, ob_list in molecules.items():
            mol.remove_objects(ob_list)
        
        bpy.ops.object.delete()
        return {'FINISHED'}
        
    
    
class MB_OT_center_mol_parent(Operator):
    '''
    Custom delete function
    find out what is deleted and update data accordingly
    '''
    bl_idname = "mb.center_mol_parent"
    bl_description = "Put origin to geometric center"
    bl_label = "Center"
    bl_options = {'UNDO', 'REGISTER'}
    
    molecule_name = StringProperty()
    
    @classmethod
    def poll(cls, context):
        try:
            context.object.mb.get_molecule().name
            return True
        except AttributeError:
            return False
    
    def invoke(self, context, event):
        # get molecule from active object
        self.molecule_name = context.object.mb.get_molecule().name
        return self.execute(context)
        
    def execute(self, context):
        if self.molecule_name:
            molecule = context.scene.mb.molecules.get(self.molecule_name)
            if not molecule:
                debug_print("ERROR in mb.center_mol_parent: Molecule {} not found in scene.".format(self.molecule_name))
                return {'CANCELLED'}
            else:
                origin = Vector((0.0,0.0,0.0))
                center = sum([atom.get_object().location for atom in molecule.objects.atoms], origin) / len(molecule.objects.atoms)
                for atom in molecule.objects.atoms:
                    atom.get_object().location -= center
                molecule.parent.get_object().location = center
        return {'FINISHED'}


class MB_MT_group_menu(Menu):
    bl_idname = "MB_MT_group_menu"
    bl_label = "Group list"
    
    def draw(self, context):
        group_names = set()
        for o in context.selected_objects:
            for g in range(len(o.users_group)):
                group_names.add(o.users_group[g].name)
        
        layout = self.layout
        for name in group_names:
            layout.operator("mb.group_selected", text=name, icon='GROUP').active_group_name=name
        layout.operator("mb.group_selected", text="all", icon='GROUP').all_groups=True

class MB_MT_right_click_menu(Menu):
    bl_idname = "MB_MT_right_click_menu"
    bl_label = "Right Click"
    
    def draw(self, context):
        layout = self.layout
        # get active object (or object that mouse hovers over) to display different things
        
        if context.object and context.object.users_group:
            op = layout.operator("mb.select", text="Select group")
            op.group = True

class MB_OT_group_selected(Operator):
    '''
    Adds an atom at the current mouse position
    '''
    bl_idname = "mb.group_selected"
    bl_label = "Group selected"
    bl_options = {'UNDO', 'REGISTER'}
    
    group_name = StringProperty(name="Name", default="Group", description="Group name")
    active_group_name = StringProperty(default="", description="Name of group to add selection to")
    all_groups = BoolProperty(default=False, description="Add selection to all groups")
    
    @classmethod
    def poll(cls, context):
        if context.selected_objects:
            return True
        else:
            #print("Nothing selected.")
            return False
    
    def execute(self, context):
        if self.active_group_name:
            group = bpy.data.groups.get(self.active_group_name)
            for o in context.selected_objects:
                if not o.name in group.objects:
                    group.objects.link(o)
            self.active_group_name = ""
        
        elif self.all_groups:
            group_names = set()
            for o in context.selected_objects:
                for g in range(len(o.users_group)):
                    group_names.add(o.users_group[g].name)
            for name in group_names:
                group = bpy.data.groups.get(name)
                for o in context.selected_objects:
                    if not o.name in group.objects:
                        group.objects.link(o)
            self.all_groups = False
        
        elif not context.window_manager.mb.globals.group_selected_extend:
            group = bpy.data.groups.new(self.group_name)
            for o in context.selected_objects:
                group.objects.link(o)
        
        else:
            last_name = ''
            for o in context.selected_objects:
                if not o.users_group:
                    continue
                elif len(o.users_group) == 1:
                    if last_name and last_name != o.users_group[0]:
                        bpy.ops.wm.call_menu(name="MB_MT_group_menu")
                    else:
                        last_name = o.users_group[0]
                else:
                    bpy.ops.wm.call_menu(name="MB_MT_group_menu")
        return {'FINISHED'}
    
    #def invoke(self, context, event):
        
        
        
        
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Group name: ")
        row = layout.row()
        row.prop(self.properties, "group_name")
        row = layout.row()
        row.prop(context.window_manager.mb.globals, "group_selected_extend")
        
    
class MB_OT_select(Operator):
    '''
    Operator for extending default selection operator
    '''
    bl_idname = "mb.select"
    bl_label = "Select an object"
    bl_options = {'REGISTER', 'UNDO'}
    
    shift = BoolProperty(default=False, description="Extend current selection")
    ctrl = BoolProperty(default=False, description="Select whole molecule")

    group = BoolProperty(default=False, description="Select group selected object belongs to")
    
    @classmethod
    def poll(cls, context):
        return context.area.type == 'VIEW_3D'
    
    def execute(self, context):
        # remember last selection
        selected = context.selected_objects
        last_active = context.object
        if last_active:
            last_active_select = last_active.select
        else:
            last_active_select = None
        
        # if shift wasn't held, deselect all
        if not self.shift:
            bpy.ops.object.select_all(action='DESELECT')
        
        context.scene.objects.active = None
        # invoke selection
        bpy.ops.view3d.select('INVOKE_DEFAULT', extend=self.shift)
        ob = context.object
        
        if self.ctrl:
            # select whole molecule
            pass
            self.ctrl = False
        # TODO group selection with right click selects anything under the cursor.
        elif self.group:
            if ob and ob.users_group:
                bpy.ops.object.select_grouped(extend=self.shift, type='GROUP')
            self.group = False
        
        if ob and self.shift:
            # deselect if it was already selected and make last active object 
            # active again
            if ob in selected:
                ob.select = False
                #context.scene.objects.active = last_active
                if self.ctrl:
                    # deselect whole molecule
                    pass
                    self.ctrl = False
                #elif self.group:
                    #if ob.users_group:
                        #for g in range(len(ob.users_group)):
                            #for o in ob.users_group[g].objects:
                                #o.select = False
                
        return {'FINISHED'}
    
    def invoke(self, context, event):
        return self.execute(context)

#class MB_OT_unify_molecule(Operator):
    #bl_idname = "mb.unify_molecule"
    #bl_label = "Unify the molecular properties of the selected atoms and bonds, using the active object's properties"
    #bl_options = {'REGISTER', 'UNDO'}
    
    #@classmethod
    #def poll(cls, context):
        #if context.object and context.object.mb.type in ('ATOM', 'BOND'):
            #return True
        #else:
            #return False
    
    #def execute(self, context):
        #mol = context.object.mb.get_molecule()
        ## update all atoms first
        #for ob in context.selected_objects:
            #if ob.mb.type == 'ATOM':
                
                
        # then update the bonds
        
    
            
        return {'FINISHED'}
    
class MB_OT_hover(Operator):
    '''
    Operator for extending default selection operator
    '''
    bl_idname = "mb.hover"
    bl_label = "Hover selection"
    bl_options = {'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        return True
    
    def modal(self, context, event):
        #if not context.window_manager.mb.is_running:
            #return {'FINISHED'}
        #else:
            bpy.ops.object.select_all(action='DESELECT')
            hover_ob = mb_utils.return_cursor_object(context, event)
            if hover_ob is not None:
                hover_ob.select = True
            context.scene.objects.active = hover_ob
            #if event.type in ['MOUSEMOVE', 'INBETWEEN_MOUSEMOVE', 'WHEELUPMOUSE',
                #'WHEELDOWNMOUSE', 'WHEELINMOUSE', 'WHEELOUTMOUSE']:
                #pass
            if event.type in ['ESC', 'LEFTMOUSE']:
                return {'FINISHED'}
            
            return {'RUNNING_MODAL'}
    
    def invoke(self, context, event):
        #print("invoke")
        return self.execute(context)
    
    def execute(self, context):
        if context.space_data.type == 'VIEW_3D':
            #print('execute')
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            return {'CANCELLED'}

# This is the class for the file dialog.
class MB_OT_import_molecule(Operator):
    bl_idname = "mb.import_molecule"
    bl_label  = "Import XYZ/PDB (*.xyz,*.pdb)"
    bl_options = {'PRESET', 'UNDO'}
    
    filename_ext = "*.pdb;*.xyz"
    #filter_glob  = StringProperty(default=filename_ext, options={'HIDDEN'},)
    filepath = StringProperty(
            name="File Path",
            description="Filepath used for importing one file",
            maxlen=1024,
            subtype='FILE_PATH',
            )
    #modepath = StringProperty(
            #name="Mode Path",
            #description="Filepath to normal mode file",
            #maxlen=1024,
            #subtype='FILE_PATH',
            #)
    directory = StringProperty(
            name="Directory",
            description="Directory used for importing the file",
            maxlen=1024,
            subtype='FILE_PATH',
            )
    files = CollectionProperty(
        name="File Path",
        description="List with file names used for importing",
        type=bpy.types.OperatorFileListElement,
        )
    
    #--- molecule properties --------------------------------------------------#
    name_mol = StringProperty(name="Molecule Name", default="Molecule",
                              description="Name of imported molecule") # human readable name
    bond_material = EnumProperty(name="Bond material",
                                description="Choose bond material",
                                items=mb_utils.enums.bond_material, default='ATOMS')
    bond_color = FloatVectorProperty(name='Bond color', default=(0.8, 0.8, 0.8),
                                     subtype='COLOR')
    # dito
    draw_style = EnumProperty(name="Display style", items=mb_utils.enums.molecule_styles,
                         description="Style to draw atoms and bonds",
                         default='BAS')
    radius_type = EnumProperty(name="Radius type", items=mb_utils.enums.radius_types,
                               default='covalent')
    bond_radius = FloatProperty(name="Bond radius", default=0.15, min=0.0, max=3.0,
                                description="Radius of bonds for Sticks, and Ball and Sticks")

    # this is a duplicate class from mb_datastructure for 
    class atom_scale(PropertyGroup):
        name = StringProperty()
        val = FloatProperty(name="Atom scale", default=0.3, min=0.0, max=5.0,
                            precision=2)
    
    atom_scales = CollectionProperty(type=atom_scale)
    #ball = EnumProperty(
        #name="Type of ball",
        #description="Choose ball",
        #items=(('NURBS', "NURBS", "NURBS balls"),
               #('MESH', "Mesh" , "Mesh balls"),
               #('META', "Meta" , "Metaballs")),
               #default='NURBS',) 
    #stick = EnumProperty(
        #name="Type of stick",
        #description="Choose ball",
        #items=(('NURBS', "NURBS", "NURBS balls"),
               #('MESH', "Mesh" , "Mesh balls"),
               ##('2', "Meta" , "Metaballs")
               #),
               #default='NURBS',) 
    refine_atoms = IntProperty(name="Refine atoms", default=8, min=3, max=64,
        description="Refine value for atom meshes")
    refine_bonds = IntProperty(name="Refine bonds", default=8, min=3, max=64,
        description="Refine value for atom meshes")
    # TODO this might be handy for different units in files
    #scale_distances = FloatProperty (
        #name = "Distances", default=1.0, min=0.0001,
        #description = "Scale factor for all distances")
    length_unit = EnumProperty(
        name="Unit", default='1.0', items=mb_utils.enums.angstrom_per_unit,
        description="Unit in input file, will be converted to Angstrom")
    length_unit_other = FloatProperty(
        name="Custom Unit", default=1.0, min=0.000001,
        description="Enter conversion factor in Angstrom/unit in file")
    bond_guess = BoolProperty(
        name = "Guess bonds", default=True,
        description = "Guess bonds that are not in the file.")
    use_mask = StringProperty(
        name="Masking object",
        description="Select object that sets boundaries to imported strucure.")
    mask_flip = BoolProperty(
        name="Mask flip",
        description="Invert masking effect (only atoms outside of mask object are imported).")
    supercell = IntVectorProperty(
        name="Supercell", size=3, default=(2,2,2), min=1, subtype='XYZ',
        description="Import a supercell")
    use_center = BoolProperty(
        name = "Object to origin (first frame)", default=False,
        description = "Put the object into the global origin, the first frame only")
    #use_center_all = BoolProperty(
        #name = "Object to origin (all frames)", default=False,
        #description = "Put the object into the global origin, all frames") 
    #use_all_frames = BoolProperty(
        #name = "Load all frames", default=False,
        #description = "Do you want to load all frames?")
    #use_select_frames = BoolProperty(
        #name = "Load select frames", default=False,
        #description = "Load specified frames only.")
    #select_frames = StringProperty(
        #name = "", description="Specify which frames you want to use (e.g \"1, 3-5, 6+7, 3-5\"). If empty, first frame is used. Frames will be processed in the order specified.",
        #maxlen = 256, default = "")
    #skip_frames = IntProperty(
        #name="", default=0, min=0,
        #description="Only use every xth frame.")
    #images_per_key = IntProperty(
        #name="", default=1, min=1,
        #description="Choose the number of images between 2 keys.")
    #interpolation = EnumProperty(
        #name="Interpolation Mode",
        #description="Choose interpolation between keyframes",
        #items=(('BEZIER', "Bezier" , "Smooth interpolation between keyframes"),
               #('LINEAR', "Linear", "Linear interpolation between keyframes"),
               #('CONSTANT', "Constant", "Step-function like interpolation")),
               #default='BEZIER',)
    
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        col = row.column()
        col.prop(self, "name_mol")
        row = layout.row()
        col = row.column()
        col.prop(self, "draw_style")
        
        box = layout.box()
        row = box.row()
        row.label(text="Atoms")
        #row = box.row()
        #col = row.column()
        #col.prop(self, "ball")
        #row.active = (self.ball == "MESH")
        #col = row.column(align=True)        
        #col.prop(self, "mesh_azimuth")
        #col.prop(self, "mesh_zenith")
        row = box.row()
        row.prop(self, "refine_atoms")
        row = box.row()
        row.prop(self, "refine_bonds")
        row = box.row()
        row.prop(self, "radius_type")
        row = box.row()
        row.label(text="Scaling factors")
        row = box.row()
        row.active = (self.draw_style != "STICKS")
        row.prop(self.atom_scales[self.draw_style], "val", text="")
        #row = box.row()
        #row.prop(self, "scale_distances")
        
        box = layout.box()
        box.active = (self.draw_style != "BALLS") # no bonds in balls style
        
        row = box.row()
        row.label(text="Bonds")
        #row = box.row()
        #col = row.column()
        #col.prop(self, "stick")
        row = box.row()
        col = row.column()
        col.prop(self, "bond_radius")
        row = box.row()
        #row.active = (self.stick == 'MESH') # bond sectors only for mesh
        #col = row.column()
        #col.prop(self, "bond_sectors")
        row = box.row()
        col = row.column()
        col.prop(self, "bond_guess")
        row = box.row()
        col = row.column()
        col.prop(self, "bond_material")
        row = box.row()
        row.active = (self.bond_material == 'GENERIC')
        col = row.column()
        col.prop(self, "bond_color")
        
        #box = layout.box()
        #row = box.row()
        #row.label(text="Normal modes")
        #row = box.row()
        #row.prop(self, "modepath")
        
        box = layout.box()
        row = box.row()
        row.label(text="Masking object")
        
        row = box.row()
        row.prop_search(self, "use_mask", bpy.data, "objects", text="")
        row = box.row()
        row.prop(self, "mask_flip")
        
        box = layout.box()
        row = box.row()
        row.label(text="Scene settings")
        
        row = box.row()
        row.prop(self, "use_center")
        
        row = box.row()
        row.prop(self, "length_unit")
        row = box.row()
        row.active = (self.length_unit == 'OTHER')
        row.prop(self, "length_unit_other")
        row = box.row()
        row.prop(self, "supercell")
        
        #row = box.row()
        #row.active = (self.use_all_frames or self.use_select_frames)
        #row.prop(self, "use_center_all")
        
        #box = layout.box()
        #box.active = (self.use_parenting or self.use_stretch_to)
        #row = box.row()
        #row.label(text="Animation/Frames")
        #row = box.row()
        #row.prop(self, "use_all_frames")
        #row = box.row()
        #row.active = self.use_all_frames
        #col = row.column()
        #col.label(text="Skip frames")
        #col = row.column()
        #col.prop(self, "skip_frames")
        #row = box.row()
        #row.prop(self, "use_select_frames")
        #row = box.row()
        #row.active = self.use_select_frames
        #row.prop(self, "select_frames")
        #row = box.row()
        #row.active = (self.use_all_frames or self.use_select_frames)
        #col = row.column()
        #col.label(text="Frames/key")
        #col = row.column()
        #col.prop(self, "images_per_key")
        #row = box.row()
        #row.active = (self.use_all_frames or self.use_select_frames)
        #row.prop(self, "interpolation")
    
    def invoke(self, context, event):
        # before import dialog is opened, initialize atom scales
        if not len(context.scene.mb.globals.atom_scales):
            default_scales = {'BALLS': 1.0, 'BAS': 0.3, 'STICKS': 0.001}
        else:
            default_scales = {}
            for style in ('BALLS', 'BAS', 'STICKS'):
                default_scales[style] = context.scene.mb.globals.atom_scales[style].val
        for style in ('BALLS', 'BAS', 'STICKS'):
            atom_scale = self.atom_scales.add()
            atom_scale.name = style
            atom_scale.val = default_scales[style]
        
        return context.window_manager.invoke_props_dialog(self)
                
    def execute(self, context):
        ##import time
        ##start = time.time()
        filepath = context.window_manager.mb.globals.import_path
        if context.window_manager.mb.globals.import_modes:
            modes_path = context.window_manager.mb.globals.modes_path
        else:
            modes_path = ''
        
        new_molecule = context.scene.mb.new_molecule()
        
        new_molecule.name_mol = self.name_mol
        new_molecule.draw_style = self.draw_style
        new_molecule.radius_type = self.radius_type
        new_molecule.bond_radius = self.bond_radius
        new_molecule.refine_atoms = self.refine_atoms
        new_molecule.refine_bonds = self.refine_bonds
        for scale in self.atom_scales:
            new_molecule.atom_scales[scale.name].val = scale.val
        
        ## check if select_frames is ok, otherwise import first frame only
        error_list = []
        #frame_list = []
        #if self.use_select_frames:
            #try:
                #for item in self.select_frames.split(','):
                    #if '-' in item:
                        #start, end = map(int, item.split('-'))
                        #frame_list.extend(range(start, end+1))
                    #elif '+' in item:
                        #frame_list.extend(map(int, item.split('+')))
                    #else:
                        #frame_list.append(int(item))
            #except (ValueError, TypeError):
                #error_list.append('Format error in the frame list. Using first frame only.')
                #frame_list.append(1)
            #if not frame_list:
                #error_list.append('Frame list was empty. Using first frame only.')
                #frame_list.append(1)
        
        mask = bpy.data.objects.get(self.use_mask)
        mask_planes = []
        if not mask and self.use_mask:
            error_list.append('Mask object not found. Not using mask.')
        elif mask:
            world_mat = mask.matrix_world
            # get all planes/faces (normal vector, point on plane) from mask object
            mask_planes = [(world_mat*pg.normal.copy(), world_mat*mask.data.vertices[pg.vertices[0]].co) for pg in mask.data.polygons]
        
        if error_list:
            debug_print('\n'.join(error_list), level=1)
        
        if self.length_unit == 'OTHER':
            scale_distances = self.length_unit_other
        else:
            scale_distances = float(self.length_unit)
        # Execute main routine
        mb_import_export.import_molecule(
                                filepath,
                                modes_path,
                                new_molecule,
                                #self.ball,
                                self.refine_atoms,
                                self.refine_bonds,
                                scale_distances,
                                #self.stick,
                                #self.bond_sectors,
                                self.bond_guess,
                                self.use_center,
                                #self.use_center_all,
                                #self.use_camera,
                                #self.use_lamp,
                                mask_planes,
                                self.mask_flip,
                                self.supercell,
                                )
            #if self.use_all_frames:
                #frame_list.extend(range(1, len(import_molecule.ALL_FRAMES)+1))
            ## Load frames
            #if len(import_molecule.ALL_FRAMES) > 1 and frame_list and frame_list != [1]:
                
                #import_molecule.build_frames(self.images_per_key, self.skip_frames,
                                        #frame_list, self.interpolation)
        #print("Imported in {} s".format(time.time()-start))
        return {'FINISHED'}

class ExportPDB(Operator, ExportHelper):
    bl_idname = "mb.export_molecule"
    bl_label  = "Export molecule structure to file(*.pdb, *.xyz)"
    filename_ext = "*.pdb;*.xyz"
    filter_glob  = StringProperty(default=filename_ext, options={'HIDDEN'},)
    
    #atom_pdb_export_type = EnumProperty(
        #name="Type of Objects",
        #description="Choose type of objects",
        #items=(('0', "All", "Export all active objects"),
               #('1', "Elements", "Export only those active objects which have"
                                 #" a proper element name")),
               #default='1',) 
    file_type = EnumProperty(
        name="File type", default="XYZ", items=mb_utils.enums.file_types,
        description="File format to export to")
    length_unit = EnumProperty(
        name="Unit", default='1.0', items=mb_utils.enums.angstrom_per_unit,
        description="Unit in input file, will be converted to Angstrom")
    length_unit_other = FloatProperty(
        name="Custom Unit", default=1.0, min=0.000001,
        description="Enter unit of input file as Angstrom/unit")
    selection_only = BoolProperty(name="Selected Objects", default=True,
        description="Only export selected objects")
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, "file_type")
        row = layout.row()
        row.prop(self, "selection_only")
        row = layout.row()
        row.prop(self, "length_unit")
        row = layout.row()
        row.active = (self.length_unit == 'OTHER')
        row.prop(self, "length_unit")

    def execute(self, context):
        if self.length_unit == 'OTHER':
            scale_distances = self.length_unit_other
        else:
            scale_distances = float(self.length_unit)
        
        if self.file_type == 'PDB':
            mb_import_export.export_pdb(bpy.path.abspath(self.filepath),
                                        scale_distances,
                                        self.selection_only)
        elif self.file_type == 'XYZ':
            mb_import_export.export_xyz(bpy.path.abspath(self.filepath),
                                        scale_distances,
                                        self.selection_only)
        return {'FINISHED'}