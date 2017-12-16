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
    imp.reload(mb_geometry)
    imp.reload(mb_import_export)
else:
    from molblend import mb_utils
    from molblend import mb_geometry
    from molblend import mb_import_export

import os
import sys
from operator import methodcaller

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
from mathutils import Vector, Matrix

from molblend.mb_helper import debug_print


class MB_OT_initialize(Operator):
    bl_idname = 'mb.initialize'
    bl_label = 'Initialize MolBlend'
    bl_description = 'Make sure that drivers work and elements are loaded.'
    bl_options = {'REGISTER'}
    
    def draw(self, context):
        layout = self.layout
        row.label("Python scripts auto execute needs "+
                  "to be enabled in order for this "+
                  "script to run.")
        layout.prop(context.user_preferences.system, "use_scripts_auto_execute")
    
    def invoke(self, context, event):
        debug_print("mb.initialize.invoke", level=6)
        # check if python scripts can be executed. Needed for drivers
        if not context.user_preferences.system.use_scripts_auto_execute:
            return context.window_manager.invoke_props_dialog(self)
        else:
            return self.execute(context)
    
    def execute(self, context):
        if not context.user_preferences.system.use_scripts_auto_execute:
            self.report({'ERROR'}, "Python scripts auto execute not enabled")
            return {'CANCELLED'}
        
        debug_print('Initialize MolBlend', level=2)
        wm = context.window_manager
        
        # initialize elements
        mb_utils.initialize_elements(context)
        # initialize atom scales
        default_scales = {'BALLS': 1.0, 'BAS': 0.4, 'STICKS': 0.001}
        if not len(context.scene.mb.globals.atom_scales):
            for style in ('BALLS', 'BAS', 'STICKS'):
                atom_scale = context.scene.mb.globals.atom_scales.add()
                atom_scale.name = style
                atom_scale.val = default_scales[style]
        
        # don't show parent lines
        context.space_data.show_relationship_lines = False
        context.scene.mb.is_initialized = True
        return {'FINISHED'}


class MB_OT_modal_add(Operator):
    bl_idname = 'mb.modal_add'
    bl_label = 'activate MolBlend'
    bl_options = {'REGISTER'}
    
    is_running_bool = BoolProperty(
        name="Modal_is_running", 
        description="Knows if main modal operator is running",
        default=False)
    
    @classmethod
    def is_running(cls):
        return cls.is_running_bool
    
    def modal(self, context, event):
        if event.type in ('ESC',) or not type(self).is_running_bool:
            return self.cancel(context)
        #print("modal")
        # get 3D Window region
        
        if context.object and context.object.mode == "EDIT":
            self.report({'ERROR'}, 
                        "MolBlend modal operator doesn't work in edit mode.")
            return self.cancel(context)
        
        min_max_lst = []
        for area in context.screen.areas:
            if area.type == "VIEW_3D":
                min_x, min_y = (10000, 10000)
                max_x, max_y = (-10000, -10000)
                for region in area.regions:
                    if region.type == "WINDOW":
                        if region.x < min_x:
                            min_x = region.x
                        if region.y < min_y:
                            min_y = region.y
                        if region.x+region.width > max_x:
                            max_x = region.x+region.width
                        if region.y+region.height > max_y:
                            max_y = region.y+region.height
                min_max_lst.append((min_x, min_y, max_x, max_y))
        x, y = event.mouse_x, event.mouse_y
        for min_max in min_max_lst:
            if (min_max[0] < x < min_max[2] and
                min_max[1] < y < min_max[3]):
                break
        else:
            context.window.cursor_modal_restore()
            return {'PASS_THROUGH'}
        
        # cursor in View3D Window, continue
        context.window.cursor_modal_set("CROSSHAIR")
        
        bpy.ops.object.select_all(action='DESELECT')
        hover_ob = mb_utils.return_cursor_object(context, event,
                                                mb_type='ATOM')
        if hover_ob is not None:
            hover_ob.select = True
            context.scene.mb.modal_last_active = hover_ob
        context.scene.objects.active = hover_ob
    
        if (event.type == 'LEFTMOUSE' and event.value == 'PRESS'):
            bpy.ops.mb.add_atom('INVOKE_DEFAULT',
                                shift=event.shift,
                                ctrl=event.ctrl,
                                alt=event.alt)
            context.scene.mb.modal_last_active = context.scene.objects.active
            return {'RUNNING_MODAL'}
        return {'PASS_THROUGH'}
        
        
    def invoke(self, context, event):
        if not context.scene.mb.is_initialized:
            bpy.ops.mb.initialize("INVOKE_DEFAULT")
        
        if context.area.type == 'VIEW_3D':
            if type(self).is_running_bool == True:
                self.cancel(context)
            else:
                type(self).is_running_bool = True
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            return {'CANCELLED'}

    def cancel(self, context):
        debug_print("cancel modal", level=6)
        type(self).is_running_bool = False
        context.window.cursor_modal_set('DEFAULT')
        context.scene.objects.active = context.scene.mb.modal_last_active
        context.scene.mb.modal_last_active = None
        return {'CANCELLED'}

class MB_OT_make_static(Operator):
    '''
    Apply and remove constraints of bonds
    '''
    bl_idname = "mb.make_static"
    bl_label = "Make static"
    bl_options = {'UNDO', 'REGISTER'}
    bl_description = "apply and remove bond constraints"

    def invoke(self, context, event):
        return self.execute(context)
    
    def remove_constraints(self, ob):
        mat = ob.matrix_world.copy()
        for cname in ("mb.stretch", "mb.parent"):
            c = ob.constraints.get(cname, None)
            if c:
                ob.constraints.remove(c)
        ob.parent = ob.mb.get_molecule().objects.parent
        ob.matrix_world = mat
    
    def execute(self, context):
        if not context.selected_editable_objects:
            debug_print("No objects selected", level=2)
        for ob in context.selected_editable_objects:
            if ob.mb.type == 'BOND':
                self.remove_constraints(ob)
            elif ob.mb.type == 'ATOM':
                for bond in ob.mb.bonds:
                    self.remove_constraints(bond)
                
        return {'FINISHED'}

class MB_OT_add_atom(Operator):
    '''
    Adds an atom at the current mouse position
    '''
    bl_idname = "mb.add_atom"
    bl_label = "Add atom"
    bl_options = {'UNDO', 'REGISTER'}
    
    element = StringProperty(name="Element", default="C")
    coord_3d = FloatVectorProperty(
        name="3D position", description="3D position of new atom",
        size=3, default=(0.0,0.0,0.0), subtype='XYZ')
    parent_coord_3d = FloatVectorProperty(
        name="3D position", description="3D position of parent molecule",
        size=3, default=(0.0,0.0,0.0), subtype='XYZ')
    depth_location = FloatVectorProperty(
        name="Depth", description="Depth of the new atom",
        size=3, default=(0.0,0.0,0.0), subtype='XYZ')
    
    new_bond_name = StringProperty()
    new_atom_name = StringProperty()
    
    shift = BoolProperty(default=False)
    ctrl = BoolProperty(default=False)
    alt = BoolProperty(default=False)
    
    geometry = EnumProperty(
        name="Geometry",
        description="Geometry the new bond should be in relative to "
                    "existing bonds. Press CTRL to activate.",
        items=mb_utils.enums.geometries, default='NONE')
    
    def mb_atom_objects(self, context):
        items = [(" ", " ", "no bond")]
        items.extend(
            [(ob.name, ob.name, "") for ob in context.scene.objects
             if ob.mb.type == 'ATOM' and not ob.name == self.new_atom_name]
            )
        return items
    
    first_atom_name = EnumProperty(
        name="Atom name", description="Name of atom to bond the new atom to",
        items=mb_atom_objects)
    
    def draw(self, context):
        layout = self.layout
        layout.prop(self, "element")
        col = layout.column()
        col.prop(self, "coord_3d", text="Location")
        col = layout.column()
        col.prop(self, "first_atom_name", text="Bond to")
        
    def modal(self, context, event):
        mouse_2d = event.mouse_x, event.mouse_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(
            context, mouse_2d, region=self.region, 
            rv3d=self.rv3d, depth=self.depth_location)
        
        if event.type == 'MOUSEMOVE':
            new_atom = context.scene.objects.get(self.new_atom_name)
            context.scene.objects.active = new_atom
            new_atom.select = True
            new_bond = context.scene.objects.get(self.new_bond_name)
            first_atom = context.scene.objects.get(self.first_atom_name)
            
            hover_ob = mb_utils.return_cursor_object(
                context, event, exclude=[new_atom], mb_type='ATOM')
            if hover_ob:
                new_atom.draw_bounds_type = 'SPHERE'
                new_atom.draw_type = 'BOUNDS'
                if new_bond:
                    new_bond.constraints["mb.stretch"].target = hover_ob
            else:
                new_atom.draw_type = 'SOLID'
                if new_bond:
                    new_bond.constraints["mb.stretch"].target = new_atom
                    if not event.alt:
                        self.coord_3d = mb_geometry.get_fixed_geometry(
                            context, first_atom, new_atom, self.coord_3d,
                            self.geometry)
                    
                    if event.ctrl:
                         # constrain length
                        self.coord_3d = mb_geometry.get_fixed_length(
                            context, first_atom, new_atom, self.coord_3d,
                            length=-1)
            
            new_atom.location = self.coord_3d - self.parent_coord_3d
            # sometimes, when bond is exactly along axis, the dimension goes
            # to zero due to the stretch constraint
            # check for this case and fix it
            if new_bond:
                mb_utils.check_ob_dimensions(new_bond)
        
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            # check if new atom is above already existing atom
            new_atom = context.object
            first_atom = context.scene.objects.get(self.first_atom_name)
            hover_ob = mb_utils.return_cursor_object(
                context, event, exclude=[new_atom], mb_type='ATOM')
            if hover_ob:
                mol = new_atom.mb.get_molecule()
                mol.remove_object(new_atom)
                mol.atom_index -= 1
                context.scene.objects.unlink(new_atom)
                bpy.data.objects.remove(new_atom)
                new_bond = context.scene.objects.get(self.new_bond_name)
                if new_bond:
                    mol.remove_object(new_bond)
                    context.scene.objects.unlink(new_bond)
                    bpy.data.objects.remove(new_bond)
                    mb_utils.add_bond(context, first_atom, hover_ob)
            return {'FINISHED'}
        
        return {'RUNNING_MODAL'}
    
    def invoke(self, context, event):
        self.element = context.scene.mb.globals.element_to_add
        self.geometry = context.scene.mb.globals.geometry_to_add
        hover_ob = mb_utils.return_cursor_object(context, event,
                                                 mb_type='ATOM')
        
        self.region, self.rv3d = mb_utils.get_region_data(
            context, event.mouse_x, event.mouse_y
            )
        
        if hover_ob:
            self.first_atom_name = hover_ob.name
            molecule = hover_ob.mb.get_molecule()
            self.parent_coord_3d = molecule.objects.parent.location
            self.depth_location = hover_ob.location + self.parent_coord_3d
        else:
            self.first_atom_name = " "
            self.depth_location = context.scene.cursor_location.copy()
        mouse_2d = event.mouse_x, event.mouse_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(
            context, mouse_2d, depth=self.depth_location)
        self.coord_3d -= self.parent_coord_3d
        
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

        if first_atom and first_atom.mb.type != 'ATOM':
            self.first_atom_name = " "
            first_atom = None
        
        if self.first_atom_name.strip() and not first_atom:
            debug_print('Object "{}" not found.'.format(self.first_atom_name),
                        level=1)
            return {'CANCELLED'}
        
        if first_atom:
            molecule = first_atom.mb.get_molecule()
        else:
            molecule = context.scene.mb.new_molecule()
        
        # create a new atom object with the molecule's properties
        new_atom = mb_utils.add_atom(context, self.coord_3d, self.element,
                                     self.element, molecule)
        self.new_atom_name = new_atom.name
        
        # add a bond if atom is added to existing molecule
        if first_atom:
            new_bond = mb_utils.add_bond(context, first_atom, new_atom)
            self.new_bond_name = new_bond.name
            new_bond.select = True
        
        context.scene.objects.active = new_atom
        new_atom.select = True
        return {'FINISHED'}


class MB_OT_select_bonded(Operator):
    '''
    Select connected molecule based on mb data
    '''
    bl_idname = "mb.select_bonded"
    bl_description = "Select connected molecule based on atom bonds"
    bl_label = "Select bonded"
    bl_options = {'UNDO', 'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        return context.object and context.object.mb.type in ('ATOM', 'BOND')
    
    def execute(self, context):
        # recursive functions
        def atom(ob):
            for b in ob.mb.bonds:
                if b not in objects:
                    objects.append(b)
                    b.select = True
                    bond(b)
            return {'FINISHED'}
        
        def bond(ob):
            for a in ob.mb.bonded_atoms:
                if a not in objects:
                    objects.append(a)
                    a.select = True
                    atom(a)
            return {'FINISHED'}
        
        objects = []
        
        for ob in context.selected_objects:
            if ob.mb.type in ('ATOM', 'BOND'):
                if ob.mb.type == 'ATOM':
                    return atom(ob)
                elif ob.mb.type == 'BOND':
                    return bond(ob)
            else:
               debug_print('mb.type {} not compatible'.format(ob.mb.type),
                           level=2)


class MB_OT_select_molecule(Operator):
    '''
    Select all objects belonging to molecule based on mb data
    '''
    bl_idname = "mb.select_molecule"
    bl_description = "Select all objects belonging to molecule"
    bl_label = "Select molecule"
    bl_options = {'UNDO', 'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        return context.object
    
    def execute(self, context):
        mol = context.object.mb.get_molecule()
        if not mol:
            self.report("Active object is not part of a molecule.")
            return {'CANCELLED'}
        
        for ob in mol.objects.get_all_objects():
            ob.select = True
        return {'FINISHED'}


class MB_OT_combine_molecules(Operator):
    '''
    Select all objects belonging to molecule based on mb data
    '''
    bl_idname = "mb.combine_molecules"
    bl_description = "Combine selected molecules into one."
    bl_label = "Combine molecules"
    bl_options = {'UNDO', 'REGISTER'}
    
    molecule_id = EnumProperty(name="Molecules", items=get_molecules)
    
    @classmethod
    def poll(cls, context):
        return context.object and context.object.mb.get_molecule()
    
    def draw(self, context):
        mol = context.scene.mb.molecules.get(self.molecule_id)
        layout = self.layout
        label = "Do you want to add all selected atoms and bonds to"
        label += " molecule {} ({})?".format(mol.name_mol, mol.name)
        layout.label(label)
    
    def invoke(self, context, event):
        mol = context.object.mb.get_molecule()
        self.molecule_id = mol.name
        return context.window_manager.invoke_props_dialog(self)
    
    def execute(self, context):
        mol = context.scene.mb.molecules.get(self.molecule_id)
        
        for ob in context.selected_objects:
            if ob.mb.type in ("ATOM", "BOND"):
                if not ob.mb.get_molecule() == mol:
                    ob.mb.get_molecule().remove_object(ob)
                    mol.add_object(ob)
        
        return {'FINISHED'}


class MB_OT_center_mol_parent(Operator):
    '''
    Set molecule parent into center of mass of atoms
    '''
    bl_idname = "mb.center_mol_parent"
    bl_description = "Put parent to center of mass of molecule"
    bl_label = "Parent to CoM"
    bl_options = {'UNDO', 'REGISTER'}
    
    molecule_id = EnumProperty(name="Molecules", items=get_molecules)
    
    @classmethod
    def poll(cls, context):
        return context.object
    
    def execute(self, context):
        molecule = context.object.mb.get_molecule()
        if not molecule:
        
            debug_print(
                "ERROR in mb.center_mol_parent: Molecule "
                "{} not found in scene.".format(self.molecule_id),
                level=0)
            return {'CANCELLED'}
        origin = Vector((0.0,0.0,0.0))
        atoms = molecule.objects.atoms
        locs = [atom.location for atom in atoms]
        center = sum(locs, origin) / len(molecule.objects.atoms)
        for atom in molecule.objects.atoms:
            atom.location -= center
        molecule.objects.parent.location = center
        return {'FINISHED'}
        

def get_molecules(self, context):
        lst = []
        for mol in context.scene.mb.molecules:
            lst.append((mol.name, mol.name_mol, ""))
        return lst


class MB_OT_draw_unit_cell(Operator):
    bl_idname = "mb.draw_unit_cell"
    bl_label = "Draw unit cell of structure"
    bl_options = {'REGISTER', 'UNDO'}
    
    molecule_id = EnumProperty(name="Molecules", items=get_molecules)
    
    @classmethod
    def poll(self, context):
        return context.object and context.object.mb.get_molecule()
    
    def invoke(self, context, event):
        mol = context.object.mb.get_molecule()
        self.molecule_id = mol.name
        return self.execute(context)
    
    def execute(self, context):
        mol = context.scene.mb.molecules.get(self.molecule_id)
        if not "unit_cells" in mol or not mol["unit_cells"]:
            mol["unit_cells"] = [[Vector((5,0,0)),
                                 Vector((0,5,0)),
                                 Vector((0,0,5))]]
        obs = mb_utils.draw_unit_cell(mol, context)
        if obs:
            for ob in obs:
                ob.select = True
            context.scene.objects.active = obs[0]
        
        return {'FINISHED'}

class MB_OT_draw_dipole(Operator):
    bl_idname = "mb.draw_dipole"
    bl_label = "Draw dipole of molecule"
    bl_options = {'REGISTER', 'UNDO'}
    
    dipole_vec = FloatVectorProperty(name="Dipole vector", size=3)
    molecule_id = EnumProperty(name="Molecules", items=get_molecules)
    
    def draw(self, context):
        layout = self.layout
        col = layout.column()
        col.prop(self, "dipole_vec")
        
    def invoke(self, context, event):
        if context.object and context.object.mb.get_molecule():
            self.molecule_id = context.object.mb.get_molecule().name
        return context.window_manager.invoke_props_dialog(self)
    
    def execute(self, context):
        mol = context.scene.mb.molecules.get(self.molecule_id)
        obs = mb_utils.draw_dipole(mol, self.dipole_vec, context)
        for ob in obs:
            ob.select = True
        context.scene.objects.active = obs[0]
        return {'FINISHED'}


class MB_OT_remove_dipole(Operator):
    bl_idname = "mb.remove_dipole"
    bl_label = "Remove dipole"
    bl_options = {'REGISTER', 'UNDO'}
    
    @classmethod
    def poll(self, context):
        return context.object and context.object.mb.get_molecule()
    
    def execute(self, context):
        mol = context.object.mb.get_molecule()
        mb_utils.remove_dipole(mol, context)
        if not context.object:
            context.scene.objects.active = mol.objects.parent
        return {'FINISHED'}


class MB_OT_remove_unit_cell(Operator):
    bl_idname = "mb.remove_unit_cell"
    bl_label = "Remove unit cell"
    bl_options = {'REGISTER', 'UNDO'}
    
    @classmethod
    def poll(self, context):
        return context.object and context.object.mb.get_molecule()
    
    def execute(self, context):
        mol = context.object.mb.get_molecule()
        mb_utils.remove_unit_cell(mol, context)
        if not context.object:
            context.scene.objects.active = mol.objects.parent
        return {'FINISHED'}


class MD_OT_import_modes(bpy.types.Operator):
    bl_idname = "mb.import_modes"
    bl_label = "Import vibrational modes for molecule"
    bl_options = {'REGISTER', 'UNDO'}
    
    filepath = StringProperty(
        name="File Path", description="Filepath used for importing one file",
        maxlen=1024, subtype='FILE_PATH')
    
    file_format = EnumProperty(
        name="File format", description="Choose file format of mode file",
        items=mb_utils.enums.mode_file_format, default='XYZ')
    
    def draw(self, context):
        layout = self.layout
        layout.prop(self, "file_format")
    
    @classmethod
    def poll(cls, context):
        return context.object and context.object.mb.get_molecule()
    
    def invoke(self, context, event):
        molecule = context.object.mb.get_molecule()
        if molecule.objects.atoms[0].animation_data:
            # TODO allow user a choice
            self.report({'WARNING'}, "Atoms already contain animation data. Will overwrite")
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
    def execute(self, context):
        ret = mb_import_export.import_modes(
            context,
            self.report,
            self.filepath,
            self.file_format,
            context.object.mb.get_molecule()
            )
        if ret:
            return {'FINISHED'}
        else:
            return {'CANCELLED'}
        

class MD_OT_import_molecules(bpy.types.Operator):
    bl_idname = "mb.import_molecules"
    bl_label = "Import files"
    __doc__ = ""

    directory = StringProperty(
        name="Directory", description="Directory used for importing the file",
        maxlen=1024, subtype='FILE_PATH')
    filepath = StringProperty(
        name="File Path", description="Filepath used for importing one file",
        maxlen=1024, subtype='FILE_PATH')
    files = CollectionProperty(
        name="File Path",
        description="List with file names used for importing",
        type=bpy.types.OperatorFileListElement)
    
    #--- molecule properties -------------------------------------------------#
    name_mol = StringProperty(
        name="Molecule Name", description="Name of imported molecule",
        default="Molecule") # human readable name
    bond_material = EnumProperty(
        name="Bond material", description="Choose bond material",
        items=mb_utils.enums.bond_material, default='ATOMS')
    bond_color = FloatVectorProperty(
        name='Bond color',
        default=(0.8, 0.8, 0.8), subtype='COLOR')
    draw_style = EnumProperty(
        name="Display style", description="Style to draw atoms and bonds",
        items=mb_utils.enums.molecule_styles, default='BAS')
    radius_type = EnumProperty(
        name="Radius type",
        items=mb_utils.enums.radius_types, default='covalent')
    bond_radius = FloatProperty(
        name="Bond radius",
        description="Radius of bonds for Sticks, and Ball and Sticks",
        default=0.1, min=0.0, max=3.0)
    
    # this is a duplicate class from mb_datastructure for
    class atom_scale(PropertyGroup):
        name = StringProperty()
        val = FloatProperty(name="Atom scale", default=0.4, min=0.0, max=5.0,
                            precision=2)
    
    atom_scales = CollectionProperty(type=atom_scale)
    refine_atoms = IntProperty(
        name="Refine atoms", description="Refine value for atom meshes",
        default=8, min=2, max=64)
    refine_bonds = IntProperty(
        name="Refine bonds", description="Refine value for atom meshes",
        default=8, min=2, max=64)
    bond_type = EnumProperty(
        name="Bond type", description="Select how bonds should behave",
        items=mb_utils.enums.bond_types, default='CONSTRAINT')
    # TODO this might be handy for different units in files
    #scale_distances = FloatProperty (
        #name = "Distances", default=1.0, min=0.0001,
        #description = "Scale factor for all distances")
    length_unit = EnumProperty(
        name="Unit",
        description="Unit in input file, will be converted to Angstrom",
        items=mb_utils.enums.angstrom_per_unit, default='1.0')
    length_unit_other = FloatProperty(
        name="Custom Unit",
        description="Enter conversion factor in Angstrom/unit in file",
        default=1.0, min=0.000000001)
    bond_guess = BoolProperty(
       name="Guess bonds", description="Guess bonds that are not in the file.",
       default=True)
    use_mask = StringProperty(
        name="Masking object",
        description="Select object that sets boundaries to imported strucure.")
    mask_flip = BoolProperty(
        name="Mask flip",
        description="Invert masking effect (only atoms outside of mask "
                    "object are imported).")
    draw_unit_cell = BoolProperty(
       name="Draw unit cell", description="Draw the unit cell if applicable.",
       default=False)
    supercell = IntVectorProperty(
        name="Supercell", description="Specify supercell dimensions",
        size=3, default=(1,1,1), min=1, subtype='XYZ')
    put_origin = BoolProperty(
        name="Parent to origin",
        description="Put the parent object into the global origin",
        default=False)
    parent_center = BoolProperty(
        name="Center parent",
        description="Put the parent object into the center of mass",
        default=False)
    
    @classmethod
    def poll(cls, context):
        return context.scene.mb.is_initialized
    
    def invoke(self, context, event):
        if not len(self.atom_scales):
            for style in ('BALLS', 'BAS', 'STICKS'):
                atom_scale = self.atom_scales.add()
                atom_scale.name = style
                val = context.scene.mb.globals.atom_scales[style].val
                atom_scale.val = val
        
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def execute(self, context):
        if not len(self.atom_scales):
            for style in ('BALLS', 'BAS', 'STICKS'):
                atom_scale = self.atom_scales.add()
                atom_scale.name = style
                val = context.scene.mb.globals.atom_scales[style].val
                atom_scale.val = val
        
        files = set([os.path.join(self.directory, f.name) for f in self.files])
        files.add(bpy.path.abspath(self.filepath))
        
        for filepath in files:
            if not os.path.exists(filepath):
                debug_print("ERROR: {} not found".format(filepath), level=0)
                return {'CANCELLED'}
        
            new_molecule = context.scene.mb.new_molecule(
                                name_mol=self.name_mol,
                                draw_style=self.draw_style,
                                radius_type=self.radius_type,
                                bond_radius=self.bond_radius,
                                refine_atoms=self.refine_atoms,
                                refine_bonds=self.refine_bonds,
                                atom_scales=self.atom_scales)
            
            ## check if select_frames is ok, otherwise import first frame only
            error_list = []
            
            mask = bpy.data.objects.get(self.use_mask)
            mask_planes = []
            if not mask and self.use_mask:
                error_list.append('Mask object not found. Not using mask.')
            elif mask:
                world_mat = mask.matrix_world
                # only rotate normal vectors
                rot = world_mat.to_3x3().normalized()
                # get all faces (normal vector, point on plane) from mask object
                mask_planes = [(rot*pg.normal.copy(),
                                world_mat*mask.data.vertices[pg.vertices[0]].co)
                                for pg in mask.data.polygons]
            if error_list:
                debug_print('\n'.join(error_list), level=1)
            
            if self.length_unit == 'OTHER':
                scale_distances = self.length_unit_other
            else:
                scale_distances = float(self.length_unit)
            # Execute main routine
            try:
                worked = mb_import_export.import_molecule(
                    context,
                    self.report,
                    filepath,
                    new_molecule,
                    self.refine_atoms,
                    self.refine_bonds,
                    self.bond_material,
                    self.bond_type,
                    scale_distances,
                    self.bond_guess,
                    self.put_origin,
                    self.parent_center,
                    mask_planes,
                    self.mask_flip,
                    self.draw_unit_cell,
                    self.supercell,
                    )
            except:
                worked = False
                raise
            finally:
                if not worked:
                    context.scene.mb.remove_molecule(new_molecule)
        return {'FINISHED'}



    def draw(self, context):
        layout = self.layout
        layout.label("Molecule name")
        layout.prop(self, "name_mol", text="")
        
        layout.separator()
        layout.label("Units")
        layout.prop(self, "length_unit", text="")
        row = layout.row()
        row.active = (self.length_unit == 'OTHER')
        row.prop(self, "length_unit_other")
        
        layout.separator()
        layout.label("Periodic systems")
        layout.prop(self, "draw_unit_cell")
        layout.prop(self, "supercell")
        
        layout.separator()
        layout.label("Parent location")
        row = layout.row()
        row.prop(self, "put_origin")
        row.prop(self, "parent_center")
        
        layout.separator()
        layout.label(text="Masking object")
        row = layout.row()
        row.prop_search(self, "use_mask", bpy.data, "objects", text="")
        row.prop(self, "mask_flip")
        
        layout.separator()
        layout.label("Draw style")
        layout.prop(self, "draw_style", text="")
        
        split = layout.split()
        col = split.column()
        col.prop(self, "refine_atoms")
        col.prop(self.atom_scales[self.draw_style], "val", text="Atom scaling")
        col.label(text="Atom radius")
        col.prop(self, "radius_type", text="")
        
        col = split.column()
        col.prop(self, "refine_bonds")
        col.prop(self, "bond_radius")
        col.label(text="Bond material")
        col.prop(self, "bond_material", text="")
        col.prop(self, "bond_color")
        col.prop(self, "bond_guess")
        col.prop(self, "bond_type")
        

class MB_OT_draw_bond_lengths(bpy.types.Operator):
    """Draw a line with the mouse"""
    bl_idname = "mb.show_bond_lengths"
    bl_label = "Show bond lengths"
    
    def modal(self, context, event):
        context.area.tag_redraw()
        if not context.scene.mb.globals.show_bond_lengths:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'CANCELLED'}
        return {'PASS_THROUGH'}
    
    def execute(self, context):
        if context.area.type == 'VIEW_3D':
            # the arguments we pass to the callback
            args = (self, context)
            # Add the region OpenGL drawing callback
            # draw in view space with 'POST_VIEW' and 'PRE_VIEW'
            self._handle = bpy.types.SpaceView3D.draw_handler_add(
                mb_utils.callback_draw_length, args, 'WINDOW', 'POST_PIXEL'
                )

            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "View3D not found, cannot run operator")
            return {'CANCELLED'}