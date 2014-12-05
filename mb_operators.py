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


class MB_OT_initialize(Operator):
    bl_idname = 'mb.initialize'
    bl_label = 'Initialize MolBlend'
    bl_description = 'Make sure that drivers work and elements are loaded.'
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
        return {'FINISHED'}


class MB_OT_modal_add(Operator):
    bl_idname = 'mb.modal_add'
    bl_label = 'test MolBlend'
    bl_options = {'REGISTER'}
    
    @classmethod
    def poll(cls, context):
        return mb_utils.is_initialized(context)
    
    def modal(self, context, event):
        #print("modal")
        if bpy.ops.object.select_all.poll():
            bpy.ops.object.select_all(action='DESELECT')
            hover_ob = mb_utils.return_cursor_object(context, event, 
                                                     mb_type='ATOM')
            if hover_ob is not None:
                hover_ob.select = True
            context.scene.objects.active = hover_ob
        if (event.type == 'LEFTMOUSE' and event.value == 'PRESS' and
                context.area.type == 'VIEW_3D'):
            # TODO make sure we're in the 3D View window
            bpy.ops.mb.add_atom('INVOKE_DEFAULT',
                                shift=event.shift,
                                ctrl=event.ctrl,
                                alt=event.alt)
            return {'RUNNING_MODAL'}
        elif event.type in ('RIGHTMOUSE', 'ESC'):
            context.window.cursor_modal_restore()
            return {'CANCELLED'}
        return {'PASS_THROUGH'}
    
    def invoke(self, context, event):
        if context.area.type == 'VIEW_3D':
            context.window.cursor_modal_set("CROSSHAIR")
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            return {'CANCELLED'}


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
        items=mb_utils.enums.geometries, default='SINGLE')
    
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
        mouse_2d = event.mouse_region_x, event.mouse_region_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(
            context, mouse_2d, depth=self.depth_location)
        
        if event.type == 'MOUSEMOVE':
            new_atom = context.scene.objects.get(self.new_atom_name)
            context.scene.objects.active = new_atom
            new_atom.select = True
            new_bond = context.scene.objects.get(self.new_bond_name)
            first_atom = context.scene.objects.get(self.first_atom_name)
            
            if not event.ctrl and not event.alt:
                hover_ob = mb_utils.return_cursor_object(
                    context, event, exclude=[new_atom], mb_type='ATOM')
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
                    self.coord_3d = mb_utils.get_fixed_geometry(
                        context, first_atom, new_atom, self.coord_3d,
                        self.geometry)
                if event.alt and new_bond:
                    # constrain length
                    length = 1.0
                    self.coord_3d = mb_utils.get_fixed_length(
                        context, first_atom, new_atom, self.coord_3d, 
                        length=-1)
            
            new_atom.location = self.coord_3d
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
                if new_bond and hover_ob:
                    context.scene.mb.remove_object(new_bond)
                    mb_utils.add_bond(context, first_atom, hover_ob)
            return {'FINISHED'}
        
        return {'RUNNING_MODAL'}
    
    def invoke(self, context, event):
        self.element = context.window_manager.mb.globals.element_to_add
        self.geometry = context.window_manager.mb.globals.geometry_to_add
        hover_ob = mb_utils.return_cursor_object(context, event, 
                                                 mb_type='ATOM')
        
        if hover_ob:
            self.first_atom_name = hover_ob.name
            self.depth_location = hover_ob.location
        else:
            self.first_atom_name = " "
            self.depth_location = context.scene.cursor_location.copy()
        mouse_2d = event.mouse_region_x, event.mouse_region_y
        self.coord_3d = mb_utils.mouse_2d_to_location_3d(
            context, mouse_2d, depth=self.depth_location)
        
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
        # so throw a warning and stop operator
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
        
        context.scene.objects.active = new_atom
        new_atom.select = True
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
                debug_print(
                    "ERROR in mb.center_mol_parent: Molecule "
                    "{} not found in scene.".format(self.molecule_name),
                    level=0)
                return {'CANCELLED'}
            else:
                origin = Vector((0.0,0.0,0.0))
                atoms = molecule.objects.atoms
                locs = [atom.get_object().location for atom in atoms]
                center = sum(locs, origin) / len(molecule.objects.atoms)
                for atom in molecule.objects.atoms:
                    atom.get_object().location -= center
                molecule.objects.parent.get_object().location = center
        return {'FINISHED'}


class MB_OT_draw_dipole(Operator):
    bl_idname = "mb.draw_dipole"
    bl_label = "Draw dipole of molecule"
    bl_options = {'REGISTER', 'UNDO'}
    
    dipole_vec = FloatVectorProperty(name="Dipole vector", size=3)
    molecule_id = StringProperty(
        name="Molecule identifier", 
        update=mb_utils.update_molecule_selection)
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop_search(self, "molecule_id", context.scene.mb, 
                        "molecules", text="")
        col = layout.column()
        col.prop(self, "dipole_vec")
        
    def invoke(self, context, event):
        if context.object and context.object.mb.get_molecule():
            self.molecule_id = context.object.mb.get_molecule().name
        return context.window_manager.invoke_props_dialog(self)
    
    def execute(self, context):
        mol = context.scene.mb.molecules.get(self.molecule_id)
        if not mol:
            debug_print(
                "ERROR: draw_dipole: {} not found.".format(self.molecule_id),
                level=0)
            return {'CANCELLED'}
        
        # add empty as stretch target
        dipole_ob = bpy.data.objects.new(
            "{}_dipole_target".format(mol.name_mol), None)
        dipole_ob.empty_draw_type = 'SINGLE_ARROW'
        dipole_ob.empty_draw_size = 0.5
        bpy.context.scene.objects.link(dipole_ob)
        mol.objects.dipole.name = dipole_ob.name
        dipole_ob.location = self.dipole_vec
        dipole_ob.parent = mol.objects.parent.get_object()
        dipole_ob.mb.molecule_name = mol.name
        
        # add arrow object
        arrow_mesh = mb_utils.get_arrow_data()
        arrow_ob = bpy.data.objects.new("{}_dipole".format(mol.name_mol), 
                                        arrow_mesh)
        arrow_ob.parent = mol.objects.parent.get_object()
        bpy.context.scene.objects.link(arrow_ob)
        bpy.context.scene.objects.active = arrow_ob
        arrow_ob.mb.molecule_name = mol.name

        c = arrow_ob.constraints.new('STRETCH_TO')
        c.name = "stretch"
        c.rest_length = 1.0
        c.volume = 'NO_VOLUME'
        c.target = dipole_ob
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
        bpy.ops.object.select_all(action='DESELECT')
        hover_ob = mb_utils.return_cursor_object(context, event)
        if hover_ob is not None:
            hover_ob.select = True
        context.scene.objects.active = hover_ob
        if event.type in ['ESC', 'LEFTMOUSE']:
            return {'FINISHED'}
        
        return {'RUNNING_MODAL'}
    
    def invoke(self, context, event):
        return self.execute(context)
    
    def execute(self, context):
        if context.space_data.type == 'VIEW_3D':
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            return {'CANCELLED'}


class MB_OT_import_molecule(Operator):
    bl_idname = "mb.import_molecule"
    bl_label  = "Import XYZ/PDB (*.xyz,*.pdb)"
    bl_options = {'PRESET', 'UNDO'}
    
    filename_ext = "*.pdb;*.xyz"
    filepath = StringProperty(
        name="File Path", description="Filepath used for importing one file",
        maxlen=1024, subtype='FILE_PATH')
    directory = StringProperty(
        name="Directory", description="Directory used for importing the file",
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
        default=8, min=3, max=64)
    refine_bonds = IntProperty(
        name="Refine bonds", description="Refine value for atom meshes",
        default=8, min=3, max=64)
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
        default=1.0, min=0.000001)
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
    use_center = BoolProperty(
        name="Object to origin (first frame)",
        description="Put the object into the global origin, "
                    "the first frame only",
        default=False)    
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, "name_mol")
        
        layout.separator()
        row = layout.row()
        row.prop(self, "draw_style")
        row = layout.row()
        # Atom props
        col = row.column()
        col.prop(self, "refine_atoms")
        col.prop(self.atom_scales[self.draw_style], "val", text="Atom scaling")
        col.label(text="Atom radius")
        col.prop(self, "radius_type", text="")
        
        col = row.column()        
        col.prop(self, "refine_bonds")
        row = layout.row()
        col.prop(self, "bond_radius")
        col.label(text="Bond material")
        col.prop(self, "bond_material", text="")
        col.prop(self, "bond_color")
        col.prop(self, "bond_guess")
        
        layout.separator()
        row = layout.row()
        row.label(text="Masking object")
        row.prop_search(self, "use_mask", bpy.data, "objects", text="")
        row = layout.row()
        row.prop(self, "mask_flip")
        
        layout.separator()
        row = layout.row()
        row.prop(self, "use_center")
        
        row = layout.row()
        row.prop(self, "length_unit")
        row = layout.row()
        row.active = (self.length_unit == 'OTHER')
        row.prop(self, "length_unit_other")
        row = layout.row()
        row.prop(self, "draw_unit_cell")
        row = layout.row()
        row.prop(self, "supercell")
    
    def invoke(self, context, event):
        # before import dialog is opened, initialize atom scales
        if not len(self.atom_scales):
            if not len(context.scene.mb.globals.atom_scales):
                default_scales = {'BALLS': 1.0, 'BAS': 0.4, 'STICKS': 0.001}
            else:
                default_scales = {}
                for style in ('BALLS', 'BAS', 'STICKS'):
                    val = context.scene.mb.globals.atom_scales[style].val
                    default_scales[style] = val
            for style in default_scales:
                atom_scale = self.atom_scales.add()
                atom_scale.name = style
                atom_scale.val = default_scales[style]
        
        return context.window_manager.invoke_props_dialog(self)
    
    def execute(self, context):
        if not len(self.atom_scales):
            if not len(context.scene.mb.globals.atom_scales):
                default_scales = {'BALLS': 1.0, 'BAS': 0.4, 'STICKS': 0.001}
            else:
                default_scales = {}
                for style in ('BALLS', 'BAS', 'STICKS'):
                    val = context.scene.mb.globals.atom_scales[style].val
                    default_scales[style] = val
            for style in default_scales:
                atom_scale = self.atom_scales.add()
                atom_scale.name = style
                atom_scale.val = default_scales[style]
        
        import_props = context.scene.mb.globals.import_props
        filepath = bpy.path.abspath(import_props.filepath)
        if not os.path.exists(filepath):
            debug_print("ERROR: {} not found".format(filepath), level=0)
            return {'CANCELLED'}
        if import_props.modes:
            modes_path = import_props.modes_path
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
        
        mask = bpy.data.objects.get(self.use_mask)
        mask_planes = []
        if not mask and self.use_mask:
            error_list.append('Mask object not found. Not using mask.')
        elif mask:
            world_mat = mask.matrix_world
            # get all faces (normal vector, point on plane) from mask object
            mask_planes = [(world_mat*pg.normal.copy(), 
                            world_mat*mask.data.vertices[pg.vertices[0]].co) 
                            for pg in mask.data.polygons]
        
        if error_list:
            debug_print('\n'.join(error_list), level=1)
        
        if self.length_unit == 'OTHER':
            scale_distances = self.length_unit_other
        else:
            scale_distances = float(self.length_unit)
        # Execute main routine
        worked = mb_import_export.import_molecule(
                                filepath,
                                modes_path,
                                new_molecule,
                                self.refine_atoms,
                                self.refine_bonds,
                                scale_distances,
                                self.bond_guess,
                                self.use_center,
                                mask_planes,
                                self.mask_flip,
                                self.draw_unit_cell,
                                self.supercell,
                                )
        if not worked:
            context.scene.mb.remove_molecule(new_molecule)
        return {'FINISHED'}

#class ExportMolecule(Operator):
    #bl_idname = "mb.export_molecule"
    #bl_label  = "Export"
    ##filename_ext = "*.pdb;*.xyz"
    ##filter_glob  = StringProperty(default=filename_ext, options={'HIDDEN'},)
    #bl_description = "Export molecule structure to file (*.pdb, *.xyz)"
    
    ##atom_pdb_export_type = EnumProperty(
        ##name="Type of Objects",
        ##description="Choose type of objects",
        ##items=(('0', "All", "Export all active objects"),
               ##('1', "Elements", "Export only those active objects which have"
                                 ##" a proper element name")),
               ##default='1',) 
    ##file_type = EnumProperty(
        ##name="File type", default="XYZ", items=mb_utils.enums.file_types,
        ##description="File format to export to")
    ##length_unit = EnumProperty(
        ##name="Unit", default='1.0', items=mb_utils.enums.angstrom_per_unit,
        ##description="Unit in output file (to convert to from Angstrom)")
    ##length_unit_other = FloatProperty(
        ##name="Custom Unit", default=1.0, min=0.000001,
        ##description="Enter unit of export file as Angstrom/unit")
    ##selection_only = BoolProperty(name="Selected Objects", default=True,
        ##description="Only export selected objects")
    
    ##def draw(self, context):
        ##layout = self.layout
        ##row = layout.row()
        ##row.prop(self, "file_type")
        ##row = layout.row()
        ##row.prop(self, "selection_only")
        ##row = layout.row()
        ##row.prop(self, "length_unit")
        ##row = layout.row()
        ##row.active = (self.length_unit == 'OTHER')
        ##row.prop(self, "length_unit")
    
    ##def invoke(self, context, event):
        ##active_ob = context.object
        ##if active_ob and hasattr(active_ob, 'mb') and active_ob.mb.get_molecule():
            ##self.filepath = active_ob.mb.get_molecule().name_mol
        
    
    #def execute(self, context):
        #export_props = context.scene.mb.globals.export_props
        #if export_props.length_unit == 'OTHER':
            #scale_distances = export_props.length_unit_other
        #else:
            #scale_distances = float(export_props.length_unit)
        
        ## convert into relative path if only filename is given (otherwise writes to homedirectory)
        #filepath = export_props.filepath
        #if not os.path.isabs(filepath) and not filepath.startswith("//"):
            #filepath = "//" + filepath
        
        #if export_props.file_type == 'PDB':
            #if not export_props.filepath[-4:] == ".pdb":
                #export_props.filepath += ".pdb"
            #mb_import_export.export_pdb(bpy.path.abspath(export_props.filepath),
                                        #export_props.selection_only,
                                        #scale_distances)
        #elif export_props.file_type == 'XYZ':
            #if not export_props.filepath[-4:] == ".xyz":
                #export_props.filepath += ".xyz"
            #mb_import_export.export_xyz(bpy.path.abspath(export_props.filepath),
                                        #export_props.selection_only,
                                        #scale_distances)
        #elif export_props.file_type == 'B4W':
            #if not export_props.filepath[-4:] == ".xyz":
                #export_props.filepath += ".xyz"
            #mb_import_export.export_b4w(bpy.path.abspath(export_props.filepath))
        #return {'FINISHED'}
