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

if "bpy" in locals():
    import imp
    imp.reload(mb_utils)
else:
    from molblend import mb_utils

import bpy
import os

from bpy.types import PropertyGroup

from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolVectorProperty,
                       PointerProperty,
                       CollectionProperty,
                       EnumProperty)

from .helper import debug_print


class mb_object_pointer(PropertyGroup):
    name = StringProperty(name="Object name")
    
    def get_object(self):
        return bpy.data.objects.get(self.name)


class mb_mesh_pointer(PropertyGroup):
    name = StringProperty(name="Mesh name")
    
    def get_data(self):
        return bpy.data.meshes.get(self.name)


class atom_scale(PropertyGroup):
    name = StringProperty()
    val = FloatProperty(name="Atom scale", default=0.4, min=0.0, max=5.0,
                        precision=2, update=mb_utils.update_all_meshes)


class mb_mesh(PropertyGroup):
    type = EnumProperty(name="type", items=[('MESH', 'MESH', 'MESH'),],
                        description="Do not change", default='MESH')


class mb_atom_mode(PropertyGroup):
    name = IntProperty(name="index")
    index = IntProperty(name="index")
    freq = FloatProperty(name="frequency")
    vec = FloatVectorProperty(name="vector", subtype="XYZ")


class mb_object(PropertyGroup):
    index = IntProperty(name="Index")
    name = StringProperty(name="Object name")
    type = EnumProperty(
        name="type", description="Select the object type",
        items=mb_utils.enums.object_types, default='NONE')
    molecule_name = StringProperty(name="Molecule identifier")
    
    # used by type == 'ATOM'
    bonds = CollectionProperty(type=mb_object_pointer)
    atom_name = StringProperty(name="Atom name")
    element = StringProperty(
        name="Element", description="Element Symbol",
        update=mb_utils.update_atom_element)
    element_long = StringProperty(
        name="Element name", description="Full element name")
    
    # used by type == 'BOND'
    bonded_atoms = CollectionProperty(type=mb_object_pointer)
    modes = CollectionProperty(type=mb_atom_mode)
    
    def get_molecule(self):
        return bpy.context.scene.mb.molecules.get(self.molecule_name)
    
    def get_object(self):
        return bpy.data.objects.get(self.name)
    
    def add_bond(self, ob, replace_name=""):
        if not self.type == 'ATOM':
            debug_print("WARNING: Something is trying to add bond to "
                        "non-ATOM type object", level=1)
            return
        
        ob_name = ob.name
        if ob_name in self.bonds:
            item = self.bonds[ob_name]
        elif replace_name:
            item = self.bonds.get(replace_name)
            if item:
                # if object is already in self.bonds, just delete replace_item
                if ob_name in self.bonds:
                    self.bonds.remove_bond(item.get_object())
                else:
                    item.name = ob_name
            else:
                debug_print("WARNING: {} not found in {}".format(
                            replace_name, self.bonds), level=3)
        else:
            item = self.bonds.add()
            item.name = ob_name
        return item
    
    def remove_bond(self, ob):
        if not self.type == 'ATOM':
            debug_print("WARNING: Something is trying to remove bond from "
                        "non-ATOM type object", level=1)
            return
        for i, b in enumerate(self.bonds):
            if b.name == ob.name:
                self.bonds.remove(i)
                return
    
    def add_bonded_atom(self, ob, replace_name=""):
        if not self.type == 'BOND':
            debug_print("WARNING: Something is trying to add bonded_atom to "
                        "non-BOND type object", level=1)
            return
        
        ob_name = ob.name
        if ob_name in self.bonded_atoms:
            item = self.bonded_atoms[ob_name]
        elif replace_name:
            item = self.bonded_atoms.get(replace_name)
            if item:
                # if object is in self.bonded_atoms, delete replace_item
                if ob_name in self.bonded_atoms:
                    self.bonded_atoms.remove_bonded_atom(item.get_object())
                else:
                    item.name = ob_name
            else:
                debug_print("WARNING: {} not found in {}".format(
                            replace_name, self.bonded_atoms), level=3)
        else:
            item = self.bonded_atoms.add()
            item.name = ob_name
        return item
    
    def remove_bonded_atom(self, ob):
        if not self.type == 'BOND':
            debug_print("WARNING: Something is trying to remove bonded_atom "
                        "{} ({}) from non-BOND type object {} ({})".format(
                        ob.name, ob.mb.type, self.name, self.type), level=1)
            return
        
        for i, a in enumerate(self.bonded_atoms):
            if a.name == ob.name:
                self.bonded_atoms.remove(i)
                return
    
    def draw_properties(self, context, layout, ob):
        element = context.scene.mb.elements[self.element]
        props = {
            "Element": [self, "element", 10],
            "Name": [self, "atom_name", 20],
            "Atom radius": [element, "covalent", 30],
            "Atom color": [ob.material_slots[0].material, "diffuse_color", 40],
            }
        for label, (data, prop, i) in sorted(props.items(), 
                                             key=lambda t: t[-1][-1]):
            row = layout.row()
            row.prop(data, prop, text=label)


class mb_molecule_objects(PropertyGroup):
    atoms = CollectionProperty(name="Atoms", type=mb_object_pointer)
    bonds = CollectionProperty(name="Bonds", type=mb_object_pointer)
    bond_curve = PointerProperty(name="Parent", type=mb_object_pointer)
    meshes = CollectionProperty(name="Molecule meshes", type=mb_mesh_pointer)
    parent = PointerProperty(name="Parent", type=mb_object_pointer)
    dipole = PointerProperty(name="Dipole", type=mb_object_pointer)
    other = CollectionProperty(name="Bonds", type=mb_object_pointer)


class mb_molecule(PropertyGroup):
    index = IntProperty(name="Molecule index")
    name = StringProperty(name="Molecule identifier")
    name_mol = StringProperty(name="Molecule Name")
    # index that increases with each added atom in the molecule, but doesn't 
    # decrease when atom is deleted. => Not an indicator of size of molecule! 
    # Only guarantees uniqueness for atom names
    atom_index = IntProperty()
    objects = PointerProperty(name="Molecule objects", 
        type=mb_molecule_objects)
    
    # display properties
    bond_material = EnumProperty(
        name="Bond material", description="Choose bond material",
        items=mb_utils.enums.bond_material, default='ATOMS',
        update=mb_utils.update_bond_material)
    bond_color = FloatVectorProperty(
        name='Bond color', description="Color for generic bond",
        default=(0.8, 0.8, 0.8), min=0.0, max=1.0, subtype='COLOR',
        update=mb_utils.update_all_meshes)
    draw_style = EnumProperty(
        name="Display style", description="Style to draw atoms and bonds",
        items=mb_utils.enums.molecule_styles, default='BAS',
        update=mb_utils.set_draw_style)
    radius_type = EnumProperty(
        name="Radius type", description="Type of radius to use as reference",
        items=mb_utils.enums.radius_types, default='covalent', 
        update=mb_utils.update_radius_type)
    bond_radius = FloatProperty(
        name="Bond radius", description="Radius for bond objects",
        default=0.1, min=0.0, max=3.0,
        update=mb_utils.update_all_meshes)
    atom_scales = CollectionProperty(type=atom_scale)
    refine_atoms = IntProperty(
        name="Refine atoms", description="Refine value for atom meshes", 
        default=8, min=3, max=64,
        update=mb_utils.update_refine_atoms)
    refine_bonds = IntProperty(
        name="Refine bonds", description="Refine value for atom meshes", 
        default=8, min=3, max=64,
        update=mb_utils.update_refine_bonds)
    
    # vibrational modes
    max_mode = IntProperty(
        name="Number of modes", 
        description="Number of vibrational modes of molecule",
        default=0, min=0,
        )
    active_mode = IntProperty(
        name="Active Mode",
        description="Active Mode to display. 0 = equilibrium position",
        default=0, min=0,
        update=mb_utils.update_active_mode)
    mode_scale = FloatProperty(
        name="Mode Scale", description="Scale of normal mode displacement",
        default=1.0, min=-10.0, max=10.0, 
        update=mb_utils.update_active_mode)
    
    def draw_properties(self, layout):
        #layout.label("Molecule properties")
        row = layout.row()
        row.prop(self, "name_mol", text="")
        row = layout.row()
        row.label("(id: '{}')".format(self.objects.parent.name))
        
        row = layout.row()
        row.active = bool(self.max_mode)
        row.prop(self, "active_mode")
        row = layout.row()
        row.active = bool(self.max_mode)
        row.prop(self, "mode_scale")
        
        if self.objects.parent.get_object():
            col = layout.column()
            col.prop(self.objects.parent.get_object(), "location", 
                     text="Center of mass")
        row = layout.row()
        row.operator("mb.center_mol_parent")
        
        if self.objects.dipole.get_object():
            col = layout.column()
            col.prop(self.objects.dipole.get_object(), "location", 
                     text="Dipole")
        else:
            row = layout.row()
            row.operator("mb.draw_dipole")
    
    def draw_styles(self, layout):
        #layout.label("Molecule draw style")
        row = layout.row()
        row.label(self.name_mol)
        row = layout.row()
        row.label("(id: {}".format(self.objects.parent.name))
        
        props = {
            "Atom scale": [self.atom_scales[self.draw_style], "val", 10],
            "Bond radius": [self, "bond_radius", 20],
            "Radius type": [self, "radius_type", 30],
            "Draw style": [self, "draw_style", 40],
            "Bond material": [self, "bond_material", 50],
            "Bond color": [self, "bond_color", 60],
            "Refine atoms": [self, "refine_atoms", 70],
            "Refine bonds": [self, "refine_bonds", 80],
            }
        for label, (data, prop, i) in sorted(props.items(), 
                                             key=lambda t: t[-1][-1]):
            row = layout.row()
            row.prop(data, prop)
        
    def add_object(self, ob, replace_name=""):
        '''
        add an object's name to the molecule's atoms collection
        if replace_name is given, first atom's name is replaced by the new one
        if name is already in collection, just return the atom
        '''
        collection = {
            'ATOM': self.objects.atoms,
            'BOND': self.objects.bonds,
            'MESH': self.objects.meshes,
            'NONE': self.objects.other
            }[ob.mb.type]
        ob_name = ob.name
        if ob_name in collection:
            item = collection[ob_name]
        elif replace_name:
            item = collection.get(replace_name)
            if item:
                # if object is already in collection, just delete replace_item
                if ob_name in collection:
                    collection.remove_object(item.get_object())
                else:
                    item.name = ob_name
            else:
                debug_print(
                    "WARNING: "
                    "{} not found in {}".format(replace_name, collection), 
                    level=3)
        else:
            item = collection.add()
            item.name = ob_name
        return item
    
    def remove_object(self, ob):
        collection = {
            'ATOM': self.objects.atoms,
            'BOND': self.objects.bonds,
            'MESH': self.objects.meshes,
            'NONE': self.objects.other
            }[ob.mb.type]
        for i, a in enumerate(collection):
            if a.name == ob.name:
                collection.remove(i)
                return
    
    def remove_objects(self, ob_list=None):
        ob_list = ob_list or []
        objects = {'ATOM': [ob.name for ob in ob_list if ob.mb.type == 'ATOM'],
                   'BOND': [ob.name for ob in ob_list if ob.mb.type == 'BOND'],
                   'MESH': [ob.name for ob in ob_list if ob.mb.type == 'MESH'],
                   'NONE': [ob.name for ob in ob_list if ob.mb.type == 'NONE']}
        
        collections = {'ATOM': self.objects.atoms,
                      'BOND': self.objects.bonds,
                      'MESH': self.objects.meshes,
                      'NONE': self.objects.other}
        for ob_type, ob_name_list in objects.items():
            # find all the indices of the objects to delete
            indeces = []
            for i, a in enumerate(collections[ob_type]):
                if a.name in ob_name_list:
                    indeces.append(i)
            # delete higher numbers first to not mess up the order of the coll.
            for i in reversed(indeces):
                collections[ob_type].remove(i)
        return


class mb_action(PropertyGroup):
    name = StringProperty(name="Action name")
    mode_vector = FloatVectorProperty(
        name="Mode vector",
        description="Original mode vector for atom as read from file")


class mb_element(PropertyGroup):
    name = StringProperty(name="Element")
    element = StringProperty(name="Element")
    element_name = StringProperty(name="Element name")
    atomic_number = IntProperty(name="Atomic number")
    color = FloatVectorProperty(name="Color", subtype='COLOR', size=3)
    covalent = FloatProperty(
        name="Covalent radius",
        min=0.0, max=5.0, 
        update=mb_utils.update_all_meshes)
    vdw = FloatProperty(
        name="vdW radius", 
        min=0.0, max=5.0, 
        update=mb_utils.update_all_meshes)
    constant = FloatProperty(
        name="Constant radius",
        min=0.0, max=5.0,
        update=mb_utils.update_all_meshes)


class mb_scn_import(PropertyGroup):
    filepath = StringProperty(
        name="Import file", 
        description="Filepath to molecule file to import (.xyz, .pdb)",
        default=os.getcwd() + "/", subtype="FILE_PATH")
    modes = BoolProperty(
        name="Modes", 
        description="Import normal modes of molecule as keyframes.",
        default=False)
    modes_path = StringProperty(
        name="Modes file",
        description="Filepath to modes file to import "
                    "(In Quantum Espresso: dynmat.out)", 
        default="", subtype="FILE_PATH")


class mb_scn_export(PropertyGroup):
    filepath = StringProperty(name="Export file", default="", 
        subtype="FILE_PATH",
        description="Filepath to exported molecule file (.xyz, .pdb)")
    selection_only = BoolProperty(name="Selected Objects", default=True,
        description="Only export selected objects")
    file_type = EnumProperty(
        name="File type", default="XYZ", items=mb_utils.enums.file_types,
        description="File format to export to", 
        update=mb_utils.update_export_file_type)
    length_unit = EnumProperty(
        name="Unit", default='1.0', items=mb_utils.enums.angstrom_per_unit,
        description="Unit in output file (to convert to from Angstrom)")
    length_unit_other = FloatProperty(
        name="Custom Unit", default=1.0, min=0.000001,
        description="Enter unit of export file as Angstrom/unit")


class mb_scn_globals(PropertyGroup):
    draw_style = EnumProperty(
        name="Draw style", description="Style to draw atoms and bonds",
        items=mb_utils.enums.molecule_styles, default='BAS',
        )
    radius_type = EnumProperty(
        name="Radius type", 
        items=mb_utils.enums.radius_types, default='covalent')
    bond_radius = FloatProperty(
        name="Bond radius",
        description="Radius of bonds for Sticks, and Ball and Sticks", 
        default=0.1, min=0.0, max=3.0)
    atom_scales = CollectionProperty(type=atom_scale)
    import_props = PointerProperty(type=mb_scn_import)
    export_props = PointerProperty(type=mb_scn_export)


class mb_scene(PropertyGroup):
    elements = CollectionProperty(type=mb_element)
    molecules = CollectionProperty(type=mb_molecule)
    # index that increases with each added molecule, but doesn't decrease when
    # molecule is deleted. Guarantees uniqueness for molecule names
    molecule_index = IntProperty()
    globals = PointerProperty(type=mb_scn_globals)
    
    def new_molecule(self):
        mol = self.molecules.add()
        
        mol.index = self.molecule_index
        self.molecule_index += 1
        mol.name = "molecule_{}".format(mol.index)
        mol.name_mol = "Molecule {}".format(mol.index)
        mol.draw_style = self.globals.draw_style
        mol.radius_type = self.globals.radius_type
        mol.bond_radius = self.globals.bond_radius
        for scale in self.globals.atom_scales:
            new_scale = mol.atom_scales.add()
            new_scale.name = scale.name
            new_scale.val = scale.val
        # create new empty that will be the parent for the molecule
        parent_ob = bpy.data.objects.new(mol.name_mol, None)
        parent_ob.empty_draw_type = 'SPHERE'
        parent_ob.empty_draw_size = 0.3
        bpy.context.scene.objects.link(parent_ob)
        mol.objects.parent.name = parent_ob.name
        parent_ob.mb.molecule_name = mol.name
        return mol
    
    def remove_molecule(self, mol):
        # Make sure all objects are deleted first
        if (0 == len(mol.objects.atoms) and
            0 == len(mol.objects.other) == len(mol.objects.bonds)):
            # delete all meshes
            for me in mol.objects.meshes:
                bpy.data.meshes.remove(me.get_data())
            # delete parent
            bpy.context.scene.objects.unlink(mol.objects.parent.get_object())
            bpy.data.objects.remove(mol.objects.parent.get_object())
            
            # Finally delete molecule from scene collection
            for i, mol_item in enumerate(self.molecules):
                if mol == mol_item:
                    mol_id = mol.name
                    name = mol.name_mol
                    self.molecules.remove(i)
                    debug_print("Deleted {} ({}).".format(mol_id, name), 
                                level=3)
                    return
        else:
            ob_list = ([ob.get_object() for ob in mol.objects.atoms] + 
                       [ob.get_object() for ob in mol.objects.other] + 
                       [ob.get_object() for ob in mol.objects.bonds])
            self.remove_objects(ob_list) # will call remove_molecule again
            
    def remove_object(self, ob):
        if ob.mb.type == 'MESH':
            debug_print("WARNING: Deleting individual mesh not necessary.",
                        level=1)
        elif ob.mb.type == 'ATOM':
            # remove link from bonds
            for b in ob.mb.bonds:
                b_ob = b.get_object()
                if b_ob:
                    b_ob.mb.remove_bonded_atom(ob)
        elif ob.mb.type == 'BOND':
            # remove link from atoms
            for a in ob.mb.bonded_atoms:
                a_ob = a.get_object()
                if a_ob:
                    a_ob.mb.remove_bond(ob)
        # remove link from molecule
        mol = ob.mb.get_molecule()
        if mol:
            mol.remove_object(ob)
        # delete object from blender
        bpy.context.scene.objects.unlink(ob)
        bpy.data.objects.remove(ob)
        # check if object was last object in molecule and delete molecule if so
        if (mol and 0 == len(mol.objects.atoms) == len(mol.objects.other) and
                0 == len(mol.objects.bonds)):
            self.remove_molecule(mol)
    
    def remove_objects(self, ob_list):
        # now collect all objects per molecule
        molecules = {}
        found_mesh = False
        for ob in ob_list:
            if ob.mb.type == 'MESH':
                found_mesh = True
            if ob.mb.type == 'ATOM':
                # remove link from bonds
                for b in ob.mb.bonds:
                    b_ob = b.get_object()
                    if b_ob:
                        b_ob.mb.remove_bonded_atom(ob)
            elif ob.mb.type == 'BOND':
                # remove link from atoms
                for a in ob.mb.bonded_atoms:
                    a_ob = a.get_object()
                    if a_ob:
                        a_ob.mb.remove_bond(ob)
            try:
                molecules[ob.mb.get_molecule()].append(ob)
            except KeyError:
                molecules[ob.mb.get_molecule()] = [ob]
        if found_mesh:
            debug_print("WARNING: Deleting individual meshes not necessary.", 
                        level=1)
        
        for mol, ob_list in molecules.items():
            if mol:
                mol.remove_objects(ob_list)
            for ob in ob_list:
                bpy.context.scene.objects.unlink(ob)
                bpy.data.objects.remove(ob)
        if (mol and 0 == len(mol.objects.atoms) == len(mol.objects.other) and
                0 == len(mol.objects.bonds)):
            self.remove_molecule(mol)


class mb_wm_globals(PropertyGroup):
    group_selected_extend = BoolProperty(
        name="Extend", description="Extend selected group to selection",
        default=False)
    element_to_add = StringProperty(
        name="Element", description="Element to add to scene", 
        default="C")
    geometry_to_add = EnumProperty(
        name="Geometry",
        description="Geometry the new bond should be in relative to "
                    "existing bonds. Press ALT to activate.", 
        items=mb_utils.enums.geometries, default='SINGLE')    


class mb_window_manager(PropertyGroup):
    globals = PointerProperty(type=mb_wm_globals)


def register():
    bpy.types.Object.mb = PointerProperty(type=mb_object)
    bpy.types.Mesh.mb = PointerProperty(type=mb_mesh)
    bpy.types.Action.mb = PointerProperty(type=mb_action)
    bpy.types.Scene.mb = PointerProperty(type=mb_scene)
    bpy.types.WindowManager.mb = PointerProperty(type=mb_window_manager)


def unregister():
    del bpy.types.Object.mb
    del bpy.types.Mesh.mb
    del bpy.types.Action.mb
    del bpy.types.Scene.mb
    del bpy.types.WindowManager.mb
