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

import os
import string
import random

import bpy
from bpy.types import (PropertyGroup,
                       UIList)
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolVectorProperty,
                       PointerProperty,
                       CollectionProperty,
                       EnumProperty)

from molblend.mb_helper import debug_print

class mb_object_pointer(PropertyGroup):
    #name = StringProperty(name="Object name")
    object = PointerProperty(name="Object", type=bpy.types.Object)
    @property
    def name(self):
        return self.object.name

class mb_element_mesh(PropertyGroup):
    name = StringProperty(name="Element")
    data = PointerProperty(name="Mesh", type=bpy.types.Mesh)
    #def get_data(self):
        #return bpy.data.meshes.get(self.name)

class atom_scale(PropertyGroup):
    name = StringProperty()
    val = FloatProperty(name="Atom scale", default=0.4, min=0.0, max=5.0,
                        precision=2, update=mb_utils.update_all_meshes)

class mb_mesh(PropertyGroup):
    type = EnumProperty(
        name="type", description="Identifies the mesh type",
        items=mb_utils.enums.mesh_types, default='NONE')

class mb_atom_mode(PropertyGroup):
    name = IntProperty(name="index")
    index = IntProperty(name="index")
    freq = FloatProperty(name="frequency")
    vec = FloatVectorProperty(name="vector", subtype="XYZ")

class mb_dipole(PropertyGroup):
    empty = PointerProperty(name="Dipole", type=mb_object_pointer)
    arrow = PointerProperty(name="Vector", type=mb_object_pointer)

class mb_unit_cell(PropertyGroup):
    uc_cube = PointerProperty(name="Unit cell base", type=mb_object_pointer)
    objects = CollectionProperty(name="Unit cell objects", type=mb_object_pointer)

class mb_molecule_objects(PropertyGroup):
    atoms = CollectionProperty(name="Atoms", type=mb_object_pointer)
    bonds = CollectionProperty(name="Bonds", type=mb_object_pointer)
    parent = PointerProperty(name="Parent", type=mb_object_pointer)
    dipole = PointerProperty(name="Dipole", type=mb_dipole)
    unit_cell = PointerProperty(name="Unit cell objects", type=mb_unit_cell)
    other = CollectionProperty(name="Bonds", type=mb_object_pointer)

class mb_mode_displacement(PropertyGroup):
    real = FloatVectorProperty(name="real", size=3)
    imag = FloatVectorProperty(name="imag", size=3)
    action = PointerProperty(name="action", type=bpy.types.Action)

class mb_mode(PropertyGroup):
    name = StringProperty(name="name")
    index = IntProperty(name="index")
    freq = StringProperty(
        name="Mode frequency",
        default="")
    symmetry = StringProperty(
        name="Mode symmetry", description="Symmetry of the active mode",
        default="")
    evecs = CollectionProperty(name="Displacements", type=mb_mode_displacement)

class mb_qmode(PropertyGroup):
    name = StringProperty(name="name")
    nqpt = IntProperty(name="nqpt", default=1)
    qvec = FloatVectorProperty(name="q vector", default=(0,0,0))
    modes = CollectionProperty(type=mb_mode)

class MB_UL_modes(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, 
                  active_propname, index):
        split = layout.split(0.1)
        col = split.column()
        col.label(str(item.nqpt))
        col = split.column()
        col.label("q=({:5.3f}, {:5.3f}, {:5.3f})".format(*item.qvec))

class mb_object(PropertyGroup):
    index = IntProperty(name="Index")
    name = StringProperty(name="Object name")
    type = EnumProperty(
        name="type", description="Select the object type",
        items=mb_utils.enums.object_types, default='NONE')
    
    molecule_ident = StringProperty(name="Molecule identifier")
    
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
    #modes = CollectionProperty(type=mb_atom_mode)
    
    supercell = IntVectorProperty(name="supercell", default=(0,0,0),
                                  size=3)
    
    @property
    def object(self):
        return self.id_data
    @property
    def world_location(self):
        return self.id_data.matrix_world.to_translation()
    
    def get_molecule(self):
        return bpy.context.scene.mb.molecules.get(self.molecule_ident)
    
    def add_bond(self, ob):
        """Add object to bond collection and return new collection item."""
        if not self.type == 'ATOM':
            debug_print("WARNING: Something is trying to add bond to "
                        "non-ATOM type object", level=1)
            return None
        
        bond = None
        for existing in self.bonds:
            if existing.object == ob:
                bond = existing
                break
        else:
            bond = self.bonds.add()
            bond.object = ob
        return bond
    
    def remove_bond(self, ob):
        if not self.type == 'ATOM':
            debug_print("WARNING: Something is trying to remove bond from "
                        "non-ATOM type object", level=1)
            return None
        for i, b in enumerate(self.bonds):
            if b.object == ob:
                self.bonds.remove(i)
                return
    
    def add_bonded_atom(self, ob):
        if not self.type == 'BOND':
            debug_print("WARNING: Something is trying to add bonded_atom to "
                        "non-BOND type object", level=1)
            return
        
        atom = None
        for existing in self.bonded_atoms:
            if existing.object == ob:
                atom = existing
                break
        else:
            atom = self.bonded_atoms.add()
            atom.object = ob
        return atom
    
    def remove_bonded_atom(self, ob):
        if not self.type == 'BOND':
            debug_print("WARNING: Something is trying to remove bonded_atom "
                        "{} ({}) from non-BOND type object {} ({})".format(
                        ob.name, ob.mb.type, self.name, self.type), level=1)
            return
        
        for i, a in enumerate(self.bonded_atoms):
            if a.object == ob:
                self.bonded_atoms.remove(i)
                return
    
    def draw_properties(self, context, layout, ob):
        element = context.scene.mb.elements[self.element]
        mat = ob.material_slots[0].material
        if context.scene.render.engine == 'CYCLES':
            atom_color = [mat.node_tree.nodes['Diffuse BSDF'].inputs[0],
                          "default_value", 40]
        else:
            atom_color = [mat, "diffuse_color", 40]
            
        props = {
            "Element": [self, "element", 10],
            "Name": [self, "atom_name", 20],
            "Atom radius": [element, "covalent", 30],
            "Atom color": atom_color,
            }
        for label, (data, prop, i) in sorted(props.items(), 
                                             key=lambda t: t[-1][-1]):
            row = layout.row()
            row.prop(data, prop, text=label)


#class mb_action(PropertyGroup):
    #name = StringProperty(name="Action name")
    #mode_vector = FloatVectorProperty(
        #name="Mode vector",
        #description="Original mode vector for atom as read from file")


class mb_molecule(PropertyGroup):
    index = IntProperty(name="Molecule index")
    name = StringProperty(name="Molecule identifier")
    name_mol = StringProperty(name="Molecule Name", 
                              update=mb_utils.update_molecule_name)
    # index that increases with each added atom in the molecule, but doesn't 
    # decrease when atom is deleted. => Not an indicator of size of molecule! 
    # Only guarantees uniqueness for atom names
    atom_index = IntProperty(name="Atom counter")
    objects = PointerProperty(name="Molecule objects", 
                              type=mb_molecule_objects)
    meshes = CollectionProperty(name="Meshes", type=mb_element_mesh)
    
    # display properties
    bond_material = EnumProperty(
        name="Bond material", description="Choose bond material",
        items=mb_utils.enums.bond_material, default='ATOMS',
        update=mb_utils.update_bond_material
        )
    bond_color = FloatVectorProperty(
        name='Bond color', description="Color for generic bond",
        default=(0.8, 0.8, 0.8), min=0.0, max=1.0, subtype='COLOR',
        update=mb_utils.update_all_meshes
        )
    draw_style = EnumProperty(
        name="Display style", description="Style to draw atoms and bonds",
        items=mb_utils.enums.molecule_styles, default='BAS',
        update=mb_utils.set_draw_style
        )
    radius_type = EnumProperty(
        name="Radius type", description="Type of radius to use as reference",
        items=mb_utils.enums.radius_types, default='covalent', 
        update=mb_utils.update_radius_type
        )
    bond_radius = FloatProperty(
        name="Bond radius", description="Radius for bond objects",
        default=0.1, min=0.0, max=3.0,
        update=mb_utils.update_all_meshes
        )
    atom_scales = CollectionProperty(type=atom_scale)
    refine_atoms = IntProperty(
        name="Refine atoms", description="Refine value for atom meshes", 
        default=8, min=2, max=64,
        update=mb_utils.update_refine_atoms
        )
    refine_bonds = IntProperty(
        name="Refine bonds", description="Refine value for atom meshes", 
        default=8, min=2, max=64,
        update=mb_utils.update_refine_bonds
        )
    
    # vibrational modes
    qpts = CollectionProperty(name="Modes", type=mb_qmode)
    active_nqpt = IntProperty(
        name="Active q-point",
        description="Active q-point",
        default=0, min=0,
        update=mb_utils.update_active_mode
        )
    max_mode = IntProperty(
        name="Number of modes", 
        description="Number of vibrational modes of molecule",
        default=0, min=0,
        )
    active_mode = IntProperty(
        name="Active Mode",
        description="Active Mode to display. 0 = equilibrium position",
        default=0, min=0,
        update=mb_utils.update_active_mode
        )
    mode_scale = FloatProperty(
        name="Mode Scale", description="Scale of normal mode displacement",
        default=1.0, min=-1000.0, max=1000.0, 
        update=mb_utils.update_active_mode
        )
    
    def draw_vibrations(self, layout):
        
        row = layout.row()
        row.operator("mb.import_modes")
        if self.qpts:
            row = layout.row()
            row.template_list("MB_UL_modes", "", self, "qpts", self, 
                              "active_nqpt",rows=1)
            row = layout.row()
            row.prop(self, "active_mode")
            row = layout.row()
            row.prop(self, "mode_scale")
            row = layout.row()
            row.prop(self.qpts[self.active_nqpt].modes[self.active_mode], 
                     "freq", text="Frequency")
            row = layout.row()
            row.prop(self.qpts[self.active_nqpt].modes[self.active_mode], 
                     "symmetry", text="Symmetry")
        
    
    def draw_properties(self, layout):
        #layout.label("Molecule properties")
        row = layout.row()
        row.prop(self, "name_mol", text="")
        row = layout.row()
        row.label("(parent: {}".format(self.objects.parent.name))
        row = layout.row()
        row.label("(parent_ob: {}".format(self.objects.parent.object.name))
        row = layout.row()
        row.label("(ident: '{}')".format(self.name))
        
        #row = layout.row()
        #row.active = bool(self.vibrations.max_mode)
        #row.prop(self.vibrations, "active_mode")
        #row = layout.row()
        #row.active = bool(self.vibrations.max_mode)
        #row.prop(self.vibrations, "mode_scale")
        
        if self.objects.parent.object:
            col = layout.column()
            col.prop(self.objects.parent.object, "location", 
                     text="Center of mass")
        row = layout.row()
        row.operator("mb.center_mol_parent")
        
        if self.objects.dipole.empty.object:
            col = layout.column()
            col.prop(self.objects.dipole.empty.object, "location", 
                     text="Dipole")
        else:
            row = layout.row()
            row.operator("mb.draw_dipole")
        #if not self.objects.unit_cell.uc_cube.object:
        row = layout.row()
        row.operator("mb.draw_unit_cell")
        #else:
            #uc_cube =self.objects.unit_cell.uc_cube.object
            #vgs = [vg for vg in uc_cube.vertex_groups 
                   #if vg.name in ("a", "b", "c")]
            #vg_ids = [vg.index for vg in sorted(vgs, key=lambda vg: vg.name)]
            #verts = []
            #for vg_id in vg_ids:
                #vs = [v for v in uc_cube.data.vertices 
                      #if vg_id in [vg.group for vg in v.groups]]
                #if len(vs) != 1:
                    #msg = "ERROR: Something wrong with uc_cube vertext groups"
                    #debug_print(msg, level=1)
                    #return
                #verts.append(vs[0])
            #for v in verts:
                #row = layout.row()
                #row.prop(v, "co")
            
    def draw_styles(self, layout):
        #layout.label("Molecule draw style")
        row = layout.row()
        row.label(self.name_mol)
        row = layout.row()
        row.label("(parent: {}".format(self.objects.parent.name))
        row = layout.row()
        row.label("(ident: {}".format(self.name))
        
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
        
    def add_object(self, ob):
        '''
        Add an object to the molecule's atoms collection and return the
        collection item. If object is already in collection, just return the
        collection item.
        '''
        collection = {
            'ATOM': self.objects.atoms,
            'BOND': self.objects.bonds,
            'NONE': self.objects.other
            }
        objects = collection[ob.mb.type]
        
        for existing in objects:
            if existing.object == ob:
                item = existing
                break
        else:
            item = objects.add()
            item.object = ob
        return item
    
    def remove_object(self, ob):
        collection = {
            'ATOM': self.objects.atoms,
            'BOND': self.objects.bonds,
            'NONE': self.objects.other
            }
        objects = collection[ob.mb.type]
        
        for i, a in enumerate(objects):
            if a == ob:
                objects.remove(i)
                return
    
    def remove_objects(self, ob_list=None):
        ob_list = ob_list or []
        objects = {'ATOM': [ob for ob in ob_list if ob.mb.type == 'ATOM'],
                   'BOND': [ob for ob in ob_list if ob.mb.type == 'BOND'],
                   'NONE': [ob for ob in ob_list if ob.mb.type == 'NONE']}
        
        collections = {'ATOM': self.objects.atoms,
                      'BOND': self.objects.bonds,
                      'NONE': self.objects.other}
        for ob_type, obs in objects.items():
            # find all the indices of the objects to delete
            indeces = []
            for i, a in enumerate(collections[ob_type]):
                if a in obs:
                    indeces.append(i)
            # delete higher numbers first to not mess up the order of the coll.
            for i in reversed(indeces):
                collections[ob_type].remove(i)
        return


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
    n_q = IntProperty(name="q point", default=1,
        min=1, description="Import modes of 'n_q'th q point in file")
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
    show_bond_lengths = BoolProperty(
        name="Show bond lengths", default=False,
        description="Display bond length of selected bonds",
        update=mb_utils.update_show_bond_lengths)
    #show_bond_angles = BoolProperty(
        #name="Show bond angles", default=False, 
        #description="Display bond angle of selected bonds",
        #update=mb_utils.update_show_bond_angles)
    element_to_add = StringProperty(
        name="Element", description="Element to add to scene", 
        default="C")
    geometry_to_add = EnumProperty(
        name="Geometry",
        description="Geometry the new bond should be in relative to "
                    "existing bonds.", 
        items=mb_utils.enums.geometries, default='NONE')


class mb_scene(PropertyGroup):
    is_initialized = BoolProperty(default=False)
    elements = CollectionProperty(type=mb_element)
    molecules = CollectionProperty(type=mb_molecule)
    # index that increases with each added molecule, but doesn't decrease when
    # molecule is deleted.
    molecule_count = IntProperty(name="Molecule counter")
    globals = PointerProperty(type=mb_scn_globals)
    # store last active object for modal operator
    modal_last_active = PointerProperty(name="Last active", 
                                        type=bpy.types.Object)
    
    #@class
    def id_generator(self, size=6,
                     chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))
    
    def new_molecule(self,
                     name_mol=None,
                     draw_style=None,
                     radius_type=None,
                     bond_radius=None,
                     refine_atoms=None,
                     refine_bonds=None,
                     atom_scales=None):
        mol = self.molecules.add()
        
        new_id = self.id_generator()
        while new_id in self.molecules:
            new_id = self.id_generator()
        mol.name = new_id
        
        mol.index = self.molecule_count
        self.molecule_count += 1
        mol.draw_style = draw_style or self.globals.draw_style
        mol.radius_type = radius_type or self.globals.radius_type
        mol.bond_radius = bond_radius or self.globals.bond_radius
        if refine_atoms:
            mol.refine_atoms = refine_atoms
        if refine_bonds:
            mol.refine_bonds = refine_bonds
        for scale in (atom_scales or self.globals.atom_scales):
            new_scale = mol.atom_scales.add()
            new_scale.name = scale.name
            new_scale.val = scale.val
        # create new empty that will be the parent for the molecule
        parent_ob = bpy.data.objects.new(mol.name_mol, None)
        parent_ob.empty_draw_type = 'SPHERE'
        parent_ob.empty_draw_size = 0.3
        bpy.context.scene.objects.link(parent_ob)
        mol.objects.parent.object = parent_ob
        parent_ob.mb.molecule_ident = mol.name
        parent_ob.mb.type = 'PARENT'
        mol.name_mol = name_mol or "Molecule"
        return mol
    
    def remove_molecule(self, mol):
        # Make sure all objects are deleted first
        if (0 == len(mol.objects.atoms) == len(mol.objects.bonds)
              == len(mol.objects.unit_cell) == len(mol.objects.other)):
            # delete all meshes
            for me in mol.meshes:
                bpy.data.meshes.remove(me.data)
            # delete parent
            parent = mol.objects.parent.object
            if parent:
                if parent.name in bpy.context.scene.objects:
                    bpy.context.scene.objects.unlink(parent)
                else:
                    debug_print("Object {} not in scene {}.".format(
                        parent.name, bpy.context.scene.name), level=1)
                try:
                    bpy.data.objects.remove(parent)
                except RuntimeError:
                    debug_print(
                        "Object {} ".format(parent.name)
                        + "can not be deleted from bpy.data.objects.",
                        level=1)
            # delete dipole
            for ob in (mol.objects.dipole.empty.object,
                       mol.objects.dipole.arrow.object):
                if ob:
                    if dipole.name in bpy.context.scene.objects:
                        bpy.context.scene.objects.unlink(dipole)
                    else:
                        debug_print("Object {} not in scene {}.".format(
                            dipole.name, bpy.context.scene.name), level=1)
                    try:
                        bpy.data.objects.remove(dipole)
                    except RuntimeError:
                        debug_print(
                            "Object {} ".format(dipole.name)
                            + "can not be deleted from bpy.data.objects.",
                            level=1)
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
            ob_list = ([item.object for item in mol.objects.atoms] +
                       [item.object for item in mol.objects.bonds] +
                       [item.object for item in mol.objects.unit_cell] +
                       [item.object for item in mol.objects.other]
                       )
            self.remove_objects(ob_list) # will call remove_molecule again
            
    def remove_object(self, ob):
        if ob.mb.type == 'ATOM':
            # remove link from bonds
            for b in ob.mb.bonds:
                b_ob = b.object
                if b_ob:
                    b_ob.mb.remove_bonded_atom(ob)
        elif ob.mb.type == 'BOND':
            # remove link from atoms
            for a in ob.mb.bonded_atoms:
                a_ob = a.object
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
        if (mol and 0 == len(mol.objects.atoms) == len(mol.objects.other)
                      == len(mol.objects.bonds)):
            self.remove_molecule(mol)
    
    def remove_objects(self, ob_list):
        # now collect all objects per molecule
        molecules = {}
        for ob in ob_list:
            if ob.mb.type == 'ATOM':
                # remove link from bonds
                for b in ob.mb.bonds:
                    b_ob = b.object
                    if b_ob:
                        b_ob.mb.remove_bonded_atom(ob)
            elif ob.mb.type == 'BOND':
                # remove link from atoms
                for a in ob.mb.bonded_atoms:
                    a_ob = a.object
                    if a_ob:
                        a_ob.mb.remove_bond(ob)
            try:
                molecules[ob.mb.get_molecule()].append(ob)
            except KeyError:
                molecules[ob.mb.get_molecule()] = [ob]
        
        for mol, ob_list in molecules.items():
            if mol:
                mol.remove_objects(ob_list)
            for ob in ob_list:
                bpy.context.scene.objects.unlink(ob)
                bpy.data.objects.remove(ob)
        if (mol and 0 == len(mol.objects.atoms) == len(mol.objects.other)
                      == len(mol.objects.bonds)):
            self.remove_molecule(mol)


def register():
    bpy.types.Object.mb = PointerProperty(type=mb_object)
    bpy.types.Mesh.mb = PointerProperty(type=mb_mesh)
    bpy.types.Scene.mb = PointerProperty(type=mb_scene)


def unregister():
    del bpy.types.Object.mb
    del bpy.types.Mesh.mb
    del bpy.types.Scene.mb
