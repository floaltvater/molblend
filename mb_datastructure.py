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
#------IMPORTS

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
    
    def get_object(self):
        return bpy.data.meshes.get(self.name)

class atom_scale(PropertyGroup):
    name = StringProperty()
    val = FloatProperty(name="Atom scale", default=0.3, min=0.0, max=5.0,
                        precision=2, update=mb_utils.update_all_meshes)

class mb_object(PropertyGroup):
    
    index = IntProperty(name="Index")
    name = StringProperty(name="Object name")
    type = EnumProperty(name="type", items=mb_utils.enums.object_types,
                        description="Select the object type",
                        default='NONE')
    id_molecule = StringProperty(name="Molecule ID")
    bonds = CollectionProperty(type=mb_object_pointer)
    bonded_atoms = CollectionProperty(type=mb_object_pointer)
    atom_name = StringProperty(name="Atom name")
    element = StringProperty(name="Element",
                             description="Element Symbol",
                             update=mb_utils.update_atom_element)
    element_long = StringProperty(name="Element name", 
                                  description="Full element name")
    
    def get_molecule(self):
        return bpy.context.scene.mb.molecules.get(self.id_molecule)
    
    def get_object(self):
        return bpy.data.objects.get(self.name)
        
    def draw_properties(self, layout, name):
        row = layout.row()
        row.prop(self, "element")
        row = layout.row()
        row.prop(self, "atom_name")
        row = layout.row()
        row.prop(bpy.context.scene.objects[name], "name")

class mb_molecule(PropertyGroup):
    index = IntProperty(name="Molecule index")
    name = StringProperty(name="Molecule identifier")
    name_mol = StringProperty(name="Molecule Name")
    atoms = CollectionProperty(name="Atoms", type=mb_object_pointer)
    # index that increases with each added atom in the molecule, but doesn't decrease when
    # atom is deleted. => Not an indicator of size of molecule! Only guarantees uniqueness for atom names
    atom_index = IntProperty()
    bonds = CollectionProperty(name="Bonds", type=mb_object_pointer)
    bond_material = EnumProperty(name="Bond material",
                                description="Choose bond material",
                                items=mb_utils.enums.bond_material, default='ATOMS',
                                update=mb_utils.update_bond_material)
    bond_color = FloatVectorProperty(name='Bond color', default=(0.8, 0.8, 0.8),
                                     min=0.0, max=1.0, subtype='COLOR',
                                     update=mb_utils.update_all_meshes)
    # dito
    meshes = CollectionProperty(name="Molecule meshes", type=mb_mesh_pointer)
    draw_style = EnumProperty(name="Display style", items=mb_utils.enums.molecule_styles,
                         description="Style to draw atoms and bonds",
                         default='BAS',
                         update=mb_utils.set_draw_style)
    radius_type = EnumProperty(name="Radius type", items=mb_utils.enums.radius_types,
                               default='covalent', update=mb_utils.update_radius_type)
    bond_radius = FloatProperty(name="Bond radius", default=0.15, min=0.0, max=3.0,
                                description="Radius of bonds for Sticks, and Ball and Sticks",
                                update=mb_utils.update_all_meshes)
    atom_scales = CollectionProperty(type=atom_scale)
    parent = PointerProperty(type=mb_object_pointer)
    
    def add_atom(self, ob, replace_name=""):
        '''
        add an atom's object name to the molecule's atoms collection
        if replace_name is given, the first atom's name is replaced by the new one
        if name is already in collection, print a Warning and just return the atom
        '''
        if ob.name in self.atoms:
            debug_print("WARNING: {} already part of molecule.".format(ob.name), 1)
            return self.atoms[ob.name]
        
        if replace_name:
            for atom in self.atoms:
                if atom.name == replace_name:
                    atom.name = ob.name
                    return atom
        else:
            atom = self.atoms.add()
            atom.name = ob.name
            return atom
        
    def add_bond(self, ob):
        bond = self.bonds.add()
        bond.name = ob.name
        
    def remove_atom(self, ob):
        for index, a in enumerate(self.atoms):
            if a.name == ob.name:
                self.atoms.remove(index)
    
class mb_element(PropertyGroup):
    name = StringProperty(name="Element")
    element = StringProperty(name="Element")
    element_name = StringProperty(name="Element name")
    atomic_number = IntProperty(name="Atomic number")
    color = FloatVectorProperty(name="Color", subtype='COLOR', size=3)
    covalent = FloatProperty(name="Covalent radius", min=0.0, max=5.0, 
                             update=mb_utils.update_all_meshes)
    vdw = FloatProperty(name="vdW radius", min=0.0, max=5.0, 
                        update=mb_utils.update_all_meshes)
    constant = FloatProperty(name="Constant radius", min=0.0, max=5.0,
                             update=mb_utils.update_all_meshes)

class mb_scn_globals(PropertyGroup):
    draw_style = EnumProperty(name="Draw style", items=mb_utils.enums.molecule_styles,
                         description="Style to draw atoms and bonds", default='BAS')
    radius_type = EnumProperty(name="Radius type", items=mb_utils.enums.radius_types,
                               default='covalent')
    bond_radius = FloatProperty(name="Bond radius", default=0.15, min=0.0, max=3.0,
                                description="Radius of bonds for Sticks, and Ball and Sticks")
    atom_scales = CollectionProperty(type=atom_scale)

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
        parent_ob = bpy.data.objects.new("molecule_{}".format(mol.index), None)
        parent_ob.empty_draw_type = 'SPHERE'
        parent_ob.empty_draw_size = 0.3
        bpy.context.scene.objects.link(parent_ob)
        mol.parent.name = parent_ob.name
        return mol
    #view = bpy.props.CollectionProperty(type=mb_view)
    #display = bpy.props.PointerProperty(type=mb_display_properties)

class mb_wm_globals(PropertyGroup):
    group_selected_extend = BoolProperty(name="Extend", default=False, description="Extend selected group to selection")
    element_to_add = StringProperty(name="Element", default="C", description="Element to add to scene")
    geometry_to_add = EnumProperty(name="Geometry", default='VIEW', items=mb_utils.enums.geometries,
                      description="Geometry the new bond should be in relative to existing bonds. Press ALT to activate.")
    active_tool = EnumProperty(name="Tool", items=mb_utils.enums.mb_tools,
                  description="Select active tool")

class mb_window_manager(PropertyGroup):
    is_running = BoolProperty(default=False, description="MolBlend is running.")
    globals = PointerProperty(type=mb_wm_globals)
    last_keyconfig = StringProperty(name="Last keyconfig")
    
def register():
    bpy.types.Object.mb = PointerProperty(type=mb_object)
    #bpy.types.Group.mb = PointerProperty(type=mb_group)
    #bpy.types.World.mb = PointerProperty(type=mb_world)
    bpy.types.Scene.mb = PointerProperty(type=mb_scene)
    bpy.types.WindowManager.mb = PointerProperty(type=mb_window_manager)

def unregister():
    #del bpy.types.Group.mb
    del bpy.types.Object.mb
    del bpy.types.Scene.mb
    #del bpy.types.World.mb
    #del bpy.types.WindowManager.mb


#class extend_blender_data():
    #bpy.types.Material.mb = PointerProperty(type = mb_material)
    #bpy.types.Object.mb = PointerProperty(type = mb_object)
    #bpy.types.Group.mb = PointerProperty(type = mb_group)
    #bpy.types.Scene.mb = PointerProperty(type = mb_scene)
    #bpy.types.World.mb = PointerProperty(type = mb_world)
    #bpy.types.WindowManager.mb = PointerProperty(type = mb_window_manager)