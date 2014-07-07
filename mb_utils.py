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

import bpy
import bmesh
import math
from bpy.types import Operator

from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolVectorProperty,
                       PointerProperty,
                       CollectionProperty,
                       EnumProperty)

from bpy_extras import view3d_utils

from mathutils import Vector, Matrix

import math

import re

from bisect import bisect_left

from .elements_default import ELEMENTS as ELEMENTS_DEFAULT

from .helper import debug_print

import time

class enums():
    object_types = [('NONE', "None", "None"),
                    ('ATOM', "Atom", "Atom"),
                    #('MOLECULE', "Molecule", "Molecule"),
                    ('BOND', "Bond", "Bond"),
                    #('CHARGE', "Charge", "Charge"),
                   ]
    #group_types = [('NONE', "None", "None"),
                    #('MOLECULE', "Molecule", "Molecule"),
                   #]
    mb_tools = [('BLENDER', "Blender", "Use normal Blender input"),
                ('ADD_ATOM', "Add atom", "Add atom"),
                ('SELECT', "Select", "Select"),
                #('MOVE', "Move", "Move"),
                #('ROTATE', "Rotate", "Rotate")
               ]
    radius_types = [('covalent', 'covalent', 'covalent'),
                    ('vdw', 'van der Waals', 'van der Waals'),
                    ('constant', 'constant', 'constant'),
                   ]
    molecule_styles = [('BALLS', 'Balls', 'Space filling'),
                       ('BAS', 'Balls and Sticks', 'Balls and Sticks'),
                       ('STICKS', 'Sticks', 'Sticks'),
                      ]
    bond_material = [('ATOMS', "Atoms", "Same as atoms"),
                     ('GENERIC', "Generic" , "Single bond color"),
                    ]
    #geometries = [('VIEW', "Planar angles", "Angles are multiples of 30 and 45 deg. in the view plane"),
                  #('OCTAHEDRAL', "Octahedral", "Octahedral"),
                  #('TETRAHEDRAL', "Tetrahedral", "Tetrahedral or sp3"),
                  #('TRIGONAL', "Trig. planar", "Trigonal planar or sp2"),
                  #('LINEAR', "Linear", "Linear or sp")]
    geometries = [('SINGLE', "Single atom", "Angles are multiples of 30 and 45 deg. in the view plane"),
                  ('LINEAR', "Linear", "Linear or sp"),
                  ('TRIGONAL', "Trig. planar", "Trigonal planar or sp2"),
                  ('TETRAHEDRAL', "Tetrahedral", "Tetrahedral or sp3"),
                  ('OCTAHEDRAL', "Octahedral", "Octahedral"),
                 ]
    angstrom_per_unit = [('1.0', "Angstrom", "Angstrom"),
                         ('0.529177249', "Bohr", "Bohr"),
                         ('0.01', "pm", "Picometer"),
                         ('OTHER', "Other", "Custom Unit"),
                        ]
    file_types = [('XYZ', "*.xyz", "xyz format"),
                  ('PDB', "*.pdb", "Protein Databank format"),
                 ]

def create_new_keyconfig(context):
    
    if not "MolBlend" in context.window_manager.keyconfigs:
        debug_print('Creating new keyconfiguration "MolBlend".', 3)
        from molblend import mb_keyconfig
        #keyconfig = mb_keyconfig.kc
    else:
        debug_print('Keyconfiguration "MolBlend" already exists.', 3)
        
    #print(keyconfig)
    return context.window_manager.keyconfigs.get("MolBlend")

def add_element(context, element, element_dict):
    '''
    add element data to scene
    '''
    default = ELEMENTS_DEFAULT["Default"]
    
    new = context.scene.mb.elements.add()
    new.name = element
    new.element = element
    new.element_name = element_dict.get("element name", default["element name"])
    new.atomic_number = element_dict.get("atomic number", default["atomic number"])
    new.color = element_dict.get("color", default["color"])
    new.covalent = element_dict.get("covalent", default["covalent"])
    if "vdw1" in element_dict or "vdw2" in element_dict:
        new.vdw = (element_dict["vdw1"] if element_dict["vdw1"] != 0.0 
                    else (element_dict["vdw2"] if element_dict["vdw2"] != 0.0 
                    else element_dict["covalent"]))
    else:
        new.vdw = element_dict.get("covalent", default["covalent"])
    new.constant = 1.0
    return new
    
def initialize_elements(context):
    for element, data in ELEMENTS_DEFAULT.items():
        add_element(context, element, data)
        
def update_all_meshes(self, context):
    # TODO this callback might be to heavy for scenes with lots of meshes
    for me in bpy.data.meshes:
        me.update()

def update_active_mode(self, context):
    for atom in self.atoms:
        atom_ob = atom.get_object()
        atom_id = '{}.{}'.format(self.index, atom_ob.mb.index)
        if self.active_mode == 0:
            action_name = "equilibrium_{}".format(atom_id)
        else:
            action_name = "mode_{}_{}".format(self.active_mode, atom_id)
        action = bpy.data.actions.get(action_name)
        if not action:
            debug_print("Mode {} not found for molecule '{}'.".format(self.active_mode, self.name_mol), 1)
            return
        anim_data = atom_ob.animation_data
        if not anim_data:
            debug_print("No animation data for molecule '{}'.".format(self.name_mol), 1)
            return
        anim_data.action = action
            
    
def mouse_2d_to_location_3d(context, coord, depth=Vector((0, 0, 0))):
    region = context.region
    rv3d = context.region_data
    if depth:
        depth_location = depth
    else:
        depth_location = context.scene.cursor_location.copy()
    
    return view3d_utils.region_2d_to_location_3d(region, rv3d, coord, depth_location)

def get_atom_data(element, molecule, type='MESH', refine=8):
    if type == 'MESH':
        mesh_name = "atom_mesh_{}.{}".format(molecule.index, element)
        me = bpy.context.blend_data.meshes.get(mesh_name)
        if not me:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            debug_print("Create {}.".format(mesh_name), 3)
            
            # create uv sphere and get mesh data
            bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0), segments=refine*2, ring_count=refine)
            new_atom = bpy.context.object
            bpy.ops.object.shade_smooth()
            me = new_atom.data
            me.name = mesh_name
            
            # adds material slot to mesh, but don't assign material yet
            new_atom.data.materials.append(None)
            
            # finally delete object and return mesh data
            bpy.context.scene.objects.unlink(new_atom)
            
            # restore old selection
            for o in selected:
                o.select = True
            bpy.context.scene.objects.active = last_active
        return me

def get_bond_data(type='MESH', refine=8):
    # TODO maybe have one bond mesh per molecule. Problem: Update meshes when joining molecules?
    new_bond = None
    if type == 'MESH':
        data_name = "bond_mesh"
        data = bpy.context.blend_data.meshes.get(data_name)
        if not data:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            debug_print("Create {}.".format(data_name), 3)
            
            bpy.ops.mesh.primitive_cylinder_add(location=(0,0,0), vertices=refine*2, depth=1, radius=1.0,
                        end_fill_type="NOTHING")
            new_bond = bpy.context.object
            for i in range(2):
                new_bond.data.materials.append(None)
            
            data = new_bond.data
            data.name = data_name
            bm = bmesh.new()
            bm.from_mesh(data)
            
            # rotate and shrink first, then add another row of vertices
            for vert in bm.verts:
                # rotate 90 degrees around x, and shift along y axis
                tmp_co = vert.co.copy()
                vert.co.y = -tmp_co.z + .5
                vert.co.z = -tmp_co.y
                if vert.co.y > 0.01:
                    vert.select = False
            new_verts = []
            for edge in bm.edges:
                if edge.calc_length() == 1.0:
                    e, v = bmesh.utils.edge_split(edge, edge.verts[0], 0.5)
                    new_verts.append(v)
            n_verts = len(new_verts)
            #for i in range(n_verts):
                #bm.edges.new((new_verts[i], new_verts[(i+1)%n_verts]))
            
            # little hacky, but don't understand how bmesh.utils.face_split works
            # remove all faces
            for f in bm.faces:
                bm.faces.remove(f)
            # now sort bm.verts
            # v.co.y is either 0, 0.5, or 1.0. So multiply y with at least 4pi to sort by y value first
            key = lambda v: v.co.y * 15 + math.atan2(v.co.x, v.co.z)
            verts_sorted = sorted((v for v in bm.verts), key=key)
            for i, v in enumerate(verts_sorted):
                v.index = i
            
            # add new faces and assign the two different material slots
            for i in range(2*n_verts):
                v1 = verts_sorted[i]
                v2 = verts_sorted[(i + 1)%n_verts + n_verts*(i//n_verts)]
                v3 = verts_sorted[(i + 1)%n_verts + n_verts*(i//n_verts + 1)]
                v4 = verts_sorted[i + n_verts]
                f = bm.faces.new((v1, v2, v3, v4))
                f.material_index = i//n_verts # gives 0 or 1
                f.smooth = True
            
            # again, sort by center.y first, than angle
            key = lambda f: f.calc_center_median().y * 15 + math.atan2(f.normal.x, f.normal.z)
            half_faces = len(bm.faces)/2
            for i, f in enumerate(sorted((f for f in bm.faces), key=key)):
                f.index = i
            
            bm.to_mesh(data)
            bm.free()
            data.update()
            #bpy.ops.object.shade_smooth()

            
    if type == 'CURVE':
        # save last selection to restore later
        selected = bpy.context.selected_objects
        last_active = bpy.context.object
        
        # get bevel object
        bond_bevel = bpy.context.blend_data.objects.get('bond_bevel')
        if not bond_bevel:
            debug_print("Create bond_bevel.", 3)
            bpy.ops.curve.primitive_bezier_circle_add(radius=0.15, location=(0,0,0))
            bond_bevel = bpy.context.object
            bond_bevel.name = 'bond_bevel'
            
        bpy.ops.curve.primitive_bezier_curve_add(location=(0,0,0))
        new_bond = bpy.context.object
        
        data = new_bond.data
        data.show_handles = False
        
        bp = data.splines[0].bezier_points
        for i in range(2):
            bp[i].handle_left_type = 'VECTOR'
            bp[i].handle_right_type = 'VECTOR'
        
        data.bevel_object = bond_bevel
        #bp[0].co = (0,0,0)
        #bp[1].co = (0,1,0)
    
    if new_bond:
        # finally delete object and reselect old selection
        bpy.context.scene.objects.unlink(new_bond)
        
        for o in selected:
            o.select = True
        bpy.context.scene.objects.active = last_active
    
    return data

def set_atom_drivers(context, atom, molecule):
    fc_list = atom.driver_add('scale', -1) # add new driver
    for fcurve in fc_list:
        drv = fcurve.driver
        drv.type = 'SCRIPTED'
        drv.show_debug_info = True
        
        var = drv.variables.get('atom_radius')
        if not var:
            var = drv.variables.new()
            var.name = 'atom_radius' # name to use in scripting
            var.type = 'SINGLE_PROP'
        targ = var.targets[0]
        targ.id_type = 'SCENE'
        targ.id = context.scene
        targ.data_path = 'mb.elements["{}"].{}'.format(atom.mb.element, molecule.radius_type)
        
        var = drv.variables.get('atom_scale')
        if not var:
            var = drv.variables.new()
            var.name = 'atom_scale' # name to use in scripting
            var.type = 'SINGLE_PROP'
        targ = var.targets[0]
        targ.id_type = 'SCENE'
        targ.id = context.scene
        targ.data_path = 'mb.molecules["{}"].atom_scales["{}"].val'.format(molecule.name, molecule.draw_style)
        
        var = drv.variables.get('bond_radius')
        if not var:
            var = drv.variables.new()
            var.name = 'bond_radius' # name to use in scripting
            var.type = 'SINGLE_PROP'
        targ = var.targets[0]
        targ.id_type = 'SCENE'
        targ.id = context.scene
        targ.data_path = 'mb.molecules["{}"].bond_radius'.format(molecule.name)
        
        if molecule.draw_style in ('BALLS', 'BAS'):
            drv.expression = "max(atom_radius * atom_scale, bond_radius)"
        elif molecule.draw_style == 'STICKS':
            drv.expression = "bond_radius"
            

def set_bond_drivers(context, bond, molecule, type='MESH'):
    if type == 'MESH':
        fc_x = bond.driver_add('scale', 0)
        fc_z = bond.driver_add('scale', 2)
        for fcurve in (fc_x, fc_z):
            drv = fcurve.driver
            drv.type = 'AVERAGE'
            drv.show_debug_info = True
            
            var = drv.variables.get('bond_radius')
            if not var:
                var = drv.variables.new()
                var.name = 'bond_radius' # name to use in scripting
                var.type = 'SINGLE_PROP'
            targ = var.targets[0]
            targ.id_type = 'SCENE'
            targ.id = context.scene
            targ.data_path = 'mb.molecules["{}"].bond_radius'.format(molecule.name)

    elif type == 'CURVE':
        bp = bond.data.splines[0].bezier_points
        for i in range(2):
            fc_list = bp[i].driver_add('co', -1)
            for dim, fcurve in zip(('X', 'Y', 'Z'), fc_list):
                drv = fcurve.driver
                drv.type = 'AVERAGE'
                drv.show_debug_info = True
            
                var = drv.variables.get('atom_location')
                if not var:
                    var = drv.variables.new()
                    var.name = 'atom_location' # name to use in scripting
                    var.type = 'TRANSFORMS'
                targ = var.targets[0]
                #targ.id_type = 'OBJECT'
                targ.id = bond.mb.bonded_atoms[i].get_object()
                #targ.data_path = 'location'
                targ.transform_type = 'LOC_{}'.format(dim)
                targ.transform_space = 'TRANSFORM_SPACE'

def new_material(name, color=(0.8, 0.8, 0.8), molecule=None):
    '''
    creates new material. If molecule is given, the molecule.mb.bond_color property will be added as driver.
    '''
    material = bpy.data.materials.new(name)
    #scn_elements = bpy.context.scene.mb.elements
    #color = scn_elements.get(element, scn_elements["Default"]).color
    if molecule == None:
        material.diffuse_color = color
    else:
        for i in range(3):
            fcurve =  material.driver_add('diffuse_color', i) # add new driver
            drv = fcurve.driver
            drv.type = 'AVERAGE'
            drv.show_debug_info = True
            
            var = drv.variables.new()
            var.name = 'diffuse_color' # name to use in scripting
            var.type = 'SINGLE_PROP'
            targ = var.targets[0]
            targ.id_type = 'SCENE'
            targ.id = bpy.context.scene
            targ.data_path = 'mb.molecules["{}"].bond_color[{}]'.format(molecule.name, i)
    
    if bpy.context.scene.render.engine == 'CYCLES':
        material.use_nodes = True
        # add alpha value to atom.color
        #color_alpha = list(color) + [1.0]
        #material.node_tree.nodes['Diffuse BSDF'].inputs[0].default_value = color_alpha
        
        # add driver to rendered color to be the same as display color
        nodesocketcolor = material.node_tree.nodes['Diffuse BSDF'].inputs[0]
        for i in range(3): # not for alpha channel
            fcurve =  nodesocketcolor.driver_add('default_value', i) # add new driver
            drv = fcurve.driver
            drv.type = 'AVERAGE'
            drv.show_debug_info = True
            
            var = drv.variables.new()
            var.name = 'diffuse_color' # name to use in scripting
            var.type = 'SINGLE_PROP'
            targ = var.targets[0]
            targ.id_type = 'MATERIAL'
            targ.id = material
            targ.data_path = 'diffuse_color[{}]'.format(i)
    return material

def assign_atom_material(ob, molecule):
    # make sure there is at least one material slot
    if len(ob.material_slots) < 1:
        ob.data.materials.append(None)
    
    # get element and molecule index of atom to make unique material name
    element = ob.mb.element
    material_name = "mat_{}_{}".format(element, molecule.index)
    
    # get or create element material per molecule
    material = bpy.data.materials.get(material_name)
    if not material:
        scn_elements = bpy.context.scene.mb.elements
        color = scn_elements.get(element, scn_elements["Default"]).color
        # get material color from elements list, and Default if not an element
        material = new_material(material_name, color=color)
    
    # finally, assign material to first slot.
    ob.material_slots[0].link = 'DATA'
    ob.material_slots[0].material = material

def assign_bond_material(ob):
    
    #bond_type = 'CURVE'
    bond_type = 'MESH'
    
    bond_mol = ob.mb.get_molecule()
    
    first_atom = ob.mb.bonded_atoms[0].get_object()
    second_atom = ob.mb.bonded_atoms[1].get_object()
    
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    
    first_mat = None
    second_mat = None
    
    if bond_type != 'CURVE':
        # the molecule properties of the two bonded atoms are used.
        # to use the bond_mol properties change accordingly
        if first_mol.bond_material == 'ATOMS':
            first_mat = first_atom.material_slots[0].material
        elif first_mol.bond_material == 'GENERIC':
            first_mat_name = "mat_bond_{}".format(first_mol.index)
            first_mat = bpy.data.materials.get(first_mat_name)
            if not first_mat:
                first_mat = new_material(first_mat_name, molecule=first_mol)
        
        if second_mol.bond_material == 'ATOMS':
            second_mat = second_atom.material_slots[0].material
        elif second_mol.bond_material == 'GENERIC':
            second_mat_name = "mat_bond_{}".format(second_mol.index)
            second_mat = bpy.data.materials.get(second_mat_name)
            if not second_mat:
                second_mat = new_material(second_mat_name, molecule=second_mol)
                
        # make sure to have at least two material slots
        for i in range(2 - len(ob.material_slots)):
            ob.data.materials.append(None)
        
        ob.material_slots[0].link = 'OBJECT'
        ob.material_slots[1].link = 'OBJECT'
        ob.material_slots[0].material = first_mat
        ob.material_slots[1].material = second_mat

    elif bond_type == 'CURVE':
        # the bond_mol properties are used
        if bond_mol.bond_material == 'ATOMS':
            debug_print("Atom colored bonds are not yet supported with Curve bonds.", level=2)
            bond_mol.bond_material == 'GENERIC'
            # changing bond_material will call this function again
            return
        elif bond_mol.bond_material == 'GENERIC':
            mat_name = "mat_bond_{}".format(bond_mol.index)
            mat = bpy.data.materials.get(mat_name)
            if not mat:
                mat = new_material(mat_name, molecule=bond_mol)
        
        # make sure to have at least one material slot
        if len(ob.material_slots) == 0:
            ob.data.materials.append(None)
        ob.material_slots[0].link = 'OBJECT'
        ob.material_slots[0].material = mat
        
        
        ### The following is the beginning of an attempt to create a cycles material that features both atom colors (or generic)
            #material = bpy.data.materials.new(ob.name)
            
            #if molecule.bond_material == 'GENERIC':
                #color = first_atom.mb.get_molecule().bond_color
                #material.diffuse_color = color
                #if bpy.context.scene.render.engine == 'CYCLES':
                    #material.use_nodes = True
                    ## add driver to rendered color to be the same as display color
                    #nodesocketcolor = material.node_tree.nodes['Diffuse BSDF'].inputs[0]
                    #for i in range(3): # not for alpha channel
                        #fcurve =  nodesocketcolor.driver_add('default_value', i) # add new driver
                        #drv = fcurve.driver
                        #drv.type = 'AVERAGE'
                        #drv.show_debug_info = True
                        
                        #var = drv.variables.new()
                        #var.name = 'diffuse_color' # name to use in scripting
                        #var.type = 'SINGLE_PROP'
                        #targ = var.targets[0]
                        #targ.id_type = 'MATERIAL'
                        #targ.id = material
                        #targ.data_path = 'diffuse_color[{}]'.format(i)
            #else:
                #material.use_nodes = True
                #tree = material.node_tree
                #links = tree.links
                #for n in tree.nodes:
                    #tree.nodes.remove(n)
                #coords = tree.nodes.new('ShaderNodeTexCoord')
                #coords.location = 0, 200
                #mapping = tree.nodes.new('ShaderNodeMapping')
                #mapping.location = 200, 200
                #script = tree.nodes.new('ShaderNodeScript')
                #script.location = 600, 200
                #script.mode = 'EXTERNAL'
                #script.filepath = "//bond_shader.osl"
                #script.shader_script_update()
                #output = tree.nodes.new('ShaderNodeOutputMaterial')
                #output.location = 800, 200
                
                #links.new(coords.outputs[3], mapping.inputs[0])
                #links.new(mapping.outputs[0], script.inputs[0])
                #links.new(script.outputs[0], output.inputs[0])
                
                ## add driver to mapping node to rotate texture coordinates
                
                ## add driver to colors 1 and 2

def update_atom_element(self, context):
    '''
    assign a mesh, give a new object name etc.
    '''
    # remove all spaces
    if self.element.strip() != self.element:
        self.element = self.element.strip()
        return
    # Comparison is case insensitive
    # get dictionary that maps all lower case elements to case style as it is in the list
    elements_map = dict([(element.name.lower(), element.name) for element in context.scene.mb.elements])
    # if element is not yet in scene elements list, add new element. Use default settings
    if self.element.lower() not in elements_map:
        add_element(context, self.element, ELEMENTS_DEFAULT["Default"])
    # adjust case of entered element to match already existing element
    elif self.element not in elements_map.values():
        self.element = elements_map[self.element.lower()]
        # because this assignment calls this function again, just return
        return
        
    # get object and molecule to change all element specific properties
    atom_ob = self.get_object()
    molecule = self.get_molecule()
    
    # update mesh and material
    me = get_atom_data(self.element, molecule)
    atom_ob.data = me
    assign_atom_material(atom_ob, molecule)
    
    set_atom_drivers(context, atom_ob, molecule)
    
    # update bond materials
    for bond in self.bonds:
        assign_bond_material(bond.get_object())
    
    # assign type last, to be able to check if element is newly assigned or just updated
    self.type = 'ATOM'

def update_bond_material(self, context):
    for bond in self.bonds:
        assign_bond_material(bond.get_object())

        
# TODO this could be a very central function, and needs a lot of time to think about
#def add_object_to_molecule(context, ob, molecule):
    ## check if object is already part of molecule
    #if ob.mb.id_molecule == molecule.name:
        #return
    #else:
        ## update name, delete from old molecule and add to new one
        #if ob.type == 'ATOM':
            #old_name = ob.name
        
        #assign_atom_material(ob, molecule)
        #set_atom_drivers(context, ob, molecule)
    
    
def add_atom(context, location, element, atom_name, molecule):
    # get new unique name for object
    name = "atom_{}.{}".format(molecule.index, molecule.atom_index)
    
    mesh_data = get_atom_data(element, molecule, type='MESH', refine=8)
    new_atom = bpy.data.objects.new(name, mesh_data)
    context.scene.objects.link(new_atom)
    #bpy.ops.object.add(type='MESH')
    #new_atom = context.object
    new_atom.location = location
    
    #new_atom.name = name
    # set mb properties
    new_atom.mb.name = new_atom.name
    new_atom.mb.id_molecule = molecule.name
    new_atom.mb.index = molecule.atom_index
    new_atom.mb.atom_name = atom_name
    
    # add to molecule
    molecule.atom_index += 1
    molecule.add_atom(new_atom)
    
    # parent to molecule origin
    new_atom.parent = molecule.parent.get_object()
    
    # updating the element will call update_atom_element, which assigns a mesh, and sets all the drivers
    new_atom.mb.element = element
    
    return new_atom

def add_bond(context, first_atom, second_atom, refine=8):
    if first_atom == second_atom:
        debug_print('WARNING: add_bond: first_atom == second_atom', 3)
        return None
    for b in first_atom.mb.bonds:
        for ba in b.get_object().mb.bonded_atoms:
            if ba.get_object() == second_atom:
                debug_print('WARNING: add_bond: Bond already exists', 3)
                return None
    #if first_atom.mb.id_molecule != second_atom.mb.id_molecule:
        ## TODO join molecule properties if connecting different molecules
        #pass
    # get new unique name for bond
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    name = "bond_{}.{}-{}.{}".format(first_mol.index, first_atom.mb.index, second_mol.index, second_atom.mb.index)
    
    bond_mesh = get_bond_data(type='MESH', refine=8)
    #bond_mesh = get_bond_data(type='CURVE')
    #if True:
        #bond_mesh.name = name
    new_bond = bpy.data.objects.new(name, bond_mesh)
    context.scene.objects.link(new_bond)
    new_bond.hide = (first_mol.draw_style == 'BALLS')
    
    # set mb properties
    new_bond.mb.type = 'BOND'
    new_bond.mb.name = new_bond.name
    new_bond.mb.id_molecule = first_mol.name
    first = new_bond.mb.bonded_atoms.add()
    first.name = first_atom.name
    second = new_bond.mb.bonded_atoms.add()
    second.name = second_atom.name
    
    # add bond to atoms mb props
    b = first_atom.mb.bonds.add()
    b.name = new_bond.name
    b = second_atom.mb.bonds.add()
    b.name = new_bond.name
    
    # add it to first molecule collection
    first_mol.add_bond(new_bond)
    
    # don't parent, as parenting also affects the scale
    c = new_bond.constraints.new('COPY_LOCATION')
    c.name = "parent"
    c.target = first_atom
    
    c = new_bond.constraints.new('STRETCH_TO')
    c.name = "stretch"
    c.rest_length = 1.0
    c.volume = 'NO_VOLUME'
    c.target = second_atom
    
    assign_bond_material(new_bond)
    set_bond_drivers(context, new_bond, new_bond.mb.get_molecule())
    new_bond.parent = first_mol.parent.get_object()
    
    return new_bond
    
def return_cursor_object(context, event, ray_max=10000.0, exclude=[], mb_type=''):
    """ This is a function that can be run from a modal operator
        to select the 3D object the mouse is hovered over.
    """
    # get the context arguments
    scene = context.scene
    region = context.region
    rv3d = context.region_data
    if not rv3d:
        return None
    coord = event.mouse_region_x, event.mouse_region_y
    # get the ray from the viewport and mouse
    view_vector = view3d_utils.region_2d_to_vector_3d(region, rv3d, coord)
    
    # have ray origin a little in front of actual origin, otherwise it might not get all of the objects
    ray_origin = view3d_utils.region_2d_to_origin_3d(region, rv3d, coord) + (0.5 * view_vector * ray_max)
    ray_target = ray_origin - (view_vector * ray_max)
    
    def visible_objects_and_duplis():
        """Loop over (object, matrix) pairs (mesh only)"""
        for obj in context.visible_objects:
            if obj.type == 'MESH':
                if (mb_type and mb_type == obj.mb.type) or mb_type == '':
                    yield (obj, obj.matrix_world.copy())
            #if obj.dupli_type != 'NONE':
                #obj.dupli_list_create(scene)
                #for dob in obj.dupli_list:
                    #obj_dupli = dob.object
                    #if obj_dupli.type == 'MESH':
                        #yield (obj_dupli, dob.matrix.copy())

            #obj.dupli_list_clear()

    def obj_ray_cast(obj, matrix):
        """Wrapper for ray casting that moves the ray into object space"""
        
        try:
            # get the ray relative to the object
            matrix_inv = matrix.inverted()
            ray_origin_obj = matrix_inv * ray_origin
            ray_target_obj = matrix_inv * ray_target
            
            # cast the ray
            hit, normal, face_index = obj.ray_cast(ray_origin_obj, ray_target_obj)
            
            if face_index != -1:
                return hit, normal, face_index
            else:
                return None, None, None
        except ValueError as e:
            debug_print("ERROR: {}: {}".format(obj.name, e), 1)
            return None, None, None
        finally:
            pass


    # cast rays and find the closest object
    best_length_squared = ray_max * ray_max
    best_obj = None
    for obj, matrix in visible_objects_and_duplis():
        if obj not in exclude and obj.type == 'MESH':
            if len(obj.data.vertices) >= 4:
                hit, normal, face_index = obj_ray_cast(obj, matrix)
                if hit is not None:
                    hit_world = matrix * hit
                    length_squared = (hit_world - ray_origin).length_squared
                    if length_squared < best_length_squared:
                        best_length_squared = length_squared
                        best_obj = obj
    # now we have the object under the mouse cursor,
    # we could do lots of stuff but for the example just select.
    return best_obj

def get_fixed_angle(context, first_atom, coord_3d, angle_list=[]):
    # get current vector between first_atom and the mouse pointer
    bond_vector = coord_3d - first_atom.location
    
    basis_change_matrix = context.region_data.view_matrix.to_3x3()
    basis_change_matrix.transpose()
    # bond vector in viewing coordinates
    transformed = basis_change_matrix.inverted() * bond_vector
    # projection of coordinates into screen plane
    xy_vec = Vector((transformed.x, transformed.y))
    # unambiguous angle between screen-x-axis and bond vector
    angle = math.atan2(xy_vec.y, xy_vec.x)
    
    if not angle_list:
        angles_30 = list(range(-180, 181, 30))
        angles_45 = list(range(-180, 181, 45))
        angle_list = set(angles_30 + angles_45)
    fixed_angles = sorted(list(map(math.radians, angle_list)))
    # find closest fixed angle in list
    pos = bisect_left(fixed_angles, angle)
    if pos == 0:
        fixed = fixed_angles[0]
    elif pos == len(fixed_angles):
        fixed = fixed_angles[-1]
    else:
        before = fixed_angles[pos - 1]
        after = fixed_angles[pos]
        if after - angle < angle - before:
            fixed = after
        else:
            fixed = before
    
    # get fixed angle in the screen plane with same projected length as before
    transformed.x = math.cos(fixed) * xy_vec.length
    transformed.y = math.sin(fixed) * xy_vec.length
    
    # transform back to world coordinates
    new_bond_vector = basis_change_matrix * transformed
    # calculate new coordinates of second atom
    fixed_vector = first_atom.location + new_bond_vector
    return fixed_vector

def get_fixed_length(context, first_atom, coord_3d, length):
    bond_vector = (coord_3d - first_atom.location).normalized()
    return first_atom.location + bond_vector * length

def get_fixed_geometry(context, first_atom, new_atom, coord_3d, geometry):
    '''
    use existing bonds of first_atom to calculate new possible bond vectors for new bond
    based on geometry return the position that is closest to the mouse pointer (coord_3d)
    '''
    # get current vector between first_atom and the mouse pointer
    #bond_vector = coord_3d - first_atom.location
    
    #basis_change_matrix = context.region_data.view_matrix.to_3x3()
    #basis_change_matrix.transpose()
    ## bond vector in viewing coordinates
    #transformed = basis_change_matrix.inverted() * bond_vector
    ## and projected into screen plane
    #xy_vec = Vector((transformed.x, transformed.y))
    
    bond_vecs = []
    for bond in first_atom.mb.bonds:
        bond_ob = bond.get_object()
        for atom in bond_ob.mb.bonded_atoms:
            if not atom.name == first_atom.name and not atom.name == new_atom.name:
                bond_vecs.append((atom.get_object().location - first_atom.location).normalized())
                break
    
    # existing number of bonds
    n_bonds = len(bond_vecs)
    #if n_bonds == 0:
        #return coord_3d
    
    if geometry == 'SINGLE' or n_bonds == 0:
        return get_fixed_angle(context, first_atom, coord_3d)
    
    #elif n_bonds == 0:
        #return coord_3d
    
    elif geometry == 'LINEAR':
        # matrix to change between view and world coordinates
        basis_change_matrix = context.region_data.view_matrix.to_3x3()
        basis_change_matrix.transpose()
        # get all possible angles projected into the view plane
        fixed_xy = []
        for bond_vec in bond_vecs:
            # bond vector in viewing coordinates
            transformed = basis_change_matrix.inverted() * bond_vec
            # projection of point mirrored coordinates into screen plane
            fixed_xy.append(-Vector((transformed.x, transformed.y)))
            ## unambiguous angle between screen-x-axis and bond vector
            #fixed_angles.append(math.atan2(xy_vec.y, xy_vec.x))
        
        # get current vector between first_atom and the mouse pointer
        bond_vector = coord_3d - first_atom.location
        
        # bond vector in viewing coordinates
        transformed = basis_change_matrix.inverted() * bond_vector
        # projection of coordinates into screen plane
        xy_vec = Vector((transformed.x, transformed.y))
        ## unambiguous angle between screen-x-axis and bond vector
        #angle = math.atan2(xy_vec.y, xy_vec.x)
        
        # find index of closest fixed xy in list
        pos, min_angle = min(enumerate(xy_vec.angle(other) for other in fixed_xy), key=lambda p: p[1])
        
        fixed_bond_vector = -bond_vecs[pos]
        
        # calculate new coordinates of second atom
        fixed_vector = first_atom.location + fixed_bond_vector * xy_vec.length
        return fixed_vector

    elif geometry == 'TRIGONAL':
        
        basis_change_matrix = context.region_data.view_matrix.to_3x3()
        basis_change_matrix.transpose()
        
        # get current vector between first_atom and the mouse pointer
        bond_vector = coord_3d - first_atom.location
        
        # bond vector in viewing coordinates
        transformed = basis_change_matrix.inverted() * bond_vector
        # projection of coordinates into screen plane
        xy_vec = Vector((transformed.x, transformed.y))
        
        if n_bonds >= 2:
            #TODO only gets one vector!
            # get the bisecting vector on the larger angle side of any two existing bonds
            fixed_xy = []
            bisect_vectors = []
            for i, b1 in enumerate(bond_vecs[:-1]):
                for b2 in bond_vecs[i+1:]:
                    bisect_vector = -(b1 + b2).normalized()
                    if (b1+b2).cross(b1).length < 1e-06:
                        # if the two existing bonds are exactly linear, the new vector will have zero length
                        # and the new bond can lie on a circle around the middle atom.
                        # return the two fixed angles that are orthogonal to the two bonds and in the view plane
                        view = basis_change_matrix.col[2]
                        bisect_vector = b1.cross(view).normalized()
                        # need to append an extra bisecting vector
                        bisect_vectors.append(bisect_vector)
                        transformed = basis_change_matrix.inverted() * bisect_vector
                        fixed_xy.append(Vector((transformed.x, transformed.y)))
                        # ... and the other one
                        bisect_vector = -bisect_vector
                    
                    bisect_vectors.append(bisect_vector)
                    transformed = basis_change_matrix.inverted() * bisect_vector
                    fixed_xy.append(Vector((transformed.x, transformed.y)))
            
            # find index of closest fixed xy in list
            pos, min_angle = min(enumerate(xy_vec.angle(other) for other in fixed_xy), key=lambda p: p[1])
            
            fixed_bond_vector = bisect_vectors[pos]
            
            # calculate new coordinates of second atom
            fixed_vector = first_atom.location + fixed_bond_vector * xy_vec.length
            return fixed_vector
        
        elif n_bonds == 1:
            # convert bond_vec to polar coordinates
            x, y, z = bond_vecs[0].xyz
            r = math.sqrt(x*x + y*y + z*z)
            theta = math.acos(z/r) # 0 <= theta <= pi
            phi = math.atan2(y, x) # -pi < phi <= pi
            
            fixed_xy = []
            bisect_vectors = []
            # increase and decrease theta by 30
            for angle in (-60, 60):
                new_theta = theta + angle * math.pi/180
                # convert back to cartesian
                x = r * math.sin(new_theta) * math.cos(phi)
                y = r * math.sin(new_theta) * math.sin(phi)
                z = r * math.cos(new_theta)
                #print(angle, x)
                bisect_vector = -Vector((x, y, z))
                bisect_vectors.append(bisect_vector)
                transformed = basis_change_matrix.inverted() * bisect_vector
                fixed_xy.append(Vector((transformed.x, transformed.y)))
            
            # find index of closest fixed xy in list
            pos, min_angle = min(enumerate(xy_vec.angle(other) for other in fixed_xy), key=lambda p: p[1])
            
            fixed_bond_vector = bisect_vectors[pos]
            
            # calculate new coordinates of second atom
            fixed_vector = first_atom.location + fixed_bond_vector * xy_vec.length
            return fixed_vector
    #elif geometry == 'OCTAHEDRAL':
        ##TODO cover all possible 90 degree angles and check against other existing bonds
        ## convert bond_vec to polar coordinates
        #x, y, z = bond_vecs[0].xyz
        #r = math.sqrt(x*x + y*y + z*z)
        #theta = math.acos(z/r) # 0 <= theta <= pi
        #phi = math.atan2(y, x) # -pi < phi <= pi
        
        #fixed_xy = []
        #bisect_vectors = []
        ## increase and decrease theta by 30
        #for angle in (-90, 0, 90):
            #new_theta = theta + angle * math.pi/180
            ## convert back to cartesian
            #x = r * math.sin(new_theta) * math.cos(phi)
            #y = r * math.sin(new_theta) * math.sin(phi)
            #z = r * math.cos(new_theta)
            ##print(angle, x)
            #bisect_vector = -Vector((x, y, z))
            #bisect_vectors.append(bisect_vector)
            #transformed = basis_change_matrix.inverted() * bisect_vector
            #fixed_xy.append(Vector((transformed.x, transformed.y)))
        
        ## find index of closest fixed xy in list
        #pos, min_angle = min(enumerate(xy_vec.angle(other) for other in fixed_xy), key=lambda p: p[1])
        
        #fixed_bond_vector = bisect_vectors[pos]
        
        ## calculate new coordinates of second atom
        #fixed_vector = first_atom.location + fixed_bond_vector * xy_vec.length
        #return fixed_vector
    # TODO implement all other geometries
    
    else:
        return coord_3d
    
def check_bond_dimension(new_bond):
    if new_bond.dimensions.x < 0.0001: #< new_bond.mb.get_molecule().bond_radius:
        toggle = {'PLANE_X': 'PLANE_Z', 'PLANE_Z': 'PLANE_X'}
        c = new_bond.constraints["stretch"]
        c.keep_axis = toggle[c.keep_axis]
    
def update_radius_type(self, context):
    for atom in self.atoms:
        set_atom_drivers(context, atom.get_object(), self)

def set_draw_style(self, context):
    for atom in self.atoms:
        set_atom_drivers(context, atom.get_object(), self)
    
    hide = (self.draw_style == 'BALLS')
    for bond in self.bonds:
        bond_ob = bond.get_object()
        bond_ob.hide = hide

#******************************************************************************#
#*** Import Export utils ******************************************************#
#******************************************************************************#

# filepath_xyz: path to xyz file
def read_xyz_file(filepath_xyz, scale_distances=1.0, mask_planes=[], mask_flip=False):
    '''
    read xyz file and return all atoms in the format
    [{index1: [element, atom_name, coords], index2: [...], ...}, [second frame], ...]
    '''
    number_frames = 0
    all_frames = []
    # Open the file ...
    with open(filepath_xyz, "r") as fin:
        
        #Go through the whole file.
        number_atoms = -1
        for line in fin:
            
            if line == "":
                continue
            
            split_list = line.rsplit()
            
            if len(split_list) == 1:
                try:
                    number_atoms = int(split_list[0])
                except ValueError:
                    pass
            if number_atoms > 0:
                # dump comment line
                line = fin.readline()
                
                all_atoms = {}
                for i in range(number_atoms):
                    
                    line = fin.readline()
                    line = line.strip()
                    split_list = line.split()
                    
                    location = Vector(list(map(float, split_list[1:4]))) * scale_distances
                    # skip processing this coordinate if outside of mask
                    if mask_planes and not is_inside_of_planes(mask_planes, location, flip=mask_flip):
                        continue
                    
                    element = split_list[0]
                    atom_name = element + str(i)
                    # If element is an 'X' then it is a vacancy.
                    if "X" in element:
                        element = "Vac"
                        atom_name = "Vacancy"
                    
                    all_atoms[i] = [element, atom_name, location, i]
                
                all_frames.append(all_atoms)
                number_frames += 1
                number_atoms = -1
            
    
    return all_frames

def read_pdb_file(filepath_pdb, scale_distances=1.0, mask_planes=[], mask_flip=False):
    '''
    read xyz file and return all atoms in the format
    [{index1: [element, atom_name, coords], index2: [...], ...}, [second frame], ...]
    and bonds (index1 < index2 < index 3 < ...)
    {index1: [index2, index3], index2: [index4, index5], ...}
    '''
    
    number_frames = 0
    
    all_frames = []
    all_atoms = {}
    bonds = {}
    # TODO include as GUI option
    # if double == False: exclude double bonds
    double = False
    with open(filepath_pdb, "r") as fin:
        i = -1
        for line in fin:
            
            if line == "":
                continue
            
            # get atom information
            if line[:6] == 'HETATM' or line[:4] == 'ATOM':
                i += 1
                location = Vector(list(map(float, (line[30:38], line[38:46], line[46:54])))) * scale_distances
                #if len(location) != 3:
                    #print(location)
                # skip processing this coordinate if outside of mask
                if mask_planes and not is_inside_of_planes(mask_planes, location, flip=mask_flip):
                    continue
                
                # get atom number
                index = int(line[6:11])
                # get atom type, remove whitespace and capitalize first letter
                element = line[76:78].strip().capitalize()
                element = ''.join(i for i in element if not i.isdigit() and not i in ('+', '-'))
                #charge = line[78:80].strip()
                atom_name = line[12:16].strip()
                if "X" in element:
                    element = "Vac"
                    atom_name = "Vacancy"
                
                all_atoms[index] = [element, atom_name, location, i]
            
            # get bond information
            elif line[:6] == 'CONECT':
                # make sure index is greater than 0
                if int(line[6:11]) > 0:
                    # base atom
                    atomID1 = int(line[6:11])
                    if not atomID1 in bonds:
                        bonds[atomID1] = []
                    # loop over conect entries
                    for i in range((len(line) - 11) // 5): # // floor division: 1//2=0 
                        # get connected atomID
                        atomID2 = int(line[11 + i*5 : 16 + i*5])
                        # only store bond once (don't have a 1-2 and a 2-1 bond)
                        if atomID2 > 0 and atomID2 > atomID1:
                            if double or not atomID2 in bonds[atomID1]:
                                # append to list
                                bonds[atomID1].append(atomID2)
            
            # if the model ends
            elif line[:6] == 'ENDMDL':
                debug_print("ENDMDL found.", level=3)
                
                all_frames.append(all_atoms)
                number_frames += 1
                all_atoms = {}
        
        # if there was only one model and no ENDMDL present
        debug_print("Finished reading file. Cleaning up.", level=3, end='')
        if all_atoms:
            all_frames.append(all_atoms)
        
    return all_frames, bonds

def read_qe_file(filepath, scale_distances=1.0, mask_planes=[], mask_flip=False):
    
    all_atoms = {}
    with open(filepath, 'r') as fin:
        for line in fin:
            if "nat" in line:
                n_atoms = int(line.split('=')[-1].strip().strip(","))
                break
        unit_vectors = []
        for line in fin:
            if "ATOMIC_POSITIONS" in line:
                # double check units
                if "angstrom" in line.lower() and scale_distances != '1.0':
                    scale_distances = 1.0
                    debug_print("Found coordinates in Angstrom. Overriding unit.", level=1)
                elif "bohr" in line.lower() and scale_distances != '0.529177249':
                    scale_distances = 0.529177249
                    debug_print("Found coordinates in Bohr. Overriding unit.", level=1)
                elif "crystal" in line.lower():
                    # get unit vectors
                    fin.seek(0)
                    for line in fin:
                        if "CELL_PARAMETERS" in line:
                            # TODO include alat units etc.
                            unit = 0.529177249 if "bohr" in line else 1.0
                            for i in range(3):
                                line = fin.readline()
                                unit_vectors.append(Vector(list(map(float, line.split()))) * unit)
                            break
                    fin.seek(0)
                    for line in fin:
                        if "ATOMIC_POSITIONS" in line:
                            break
                else:
                    # TODO include more checks (alat etc.)
                    pass
                break
        for i, line in enumerate(fin):
            if i < n_atoms:
                line = line.strip()
                split_list = line.split()
                
                location = Vector(list(map(float, split_list[1:4])))
                if unit_vectors:
                    location = (location[0] * unit_vectors[0] +
                                location[1] * unit_vectors[1] +
                                location[2] * unit_vectors[2])
                else:
                    location *= scale_distances
                # skip processing this coordinate if outside of mask
                if mask_planes and not is_inside_of_planes(mask_planes, location, flip=mask_flip):
                    continue
                    
                element = split_list[0]
                atom_name = element + str(i)
                # If element is an 'X' then it is a vacancy.
                if "X" in element:
                    element = "Vac"
                    atom_name = "Vacancy"
                
                all_atoms[i] = [element, atom_name, location, i]
            else:
                break
    return [all_atoms]

def read_qe_unit_cell(filepath):
    unit_vectors = []
    with open(filepath, 'r') as fin:
        for line in fin:
            if "CELL_PARAMETERS" in line:
                # TODO include alat units etc.
                unit = 0.529177249 if "bohr" in line else 1.0
                
                for i in range(3):
                    line = fin.readline()
                    unit_vectors.append(Vector(list(map(float, line.split()))) * unit)
                break
    return unit_vectors

def draw_unit_cell(unit_vectors, draw_style='ARROWS'):
    # TODO implement different drawing styles
    #if 'CUBE' in draw_style:
    bpy.ops.mesh.primitive_cube_add(location=(0,0,0))
    uc_cube = bpy.context.object
    uc_cube.name = "unit_cell"
    bm = bmesh.new()
    bm.from_mesh(uc_cube.data)
    verts = bm.verts
    verts[0].co = (0,0,0)
    verts[1].co = unit_vectors[0]
    verts[2].co = unit_vectors[0] + unit_vectors[1]
    verts[3].co = unit_vectors[1]
    verts[4].co = unit_vectors[2]
    verts[5].co = unit_vectors[0] + unit_vectors[2]
    verts[7].co = unit_vectors[1] + unit_vectors[2]
    verts[6].co = unit_vectors[0] + unit_vectors[1] + unit_vectors[2]
    bm.to_mesh(uc_cube.data)
    bm.free()
    vg = []
    vg.append(uc_cube.vertex_groups.new('a'))
    vg[-1].add([1], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('b'))
    vg[-1].add([3], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('c'))
    vg[-1].add([4], 1, 'REPLACE')
    
    if 'ARROWS' in draw_style:
        radius = 0.1
        ring_y = 0.9
        ring_scale = 2
        
        # get material
        material = bpy.data.materials.get('axes')
        if not material:
            material = new_material('axes', color=(1,0,0))
        
        # add sphere at origin
        bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0), size=radius, segments=8, ring_count=8)
        ob = bpy.context.object
        ob.name = "unit_cell_origin"
        ob.parent = uc_cube
        ob.parent_type = 'VERTEX'
        ob.data.materials.append(None)
        ob.material_slots[0].material = material
        
        bpy.ops.mesh.primitive_cylinder_add(location=(0,0,0), radius=radius, vertices=8, 
                depth=1, end_fill_type="TRIFAN")
        ob = bpy.context.object
        ob.data.materials.append(None)
        ob.material_slots[0].material = material
        
        # convert cylinder to arrow
        bm = bmesh.new()
        bm.from_mesh(ob.data)
        
        # rotate and shrink first, then add another row of vertices
        for vert in bm.verts:
            # rotate 90 degrees around x, and shift along y axis
            tmp_co = vert.co.copy()
            vert.co.y = -tmp_co.z + .5
            vert.co.z = -tmp_co.y
            if vert.co.y > 0.01:
                vert.select = False
        new_verts = []
        for edge in bm.edges:
            if edge.calc_length() == 1.0:
                e, v = bmesh.utils.edge_split(edge, edge.verts[0], 0.5)
                new_verts.append(v)
        n_verts = len(new_verts)
        #for i in range(n_verts):
            #bm.edges.new((new_verts[i], new_verts[(i+1)%n_verts]))
        
        # little hacky, but don't understand how bmesh.utils.face_split works
        # remove all faces
        for f in bm.faces:
            if len(f.verts) == 6:
                bm.faces.remove(f)
        
        # now sort bm.verts
        # v.co.y is either 0, 0.5, or 1.0. So multiply y with at least 4pi to sort by y value first
        key = lambda v: v.co.y * 15 + math.atan2(v.co.x, v.co.z)
        verts_sorted = sorted((v for v in bm.verts if (0 < v.co.length and v.co.length != 1.0)), key=key)
        #for i, v in enumerate(verts_sorted):
        #    v.index = i
        
        # add new faces and assign the two different material slots
        for i in range(2*n_verts):
            v1 = verts_sorted[i]
            v2 = verts_sorted[(i + 1)%n_verts + n_verts*(i//n_verts)]
            v3 = verts_sorted[(i + 1)%n_verts + n_verts*(i//n_verts + 1)]
            v4 = verts_sorted[i + n_verts]
            f = bm.faces.new((v1, v2, v3, v4))
        
        # now shape the arrow head
        for vert in bm.verts:
            if vert.co.y == 1.0 and not vert.co.length == 1.0:
                vert.co.y = ring_y
                vert.co.x = vert.co.x * ring_scale
                vert.co.z = vert.co.z * ring_scale
            elif vert.co.y == 0.5:
                vert.co.y = ring_y
        bm.to_mesh(ob.data)
        bm.free()
        arrow_mesh = ob.data
        bpy.context.scene.objects.unlink(ob)
        
        for di, vec in enumerate(unit_vectors):
            ob = bpy.data.objects.new({0: 'a', 1: 'b', 2: 'c'}[di], arrow_mesh)
            ob.parent = uc_cube
            ob.parent_type = 'VERTEX'
            bpy.context.scene.objects.link(ob)
            bpy.context.scene.objects.active = ob
            bpy.ops.object.constraint_add(type='STRETCH_TO')
            c = ob.constraints[-1]
            c.name = "stretch"
            c.rest_length = 1.0
            c.volume = 'NO_VOLUME'
            c.target = uc_cube
            c.subtarget = vg[di].name
            
            
    return unit_vectors

def read_modes(filepath):
    # read mode file, modes need to be in same order as atoms in input file
    # currently only supports dynmat.out
    with open(filepath, 'r') as fin:
        all_evecs = []
        freqs = []
        for line in fin:
            lstrip = line.strip()
            # new mode
            if lstrip.startswith('omega('):
                m = re.search('omega\(([ 0-9]+)\).+ ([-.0-9]+)(?= \[cm-1\])', lstrip)
                i = int(m.group(1))
                freq = float(m.group(2))
                freqs.append(freq)
                current = []
                all_evecs.append(current) # links current to all_evecs
            elif lstrip.startswith('('):
                lsplit = lstrip[1:-1].split()
                current.append(list(map(float, lsplit[::2])))
    return freqs, all_evecs

def find_bonds(all_atoms, bonds={}):
    '''
    loops through all atoms and returns a list of bonds between them in the format
    {index1: [index2, index3], index2: [index4, index5], ...}
    If bond argument is given, additional bonds are added
    TODO for now, this only checks in the unit of Angstroms!
    '''
    for index1, data1 in all_atoms.items():
        for index2, data2 in all_atoms.items():
            if index1 < index2:
                distance = (data2[2] - data1[2]).length
                # get covalent radii
                try:
                    cov1 = ELEMENTS_DEFAULT[data1[0]]['covalent']
                except KeyError:
                    cov1 = 1
                try:
                    cov2 = ELEMENTS_DEFAULT[data2[0]]['covalent']
                except KeyError:
                    cov2 = 1
                # the bond length should just be the sum of covalent radii (plus some lee way)
                max_dist = cov1 + cov2 + 0.2 # TODO this value might have to be adjusted
                if distance < max_dist:
                    # append the bond
                    try:
                        if index2 not in bonds[index1]:
                            bonds[index1].append(index2)
                    except KeyError:
                        bonds[index1] = []
                        bonds[index1].append(index2)
    return bonds

def import_molecule(filepath,
                    modepath,
                    molecule,
                    #ball,
                    refine,
                    scale_distances,
                    #stick,
                    #bond_sectors,
                    bond_guess,
                    use_center,
                    #use_center_all,
                    #use_camera,
                    #use_lamp,
                    mask_planes,
                    mask_flip,
                    supercell,
                    ):
    
    start = time.time()
    
    all_frames = []
    bonds = {}
    axes = []
    
    if modepath:
        frequencies, modes = read_modes(modepath)
        if not modes:
            debug_print("WARNING: Couldn't read any normal modes in file {}".format(modepath), 1)
    
    if filepath.rsplit('.')[-1] == 'xyz':
        all_frames = read_xyz_file(filepath, scale_distances, mask_planes, mask_flip)
        bonds = {}
    elif filepath.rsplit('.')[-1] == 'pdb':
        all_frames, bonds = read_pdb_file(filepath, scale_distances, mask_planes, mask_flip)
    #elif filepath.rsplit('.')[-1] in ('axsf', 'xsf'):
        #all_frames, axes = read_xsf_file(filepath, scale_distances, mask_planes, mask_flip)
    else:
        # read first lines of file to determine file format
        with open(filepath, 'r') as fin:
            lines = fin.read()
            if "&control" in lines:
                all_frames = read_qe_file(filepath, scale_distances, mask_planes, mask_flip)
                bonds = {}
                # read unit cell and create cube
                axes = read_qe_unit_cell(filepath)
                draw_unit_cell(axes)
    # replace the index i with the corresponding mode
    
    if modepath and modes:
        for all_atoms in all_frames:
            for index, data in all_atoms.items():
                #try:
                    data[-1] = [mode[data[-1]] for mode in modes]
                #except (IndexError, TypeError):
                    #debug_print("WARNING: Modes couldn't be matched with atoms.", 1)
                    #modes = None
                    #break
            if not modes:
                break
    
    if sum(supercell) > 3 and axes:
        debug_print("Creating supercell.", level=3)
        for frame in range(len(all_frames)):
            max_index = max(all_frames[frame]) + 1
            n_unit_cell = 0
            all_atoms = {}
            for i in range(supercell[0]):
                for j in range(supercell[1]):
                    for k in range(supercell[2]):
                        for index, (element, atom_name, location, modes) in sorted(all_frames[frame].items()):
                            location = location.copy() + i*axes[0] + j*axes[1] + k*axes[2]
                            all_atoms[index + max_index*n_unit_cell] = [element, atom_name, location, modes]
                        n_unit_cell += 1
            all_frames[frame] = all_atoms
    
    debug_print("read", 4)
    debug_print(time.time() - start, 4)
    
    if bond_guess:
        find_bonds(all_frames[0], bonds)
        debug_print("guess", 4)
        debug_print(time.time() - start, 4)
    
    origin = Vector((0.0,0.0,0.0))
    center_of_mass = sum([data[2] for data in all_frames[0].values()], origin) / len(all_frames[0])
    # set molecule parent in center of mass
    molecule.parent.get_object().location = center_of_mass
    debug_print("center", 4)
    debug_print(time.time() - start, 4)
    
    # add all atoms to scene
    atom_obs = {}
    error = set()
    debug_print("", 4)
    for index, data in all_frames[0].items():
        debug_print("\ratom {}".format(index), level=4, end='')
        new_atom = add_atom(bpy.context, data[2]-center_of_mass, data[0], data[1], molecule)
        
        # adjust index to be the same as in imported file
        if new_atom.mb.index < index:
            new_atom.mb.index = index
            molecule.atom_index = index + 1
        elif new_atom.mb.index > index:
            error.add("WARNING: Indeces will not be the same as imported.")
        
        if modes:
            atom_id = '{}.{}'.format(new_atom.mb.get_molecule().index, new_atom.mb.index)
            # first create default action where the atom is in it's rest position
            anim_data = new_atom.animation_data_create()
            anim_data.action = bpy.data.actions.new(name="equilibrium_{}".format(atom_id))
            anim_data.action.use_fake_user = True
            for dim in range(3):
                fcu = anim_data.action.fcurves.new(data_path="location", index=dim)
                fcu.keyframe_points.add(1)
                fcu.keyframe_points[0].co = 1.0, new_atom.location[dim]
            # then add one action per mode
            for i, mode_vec in enumerate(data[3]):
                action = bpy.data.actions.new(name="mode_{}_{}".format(i, atom_id))
                action.use_fake_user = True
                for dim in range(3):
                    fcu = action.fcurves.new(data_path="location", index=dim)
                    fcu.keyframe_points.add(3)
                    fcu.keyframe_points[0].co = 1.0, new_atom.location[dim]+mode_vec[dim]
                    fcu.keyframe_points[0].interpolation = 'SINE'
                    fcu.keyframe_points[1].co = 11.0, new_atom.location[dim]-mode_vec[dim]
                    fcu.keyframe_points[1].interpolation = 'SINE'
                    fcu.keyframe_points[2].co = 21.0, new_atom.location[dim]+mode_vec[dim]
                    fcu.keyframe_points[2].interpolation = 'SINE'
                    fcu.update()
                    
        atom_obs[index] = new_atom
    debug_print("", 4)

    debug_print("atoms", 4)
    debug_print(time.time() - start, 4)
    
    # add bonds to scene
    debug_print("", 4)
    
    for index1, other in bonds.items():
        for index2 in other:
            debug_print("\rbond {}-{}".format(index1, index2), level=4, end='')
            new_bond = add_bond(bpy.context, atom_obs[index1], atom_obs[index2])
    debug_print("", 4)
    
    debug_print("bonds", 4)
    debug_print(time.time() - start, 4)
    
    if error:
        debug_print('\n'.join(error), 1)
    
    if use_center:
        molecule.parent.get_object().location -= center_of_mass

def export_xyz(filepath, selection_only, scale_distances):

    all_atoms = {}
    if selection_only:
        objects = bpy.context.selected_objects
    else:
        objects = bpy.context.scene.objects
    
    n_atoms = 0
    for ob in objects:
        if ob.mb.type == 'ATOM':
            try:
                all_atoms[ob.mb.id_molecule].append((ob.mb.index, ob.mb.element, ob.location))
            except KeyError:
                all_atoms[ob.mb.id_molecule] = [(ob.mb.index, ob.mb.element, ob.location)]
            n_atoms += 1
        
    with open(filepath, "w") as fout:
        #fout.write("REMARK This pdb file has been created with Blender "
                     #"and the addon MolBlend by Florian Altvater\n"
                     #"REMARK\n"
                     #"REMARK\n")
        fout.write(n_atoms)
        fout.write("This file has been created with Blender and the MolBlend addon")
    
        for mol_id in sorted(all_atoms):
            for i, element, location in sorted(all_atoms[mol_id]):
                fout.write("{}   {:>10.6f}   {:>10.6f}   {:>10.6f}\n".format(element, *location))
    return True

def export_pdb(filepath, selection_only, scale_distances):

    all_atoms = {}
    if selection_only:
        objects = bpy.context.selected_objects
    else:
        objects = bpy.context.scene.objects
    
    n_atoms = 0
    for ob in objects:
        if ob.mb.type == 'ATOM':
            try:
                all_atoms[ob.mb.id_molecule].append((
                    ob.mb.index,
                    ob.mb.atom_name if len(ob.mb.element) != 1 else " "+ob.mb.atom_name,
                    ob.mb.get_molecule().name[:3].upper(),
                    ob.location,
                    ob.mb.ob.mb.element
                    ))
            except KeyError:
                all_atoms[ob.mb.id_molecule] = [(
                    ob.mb.index,
                    ob.mb.atom_name if len(ob.mb.element) != 1 else " "+ob.mb.atom_name,
                    ob.mb.get_molecule().name[:3].upper(),
                    ob.location,
                    ob.mb.ob.mb.element
                    )]
            n_atoms += 1
        
    with open(filepath, "w") as fout:
        fout.write("REMARK This pdb file has been created with Blender "
                  "and the addon MolBlend by Florian Altvater\n"
                  "REMARK\n"
                  "REMARK\n")
        index = 1
        for mol_id in sorted(all_atoms):
            for i, name, res_name, coords, element in sorted(all_atoms[mol_id]):
                fout.w('HETATM{ID:>5} {name:<4} {res}     1    {c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00          {element:>2}\n'.format(ID=index, name=name[:4], res=res_name, c=coords, element=element))
                index += 1
                if index > 99999:
                    index = 1
    return True

# import unit cell
#import bpy
#import bmesh
#from mathutils import Vector

#molecule = "anthracene"

#filename = "/home/flo/Eigene/a_research/calculations/mobilities/structures/{}/in".format(molecule)

#unit_vectors = []
#with open(filename, 'r') as fin:
    #for line in fin:
        #if "CELL_PARAMETERS" in line:
            #unit = 0.529177249 if "bohr" in line else 1
            
            #for i in range(3):
                #line = fin.readline()
                #unit_vectors.append(Vector(list(map(float, line.split())))*unit)
            #break

#bpy.ops.mesh.primitive_cube_add()
#uc = bpy.context.object
#bm = bmesh.new()
#bm.from_mesh(uc.data)
#verts = bm.verts
#verts[0].co = (0,0,0)
#verts[1].co = unit_vectors[0]
#verts[2].co = unit_vectors[0] + unit_vectors[1]
#verts[3].co = unit_vectors[1]
#verts[4].co = unit_vectors[2]
#verts[5].co = unit_vectors[0] + unit_vectors[2]
#verts[7].co = unit_vectors[1] + unit_vectors[2]
#verts[6].co = unit_vectors[0] + unit_vectors[1] + unit_vectors[2]
#bm.to_mesh(uc.data)
#bm.free()

#--------------------------------------------------------------
#import bpy
#import time

#bpy.ops.mesh.primitive_uv_sphere_add()
#mesh = bpy.context.object.data

#start = time.time()
#with open('test_blender', 'w') as fin:
    #for i in range(1000):
        #last = time.time()
        #print('\r{}'.format(i), end='')
        
        ##bpy.ops.object.add(type='MESH')
        ##new_atom = bpy.context.object
        ##new_atom.name = 'test{}'.format(i)
        ##new_atom.data = mesh
        
        #new_atom = bpy.data.objects.new('test{}'.format(i), mesh)
        #bpy.context.scene.objects.link(new_atom)
        #fin.write('{}\t{}\n'.format(i, time.time()-last))
#print()
#print(time.time()-start)