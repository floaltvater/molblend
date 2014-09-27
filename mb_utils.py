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
    file_types = [('XYZ', "xyz", "xyz format"),
                  ('PDB', "pdb", "Protein Databank format"),
                  ('B4W', "b4w", "standalone HTML with embedded 3D viewer")
                 ]

#--- Update functions ---------------------------------------------------------#

def update_all_meshes(self, context):
    debug_print("mb_utils.update_all_meshes", 6)
    # TODO this callback might be too heavy for scenes with lots of meshes
    for me in bpy.data.meshes:
        me.update()

def update_active_mode(self, context):
    debug_print("mb_utils.update_active_mode", 6)
    if self.max_mode == 0:
        self.active_mode = 0
        return
    elif self.active_mode > self.max_mode:
        self.active_mode = self.max_mode
        return
    
    for atom in self.objects.atoms:
        atom_ob = atom.get_object()
        atom_id = '{}.{}'.format(self.index, atom_ob.mb.index)
        anim_data = atom_ob.animation_data
        if anim_data:
            action = anim_data.action
            if action:
                atom_vec = atom_ob.mb.modes[self.active_mode].vec * self.mode_scale
                for dim in range(3):
                    fcu = action.fcurves[dim]
                    
                    kf1 = fcu.keyframe_points[0]
                    kf2 = fcu.keyframe_points[1]
                    #kf3 = fcu.keyframe_points[2]
                    middle = (kf1.co[1] + kf2.co[1]) / 2.0
                    
                    for p in range(3):
                        loc = middle + pow(-1, p) * atom_vec[dim]
                        fcu.keyframe_points[p].co = 1.0 + 10*p, loc
                        fcu.keyframe_points[p].interpolation = 'BEZIER'
                    fcu.update()
            else:
                debug_print("No action for atom '{}'.".format(atom_id), 1)
        else:
            debug_print("No animation data for atom '{}'.".format(atom_id), 1)
    
    if self.active_mode == 0:
        # stop animation
        if context.screen.is_animation_playing:
            bpy.ops.screen.animation_play()
            context.scene.frame_current = 1
    else:
        # start animation
        if not context.screen.is_animation_playing:
            context.scene.frame_end = 20
            bpy.ops.screen.animation_play()

###########################################
### old update_active_mode

    #debug_print("mb_utils.update_active_mode", 6)
    #if self.max_mode == 0:
        #self.active_mode = 0
        #return

    #for atom in self.objects.atoms:
        #atom_ob = atom.get_object()
        
        #if self.active_mode > self.max_mode:
            #self.active_mode = self.max_mode
            #return
        #elif self.active_mode == 0:
            #action_name = "equilibrium_{}".format(atom_id)
        #else:
            #action_name = "mode_{}_{}".format(self.active_mode, atom_id)
        
        #action = bpy.data.actions.get(action_name)
        #anim_data = atom_ob.animation_data
        #if action and anim_data:
            #anim_data.action = action
            #update_mode_scale(self, context)
        #elif not action:
            #debug_print("Mode {} not found for atom '{}'. ({})".format(self.active_mode, atom_id, action_name), 1)
        #elif not anim_data:
            #debug_print("No animation data for atom '{}'.".format(atom_id), 1)
        
    #if self.active_mode == 0:
        ## stop animation
        #if context.screen.is_animation_playing:
            #bpy.ops.screen.animation_play()
            #context.scene.frame_current = 1
    #else:
        ## start animation
        #if not context.screen.is_animation_playing:
            #context.scene.frame_end = 20
            #bpy.ops.screen.animation_play()

##################################################
#def update_mode_scale(self, context):
    #debug_print("mb_utils.update_mode_scale", 6)
    #for atom in self.objects.atoms:
        #atom_ob = atom.get_object()
        #try:
            #action = atom_ob.animation_data.action
        #except AttributeError:
            #debug_print("WARNING: No animation data when updating mode scale.", 3)
            #return
        #mode_vec = action.mb.mode_vector
        #for dim, fcu in enumerate(action.fcurves):
            #if len(fcu.keyframe_points) > 2:
                #kf1 = fcu.keyframe_points[0]
                #kf2 = fcu.keyframe_points[1]
                #kf3 = fcu.keyframe_points[2]
                #middle = (kf1.co[1] + kf2.co[1]) / 2.0
                #scaled = self.mode_scale * mode_vec[dim]
                #kf1.co[1] = middle + scaled
                #kf2.co[1] = middle - scaled
                #kf3.co[1] = middle + scaled
                #fcu.update()

def update_atom_element(self, context):
    debug_print("mb_utils.update_atom_element", 6)
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
    debug_print("mb_utils.update_bond_material", 6)
    for bond in self.objects.bonds:
        assign_bond_material(bond.get_object())

def update_refine_atoms(self, context):
    debug_print("mb_utils.update_refine_atoms", 6)
    if self.refine_atoms < 3:
        self.refine_atoms = 3
    replaced_elements = set()
    for mesh in self.objects.meshes:
        if "atom_mesh" in mesh.name:
            element = mesh.name.split(".")[-1] #mesh.name = "atom_mesh_{}.{}".format(molecule.index, element)
            if not element in replaced_elements:
                data = mesh.get_data()
                # get new temporary atom mesh with new refine value
                new_data = get_atom_data(element, self, type='MESH', name="tmp_mesh")
                # replace mesh data
                bm = bmesh.new()
                bm.from_mesh(new_data)
                bm.to_mesh(data)
                bm.free()
                replaced_elements.add(element)
                # delete temporary mesh
                bpy.data.meshes.remove(new_data)

def update_refine_bonds(self, context):
    debug_print("mb_utils.update_refine_bonds", 6)
    if self.refine_bonds < 3:
        self.refine_bonds = 3
    mesh = self.objects.meshes.get("bond_mesh_{}".format(self.index))
    if mesh:
        data = mesh.get_data()
        # get new temporary bond mesh with new refine value
        name = "tmp_mesh"
        new_data = get_bond_data(self, type='MESH', name="tmp_mesh")
        # replace mesh data
        bm = bmesh.new()
        bm.from_mesh(new_data)
        bm.to_mesh(data)
        bm.free()
        # delete temporary mesh
        bpy.data.meshes.remove(new_data)

def update_export_file_type(self, context):
    debug_print("mb_utils.update_export_file_type", 6)
    """
    Change file extension when filetype is changed. Replace known extensions.
    Append to everything else.
    """
    if self.filepath:
        
        #bpy.path.ensure_ext(filepath, ext, case_sensitive=False)
        filetypes = {'XYZ': ".xyz",
                     'PDB': ".pdb"}
        ext = filetypes[self.file_type]
        other_ext = [f for f in filetypes.values() if f != ext]
        
        if self.filepath[-4:] in other_ext:
            self.filepath = self.filepath[:-4] + ext
        else:
            self.filepath = self.filepath + ext

def update_molecule_selection(self, context):
    debug_print("mb_utils.update_molecule_selection", 6)
    mol = context.scene.mb.molecules.get(self.molecule_id)
    if mol:
        for ob in context.selected_objects:
            ob.select = False
        for col in (mol.objects.atoms, mol.objects.bonds, mol.objects.other):
            for ob in col:
                ob.get_object().select = True
        parent = mol.objects.parent.get_object()
        parent.select = True
        context.scene.objects.active = parent
#--- General functions --------------------------------------------------------#

def create_new_keyconfig(name, context):
    debug_print("mb_utils.create_new_keyconfig", 6)
    
    if not name in context.window_manager.keyconfigs:
        debug_print('Creating new keyconfiguration "{}".'.format(name), 3)
        from molblend import mb_keyconfig
        #keyconfig = mb_keyconfig.kc
    else:
        debug_print('Keyconfiguration "{}" already exists.'.format(name), 3)
        
    #print(keyconfig)
    return context.window_manager.keyconfigs.get("MolBlend")

def add_element(context, element, element_dict):
    debug_print("mb_utils.add_element", 6)
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
    debug_print("mb_utils.initialize_elements", 6)
    for element, data in ELEMENTS_DEFAULT.items():
        add_element(context, element, data)


#--- Viewport functions -------------------------------------------------------#

def mouse_2d_to_location_3d(context, coord, depth=Vector((0, 0, 0))):
    debug_print("mb_utils.mouse_2d_to_location_3d", 6)
    region = context.region
    rv3d = context.region_data
    if depth:
        depth_location = depth
    else:
        depth_location = context.scene.cursor_location.copy()
    
    return view3d_utils.region_2d_to_location_3d(region, rv3d, coord, depth_location)

def return_cursor_object(context, event, ray_max=10000.0, exclude=None, mb_type=''):
    debug_print("mb_utils.return_cursor_object", 6)
    """ This is a function that can be run from a modal operator
        to select the 3D object the mouse is hovered over.
    """
    exclude = exclude or []
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


#--- Geometry functions -------------------------------------------------------#

def get_fixed_angle(context, first_atom, coord_3d, angle_list=None):
    debug_print("mb_utils.get_fixed_angle", 6)
    angle_list = angle_list or []
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

def get_fixed_length(context, first_atom, second_atom, coord_3d, length=-1):
    debug_print("mb_utils.get_fixed_length", 6)
    if length < 0:
        r1 = context.scene.mb.elements[first_atom.mb.element].covalent
        r2 = context.scene.mb.elements[second_atom.mb.element].covalent
        length = r1 + r2
    bond_vector = (coord_3d - first_atom.location).normalized()
    return first_atom.location + bond_vector * length

def get_fixed_geometry(context, first_atom, new_atom, coord_3d, geometry):
    debug_print("mb_utils.get_fixed_geometry", 6)
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
    
def check_ob_dimensions(ob):
    debug_print("mb_utils.check_ob_dimensions", 6)
    if ob.dimensions.x < 0.0001: #< ob.mb.get_molecule().bond_radius:
        toggle = {'PLANE_X': 'PLANE_Z', 'PLANE_Z': 'PLANE_X'}
        c = ob.constraints["stretch"]
        if c:
            c.keep_axis = toggle[c.keep_axis]


#--- Add object functions -----------------------------------------------------#

def add_atom(context, location, element, atom_name, molecule):
    debug_print("mb_utils.add_atom", 6)
    # get new unique name for object
    name = "atom_{}.{}".format(molecule.index, molecule.atom_index)
    
    mesh_data = get_atom_data(element, molecule, type='MESH', name="")
    new_atom = bpy.data.objects.new(name, mesh_data)
    context.scene.objects.link(new_atom)
    new_atom.location = location
    
    # set mb properties
    new_atom.mb.name = new_atom.name
    new_atom.mb.molecule_name = molecule.name
    new_atom.mb.index = molecule.atom_index
    new_atom.mb.atom_name = atom_name
    
    # add to molecule
    molecule.atom_index += 1
    
    # parent to molecule origin
    new_atom.parent = molecule.objects.parent.get_object()
    
    # updating the element will call update_atom_element, which assigns a mesh, and sets all the drivers
    new_atom.mb.element = element
    # add atom object and mesh to molecule collections
    molecule.add_object(new_atom)
    molecule.add_object(new_atom.data)
    
    return new_atom

def add_bond(context, first_atom, second_atom):
    #bond_type = "curve_modifier"
    bond_type = "constraint"
    
    debug_print("mb_utils.add_bond", 6)
    if first_atom == second_atom:
        debug_print('WARNING: add_bond: first_atom == second_atom', 3)
        return None
    for b in first_atom.mb.bonds:
        if second_atom.name in b.get_object().mb.bonded_atoms:
            debug_print('WARNING: add_bond: Bond {}-{} already exists'.format(first_atom.mb.index, second_atom.mb.index), 3)
            return None
    #if first_atom.mb.molecule_name != second_atom.mb.molecule_name:
        ## TODO join molecule properties if connecting different molecules
        #pass
    # get new unique name for bond
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    name = "bond_{}.{}-{}.{}".format(first_mol.index, first_atom.mb.index, second_mol.index, second_atom.mb.index)
    
    bond_mesh = get_bond_data(first_mol, type='MESH')
    #bond_mesh = get_bond_data(type='CURVE')
    #if True:
        #bond_mesh.name = name
    new_bond = bpy.data.objects.new(name, bond_mesh)
    context.scene.objects.link(new_bond)
    new_bond.hide = (first_mol.draw_style == 'BALLS')
    
    # set mb properties
    new_bond.mb.type = 'BOND'
    new_bond.mb.name = new_bond.name
    new_bond.mb.molecule_name = first_mol.name
    new_bond.mb.add_bonded_atom(first_atom)
    new_bond.mb.add_bonded_atom(second_atom)
    
    # add bond to atoms mb props
    first_atom.mb.add_bond(new_bond)
    second_atom.mb.add_bond(new_bond)
    
    # add it to first molecule collection
    first_mol.add_object(new_bond)
    
    if bond_type == "constraint":
        # don't parent, as parenting also affects the scale
        c = new_bond.constraints.new('COPY_LOCATION')
        c.name = "parent"
        c.target = first_atom
        
        c = new_bond.constraints.new('STRETCH_TO')
        c.name = "stretch"
        c.rest_length = 1.0
        c.volume = 'NO_VOLUME'
        c.target = second_atom
    elif bond_type == "curve_modifier":
        # get new bezier curve that hooks two first and second atom
        curve_ob = first_mol.objects.bond_curve.get_object()
        molcenter = first_mol.objects.parent.get_object().location
        bc_data = curve_ob.data
        # add spline to curve
        sp = bc_data.splines.new('BEZIER')
        sp.bezier_points.add(2)
        for bp, atom in zip(sp.bezier_points[-2:], (first_atom, second_atom)):
            bp.co = atom.location - molcenter
            bp.handle_left_type = "VECTOR"
            bp.handle_right_type = "VECTOR"
        
        
        
        
    assign_bond_material(new_bond)
    set_bond_drivers(context, new_bond, new_bond.mb.get_molecule())
    new_bond.parent = first_mol.objects.parent.get_object()
    
    return new_bond


#--- Get Mesh functions -------------------------------------------------------#

def get_atom_data(element, molecule, type='MESH', name=""):
    debug_print("mb_utils.get_atom_data", 6)
    if type == 'MESH':
        if name:
            mesh_name = name
        else:
            mesh_name = "atom_mesh_{}.{}".format(molecule.index, element)
            debug_print("Create {}.".format(mesh_name), 5)
        me = bpy.context.blend_data.meshes.get(mesh_name)
        if not me:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            refine = molecule.refine_atoms
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
            bpy.data.objects.remove(new_atom)
            
            # restore old selection
            for o in selected:
                o.select = True
            bpy.context.scene.objects.active = last_active
        return me

def get_bond_curve(molecule, name=""):
    debug_print("mb_utils.get_bond_curve", 6)
    if name:
        curve_name = name
    else:
        curve_name = "bond_curve_{}".format(molecule.index)
    curve = bpy.context.blend_data.curves.get(curve_name)
    if not curve:
        curve = bpy.data.curves.new(curve_name, "CURVE")
        curve.show_handles = False
        
def get_bond_data(molecule, type='MESH', name=""):
    debug_print("mb_utils.get_bond_data", 6)
    new_bond = None
    if type == 'MESH':
        if name:
            data_name = name
        else:
            data_name = "bond_mesh_{}".format(molecule.index)
            debug_print("Create {}.".format(data_name), 5)
        data = bpy.context.blend_data.meshes.get(data_name)
        if not data:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            bpy.ops.mesh.primitive_cylinder_add(location=(0,0,0), 
                vertices=molecule.refine_bonds*2, depth=1, radius=1.0,
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
                if len(edge.link_faces) == 2 and edge.calc_length() == 1.0:
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

    # TODO need to fix it so that each molecule has it's own bond.
    #if type == 'CURVE':
        ## save last selection to restore later
        #selected = bpy.context.selected_objects
        #last_active = bpy.context.object
        
        ## get bevel object
        #bevel_name = 'bond_bevel_{}'.format(molecule)
        #bond_bevel = bpy.context.blend_data.objects.get(bevel_name)
        #if not bond_bevel:
            #debug_print("Create bond_bevel.", 3)
            #bpy.ops.curve.primitive_bezier_circle_add(radius=0.15, location=(0,0,0))
            #bond_bevel = bpy.context.object
            #bond_bevel.name = bevel_name
            
        #bpy.ops.curve.primitive_bezier_curve_add(location=(0,0,0))
        #new_bond = bpy.context.object
        
        #data = new_bond.data
        #data.show_handles = False
        
        #bp = data.splines[0].bezier_points
        #for i in range(2):
            #bp[i].handle_left_type = 'VECTOR'
            #bp[i].handle_right_type = 'VECTOR'
        
        #data.bevel_object = bond_bevel
        ##bp[0].co = (0,0,0)
        ##bp[1].co = (0,1,0)
    
    if new_bond:
        # finally delete object and reselect old selection
        bpy.context.scene.objects.unlink(new_bond)
        bpy.data.objects.remove(new_bond)
        for o in selected:
            o.select = True
        bpy.context.scene.objects.active = last_active
    
    molecule.add_object(data)
    return data

def get_arrow_data(type='MESH', name="arrow", material=None, 
        radius = 0.1, ring_y = 0.9, ring_scale = 2):
    data = bpy.context.blend_data.meshes.get(name)
    if not data:
        debug_print("Create {} mesh.".format(name), 4)
        # Make arrow mesh
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
            vert.co.z = tmp_co.y
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
        # remove faces with 6 verts
        for f in bm.faces:
            if len(f.verts) == 6:
                bm.faces.remove(f)
        
        # now sort bm.verts
        # v.co.y is either 0, 0.5, or 1.0. So multiply y with at least 4pi to sort by y value first
        key = lambda v: v.co.y * 15 + math.atan2(v.co.x, v.co.z)
        verts_sorted = sorted((v for v in bm.verts if (0 < v.co.length and v.co.length != 1.0)), key=key)
        #for i, v in enumerate(verts_sorted):
        #    v.index = i
        
        # add new faces
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
        
        # make everything smooth
        for f in bm.faces:
            f.smooth = True
        
        bm.to_mesh(ob.data)
        bm.free()
        data = ob.data
        data.name = name
        bpy.context.scene.objects.unlink(ob)
    return data
#--- Driver setting functions -------------------------------------------------#

def set_atom_drivers(context, atom, molecule):
    debug_print("mb_utils.set_atom_drivers", 6)
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
    debug_print("mb_utils.set_bond_drivers", 6)
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


#--- Material functions -------------------------------------------------------#

def new_material(name, color=(0.8, 0.8, 0.8), molecule=None):
    debug_print("mb_utils.new_material", 6)
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
    debug_print("mb_utils.assign_atom_material", 6)
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
    debug_print("mb_utils.assign_bond_material", 6)
    
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

        
# TODO this could be a very central function, and needs a lot of time to think about
#def add_object_to_molecule(context, ob, molecule):
    debug_print("mb_utils.add_object_to_molecule", 6)
    ## check if object is already part of molecule
    #if ob.mb.molecule_name == molecule.name:
        #return
    #else:
        ## update name, delete from old molecule and add to new one
        #if ob.type == 'ATOM':
            #old_name = ob.name
        
        #assign_atom_material(ob, molecule)
        #set_atom_drivers(context, ob, molecule)


    
def update_radius_type(self, context):
    debug_print("mb_utils.update_radius_type", 6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.get_object(), self)

def set_draw_style(self, context):
    debug_print("mb_utils.set_draw_style", 6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.get_object(), self)
    
    hide = (self.draw_style == 'BALLS')
    for bond in self.objects.bonds:
        bond_ob = bond.get_object()
        bond_ob.hide = hide

#******************************************************************************#
#*** Import Export utils ******************************************************#
#******************************************************************************#



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