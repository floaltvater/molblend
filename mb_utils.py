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

import math
from bisect import bisect_left

import bpy
import bmesh
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

from molblend.elements_default import ELEMENTS as ELEMENTS_DEFAULT
from molblend.mb_helper import debug_print

class enums():
    object_types = [
        ('NONE', "None", "None"),
        ('ATOM', "Atom", "Atom"),
        ('BOND', "Bond", "Bond"),
        ('PARENT', "Parent", "Parent"),
        ]
    radius_types = [
        ('covalent', 'covalent', 'covalent'),
        ('vdw', 'van der Waals', 'van der Waals'),
        ('constant', 'constant', 'constant'),
        ]
    molecule_styles = [
        ('BALLS', 'Balls', 'Space filling'),
        ('BAS', 'Balls and Sticks', 'Balls and Sticks'),
        ('STICKS', 'Sticks', 'Sticks'),
        ]
    bond_material = [
        ('ATOMS', "Atoms", "Same as atoms"),
        ('GENERIC', "Generic" , "Single bond color"),
        ]
    bond_types = [
        ('CONSTRAINT', "constraint", "constrained by the two bonded atoms"),
        ('STATIC', "static", "independent bonds, don't move with atoms"),
        ]
    geometries = [
        ('NONE', "None", "No geometry constraints"),
        ('GENERAL', "General", 
         "Angles are multiples of 30 and 45 deg. in the view plane"),
        ('LINEAR', "Linear", "Linear or sp"),
        ('TRIGONAL', "Trig. planar", "Trigonal planar or sp2"),
        ]
    angstrom_per_unit = [
        ('1.0', "Angstrom", "Angstrom"),
        ('0.529177249', "Bohr", "Bohr"),
        ('0.01', "pm", "Picometer"),
        ('OTHER', "Other", "Custom Unit"),
        ]
    file_types = [
        ('XYZ', "xyz", "xyz format"),
        ('PDB', "pdb", "Protein Databank format"),
        ('B4W', "b4w", "standalone HTML with embedded 3D viewer")
        ]


#--- Update functions --------------------------------------------------------#

def update_all_meshes(self, context):
    debug_print("mb_utils.update_all_meshes", level=6)
    # TODO this callback might be too heavy for scenes with lots of meshes
    for me in bpy.data.meshes:
        me.update()


def update_active_mode(self, context):
    debug_print("mb_utils.update_active_mode", level=6)
    if self.max_mode == 0:
        if self.active_mode == 0:
            return
        else:
            self.active_mode = 0
            return
    elif self.active_mode > self.max_mode:
        self.active_mode = self.max_mode
        return
    
    for atom in self.objects.atoms:
        atom_ob = atom.get_object()
        atom_id = get_atom_id(self.index, atom_ob.mb.index)
        anim_data = atom_ob.animation_data
        if anim_data:
            action = anim_data.action
            if action:
                avec = atom_ob.mb.modes[self.active_mode].vec * self.mode_scale
                for dim in range(3):
                    fcu = action.fcurves[dim]
                    kf1 = fcu.keyframe_points[0]
                    kf2 = fcu.keyframe_points[1]
                    middle = (kf1.co[1] + kf2.co[1]) / 2.0
                    
                    for p in range(3):
                        loc = middle + pow(-1, p) * avec[dim]
                        fcu.keyframe_points[p].co = 1.0 + 10*p, loc
                        if self.active_mode == 0:
                            fcu.keyframe_points[p].interpolation = 'LINEAR'
                        else:
                            fcu.keyframe_points[p].interpolation = 'BEZIER'
                    fcu.update()
            else:
                debug_print("No action for atom '{}'.".format(atom_id),
                            level=1)
        else:
            debug_print("No animation data for atom '{}'.".format(atom_id),
                        level=1)
    
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


def update_atom_element(self, context):
    '''
    assign a mesh, give a new object name etc.
    '''
    debug_print("mb_utils.update_atom_element", level=6)
    # remove all spaces
    if self.element.strip() != self.element:
        self.element = self.element.strip()
        return
    # Comparison is case insensitive
    # get dictionary that maps all lower case elements to exact upper/lower 
    # case as it appears in the list
    elements_map = dict([(element.name.lower(), element.name)
                         for element in context.scene.mb.elements])
    # if element is not yet in scene elements list, add new element. 
    # Use default settings
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
    
    # assign type last, to be able to check if element is newly assigned or
    # just updated
    self.type = 'ATOM'


def update_bond_material(self, context):
    debug_print("mb_utils.update_bond_material", level=6)
    for bond in self.objects.bonds:
        assign_bond_material(bond.get_object())


def update_refine_atoms(self, context):
    debug_print("mb_utils.update_refine_atoms", level=6)
    if self.refine_atoms < 3:
        self.refine_atoms = 3
    replaced_elements = set()
    for mesh in self.objects.meshes:
        if "atom_mesh" in mesh.name:
            element = mesh.name.split(".")[-1]
            if not element in replaced_elements:
                data = mesh.get_data()
                # get new temporary atom mesh with new refine value
                new_data = get_atom_data(element, self, type='MESH', 
                                         name="tmp_mesh")
                # replace mesh data
                bm = bmesh.new()
                bm.from_mesh(new_data)
                bm.to_mesh(data)
                bm.free()
                replaced_elements.add(element)
                # delete temporary mesh
                bpy.data.meshes.remove(new_data)


def update_refine_bonds(self, context):
    debug_print("mb_utils.update_refine_bonds", level=6)
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
    """
    Change file extension when filetype is changed. Replace known extensions.
    Append to everything else.
    """
    debug_print("mb_utils.update_export_file_type", level=6)
    if self.filepath:
        filetypes = {'XYZ': ".xyz",
                     'PDB': ".pdb"}
        ext = filetypes[self.file_type]
        other_ext = [f for f in filetypes.values() if f != ext]
        
        if self.filepath[-4:] in other_ext:
            self.filepath = self.filepath[:-4] + ext
        else:
            self.filepath = self.filepath + ext


def update_molecule_selection(self, context):
    debug_print("mb_utils.update_molecule_selection", level=6)
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

def update_show_bond_lengths(self, context):
    if self.show_bond_lengths:
        bpy.ops.mb.show_bond_lengths()
        
def update_show_bond_angles(self, context):
    if self.show_bond_lengths:
        bpy.ops.mb.show_bond_angles()

#--- General functions -------------------------------------------------------#

def create_new_keyconfig(name, context):
    debug_print("mb_utils.create_new_keyconfig", level=6)
    
    if not name in context.window_manager.keyconfigs:
        debug_print('Creating new keyconfiguration "{}".'.format(name),
                    level=3)
        from molblend import mb_keyconfig
        #keyconfig = mb_keyconfig.kc
    else:
        debug_print('Keyconfiguration "{}" already exists.'.format(name), 
                    level=3)
    
    return context.window_manager.keyconfigs.get("MolBlend")


def add_element(context, element, element_dict):
    '''
    add element data to scene
    '''
    debug_print("mb_utils.add_element", level=6)
    default = ELEMENTS_DEFAULT["Default"]
    
    new = context.scene.mb.elements.add()
    new.name = element
    new.element = element
    new.element_name = element_dict.get("element name", 
                                        default["element name"])
    new.atomic_number = element_dict.get("atomic number", 
                                         default["atomic number"])
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
    debug_print("mb_utils.initialize_elements", level=6)
    for element, data in ELEMENTS_DEFAULT.items():
        add_element(context, element, data)


def is_initialized(context):
    try:
        context.scene.mb.elements['Default']
        context.scene.mb.globals.atom_scales['BALLS']
        context.scene.mb.globals.atom_scales['BAS']
        context.scene.mb.globals.atom_scales['STICKS']
    except KeyError as e:
        return False
    
    return context.user_preferences.system.use_scripts_auto_execute


#--- Viewport functions ------------------------------------------------------#

def get_region_data(context, x, y):
    
    for area in context.screen.areas:
        if area.type != 'VIEW_3D':
            continue
        is_quadview = len(area.spaces.active.region_quadviews) != 0
        i = -1
        for region in area.regions:
            if region.type == 'WINDOW':
                i += 1
                if (x > region.x and
                    y > region.y and
                    x < region.width + region.x and
                    y < region.height + region.y):
                    
                    if is_quadview:
                        rv3d = area.spaces.active.region_quadviews[i]
                    else:
                        rv3d = area.spaces.active.region_3d
                    return (region, rv3d)
    return (None, None)

def mouse_2d_to_location_3d(context, mouse2d, depth=Vector((0, 0, 0)),
                            region=None, rv3d=None):
    debug_print("mb_utils.mouse_2d_to_location_3d", level=6)
    x, y = mouse2d
    
    # Get region and region data from mouse position.
    # If region is given, passed rv3d is ignored.
    if region == None:
        region, rv3d = get_region_data(context, x, y)
    
    # mouse coordinates relative to region
    coord2d = (x - region.x, y - region.y)
    if depth:
        depth_location = depth
    else:
        depth_location = context.scene.cursor_location.copy()
    return view3d_utils.region_2d_to_location_3d(region, rv3d, coord2d, 
                                                 depth_location)


def return_cursor_object(context, event, ray_max=10000.0, exclude=None, 
                         mb_type=''):
    """ This is a function that can be run from a modal operator
        to select the 3D object the mouse is hovered over.
    """
    debug_print("mb_utils.return_cursor_object", level=6)
    exclude = exclude or []
    # get the context arguments
    scene = context.scene
    x, y = event.mouse_x, event.mouse_y
    region, rv3d = get_region_data(context, x, y)
    if not rv3d:
        return None
    coord2d = (x - region.x, y - region.y)
    # get the ray from the viewport and mouse
    view_vector = view3d_utils.region_2d_to_vector_3d(region, rv3d, coord2d)
    
    # have ray origin a little in front of actual origin, otherwise it might 
    # not get all of the objects
    ray_origin = (view3d_utils.region_2d_to_origin_3d(region, rv3d, coord2d) +
                  0.5 * view_vector * ray_max)
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
            #print(obj.ray_cast(ray_origin_obj, ray_target_obj))
            result, hit, normal, face_index = obj.ray_cast(ray_origin_obj, 
                                                           ray_target_obj)
            if face_index != -1:
                return hit, normal, face_index
            else:
                return None, None, None
        except ValueError as e:
            debug_print("ERROR in obj_ray_cast: {}: {}".format(obj.name, e), 
                        level=0)
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




def check_ob_dimensions(ob):
    debug_print("mb_utils.check_ob_dimensions", level=6)
    if ob.dimensions.x < 0.0001: #< ob.mb.get_molecule().bond_radius:
        toggle = {'PLANE_X': 'PLANE_Z', 'PLANE_Z': 'PLANE_X'}
        c = ob.constraints.get("mb.stretch", None)
        if c:
            c.keep_axis = toggle[c.keep_axis]


#--- Add object functions ----------------------------------------------------#
def get_atom_id(mol_index, atom_index):
    return "{:>02d}.{:>04d}".format(mol_index, atom_index)

def add_atom(context, location, element, atom_name, molecule):
    debug_print("mb_utils.add_atom", level=6)
    # get new unique name for object
    name = "atom_{}".format(get_atom_id(molecule.index, molecule.atom_index))
    
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
    
    # updating the element will call update_atom_element, which assigns a mesh,
    # and sets all the drivers
    new_atom.mb.element = element
    # add atom object and mesh to molecule collections
    molecule.add_object(new_atom)
    molecule.add_object(new_atom.data)
    
    return new_atom


def add_bond(context, first_atom, second_atom, bond_type="CONSTRAINT"):
    
    debug_print("mb_utils.add_bond", level=6)
    if first_atom == second_atom:
        debug_print('WARNING: add_bond: first_atom == second_atom', level=3)
        return None
    for b in first_atom.mb.bonds:
        bob = b.get_object()
        if bob != None and second_atom.name in bob.mb.bonded_atoms:
            debug_print(
                "WARNING: add_bond: Bond {}-{} already exists".format(
                    first_atom.mb.index, second_atom.mb.index), level=3)
            return None
    # get new unique name for bond
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    name = "bond_{}-{}".format(
        get_atom_id(first_mol.index, first_atom.mb.index), 
        get_atom_id(second_mol.index, second_atom.mb.index)
        )
    bond_mesh = get_bond_data(first_mol, type='MESH')
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
    
    if bond_type == "CONSTRAINT":
        # don't parent, as parenting also affects the scale
        c = new_bond.constraints.new('COPY_LOCATION')
        c.name = "mb.parent"
        c.target = first_atom
        
        c = new_bond.constraints.new('STRETCH_TO')
        c.name = "mb.stretch"
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
    elif bond_type == "STATIC":
        y_axis = Vector((0,1,0))
        loc1 = first_atom.location
        loc2 = second_atom.location
        # Location
        location = loc1
        vec = (loc2 - loc1)
        angle = vec.angle(y_axis, 0)
        # vector of rotation
        axis = y_axis.cross(vec)
        new_bond.rotation_euler = Matrix.Rotation(angle, 4, axis).to_euler()
        new_bond.location = loc1
        new_bond.scale[1] = vec.length
        
    assign_bond_material(new_bond)
    set_bond_drivers(context, new_bond, new_bond.mb.get_molecule())
    new_bond.parent = first_mol.objects.parent.get_object()
    
    return new_bond


#--- Get Mesh functions ------------------------------------------------------#

def get_atom_data(element, molecule, type='MESH', name=""):
    debug_print("mb_utils.get_atom_data", level=6)
    if type == 'MESH':
        if name:
            mesh_name = name
        else:
            mesh_name = "atom_mesh_{}.{}".format(molecule.index, element)
            debug_print("Create {}.".format(mesh_name), level=5)
        me = bpy.context.blend_data.meshes.get(mesh_name)
        if not me:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            refine = molecule.refine_atoms
            # create uv sphere and get mesh data
            bpy.ops.mesh.primitive_uv_sphere_add(
                location=(0,0,0), segments=refine*2, ring_count=refine)
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
    debug_print("mb_utils.get_bond_curve", level=6)
    if name:
        curve_name = name
    else:
        curve_name = "bond_curve_{}".format(molecule.index)
    curve = bpy.context.blend_data.curves.get(curve_name)
    if not curve:
        curve = bpy.data.curves.new(curve_name, "CURVE")
        curve.show_handles = False


def get_bond_data(molecule, type='MESH', name=""):
    debug_print("mb_utils.get_bond_data", level=6)
    new_bond = None
    if type == 'MESH':
        if name:
            data_name = name
        else:
            data_name = "bond_mesh_{}".format(molecule.index)
            debug_print("Create {}.".format(data_name), level=5)
        data = bpy.context.blend_data.meshes.get(data_name)
        if not data:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            bpy.ops.mesh.primitive_cylinder_add(
                location=(0,0,0), vertices=molecule.refine_bonds*2, 
                depth=1, radius=1.0)#, end_fill_type="NOTHING")
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
                    if hasattr(edge, "ensure_lookup_table"):
                        edge.ensure_lookup_table()
                    e, v = bmesh.utils.edge_split(edge, edge.verts[0], 0.5)
                    new_verts.append(v)
            n_verts = len(new_verts)
            
            # bad hack, but don't understand how bmesh.utils.face_split works
            # remove all faces
            for f in bm.faces:
                bm.faces.remove(f)
            # now sort bm.verts
            # v.co.y is either 0, 0.5, or 1.0. 
            # So multiply y with at least 4pi to sort by y value first
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
            key = lambda f: (f.calc_center_median().y * 15 +
                             math.atan2(f.normal.x, f.normal.z))
            half_faces = len(bm.faces)/2
            for i, f in enumerate(sorted((f for f in bm.faces), key=key)):
                f.index = i
            
            bm.to_mesh(data)
            bm.free()
            data.update()
    
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
        debug_print("Create {} mesh.".format(name), level=4)
        # Make arrow mesh
        bpy.ops.mesh.primitive_cylinder_add(
            location=(0,0,0), radius=radius, vertices=8, depth=1,
            end_fill_type="TRIFAN")
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
                if hasattr(edge, "ensure_lookup_table"):
                    edge.ensure_lookup_table()
                e, v = bmesh.utils.edge_split(edge, edge.verts[0], 0.5)
                new_verts.append(v)
        n_verts = len(new_verts)
        
        # bad hack, but don't understand how bmesh.utils.face_split works
        # remove faces with 6 verts
        for f in bm.faces:
            if len(f.verts) == 6:
                bm.faces.remove(f)
        
        # now sort bm.verts
        # v.co.y is either 0, 0.5, or 1.0.
        # So multiply y with at least 4pi to sort by y value first
        key = lambda v: v.co.y * 15 + math.atan2(v.co.x, v.co.z)
        verts_sorted = sorted(
            (v for v in bm.verts if (0 < v.co.length and v.co.length != 1.0)), 
            key=key
            )
        
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


#--- Driver setting functions ------------------------------------------------#

def set_atom_drivers(context, atom, molecule):
    debug_print("mb_utils.set_atom_drivers", level=6)
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
        targ.data_path = 'mb.elements["{}"].{}'.format(atom.mb.element, 
                                                       molecule.radius_type)
        
        var = drv.variables.get('atom_scale')
        if not var:
            var = drv.variables.new()
            var.name = 'atom_scale' # name to use in scripting
            var.type = 'SINGLE_PROP'
        targ = var.targets[0]
        targ.id_type = 'SCENE'
        targ.id = context.scene
        targ.data_path = 'mb.molecules["{}"].atom_scales["{}"].val'.format(
                         molecule.name, molecule.draw_style)
        
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
    debug_print("mb_utils.set_bond_drivers", level=6)
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
            targ.data_path = 'mb.molecules["{}"].bond_radius'.format(
                             molecule.name)

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
                targ.id = bond.mb.bonded_atoms[i].get_object()
                targ.transform_type = 'LOC_{}'.format(dim)
                targ.transform_space = 'TRANSFORM_SPACE'


#--- Material functions ------------------------------------------------------#

def new_material(name, color=(0.8, 0.8, 0.8), molecule=None):
    '''
    creates new material. If molecule is given, the molecule.mb.bond_color 
    property will be added as driver.
    '''
    debug_print("mb_utils.new_material", level=6)
    material = bpy.data.materials.new(name)
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
            targ.data_path = 'mb.molecules["{}"].bond_color[{}]'.format(
                             molecule.name, i)
    
    if bpy.context.scene.render.engine == 'CYCLES':
        material.use_nodes = True
        # add driver to rendered color to be the same as display color
        nodesocketcolor = material.node_tree.nodes['Diffuse BSDF'].inputs[0]
        nodesocketcolor.default_value[:3] = color
        for i in range(3): # not for alpha channel
            fcurve =  material.driver_add('diffuse_color', i)
            drv = fcurve.driver
            drv.type = 'AVERAGE'
            drv.show_debug_info = True
            
            var = drv.variables.new()
            var.name = 'diffuse_color' # name to use in scripting
            var.type = 'SINGLE_PROP'
            targ = var.targets[0]
            targ.id_type = 'MATERIAL'
            targ.id = material
            targ.data_path = ('node_tree.nodes["Diffuse BSDF"].inputs[0]'
                              + '.default_value[{}]'.format(i))
    return material


def assign_atom_material(ob, molecule):
    debug_print("mb_utils.assign_atom_material", level=6)
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
    debug_print("mb_utils.assign_bond_material", level=6)
    
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
            debug_print("Atom colored bonds are not yet supported with "
                        "Curve bonds.", level=2)
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


def update_radius_type(self, context):
    debug_print("mb_utils.update_radius_type", level=6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.get_object(), self)


def set_draw_style(self, context):
    debug_print("mb_utils.set_draw_style", level=6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.get_object(), self)
    
    hide = (self.draw_style == 'BALLS')
    for bond in self.objects.bonds:
        bond_ob = bond.get_object()
        bond_ob.hide = hide
