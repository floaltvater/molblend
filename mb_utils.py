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
    mesh_types = [
        ('NONE', "None", "None"),
        ('ELEMENT', "Element", "Element"),
        ('BOND', "Bond", "Bond"),
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
    mode_file_format = [
        ('XYZ', "xyz", "xyz-style format"),
        ('QE_DYNMAT', "QE dynmat", "Quantum ESPRESSO output"),
        ]


#--- Update functions --------------------------------------------------------#

def update_all_meshes(self, context):
    debug_print("mb_utils.update_all_meshes", level=6)
    # TODO this callback might be too heavy for scenes with lots of meshes
    for me in bpy.data.meshes:
        me.update()

def update_active_mode(self, context):
    debug_print("mb_utils.update_active_mode", level=6)
    if len(self.qpts) == 0:
        if context.screen.is_animation_playing:
            bpy.ops.screen.animation_play()
            context.scene.frame_current = 1
        return
    
    modes = self.qpts[self.active_nqpt].modes
    if len(modes) == 0:
        if self.active_mode == 0:
            return
        else:
            self.active_mode = 0
            return
    elif self.active_mode >= len(modes):
        self.active_mode = len(modes) - 1
        return
    
    for atom in self.objects.atoms:
        #update_mode_drivers(atom.object, self)
        update_mode_action(atom.object, self)
    
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
    atom_ob = self.object
    molecule = self.get_molecule()
    
    # update mesh and material
    me = get_atom_data(self.element, molecule)
    atom_ob.data = me
    assign_atom_material(atom_ob, molecule)
    
    set_atom_drivers(context, atom_ob, molecule)
    
    # update bond materials
    for bond in self.bonds:
        assign_bond_material(bond.object)
    
    # assign type last, to be able to check if element is newly assigned or
    # just updated
    self.type = 'ATOM'


def update_bond_material(self, context):
    debug_print("mb_utils.update_bond_material", level=6)
    for bond in self.objects.bonds:
        assign_bond_material(bond.object)


def update_refine_atoms(self, context):
    debug_print("mb_utils.update_refine_atoms", level=6)
    if self.refine_atoms < 2:
        self.refine_atoms = 2
        return
    replaced_elements = set()
    for mesh in self.meshes:
        print(mesh.name, mesh.data.mb.type)
        if mesh.data.mb.type == "ELEMENT":
            element = mesh.name
            #if not element in replaced_elements:
            data = mesh.data
            # get new temporary atom mesh with new refine value
            new_data = get_atom_data(element, self, type='MESH', 
                                        mesh_name="tmp_mesh")
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
    if self.refine_bonds < 2:
        self.refine_bonds = 2
        return
    mesh = self.meshes.get("bond")
    if mesh:
        data = mesh.data
        # get new temporary bond mesh with new refine value
        name = "tmp_mesh"
        new_data = get_bond_data(self, type='MESH', mesh_name="tmp_mesh")
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
            for item in col:
                item.object.select = True
        parent = mol.objects.parent.object
        parent.select = True
        context.scene.objects.active = parent

def update_molecule_name(self, context):
    self.objects.parent.object.name = self.name_mol

def update_show_bond_lengths(self, context):
    if self.show_bond_lengths:
        bpy.ops.mb.show_bond_lengths()
        
def update_show_bond_angles(self, context):
    if self.show_bond_angles:
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

def callback_draw_length(self, context):
    try:
        font_id = 0
        blf.size(font_id, 12, 72)
        offset = 0
        
        rv3d = context.space_data.region_3d
        width = context.region.width
        height = context.region.height
        persp_mat = rv3d.perspective_matrix
        persinv = persp_mat.inverted()
        
        for ob in context.selected_objects:
            if ob.mb.type == "BOND":
                locs = [o.object.mb.world_location for o in ob.mb.bonded_atoms]
                co_3d = (locs[0] + locs[1]) / 2.
                prj = persp_mat * co_3d.to_4d()
                x = width/2 + width/2 * (prj.x / prj.w)
                y = height/2 + height/2 * (prj.y / prj.w)
                blf.position(font_id, x, y, 0)
                blf.draw(font_id, "{:6.4f}".format((locs[1]-locs[0]).length))
    except:
        debug_print(sys.exc_info()[0], level=1)
        context.scene.mb.globals.show_bond_lengths = False

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
    
    mesh_data = get_atom_data(element, molecule)
    new_atom = bpy.data.objects.new(name, mesh_data)
    context.scene.objects.link(new_atom)
    new_atom.location = location
    
    # set mb properties
    new_atom.mb.name = new_atom.name
    new_atom.mb.molecule_ident = molecule.name
    new_atom.mb.atom_name = atom_name
    
    # parent to molecule origin
    new_atom.parent = molecule.objects.parent.object
    
    # updating the element will call update_atom_element, which assigns a mesh,
    # and sets all the drivers
    new_atom.mb.element = element
    # add atom object and mesh to molecule collections
    molecule.add_object(new_atom)
    
    return new_atom


def add_bond(context, first_atom, second_atom, bond_type="CONSTRAINT"):
    
    debug_print("mb_utils.add_bond", level=6)
    if first_atom == second_atom:
        debug_print('WARNING: add_bond: first_atom == second_atom', level=3)
        return None
    for b in first_atom.mb.bonds:
        bob = b.object
        if bob != None:
            for ba in bob.mb.bonded_atoms:
                if ba.object == second_atom:
                    debug_print(
                        "WARNING: add_bond: Bond {}-{} already exists".format(
                            first_atom.mb.index, second_atom.mb.index),
                        level=3)
                    return None
    # get new unique name for bond
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    name = "bond_{}-{}".format(
        get_atom_id(first_mol.index, first_atom.mb.index), 
        get_atom_id(second_mol.index, second_atom.mb.index)
        )
    bond_mesh = get_bond_data(first_mol)
    new_bond = bpy.data.objects.new(name, bond_mesh)
    context.scene.objects.link(new_bond)
    new_bond.hide = (first_mol.draw_style == 'BALLS')
    
    # set mb properties
    new_bond.mb.type = 'BOND'
    new_bond.mb.name = new_bond.name
    new_bond.mb.molecule_ident = first_mol.name
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
    new_bond.parent = first_mol.objects.parent.object
    
    return new_bond


#--- Get Mesh functions ------------------------------------------------------#

def get_atom_data(element, molecule, type='MESH', mesh_name=""):
    """
    Retrieve mesh for a certain element. If mesh_name is given, the mesh is 
    retrieved from bpy.data.meshes if it exists, or created and assigned that
    name. Otherwise the mesh is retrieved from molecule.meshes.
    """
    debug_print("mb_utils.get_atom_data", level=6)
    if type == 'MESH':
        if mesh_name:
            me = bpy.context.blend_data.meshes.get(mesh_name)
        else:
            element = element.capitalize()
            atom_name = "atom_mesh_{}.{}".format(molecule.index, element)
            item = molecule.meshes.get(element)
            if item:
                me = item.data
            else:
                me = None
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
            me.name = mesh_name or atom_name
            me.mb.type = "ELEMENT"
            if not mesh_name:
                print("add")
                item = molecule.meshes.add()
                item.name = element
                item.data = me
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

def get_bond_data(molecule, type='MESH', mesh_name=""):
    debug_print("mb_utils.get_bond_data", level=6)
    new_bond = None
    if type == 'MESH':
        if mesh_name:
            me = bpy.context.blend_data.meshes.get(mesh_name)
        else:
            bond_name = "bond"
            item = molecule.meshes.get(bond_name)
            if item:
                me = item.data
            else:
                me = None
        if not me:
            # save last selection to restore later
            selected = bpy.context.selected_objects
            last_active = bpy.context.object
            
            bpy.ops.mesh.primitive_cylinder_add(
                location=(0,0,0), vertices=molecule.refine_bonds*2, 
                depth=1, radius=1.0)#, end_fill_type="NOTHING")
            new_bond = bpy.context.object
            for i in range(2):
                new_bond.data.materials.append(None)
            
            me = new_bond.data
            me.name = mesh_name or "bond_mesh_{}".format(molecule.index)
            me.mb.type = "BOND"
            bm = bmesh.new()
            bm.from_mesh(me)
            
            # rotate and shrink first, then add another row of vertices
            for vert in bm.verts:
                # rotate 90 degrees around x, and shift along y axis
                tmp_co = vert.co.copy()
                vert.co.y = -tmp_co.z + .5
                vert.co.z = -tmp_co.y
                if vert.co.y > 0.01:
                    vert.select = False
            new_verts = []
            bm.edges.ensure_lookup_table()
            for edge in bm.edges:
                if abs(edge.verts[0].co.y - edge.verts[1].co.y) > .5:
                    #if hasattr(edge, "ensure_lookup_table"):
                        #edge.ensure_lookup_table()
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
            
            bm.to_mesh(me)
            bm.free()
            me.update()
    
            # finally delete object and reselect old selection
            bpy.context.scene.objects.unlink(new_bond)
            bpy.data.objects.remove(new_bond)
            for o in selected:
                o.select = True
            bpy.context.scene.objects.active = last_active
            if not mesh_name:
                item = molecule.meshes.add()
                item.name = bond_name
                item.data = me
    return me


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


def set_bond_drivers(context, bond, molecule):
    debug_print("mb_utils.set_bond_drivers", level=6)

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

def clear_modes(molecule):
    #for atom in molecule.objects.atoms:
        #clear_mode_drivers(atom.object)
    while len(molecule.qpts):
        molecule.qpts.remove(0)

#def clear_mode_drivers(atom):
    #atom.driver_remove('location', -1)

#def update_mode_drivers(atom, molecule):
    #debug_print("mb_utils.set_mode_drivers", level=6)
    #fc_list = atom.driver_add('location', -1) # add new driver
    
    #nqpt = molecule.active_nqpt
    #mode = molecule.active_mode
    #nvec = len(molecule.qpts[nqpt].modes[mode].evecs)
    
    #for dim, fcurve in enumerate(fc_list):
        #drv = fcurve.driver
        #drv.type = 'SCRIPTED'
        #drv.show_debug_info = True
        
        #var = drv.variables.get('Real')
        #trg = 'mb.molecules["{}"].qpts[{}].modes[{}].evecs[{}].real[{}]'.format(molecule.name, nqpt, mode, atom.mb.index%nvec, dim)
        #var.targets[0].data_path = trg
        
        #var = drv.variables.get('Imag')
        #trg = 'mb.molecules["{}"].qpts[{}].modes[{}].evecs[{}].imag[{}]'.format(molecule.name, nqpt, mode, atom.mb.index%nvec, dim)
        #var.targets[0].data_path = trg

        #var = drv.variables.get('qvec')
        #trg = 'mb.molecules["{}"].qpts[{}].qvec'.format(molecule.name, nqpt)
        #var.targets[0].data_path = trg
        
        #loc = atom.location[dim]
        
        #qR = "+".join(["qvec[{0}]*sc[{0}]".format(j) for j in range(3)])
        #T = 20
        #arg = "2*pi*({}-(frame-1)/{})".format(qR, T)
        #expr = "{}+scale*(Real*cos({})+Imag*sin({}))".format(loc, arg, arg)
        
        #drv.expression = expr
        

#def set_mode_drivers(context, atom, molecule):
    #debug_print("mb_utils.set_mode_drivers", level=6)
    #fc_list = atom.driver_add('location', -1) # add new driver
    
    #nqpt = molecule.active_nqpt
    #mode = molecule.active_mode
    #nvec = len(molecule.qpts[nqpt].modes[mode].evecs)
    
    #for dim, fcurve in enumerate(fc_list):
        #drv = fcurve.driver
        #drv.type = 'SCRIPTED'
        #drv.show_debug_info = True
        
        #var = drv.variables.get('scale')
        #if not var:
            #var = drv.variables.new()
            #var.name = 'scale' # name to use in scripting
            #var.type = 'SINGLE_PROP'
        #targ = var.targets[0]
        #targ.id_type = 'SCENE'
        #targ.id = context.scene
        #targ.data_path = 'mb.molecules["{}"].mode_scale'.format(molecule.name)
        
        #var = drv.variables.get('sc')
        #if not var:
            #var = drv.variables.new()
            #var.name = 'sc' # name to use in scripting
            #var.type = 'SINGLE_PROP'
        #targ = var.targets[0]
        #targ.id_type = 'OBJECT'
        #targ.id = atom
        #trg = 'mb.supercell'
        #targ.data_path = trg
        
        #var = drv.variables.get('Real')
        #if not var:
            #var = drv.variables.new()
            #var.name = 'Real' # name to use in scripting
            #var.type = 'SINGLE_PROP'
        #targ = var.targets[0]
        #targ.id_type = 'SCENE'
        #targ.id = context.scene
        #trg = 'mb.molecules["{}"].qpts[{}].modes[{}].evecs[{}].real[{}]'.format(molecule.name, nqpt, mode, atom.mb.index%nvec, dim)
        #targ.data_path = trg
        
        #var = drv.variables.get('Imag')
        #if not var:
            #var = drv.variables.new()
            #var.name = 'Imag' # name to use in scripting
            #var.type = 'SINGLE_PROP'
        #targ = var.targets[0]
        #targ.id_type = 'SCENE'
        #targ.id = context.scene
        #trg = 'mb.molecules["{}"].qpts[{}].modes[{}].evecs[{}].imag[{}]'.format(molecule.name, nqpt, mode, atom.mb.index%nvec, dim)
        #targ.data_path = trg
        
        #var = drv.variables.get('qvec')
        #if not var:
            #var = drv.variables.new()
            #var.name = 'qvec' # name to use in scripting
            #var.type = 'SINGLE_PROP'
        #targ = var.targets[0]
        #targ.id_type = 'SCENE'
        #targ.id = context.scene
        #trg = 'mb.molecules["{}"].qpts[{}].qvec'.format(molecule.name, nqpt)
        #targ.data_path = trg
        
        #loc = atom.location[dim]
        
        #qR = "+".join(["qvec[{0}]*sc[{0}]".format(j) for j in range(3)])
        #T = 20
        #arg = "2*pi*({}-(frame-1)/{})".format(qR, T)
        #expr = "{}+scale*(Real*cos({})+Imag*sin({}))".format(loc, arg, arg)
        
        #drv.expression = expr
        ##amp = "scale*qpts[nqpt].modes[mode].evecs[i].real[{}]".format(dim)
        ##drv.expression = "{} +{}*sin((frame-1)*pi/10)".format(loc, amp)

def update_mode_action(atom, molecule):
    action = atom.animation_data.action
    qpt = molecule.qpts[molecule.active_nqpt]
    qvec = qpt.qvec
    sc = atom.mb.supercell
    qR = qvec[0]*sc[0] + qvec[1]*sc[1] + qvec[2]*sc[2]
    T = 20
    nevecs = len(qpt.modes[molecule.active_mode].evecs)
    evec = qpt.modes[molecule.active_mode].evecs[atom.mb.index%nevecs]
    Re = evec.real
    Im = evec.imag
    
    for dim in range(3):
        t_max = T*(qR - math.atan2(Im[dim], Re[dim])/(2*math.pi))
        t0 = (t_max - T/4)
        t0 = t0 - T*(t0//T) - T
        arg = 2*math.pi*(qR-t_max/T)
        cos_max = math.cos(arg)
        sin_max = math.sin(arg)
        fcu = action.fcurves[dim]
        loc = fcu.keyframe_points[0].co[1]
        vec = Re[dim]*cos_max - Im[dim]*sin_max
        
        for p in range(9):
            frame = t0 + 5*p
            coords = loc + pow(-1, p//2)*vec if p%2 else loc
            #if atom.mb.index == 1 and dim == 0:
                #print("{:6.3f}".format(coords),end=" ")
            fcu.keyframe_points[p].co = (frame, coords)
            if p%2 == 0:
                fcu.keyframe_points[p].handle_left = (frame, loc)
                fcu.keyframe_points[p].handle_right = (frame, loc)
        fcu.update()
    #if atom.mb.index == 1:
        #print()
            
def create_mode_action(context, atom, molecule):
    anim_data = atom.animation_data_create()
    atom_id = '{}.{}'.format(molecule.index, atom.mb.index)
    action = bpy.data.actions.new(name="mode_{}".format(atom_id))
    anim_data.action = action
    # make new group
    ag = action.groups.new("Location")
    for dim in range(3):
        fcu = action.fcurves.new(data_path="location", index=dim)
        fcu.group = ag
        fcu.keyframe_points.add(9)
        loc = atom.location[dim]
        for p in range(9):
            fcu.keyframe_points[p].co = 1.0 + 5*p, loc
            fcu.keyframe_points[p].interpolation = 'BEZIER'
        fcu.update()
        
#--- Unit cell functions -----------------------------------------------------#

def draw_unit_cell(molecule, draw_style='ARROWS'):
    # TODO implement different drawing styles
    all_obs = []
    
    unit_vectors = Matrix(molecule.unit_cell_frames[0])

    me = bpy.data.meshes.new("unit_cell_{}".format(molecule.index))
    coords = (
        (0,0,0), #0, O
        unit_vectors[0], #1, x
        unit_vectors[0] + unit_vectors[1], #2, xy
        unit_vectors[1], #3, y
        unit_vectors[2], #4, z
        unit_vectors[0] + unit_vectors[2], #5, xz
        unit_vectors[1] + unit_vectors[2], #6, yz
        unit_vectors[0] + unit_vectors[1] + unit_vectors[2], #7, xyz
        )
    faces = (
        (0, 3, 2, 1),
        (0, 4, 6, 3),
        (0, 1, 5, 4),
        (7, 6, 4, 5),
        (7, 5, 1, 2),
        (7, 2, 3, 6),
        )
    me.from_pydata(coords, [], faces)
    uc_cube = bpy.data.objects.new("unit_cell_{}".format(molecule.index), me)
    bpy.context.scene.objects.link(uc_cube)
    uc_cube.draw_type = 'WIRE'
    
    vg = []
    vg.append(uc_cube.vertex_groups.new('a'))
    vg[-1].add([1], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('b'))
    vg[-1].add([3], 1, 'REPLACE')
    vg.append(uc_cube.vertex_groups.new('c'))
    vg[-1].add([4], 1, 'REPLACE')
    
    all_obs.append(uc_cube)
    
    if 'ARROWS' in draw_style:
        radius = 0.1
        # get material
        material = bpy.data.materials.get('axes')
        if not material:
            material = mb_utils.new_material('axes', color=(1,0,0))
        
        # add sphere at origin
        bpy.ops.mesh.primitive_uv_sphere_add(location=(0,0,0), size=radius, 
                                             segments=8, ring_count=8)
        ob = bpy.context.object
        ob.name = "unit_cell_origin"
        ob.parent = uc_cube
        ob.parent_type = 'VERTEX'
        ob.data.materials.append(None)
        ob.material_slots[0].material = material
        all_obs.append(ob)
        
        arrow_mesh = mb_utils.get_arrow_data(material=material, radius=radius)
        
        for di, vec in enumerate(unit_vectors):
            ob = bpy.data.objects.new(('a', 'b', 'c')[di], arrow_mesh)
            ob.parent = uc_cube
            ob.parent_type = 'VERTEX'
            bpy.context.scene.objects.link(ob)
            bpy.context.scene.objects.active = ob
            
            c = ob.constraints.new('STRETCH_TO')
            c.name = "mb.stretch"
            c.rest_length = 1.0
            c.volume = 'NO_VOLUME'
            c.target = uc_cube
            c.subtarget = vg[di].name
            all_obs.append(ob)
    
    if len(molecule.unit_cell_frames) > 1:
        print("animating uc")
        # add one shape key per frame
        for i, uvs in enumerate(len(molecule.unit_cell_frames)):
            shape_key = uc_cube.shape_key_add("{}.{}".format(uc_cube.name, i), 
                                              from_mix=False)
            unit_vectors = Matrix(uvs)
            coords = (
                (0,0,0), #0, O
                unit_vectors[0], #1, x
                unit_vectors[0] + unit_vectors[1], #2, xy
                unit_vectors[1], #3, y
                unit_vectors[2], #4, z
                unit_vectors[0] + unit_vectors[2], #5, xz
                unit_vectors[1] + unit_vectors[2], #6, yz
                unit_vectors[0] + unit_vectors[1] + unit_vectors[2], #7, xyz
                )
            for vert, coord in zip(uc_cube.data.vertices, coords):
                shape_key.data[vert.index].co = coord
        uc_cube.data.shape_keys.use_relative = False
        
        anim_data = uc_cube.data.shape_keys.animation_data_create()
        action = bpy.data.actions.new(name="uc_{}".format(molecule.index))
        anim_data.action = action
        
        fcu = action.fcurves.new(data_path="eval_time")
        for nf, kb in enumerate(uc_cube.data.shape_keys.key_blocks):
            fcu.keyframe_points.add(1)
            fcu.keyframe_points[-1].co = nf + 1, kb.frame
            fcu.keyframe_points[-1].interpolation = 'LINEAR'
        
    return all_obs

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
    
    first_atom = ob.mb.bonded_atoms[0].object
    second_atom = ob.mb.bonded_atoms[1].object
    
    first_mol = first_atom.mb.get_molecule()
    second_mol = second_atom.mb.get_molecule()
    
    first_mat = None
    second_mat = None
    
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


def update_radius_type(self, context):
    debug_print("mb_utils.update_radius_type", level=6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.object, self)


def set_draw_style(self, context):
    debug_print("mb_utils.set_draw_style", level=6)
    for atom in self.objects.atoms:
        set_atom_drivers(context, atom.object, self)
    
    hide = (self.draw_style == 'BALLS')
    for bond in self.objects.bonds:
        bond_ob = bond.object
        bond_ob.hide = hide
