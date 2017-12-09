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

import re
import math
import itertools

from mathutils import Vector, Matrix

from molblend.mb_helper import debug_print
from molblend.elements_default import ELEMENTS

A_per_Bohr = 0.529177249

class MB_Mode():
    
    def __init__(self):
        self.displacements = []
        self.freq = []
    

class MB_Q_Modes():
    """This class is a collection of all modes of a q-point.
    """
    def __init__(self, n_q=-1, qvec=None):
        self.modes = []
        self.n_q = n_q
        self.qvec = qvec

class MB_Mode_Collection():
    
    def __init__(self):
        self.qmodes = []
        
    @classmethod
    def from_file(cls, modefilepath, file_format):
        """
        Read vibrational modes from file, agnostic of the underlying atomic
        system.
        
        Currently implemented file formats:
        - qe_dynmat: Quantum espresso (dynmat/matdyn format)
        - simple: xyz style file with Nx3 (real) or Nx6 (complex) floats
            (i.e. real_modes.xyz:
                2
                freq = 10 cm^-1
                C  0.1  0.2  0.1
                O  0.1  0.9  0.1
                
                2
                freq = 0.1 THz
                C -0.3  0.2  0.1
                O  0.2 -0.1  0.2
                
                etc.
            )
            complex: Re_x Im_x Re_y Im_y Re_z Im_z
            The first column is completely ignored!
            This format doesn't allow separate q-points.
        """
        read_funcs = {
            "qe_dynmat": cls._read_qe_dynmat_out,
            "xyz": cls._read_simple_modes
        }
        if not file_format in read_funcs:
            msg = "ERROR: File format {}".format(file_format)
            msg += " not implemented yet."
            debug_print(msg, level=1)
        
        ret = read_modes_funcs[fmt](modefilepath)
        if not ret == True:
            debug_print("Problem when reading modes")
            structure.modes = []
            structure.freqs = []
    
    @classmethod
    def _read_simple_modes(cls, modefilepath):
        
        mode_col = cls()
        qmode = MB_Q_Modes()
        
        # This is a custom file format similar to the xyz file format
        with open(modefilepath, "r") as fin:
            number_atoms = -1
            for line in fin:
                
                # Simply ignore empty lines
                if line == "":
                    continue
                
                # Read until I find the number of atoms
                split_line = line.rsplit()
                if len(split_line) == 1:
                    try:
                        number_atoms = int(split_line[0])
                    except ValueError:
                        pass
                
                # Read as many lines as there are atoms (plus one comment line)
                # and repeat the procedure for additional frames
                if number_atoms > 0:
                    new_mode = MB_Mode()
                    # comment line might contain frequency, search for float
                    line = fin.readline()
                    m = re.search("[-+]?[0-9]*\.[0-9]+([eEdD][-+]?[0-9]+)?", line)
                    if m:
                        new_mode.freqs.append(m.group(0))
                    else:
                        new_mode.freqs.append("")
                    
                    for i in range(number_atoms):
                        
                        split_line = fin.readline().strip().split()
                        disp = list(map(float, split_line[1:]))
                        if len(disp) == 3:
                            real = disp
                            imag = (0., 0., 0.)
                        elif len(disp) == 6:
                            real = disp[::2]
                            imag = disp[1::2]
                        
                        new_mode.displacements.append((Vector(real), 
                                                       Vector(imag)))
                    
                    number_atoms = -1
                    qmode.modes.append(new_mode)
        mode_col.qmodes.append(qmode)
        return mode_col

    @classmethod
    def _read_qe_dynmat_out(cls, filepath):
        # read mode file, modes need to be in same order as atoms in input file
        # currently only supports dynmat.out
        
        mode_col = cls()

        with open(filepath, 'r') as fin:
            line = next(fin)
            q_count = 0
            while line:
                if 'q =' in line:
                    q_count += 1
                    q = list(map(float, line.split()[-3:]))
                qmodes = cls(n_q=q_count, qvec=Vector(q))
                debug_print("Reading q-point {}: ({}, {}, {})".format(n_q, *q), level=3)
                
                line = next(fin) # stars
                line = next(fin)
                
                new_mode = None
                for line in fin:
                    lstrip = line.strip()
                    # new mode
                    if lstrip.startswith('omega(') or lstrip.startswith('freq ('):
                        if new_mode:
                            qmodes.modes.append(new_mode)
                        new_mode = MB_Mode()
                        m = re.search(
                            '(omega|freq )\(([ 0-9*]+)\).+ ([-.0-9]+)(?= \[cm-1\])', 
                            lstrip)
                        new_mode.freq = float(m.group(3))
                    elif lstrip.startswith('('):
                        lsplit = lstrip[1:-1].split()
                        
                        disp = list(map(float, lsplit))
                        real = disp[::2]
                        imag = disp[1::2]
                        new_mode.displacements.append((Vector(real), 
                                                       Vector(imag)))
                    elif '**********' in line:
                        mode_col.qmodes.append(qmodes)
                        break
                line = next(fin)

        return mode_col


class MB_Structure():
    
    def __init__(self):
        self.nframes = 0
        self.all_atoms = {}
        self.bonds = {}
        self.axes = []
    
    def guess_bonds(self, tol=0.2):
        '''
        loops through all atoms and adds list of bonds to structure in the 
        format {index1: [index2, index3], index2: [index4, index5], ...}
        If bond argument is given, only additional bonds are added
        TODO for now, this only checks in the unit of Angstroms!
        '''
        
        for index1, atom1 in self.all_atoms.items():
            for index2, atom2 in self.all_atoms.items():
                if (index1 < index2 and not
                        (index1 in self.bonds
                         and index2 in self.bonds[index1])):
                    distance = (atom1["coords"][0] - atom2["coords"][0]).length
                    try:
                        cov1 = ELEMENTS[atom1["element"]]['covalent']
                    except KeyError:
                        cov1 = 1
                    try:
                        cov2 = ELEMENTS[atom1["element"]]['covalent']
                    except KeyError:
                        cov2 = 1
                    # the bond length should just be the sum of covalent radii 
                    # (plus some extra just to make sure)
                    max_dist = cov1 + cov2 + tol
                    if distance < max_dist:
                        try:
                            self.bonds[index1].append(index2)
                        except KeyError:
                            self.bonds[index1] = [index2]
            if not self.bonds:
                debug_print("WARNING: No bonds found.", level=1)
    
    def create_supercell(self, supercell):
        if self.axes:
            debug_print("Creating supercell.", level=4)
            max_atom_index = max(self.all_atoms) + 1
            
            n_unit_cell = 0
            iter_3d = itertools.product(range(supercell[0]),
                                        range(supercell[1]),
                                        range(supercell[2]))
            single_cell = list(self.all_atoms.items())
            for i, j, k in iter_3d:
                for index, atom in single_cell:
                    sc_coords = []
                    for loc, uc in zip(atom["coords"], self.axes):
                        sc_coords.append(
                            loc.copy() + i*uc[0] + j*uc[1] + k*uc[2]
                            )
                    #d = [element, atom_name, location, mode_vec]
                    new_index = index + max_atom_index*n_unit_cell
                    self.all_atoms[new_index] = {
                        "element": atom["element"],
                        "name": atom["name"],
                        "coords": sc_coords,
                        "id": index, # keep old index for reference
                        }
                n_unit_cell += 1
        else:
            debug_print("ERROR: Supercell requested, but no unit cell vectors"
                        + " read.", level=0)
    
    def get_center_of_mass(self):
        origin = Vector((0.0,0.0,0.0))
        coords = [atom["coords"][0] for atom in self.all_atoms.values()]
        center_of_mass = sum(coords, origin) / len(coords)
        return center_of_mass
    
    @classmethod
    def from_file(cls, filepath,
                  unit_fac=1.0):
        read_funcs = {
            "xyz": cls._read_xyz_file,
            "pdb": cls._read_pdb_file,
            "vasp": cls._read_vasp_file,
            "guo": cls._read_guo_file,
            "abinit": cls._read_abinit_output_file,
            "qe_input": cls._read_qe_input_file,
            "qe_output": cls._read_qe_rlx_output_file
            }
        
        fmt = filepath.rsplit('.')[-1]
        # Determine file format
        if not fmt in read_funcs:
            # read first lines of file to determine file format
            with open(filepath, 'r') as fin:
                lines = fin.read()
                if "ABINIT" in lines:
                    fmt = "abinit"
                if "&control" in lines.lower():
                    fmt = "qe_input"
                elif "BFGS Geometry Optimization" in lines:
                    fmt = "qe_output"
        
        structure = read_funcs[fmt](filepath)
        
        if unit_fac != 1.0:
            for atom in structure.all_atoms.values():
                atom["coords"] = [loc*unit_fac for loc in atom["coords"]]
        return structure
    
    @classmethod
    def _read_guo_file(cls, filepath_guo):
        
        strc = cls()
        
        # Open the file ...
        with open(filepath_guo, "r") as fin:
            # first line is name or comment
            #next(fin)
            molname = next(fin).strip()
            next(fin)
            strc.axes = [Vector(list(map(float, next(fin).split())))
                         for i in range(3)]
            element_list = next(fin).split()
            n_atoms = list(map(int, next(fin).split()))
            elements = [element_list[i] 
                        for i, n in enumerate(n_atoms) for j in range(n)]
            nat = sum(n_atoms)
            coord_type = next(fin).strip()
            if coord_type[0] in "dD":
                cell_matrix = Matrix(strc.axes).transposed()
            elif coord_type[0] in "cCkK":
                cell_matrix = Matrix.Identity(3)
            pos = [cell_matrix*Vector(list(map(float, next(fin).split())))
                   for i in range(nat)]
                
            for i, (element, location) in enumerate(zip(elements, pos)):
                atom_name = "{}{}".format(element.capitalize(), i)
                strc.all_atoms[i] = {"element": element,
                                     "name": atom_name,
                                     "coords": [location],
                                     "id": i}
        if not strc.all_atoms:
            raise IOError("No atoms found in {}".format(filepath_guo))
        strc.axes = [strc.axes]
        strc.nframes += 1
        return strc

    @classmethod
    def _read_vasp_file(cls, filepath_vasp):
        
        strc = cls()
        
        with open(filepath_vasp, "r") as fin:
            # first line is name or comment
            next(fin)
            scale = float(next(fin))
            # if scale is negative, it is the total volume of the cell
            if scale < 0.0:
                scale = 1.0
            # lattice vectors
            for i in range(3):
                line = next(fin)
                coords = list(map(float, line.split()))
                strc.axes.append(Vector(coords) * scale)
            # probably elements
            elements_list = next(fin).split()
            try:
                n_atoms = map(int, elements_list)
                elements = ["Default"] * sum(n_atoms)
            except ValueError:
                n_atoms = map(int, next(fin).split())
                elements = [elements_list[i] for i, n in enumerate(n_atoms)
                            for j in range(n)]
            line = next(fin).strip()
            # dynamical switch that is ignored here
            if line[0] in "Ss":
                line = next(fin).strip()
            # The next character determines format of coordinates.
            # - in cartesian coordinates
            if line[0] in "CcKk":
                cell_matrix = Matrix.Identity(3)
            # - fractional coordinates
            elif line[0] in "Dd":
                cell_matrix = Matrix(strc.axes).transposed()
            
            for i, element in enumerate(elements):
                line = next(fin)
                coords = list(map(float, line.split()))
                location = cell_matrix * Vector(coords)
                
                atom_name = "{}{}".format(element.capitalize(), i)
                strc.all_atoms[i] = {"element": element,
                                     "name": atom_name,
                                     "coords": [location],
                                     "id": i}
        strc.axes = [strc.axes]
        strc.nframes += 1
        return strc
    
    @classmethod
    def _read_xyz_file(cls, filepath_xyz):
        
        strc = cls()
        
        all_frames = []
        
        element_by_number = dict(
            [(ELEMENTS[element]["atomic number"], element) for element in ELEMENTS]
            )
        
        # XYZ is not a well defined file format. This is an attempt to be very
        # flexible. about it.
        with open(filepath_xyz, "r") as fin:
            number_atoms = -1
            for line in fin:
                print(line)
                # Simply ignore empty lines
                if line == "":
                    continue
                
                # Read until I find the number of atoms
                split_line = line.rsplit()
                if len(split_line) == 1:
                    try:
                        number_atoms = int(split_line[0])
                    except ValueError:
                        pass
                
                # Read as many lines as there are atoms (plus one comment line)
                # and repeat the procedure for additional frames
                if number_atoms > 0:
                    # dump comment line
                    line = next(fin)
                    print(line)
                    all_atoms = []
                    for i in range(number_atoms):
                        
                        split_line = next(fin).strip().split()
                        print(split_line)
                        coords = list(map(float, split_line[1:4]))
                        print(coords)
                        location = Vector(coords)
                        
                        element = split_line[0]
                        # Apparently some formats use atomic numbers instead of
                        # element names
                        try:
                            element = element_by_number[int(element)]
                        except:
                            pass
                        atom_name = "{}{}".format(element.capitalize(), i)
                        # If element is an 'X' then it is a vacancy.
                        if "X" in element:
                            element = "Vac"
                            atom_name = "Vacancy"
                        
                        all_atoms.append([element, atom_name, location, i])
                    
                    all_frames.append(all_atoms)
                    strc.nframes += 1
                    number_atoms = -1
        
        for element, atom_name, location, i in all_frames[0]:
            strc.all_atoms[i] = {"element": element,
                                 "name": atom_name,
                                 "coords": [location],
                                 "id": i}
        # Add other frames to coordinate lists
        for frame in all_frames[1:]:
            for el, an, loc, i in frame:
                if strc.all_atoms[i]["element"] == el:
                    strc.all_atoms[i]["coords"].append(loc)
                else:
                    raise ValueError(
                        "ERROR: Mismatch across frames in file {}.\n".format(
                            filepath_pdb
                            )
                        + "atom {}: {} (frame 0)".format(
                            i, strc.all_atoms[i]["element"]
                            )
                        + "vs {} (frame {})".format(el, nf+1),
                        level=0
                        )
        return strc
    
    @classmethod
    def _read_abinit_output_file(cls, filepath_abi):
        '''
        read abinit output file
        '''
        
        strc = cls()
        
        all_unit_vectors = []
        
        # some default variables
        optcell = 0
        ionmov = 0
        
        with open(filepath_abi, "r") as fin:
            for line in fin:
                if "ionmov" in line and re.search("ionmov +[0-9]+", line):
                    ionmov = int(line.split()[-1])
                elif ("acell" in line
                      and re.search("acell +( [ .0-9E+-]+){3}.*", line)):
                    acell_split = line.split()
                    acell = Vector([float(f) for f in acell_split[1:4]])
                    acell_unit = A_per_Bohr
                    if len(acell_split) == 5:
                        if acell_split[-1].startswith("Ang"):
                            acell_unit = 1.0
                        elif acell_split[-1] == "Bohr":
                            acell_unit = A_per_Bohr
                        else:
                            msg = ("WARNING: Didn't understand acell unit in "
                                   + "{} ({})".format(filepath_abi,
                                                      acell_split[-1]))
                            debug_print(msg, level=1)
                    acell *= acell_unit
                elif "natom" in line and re.search("natom +[0-9]+", line):
                    natom = int(line.split()[-1])
                elif "ndtset" in line and re.search("ndtset +[0-9]+", line):
                    ndtset = int(line.split()[-1])
                    if ndtset > 1:
                        debug_print("WARNING: More than one dataset present. "
                            "Will probably result in unwanted behavior.",
                            level=1)
                elif "optcell" in line and re.search("optcell +[0-9]+", line):
                    optcell = int(line.split()[-1])
                elif ("rprim" in line
                      and re.search("rprim +( [ .0-9E+-]+){3} *", line)):
                    rprim = []
                    ls = line.split()
                    rprim.append(acell[0]*Vector([float(f) for f in ls[1:4]]))
                    rprim.append(acell[1]*Vector([float(f)
                                                  for f in next(fin).split()]))
                    rprim.append(acell[2]*Vector([float(f)
                                                  for f in next(fin).split()]))
                    strc.axes.append(rprim)
                
                elif ("typat" in line
                      and re.search(" typat +( [0-9]+)+ *", line)):
                    try:
                        typat = []
                        ls = line.split()
                        typat.extend([int(i) for i in ls[1:]])
                        while len(typat) < natom:
                            typat.extend([int(i) for i in next(fin).split()])
                    except NameError as e:
                        msg = "ERROR: natom should be defined before typat"
                        debug_print(msg, level=1)
                        debug_print(e, level=1)
                elif ("xangst" in line
                      and re.search("xangst +( [ .0-9E+-]+){3} *", line)):
                    try:
                        xangst = []
                        ls = line.split()
                        xangst.append(Vector([float(f) for f in ls[1:4]]))
                        while len(xangst) < natom:
                            xangst.append(Vector([float(f) 
                                                  for f in next(fin).split()]))
                    except NameError:
                        msg = "ERROR: natom should be defined before xangst"
                        debug_print(msg, level=1)
                        debug_print(e, level=1)
                elif "znucl" in line and re.search("znucl +( [.0-9]+)+", line):
                    znucl = [int(i.split('.')[0]) for i in line.split()[1:]]
                    znucl_str = [""] * len(znucl)
                    for el, vals in ELEMENTS_DEFAULT.items():
                        if vals["atomic number"] in znucl:
                            znucl_str[znucl.index(vals["atomic number"])] = el
                    elements = [znucl_str[i-1] for i in typat]
                
                elif "== DATASET" in line:
                    # compile all information
                    strc.nframes += 1
                    for i, (el, location) in enumerate(zip(elements, xangst)):
                        atom_name = "{}{}".format(el.capitalize(), i)
                        strc.all_atoms[i] = {"element": el,
                                             "name": atom_name,
                                             "coords": [location],
                                             "id": i}
                    break
            if ionmov > 0:
                for line in fin:
                    if "Iteration" in line:
                        strc.nframes += 1
                        while "---OUTPUT" not in line:
                            line = next(fin)
                        next(fin)
                        next(fin)
                        for i in range(natom):
                            ls = next(fin).split()
                            coord = Vector([float(f) for f in ls]) * A_per_Bohr
                            strc.all_atoms[i]["coords"].append(coord)
                        
                        # now find cell
                        if optcell > 0:
                            stop = "Real space primitive translations (rprimd)"
                            while stop not in line:
                                line = next(fin)
                            pattern = "\(rprimd\) \[([a-z]+)\]"
                            unit = re.search(pattern, line).group(1)
                            if unit.lower() == "bohr":
                                fac = A_per_Bohr
                            elif unit.lower().startswith("ang"):
                                fac = 1.0
                            else:
                                msg = ("WARNING: Scale of Primitive Cell unit"
                                       + " not recognized ({})".format(unit))
                                debug_print(msg, level=1)
                                fac = 1.0
                            rprimd = []
                            for i in range(3):
                                xyz = [float(f) for f in next(fin).split()]
                                rprimd.append(Vector(xyz) * fac)
                            strc.axes.append(rprimd)
                    elif "== DATASET" in line:
                        break
        if strc.nframes > 1 and len(strc.axes) == 1:
            strc.axes = [strc.axes for i in range(strc.nframe)]
        
        return strc
        
    @classmethod
    def _read_pdb_file(cls, filepath_pdb):
        '''
        read pdb file (only supports the old standard with < 100k atoms)
        '''
        
        strc = cls()
        
        
        all_frames = []
        all_atoms = {}
        double = False
        with open(filepath_pdb, "r") as fin:
            i = -1
            for line in fin:
                
                if line == "":
                    continue
                
                if line[:6] == "CRYST1":
                    a, b, c, alpha, beta, gamma = map(float, line.split()[1:7])
                    alpha, beta, gamma = map(math.radians, 
                                             (alpha, beta, gamma))
                    avec = Vector((a, 0., 0.))
                    bvec = Vector((b*math.cos(gamma), b*math.sin(gamma), 0.))
                    cx = c*math.cos(beta)
                    cy = c*(math.cos(alpha) - math.cos(beta)*math.cos(gamma))
                    cz = math.sqrt(c*c - cx*cx - cy*cy)
                    cvec = Vector((cx,cy,cz))
                    strc.axes.append([avec, bvec, cvec])
                
                elif line[:6] == 'HETATM' or line[:4] == 'ATOM':
                    i += 1
                    coords = list(map(float,
                                      (line[30:38], line[38:46], line[46:54])))
                    location = Vector(coords)
                    
                    index = int(line[6:11])
                    element = line[76:78].strip().capitalize()
                    element = ''.join(s for s in element if s.isalpha())
                    atom_name = line[12:16].strip()
                    if "X" in element:
                        element = "Vac"
                    if index in all_atoms:
                        debug_print(
                            "ERROR: duplicate atom indeces in {}".format(
                                filepath_pdb
                                ),
                            level=0
                            )
                    all_atoms[index] = [element, atom_name, location, index]
                    
                elif line[:6] == 'CONECT':
                    # Need to make sure index is greater than 0
                    if int(line[6:11]) > 0:
                        atomID1 = int(line[6:11])
                        if not atomID1 in strc.bonds:
                            strc.bonds[atomID1] = []
                        for i in range((len(line) - 11) // 5):
                            atomID2 = int(line[11 + i*5 : 16 + i*5])
                            if atomID2 > 0 and atomID2 > atomID1:
                                if atomID2 not in strc.bonds[atomID1]:
                                    strc.bonds[atomID1].append(atomID2)
                
                elif line[:6] == 'ENDMDL':
                    debug_print("ENDMDL found.", level=6)
                    
                    all_frames.append(all_atoms)
                    strc.nframes += 1
                    all_atoms = {}
            
            # If there was only one model and no 'ENDMDL' present in file
            if all_atoms:
                all_frames.append(all_atoms)
                strc.nframes += 1
        
        for i in all_frames[0]:
            element, atom_name, location, index = all_frames[0][i]
            strc.all_atoms[i] = {"element": element,
                                 "name": atom_name,
                                 "coords": [location],
                                 "id": i}
        # Add other frames to coordinate lists
        for nf, frame in enumerate(all_frames[1:]):
            for i in frame:
                el, an, loc, index = frame[i]
                if strc.all_atoms[i]["name"] == an:
                    strc.all_atoms[i]["coords"].append(loc)
                else:
                    raise ValueError(
                        "ERROR: Mismatch across frames in file {}.\n".format(
                            filepath_pdb
                            )
                        + "atom {}: {}, {} (frame 0)".format(
                            i, strc.all_atoms[i]["element"],
                            strc.all_atoms[i]["name"]
                            )
                        + "vs {}, {} (frame {})".format(el, an, nf+1),
                        level=0
                        )
        
        if strc.nframes > 1 and len(strc.axes) == 1:
            strc.axes = [strc.axes for i in range(strc.nframe)]
        
        return strc
    
    @classmethod
    def _read_qe_input_file(cls, filepath):
        
        strc = cls()
        
        crystal_units = False
        with open(filepath, 'r') as fin:
            # Read all parameters
            for line in fin:
                for keyval in line.split(","):
                    if "nat" in keyval.lower():
                        na = keyval.split('=')[-1].strip().strip(",").strip()
                        n_atoms = int(na)
                    elif "celldm(1)" in keyval.lower():
                        al = keyval.split('=')[-1].strip().strip(",").strip()
                        alat = float(al)
                    elif "CELL_PARAMETERS" in keyval.upper():
                        if re.search("bohr|cubic", line):
                            unit = A_per_Bohr
                        else:
                            unit = 1.0
                        for i in range(3):
                            line = next(fin)
                            coords = list(map(float, line.split()))
                            strc.axes.append(Vector(coords) * unit)
                    elif "ATOMIC_POSITIONS" in keyval.upper():
                        # set correct conversion matrix
                        cell_matrix = Matrix.Identity(3)
                        if "angstrom" in line.lower():
                            cell_matrix = 1.0
                            #if scale_distances != 1.0:
                                #debug_print(
                                    #"Found coordinates in Angstrom. Overriding unit.",
                                    #level=1)
                        elif "bohr" in line.lower():
                            cell_matrix = A_per_Bohr * cell_matrix
                            #if scale_distances != A_per_Bohr:
                                #debug_print(
                                    #"Found coordinates in Bohr. Overriding unit.", 
                                    #level=1)
                        elif "crystal" in line.lower():
                            # process at the end when both atomic coordinates and 
                            # vectors are read for certain
                            crystal_units = True
                            # get unit vectors
                        elif "alat" in line.lower():
                            msg = ("Unit of ATOMIC_POSITIONS in "
                                   + "{} is alat. Implement!".format(filepath))
                            raise ValueError(msg)
                            #alat = int(line.split('=')[-1].strip().strip(",").strip())
                            #cell_matrix = alat * A_per_Bohr * cell_matrix
                        else:
                            msg = ("Unit of ATOMIC_POSITIONS in "
                                   + "{} unclear.\n".format(filepath)
                                   + "Line:\n'{}'".format(line))
                            raise ValueError(msg)
                        
                        # read coordinates
                        for i in range(n_atoms):
                            line = next(fin)
                            line = line.strip()
                            split_line = line.split()
                            coords = list(map(float, split_line[1:4]))
                            location = cell_matrix * Vector(coords)

                            element = split_line[0]
                            atom_name = "{}{}".format(element.capitalize(), i)
                            # If element is an 'X' then it is a vacancy.
                            if "X" in element:
                                element = "Vac"
                                atom_name = "Vacancy"
                            
                            strc.all_atoms[i] = {"element": element,
                                                 "name": atom_name,
                                                 "coords": [location],
                                                 "id": i}
            if crystal_units:
                cell_matrix = Matrix(strc.axes).transposed()
                for i in strc.all_atoms:
                    location = cell_matrix * strc.all_atoms[i]["coords"][0]
                    strc.all_atoms[i]["coords"][0] = location

        strc.axes = [strc.axes]
        strc.nframes += 1
        return strc
    
    @classmethod
    def _read_qe_rlx_output_file(cls, filepath):
        
        strc = cls()
        
        all_atoms = {}
        all_frames = []
        all_unit_vectors = []
        with open(filepath, 'r') as fin:
            # read all parameters
            for line in fin:
                if "number of atoms/cell" in line:
                    n_atoms = int(line.split()[-1])
                #if "lattice parameter (alat)" in line:
                    #if line.split()[-1] == 'a.u.':
                        #fac = A_per_Bohr 
                    #else:
                        #raise ValueError("check units in {}".format(filepath))
                    #alat = float(line.split()[-2]) * fac
                elif "lattice parameter (alat)" in line:
                    alat = A_per_Bohr * float(line.split()[-2])
                elif "crystal axes: (cart. coord. in units of alat)" in line:
                    unit_vectors = []
                    for i in range(3):
                        line = fin.readline()
                        coords = list(map(float, line.split()[3:6]))
                        unit_vectors.append(Vector(coords) * alat)
                    all_unit_vectors.append(unit_vectors)
                elif "Cartesian axes" in line:
                    break
            
            # Read starting geometry
            next(fin) # empty line
            next(fin) # column headers
            cell_matrix = alat * Matrix.Identity(3)
            # read coordinates
            strc.nframes += 1
            for i in range(n_atoms):
                line = next(fin)
                line = line.strip()
                split_line = line.split()
                coords = list(map(float, split_line[-4:-1]))
                location = cell_matrix * Vector(coords)

                element = split_line[1]
                atom_name = "{}{}".format(element.capitalize(), i)
                # If element is an 'X' then it is a vacancy.
                if "X" in element:
                    element = "Vac"
                    atom_name = "Vacancy"
                
                strc.all_atoms[i] = {"element": element,
                                     "name": atom_name,
                                     "coords": [location],
                                     "id": i}
            
            # Now read all steps
            for line in fin:
                if "CELL_PARAMETERS" in line:
                    # vc-relax calculation
                    if "alat" in line:
                        pattern = r"CELL_PARAMETERS \(alat=(.+)\)$"
                        alat = float(re.match(pattern,
                                            line).group(1)) * A_per_Bohr
                    elif "angstrom" in line:
                        alat = 1.0
                    unit_vectors = []
                    for i in range(3):
                        line = fin.readline()
                        coords = list(map(float, line.split()))
                        unit_vectors.append(Vector(coords) * alat)
                    all_unit_vectors.append(unit_vectors)
                
                if "ATOMIC_POSITIONS" in line:
                    strc.nframes += 1
                    # double check units and set correct conversion matrix
                    if "angstrom" in line.lower():
                        cell_matrix = 1.0 #Matrix.Identity(3)
                        #if scale_distances != 1.0:
                            #debug_print(
                                #"Found coordinates in Angstrom. Overriding unit.", 
                                #level=1)
                    elif "bohr" in line.lower():
                        cell_matrix = A_per_Bohr #* Matrix.Identity(3)
                        #if scale_distances != A_per_Bohr:
                            #debug_print(
                                #"Found coordinates in Bohr. Overriding unit.",
                                #level=1)
                    elif "crystal" in line.lower():
                        # get unit vectors
                        cell_matrix = Matrix(all_unit_vectors[-1]).transposed()
                    elif "alat" in line.lower():
                        #alat = float(line.split('=')[-1].strip().strip(",").strip()) 
                        cell_matrix = alat #* Matrix.Identity(3)
                    else:
                        msg = ("Unit of ATOMIC_POSITIONS in "
                               + "{} unclear.\n".format(filepath)
                               + "Line:\n'{}'".format(line))
                        raise ValueError(msg)
                    
                    # read coordinates
                    for i in range(n_atoms):
                        line = next(fin)
                        line = line.strip()
                        split_line = line.split()
                        coords = list(map(float, split_line[1:4]))
                        location = cell_matrix * Vector(coords)

                        element = split_line[0]
                        atom_name = "{}{}".format(element.capitalize(), i)
                        # If element is an 'X' then it is a vacancy.
                        if "X" in element:
                            element = "Vac"
                            atom_name = "Vacancy"
                        
                        if strc.all_atoms[i]["element"] == element:
                            strc.all_atoms[i]["coords"].append(location)
                        else:
                            raise ValueError(
                                "ERROR: Mismatch across frames in file "
                                + "{}.\n".format(filepath)
                                + "atom {}: {} (frame 0)".format(
                                    i, strc.all_atoms[i]["element"]
                                    )
                                + "vs {} (frame {})".format(element,
                                                            strc.nframes),
                                level=0
                                )
        
        return strc



## TODO combine export functions into one and split different fileformats into
## write functions
#def export_xyz(filepath, selection_only, scale_distances):
    #debug_print("Write to file {}.".format(filepath), level=1)
    #all_atoms = {}
    #if selection_only:
        #objects = bpy.context.selected_objects
    #else:
        #objects = bpy.context.scene.objects
    
    #n_atoms = 0
    #for ob in objects:
        #if ob.mb.type == 'ATOM':
            #try:
                #all_atoms[ob.mb.molecule_name].append(
                    #(ob.mb.index, ob.mb.element, get_world_coordinates(ob)))
            #except KeyError:
                #all_atoms[ob.mb.molecule_name] = [(ob.mb.index, ob.mb.element, 
                                                  #get_world_coordinates(ob))]
            #n_atoms += 1
    
    #with open(filepath, "w") as fout:
        #fout.write("{}\n".format(n_atoms))
        #fout.write("This file has been created with Blender "
                   #"and the MolBlend addon.\n")
    
        #for mol_id in sorted(all_atoms):
            #for i, element, location in sorted(all_atoms[mol_id]):
                #l = "{}   {:>10.5f}   {:>10.5f}   {:>10.5f}\n"
                #fout.write(l.format(element, *location))
        #debug_print("Exported {} atoms.".format(n_atoms), level=1)
    #return True


#def export_pdb(filepath, selection_only, scale_distances):
    #debug_print("Write to file {}.".format(filepath), level=1)
    #all_atoms = {}
    #if selection_only:
        #objects = bpy.context.selected_objects
    #else:
        #objects = bpy.context.scene.objects
    
    #n_atoms = 0
    #for ob in objects:
        #if ob.mb.type == 'ATOM':
            #atom_name = ob.mb.atom_name
            #info = (
                #ob.mb.index,
                #atom_name if len(ob.mb.element) != 1 else " " + atom_name,
                #ob.mb.get_molecule().name[:3].upper(),
                #get_world_coordinates(ob),
                #ob.mb.element,
                #)
            #try:
                #all_atoms[ob.mb.molecule_name].append(info)
            #except KeyError:
                #all_atoms[ob.mb.molecule_name] = [info]
            #n_atoms += 1
        
    #with open(filepath, "w") as fout:
        #fout.write("REMARK This pdb file has been created with Blender "
                  #"and the addon MolBlend\n"
                  #"REMARK\n"
                  #"REMARK\n")
        #index = 1
        #for mol_id in sorted(all_atoms):
            #for data in sorted(all_atoms[mol_id]):
                #i, name, res_name, coords, element = data
                #l = ("HETATM{ID:>5} {name:<4} {res}     1    {c[0]:8.3f}" 
                     #"{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00          "
                     #"{element:>2}\n")
                #fout.write(l.format(ID=index, name=name[:4], res=res_name, 
                                    #c=coords, element=element))
                #index += 1
                #if index > 99999:
                    #index = 1
    #debug_print("Exported {} atoms.".format(n_atoms), level=1)
    #return True


#def export_b4w(filepath):
    ## need to link all materials to meshes instead of objects
    ## so first make each mesh single user
    ## then link the material slots to data
    ## export to b4w without autosaving blend file
    ## undo all the changes and revert to state before (or alternatively copy
    ## selection to new file and do everything there)
    #debug_print("Blend4Web export not implemented yet.", level=1)
