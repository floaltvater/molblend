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

from molblend.mb_helper import debug_print

def read_guo_file(filepath_guo, scale_distances=1.0, mask_planes=None, 
                  mask_flip=False):
    mask_planes = mask_planes or []
    number_frames = 0
    all_atoms = {}
    unit_vectors = []
    # Open the file ...
    with open(filepath_guo, "r") as fin:
        # first line is name or comment
        #next(fin)
        molname = next(fin).strip()
        next(fin)
        unit_vectors = [Vector(list(map(float, next(fin).split()))) * scale_distances for i in range(3)]
        element_list = next(fin).split()
        n_atoms = list(map(int, next(fin).split()))
        elements = [element_list[i] for i, n in enumerate(n_atoms) for j in range(n)]
        nat = sum(n_atoms)
        coord_type = next(fin).strip()
        if coord_type[0] in "dD":
            cell_matrix = Matrix(unit_vectors).transposed()
        elif coord_type[0] in "cCkK":
            cell_matrix = Matrix.Identity(3)
        pos = [cell_matrix*Vector(list(map(float, next(fin).split()))) for i in range(nat)]
            
        for i, (element, location) in enumerate(zip(elements, pos)):
            if (mask_planes and
                not is_inside_of_planes(mask_planes, location,
                                        flip=mask_flip)):
                continue
            atom_name = element + str(i)
            all_atoms[i] = [element, atom_name, location, i]
    if not all_atoms:
        raise IOError("no atoms found")
    return [all_atoms], unit_vectors

def read_guo_modes(filepath):
    all_evecs = []
    freqs = []
    with open(filepath, "r") as fin:
        
        #next(fin)
        molname = next(fin).strip()
        next(fin)
        [next(fin) for i in range(3)]
        next(fin)
        n_atoms = list(map(int, next(fin).split()))
        nat = sum(n_atoms)
        
        for line in fin:
            if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
                next(fin)
                break
#         for line in fin:
#             if "-----" in line:
#                 break
        
        for line in fin:
            if "cm-1" in line:
                lsplit = line.split()
                mode_n = int(lsplit[0])
                freqs.append(float(lsplit[-4]))
                current = []
                all_evecs.append(current)
                next(fin)
                for n in range(nat):
                    ls = next(fin).strip().split()
                    current.append(list(map(float, ls[3:])))
            elif "Eigenvector" in line:
                break
    return freqs, all_evecs


def read_vasp_file(filepath_vasp, scale_distances=1.0, mask_planes=None, 
                  mask_flip=False):
    mask_planes = mask_planes or []
    number_frames = 0
    all_atoms = {}
    unit_vectors = []
    # Open the file ...
    with open(filepath_vasp, "r") as fin:
        # first line is name or comment
        next(fin)
        # scaling factor
        scale = float(next(fin))
        # if scale is negative, it is the total volume of the cell
        if scale < 0.0:
            scale = 1.0
        # lattice vectors
        for i in range(3):
            line = next(fin)
            coords = list(map(float, line.split()))
            unit_vectors.append(Vector(coords) * scale)
        # probably elements
        elements_list = next(fin).split()
        try:
            n_atoms = map(int, elements_list)
            elements = ["Default"] * sum(n_atoms)
        except ValueError:
            n_atoms = map(int, next(fin).split())
            elements = [elements_list[i] for i, n in enumerate(n_atoms) for j in range(n)]
        line = next(fin).strip()
        # dynamical switch that is ignored here
        if line[0] in "Ss":
            line = next(fin).strip()
        # Cartesian coordinates
        if line[0] in "CcKk":
            cell_matrix = Matrix.Identity(3)
        # fractional coordinates
        elif line[0] in "Dd":
            cell_matrix = Matrix(unit_vectors).transposed()
        
        for i, element in enumerate(elements):
            line = next(fin)
            coords = list(map(float, line.split()))
            location = cell_matrix * Vector(coords)
            # skip processing this coordinate if outside of mask
            if (mask_planes and 
                not is_inside_of_planes(mask_planes, location,
                                        flip=mask_flip)):
                continue
            atom_name = element + str(i)
            all_atoms[i] = [element, atom_name, location, i]
    
    return [all_atoms], unit_vectors

def read_xyz_file(filepath_xyz, scale_distances=1.0, mask_planes=None, 
                  mask_flip=False):
    '''
    read xyz file and return all atoms in the format
    [{index1: [element, atom_name, coords], index2: [...], ...}, [2nd frame],
     ...]
    '''
    mask_planes = mask_planes or []
    number_frames = 0
    all_frames = []
    
    truncate_manually = False
    if truncate_manually:
        print("WARNING: Truncating structure!")
        xlim = (0, 100)
        ylim = (0, 100)
        zlim = (0, 15.5)
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
                    coords = list(map(float, split_list[1:4]))
                    location = Vector(coords) * scale_distances
                    # skip processing this coordinate if outside of mask
                    if (mask_planes and 
                        not is_inside_of_planes(mask_planes, location, 
                                                flip=mask_flip)):
                        continue
                    ## WARNING comment out
                    if truncate_manually:
                        if not (xlim[0] <= location[0] <= xlim[1] and
                                ylim[0] <= location[1] <= ylim[1] and
                                zlim[0] <= location[2] <= zlim[1]):
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

def read_abinit_output_file(filepath_abi, scale_distances=1.0, mask_planes=None, 
                     mask_flip=False):
    '''
    read abinit output file and return all atoms in the format
    [{index1: [element, atom_name, coords], index2: [...], ...}, [2nd frame],
     ...]
    '''
    mask_planes = mask_planes or []
    all_frames = []
    all_unit_vectors = []
    
    truncate_manually = False
    if truncate_manually:
        print("WARNING: Truncating structure!")
        xlim = (0, 100)
        ylim = (0, 100)
        zlim = (0, 15.5)
    # Open the file ...
    
    # some default variables
    optcell = 0
    ionmov = 0
    
    with open(filepath_abi, "r") as fin:
        for line in fin:
            if "ionmov" in line and re.search("ionmov +[0-9]+", line):
                ionmov = int(line.split()[-1])
            elif "acell" in line and re.search("acell +( [ .0-9E+-]+){3}.*", line):
                acell_split = line.split()
                acell = Vector([float(f) for f in acell_split[1:4]])
                acell_unit = A_per_Bohr
                if len(acell_split) == 5:
                    if acell_split[-1].startswith("Ang"):
                        acell_unit = 1.0
                    elif acell_split[-1] == "Bohr":
                        acell_unit = A_per_Bohr
                    else:
                        debug_print("WARNING: Didn't understand acell unit in " 
                                    "{} ({})".format(filepath_abi, acell_split[-1]),
                                    level=1)
                acell *= acell_unit
            elif "natom" in line and re.search("natom +[0-9]+", line):
                natom = int(line.split()[-1])
            elif "ndtset" in line and re.search("ndtset +[0-9]+", line):
                ndtset = int(line.split()[-1])
                if ndtset > 1:
                    debug_print("WARNING: More than one dataset present. "
                        "Will probably result in unwanted behavior.", level=1)
            elif "optcell" in line and re.search("optcell +[0-9]+", line):
                optcell = int(line.split()[-1])
            elif "rprim" in line and re.search("rprim +( [ .0-9E+-]+){3} *", line):
                rprim = []
                ls = line.split()
                rprim.append(acell[0]*Vector([float(f) for f in ls[1:4]]))
                rprim.append(acell[1]*Vector([float(f) for f in next(fin).split()]))
                rprim.append(acell[2]*Vector([float(f) for f in next(fin).split()]))
                all_unit_vectors.append(rprim)
            
            elif "typat" in line and re.search(" typat +( [0-9]+)+ *", line):
                try:
                    typat = []
                    ls = line.split()
                    typat.extend([int(i) for i in ls[1:]])
                    while len(typat) < natom:
                        typat.extend([int(i) for i in next(fin).split()])
                except NameError as e:
                    debug_print("ERROR: natom should be defined before typat", level=1)
                    debug_print(e, level=1)
            elif "xangst" in line and re.search("xangst +( [ .0-9E+-]+){3} *", line):
                try:
                    xangst = []
                    ls = line.split()
                    xangst.append(Vector([float(f) for f in ls[1:4]]))
                    while len(xangst) < natom:
                        xangst.append(Vector([float(f) for f in next(fin).split()]))
                except NameError:
                    debug_print("ERROR: natom should be defined before xangst")
                    debug_print(e, level=1)
            elif "znucl" in line and re.search("znucl +( [.0-9]+)+", line):
                znucl = [int(i.split('.')[0]) for i in line.split()[1:]]
                znucl_str = [""] * len(znucl)
                for element, vals in ELEMENTS_DEFAULT.items():
                    if vals["atomic number"] in znucl:
                        znucl_str[znucl.index(vals["atomic number"])] = element
                elements = [znucl_str[i-1] for i in typat]
            
            elif "== DATASET" in line:
                # compile all information
                all_atoms = {}
                for i, (element, location) in enumerate(zip(elements, xangst)):
                    atom_name = element + str(i)
                    all_atoms[i] = [element, atom_name, location, i]
                
                all_frames.append(all_atoms)
                break
        if ionmov > 0:
        # Dataset 1
            for line in fin:
                if "Iteration" in line:
                    frame = int(re.search("Iteration: \(([ 0-9]+)/[ 0-9]+\)", line).group(1))
                    all_atoms = {}
                elif "---OUTPUT" in line:
                    ref = all_frames[0]
                    next(fin)
                    next(fin)
                    for i in range(natom):
                        coords = Vector([float(f) for f in next(fin).split()]) * A_per_Bohr
                        all_atoms[i] = [ref[i][0], ref[i][1], coords, i]
                    all_frames.append(all_atoms)
                    
                    # now find cell
                    if optcell > 0:
                        while not "Real space primitive translations (rprimd)" in line:
                            line = next(fin)
                            
                        unit = re.search("\(rprimd\) \[([a-z]+)\]", line).group(1)
                        if unit.lower() == "bohr":
                            fac = A_per_Bohr
                        elif unit.lower().startswith("ang"):
                            fac = 1.0
                        else:
                            debug_print("WARNING: Scale of Primitive Cell unit not recognized"
                                " ({})".format(unit), level=1)
                            fac = 1.0
                        rprimd = []
                        for i in range(3):
                            rprimd.append(Vector([float(f) for f in next(fin).split()]) * fac)
                        all_unit_vectors.append(rprimd)
                elif "== DATASET" in line:
                    break
    if len(all_frames) > 1 and len(all_unit_vectors) == 1:
        all_unit_vectors = all_unit_vectors[0]
    return all_frames, all_unit_vectors


def read_pdb_file(filepath_pdb, scale_distances=1.0, mask_planes=None, 
                  mask_flip=False):
    '''
    read pdb file and return all atoms in the format
    [{index1: [element, atom_name, coords], index2: [...], ...}, [2nd frame],
     ...]
    and bonds (index1 < index2 < index 3 < ...)
    {index1: [index2, index3], index2: [index4, index5], ...}
    '''
    mask_planes = mask_planes or []
    
    number_frames = 0
    
    all_frames = []
    all_atoms = {}
    all_unit_vectors = []
    bonds = {}
    # TODO include as GUI option
    # if double == False: exclude double bonds
    double = False
    with open(filepath_pdb, "r") as fin:
        i = -1
        for line in fin:
            
            if line == "":
                continue
            
            if line[:6] == "CRYST1":
                a, b, c, alpha, beta, gamma = map(float, line.split()[1:7])
                # convert to radians
                alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
                avec = Vector((a, 0., 0.))
                bvec = Vector((b*math.cos(gamma), b*math.sin(gamma), 0.))
                cx = c*math.cos(beta)
                cy = c*(math.cos(alpha) - math.cos(beta)*math.cos(gamma))
                cz = math.sqrt(c*c - cx*cx - cy*cy)
                cvec = Vector((cx,cy,cz))
                all_unit_vectors.append([avec, bvec, cvec])
            
            # get atom information
            elif line[:6] == 'HETATM' or line[:4] == 'ATOM':
                i += 1
                coords = list(
                    map(float, (line[30:38], line[38:46], line[46:54])))
                location = Vector(coords) * scale_distances
                # skip processing this coordinate if outside of mask
                if (mask_planes and
                    not is_inside_of_planes(mask_planes, location, 
                                            flip=mask_flip)):
                    continue
                
                # get atom number
                index = int(line[6:11])
                # get atom type, remove whitespace and capitalize first letter
                element = line[76:78].strip().capitalize()
                element = ''.join(i for i in element if (not i.isdigit() and 
                                                         not i in ('+', '-')))
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
                    for i in range((len(line) - 11) // 5):
                        # get connected atomID
                        atomID2 = int(line[11 + i*5 : 16 + i*5])
                        # only store bond once
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
    
    if len(all_frames) > 1 and len(all_unit_vectors) == 1:
        all_unit_vectors = all_unit_vectors[0]
    
    return all_frames, all_unit_vectors, bonds


def read_qe_input_file(filepath, scale_distances=1.0, mask_planes=None, 
                       mask_flip=False):
    mask_planes = mask_planes or []
    all_atoms = {}
    unit_vectors = []
    crystal_units = False
    with open(filepath, 'r') as fin:
        # Read all parameters
        for line in fin:
            for keyval in line.split(","):
                if "nat" in keyval.lower():
                    n_atoms = int(keyval.split('=')[-1].strip().strip(",").strip())
                elif "celldm(1)" in keyval.lower():
                    alat = float(keyval.split('=')[-1].strip().strip(",").strip())
                elif "CELL_PARAMETERS" in keyval.upper():
                    unit = A_per_Bohr if re.search("bohr|cubic", line) else 1.0
                    for i in range(3):
                        line = next(fin)
                        coords = list(map(float, line.split()))
                        unit_vectors.append(Vector(coords) * unit)
                elif "ATOMIC_POSITIONS" in keyval.upper():
                    # set correct conversion matrix
                    cell_matrix = Matrix.Identity(3)
                    if "angstrom" in line.lower():
                        # Identity is already correct
                        if scale_distances != 1.0:
                            debug_print(
                                "Found coordinates in Angstrom. Overriding unit.",
                                level=1)
                    elif "bohr" in line.lower():
                        cell_matrix = A_per_Bohr * cell_matrix
                        if scale_distances != A_per_Bohr:
                            debug_print(
                                "Found coordinates in Bohr. Overriding unit.", 
                                level=1)
                    elif "crystal" in line.lower():
                        # process at the end when both atomic coordinates and 
                        # vectors are read for certain
                        crystal_units = True
                        # get unit vectors
                    elif "alat" in line.lower():
                        raise ValueError("Unit of ATOMIC_POSITIONS in {} is alat. Implement!".format(filepath))
                        #alat = int(line.split('=')[-1].strip().strip(",").strip())
                        #cell_matrix = alat * A_per_Bohr * cell_matrix
                    else:
                        raise ValueError("Unit of ATOMIC_POSITIONS in {} unclear.\nLine:\n'{}'".format(filepath, line))
                    
                    # read coordinates
                    for i in range(n_atoms):
                        line = next(fin)
                        line = line.strip()
                        split_list = line.split()
                        coords = list(map(float, split_list[1:4]))
                        location = cell_matrix * Vector(coords)

                        # skip processing this coordinate if outside of mask
                        if (mask_planes and 
                            not is_inside_of_planes(mask_planes, location,
                                                    flip=mask_flip)):
                            continue
                        
                        element = split_list[0]
                        atom_name = element + str(i)
                        # If element is an 'X' then it is a vacancy.
                        if "X" in element:
                            element = "Vac"
                            atom_name = "Vacancy"
                        
                        all_atoms[i] = [element, atom_name, location, i]
        if crystal_units:
            cell_matrix = Matrix(unit_vectors).transposed()
            for i, (element, atom_name, location, i) in all_atoms.items():
                location = cell_matrix * location
                all_atoms[i][2] = location
    return [all_atoms], unit_vectors


def read_qe_rlx_output_file(filepath, scale_distances=1.0, mask_planes=None, 
                            mask_flip=False):
    
    mask_planes = mask_planes or []
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
        for i in range(n_atoms):
            line = next(fin)
            line = line.strip()
            split_list = line.split()
            coords = list(map(float, split_list[-4:-1]))
            location = cell_matrix * Vector(coords)

            # skip processing this coordinate if outside of mask
            if (mask_planes and 
                not is_inside_of_planes(mask_planes, location,
                                        flip=mask_flip)):
                continue
            
            element = split_list[1]
            atom_name = element + str(i)
            # If element is an 'X' then it is a vacancy.
            if "X" in element:
                element = "Vac"
                atom_name = "Vacancy"
            
            all_atoms[i] = [element, atom_name, location, i]
        all_frames.append(all_atoms)
        all_atoms = {}
        
        # Now read all steps
        for line in fin:
            if "CELL_PARAMETERS" in line:
                # vc-relax calculation
                if "alat" in line:
                    alat = float(re.match(r"CELL_PARAMETERS \(alat=(.+)\)$", 
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
                # double check units and set correct conversion matrix
                if "angstrom" in line.lower():
                    cell_matrix = 1.0 #Matrix.Identity(3)
                    if scale_distances != 1.0:
                        debug_print(
                            "Found coordinates in Angstrom. Overriding unit.", 
                            level=1)
                elif "bohr" in line.lower():
                    cell_matrix = A_per_Bohr #* Matrix.Identity(3)
                    if scale_distances != A_per_Bohr:
                        debug_print(
                            "Found coordinates in Bohr. Overriding unit.",
                            level=1)
                elif "crystal" in line.lower():
                    # get unit vectors
                    cell_matrix = Matrix(all_unit_vectors[-1]).transposed()
                elif "alat" in line.lower():
                    #alat = float(line.split('=')[-1].strip().strip(",").strip()) 
                    cell_matrix = alat #* Matrix.Identity(3)
                else:
                    raise ValueError("Unit of ATOMIC_POSITIONS in {} unclear.\nLine:\n'{}'".format(filepath, line))
                
                # read coordinates
                for i in range(n_atoms):
                    line = next(fin)
                    line = line.strip()
                    split_list = line.split()
                    coords = list(map(float, split_list[1:4]))
                    location = cell_matrix * Vector(coords)

                    # skip processing this coordinate if outside of mask
                    if (mask_planes and 
                        not is_inside_of_planes(mask_planes, location,
                                                flip=mask_flip)):
                        continue
                    
                    element = split_list[0]
                    atom_name = element + str(i)
                    # If element is an 'X' then it is a vacancy.
                    if "X" in element:
                        element = "Vac"
                        atom_name = "Vacancy"
                    
                    all_atoms[i] = [element, atom_name, location, i]
                all_frames.append(all_atoms)
                all_atoms = {}
    debug_print("Found {} frames".format(len(all_frames)), level=1)
    return all_frames, all_unit_vectors


#def read_qe_unit_cell(filepath):
    #unit_vectors = []
    #with open(filepath, 'r') as fin:
        #for line in fin:
            #if "celldm(1)" in line.lower():
                #alat = float(line.split('=')[-1].strip().strip(",").strip())
            #if "CELL_PARAMETERS" in line:
                ## TODO include alat units etc.
                #if re.search("bohr|cubic", line.lower()):
                    #unit = A_per_Bohr
                #elif "angstrom" in line.lower():
                    #unit = 1.0
                #elif "alat" in line.lower():
                    #unit = alat * A_per_Bohr
                #for i in range(3):
                    #line = next(fin)
                    #coords = list(map(float, line.split()))
                    #unit_vectors.append(Vector(coords) * unit)
                #break
    #return unit_vectors

def read_phonopy_mode(filepath, all_frames):
    mass_dict = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "Si": 28.085,
        "S": 32.06,
        }
    masses = [mass_dict[all_frames[0][i][0]] for i in sorted(all_frames[0])]
    print(masses)
    print(len(masses))
    all_evecs = []
    freqs = []
    with open(filepath, "r") as fin:
        nqpt = int(next(fin).split()[-1])
        nat = int(next(fin).split()[-1])
        
        # reciprocal lattice
        next(fin)
        #recip_unit_vectors = []
        for i in range(3):
            line = next(fin)
            #p = "\[([- .0-9]+),([- .0-9]+),([- .0-9]+)\]"
            #coords = list(map(float, re.search(p, line).group(1,2,3)))
            #recip_unit_vectors.append(Vector(coords) * alat)
        # phonon
        next(fin)
        for nq in range(nqpt):
            q = next(fin)
            next(fin) #band
            modes = []
            for nm in range(3*nat):
                nmode = int(next(fin).split()[-1])
                freq = float(next(fin).split()[-1])
                next(fin) # eigenvector
                mode = []
                for na in range(nat):
                    next(fin)
                    coords = []
                    try:
                        for i in range(3):
                            line = next(fin)
                            p = "\[([- .0-9]+),([- .0-9]+)\]"
                            real, im = list(map(float, re.search(p, line).group(1,2)))
                            coords.append(real/math.sqrt(masses[na]))
                    except AttributeError:
                        print(line)
                        raise
                    mode.append(coords)
                freqs.append(freq)
                all_evecs.append(mode)
    for i in range(3):
        for c in all_evecs[i]:
            print("{:>8.5f} {:>8.5f} {:>8.5f}".format(*c))
    return freqs, all_evecs

def read_manual_modes(filepath):
    all_evecs = []
    with open(filepath, 'r') as fin:
        for line in fin:
            ls = list(map(float, line.split()))
            all_evecs.append([[ls[i], ls[i+1], ls[i+2]] for i in range(0, len(ls), 3)])
    with open(filepath.rsplit(".",1)[0] + ".freqs", "r") as fin:
        freqs = [float(line) for line in fin]
            
    return freqs, all_evecs

def read_modes(filepath, n_q=1):
    # read mode file, modes need to be in same order as atoms in input file
    # currently only supports dynmat.out
    with open(filepath, 'r') as fin:
        line = next(fin)
        q_count = 0
        while line and q_count != n_q:
            if 'q =' in line:
                q_count += 1
                q = list(map(float, line.split()[-3:]))
            line = next(fin)
        debug_print("Reading q-point {}: ({}, {}, {})".format(n_q, *q), level=2)
        all_evecs = []
        freqs = []
        for line in fin:
            lstrip = line.strip()
            # new mode
            if lstrip.startswith('omega(') or lstrip.startswith('freq ('):
                m = re.search(
                    '(omega|freq )\(([ 0-9*]+)\).+ ([-.0-9]+)(?= \[cm-1\])', 
                    lstrip)
                freq = float(m.group(3))
                freqs.append(freq)
                current = []
                all_evecs.append(current) # links current to all_evecs
            elif lstrip.startswith('('):
                lsplit = lstrip[1:-1].split()
                current.append(list(map(float, lsplit[::2])))
            elif 'q = ' in line:
                break
    
    return freqs, all_evecs

# TODO combine export functions into one and split different fileformats into
# write functions
def export_xyz(filepath, selection_only, scale_distances):
    debug_print("Write to file {}.".format(filepath), level=1)
    all_atoms = {}
    if selection_only:
        objects = bpy.context.selected_objects
    else:
        objects = bpy.context.scene.objects
    
    n_atoms = 0
    for ob in objects:
        if ob.mb.type == 'ATOM':
            try:
                all_atoms[ob.mb.molecule_name].append(
                    (ob.mb.index, ob.mb.element, get_world_coordinates(ob)))
            except KeyError:
                all_atoms[ob.mb.molecule_name] = [(ob.mb.index, ob.mb.element, 
                                                  get_world_coordinates(ob))]
            n_atoms += 1
    
    with open(filepath, "w") as fout:
        fout.write("{}\n".format(n_atoms))
        fout.write("This file has been created with Blender "
                   "and the MolBlend addon.\n")
    
        for mol_id in sorted(all_atoms):
            for i, element, location in sorted(all_atoms[mol_id]):
                l = "{}   {:>10.5f}   {:>10.5f}   {:>10.5f}\n"
                fout.write(l.format(element, *location))
        debug_print("Exported {} atoms.".format(n_atoms), level=1)
    return True


def export_pdb(filepath, selection_only, scale_distances):
    debug_print("Write to file {}.".format(filepath), level=1)
    all_atoms = {}
    if selection_only:
        objects = bpy.context.selected_objects
    else:
        objects = bpy.context.scene.objects
    
    n_atoms = 0
    for ob in objects:
        if ob.mb.type == 'ATOM':
            atom_name = ob.mb.atom_name
            info = (
                ob.mb.index,
                atom_name if len(ob.mb.element) != 1 else " " + atom_name,
                ob.mb.get_molecule().name[:3].upper(),
                get_world_coordinates(ob),
                ob.mb.element,
                )
            try:
                all_atoms[ob.mb.molecule_name].append(info)
            except KeyError:
                all_atoms[ob.mb.molecule_name] = [info]
            n_atoms += 1
        
    with open(filepath, "w") as fout:
        fout.write("REMARK This pdb file has been created with Blender "
                  "and the addon MolBlend\n"
                  "REMARK\n"
                  "REMARK\n")
        index = 1
        for mol_id in sorted(all_atoms):
            for data in sorted(all_atoms[mol_id]):
                i, name, res_name, coords, element = data
                l = ("HETATM{ID:>5} {name:<4} {res}     1    {c[0]:8.3f}" 
                     "{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00          "
                     "{element:>2}\n")
                fout.write(l.format(ID=index, name=name[:4], res=res_name, 
                                    c=coords, element=element))
                index += 1
                if index > 99999:
                    index = 1
    debug_print("Exported {} atoms.".format(n_atoms), level=1)
    return True


def export_b4w(filepath):
    # need to link all materials to meshes instead of objects
    # so first make each mesh single user
    # then link the material slots to data
    # export to b4w without autosaving blend file
    # undo all the changes and revert to state before (or alternatively copy
    # selection to new file and do everything there)
    debug_print("Blend4Web export not implemented yet.", level=1)
