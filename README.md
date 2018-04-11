# MolBlend

This addon for Blender (>2.79) primarily helps to import molecular structures
of (almost) any kind. In addition it adds some basic functionality to add to
imported structures, or even draw molecules from scratch.

## Import
### File formats

The currently supported file formats are PDB, XYZ, POSCAR (VASP), Abinit
output, Quantum ESPRESSO input and output files, and Gaussian cube files.

In addition, once a structure has been imported, you can import a file
containing vibrations or phonon modes. Currently Quantum ESPRESSO, Abinit
(anaddb), phonopy and a custom "XYZ" format are supported.

MolBlend can also read the volume data in cube files to create isosurfaces.
However, to use this feature, you have to install pymcubes on your computer
and change the path to its library in mb_import_export.py. Make sure that 
the python version for which you install pymcubes is the same as the one
Blender uses (for Blender 2.79 it is python 3.5).

I decided to implement all file imports by myself to avoid having to install
any external dependencies.
If you would like to see a file format implemented, please don't hesitate to
ask!

### Animation

For files that support several frames (e.g., PDB, XYZ, output from relaxation
calculations), the script imports all steps and animates the atom coordinates
1 step/frame.
This is obviously meant for analysis, not to make a pretty movie. For a movie
with a more realistic frame rate, you could adjust the time remapping value in
the Render properties panel.

After importing a vibration/mode file, you can select each individual mode to
be animated. Complex valued phonon modes and off-center q-points are also 
supported.

## Drawing

The script contains a simple modal operator, that lets you draw atoms with
left-clicks. Click-and-drag when hovering over another atom will create a bond
between that atom and the new atom. 
Simple geometry constraints allow for nicer structure drawing. There is a 
drop-down menu for angle constraints, and holding Ctrl while drawing a bonded
atom enforces a length constraint according to the covalent radii.

## Representations

Representations of elements can be changed independently for each molecule. 
So one could have small green carbon atoms in one molecule, and large red 
carbon atoms in another.

You can switch between Ball-and-Stick, Stick, or Ball drawing, with
bonds colored generically or according to the bonded elements. Atom radii are
relative to the (adjustable) covalent or van der Waals radii, or just a 
constant value.

Molecules are not necessarily bonded. All atoms imported from a single file
initially belong to a single atom, no matter if they are all connected or not.
One can split a group of atoms into separate molecules manually, again whether
they are connected or not.

## Export

Exporting the structures to different file formats is planned for the future.
Please let me know if you have any preferences.

## Known issues

- Since every atom and every bond is its own object, large structures soon
  slow down the usage of Blender. Usability depends on the machine you are on.
  To work with huge structures, look into Bioblender (bioblender.org), or 
  other addons that use Dupliverts etc.
- Guessing the bonds is not very efficient or stable, since it compares every
  atom pair to the sum of their covalent radii plus some tolerance.
- Unit handling is not very elegant. Currently the addon works exclusively in
  Angstrom (like guessing bonds), and converts, or tries to convert, imported
  files to that unit.
- Currently MolBlend doesn't support double/triple bonds, or other even less
  common bonds.
- Duplication (Shift+D, Alt+D) of a molecule has various issues, since the 
  underlying MB_Molecule Property is not updated (doesn't know of the new
  atoms, atom colors, draw styles are linked between old and new atoms, etc.).
  A custom MolBlend duplicate operator might be in the future.
- MolBlend doesn't necessarily preserve indices that are explicitly written 
  in input files (like PDB files). It does preserve the order of atoms though.
- Dipoles or unit cell axes might end up having zero dimension due to the
  stretch constraint. Manually switching the "Plane" in the corresponding 
  mb.stretch constraint should fix that. For some reason this doesn't work
  reliably from the python API.
- Logging is somewhat sparse and not functional for debugging.
