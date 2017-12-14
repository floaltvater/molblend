# MolBlend

This addon for Blender (>2.79) primarily helps to import molecular structures
of (almost) any kind. In addition it adds some basic functionality to add to
imported structures, or even draw molecules from scratch.

## Import
### File formats

The currently supported file formats are PDB and XYZ, POSCAR (VASP), Abinit
and Quantum ESPRESSO input and output files. It is fairly easy to add other
file formats. (I decided to not got the openbabel route to avoid dependencies.)

In addition, once a structure has been imported, you can import a file
containing vibrations or phonon modes. Currently only the Quantum ESPRESSO
output format and a custom "XYZ" format is supported.

If you would like to see a file format implemented, please don't hesitate to
ask!

### Animation

For files that support several frames (e.g., PDB, XYZ, output from relaxation
calculations), the script imports all steps and animates the atom coordinates
1 step/frame.

After importing a vibration/mode file, you can select each individual mode to
be animated. Complex phonon modes and off-center q-points are also supported.

## Drawing

The script contains a simple modal operator, that lets you draw atoms with
left-clicks. Click-and-drawing from one atom will create a bond between this
to the new atom. Simple geometry constraints allow for nicer structures (drop-
down menu for angle constraints, and Ctrl for length constraint according to
covalent radius).

## Representations

Elements of the atoms can easily be changed, and are independent for each 
molecule (so one could have green, the other red carbon atoms).

You can switch easily between Ball-and-Stick, Stick, or Ball drawing, with
bonds colored generically or according to the bonded elements. Atom radii are
relative to the covalent or van der Waals radii, or just a constant value.

## Export

Exporting the structures to different file formats is planned for the future.
Please let me know if you have any preferences.

## Known issues

- Since every atom and every bond is its own object, large structures soon
  slow down the usage of Blender. Usability depends on the machine you are on.
  To work with huge structures, look into Bioblender (bioblender.org), or 
  other addons that use Dupliverts etc.
- Guessing the bonds is not very efficient or stable, since it compares every
  atom pair to the sum of their covalent radii.
- Deleting objects might lead to unintended consequences, since info attached
  to still existing objects is not necessarily updated. Please report any
  issues you run into!
- Unit handling is not very elegant. Currently the addon works exclusively in
  Angstrom (like guessing bonds), and converts, or tries to convert, imported
  files to that unit.
- Currently MolBlend doesn't support double/triple bonds, or other even less
  common bonds. But that is why you're working with Blender in the first
  place, right?
- Duplication of a molecule is like a "deep-duplication" (Alt+D) and most 
  properties (like atom colors etc.) are linked between the original and the
  duplicate.
  properties are still linked to each other.
