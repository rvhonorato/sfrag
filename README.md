** What is SFrag? **
SFrag - Spheres and Fragments Protein Comparator is an algorithm that compares two conformations of the same protein. It divides the conformations in fragments (or spheres, depending on the user's choice) and compares them fragment by fragment (or sphere by sphere).

** Parameters **
The algorithm has the following parameters:
- f1, f2: names of the desired PDB files for comparation;
- c1, c2: desired chains of f1 and f2, respectively;
- sph: selects the spheres method;
- frag: selects the fragments method;
- r: sets the spheres radius (in angstrons);
- s: sets the fragments size (in number of residues) (obs: it must be an odd number);
- rmsd: selects only the RMSD metric;
- can: selects only the Canyon metric (explained below);
- both: selects both metrics;
- d: optional parameter that sets a limit distance in the Canyon metric;
- e: optional parameter that sets a constant error in the Canyon metric;
- a: in fragments method, considers all the atoms in the analysis (the default is to consider only CAs).

Only one method (fragments or spheres) is allowed each time that the algorithm is executed.
Usage example: *1b39.pdb 1fq1.pdb A B -frag -rmsd -s=7*

** RMSD **
See https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions.

** Canyon metric **
Canyon is a metric developed for this project. Its objective is to penalize with a constant error distances bigger than a limit distance. Distances smaller than it are penalized with an error related to the limit distance, calculated distance and constant error (the formula is in the code).
For each pair of atoms, the error is calculated and added to a variable. At the end of the fragment/sphere, the value is normalized by the number of pairs compared. So the maximum value is the constant error.

** Fragments **
If requested, the algorithm divides the proteins in fragments. For example, if the amino acid sequence is *1, 2, 3, 4, 5, 6, 8...* and the user sets s = 5, the algorithm will consider:
- fragment 1 = 1, 2, 3, 4, 5
- fragment 2 = 2, 3, 4, 5, 6
- fragment 3 = 3, 4, 5, 6, 8 etc
For each fragment, a translation and a rotation are done to verify (by RMSD and/or Canyon) how much they match. So every fragment has a value that indicates its similarity, and this value is associated with the middle residue of the fragment (3, 4 and 5, in the example).

** Spheres **
If requested, the algorithm divides the proteins in spheres. Each CA of each protein is the center of a sphere with radius set by the user. Then, the union of the atoms of both spheres is made and, for each sphere, RMSD and/or Canyon is calculated.

** Output **
The algorithm gives as output PDB files with the calculated errors replacing the temperature factor.

** The steps of the algorithm are commented on it. ** 

** Contact **
Bruna Inacio Trajano (brunainaciot@gmail.com)
