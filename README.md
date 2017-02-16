## What is SFrag?
**SFrag** - *Spheres and Fragments Protein Comparator* is an algorithm that compares two conformations of the same protein. It divides the conformations in **fragments** *or* **spheres** and compares them fragment by fragment or sphere by sphere.

This software was developed by **Bruna Inácio Trajano** during the [26th Summer Internship Program](http://pages.cnpem.br/bolsasdeverao/) at the [National Center of Research in Energy and Materials - CNPEM](http://cnpem.br) and was advised by [Rodrigo Honorato](https://github.com/rodrigovrgs), PhD and [José Geraldo de Carvalho](https://github.com/jgcarvalho), MsC.

## Requirements

- Python 2.7
- NumPy

## Usage examples

1) Compare using both metrics and sphere method with 5Å radius

`python sfrag.py proteinA.pdb proteinB.pdb A B -sph -r 5 -both`

2) Compare using Canyon metric and fragment method

`python sfrag.py proteinA.pdb proteinB.pdb A B -frag -s 3 -can`

3) Compare using RMSD metric and sphere method with 7Å radius

`python sfrag.py proteinA.pdb proteinB.pdb A B -sph -r 7 -rmsd`

## Parameters

**Note**: Methods are mutually exclusive

#### Input
- *-f1*, *-f2*
    - names of the desired PDB files for comparation;
- *-c1*, *-c2*
    - desired chains of f1 and f2, respectively;

#### Methods
- *-sph*
    - selects the spheres method;
- *-frag*
    - selects the fragments method;

#### Specific parameters
- *-r*
    - sets the spheres radius (in angstrons);
- *-s*
    - sets the fragments size (in number of residues) (obs: it must be an odd number);

#### Metrics
- *-rmsd*
    - selects only the RMSD metric;
- *-can*
    - selects only the Canyon metric (explained below);
- *-both*
    - selects both metrics;

#### Optional
- *-d*
    - optional parameter that sets a limit distance in the Canyon metric;
- *-e*
    - optional parameter that sets a constant error in the Canyon metric;
- *-a*
    - in fragments method, considers all the atoms in the analysis (the default is to consider only CAs).

## Metrics
### Canyon
Canyon is a metric developed for this project. Its objective is to penalize with a constant error distances bigger than a limit distance. Distances smaller than it are penalized with an error related to the limit distance, calculated distance and constant error (the formula is in the code).

For each pair of atoms, the error is calculated and added to a variable. At the end of the fragment/sphere, the value is normalized by the number of pairs compared. So the maximum value is the constant error.


![](https://www.ime.usp.br/~rvargas/sfrag_canyon_form.gif)


![](https://www.ime.usp.br/~rvargas/sfrag_canyon.png)

### RMSD
Check [Wikipedia](https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions) for more info.


## Fragments
If requested, the algorithm divides the proteins in fragments `-frag`. For example, if the amino acid sequence is *1, 2, 3, 4, 5, 6, 8...* and the user sets `s 5`, the algorithm will consider:
- frag 1 = 1, 2, 3, 4, 5
- frag 2 = 2, 3, 4, 5, 6
- frag 3 = 3, 4, 5, 6, 8
- frag 4 = 4, 5, 6, 7, 9
- ...



For each fragment, a translation and a rotation are done to verify (by RMSD and/or Canyon) how much they match. So every fragment has a value that indicates its similarity, and this value is associated with the middle residue of the fragment (3, 4, 5 and 6, in the example).

## Spheres
If requested, the algorithm divides the proteins in spheres. Each CA of each protein is the center of a sphere with radius set by the user. Then, the union of the atoms of both spheres is made and, for each sphere, RMSD and/or Canyon is calculated.

## Output
The algorithm gives as output PDB files with the calculated errors replacing the temperature factor.

PS.: The steps of the algorithm are commented on it.

## Problems? Suggestions? Let me know!
Send an e-mail to Bruna Inacio Trajano (brunainaciot@gmail.com).


___
*SFrag is an open source tool distributed under GPLv3*
