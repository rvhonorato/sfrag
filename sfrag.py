###############################################################################################################
#
#                               Bruna Inacio Trajano (brunainaciot@gmail.com)
#                                 26 Programa Bolsas de Verao (jan-feb 2017)
#                            Centro Nacional de Pesquisa em Energia e Materiais
#                                        Campinas - Sao Paulo - Brazil
#
###############################################################################################################
#
#                                      STEPS OF THE ALGORITHM - SPHERES VERSION
#
# 1. Gets the parameters: the files, the desired radius, chains and errors
# 2. Verifies if the args are ok
# 3. Opens the given files and gets only the common residues, that are stored in new-[nameoffile1].pdb and
# new-[nameoffile2].pdb
# 4. Gets the coordinates of these common atoms
# 5. Sets the convention that protein1 is the smaller one
# 6. Calculates the centroid of each protein
# 7. Translates each protein to its centroid and finds the optimal rotation between them
# 8. Creates the spheres
# 9. Calls the desired error function
# 10. Creates the output files
# 11. Centralizes and rotates the atoms of each sphere finding an optimal rotation
# 12. Calculates the errors
# 13. Writes the output files
#                                   STEPS OF THE ALGORITHM - FRAGMENTS VERSION
#
# 1. Gets the parameters: the files, the size of the fragment, the desired chainsand errors and if all
# the atoms will be used for the calculations
# 2-7. Idem previous 2-7.
# 8. Creates the fragments
# 9-10. Idem previous 9-10.
# 11. Centralizes and rotates the atoms of each fragment finding an optimal rotation
# 12-13. Idem previous 12-13.
#
# PS1.: in rmsd and canyon_metric, pr1, prot1 etc is the fragment or the sphere, not the entire protein
# PS2.: In the error functions, there are comments only in main_both
#
###############################################################################################################
from math import sqrt, pi, cos
import numpy as np
import argparse
import os
###############################################################################################################
#
#                                       INITIAL TREATMENT OF THE FILE
#
###############################################################################################################
# this is a test
# Calculates the centroid of the protein
def calc_centroid(protein):

    length = len(protein)
    sumxcentr = 0
    sumycentr = 0
    sumzcentr = 0

    for i in range(length):
        sumxcentr += protein[i][0]
        sumycentr += protein[i][1]
        sumzcentr += protein[i][2]

    xcentr = sumxcentr/length
    ycentr = sumycentr/length
    zcentr = sumzcentr/length

    return xcentr, ycentr, zcentr

# Translates the protein to the centroid previously calculated
def transl_to_centroid(prot, xcentr, ycentr, zcentr):

    length = len(prot)

    for i in range(length):
        prot[i][0] -= xcentr
        prot[i][1] -= ycentr
        prot[i][2] -= zcentr

    return prot

# Given a protein, this function returns a nx6 matrix with the CA data
def get_caprot(prot):

    caprot = []

    # caprot is a matrix nx6 (n = number of CAs)
    for i in range(len(prot)):
        if 'CA' in prot[i][5]:
            x = prot[i][0]
            y = prot[i][1]
            z = prot[i][2]
            numres = prot[i][3]
            nameres = prot[i][4]
            nameatom = prot[i][5]
            caprot.append([x, y, z, numres, nameres, nameatom])

    return caprot

# Translates and rotates the proteins
def translate_rotate(prot1, prot2):

    # calculates the centroid of each protein
    xcentr1, ycentr1, zcentr1 = calc_centroid(prot1)
    xcentr2, ycentr2, zcentr2 = calc_centroid(prot2)

    # translates each protein to its centroid
    prot1 = transl_to_centroid(prot1, xcentr1, ycentr1, zcentr1)
    prot2 = transl_to_centroid(prot2, xcentr2, ycentr2, zcentr2)

    # select the submatrices (with only coordinates) of prot1 and prot2
    mprot1 = [[prot1[i][j] for j in range(3)] for i in range(len(prot1))]
    mprot2 = [[prot2[i][j] for j in range(3)] for i in range(len(prot2))]

    # gets the rotation matrix
    mrotation = rotate(mprot1, mprot2)
    # gets a new prot1 (only this protein is rotated)
    rotprot1 = np.dot(mprot1, mrotation)

    # substitutes the coordinates of the old prot1 for the new ones
    for i in range(len(prot1)):
        prot1[i][0] = rotprot1[i][0]
        prot1[i][1] = rotprot1[i][1]
        prot1[i][2] = rotprot1[i][2]

    # gets the CA coordinates from the proteins already translated and rotated
    caprot1 = get_caprot(prot1)
    caprot2 = get_caprot(prot2)

    return prot1, prot2, caprot1, caprot2


# If the two given proteins are different, this function gets the common residues between them
def equal_protein(new_file1, new_file2, f1, f2, chain1, chain2):

    # lists with all the lines of the files
    lines_f1 = f1.readlines()
    lines_f2 = f2.readlines()

    # lists of residues of proteins 1 and 2
    list_of_res1 = []
    list_of_res2 = []

    # gets the residues of protein 1 that are in the desired chain
    for l1 in lines_f1:
        if l1[:4] == 'ATOM':
            if l1[21] == chain1:
                # appends the atom number, the residue name and its number, respectively
                list_of_res1.append((l1[7:11], l1[17:20], l1[23:26]))

    # eliminates the multiple entries
    list_of_res1 = list(set(list_of_res1))

    for l2 in lines_f2:
        if 'ATOM' in l2[:4]:
            if l2[21] == chain2:
                list_of_res2.append((l2[7:11], l2[17:20], l2[23:26]))
    list_of_res2 = list(set(list_of_res2))

    # gets the intersection between the groups of residues
    intersec_res = [line for line in list_of_res1 if line in list_of_res2]

    if not intersec_res:
        print "The chains intersection is empty! Are you sure about them?"
        exit()

    # writes two pdb files only with the common residues (they'll be the new 'input')
    for line in lines_f1:
        if (line[7:11], line[17:20], line[23:26]) in intersec_res:
            if 'ATOM' in line[:4] and line[21] == chain1:
                if 'H' not in line[77:78]:
                    new_file1.write(line)

    for line in lines_f2:
        if (line[7:11], line[17:20], line[23:26]) in intersec_res:
            if 'ATOM' in line[:4] and line[21] == chain2:
                if 'H' not in line[77:78]:
                    new_file2.write(line)

    return new_file1, new_file2


###############################################################################################################
#
#                                               METRIC FUNCTIONS
#
###############################################################################################################

# RMSD
def rmsd(pr1, pr2, matrot):

    sumdistances = 0
    #rotates the protein
    newprot1 = np.dot(pr1, matrot)

    #calculates rmsd
    for i in range(len(pr1)):
        sumdistances += dist(newprot1[i], pr2[i])**2

    rmsdresult = sqrt(sumdistances/len(pr1))
    return rmsdresult


# My metric
def canyon_metric(pr1, pr2, mrot, limitdist, errorct):

    #acumulated error
    errorac = 0
    #rotates the protein
    newprot1 = np.dot(pr1, mrot)

    for i in range(len(pr1)):
        #distance between two atoms
        dist = sqrt((newprot1[i][0]-pr2[i][0])**2 + (newprot1[i][1]-pr2[i][1])**2 + (newprot1[i][2]-pr2[i][2])**2)
        if (dist>=limitdist):
            errorac += errorct
        else:
            # my formula
            errorcos = -(errorct/2)*cos(pi*dist/limitdist) + errorct/2
            errorac = errorac + errorcos
    #normalizes the result - the maximum error is always errorct
    return errorac/len(pr1)

###############################################################################################################
#
#                                   OTHER NECESSARY FUNCTIONS
#
###############################################################################################################

# Verifies if all the arguments given are ok
def verify_args(arg_rmsd, can, both, sph, frag, atoms, size=None, rad=None):

    # if more than one type of error is selected
    if (can==True and both==True) or (can==True and arg_rmsd==True) or (both==True and arg_rmsd==True):
        print "Choose only ONE type of error!"
        exit()

    # if no one type of error is selected
    if can == False and both == False and arg_rmsd== False:
        print "Choose a type of error!"
        exit()

    # if no one type of analysis is selected
    if sph == False and frag == False:
        print "Choose a type of analysis (sphere or fragments)!"
        exit()

    # if both types of analysis are selected
    if sph == True and frag == True:
        print "Choose only ONE type of analysis!"
        exit()

    # if no one (size or radius) is given
    if size == None and rad == None:
        print "Choose a size or a radius!"
        exit()

    # if both (radius and size) are given
    if size != None and rad != None:
        print "Choose radius or size!"
        exit()

    # the spheres metric already includes all the atoms,
    # so the parameter -a is unnecessary
    if atoms == True and sph == True:
        print "The spheres metric does not accept the parameter -a!"
        exit()

    # wrong combination of parameters
    if sph == True and size != None:
        print "This is not the right parameter. You must set a radius (parameter -r)."
        exit()

    # wrong combination of parameters
    if frag == True and rad != None:
        print "This is not the right parameter. You must set a SIZE (parameter -s)."
        exit()

    # if the given size has an invalid value
    if size<=0 and rad == None:
        print "Size must be positive!"
        exit()

    # if the given radius has an invalid value
    if rad<=0 and size == None:
        print "Radius must be positive!"
        exit()

    # if the size given is odd (this is necessary because, in fragments
    # metric, the fragment will be centralized by the middle CA)
    if size != None:
        if size%2 == 0:
            print "You must set an odd number as size!"
            exit()


# Gets the coordinates of all atoms (except for H) from the files
# (the function also returns a list with the coordinates of CAs for future purposes)
def get_coord(arccoord, chain):

    coordall = []

    # goes through the lines of the file
    # coordall will be a nx6 matrix (n = number of atoms)
    for line in arccoord:
        x = float(line[31:38])
        y = float(line[39:46])
        z = float(line[47:54])
        numres = int(line[23:26])
        nameres = line[17:20]
        nameatom = line[13:16]
        coordall.append([x, y, z, numres, nameres, nameatom])

    return coordall

# Finds the optimal rotate matrix using SVD
def rotate(R1, R2):

    if len(R1) != len(R2):
        print "The proteins have different sizes! It's impossible to make the calculations!"
        exit()

    A = np.dot(np.transpose(R1), R2)
    V, S, W = np.linalg.svd(A)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    # if necessary, changes the system of coordinates to a right-handed one
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # create rotation matrix U
    U = np.dot(V, W)

    return U


# Calculates the distance between two atoms
def dist(atom1, atom2):
    return sqrt((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + (atom1[2]-atom2[2])**2)

# Gets the atoms in each sphere
def atoms_in_sphere(prot1, prot2, caprot, rad):

    atoms_sphere1 = []
    atoms_sphere2 = []
    union = []
    coords_atoms1 = []
    coords_atoms2 = []

    # atoms_sphere1, for example, is a list of lists of tuples - each tuple
    # is an atom (its residue number, residue name and name), and each list
    # of tuples is a sphere
    for i in range(len(caprot)):
        atoms_sphere1.append([])
        atoms_sphere2.append([])
        for j in range(len(prot1)):
            # checks which atoms in prot1 are inside the sphere centered at CA (of prot1)
            if dist(caprot[i], prot1[j][:3]) < rad:
                atoms_sphere1[i].append(tuple(prot1[j][3:]))
            # checks which atoms in prot2 are inside the same sphere
            if dist(caprot[i], prot2[j][:3]) < rad:
                atoms_sphere2[i].append(tuple(prot2[j][3:]))

    # gets the union between both sets of atoms, eliminating duplicates
    # and storing it in all_atoms
    for i in range(len(atoms_sphere1)):
        union.append(set(atoms_sphere1[i] + atoms_sphere2[i]))

    # gets the coordinates of the atoms in each sphere
    # (s = abbreviation for sphere)
    # coords_atoms1 is a list of lists (spheres) of lists (atoms)
    for s in range(len(union)):
        coords_atoms1.append([])
        coords_atoms2.append([])
        for i in range(len(prot1)):
            if tuple(prot1[i][3:]) in union[s]:
                coords_atoms1[s].append(prot1[i][:3])
            if tuple(prot2[i][3:]) in union[s]:
                coords_atoms2[s].append(prot2[i][:3])

    return coords_atoms1, coords_atoms2

# If 'fragments' is chosen, centralizes the fragment putting the middle ca in (0,0,0) 
# If 'spheres' is chosen, centralizes the atoms putting the CA sphere's center in (0,0,0)
def center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, coords_atoms1, coords_atoms2, pack_atm_prot1, pack_atm_prot2, all_atoms):

    if sph == True:
        # gives the sphere
        p1c = center(coords_atoms1[i], caprot1[i], True, sph)
        p2c = center(coords_atoms2[i], caprot2[i], True, sph)
    elif all_atoms == True:
        # gives the set of aminoacids
        p1c = center(pack_atm_prot1[i-size/2:i+size/2+1], caprot1[i], all_atoms, sph)
        p2c = center(pack_atm_prot2[i-size/2:i+size/2+1], caprot2[i], all_atoms, sph)
    elif all_atoms == False:
        # gives the CAs
        p1c = center(prot1[i-size/2:i+size/2+1], [], all_atoms, sph)
        p2c = center(prot2[i-size/2:i+size/2+1], [], all_atoms, sph)

    mrot = rotate(p1c, p2c)
        
    return p1c, p2c, mrot

# Translates the set of atoms
def center(atm, coord_ca, all_atoms, sph):

    v = []

    if all_atoms == True:
        for l1 in atm:
            if sph == True:
                v.append([l1[0] - coord_ca[0], l1[1] - coord_ca[1], l1[2] - coord_ca[2]])
            else:
                for l2 in l1:
                    v.append([l2[0] - coord_ca[0], l2[1] - coord_ca[1], l2[2] - coord_ca[2]])
    else:
        x = atm[len(atm)/2][0]
        y = atm[len(atm)/2][1]
        z = atm[len(atm)/2][2]

        for atom in range(len(atm)):
            v.append([atm[atom][0] - x, atm[atom][1] -y, atm[atom][2] -z])

    return v

# Creates and writes the output files 
def output(original_file, protein, caprot, answer, m, size, sph):

    reslist = []
    lines = original_file.readlines()

    # creates a list of unique residues
    for l in lines:
        residue = int(l[23:26])
        if residue not in reslist:
            reslist.append(residue)

    # creates a dictionary where the errors will be stored according to residue
    dictm = {}

    # when it's sphere metric, all the residues have an associated error
    # however, in the fragments metric, the first and last ones doesn't have it 
    # because the first/last CAs aren't center of no one fragment
    if sph == True:
        for i in range(len(reslist)):
            dictm[int(reslist[i])] = m[i]
    else:
        for i in range(len(reslist)):
            if i<size/2 or i>=len(reslist)-size/2:
                dictm[int(reslist[i])] = 0
            else:
                dictm[int(reslist[i])] = m[i-size/2]

    # writes the file
    natom = 0
    for l in lines:
        if natom != len(protein):
            x = "{:8.3f}".format(protein[natom][0])
            y = "{:8.3f}".format(protein[natom][1])
            z = "{:8.3f}".format(protein[natom][2])
            error = "{:6.2f}".format(dictm[int(l[23:26])])
            textm = l[:30] + x + y + z + l[54:60] + error + l[66:]
            answer.write(str(textm))
            natom +=1

# Creates two matrices, pack_atm_prot1 and pack_atm_prot2
# Each one of them is a list (protein) of lists (residues) of lists
# (coordinates of the atoms)
def pack_atm(prot1, prot2, rlist):
    
    pack_atm_prot1 = []
    pack_atm_prot2 = []

    for j in range(len(rlist)):
        atoms_res_prot1 = []
        atoms_res_prot2 = []
        for i in range(len(prot1)):
            if prot1[i][3] == rlist[j]:
                atoms_res_prot1.append(prot1[i][:3])
                atoms_res_prot2.append(prot2[i][:3])
        pack_atm_prot1.append(atoms_res_prot1)
        pack_atm_prot2.append(atoms_res_prot2)

    return pack_atm_prot1, pack_atm_prot2


###############################################################################################################
#
#                                   MAIN FUNCTIONS
#
# main_canyon: called when the canyon metric must be calculated
# main_rmsd: called when rmsd must be calculated
# main_both: called when the both metrics must be calculated
#
# PS.: the comments are in main_both()
###############################################################################################################

def main_canyon(new_file1, new_file2, prot1, prot2, caprot1, caprot2, dist, error, sph, coords_atoms1, coords_atoms2, size, atoms):

    c = []

    c1 = new_file1.readline()
    c2 = new_file2.readline()

    if sph == True:
        analysis = '_sph_'
    else:
        analysis = '_frag'

    if atoms == True:
        all_atoms = '_all_'
    elif sph == False:
        all_atoms = '_ca_'
    else:
        all_atoms = ''


    out1 = open('can' + analysis + all_atoms + c1[21] + '_' + new_file1.name[4:], 'w')
    out2 = open('can' + analysis + all_atoms + c2[21] + '_' + new_file2.name[4:], 'w')

    if atoms == True:
        rlist = []

        for i in range(len(prot1)):
            residue = int(prot1[i][3])
            if residue not in rlist: 
                rlist.append(residue)

        pack_atm_prot1, pack_atm_prot2 = pack_atm(prot1, prot2, rlist)

    if sph == True:
        for i in range(len(coords_atoms1)):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, 0, coords_atoms1, coords_atoms2, [], [], True)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))
    elif atoms == True:
        for i in range(size/2, len(caprot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], pack_atm_prot1, pack_atm_prot2, True)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))
    elif atoms == False:
        for i in range(size/2, len(prot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], [], [], False)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))

    new_file1.seek(0)
    new_file2.seek(0)

    output(new_file1, prot1, caprot1, out1, c, size, sph)
    output(new_file2, prot2, caprot2, out2, c, size, sph)

    out1.close()
    out2.close()

def main_rmsd(new_file1, new_file2, prot1, prot2, caprot1, caprot2, sph, coords_atoms1, coords_atoms2, size, atoms):

    r = []

    c1 = new_file1.readline()
    c2 = new_file2.readline()

    if sph == True:
        analysis = '_sph_'
    else:
        analysis = '_frag'

    if atoms == True:
        all_atoms = '_all_'
    elif sph == False:
        all_atoms = '_ca_'
    else:
        all_atoms = ''

    out1 = open('rmsd' + analysis + all_atoms + c1[21] + '_' + new_file1.name[4:], 'w')
    out2 = open('rmsd' + analysis + all_atoms + c2[21] + '_' + new_file2.name[4:], 'w')

    rlist = []

    if atoms == True:

        for i in range(len(prot1)):
            residue = int(prot1[i][3])
            if residue not in rlist: 
                rlist.append(residue)

        pack_atm_prot1, pack_atm_prot2 = pack_atm(prot1, prot2, rlist)

    if sph == True:
        for i in range(len(coords_atoms1)):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, 0, coords_atoms1, coords_atoms2, [], [], True)
            r.append(rmsd(p1c, p2c, mrot))
    elif atoms == True:
        for i in range(size/2, len(caprot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], pack_atm_prot1, pack_atm_prot2, True)
            r.append(rmsd(p1c, p2c, mrot))
    elif atoms == False:
        for i in range(size/2, len(prot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], [], [], False)
            r.append(rmsd(p1c, p2c, mrot))

    new_file1.seek(0)
    new_file2.seek(0)
    
    output(new_file1, prot1, caprot1, out1, r, size, sph)
    output(new_file2, prot2, caprot2, out2, r, size, sph)
    out1.close()
    out2.close()

def main_both(new_file1, new_file2, prot1, prot2, caprot1, caprot2, dist, error, sph, coords_atoms1, coords_atoms2, size, atoms):
    
    # creates the lists where the errors will be stored
    c = [] # canyon
    r = [] # rmsd

    # some variables that will be in the output file's name
    c1 = new_file1.readline()
    c2 = new_file2.readline()

    if sph == True:
        analysis = '_sph_'
    else:
        analysis = '_frag'

    if atoms == True:
        all_atoms = '_all_'
    elif sph == False:
        all_atoms = '_ca_'
    else:
        all_atoms = ''

    # creates, and opens, the files where the errors will be written
    out1 = open('can' + analysis + all_atoms + c1[21] + '_' + new_file1.name[4:], 'w')
    out2 = open('can' + analysis + all_atoms + c2[21] + '_' + new_file2.name[4:], 'w')
    out3 = open('rmsd' + analysis + all_atoms + c1[21] + '_' + new_file1.name[4:], 'w')
    out4 = open('rmsd' + analysis + all_atoms + c2[21] + '_' + new_file2.name[4:], 'w')

    if atoms == True:
        rlist = []

        # creates a list with the residues - just like reslist in output()
        for i in range(len(prot1)):
            residue = int(prot1[i][3])
            if residue not in rlist:
                rlist.append(residue)

        # see pack_atm explanation
        pack_atm_prot1, pack_atm_prot2 = pack_atm(prot1, prot2, rlist)

    # centralizes and rotates the set of atoms, also calculating the desired errors
    if sph == True:
        for i in range(len(coords_atoms1)):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, 0, coords_atoms1, coords_atoms2, [], [], True)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))
            r.append(rmsd(p1c, p2c, mrot))
    elif atoms == True:
        for i in range(size/2, len(caprot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], pack_atm_prot1, pack_atm_prot2, True)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))
            r.append(rmsd(p1c, p2c, mrot))
    elif atoms == False:
        for i in range(size/2, len(prot1)-size/2):
            p1c, p2c, mrot = center_and_rotate(caprot1, caprot2, prot1, prot2, sph, i, size, [], [], [], [], False)
            c.append(canyon_metric(p1c, p2c, mrot, dist, error))
            r.append(rmsd(p1c, p2c, mrot))

    new_file1.seek(0)
    new_file2.seek(0)

    # writes in output the files
    output(new_file1, prot1, caprot1, out1, c, size, sph)
    output(new_file2, prot2, caprot2, out2, c, size, sph)

    new_file1.seek(0)
    new_file2.seek(0)

    output(new_file1, prot1, caprot1, out3, r, size, sph)
    output(new_file2, prot2, caprot2, out4, r, size, sph)

    # closes all the files
    out1.close()
    out2.close()
    out3.close()
    out4.close()

def main():

    # gets the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("f1", help="Choose the first file pdb")
    parser.add_argument("f2", help="Choose the second file pdb")
    parser.add_argument("c1", help="Set the wanted chain of the first protein")
    parser.add_argument("c2", help="Set the wanted chain of the second protein")
    parser.add_argument("-sph", help="Set the type of the analysis: by the spheres method", action="store_true")
    parser.add_argument("-frag", help="Set the type of the analysis: by the fragments method", action="store_true")
    parser.add_argument("-r", help="Set the radius for the sphere of residues (in angstrons)", type=float)
    parser.add_argument("-s", help="Set the size for the fragments (in number of residues)", type=int)
    parser.add_argument("-rmsd", help="Calculate only the RMSD",action="store_true")
    parser.add_argument("-can", help="Calculate only the Canyon metric", action="store_true")
    parser.add_argument("-both", help="Calculate the RMSD and the Canyon metric", action="store_true")
    parser.add_argument("-d", help="Set a cut-off distance for the Canyon metric (default=2)", action="store", default=2)
    parser.add_argument("-e", help="Set the constant error for the Canyon metric (default=4)", action="store", default=4)
    parser.add_argument("-a", help="Calculate the selected error for all the atoms, not just CA", default=False, action="store_true")

    parser.parse_args()
    args = vars(parser.parse_args())

    file1 = args['f1']
    file2 = args['f2']
    chain1 = args['c1']
    chain2 = args['c2']
    sph = args['sph']
    frag = args['frag']
    rad = args['r']
    size = args['s']
    arg_rmsd = args['rmsd']
    can = args['can']
    both = args['both']
    dist = float(args['d'])
    error = float(args['e'])
    atoms = args['a']

    # verifies if the args are ok
    verify_args(arg_rmsd, can, both, sph, frag, atoms, size, rad)

    # opens the given files 1 and 2, and creates new-file1 and 2
    # (where the common atoms will be stored)
    f1 = open(file1, "r")
    f2 = open(file2, "r")
    new_file1 = open('new-'+ file1, 'w+')
    new_file2 = open('new-' + file2, 'w+')

    new_file1, new_file2 = equal_protein(new_file1, new_file2, f1, f2, chain1, chain2)

    # takes the file's pointer to line 0
    new_file1.seek(0)
    new_file2.seek(0)

    # gets the coordinates from the new files
    prot1 = get_coord(new_file1, chain1)
    prot2 = get_coord(new_file2, chain2)

    new_file1.seek(0)
    new_file2.seek(0)

    # uses always the smaller protein (1) as reference
    if len(prot1)>len(prot2):
        prot3 = prot1
        prot1 = prot2
        prot2 = prot3
    
    # translates and rotates the proteins
    prot1, prot2, caprot1, caprot2 = translate_rotate(prot1, prot2)

    if sph == True:
        coords_atoms1, coords_atoms2 = atoms_in_sphere(prot1, prot2, caprot1, rad)
    else:
        coords_atoms1 = []
        coords_atoms2 = []

    # calls the desired error function
    if can == True:
        main_canyon(new_file1, new_file2, prot1, prot2, caprot1, caprot2, dist, error, sph, coords_atoms1, coords_atoms2, size, atoms)

    if both == True:
        main_both(new_file1, new_file2, prot1, prot2, caprot1, caprot2, dist, error, sph, coords_atoms1, coords_atoms2, size, atoms)

    if arg_rmsd == True:
        main_rmsd(new_file1, new_file2, prot1, prot2, caprot1, caprot2, sph, coords_atoms1, coords_atoms2, size, atoms)

    # closes the files and deletes the ones that aren't necessary anymore
    f1.close()
    f2.close()
    os.remove(new_file1.name)
    os.remove(new_file2.name)

if __name__ == '__main__':
    main()
