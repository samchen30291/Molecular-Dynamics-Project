import sys
import math
import itertools

mol = []
moldata = []
molcoord = []
bonds = []
dict_parameter = {
                 "CC":1.53,"CH":1.11,
                 "HCH":109.5,"HCC":109.5,"HCC":109.5,"CCC":109.5,
                 "XCCX":3,
                 "sigma_H":1.20,"sigma_C":1.75,
                 "kb_CC":300,"kb_CH":350,
                 "ka_HCH":35,"ka_HCC":35,"ka_CCC":60,
                 "Aphi":0.3,
                 "epsilon_H":0.03,"epsilon_C":0.07
                 }

def readfile():
    global n_atom
    global n_bond
    global n_carbon
    global n_cc
    global n_inputfield
    with open(sys.argv[1]) as file:
    
        mol = file.readlines() # the return is in ["line 1","line 2","line 3","line 4",....]
        
    ## first block in the reading file
        firstline = mol[0].split() # delete the space
        for f in firstline:
            moldata.append(f) # input the first line, get to know the basic mol data
        n_atom = int(moldata[0])
        n_bond = int(moldata[1])
        n_carbon = int(moldata[2])
        n_cc = int(moldata[3])
        n_inputfield = list(moldata[4:]) # the unuse force field
        
    ## second block in the reading file
        atomcoord = mol[1:1+n_atom] #read the coord line
        for lines in atomcoord:
            word = lines.split() # words in every line
            coords = list(map(float,word[0:3])) 
            atom = word[3]
            forcefield = list(map(int, word[4:]))
            atomcoord = {"coord":coords,"atom":atom,"forcefield":forcefield}
            molcoord.append(atomcoord)# add dict into list, list can sequence the atoms
    
    ## third block in the reading file
        restline = mol[1+n_atom:1+n_atom+n_bond]
        for lines in restline:
            words = lines.split()
            bondinfo = list(map(int,words[0:3]))
            forcefield = list(map(int,words[3:]))
            bond = {"bondinfo":bondinfo,"forcefield":forcefield}
            bonds.append(bond)
def distance(a,b):
    a1,a2,a3 = a
    b1,b2,b3 = b
    return math.sqrt((a1 - b1)**2 + (a2 - b2)**2 + (a3 - b3)**2)
def bondlength():
    global dist_info
    dist_info = []
    for i in range(n_bond):
        Atom_a = bonds[i].get("bondinfo")[0]
        Atom_b = bonds[i].get("bondinfo")[1]
        dist = distance(molcoord[Atom_a-1].get("coord"),molcoord[Atom_b-1].get("coord"))
        dist_info.append(dist)

def compute_angle(a, b, c):
    ax,ay,az = molcoord[a-1].get("coord")
    bx,by,bz = molcoord[b-1].get("coord")
    cx,cy,cz = molcoord[c-1].get("coord")     
    v1 = (ax - bx, ay - by, az - bz)
    v2 = (cx - bx, cy - by, cz - bz)
    dot = sum(v1[i] * v2[i] for i in range(3))
    norm1 = math.sqrt(sum(v1[i]**2 for i in range(3)))
    norm2 = math.sqrt(sum(v2[i]**2 for i in range(3)))  
    cos_theta = dot / (norm1 * norm2)
    return math.acos(cos_theta) #output of acos in radians -> convert to degree

def angle():
    global angles
    global neighbors
    global degree_angles
    angles = []
    degree_angles = []
    bond_profile = []
    for i in range(n_bond):
        a1,a2,a3 = bonds[i].get("bondinfo",0)
        bond_profile.append(tuple([a1,a2])) # make profile for all bonds
    neighbors = {i: set() for i in range(1, n_atom+1)} # build a dict for every center(key) : {set of neighbors}
    for a1, a2 in bond_profile:
        neighbors[a1].add(a2)
        neighbors[a2].add(a1) # to add neighbors mutually, ouput {1: {2, 3, 4, 5}, 2: {1}, 3: {1}, 4: {1}, 5: {1}}
    for center, nbs in neighbors.items():
        if len(nbs) < 2: # nbs >=2 to make an angle
            continue
        for a, c in itertools.combinations(nbs,2): # make a combination in nbs C([.....],2)
            angle = compute_angle(a,center,c)
            angles.append((a,center,c,angle))
    for angle in angles:
        degree_angles.append(math.degrees(angle[3]))
        
def compute_torsion(a,b,c,d):
    r1 = molcoord[a-1].get("coord")
    r2 = molcoord[b-1].get("coord")
    r3 = molcoord[c-1].get("coord")
    r4 = molcoord[d-1].get("coord")
    def vec_sub(u,v):
        return (u[0]-v[0],u[1]-v[1],u[2]-v[2])
    def dot(u,v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u,v):
        return (
            u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0],
                )
    b1 = vec_sub(r2,r1)
    b2 = vec_sub(r3,r2)
    b3 = vec_sub(r4,r3)
    n1 = cross(b1,b2)
    n2 = cross(b2,b3)
    n1_norm = math.sqrt(sum(n1[i]**2 for i in range(3)))
    n2_norm = math.sqrt(sum(n2[i]**2 for i in range(3)))
    if n1_norm == 0 or n2_norm == 0:
        raise ValueError("Atoms are collinear; torsion undefined")
    cos_theta = dot(n1,n2) / (n1_norm*n2_norm)
    cos_theta = max(-1,min(cos_theta,1))
    return math.acos(cos_theta) # return radians

def torsion():
    global torsions
    global torsions_angles
    global torsions_degree_angles
    torsions = []
    torsions_angles = []
    bond_profile = []
    torsions_degree_angles = []

    for info in bonds:
        B = info.get("bondinfo")[0]
        C = info.get("bondinfo")[1]
        
        if (molcoord[B-1].get("atom") == "C" and molcoord[C-1].get("atom") == "C"):
            for A in neighbors[B]:
                for D in neighbors[C]:
                    if A == D: continue
                    elif (B == D or A == C): continue
                    temp = tuple(sorted([A,D])+[B,C])
                    if temp not in torsions:
                        torsions.append(temp)
    for torsion in torsions:
        torsions_angles.append(compute_torsion(torsion[0],torsion[2],torsion[3],torsion[1]))
    for torsion_angle in torsions_angles:
        torsions_degree_angles.append(math.degrees(torsion_angle))
    print(torsions_angles)

def V():
    global energies
    energies = []
    def stretch_energy(dists):
        global stretch_energies
        stretch_energies = []
        count = 0
        for dist in dists:
            Atom_a = bonds[count].get("bondinfo")[0]-1
            Atom_b = bonds[count].get("bondinfo")[1]-1
            if (molcoord[Atom_a].get("atom","")=="C" and molcoord[Atom_b].get("atom","")=="C"):
                Kb = dict_parameter["kb_CC"]
                R0 = dict_parameter["CC"]
            elif (molcoord[Atom_a].get("atom","")=="C" and molcoord[Atom_b].get("atom","")=="H"):
                Kb = dict_parameter["kb_CH"]
                R0 = dict_parameter["CH"]
            else:
                raise KeyError("Something went wrong")
            stretch_energies.append(Kb * (dist - R0)**2)
            count = count + 1
        return sum(stretch_energies)
    def bend_energy():
        global bend_energies
        bend_energies = []
        for angle in angles:
            if molcoord[angle[0]-1].get("atom") == "H" and molcoord[angle[2]-1].get("atom") == "H": # H-C-H
                bend_energies.append(dict_parameter["ka_HCH"] * (angle[3] - math.radians(dict_parameter["HCH"]))**2)
            elif molcoord[angle[0]-1].get("atom") == "C" and molcoord[angle[2]-1].get("atom") == "H": # C-C-H
                bend_energies.append(dict_parameter["ka_HCC"] * (angle[3] - math.radians(dict_parameter["HCC"]))**2)
            elif molcoord[angle[0]-1].get("atom") == "H" and molcoord[angle[2]-1].get("atom") == "C": # H-C-C
                bend_energies.append(dict_parameter["ka_HCC"] * (angle[3] - math.radians(dict_parameter["HCC"]))**2)
            elif molcoord[angle[0]-1].get("atom") == "C" and molcoord[angle[2]-1].get("atom") == "C": # C-C-C
                bend_energies.append(dict_parameter["ka_CCC"] * (angle[3] - math.radians(dict_parameter["CCC"]))**2)
            else: continue
        return sum(bend_energies)
    def torsion_energy():
        global torsion_energies
        torsion_energies = []
        for torsion in torsions_angles:
            torsion_energies.append(dict_parameter["Aphi"] * (1 + math.cos(3*torsion)))
        return sum(torsion_energies)

    def LJ():
        global LJ_terms
        global LJ_terms_energies
        LJ_terms = []
        LJ_terms_energies = []
        energy = 0
        for i in range(1,n_atom+1):
            for j in range(1,n_atom+1):
                if i == j: continue
                elif j in neighbors[i]: continue
                skip = False
                for ii,_,jj,_ in angles:
                    if (ii == i and jj == j) or (ii == j and jj == i):
                        skip = True
                        break
                if skip: continue
                if set([i,j]) not in LJ_terms: LJ_terms.append(set([i,j]))
        for LJ_term in LJ_terms:
            a,b = LJ_term
            if molcoord[a-1].get("atom") == "C" and molcoord[b-1].get("atom") == "C":
                epsilon_ij = math.sqrt(dict_parameter["epsilon_C"]**2)
                sigma_ij = 2 * math.sqrt(dict_parameter["sigma_C"]**2)
            elif molcoord[a-1].get("atom") == "H" and molcoord[b-1].get("atom") == "H":
                epsilon_ij = math.sqrt(dict_parameter["epsilon_H"]**2)
                sigma_ij = 2 * math.sqrt(dict_parameter["sigma_H"]**2)
            elif molcoord[a-1].get("atom") != molcoord[b-1].get("atom"):
                epsilon_ij = math.sqrt(dict_parameter["epsilon_H"]*dict_parameter["epsilon_C"])
                sigma_ij = 2 * math.sqrt(dict_parameter["sigma_H"]*dict_parameter["sigma_C"])
            dist = distance(molcoord[a-1].get("coord"),molcoord[b-1].get("coord"))
            LJ_terms_energies.append(4 * epsilon_ij * ((sigma_ij/dist)**12 - (sigma_ij/dist)**6))
        return sum(LJ_terms_energies)
    energies = [stretch_energy(dist_info)+bend_energy()+torsion_energy()+LJ(),stretch_energy(dist_info),bend_energy(),torsion_energy(),LJ()]
    return energies

def writefile():
    with open(sys.argv[2],'w') as f:
        f.write("The input file has: {:d} atoms\n".format(n_atom))
        f.write("Atoms and coordinates (in Å):\n")
        for mcd in molcoord:
            text = "{:s}{:>14.6f}{:>14.6f}{:>14.6f}".format(mcd["atom"],mcd["coord"][0],mcd["coord"][1],mcd["coord"][2])
            f.write(text+"\n")
        text = """Number of coordinates:
Stretching:{:>7d} Bending:{:>7d} Torsion:{:>7d}
Internal:{:>7d} Cartesian:{:>7d} 
Potential energy at input structure:
   {:.6f} kcal/mol
Stretch, Bend, Torsion, VDW components of potential energy:
{:12.6f}{:13.6f}{:13.6f}{:13.6f}""".format(n_bond,len(angles),len(torsions),n_bond+len(angles)+len(torsions),3*n_atom,energies[0],energies[1],energies[2],energies[3],energies[4])
        f.write(text+"\n")
        f.write("List of all bonds: (At1 - At2, with labels, and distance in Angstrom, energy contrib in kcal/mol)\n")
        count = 0
        for bond in bonds:
            a = bond.get("bondinfo")[0]
            b = bond.get("bondinfo")[1]
            text = """{:s}{:>4d}   -   {:s}{:>4d}:{:>13.5f}{:>14.5f}
""".format(molcoord[a-1].get("atom"),a,molcoord[b-1].get("atom"),b,dist_info[count],stretch_energies[count])
            count = count + 1
            f.write(text) 
        f.write("List of all bending angles: (At1 - At2 - At3, with labels, angle in radian then degrees, energy contribution)\n")
        count = 0
        for angle in angles:
            a = angle[0]
            b = angle[1]
            c = angle[2]
            text = """{:s}{:>4d}   -   {:s}{:>4d}   -   {:s}{:>4d}:{:>13.5f}{:>11.3f}{:>15.5f}
""".format(molcoord[a-1].get("atom"),a,molcoord[b-1].get("atom"),b,molcoord[c-1].get("atom"),c,angle[3],degree_angles[count],bend_energies[count])
            count = count + 1
            f.write(text)
        f.write("List of all torsional angles: (At1 - At2 - At3 - At4, with labels, angle in radian then degrees,  energy contrib in kcal/mol)\n")
        count = 0
        for torsion in torsions:
            a = torsion[0]
            b = torsion[2]
            c = torsion[3]
            d = torsion[1]
            text = """{:s}{:>4d}   -   {:s}{:>4d}   -   {:s}{:>4d}   -   {:s}{:>4d}:{:>13.6f}{:>11.3f}{:>14.6f}
""".format(molcoord[a-1].get("atom"),a,molcoord[b-1].get("atom"),b,molcoord[c-1].get("atom"),c,molcoord[d-1].get("atom"),d,torsions_angles[count],torsions_degree_angles[count],torsion_energies[count])
            count = count + 1
            f.write(text)
        f.write("List of all unique atom pairs: (At1 - At2, with labels, distance,  vdW energy contrib in kcal/mol)\n")
        for i in range(1,n_atom+1):
            for j in range(1,n_atom+1):
                ind = set([i,j])
                if i >= j: continue
                elif ind not in LJ_terms:
                    a1,a2,a3 = molcoord[i-1].get("coord")
                    b1,b2,b3 = molcoord[j-1].get("coord")
                    dist = math.sqrt((a1 - b1)**2 + (a2 - b2)**2 + (a3 - b3)**2)
                    text = """{:s}{:4d}   -   {:s}{:4d}:{:13.5f}{:14.5f}
""".format(molcoord[i-1].get("atom"),i,molcoord[j-1].get("atom"),j,dist,0)
                    f.write(text)
                elif ind in LJ_terms:
                    a1,a2,a3 = molcoord[i-1].get("coord")
                    b1,b2,b3 = molcoord[j-1].get("coord")
                    dist = math.sqrt((a1 - b1)**2 + (a2 - b2)**2 + (a3 - b3)**2)
                    text = """{:s}{:4d}   -   {:s}{:4d}:{:13.5f}{:14.5f}
""".format(molcoord[i-1].get("atom"),i,molcoord[j-1].get("atom"),j,dist,LJ_terms_energies[LJ_terms.index(ind)])
                    f.write(text)
def main():
    readfile()        
    bondlength()
    angle()
    torsion()
    V()
    writefile()

    print("----------------------")  
    print(moldata)
    print("----------------------")  
    for i in molcoord:
        print(i,end="\n")
    print("----------------------")  
    for i in bonds:
        print(i,end="\n")



main()