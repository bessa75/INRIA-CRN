from recherche_chemin import *
import re

def numero3(texte):
    if type(texte) is str:
        return numero(texte,liste_molecules)
    if type(texte) is list:
        return numero2(texte,liste_molecules)

"""    Lecture du fichier catalogue"""


file = "catalog.bc"
nb_réactions_max = 50

elmts = []
inputs = []
enzymes = []
regex = ["(present\()|(, [e\-0-9]+\))|\.", "(MA.*for )|\+|(=>)|\."]
blocs = {"inputs": (regex[0], inputs), "enzymes": (regex[0], enzymes), "elmts": (regex[1], elmts)}

""" Création du tableau "reaction" de l'ensemble des réactions possible données par le fichier file"""
friends = {'test': []}  # Chaque elmt et l'ensemble des elements avec lesquels il reagit
# par contre les enzymes ne sont pas comptés comme éléments, du coup ils ne font pas partie des clés du dict.

reaction = []  # Liste des réactions en version texte
with open(file) as f:
    while f.readline() != "% Inputs\n":
        continue
    bloc = 'inputs'
    while True:
        l = f.readline()
        if not l:
            break

        if l[0] == '%':
            if l.find("Enzymes") > 0:
                bloc = 'enzymes'
                continue
            if l.find('reaction') > 0:
                bloc = 'elmts'
                continue
            continue

        if len(l) < 3:
            continue

        if bloc == 'elmts':
            elmt = l.split('=>')
            reactifs = re.sub("(MA.*for )|\+", " ", elmt[0]).split()
            for r in reactifs:
                try:
                    friends[r].extend([e for e in reactifs if e not in friends[r]])
                except:
                    friends.update({r: [e for e in reactifs if e != r]})
            produits = re.sub("\+|\.", " ", elmt[1]).split()
            reaction.append((reactifs, produits))
            elts = reactifs + produits
        else:
            elts = re.sub(blocs[bloc][0], '', l).split()

        blocs[bloc][1].extend([e for e in elts if (e not in blocs[bloc][1])])

liste_molecules = elmts
for molecule in enzymes :
    if not(molecule in enzymes):
        liste_molecules.append(molecule)
liste_molecules.pop(numero3('H_2O_2'))  # Remove 'H_2O_2', 'H2O2' is in list

dictionnaire_molecules={liste_molecules[i] : i for i in range(len(liste_molecules))}

REACTIONS = []  # Liste des réactions en version numéros
for a in reaction:
    reac = ([], [])
    for m in a[0]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[0].append(numero('H2O2',liste_molecules))
        else:
            reac[0].append(numero(m,liste_molecules))
    for m in a[1]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[1].append(numero('H2O2',liste_molecules))
        else:
            reac[1].append(numero(m,liste_molecules))
    REACTIONS.append(reac)


REACPARMOL = [] # REACPARMOL = ?
for k in range(0, len(liste_molecules)):
    REACPARMOL.append([])
for k in range(0, len(REACTIONS)):
    for u in REACTIONS[k][0]:
        REACPARMOL[u].append(k)

REACPARMOL2=[]#cette liste contient pour chaque molécule, les numéros des réactions dans lesquels elle est un réactif
REACPARMOLP=[]#cette liste contient pour chaque molécules, les numéros des réactions qui la produisent
MOL=liste_molecules
ENZYMES=numero2(enzymes,MOL)
CYCLES=[]
CYCLESPARMOL=[]
BoolCycles=[False]*len(REACTIONS)

for k in range (len(MOL)):
    REACPARMOL2.append([])
    CYCLESPARMOL.append([])
    REACPARMOLP.append([])
for k in range (0,len(REACTIONS)):
    for i in range (0,len(REACTIONS)):
        reac1=REACTIONS[k]
        reac2=REACTIONS[i]
        bol1=True
        if len(reac1[0])==len(reac2[1]) and len(reac1[1])==len(reac2[0]):
            for a in reac2[0]:
                if a not in reac1[1]:
                    bol1=False
            for a in reac2[1]:
                if a not in reac1[0]:
                    bol1=False
            if bol1 and (reac1 not in CYCLES):
                BoolCycles[k]=True
                BoolCycles[i]=True
                CYCLES.append(reac1)
                for a in reac1[0]:
                    CYCLESPARMOL[a].append(k)
for k in range (0,len(REACTIONS)):
    reac=REACTIONS[k]
    reactifs=reac[0]
    REACPARMOL2[reactifs[0]].append(k)
    if len(reactifs)==2:
        REACPARMOL2[reactifs[1]].append(k)
    produits=reac[1]
    REACPARMOLP[produits[0]].append(k)
    if len(produits)==2:
        REACPARMOLP[produits[1]].append(k)
#cas 1 :
ENZ1=numero2(['AO', 'ADH', 'G_1DH', 'resazurin', 'HRP', 'H2O2'],MOL)
#cas 2 :
ENZ2=numero2(['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO'],MOL)
#cas 3 :
ENZ3=numero2(['ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'G_1DH', 'O2', 'DAF', 'NAD'],MOL)

def algo_negation(A,B,C,ENZ,CYCLES,CYCLESPARMOL,REACPARMOL,reaction,MOL,REACTIONS,n):
    print('********************************************')
    PRODINT=[]
    REACT=[]
    ACONSOMMER=[C]
    CYCLESa=[False]
    BOOL=[]
    L=[]
    CYCLEC=CYCLESPARMOL[C]
    RES=[]
    for k in range (0,len(CYCLEC)):
        #print(k)
        #print(CYCLEC[k])
        bool=False
        m1,m2=REACTIONS[CYCLEC[k]]
        if len(m1)==1 and C in m2 and m1[0] not in ACONSOMMER:
            ACONSOMMER.append(m1[0])
            bool=True
            CYCLESa.append(CYCLEC[k])
        if len(m2)==1 and C in m1 and m2[0] not in ACONSOMMER:
            ACONSOMMER.append(m2[0])
            if bool==False:
                CYCLESa.append(CYCLEC[k])
    print("Molécules dont la consommation permettrait de consommer " + str(MOL[C])+" : "+str([MOL[i] for i in ACONSOMMER]))
    for i in range (0,len(ACONSOMMER)):
        M=ACONSOMMER[i]
        REACC=REACPARMOL[M]
        CYCLEC=CYCLESPARMOL[M]
        for k in range (0,len(REACC)):
            r,p=REACTIONS[REACC[k]]
            if len(r)==1:
                if (p,r) not in CYCLEC and (r,p) not in CYCLEC and M==C:
                    return("C est une molecule instable")
            elif len(r)==2: #on considere la molecularite inferieure ou egale a 2
                im=min([i for i in range (0,len(r)) if r[i]!=M])
                prodint=r[im]
                if r[im]==B:
                    print("Fin de l'algorithme, voici le résultat :'")
                    return(([[REACC[k],CYCLESa[i]],'a',[]]))
                PRODINT.append(prodint)
                REACT.append(REACC[k])
                if M==C:
                    BOOL.append(False)
                else:
                    BOOL.append((True,i))

    set1= {}
    set2= {}
    PRESENCE=pluscourtchemin(ENZ,[B],0,set1,set2,n,False,reac,MOL,REACTIONS,REACPARMOL,bool=True,bool2=True)
    ENZint=[]
    for k in range (0,len(PRODINT)):
        m=PRODINT[k]
        #print(MOL[m])
        PRESm=PRESENCE[m]
        #print(PRESm)
        #print(PRESm[0])
        if PRESm[0] and 'a' in PRESm[2]:
            print("produit intermédiaire considéré : "+MOL[m])
            P2=pluscourtchemin(ENZ,[B,0],m,set1,set2,n,False,reac,MOL,REACTIONS,REACPARMOL,bool2=True)
            if P2!=False:
                for mec in P2:
                    if mec[1]=='a':
                        meca=mec[0]
                        meca.append([REACT[k]])
                        if BOOL[k]!=False:
                            meca.append([CYCLESa[BOOL[k][1]]])
                        Res=(meca,'a',mec[2])
                        RES.append(Res)
                    if mec[1]=='e':
                        meca=mec[0]
                        meca.append([REACT[k]])
                        if BOOL[k]!=False:
                            meca.append([CYCLESa[BOOL[k][1]]])
                        Res=(meca,'e',mec[2])
                        RES.append(Res)
                        ENZint.append([MOL[i] for i in Res[2]])
    """ PRESENCE : liste des molécules présente.
    Evolue au cours de l'algorithme.
    Chaque case coorespond à une molécule et est de la forme :
    [boolean marquant la présence,
    liste des réaction l'ayant produit selon le triplet (numéro de réaction, numéro d'étape, étiquette),
    liste des étiquettes avec lesquelles la molécule a été produite,
    liste des doublets (numéro d'étape, étiquette) avec lesquelles la molécule a été produite]"""
    print("Jeux d'enzymes interdites : "+str(ENZint))
    print("Fin de l'algorithme: voici le résultat :")
    print("")
    return(RES)

#print(numero('ABTSOX',MOL))

print(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ1,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0],reaction))
print(algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
print(algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ3,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
#print(reaction[a])



