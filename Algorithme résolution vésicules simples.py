import recherche_chemin
import algoresolution_système
import algonegation
#bibliothéque python
import re
import time

start_time = time.time()

def numero(texte):
    if type(texte) is str:
        return recherche_chemin.numero(texte,liste_molecules)
    if type(texte) is list:
        return recherche_chemin.numero2(texte,liste_molecules)

"""    Déclaration variables et constantes"""


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
liste_molecules.pop(numero('H_2O_2'))  # Remove 'H_2O_2', 'H2O2' is in list

dictionnaire_molecules={liste_molecules[i] : i for i in range(len(liste_molecules))}

REACTIONS = []  # Liste des réactions en version numéros
for a in reaction:
    reac = ([], [])
    for m in a[0]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[0].append(recherche_chemin.numero('H2O2',liste_molecules))
        else:
            reac[0].append(recherche_chemin.numero(m,liste_molecules))
    for m in a[1]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[1].append(recherche_chemin.numero('H2O2',liste_molecules))
        else:
            reac[1].append(recherche_chemin.numero(m,liste_molecules))
    REACTIONS.append(reac)


REACPARMOL = [] # REACPARMOL = ?
for k in range(0, len(liste_molecules)):
    REACPARMOL.append([])
for k in range(0, len(REACTIONS)):
    u = REACTIONS[k][0][0]
    REACPARMOL[u].append(k)
recherche_chemin.REACPARMOL=REACPARMOL

REACPARMOL2=[]#cette liste contient pour chaque molécule, les numéros des réactions dans lesquels elle est un réactif
REACPARMOLP=[]#cette liste contient pour chaque molécules, les numéros des réactions qui la produisent
ENZYMES=recherche_chemin.numero2(enzymes,liste_molecules)
CYCLES=[]
CYCLESPARMOL=[]
BoolCycles=[False]*len(REACTIONS)
for k in range (len(liste_molecules)):
    REACPARMOL2.append([])
    CYCLESPARMOL.append([])
    REACPARMOLP.append([])
for k in range (0,len(REACTIONS)):
    for i in range (0,len(REACTIONS)):
        reac1=REACTIONS[k]
        reac2=REACTIONS[i]
        if (reac2[0],reac2[1])==(reac1[1],reac1[0]) and (reac1 not in CYCLES):
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
ENZ1=recherche_chemin.numero2(['ABTS','NAD', 'resazurin', 'HRP','NR', 'AO', 'POD', 'G_1DH', 'O2', 'DAF'],liste_molecules)
#cas 2 :
ENZ2=recherche_chemin.numero2(['ABTS','NAD', 'resazurin', 'HRP','NR', 'AO', 'POD', 'G_1DH', 'O2', 'DAF','ADH'],liste_molecules)
#cas 3 :
ENZ3=recherche_chemin.numero2(['ABTS','NAD', 'resazurin', 'HRP','NR', 'AO', 'POD', 'G_1DH', 'O2', 'DAF'],liste_molecules)
print("Temps d'initialisation : ",time.time()-start_time)

"""    Résolution des 3 exemples"""

'''
print(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ1,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ1,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0],reaction))
print(algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
print(algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ3,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ3,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
#print(reaction[a])
'''

def test_1():
##glucose et acetone donnent resorufin
    print("Test 1 :")
    ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    MECAS = recherche_chemin.résolution_équation(ENZ,"acetoneext + glucoseext => resorufin",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)
    print()
    print('------------------------------------------------------------------')
    print()


def test_2():
## Test lescture équation
    print("Test 2 :")
    ENZ=['ABTS','ADH', 'NADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF','NAD'] ##rajouter NAD pour fausser le résultat
    RE=['acetoneext','glucoseext']
    re=numero(RE)
    enz=numero(ENZ)

    solution=algoresolution_système.res([re],[numero('gluconolacrone')],['ab'],20,enz,reaction,liste_molecules,REACTIONS,REACPARMOL,reac,CYCLES,CYCLESPARMOL)
    print(solution)
    mt = recherche_chemin.mecatexte(solution[0],reaction)
    for d in mt[0]:
        print(d)
    print("")
    print("")
    print()
    print('------------------------------------------------------------------')
    print()


def test_3():
##NO et glucose donnent gluconolacrone avec NO3 en réactif annexe
    print("Test 3 :")
    ENZ = ['ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'NAD']
    recherche_chemin.résolution_équation(ENZ,"NO2 + glucoseext => DAFF",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)
    print()
    print('------------------------------------------------------------------')
    print()


def test_4():
## Ancien OU logique à changer sur le fonctionnement pluscourtchemin
    print("Test 4 :")
    ENZ=['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO']
    #recherche_chemin.ENZ = ['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO']
    RE = numero(['Lactateext', 'EtOHext'])

    recherche_chemin.résolution_équation(ENZ,"Lactateext + EtOHext => ABTSOX",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)

    MECAS = recherche_chemin.pluscourtchemin(numero(ENZ), RE, numero('ABTSOX'), nb_réactions_max, True,reaction,liste_molecules,REACTIONS,REACPARMOL)
    mt = recherche_chemin.mecatexte(MECAS[0][0],reaction)
    for d in mt:
        print(d)
    print()
    print('------------------------------------------------------------------')
    print()


def test_5():
##glucose et Non(acetone) donnent gluconolacrone
    print("Test 5 :")
    ENZ = ['AO', 'ADH', 'G_1DH', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    RE=['acetoneext','glucoseext']
    re=numero(RE)
    enz=numero(ENZ)

    solution=algoresolution_système.res([[numero('glucose'),numero('acetone')]],[numero('NADH')],['anb'],20,ENZ1,reaction,liste_molecules,REACTIONS,REACPARMOL,reac,CYCLES,CYCLESPARMOL)
    print(solution)
    '''
    mt = recherche_chemin.mecatexte(solution[0],reaction)
    for d in mt[0]:
        print(d)
    '''
    print()
    print('------------------------------------------------------------------')
    print()


start_time = time.time()
test_1()
#test_2()
#test_3()
#test_4()
#test_5()

#print(numero('ABTSOX',MOL))

#print(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ1,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5))
#print(mecatexte(algo_negation(numero('glucose',MOL),numero('acetone',MOL),numero('NADH',MOL),ENZ1,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0],reaction))
#print(algonegation.algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5,reac))
#print(mecatexte(algo_negation(numero('Lactateext',MOL),numero('EtOHext',MOL),numero('ABTSOX',MOL),ENZ2,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
#print(algonegation.algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ3,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5,reac))
#print(mecatexte(algo_negation(numero('glucoseext',MOL),numero('NO3ext',MOL),numero('NADH',MOL),ENZ3,CYCLES,CYCLESPARMOL,REACPARMOL2,reaction,MOL,REACTIONS,5)[0][0],reaction))
#print(reaction[a])

print("Temps de calcul des exemples : ",time.time()-start_time)
