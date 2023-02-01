import recherche_chemin
import algoresolution_système

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


print("Temps d'initialisation : ",time.time()-start_time)
    
"""    Résolution des 3 exemples"""

def test_1():
    ##glucose et acetone donnent gluconolacrone
    ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    MECAS = recherche_chemin.résolution_équation(ENZ,"acetoneext + glucoseext => resorufin",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)


def test_2():
## Test lescture équation
    ENZ=['ABTS','ADH', 'NADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF','NAD'] ##rajouter NAD pour fausser le résultat
    RE=['acetoneext','glucoseext']
    re=numero(RE)
    enz=numero(ENZ)

    solution=algoresolution_système.ressystem([re],[numero('gluconolacrone')],['ab'],20,enz,reaction,liste_molecules,REACTIONS,REACPARMOL)

    mt = recherche_chemin.mecatexte(solution[0],reaction)
    for d in mt[0]:
        print(d)
    print("")
    print("")

    
def test_3():
##NO et glucose donnent gluconolacrone avec NO3 en réactif annexe
    ENZ = ['ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'NAD']
    recherche_chemin.résolution_équation(ENZ,"NO2 + glucoseext => DAFF",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)


def test_4():
## Ancien OU logique à changer sur le fonctionnement pluscourtchemin
    ENZ=numero(['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO'])
    #recherche_chemin.ENZ = ['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO']
    RE = numero(['Lactateext', 'EtOHext'])

    recherche_chemin.résolution_équation(ENZ,"Lactateext + EtOHext => ABTSOX",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)

    MECAS = recherche_chemin.pluscourtchemin(ENZ, RE, numero('ABTSOX'), nb_réactions_max, True,reaction,liste_molecules,REACTIONS,REACPARMOL)  # Pourquoi tag a ? OU logique ?
    mt = recherche_chemin.mecatexte(MECAS[0][0],reaction)
    for d in mt:
        print(d)


def test_5():
##glucose et Non(acetone) donnent gluconolacrone
    ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    recherche_chemin.résolution_équation(ENZ,"!acetoneext + glucoseext => NADH",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)

    
def test_6():
##glucose et acetone donnent gluconolacrone
    ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    recherche_chemin.résolution_équation(ENZ,"acetoneext + glucoseext => resorufin",nb_réactions_max,reaction,liste_molecules,REACTIONS,REACPARMOL)


start_time = time.time()

test_1()
print()
print('------------------------------------------------------------------')
print()
test_2()
print()
print('------------------------------------------------------------------')
print()
test_3()
print()
print('------------------------------------------------------------------')
print()
test_5()
print()
print('------------------------------------------------------------------')
print()
test_6()
print()
print('------------------------------------------------------------------')
print()

print("Temps de calcul des exemples : ",time.time()-start_time)
