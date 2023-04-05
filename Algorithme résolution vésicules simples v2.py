import recherche_chemin
import algoresolution_système
from creation_CRN_v2 import *
from brenda import get_data
#bibliothéque python
import re
import time

def numero(texte):
    if type(texte) is str:
        return recherche_chemin.numero(texte,listeMoleculeTexte)
    if type(texte) is list:
        return recherche_chemin.numero2(texte,listeMoleculeTexte)

"""    Déclaration variables et constantes"""


def get_data_from_file(file):
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

    return elmts, enzymes, reaction


"""    Déclaration variables et constantes"""
nb_réactions_max = 50


# Utilisation de Brenda
data = get_data()
listeReactionsBrut = data['reactions']
listeMoleculeTexte = data['molecules'] # inclus molecules et enzymes
# Fin brenda


# Utilisation du catalogue
file = "catalog.bc"
enzymes, elmts, reaction = get_data_from_file(file)
listeMoleculeTexte = enzymes + elmts
listeMoleculeTexte.pop(26)

N=len(listeMoleculeTexte) # Nombre de molecules


listeReactionsBrut = []  # Liste des réactions en version numéros
for a in reaction:
    reac = ([], [])
    for m in a[0]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[0].append(recherche_chemin.numero('H2O2',listeMoleculeTexte))
        else:
            reac[0].append(recherche_chemin.numero(m,listeMoleculeTexte))
    for m in a[1]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[1].append(recherche_chemin.numero('H2O2',listeMoleculeTexte))
        else:
            reac[1].append(recherche_chemin.numero(m,listeMoleculeTexte))
    listeReactionsBrut.append(reac)

# Fin  preprocessing du Fichier

listeReactions = [[[] for j in range(N)] for i in range(N)] # liste où chaque case i,j est (numéro de réaction, listeProduits)

for num_reaction, reaction in enumerate(listeReactionsBrut):
    listeReactifs=reaction[0]
    listeProduits=reaction[1]
    if len(listeReactifs)==1:
        listeReactions[listeReactifs[0]][listeReactifs[0]]=listeReactions[listeReactifs[0]][listeReactifs[0]]+listeProduits
    elif len(listeReactifs)==2:
        listeReactions[listeReactifs[0]][listeReactifs[1]]=listeReactions[listeReactifs[0]][listeReactifs[1]]+listeProduits
        listeReactions[listeReactifs[1]][listeReactifs[0]]=listeReactions[listeReactifs[1]][listeReactifs[0]]+listeProduits
    else:
        print('--- Plus que 2 réactifs---')
    
    

"""    Résolution des 3 exemples"""


def test_1():
    start_time = time.time()
    ##glucose et acetone donnent gluconolacrone

    EspecesInitiales=numero(['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2','acetoneext', 'glucoseext'])
    nbEtapeMax=50
    A=numero('acetoneext')
    B=numero('glucoseext')
    C=numero('resorufin')

    MelangesInitiaux,melange_C=recherche_melanges_initiaux(EspecesInitiales,nbEtapeMax,C,listeReactions)
    #print(déterminationCRN (A,B,C,nbEtapeMax,"ET",EspecesInitiales,listeReactions))

    affichage_mélanges_i(MelangesInitiaux,C,listeMoleculeTexte)
    print("Temps de calcul des exemples : ",time.time()-start_time)
    print()

def test_2():
    ## Test lescture équation
    ENZ=['ABTS','ADH', 'NADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF','NAD']
    RE=['acetoneext','glucoseext']
    re=numero(RE)
    enz=numero(ENZ)

    solution=algoresolution_système.res([re],[numero('gluconolacrone')],['ab'],20,enz,reaction,listeMoleculeTexte,REACTIONS,REACPARMOL,reac,CYCLES,CYCLESPARMOL)
    print(solution)
    mt = recherche_chemin.mecatexte(solution[0],reaction)
    for d in mt[0]:
        print(d)

def test_3():
    start_time = time.time()
    ##NO et glucose donnent gluconolacrone avec NO3 en réactif annexe
    EspecesInitiales=numero(['H2O2','acetoneext', 'ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'NAD'])
    nbEtapeMax=50
    A=numero('NO2')
    B=numero('glucoseext')
    C=numero('DAFF')

    MelangesInitiaux,res=déterminationCRN (A,B,C,nbEtapeMax,'ET',EspecesInitiales,listeReactions)
    affichage_mélanges_i(MelangesInitiaux,C,listeMoleculeTexte)
    print("Temps de calcul des exemples : ",time.time()-start_time)
    print()

def test_4():
    start_time = time.time()
    ## Ancien OU logique à changer sur le fonctionnement pluscourtchemin
    EspecesInitiales=numero(['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO'])
    nbEtapeMax=50
    A=numero('Lactateext')
    B=numero('EtOHext')
    C=numero('ABTSOX')

    MelangesInitiaux,res=déterminationCRN (A,B,C,nbEtapeMax,'ET',EspecesInitiales,listeReactions)
    affichage_mélanges_i(MelangesInitiaux,C,listeMoleculeTexte)
    print("Temps de calcul des exemples : ",time.time()-start_time)
    print()

def test_5():
    ##glucose et Non(acetone) donnent gluconolacrone
    ENZ = ['AO', 'ADH', 'G_1DH', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
    RE=['acetoneext','glucoseext']
    re=numero(RE)
    enz=numero(ENZ)

    solution=algoresolution_système.res([[numero('glucose'),numero('acetone')]],[numero('NADH')],['anb'],20,ENZ1,reaction,listeMoleculeTexte,REACTIONS,REACPARMOL,reac,CYCLES,CYCLESPARMOL)
    print(solution)


start_time = time.time()
print("Test 1 :")
test_1()
print()
print('------------------------------------------------------------------')
print()
"""
print("Test 2 :")
test_2()
print()
print('------------------------------------------------------------------')
print()
"""
print("Test 3 :")
test_3()
print()
print('------------------------------------------------------------------')
print()
print("Test 4 :")
test_4()
print()
print('------------------------------------------------------------------')
print()
"""
print("Test 5 :")
test_5()
print()
print('------------------------------------------------------------------')
print()
"""