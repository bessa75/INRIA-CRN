import algoPSC2_0_1
import algoresolution_système
from creation_CRN_v2 import *
from brenda import get_data
import re #bibliothéque python

def numero(texte):
    if type(texte) is str:
        return algoPSC2_0_1.numero(texte,listeMoleculeTexte)
    if type(texte) is list:
        return algoPSC2_0_1.numero2(texte,listeMoleculeTexte)


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
            reac[0].append(algoPSC2_0_1.numero('H2O2',listeMoleculeTexte))
        else:
            reac[0].append(algoPSC2_0_1.numero(m,listeMoleculeTexte))
    for m in a[1]:
        if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
            reac[1].append(algoPSC2_0_1.numero('H2O2',listeMoleculeTexte))
        else:
            reac[1].append(algoPSC2_0_1.numero(m,listeMoleculeTexte))
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

"""
##glucose et acetone donnent gluconolacrone
algoPSC2_0_1.ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
MECAS = algoPSC2_0_1.résolution_équation("acetoneext + glucoseext => resorufin",nb_réactions_max,reaction,MOL,REACTIONS,REACPARMOL)
"""


EspecesInitiales=numero(['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2','acetoneext', 'glucoseext'])
nbEtapeMax=50
A=numero('acetoneext')
B=numero('glucoseext')
C=numero('resorufin')
MelangesInitiaux,melange_C=recherche_melanges_initiaux(EspecesInitiales,nbEtapeMax,C,listeReactions)
#print(déterminationCRN (A,B,C,nbEtapeMax,"ET",EspecesInitiales,listeReactions))

#affichage_tous_mélanges(MelangesInitiaux,listeMoleculeTexte)
affichage_mélanges_i(MelangesInitiaux,C,listeMoleculeTexte)
print()

"""
## Test lescture équation
ENZ=['ABTS','ADH', 'NADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF','NAD'] ##rajouter NAD pour fausser le résultat
RE=['acetoneext','glucoseext']
re=numero(RE)
enz=numero(ENZ)

solution=algoresolution_système.ressystem([re],[numero('gluconolacrone')],['ab'],20,enz,reaction,MOL,REACTIONS,REACPARMOL)

mt = algoPSC2_0_1.mecatexte(solution[0],reaction)
for d in mt[0]:
    print(d)
print("")
print("")
"""


"""
##NO et glucose donnent gluconolacrone avec NO3 en réactif annexe
algoPSC2_0_1.ENZ = ['ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'NAD']
algoPSC2_0_1.résolution_équation("NO2 + glucoseext => DAFF",nb_réactions_max,reaction,MOL,REACTIONS,REACPARMOL)
"""


"""
## Ancien OU logique à changer sur le fonctionnement pluscourtchemin

ENZ=numero(['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO'])
algoPSC2_0_1.ENZ = ['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO']
RE = numero(['Lactateext', 'EtOHext'])

algoPSC2_0_1.résolution_équation("Lactateext + EtOHext => ABTSOX",nb_réactions_max,reaction,MOL,REACTIONS,REACPARMOL)

MECAS = algoPSC2_0_1.pluscourtchemin(ENZ, RE, numero('ABTSOX'), nb_réactions_max, True,reaction,MOL,REACTIONS,REACPARMOL)  # Pourquoi tag a ? OU logique ?
mt = algoPSC2_0_1.mecatexte(MECAS[0][0],reaction)
for d in mt:
    print(d)
"""


"""
##glucose et Non(acetone) donnent gluconolacrone
algoPSC2_0_1.ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
algoPSC2_0_1.résolution_équation("!acetoneext + glucoseext => NADH",nb_réactions_max,reaction,MOL,REACTIONS,REACPARMOL)
"""

"""
##glucose et acetone donnent gluconolacrone
algoPSC2_0_1.ENZ = ['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2']  ##rajouter NAD pour fausser le résultat
algoPSC2_0_1.résolution_équation("acetoneext + glucoseext => resorufin",nb_réactions_max,reaction,MOL,REACTIONS,REACPARMOL)
"""