#bibliothéque python
import re
import time
import pandas as pd
import numpy as np
#importation code perso
import recherche_chemin
import algoresolution_système
from creation_CRN_v2 import *
from brenda import read_data

def numero(texte):
    if type(texte) is str:
        return recherche_chemin.numero(texte,listeNomsMolecules)
    if type(texte) is list:
        return recherche_chemin.numero2(texte,listeNomsMolecules)

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

def base_cataloge():
    # Utilisation du catalogue
    file = "catalog.bc"
    enzymes, elmts, reaction = get_data_from_file(file)
    listeNomsMolecules = enzymes + elmts
    listeNomsMolecules.pop(26)

    listeReactionsBrut = []  # Liste des réactions en version numéros
    for a in reaction:
        reac = ([], [])
        for m in a[0]:
            if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
                reac[0].append(recherche_chemin.numero('H2O2',listeNomsMolecules))
            else:
                reac[0].append(recherche_chemin.numero(m,listeNomsMolecules))
        for m in a[1]:
            if m == 'H_20_2' or m == 'H202' or m == 'H_2O_2':
                reac[1].append(recherche_chemin.numero('H2O2',listeNomsMolecules))
            else:
                reac[1].append(recherche_chemin.numero(m,listeNomsMolecules))
        listeReactionsBrut.append(reac)
    # Fin  preprocessing du Fichier
            
    return listeReactionsBrut,listeNomsMolecules

def base_brenda():
    # Utilisation de Brenda
    print("Start read data")
    data = read_data('data.json')
    listeReactionsBrut = data['reactions'] # [ [[liste des reactifs],[liste des produits],[kms],[kcats],[phs]] autant de fois qu'il y a de réactions]
    listeNomsMolecules = data['molecules_list'] # inclus molecules et enzymes
    dictionnaireNomsMolecules = data['dic_name_to_idx'] # Pour avoir l'indice à partir du nom
    #listeReactionParMolecules = data['reactions_per_molecule']
    print("End read data")
    return listeReactionsBrut,listeNomsMolecules,dictionnaireNomsMolecules


"""    Déclaration tableaux de réactions et molécules"""
listeReactionsBrut,listeNomsMolecules = base_cataloge()
#listeReactionsBrut,listeNomsMolecules,dictionnaireNomsMolecules = base_brenda()

listeReactionParMolecules = [[] for k in range(len(listeNomsMolecules))]
for numReaction, reaction in enumerate(listeReactionsBrut):
    for reactif in reaction[0]:
        listeReactionParMolecules[reactif].append(numReaction)


"""    Résolution des 3 exemples"""

def test_1():
##glucose et acetone donnent resorufin
    print("Test 1 :")
    start_time = time.time()

    EspecesInitiales=numero(['AO', 'ADH', 'G_1DH', 'NAD', 'resazurin', 'HRP', 'H2O2'])
    #EspecesInitiales=numero(['antidiuretic hormone', 'beta-d-glucose', 'NAD+', 'resazurin', 'peroxidase', 'H2O2']) #Verrsion brenda
    
    nbEtapeMax=15

    A=numero('acetoneext')
    B=numero('glucoseext')
    #Brenda name version
    #A=numero('2-propanone')
    #B=numero('d-glucose')

    C=numero('resorufin')

    MelangesInitiaux,res=déterminationCRN (A,B,C,nbEtapeMax,"ET",EspecesInitiales,listeReactionsBrut,listeReactionParMolecules)

    #Affichage des résultats
    affichage_mélanges_i(MelangesInitiaux,C,listeNomsMolecules)
    print()
    print("Affichage du CRN trouvé")
    affichage_CRN(res[1],listeReactionsBrut,listeReactionParMolecules,listeNomsMolecules)

    print("Temps de calcul : ",time.time()-start_time)
    print()
    print('------------------------------------------------------------------')
    print()

def test_2():
##NO et glucose donnent gluconolacrone avec NO3 en réactif annexe
    print("Test 2 :")
    start_time = time.time()

    EspecesInitiales=numero(['H2O2','acetoneext', 'ABTS', 'ADH', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'NAD'])
    nbEtapeMax=50
    A=numero('NO2')
    B=numero('glucoseext')
    C=numero('DAFF')

    MelangesInitiaux,res=déterminationCRN (A,B,C,nbEtapeMax,'ET',EspecesInitiales,listeReactionsBrut,listeReactionParMolecules)

    #Affichage des résultats
    affichage_mélanges_i(MelangesInitiaux,C,listeNomsMolecules)
    print()
    print("Affichage du CRN trouvé")
    affichage_CRN(res[1],listeReactionsBrut,listeReactionParMolecules,listeNomsMolecules)

    print("Temps de calcul : ",time.time()-start_time)
    print()
    print('------------------------------------------------------------------')
    print()

def test_3():
## Lactateext ET EtOHext donnent ABTSOX
    print("Test 3 :")
    start_time = time.time()
    
    EspecesInitiales=numero(['ABTS', 'ADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'LO'])
    nbEtapeMax=50
    A=numero('Lactateext')
    B=numero('EtOHext')
    C=numero('ABTSOX')

    MelangesInitiaux,res=déterminationCRN (A,B,C,nbEtapeMax,'ET',EspecesInitiales,listeReactionsBrut,listeReactionParMolecules)

    #Affichage des résultats
    affichage_mélanges_i(MelangesInitiaux,C,listeNomsMolecules)
    print()
    print("Affichage du CRN trouvé")
    affichage_CRN(MelangesInitiaux[C][0],listeReactionsBrut,listeReactionParMolecules,listeNomsMolecules)

    print("Temps de calcul : ",time.time()-start_time)
    print()
    print('------------------------------------------------------------------')
    print()

test_1()
#test_2()
#test_3()
