def recherche_melanges_initiaux(EspecesInitiales,nbEtapeMax,moleculeCherche,listeReactionsBrut,listeReactionParMolecules):

    # I) Initialisation
    N=len(listeReactionParMolecules)
    nbEtape=0 # Numero de l'etape en cours
    
    moleculesPresentes=[numMol for numMol in EspecesInitiales] # Liste de toutes les molecules crees depuis le debut
    moleculesCreees=[numMol for numMol in EspecesInitiales] # Liste des molecules cree à l'etape precedente
    nouvellesMolecule=[] # Tableau tampon pour la boucle

    melangesInitiaux=[[] for numMol in range(N)] # Tableau où chaque case i est la liste des doublets (n, combinaisons d'especes initiales permettant d'obtenir la molecule i en n etape) 
    melangesInitiauxCrees=[[] for numMol in range(N)] # Idem melangesInitiaux mais que ceux creees à l'etape precedente
    for numMol in EspecesInitiales:
        melangesInitiaux[numMol]=[(nbEtape,[numMol])]
        melangesInitiauxCrees[numMol]=[(nbEtape,[numMol])]
    nouveauxMelanges=[[] for numMol in range(N)] # Tableau tampon pour la boucle


    # II) Calcul des etapes
    while nbEtape<nbEtapeMax:
        nbEtape+=1
        # 1. Reset variables tampon de la boucle
        nouvellesMolecule=[] # Liste des molecules creees à cette etape
        for numMol in range(N): # Liste des melanges initiaux crees à cette etape
            nouveauxMelanges[numMol]=[]

        # 2. Tests pour chaque molecule 
        for reactif1 in moleculesPresentes:
            # Pour chacun des mélanges créés avant
            for melange_initial in melangesInitiauxCrees[reactif1]:
                reactions_avec_1(reactif1,melange_initial,nouvellesMolecule,moleculesPresentes,nouveauxMelanges,melangesInitiaux,melangesInitiauxCrees,nbEtape,listeReactionsBrut,listeReactionParMolecules)
            

        # 3. Màj molécules présente et les mélanges les produisants
        moleculesPresentes=list(set(moleculesPresentes+nouvellesMolecule)) # Ajout des molécules crées à cette étape

        for numMol,liste_melange in enumerate(nouveauxMelanges): # Ajout des mélange initiaux donnat la molécule numMol à cette étape
            for melange in liste_melange:
                melangesInitiaux[numMol].append(melange)
        
        moleculesCreees=[numMol for numMol in nouvellesMolecule]
        melangesInitiauxCrees=[[melange for melange in nouveauxMelanges[numMol]] for numMol in range(N)]
        
        #Fin d'étape

    return melangesInitiaux,melangesInitiaux[moleculeCherche]


# Recherches des nouveaux mélange possible avec reactif1 partant du mélange melange_initial
def reactions_avec_1(reactif1,melange_initial,nouvellesMolecule,moleculesPresentes,nouveauxMelanges,melangesInitiaux,melangesInitiauxCrees,nbEtape,listeReactionsBrut,listeReactionParMolecules):
    for reactionsDe1 in listeReactionParMolecules[reactif1]: # On test les réactions dans lesquelles reactionsDe1 intervient. Faible nombre donc on peut tout tester 

        # On vérifie que les réactifs sont présents
        reactifsPresents = True
        for reactif2 in listeReactionsBrut[reactionsDe1][0]: 
            if not(reactif2 in moleculesPresentes):
                reactifsPresents = False
                break

        #Si les réactifs sont présent ont créer les nouveaux mélanges initiaux des produits
        if reactifsPresents:
            if len(listeReactionsBrut[reactionsDe1][0])==1:
                reactif2=reactif1
            elif len(listeReactionsBrut[reactionsDe1][0])==2:
                if listeReactionsBrut[reactionsDe1][0][0]==reactif1:
                    reactif2=listeReactionsBrut[reactionsDe1][0][1]
                else:
                    reactif2=listeReactionsBrut[reactionsDe1][0][0]
            else:
                print("Plus de 2 réactifs")
                break

            listeMelangesTrouves = unionMelanges(reactif2,melange_initial,melangesInitiaux,nbEtape)
            produitsReaction = listeReactionsBrut[reactionsDe1][1]

            for produit in produitsReaction:
                if (produit not in nouvellesMolecule):
                    nouvellesMolecule.append(produit)

                for melangeTrouve in listeMelangesTrouves:
                    if not(test_melange_existant(produit,melangeTrouve,melangesInitiaux) or (melangeTrouve in nouveauxMelanges[produit])):
                        nouveauxMelanges[produit].append(melangeTrouve)


# Permet d'obtenir les mélanges initiaux entre reactif2 et reactif1 partant du mélange initial melange_initial1
def unionMelanges(reactif2,melange_initial1,melangesInitiaux,nbEtape):
    resultat=[]
    liste_réactifs_initiaux_R1=melange_initial1[1]
    for melangeR2 in melangesInitiaux[reactif2]:
        liste_réactifs_initiaux_R2=melangeR2[1]

        melangePossible=sorted(list(set(liste_réactifs_initiaux_R1+liste_réactifs_initiaux_R2))) # List(Set()) pour supprimer les doublons
        resultat.append((nbEtape,melangePossible))
    
    return resultat


# Pour regarder si un mélange initial a déjà été fait et éviterr des doublons
def test_melange_existant(num_molecule,melangeTest,melangesInitiaux):
    for melange_existant in melangesInitiaux[num_molecule]:
        if melangeTest[1]==melange_existant[1]:
            return True
    return False


# L'objectif de cette fonction est de créer un CRN à partir d'un ensemble d'espèces initiales stocké dans espIni
def creationCRN(espIni,listeReactionsBrut,listeReactionParMolecules) : 
    
    #initialisation des variables nécessaires au fonctionnement
    present = espIni #ensemble des espèces présentes dans le système 
    cree = espIni #ensemble des espèces crées à l'étape précédente
    CRN = [] #ensemble des récations obtenues et correspondant au CRN 
    liste_reaction_passe=[] #[reaction[0] for reaction in CRN]

    #corps de la fonction de création du CRN
    while cree != [] : 
        newEsp = [] #va servir à stocker les espèces crées à cette étape 
        for reactif1 in cree :
            for reactionsDe1 in listeReactionParMolecules[reactif1]:
                # On vérifie que les réactifs sont présents
                reactifsPresents = True
                for reactif2 in listeReactionsBrut[reactionsDe1][0]: 
                    if not(reactif2 in present):
                        reactifsPresents = False
                        break
                
                #Si les réactifs sont présent ont créer les nouveaux mélanges initiaux des produits
                if reactifsPresents:
                    produitsReaction = listeReactionsBrut[reactionsDe1][1]

                    if len(listeReactionsBrut[reactionsDe1][0])==1:
                        CRN.append(([reactif1],produitsReaction))
                    elif len(listeReactionsBrut[reactionsDe1][0])==2:
                        if listeReactionsBrut[reactionsDe1][0][0]==reactif1:
                            reactif2=listeReactionsBrut[reactionsDe1][0][1]
                        else:
                            reactif2=listeReactionsBrut[reactionsDe1][0][0]
                        if(not([reactif1,reactif2] in liste_reaction_passe or [reactif2,reactif1] in liste_reaction_passe)):
                            liste_reaction_passe.append([reactif1,reactif2])
                            CRN.append(([reactif1,reactif2],produitsReaction))

                    else:
                        print("Plus de 2 réactifs")
                        break

                    for produit in produitsReaction: 
                        if not(produit in present or produit in newEsp): 
                            newEsp.append(produit)
        cree = newEsp
        present = present + newEsp
    return CRN 


# L'objectif de cette fonction est de trouver une combinaison initiale d'espèces permettant d'obtenir de CRN de A+B->C
"""Attention il s'agit de la version sans prise en compte des concentrations et de la potentielle disparition de C !!! """
# A et B sont les réactifs et C le produit, listeIni est l'ensemble des espèces initiales permettenat d'obtenir C (cf. la descente)
def remonteeET(A,B,C,listeIni,listeReactionsBrut,listeReactionParMolecules) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 

    #corps de la fonction : 
    for ini in listeIni : 
        nbEtape=ini[0]
        ini=ini[1]
        #Verification que A et B sont présent (condition nécessaire)
        if not(A in ini) or not(B in ini) : 
            None
        else : 
            #Création du CRN sans A
            iniSansA = ini.copy()
            iniSansA.remove(A)
            CRNsansA = creationCRN(iniSansA,listeReactionsBrut,listeReactionParMolecules)

            listeProduits_CRNsansA=set() # pour avoir les liste des molécules présentes dans CRNsansA
            for listeProduits in [reaction[1] for reaction in CRNsansA]: # On récupére la liste des produits de chaque réaction
                for produit in listeProduits: # On ajoute chacun des produits à la liste globale
                    listeProduits_CRNsansA.add(produit)


            #Création du CRN sans B
            iniSansB = ini.copy()
            iniSansB.remove(B)
            CRNsansB = creationCRN(iniSansB,listeReactionsBrut,listeReactionParMolecules) 

            listeProduits_CRNsansB=set() # pour avoir les liste des molécules présentes dans CRNsansB
            for listeProduits in [reaction[1] for reaction in CRNsansB]: # On récupére la liste des produits de chaque réaction
                for produit in listeProduits: # On ajoute chacun des produits à la liste globale
                    listeProduits_CRNsansB.add(produit)


            #Vérifie que C n'apparait pas si A ou B est absent 
            if (C in listeProduits_CRNsansA) or (C in listeProduits_CRNsansB) :
                None

            #Si ça n'est pas le cas on a trouvé les bonnes espèces initiales car si A ou B manque on a pas de C et si A et B sont présents on a C
            #Le dernier point étant dû à l'hypothèse de construction sur listeIni 
            else :
                return (creationCRN(ini,listeReactionsBrut,listeReactionParMolecules),[nbEtape,ini])


# L'objectif de cette fonction est de trouver une combinaison initiale d'espèces permettant d'obtenir le CRN de A ou B -> C
"""Attention il s'agit de la version sans prise en compte des concentrations et de la potentielle disparition de C !!! """
# A et B sont les réactifs et C le produit, listeIni est l'ensemble des espèces initiales permettenat d'obtenir C (cf. la descente)
def remonteeOU(A,B,C,listeIni,listeReactionsBrut,listeReactionParMolecules) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 

    #corpsdef remonteeOU(A,B,C,listeIni) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 
    print("A : ",A)
    print("B : ",B)
    print("C : ",C)
    #corps de la fonction : 
    for ini1 in listeIni : 
        ini1=ini1[1]
        for ini2 in listeIni : 
            ini2=ini2[1]
            #Vérification du fait que A ou B est bien présent dans l'un des deux. 
            if not((A in ini1 and B in ini2) or not(B in ini1 and A in ini2)):  # a checker
                #print("  pas A ou B")
                None
            else : 

                #Création de l'ensemble des espèces initiales sans A ni B pour vérifier que C n'est plus dans le CRN associé
                ini1SansRien = ini1.copy()
                if A in ini1SansRien : 
                    ini1SansRien.remove(A)
                if B in ini1SansRien : 
                    ini1SansRien.remove(B)

                ini2SansRien = ini2.copy()
                if A in ini2SansRien : 
                    ini2SansRien.remove(A)
                if B in ini2SansRien : 
                    ini2SansRien.remove(B)

                iniSansRien = ini1SansRien + ini2SansRien
                #print("iniSansRien ",iniSansRien)
                CRNsansRien = creationCRN(iniSansRien,listeReactionsBrut,listeReactionParMolecules) 

                #Vérifie que C n'apparait pas si A et B sont absents 
                if (C in [triplet[2] for triplet in CRNsansRien]) : 
                    #print("  Non")
                    None

                #Si ça n'est pas le cas on a trouvé les bonnes espèces initiales car si A ou B manque on a pas de C et si A et B sont présents on a C
                #Le dernier point étant dû à l'hypothèse de construction sur listeIni 
                else : 
                    return creationCRN(ini1 + ini2,listeReactionsBrut,listeReactionParMolecules),ini1,ini2 
    

#Cette fonction renvoie le CRN entre A,B et C pour une relation choisie.
#Notez qu'il faut que le CRN existe avec moins de n étapes pour arriver à C. 
#ini correspond aux espèces initialement présentes
def déterminationCRN (A,B,C,n,relation,ini,listeReactionsBrut,listeReactionParMolecules) : 
    if relation == "ET" : 
        melangesInitiaux,melange_C=recherche_melanges_initiaux(ini+[A,B],n,C,listeReactionsBrut,listeReactionParMolecules)
        res=remonteeET(A,B,C,melange_C,listeReactionsBrut,listeReactionParMolecules)
    elif relation == "OU" :
        melangesInitiaux,melange_C=recherche_melanges_initiaux(ini+[A,B],n,C,listeReactionsBrut,listeReactionParMolecules)
        res=remonteeOU(A,B,C,melange_C,listeReactionsBrut,listeReactionParMolecules) 
    else : 
        melangesInitiaux,melange_C=recherche_melanges_initiaux(ini+[A,B],n,C,listeReactionsBrut,listeReactionParMolecules)
        res=remonteeNON(A,B,C,melange_C,listeReactionsBrut,listeReactionParMolecules)
    return melangesInitiaux,res


def remonteeNON() : 
    return None 


def affichage_mélanges_i(melangesInitiaux,numMol,listeMoleculeTexte):
    if len(melangesInitiaux[numMol])==0:
        print(f"Pas de mélange possible pour {listeMoleculeTexte[numMol]}")
        return
    print("Pour obtenir la molécule ",listeMoleculeTexte[numMol])
    for melange in melangesInitiaux[numMol]:
        print(f"   En {melange[0]} étape avec le mélange : ",end="")
        print(",".join("{0}".format(listeMoleculeTexte[mol]) for mol in melange[1]))


def affichage_tous_mélanges(melangesInitiaux,listeMoleculeTexte):
    for mol_i in range(len(listeMoleculeTexte)):
        affichage_mélanges_i(melangesInitiaux,mol_i,listeMoleculeTexte)


def affichage_CRN(melangeInitial,listeReactionsBrut,listeReactionParMolecules,listeMoleculeTexte):
    #melangeInitial liste des molécules présentes au départ
    nb_molécules=len(melangeInitial[1])
    CRN=creationCRN(melangeInitial[1],listeReactionsBrut,listeReactionParMolecules)

    print(f"   En {melangeInitial[0]} étape avec le mélange initial : ",end="")
    print(", ".join("{0}".format(listeMoleculeTexte[mol]) for mol in melangeInitial[1]))
    print()
    listeMoleculesCRN = set()
    print("   Réactions du CRN :")
    for reaction in CRN:
        print(" + ".join("{0}".format(listeMoleculeTexte[mol]) for mol in reaction[0]), '->', " + ".join("{0}".format(listeMoleculeTexte[mol]) for mol in reaction[1]))
        listeMoleculesCRN=listeMoleculesCRN.union(set([mol for mol in reaction[0]]))
        listeMoleculesCRN=listeMoleculesCRN.union(set([mol for mol in reaction[1]]))
    print()
    print("   Liste des molécules du CRN : ")
    print(", ".join("{0}".format(listeMoleculeTexte[mol]) for mol in listeMoleculesCRN))
    print()