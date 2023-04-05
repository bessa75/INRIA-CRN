def recherche_melanges_initiaux(EspecesInitiales,nbEtapeMax,moleculeCherche,listeReactions):
    N=len(listeReactions)

    # I) Initialisation
    nbEtape=0 # Numero de l'etape en cours
    
    MoleculesPresentes=[num_mol for num_mol in EspecesInitiales] # Liste de toutes les molecules crees depuis le debut
    MoleculesCreees=[num_mol for num_mol in EspecesInitiales] # Liste des molecules cree à l'etape precedente
    NouvellesMolecule=[] # Tableau tampon pour la boucle

    MelangesInitiaux=[[] for num_mol in range(N)] # Tableau où chaque case i est la liste des doublets (n, combinaisons d'especes initiales permettant d'obtenir la molecule i en n etape) 
    MelangesInitiauxCrees=[[] for num_mol in range(N)] # Idem MelangesInitiaux mais que ceux creees à l'etape precedente
    for num_mol in EspecesInitiales:
        MelangesInitiaux[num_mol]=[(nbEtape,[num_mol])]
        MelangesInitiauxCrees[num_mol]=[(nbEtape,[num_mol])]
    NouveauxMelanges=[[] for num_mol in range(N)] # Tableau tampon pour la boucle


    # II) Calcul des etapes
    while nbEtape<nbEtapeMax:
        nbEtape+=1
        # 1. Reset variables tampon de la boucle
        NouvellesMolecule=[] # Liste des molecules creees à cette etape
        for num_mol in range(N): # Liste des melanges initiaux crees à cette etape
            NouveauxMelanges[num_mol]=[]

        # 2. Tests pour chaque molecule 
        for reactif1 in MoleculesPresentes:
            # Pour chacun des mélanges créés avant
            for melange_initial in MelangesInitiauxCrees[reactif1]:
                reactions_avec_1(reactif1,melange_initial,NouvellesMolecule,MoleculesPresentes,NouveauxMelanges,MelangesInitiaux,MelangesInitiauxCrees,nbEtape,listeReactions)
            

        # 3. Màj molécules présente et les mélanges les produisants
        MoleculesPresentes=list(set(MoleculesPresentes+NouvellesMolecule)) # Ajout des molécules crées à cette étape

        for num_mol,liste_melange in enumerate(NouveauxMelanges): # Ajout des mélange initiaux donnat la molécule num_mol à cette étape
            for melange in liste_melange:
                MelangesInitiaux[num_mol].append(melange)
        
        MoleculesCreees=[num_mol for num_mol in NouvellesMolecule]
        MelangesInitiauxCrees=[[melange for melange in NouveauxMelanges[num_mol]] for num_mol in range(N)]
        
        #Fin d'étape

    return MelangesInitiaux,MelangesInitiaux[moleculeCherche]

# Recherches des nouveaux mélange possible avec reactif1 partant du mélange melange_initial
def reactions_avec_1(reactif1,melange_initial,NouvellesMolecule,MoleculesPresentes,NouveauxMelanges,MelangesInitiaux,MelangesInitiauxCrees,nbEtape,listeReactions):
    for reactif2 in MoleculesPresentes:
        produits_reaction=listeReactions[reactif1][reactif2]
        if produits_reaction!=[]:
            listeMelangesTrouves=unionMelanges(reactif1,reactif2,melange_initial,MelangesInitiaux,MelangesInitiauxCrees,nbEtape)

            for produit in produits_reaction:
                if (produit not in NouvellesMolecule):
                    NouvellesMolecule.append(produit)

                for melangeTrouve in listeMelangesTrouves:
                    if not(test_melange_existant(produit,melangeTrouve,MelangesInitiaux) or (melangeTrouve in NouveauxMelanges[produit])):
                        NouveauxMelanges[produit].append(melangeTrouve)


def unionMelanges(reactif1,reactif2,melange_initial1,MelangesInitiaux,MelangesInitiauxCrees,nbEtape):
    resultat=[]
    liste_réactifs_initiaux_R1=melange_initial1[1]
    for melangeR2 in MelangesInitiaux[reactif2]:
        liste_réactifs_initiaux_R2=melangeR2[1]

        melangePossible=sorted(list(set(liste_réactifs_initiaux_R1+liste_réactifs_initiaux_R2))) # List(Set()) pour supprimer les doublons
        resultat.append((nbEtape,melangePossible))
    
    return resultat


def test_melange_existant(num_molecule,melangeTest,MelangesInitiaux):
    for melange_existant in MelangesInitiaux[num_molecule]:
        if melangeTest[1]==melange_existant[1]:
            return True
    return False


#L'objectif de cette fonction est de créer un CRN à partir d'un ensemble d'espèces initiales stocké dans espIni
def creationCRN(espIni,R) : 
    
    #initialisation des variables nécessaires au fonctionnement
    present = espIni #ensemble des espèces présentes dans le système 
    cree = espIni #ensemble des espèces crées à l'étape précédente
    CRN = [] #ensemble des récations obtenues et correspondant au CRN 

    #corps de la fonction de création du CRN
    while cree != [] : 
        newEsp = [] #va servir à stocker les espèces crées à cette étape 
        for mol1 in cree : 
            for mol2 in present : 
                for mol3 in R[mol1-1][mol2-1]: 
                    if  not(((mol1,mol2,mol3) in CRN)or (mol2,mol1,mol3) in CRN)  : 
                        newEsp.append(mol3)
                        CRN.append((mol1,mol2,mol3))
        cree = newEsp
        present = present + newEsp

    return CRN 


#L'objectif de cette fonction est de trouver une combinaison initiale d'espèces permettant d'obtenir de CRN de A+B->C
"""Attention il s'agit de la version sans prise en compte des concentrations et de la potentielle disparition de C !!! """
#A et B sont les réactifs et C le produit, listeIni est l'ensemble des espèces initiales permettenat d'obtenir C (cf. la descente)
def remonteeET(A,B,C,listeIni,R) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 

    #corps de la fonction : 
    for ini in listeIni : 
        ini=ini[1]
        #Verification que A et B sont présent (condition nécessaire)
        if not(A in ini) or not(B in ini) : 
            None
        else : 
            #Création du CRN sans A
            iniSansA = ini.copy()
            iniSansA.remove(A)
            CRNsansA = creationCRN(iniSansA,R) 
            #print(CRNsansA)
            #Création du CRN sans B
            iniSansB = ini.copy()
            iniSansB.remove(B)
            CRNsansB = creationCRN(iniSansB,R) 
            #print(CRNsansB)

            #Vérifie que C n'apparait pas si A ou B est absent 
            if (C in [triplet[2] for triplet in CRNsansA]) or (C in [triplet[2] for triplet in CRNsansA]) : 
                None

            #Si ça n'est pas le cas on a trouvé les bonnes espèces initiales car si A ou B manque on a pas de C et si A et B sont présents on a C
            #Le dernier point étant dû à l'hypothèse de construction sur listeIni 
            else : 
                #print(ini)
                return creationCRN(ini,R),ini


#L'objectif de cette fonction est de trouver une combinaison initiale d'espèces permettant d'obtenir le CRN de A ou B -> C
"""Attention il s'agit de la version sans prise en compte des concentrations et de la potentielle disparition de C !!! """
#A et B sont les réactifs et C le produit, listeIni est l'ensemble des espèces initiales permettenat d'obtenir C (cf. la descente)
def remonteeOU(A,B,C,listeIni) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 

    #corpsdef remonteeOU(A,B,C,listeIni) : 

    #on peut ajouter une initialisation que consisterait à trier listeIni selon un ordre défénie 

    #corps de la fonction : 
    for ini1 in listeIni : 
        for ini2 in listeIni : 
            #Vérification du fait que A ou B est bien présent dans l'un des deux. 
            if not((A in ini1 and B in ini2) or not(B in ini1 and A in ini2)):  
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
                
                CRNsansRien = creationCRN(iniSansRien) 

                #Vérifie que C n'apparait pas si A et B sont absents 
                if (C in [triplet[2] for triplet in CRNsansRien]) : 
                    None

                #Si ça n'est pas le cas on a trouvé les bonnes espèces initiales car si A ou B manque on a pas de C et si A et B sont présents on a C
                #Le dernier point étant dû à l'hypothèse de construction sur listeIni 
                else : 
                    return creationCRN(ini1 + ini2),ini1,ini2 
    

#Cette fonction renvoie le CRN entre A,B et C pour une relation choisie.
#Notez qu'il faut que le CRN existe avec moins de n étapes pour arriver à C. 
#ini correspond aux espèces initialement présentes

def déterminationCRN (A,B,C,n,relation,ini,listeReactions) : 
    if relation == "ET" : 
        MelangesInitiaux,melange_C=recherche_melanges_initiaux(ini,n,C,listeReactions)
        res=remonteeET(A,B,C,melange_C,listeReactions)
    elif relation == "OU" :
        MelangesInitiaux,melange_C=recherche_melanges_initiaux(ini,n,C,listeReactions)
        res=remonteeOU(A,B,C,melange_C,listeReactions) 
    else : 
        MelangesInitiaux,melange_C=recherche_melanges_initiaux(ini,n,C,listeReactions)
        res=remonteeNON(A,B,C,melange_C,listeReactions)
    return res


def remonteeNON() : 
    return None 



def affichage_mélanges_i(MelangesInitiaux,numMol,listeMoleculeTexte):
    if len(MelangesInitiaux[numMol])==0:
        print(f"Pas de mélange possible pour {listeMoleculeTexte[numMol]}")
        return
    print("Pour obtenir la molécule ",listeMoleculeTexte[numMol])
    for melange in MelangesInitiaux[numMol]:
        print(f"   En {melange[0]} étape avec le mélange : ",end="")
        for réac_initiaux in melange[1]:
            print(listeMoleculeTexte[réac_initiaux],end=",")
        print()


def affichage_tous_mélanges(MelangesInitiaux,listeMoleculeTexte):
    for mol_i in range(len(listeMoleculeTexte)):
        affichage_mélanges_i(MelangesInitiaux,mol_i,listeMoleculeTexte)
