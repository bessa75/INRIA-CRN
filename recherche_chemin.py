# ENZ  : Liste des enzymes utilisables (indexé par un numéro du tableau MOL)
# REAC : Liste des entrées (indexé par un numéro du tableau MOL)
# prod : produit à obtenir (indexé par un numéro du tableau MOL)
# etqt : étiquette à avoir ?
# n    : nb max d'étaps ?

def pluscourtchemin(ENZ,REAC,prod,SETM,SETR,n,imprime,liste_reaction_texte,MOL,REACTIONS,REACPARMOL,bool=False,bool2=False):

    """ PRESENCE : liste des molécules présente.
    Evolue au cours de l'algorithme.
    Chaque case coorespond à une molécule et est de la forme :
    [boolean marquant la présence,
    liste des réaction l'ayant produit selon le triplet (numéro de réaction, numéro d'étape, étiquette),
    liste des étiquettes avec lesquelles la molécule a été produite,
    liste des doublets (numéro d'étape, étiquette) avec lesquelles la molécule a été produite]"""
    # print(SETM)
    """     Initialisation des variables """
    # Initialisation de chaque molécule comme absente
    PRESENCE = []
    NbPresence=[]
    for k in range(0, len(MOL)):
        PRESENCE.append([False, [], [], []])
    nbetape = 0  # Nombre d'étapes réactionnels pour obtenir le produis cherché (la plus longue chaine d'étapes)
    nbmol = 0  # Nombre total de molécules présente (au sens qui ont été produite à un moment)
    nbmolbis = 0
    print('REAC :'+str([MOL[i] for i in REAC]))
    ##initialisation de la liste de présence pour les réactifs de REAC (formés à l'étape -1, par la réaction 0 qui n'existe pas)
    for k in range(len(REAC)):
        if k == 0:
            PRESENCE[REAC[k]] = [True, [(-1, 0, 'a')], ['a'], [(-1, 'a')]]
            NbPresence.append((REAC[k],'a'))
        if k == 1:
            PRESENCE[REAC[k]] = [True, [(-1, 0, 'b')], ['b'], [(-1, 'b')]]
            NbPresence.append((REAC[k],'b'))
        nbmol += 1
    print(nbmol)
    for a in ENZ:  ##initialisation de la liste de présence pour ajouter les enzymes
        PRESENCE[a] = [True, [(-1, 0, 'e')], ['e'], [(-1, 'e')]]
        nbmol += 1
        NbPresence.append((a,'e'))

        """ Boucle de recherche descendante """
    ##exploration des différents chemins réactionnels par itérations successives
    while (nbetape < n):
        print('etape :'+str(nbetape))
        print('nombre de molecules :'+str(nbmol))
        nbetape += 1
        nbmolbis = nbmol
        ## on veut que les molécules produites soient notées présentes uniquement à la fin de l'étape pour ne pas mélanger les étapes. On ne met donc pas à jour directement PRESENCE, mais d'abord PRESENCEBIS.
        PRESENCEBIS = []

        ## itération sur les molécules présences
        for doublet in NbPresence:
            num_molecule=doublet[0]
            #print(MOL[num_molecule])
            eti_molecule=doublet[1]
            if PRESENCE[num_molecule][0]:
                # On récupére la liste des réactions où la molécule intervient
                REACPOT = REACPARMOL[num_molecule]

                ## itération sur les réactions impliquant la molécule en tant que réactif
                for a in REACPOT:
                    if a in SETR or bool2:
                        reactifs_presents = True
                        for b in REACTIONS[a][0]:
                            ## on teste si tout les réactifs de la réaction en question sont présents
                            if PRESENCE[b][0] == False:
                                #if 'L-ascorbic acid alpha-D-glucoside' in [MOL[i] for i in REACTIONS[a][0]]:
                                    #print(MOL[b])
                                reactifs_presents = False
                                break
                        if 'L-ascorbic acid alpha-D-glucoside' in [MOL[i] for i in REACTIONS[a][0]]:
                            print([[MOL[i] for i in REACTIONS[a][0]],[MOL[i] for i in REACTIONS[a][1]]])
                            print([doublet for doublet in NbPresence if doublet==(68051, 'ab')])
                        if reactifs_presents:  ##si oui on met à jour la liste de présence
                            ## cas où il y a un seul réactif (marginal) / REACTIONS[a][0] est la liste des réactifs)
                            
                            if len(REACTIONS[a][0]) == 1:
                                ##étape de mise à jour de la liste de présence (molécule produite à telle étape, avec telle étiquette, par telle réaction)
                                for produit in REACTIONS[a][1]:
                                    if produit in SETM or bool2:
                                        if eti_molecule not in PRESENCE[produit][2]:
                                            PRESENCEBIS.append((produit, eti_molecule))
                                        if (a, eti_molecule) not in PRESENCE[produit][3]:
                                            PRESENCE[produit][1].append((a, nbetape, eti_molecule))  ## ici on se permet de mettre à jour PRESENCE et pas PRESENCEBIS car cela n'a pas d'influence.
                                            PRESENCE[produit][3].append((a, eti_molecule))
                                            if PRESENCE[produit][0] == False:
                                                nbmol += 1
    
                            if len(REACTIONS[a][0]) == 2:  ## cas où il y a deux réactifs (cas commun)
                                if REACTIONS[a][0][0]==num_molecule:
                                    mol1 = REACTIONS[a][0][0]
                                    mol2 = REACTIONS[a][0][1]
                                else :
                                    mol2 = REACTIONS[a][0][0]
                                    mol1 = REACTIONS[a][0][1]
                                for e2 in PRESENCE[mol2][2]:
                                    e = bin(eti_molecule,e2)  ##utilisation de la relation binaire pour la propagation des étiquettes
                                    for produit in REACTIONS[a][1]:  ##REACTIONS [a][1] correspond aux produits de la réaction a
                                        if e not in PRESENCE[produit][2]:
                                            PRESENCEBIS.append((produit, e))
                                        if (a, e) not in PRESENCE[produit][3]:  ##contingent si on veut juste le plus court chemin
                                            PRESENCE[produit][3].append((a, e))
                                            PRESENCE[produit][1].append((a, nbetape, e))
                                        if PRESENCE[produit][0] == False:
                                            nbmol += 1
                            if len(REACTIONS[a][0]) == 3:
                                if REACTIONS[a][0][0]==num_molecule:
                                    mol1 = REACTIONS[a][0][0]
                                    mol2 = REACTIONS[a][0][1]
                                    mol3 = REACTIONS[a][0][2]
                                elif REACTIONS[a][0][1]==num_molecule:
                                    mol2 = REACTIONS[a][0][0]
                                    mol1 = REACTIONS[a][0][1]
                                    mol3=REACTIONS[a][0][2]
                                else:
                                    mol2 = REACTIONS[a][0][0]
                                    mol3 = REACTIONS[a][0][1]
                                    mol1=REACTIONS[a][0][2]
                                for e2 in PRESENCE[mol2][2]:
                                    e2b=bin(eti_molecule,e2)
                                    for e3 in PRESENCE[mol3][2]:
                                        e=bin(e2b,e3)
                                        for produit in REACTIONS[a][1]:  ##REACTIONS [a][1] correspond aux produits de la réaction a
                                            if e not in PRESENCE[produit][2]:
                                                PRESENCEBIS.append((produit, e))
                                            if (a, e) not in PRESENCE[produit][3]:  ##contingent si on veut juste le plus court chemin
                                                PRESENCE[produit][3].append((a, e))
                                                PRESENCE[produit][1].append((a, nbetape, e))
                                            if PRESENCE[produit][0] == False:
                                                nbmol += 1
        # mise à jour de PRESENCE à partir de PRESENCE BIS
        NbPresence=[]
        for a in PRESENCEBIS:
            if a[1] not in PRESENCE[a[0]][2]:
                PRESENCE[a[0]][2].append(a[1])
            PRESENCE[a[0]][0] = True
            if a not in NbPresence:
                NbPresence.append(a)
        """ Fin boucle de recherche  """

        """ Boucle de recherche ascendante """
    if bool:
        return(PRESENCE)
    #print("Mécanismes réactionnels obtenus pour le produit en " + str(nbetape) + " étapes maximum avec l'étiquette "+etqt+" :")
    #print(PRESENCE[prod][0:2])
    #print("")
    # On vérifie que le produit voulu a été créé
    if (PRESENCE[prod][0] == False):
        print("impossible d'arriver au produit")
        print("nombre d'étapes=" + str(nbetape))
        print("nombre de molecules=" + str(nbmol))
        print('********************************************')
        print('')
        return ([[[]]])

    MECANISMES = []  ## liste des mécanismes sous forme de triplet (MECANISME, étiquette, enzymes)
    MECANISME = []  ## mécanisme sous forme d'une liste d'étapes, chaque étape étant une liste de réactions
    
    if imprime:
        print("PRESENCE : ",PRESENCE[prod])
        print()

    # une fois le produit trouvé, on remonte la chaîne réactionnelle pour écrire le mécanisme par étapes
    for meca in PRESENCE[prod][1]: #Boucle sur chaque reaction ayant produit la molecule recherche (prod)
        #print("---")
        #print("meca finissant par : ",meca)
        #print()
        
        MECANISME=[]
        presence_cycles=False ## booléen pour éviter que des cycles ne se répètent
        nbetape = meca[1]
        PROD = [(prod,meca[2])] #double (molecule,etiquette)
        
        while nbetape > 0:
            #print()
            #print("PROD : ",PROD)
            # print(nbetape)
            ETAPE = []
            PRODBIS = []
            #print("liste des molecules : ",[MOL[prod[0]] for prod in PROD])
            

            # pour chaque produit on retrouve les réactifs qui l'ont formé et on leur associe l'étiquette correspondante
            for a in PROD:
                # recherche du numéro r de la réaction ayant permis la production de la molécule a en nbetape étapes
                reaction_num = -1
                #print(f"  pour le produit {MOL[a[0]]} avec l'etiquette {a[1]}")
                for reac in PRESENCE[a[0]][1]: #On recherche la reaction ayant donne le couple a=(molecule,etiquette) a l'etape actuel nbetape
                    if reac[1] == nbetape and reac[2] == a[1] and reac[1] == min([reac2[1] for reac2 in PRESENCE[a[0]][1] if reac2[2] == a[1]]):
                        reaction_num = reac[0]
                        break
                
                if reaction_num == -1: ## potentiellement contingent
                    #print(nbetape,MOL[a[0]],a[1])
                    PRODBIS.append(a)
                    
                else:
                    reactifs,produits = REACTIONS[reaction_num][0],REACTIONS[reaction_num][1]
                    
                    #print("    pour la reaction \"",end="")
                    #print(' + '.join([MOL[molecule_i] for molecule_i in reactifs]),'->',' + '.join([MOL[molecule_i] for molecule_i in produits]),end="")
                    #print(f"\" a l'etape {reac[1]} avec l'etiquette {a[1]}={reac[2]}")
                    if reaction_num not in ETAPE:
                        ETAPE.append(reaction_num)
                    for p in produits:
                        if p == prod and nbetape < meca[1]-1:
                            presence_cycles=True #détection d'un cycle rendant le mécanisme invalide : un mécanisme invalide est un mécanisme dans le quel le produit final n'apparait pas seulement dans la dernière réaction
                            break ##on veut sortir de 3 boucles imbriquées donc on va rappeler ce break
                    
                    if presence_cycles: break
                    
                    if len(reactifs) == 1:
                        PRODBIS.append((reactifs[0], a[1]))
                    else:  ##on suppose le nombre de reactifs inférieure ou égale à 2
                        determination_etiquette(reactifs,a[1],[PRESENCE[reactif][1] for reactif in reactifs],nbetape,PRODBIS)
                    
                if presence_cycles: break ##break prolongeant un autre pour sortir de 3 boucles successives
            if presence_cycles: break ##break prolongeant un autre pour sortir de 3 boucles successives
            # fin de la boucle pour remonter les reactifs
            
            if ETAPE != []:
                MECANISME = [ETAPE] + MECANISME ##interet des crochets autour d'ETAPE ?
            # print([MOL[a[0]] for a in PROD])
            # print([MOL[a[0]] for a in PRODBIS])
            PROD = [a for a in PRODBIS]
            nbetape=nbetape-1
        
        #print(presence_cycles)
        if presence_cycles==False:
            Enzs = []
            for mol in PROD:
                if mol[0] != REAC[0] and mol[0] != REAC[1]:
                    Enzs.append(mol[0])
            MECANISMES.append((MECANISME,meca[2],Enzs))




    print('Les mécanismes après sélection sont les suivants')
    print(MECANISMES)
    print(" ")
    if imprime==True: ##Si on décide de print les mécanismes sous forme de texte, on les affiche
        i=0
        for meca in MECANISMES:
            i+=1
            print("mécanisme "+str(i)+":")
            txt=mecatexte(meca[0],liste_reaction_texte)
            for ligne in txt:
                print(' --- '.join(ligne))
            print(" ")
    print('********************************************')
    return MECANISMES


def mecatexte(MECANISME,liste_reaction_texte):  ##simple fonction qui convertit le mécanisme réactionnel renvoyé par l'algorithme en un texte lisible pour l'utilisateur.
    # MECANISME : liste des etapes reactionelles
    # Chaque etape (case) dans MECANISME contient
    MT = []
    for k in range(len(MECANISME)):
        ET = []
        for a in MECANISME[k]:
            reactifs,produits = liste_reaction_texte[a]
            ET.append(' + '.join(reactifs)+' -> '+' + '.join(produits))
        MT.append(ET)
    return (MT)


def determination_etiquette(R,eti,P,nbetape,PRODBIS):
    """
    r_i les deux reactifs dont il faut determiner les etiquettes lors de leur production
    eti l'etiquette obtenue par la reaction entre les deux reactifs
    P_i réactions ayant produit r_i
    PRODBIS le tableau des (molecules,etiquette) a mettre a jour
    """
    
    # Création de L0 et L1, liste des étiquettes avec lesquels r_0 a pu être formé avant nbetape
    r_0=R[0]
    r_1=R[1]
    P_0=P[0]
    P_1=P[1]
    L0 = [(d[1],d[2]) for d in P_0 if d[1] <= nbetape - 1]
    L1 = [(d[1],d[2]) for d in P_1 if d[1] <= nbetape - 1]
    # print(L0)
    # print(L1)
    if eti == 'e':
        for r in R:
            PRODBIS.append((r,'e'))
        
    if eti == 'a':
        for i in range (0,len(R)):
            r=R[i]
            p=P[i]
            L=[d[2] for d in p if d[1] <= nbetape - 1]
            if 'a' in L:
                PRODBIS.append((r, 'a'))
            else:
                PRODBIS.append((r, 'e'))
            
    if eti == 'b':
        for i in range (0,len(R)):
            r=R[i]
            p=P[i]
            L=[d[2] for d in p if d[1] <= nbetape - 1]
            if 'b' in L:
                PRODBIS.append((r, 'b'))
            else:
                PRODBIS.append((r, 'e'))
        '''
    if eti == 'ab' and len(R)==2:
        sel = selec2ab(L0, L1)
        if sel[0] != 'o':
            PRODBIS.append((r_0, sel[0]))
        if sel[1] != 'o':
            PRODBIS.append((r_1, sel[1]))
            
        if sel[0] == 'o' and sel[1] != 'o':
            for reaction in P_0:
                if (reaction[2] == 'a') or (reaction[2] == 'b'):
                    PRODBIS.append((r_0, reaction[2]))
                    break
                    
        if sel[1] == 'o' and sel[0] != 'o':
            for reaction in P_1:
                if (reaction[2] == 'a') or (reaction[2] == 'b'):
                    PRODBIS.append((r_1, reaction[2]))
                    break
                    
        if sel[0] == 'o' and sel[1] == 'o':
            for reaction in P_0:
                if (reaction[2] == 'a') or (reaction[2] == 'b'):
                    PRODBIS.append((r_0, reaction[2]))
                    etq = reaction[2]
                    break
            for reaction in P_1:
                if ((reaction[2] == 'a') or (reaction[2] == 'b')) and reaction[2] != etq:
                    PRODBIS.append((r_1, reaction[2]))
'''
    LL=[]
    for Pi in P:
        Pb=[(d[1],d[2]) for d in Pi if d[1] <= nbetape - 1]
        LL.append(Pb)
    if eti=='ab':
        sel=selecXab(LL)
        for i in range (0,len(R)):
            PRODBIS.append((R[i],sel[i]))

def bin(c, d):  ##relation binaire de propagation des étiquettes
    if c == 'a' and d == 'b':
        return('ab')
    if c == 'b' and d == 'a':
        return('ab')
    if c == 'ab' or d == 'ab':
        return('ab')
    if c == 'a' or d == 'a':
        return('a')
    if c == 'b' or d == 'b':
        return('b')
    else:
        return('e')


def selec(L):  ##fonction utile pour selec2ab
    L1=[l[1] for l in L]
    if 'ab' in L1:
        return('ab')
    if 'a' in L1 and 'b' in L1:
        mina=1000
        minb=1000
        for i in range (0,len(L)):
            if L1[i]=='a':
                mina=min([mina,L[i][0]])
            if L1[i]=='b':
                minb=min([minb,L[i][0]])
        if mina<minb:
            return('oa')
        else:
            return('ob')
    if 'b' in L1:
        return('b')
    if 'a' in L1:
        return('a')
    else:
        return('e')


def selec2ab(L1,L2):  ##cette fonction aide l'algorithme de remontée. Le but est, quand le produit d'une réaction est avec l'étiquette 'ab', de choisir avec quelles étiquettes on va considérer les réactifs qui ont mené à ce produit. Si on a 'o' on peut choisir 'a' ou 'b' indifféremment.
    l1 = selec(L1)
    l2 = selec(L2)
    if l2[0] == 'o':
        if l1 == 'ab':
            return('ab', l2[1])
        if l1 == 'b':
            return('b', 'a')
        if l1 == 'a':
            return('a', 'b')
        if l1[0] == 'o':
            if l1[1]==l2[1]:
                return('a', 'b')
            else:
                return(l1[1],l2[1])
    if l1[0] == 'o':
        if l2[0] == 'o':
            if l1[1]==l2[1]:
                return('a', 'b')
            else:
                return(l1[1],l2[1])
        if l2 == 'b':
            return('a', 'b')
        if l2 == 'a':
            return('b', 'a')
        if l2 == 'ab':
            return(l1[1], 'ab')
    return(l1, l2)

def selecXab(LL):
    SEL=[]
    ETI=[]
    BOOL=[False,False,False]
    CORR=['ab','a','b']
    for i in range (0,len(LL)):
        sel=selec(LL[i])
        SEL.append(sel)
        ETI.append(sel)
        if sel=='ab':
            BOOL[0]=True
        if sel=='a':
            BOOL[1]=True
        if sel=='b':
            BOOL[2]=True
    if BOOL[0]:
        for eti in ETI:
            if eti[0]=='o':
                eti=eti[1]
        return(ETI)
    if BOOL[1]:
        for eti in ETI:
            if eti=='ob':
                eti='b' 
                BOOL[2]=True
        if BOOL[2]:
            for eti in ETI:
                if eti=='oa':
                    eti='a'
        else:
            b=True
            for eti in ETI:
                if eti=='oa' and b:
                    b=False
                    eti='b'
                elif eti=='oa' and b==False:
                    eti='a'
        return(ETI)
    if BOOL[2]:
        for eti in ETI:
            if eti=='oa':
                eti='a' 
                BOOL[1]=True
        if BOOL[1]:
            for eti in ETI:
                if eti=='ob':
                    eti='b'
        else:
            b=True
            for eti in ETI:
                if eti=='ob' and b:
                    b=False
                    eti='a'
                elif eti=='ob' and b==False:
                    eti='b'
        return(ETI)
    
    
    
def numero(nom_molecule,MOL):  ##retourne le numéro correspondant à un nom d'enzyme
    for i in range(0, len(MOL)):
        if MOL[i] == nom_molecule:
            return (i)


def numero2(L,MOL):  ##retourne les numéros correspondant aux noms d'enzymes d'une liste d'enzymes
    L2 = []
    for k in range(0, len(L)):
        L2.append(numero(L[k],MOL))
    return (L2)


# Prend un numéro de molécule en paramétre et renvoie le numéro de molécule faisant la logique "non" ainsi que le numéro de réaction
def molécule_non_v1(numéro_molécule,MOL,REACTIONS):
    if(MOL[numéro_molécule].endswith('ext')): #Il faut prendre la molécule une fois dans la cellule pour chercher un non
        nom_molécule=MOL[numéro_molécule]
        if MOL[numéro_molécule+1]==nom_molécule[:-3]:
            numéro_molécule=numéro_molécule+1
        else:
            for num_molécule_boucle,nom_molécule_boucle in enumerate(MOL):
                if nom_molécule_boucle==nom_molécule[:-3]:
                    numéro_molécule=num_molécule_boucle
                    break
        

    for num_reac, réacion_i in enumerate(REACTIONS): # La molécule une fois rentré ne peut pas sortir, donc pas besoin d'exclure la réaction retour (sauf s'il la réaction A=>Aext devient possible)
        if numéro_molécule in réacion_i[0]:
            if réacion_i[0][0] == numéro_molécule:
                molécule_non = réacion_i[1][0]
            else:
                molécule_non = réacion_i[0][0]
            
            print(f"Molécule non : {MOL[molécule_non]}")
            réaction_négation=REACTIONS[num_reac]
            print(f"Réaction de négation : {' + '.join([MOL[indice_réactif] for indice_réactif in réaction_négation[0]])} => {MOL[réaction_négation[1][0]]}")

            return molécule_non, num_reac
    raise Exception("Aucune molécule 'non' trouvée")


# ET : '+'
# OU : '|'
# NON : '!'
# OU EXCLUSIF : '-'
# symbole réaction : '=>'
def résolution_équation(ENZ,équation_logique,nb_réactions_max,liste_reaction_texte,MOL,REACTIONS,REACPARMOL):
    liste_mots = équation_logique.split(" ")

    # Récupération du numéro du 1er réactif
    if liste_mots[0][0] == '!':
        réactif_1, num_reac_non1 = molécule_non_v1(numero(liste_mots[0][1:],MOL),MOL,REACTIONS)
    else:
        réactif_1 = numero(liste_mots[0],MOL)

    # Récupération du numéro du 2eme réactif
    if liste_mots[2][0] == '!':
        réactif_2, num_reac_non2 = molécule_non_v1(numero(liste_mots[2][1:],MOL),MOL,REACTIONS)
    else:
        réactif_2 = numero(liste_mots[2],MOL)

    produit = numero(liste_mots[4],MOL)

    if liste_mots[1] == '+':
        MECAS = pluscourtchemin(numero2(ENZ,MOL), [réactif_1, réactif_2], produit,set(),set(), nb_réactions_max, True,liste_reaction_texte,MOL,REACTIONS,REACPARMOL)
        mt = mecatexte(MECAS[0][0],liste_reaction_texte)
        for d in mt:
            print(d)
        print("")
        print("")
        #return MECAS


# def afficher (mecanisme):
#     for i in range (0,len(mecanisme)):
#         print(mecanisme[i])
#         print("")