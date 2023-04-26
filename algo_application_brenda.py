
import recherche_chemin


def arborescence_produit(prod,n,MOL,REACTIONS,REACPARMOLP):
    # Initialisation de chaque molécule comme absente
    PRESENCE = []
    REAC=[]
    NbPresence=[]
    for k in range(0, len(MOL)):
        PRESENCE.append(False)
    for k in range(0, len(REACTIONS)):
        REAC.append(False)

    ##initialisation de la liste de présence pour les réactifs de REAC (formés à l'étape -1, par la réaction 0 qui n'existe pas)
    PRESENCE[prod]=True
    NbPresence.append(prod)
    NbPresencebis=NbPresence
    NbREAC=[]
    NbPresenceTot=NbPresence
    nbetape=0

    """ Boucle de recherche descendante """
    ##exploration des différents chemins réactionnels par itérations successives
    while (nbetape < n):
        nbetape += 1
        ## on veut que les molécules produites soient notées présentes uniquement à la fin de l'étape pour ne pas mélanger les étapes. On ne met donc pas à jour directement PRESENCE, mais d'abord PRESENCEBIS.
        NbPresence=NbPresencebis
        NbPresencebis = []

        ## itération sur les molécules présences
        for num_molecule in NbPresence:
            # On récupére la liste des réactions où la molécule intervient
            REACPOT = REACPARMOLP[num_molecule]
            for r in REACPOT:
                if REAC[r]==False:
                    REAC[r]=True
                    NbREAC.append(r)
                    reac=REACTIONS[r]
                    for re in reac[0]:
                        if PRESENCE[re]==False:
                            NbPresencebis.append(re)
                            PRESENCE[re]=True
        NbPresenceTot=NbPresenceTot+NbPresencebis

    return(NbPresenceTot,NbREAC)

def arborescence_reactif(ENZ,Reac,n,MOL,REACTIONS,REACPARMOL):
    # Initialisation de chaque molécule comme absente
    PRESENCE = []
    REAC=[]
    NbPresenceTot=[]
    #print(REAC)
    for k in range(0, len(MOL)):
        PRESENCE.append(False)
    for k in range(0, len(REACTIONS)):
        REAC.append(False)
    for r in Reac:
        PRESENCE[r]=True
        NbPresenceTot.append(r)
    for e in ENZ:
        PRESENCE[e]=True
        NbPresenceTot.append(e)
    ##initialisation de la liste de présence pour les réactifs de REAC (formés à l'étape -1, par la réaction 0 qui n'existe pas)
    NbREAC=[]
    nbetape=0
    NbPresencebis=NbPresenceTot
    print('***********')
    """ Boucle de recherche descendante """
    ##exploration des différents chemins réactionnels par itérations successives
    while (nbetape < n):
        nbetape += 1
        ## on veut que les molécules produites soient notées présentes uniquement à la fin de l'étape pour ne pas mélanger les étapes. On ne met donc pas à jour directement PRESENCE, mais d'abord PRESENCEBIS.
        NbPresence=NbPresencebis
        NbPresencebis = []
        ## itération sur les molécules présences
        for num_molecule in NbPresence:
            # On récupére la liste des réactions où la molécule intervient
            REACPOT = REACPARMOL[num_molecule]
            for r in REACPOT:
                rreac=REACTIONS[r]
                rreactifs=rreac[0]
                rproduits=rreac[1]
                bool=True
                for r in rreactifs:
                    if PRESENCE[r]==False:
                        bool=False
                if bool and REAC[r]==False:
                    NbREAC.append(r)
                    REAC[r]=True
                    #print(r)
                    #print(rreactifs,rproduits)
                    #print([MOL[r] for r in rreactifs],[MOL[r] for r in rproduits])
                    for p in rproduits:
                        if PRESENCE[p]==False:
                            PRESENCE[p]=True
                            NbPresencebis.append(p)
        NbPresenceTot+=NbPresencebis
    return(NbPresenceTot,NbREAC)

def brenda(ENZ,Reac,prod,n,imprime,liste_reaction_texte,MOL,REACTIONS,REACPARMOL,REACPARMOLP,bool=False):
    (ArbMR,ArbReacR)=arborescence_reactif(ENZ,Reac,n,MOL,REACTIONS,REACPARMOL)
    (ArbMP,ArbReacP)=arborescence_produit(prod,n,MOL,REACTIONS,REACPARMOLP)
    Reac=ArbReacR+ArbReacP
    SETM=set()
    for m in ArbMR:
        SETM.add(m)
    for m in ArbMP:
        SETM.add(m)
    SETR=set()
    for r in ArbReacR:
        SETR.add(r)
    for r in ArbReacP:
        SETR.add(r)
    return(recherche_chemin.pluscourtchemin(ENZ,Reac,prod,SETM,SETR,n,imprime,liste_reaction_texte,MOL,REACTIONS,REACPARMOL,bool=False))