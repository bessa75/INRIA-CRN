#test push3 commit without another branch
MOLTEST=["molecule1"]
MOL=["molecule1","molecule2","molecule3","molecule4","molecule5","molecule6","enzyme1","enzyme2","enzyme3"]## liste des molécules associées à leur numéro


REACTIONS=[([0],[0]),([2,5],[3,4]),([1,7],[4,6]),([7,3],[2,8]),([1,6],[4,7])] ## liste des réactions sous forme de doublet avec les réactifs à gauche et les produits à droite


REACPARMOL=[[],[2,4],[1],[3],[],[1],[4],[2],[3]] ## pour chaque molécule, liste des réactions où elle est un réactif

PRESENCE=[[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]]]

def pluscourtchemin(REAC,prod):
    PRESENCE=[[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]],[False,[]]]
    nbetape=0
    nbmol=0
    nbmolbis=0
    for a in REAC: ##initialisation
        PRESENCE[a]=[True,[(-1,0)]]
        nbmol+=1
    while (PRESENCE[prod][0]==False) and (nbmolbis!=nbmol): ##exploration des différents chemins réactionnels
        nbetape+=1
        nbmolbis=nbmol
        PRESENCEBIS=[]
        for i in range (len(PRESENCE)):
            if PRESENCE[i][0]:
                REACPOT=REACPARMOL[i]
                for a in REACPOT:
                    bool=True
                    for b in REACTIONS[a][0]:
                        if PRESENCE[b][0]==False:
                            bool=False
                    if bool:
                        for c in REACTIONS[a][1]:
                            if PRESENCE[c][0]==False:
                                PRESENCEBIS.append(c)
                                nbmol+=1
                            if notin(PRESENCE[c][1],(a,nbetape)):
                                PRESENCE[c][1].append((a,nbetape)) ##problème : doublon
        for a in PRESENCEBIS:
            PRESENCE[a][0]=True
    if (nbmolbis==nbmol) and (PRESENCE[prod][0]==False):
        return(False)
    MECANISME=[]
    PROD=[prod]
    while nbetape>0: ##une fois le produit trouvé, on remonte la chaîne réactionnelle pour écrire le mécanisme par étapes
        ETAPE=[]
        PRODBIS=[]
        for a in PROD:
            if PRESENCE[a][1][0][1]==nbetape:
                reac=PRESENCE[a][1][0][0]
                if notin(ETAPE,reac):
                    ETAPE.append(reac)
                for b in REACTIONS[reac][0]:
                    if notin(PRODBIS,b):
                        PRODBIS.append(b)
            else:
                PRODBIS.append(a)
        MECANISME=[ETAPE]+MECANISME
        print(PROD)
        PROD=PRODBIS
        PRODBIS=[]
        nbetape-=1
    print(PROD)
    return(MECANISME)

def notin(L,a):
    for b in L:
        if b==a:
            return False
    return True


