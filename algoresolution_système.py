from recherche_chemin import pluscourtchemin
import itertools
import algonegation

def check(Current,LETIQ,nmax,LREAC,LPROD,reaction,MOL,REACTIONS,REACPARMOL,n,reac,CYCLES,CYCLESPARMOL): #permet de vérifier que les différents chemins choisis n'interfèrent pas sur leurs étiquettes choisies
    ENZYME=[]
    n=len(LPROD)
    for ligne in Current:
        enzliste=ligne[2]
        for enz in enzliste:
            ENZYME.append(enz)  #on agrège toutes les enzymes pour obtenir l'ensemble des enzymes qu'on introduit dans le milieu. Il reste à retirer les doublons pour rendre plus efficace les recherches de chemins suivantes
    for k in range(n):
        chemins=pluscourtchemin(ENZYME,LREAC[k],LPROD[k],nmax,False,reaction,MOL,REACTIONS,REACPARMOL) # on lance plus court chemin avec l'ensemble des enzymes pour voir si les étiquettes sont respectées, dès que ce n'est plus vérifié on renvoie false
        if LETIQ[k]=='ab':  #si on cherche un a et b, alors tout chemin avec seulement a ou seulement b entraine une invalidité de la fonction logique donc l'étiquette n'est plus vérifiée
            for chemin in chemins:
                if chemin[1]!='ab':
                    return False
        if LETIQ[k]=='aob': #il faut trouver un cehmin avec seulement a et un avec seulement b pour vérifier l'étiquette, peu importe si il y a un et (on pourrait implémenter le ou exclusif en vérifiant qu'il n'y a pas de et logique)
            inda=False
            indb=False
            for chemin in chemins:
                if chemin[1]=='a':
                    inda=True
                if chemin[1]=='b':
                    indb=True
                if chemin[1]=='e':
                    return False
            if not inda or not indb:
                return False
        if LETIQ[k]=='anb': #on a déjà créé une ligne avec la negation, on a donc séparé le A+nonB-->C en A-->C et une autre ligne non B-->C, ce qui correspond à un 'et' entre les deux lignes
            inda=False
            for chemin in chemins:
                if chemin[1]=='a':
                    inda=True
                if chemin[1]=='e':
                    return False
            if not inda:
                return False
        if LETIQ[k]=='a': #il faut un chemin 'a' et le seul chemin interdit est tout chemin 'e'
            inda=False
            for chemin in chemins:
                if chemin[1]=='a':
                    inda=True
                if chemin[1]=='e':
                    return False
            if not inda:
                return False
    return True # si à la fin du parcours on a pas trouvé d'invalidation d'un ligne, le système fonctionne et on renvoie true


def research(Allreactions,LETIQ,nmax,LREAC,LPROD,reaction,MOL,REACTIONS,REACPARMOL,n,reac,CYCLES,CYCLESPARMOL): # permet de rechercher si il existe une combinaison des lots d'enzyme d'un chemin par ligne des différents chemins trouvés pour chaque ligne du système tel que les étiquettes de chaque ligne sont préservées
    n=len(Allreactions)
    i=[0]*(n+1)
    Borne = [0]*n+[1]
    ALL=list(itertools.product(*Allreactions))
    for Current in ALL:
        if (check(Current,LETIQ,nmax,LREAC,LPROD,reaction,MOL,REACTIONS,REACPARMOL,n,reac,CYCLES,CYCLESPARMOL)): #on fait appel à la fonction check pour vérifier que les étiquettes sont bien préservées sur chaque ligne
            return Current
    return 0 #on renvoie 0 pour signifier qu'aucune combinaison de fonctionne

def res(LREAC,LPROD,LETIQ,nmax,ENZ,reaction,MOL,REACTIONS,REACPARMOL,reac,CYCLES,CYCLESPARMOL):  #programme à apeler pour la résolution du système, avec des listes pour les produits, réactifs et étiquettes où l'indice i correspond à la ième ligne du système, 'ab' pour et, 'aob' pour a ou b, 'a' pour une simple présence, 'anb' pour a et non b, 'aonb' pour a ou non b, les négations restent à implémenter
    n=len(LPROD)
    nneg=0
    for k in range(n):
        if LETIQ[k]=='anb': #on créée des lignes suplémentaires pour traiter les négations sous la forme de deux lignes, cf explication dans check
            nneg+=1
            LREAC.append([LREAC[k][1]])
            LPROD.append(LPROD[k])
            LETIQ.append('nb')
    ENZYME=[]
    MECANISMES=[]
    Allreactions=[] # tableau de tous les chemins trouvés pour chaque ligne du système de manière indépendante
    for k in range(n):
        Allreactions.append(pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False,reaction,MOL,REACTIONS,REACPARMOL)) # on remplit le tableau pour chaque ligne du système
    for k in range(n,nneg):
        Allreactions.append(algonegation.algo_negation(0,LREAC[k],LPROD[k],CYCLES,CYCLESPARMOL,REACPARMOL,reaction,MOL,REACTIONS,nmax,reac)) # on remplit le tableau pour chaque ligne du système
    result=research(Allreactions,LETIQ,nmax, LREAC,LPROD,reaction,MOL,REACTIONS,REACPARMOL,n,reac,CYCLES,CYCLESPARMOL) # on rentre dans le programme secondaire (qui pourrait être à la suite sans problème)
    if result==0:
        return ("pas de solution au système")
    else:
        for k in range(n):
            ENZYME.append(result[k][2])
            MECANISMES.append(result[k][0])
        return [ENZYME,MECANISMES] # on récupère les mécanismes et enzyme de current qui correspond ici à une combinaison qui fonctionne, il faudrait aussi retirer les doublons dans ENZYME pour l'affichage

'''

def ressystem(LREAC,LPROD,LETIQ,nmax,ENZ,reaction,MOL,REACTIONS,REACPARMOL):
    n=len(LPROD)
    ENZYME=[]
    MECANISMES=[]
    for k in range(n):
        if LETIQ[k]=='ab':
            trouve=False
            for chemin in pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False,reaction,MOL,REACTIONS,REACPARMOL):
                if chemin[1]=='ab':
                    ENZYME.append(chemin[2])
                    MECANISMES.append(chemin[0])
                    trouve=True
                    break
            if not trouve:
                return ('erreur pour le chemin'+str(k))

        if LETIQ[k]=='aob':
            trouvea=False
            trouveb=False
            LENZYME=[]
            MECANISMEa=[]
            MECANISMEb=[]
            for chemin in pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False,reaction,MOL,REACTIONS,REACPARMOL):
                if not trouvea:
                    if not trouveb:
                        if chemin[1]=='a':
                            for e in chemin[2]:
                                LENZYME.append(e)
                            MECANISMEa.append(chemin[0])
                            trouvea=True
                        if chemin[1]=='b':
                            for e in chemin[2]:
                                LENZYME.append(e)
                            MECANISMEb.append(chemin[0])
                            trouveb=True
                    if chemin[1]=='a':
                        for e in chemin[2]:
                            LENZYME.append(e)
                        trouvea=True
                        MECANISMEa.append(chemin[0])
                if not trouveb:
                    if chemin[1]=='b':
                        for e in chemin[2]:
                            LENZYME.append(e)
                        trouveb=True
                        MECANISMEb.append(chemin[0])
                if trouvea and trouveb :
                    break
            if not trouvea and not trouveb:
                return ('erreur pour le chemin'+str(k))
            else:
                ENZYME.append(LENZYME)
                MECANISMES.append((MECANISMEa,MECANISMEb))

        if LETIQ[k]=='a':
            chemin= pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False,reaction,MOL,REACTIONS,REACPARMOL)
            ENZYME.append(chemin[2])
            MECANISME.append(chemin[0])
        if LETIQ[k]=='anb':
            LNEGATION=recherchetemoin(LREAC[k])
            i=0
            nmin=nmax
            NNEG=LNEGATION.length()
            for j in range(NNEG):
                ENZYMENEG=pluscourtcheminneg(enz,numleneg,numero(LPROD[k]),'ab',nmax,LNEGATION[j])
                if ENZYMENEG[1]<nmin:
                    i=j
                    nmin=ENZYMENEG[1]
            ENZYME.append(pluscourtcheminneg(enz,numero2(LREAC[k]),numero(LPROD[k]),'ab',nmin,LNEGATION[j]))
    return [ENZYME,MECANISMES]
    '''
