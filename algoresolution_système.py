from algoPSC2_0_1 import pluscourtchemin

def check(Current,LETIQ,nmax,LREAC,LPROD):
    ENZYME=[]
    n=len(LPROD)
    for ligne in Current:
        for enz in ligne[2]:
            ENZYME.append(enz)
    for k in range(n):
        chemins=pluscourtchemin(ENZYME,LREAC[k],LPROD[k],nmax)
        if LETIQ[k]=='ab':
            for chemin in chemins:
                if chemin[1]!='ab':
                    return False
        if LETIQ[k]=='aob':
            inda=False
            indb=False
            for chemin in chemins:
                if chemin[1]=='a':
                    inda=True
                if chemin[1]=='b':
                    indb=True
            if not inda or not indb:
                return False
    return True


def research(Allreactions,LETIQ,nmax,LREAC,LPROD):
    n=len(Allreactions)
    i=[0]*(n+1)
    Borne = [0]*n
    Current=[[]]*n
    for k in range(n):
        Borne[k]=len(Allreactions[k])
    p = 0
    while i[n]==0:
        for k in range(n):
            Current[k]=Allreactions[i[k]]
        if (check(Current,LETIQ,nmax,LREAC,LPROD)):
            return Current
        i[0]+=1
        while i[p]==MAX[p]:
            i[p]=0
            p+=1
            i[p]+=1
            if(i[p]!=MAX[p]):
                p=0
    return 0

def res(LREAC,LPROD,LETIQ,nmax, ENZ):
    n=len(LPROD)
    ENZYME=[]
    MECANISMES=[]
    Allreactions=[]
    for k in range(n):
        Allreactions.append(pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax))

    result=research(Allreactions,LETIQ,nmax, LREAC,LPROD)
    if result==0:
        return ("pas de solution au système")
    else:
        for k in range(n):
            ENZYME.append(result[k][2])
            MECANISMES.append(result[k][0])
        return [ENZYME,MECANISMES]



def ressystem(LREAC,LPROD,LETIQ,nmax, ENZ):
    n=len(LPROD)
    ENZYME=[]
    MECANISMES=[]
    for k in range(n):
        if LETIQ[k]=='ab':
            trouve=False
            for chemin in pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False):
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
            for chemin in pluscourtchemin(ENZ,LREAC[k],LPROD[k],nmax,False):
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
            chemin= pluscourtchemin(ENZ,LREAC[k],LPROD[k],'a',nmax)
            ENZYME.append(chemin[2])
            MECANISME.append(chemin[0])
    """
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
    """
    return [ENZYME,MECANISMES]
