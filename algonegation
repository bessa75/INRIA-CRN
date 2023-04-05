from recherche_chemin import pluscourtchemin
def algo_negation(A,B,C,ENZ,CYCLES,CYCLESPARMOL,REACPARMOL,reaction,MOL,REACTIONS,n,reac):
    PRODINT=[]
    REACT=[]
    ACONSOMMER=[C]
    CYCLESa=[False]
    BOOL=[]
    L=[]
    CYCLEC=CYCLESPARMOL[C]
    RES=[]
    for k in range (0,len(CYCLEC)):
        #print(k)
        #print(CYCLEC[k])
        bool=False
        m1,m2=REACTIONS[CYCLEC[k]]
        if len(m1)==1 and C in m2 and m1[0] not in ACONSOMMER:
            ACONSOMMER.append(m1[0])
            bool=True
            CYCLESa.append(CYCLEC[k])
        if len(m2)==1 and C in m1 and m2[0] not in ACONSOMMER:
            ACONSOMMER.append(m2[0])
            if bool==False:
                CYCLESa.append(CYCLEC[k])
    print("Molécules dont la consommation permettrait de consommer " + str(MOL[C])+" : "+str([MOL[i] for i in ACONSOMMER]))
    for i in range (0,len(ACONSOMMER)):
        M=ACONSOMMER[i]
        REACC=REACPARMOL[M]
        CYCLEC=CYCLESPARMOL[M]
        for k in range (0,len(REACC)):
            r,p=REACTIONS[REACC[k]]
            if len(r)==1:
                if (p,r) not in CYCLEC and (r,p) not in CYCLEC and M==C:
                    return("C est une molecule instable")
            elif len(r)==2: #on considere la molecularite inferieure ou egale a 2
                im=min([i for i in range (0,len(r)) if r[i]!=M])
                prodint=r[im]
                if r[im]==B:
                    print("Fin de l'algorithme, voici le résultat :'")
                    return(([[REACC[k],CYCLESa[i]],'a',[]]))
                PRODINT.append(prodint)
                REACT.append(REACC[k])
                if M==C:
                    BOOL.append(False)
                else:
                    BOOL.append((True,i))
    print('********************************************')
    PRESENCE=pluscourtchemin(ENZ,[B],0,n,False,reac,MOL,REACTIONS,REACPARMOL,bool=True)
    ENZint=[]
    for k in range (0,len(PRODINT)):
        m=PRODINT[k]
        #print(MOL[m])
        PRESm=PRESENCE[m]
        #print(PRESm)
        #print(PRESm[0])
        if PRESm[0] and 'a' in PRESm[2]:
            print("produit intermédiaire considéré : "+MOL[m])
            P2=pluscourtchemin(ENZ,[B,0],m,n,False,reac,MOL,REACTIONS,REACPARMOL,bool=False)
            for mec in P2:
                if mec[1]=='a':
                    meca=mec[0]
                    meca.append([REACT[k]])
                    if BOOL[k]!=False:
                        meca.append([CYCLESa[BOOL[k][1]]])
                    Res=(meca,'a',mec[2])
                    RES.append(Res)
                if mec[1]=='e':
                    meca=mec[0]
                    meca.append([REACT[k]])
                    if BOOL[k]!=False:
                        meca.append([CYCLESa[BOOL[k][1]]])
                    Res=(meca,'e',mec[2])
                    RES.append(Res)
                    ENZint.append([MOL[i] for i in Res[2]])
    """ PRESENCE : liste des molécules présente.
    Evolue au cours de l'algorithme.
    Chaque case coorespond à une molécule et est de la forme :
    [boolean marquant la présence,
    liste des réaction l'ayant produit selon le triplet (numéro de réaction, numéro d'étape, étiquette),
    liste des étiquettes avec lesquelles la molécule a été produite,
    liste des doublets (numéro d'étape, étiquette) avec lesquelles la molécule a été produite]"""
    print("Jeux d'enzymes interdites : "+str(ENZint))
    print("Fin de l'algorithme: voici le résultat :")
    print("")
    return(RES)


