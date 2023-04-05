import copy

Molecule_Lettre = ['ABTS', 'ADH', 'NADH', 'NAD', 'resazurin', 'HRP', 'AO', 'HRP2', 'POD', 'NR', 'G_1DH', 'O2', 'DAF', 'glucoseext', 'glucose', 'acetoneext', 'acetone', 'Lactateext', 'Lactate', 'EtOHext', 'EtOH', 'NO3ext', 'NO3', 'NO2ext', 'NO2', 'HRP', 'CCia5', 'resazurin', 'CCib5', 'resorufin', 'HRP2', 'NADH', 'CCf4', 'NADN', 'AO', 'isopropanol', 'CCf3', 'CCio3', 'ADH', 'CCia2', 'CCib2', 'CCfa2', 'NAD', 'G_1DH', 'CCia1', 'CCib1', 'CCfa1', 'gluconolacrone', 'NO', 'volatNO', 'O2', 'Cf6', 'NO2b', 'Cf5', 'N2O3', 'DAF', 'Cf4', 'DAFF', 'NR', 'Cia3', 'Cib3', 'Cfa3', 'Cfb3', 'Cio3', 'Cia2', 'Cfa2', 'Cfb2', 'Cio2', 'Cia1', 'Cib1', 'Cfa1', 'ABTSOX', 'DDf3', 'ABTS', 'LO', 'DDf2', 'DDio2', 'H2O2', 'Pyruvate', 'DDia1', 'DDib1', 'DDfa1', 'DDfb1', 'DDio1', 'acetaldehyde', 'POD', 'DDia5', 'DDib5']
MOLECULESS = [[False, 1] for i in range(len(Molecule_Lettre))]

def nombre_lettre(k): #Ne sert a rien à part me donner le nom d'une molécule associé à un nombre
    print(Molecule_Lettre[k])
    return
for k in range(len(Molecule_Lettre)): #ne sert rien à part me donner le numéro correspondant à molécule particulière.
    if 'resorufin' == Molecule_Lettre[k]:
        #print(k)
        break



#MOLECULESS[0] = [True, 630676000] #ABTS
MOLECULESS[1] = [True, 92078816] #ADH
#MOLECULESS[2] = [True, 31533] #NADH
MOLECULESS[3] = [True, 1e9] #NAD
MOLECULESS[4] = [True, 315338394] #resazurin
MOLECULESS[5] = [True, 6565] #HRP
MOLECULESS[6] = [True, 175012] #AO
#MOLECULESS[7] = [True, 6565] #HRP2
#MOLECULESS[8] = [True, 21884] #POD
#MOLECULESS[9] = [True, 26488392] #NR
MOLECULESS[10] = [True, 36295404] #G_1DH
#MOLECULESS[11] = [True, 79780608] #O2
#MOLECULESS[12] = [True, 63067600] #DAF
MOLECULESS[13] = [True, 38319820] #gluxoseext
MOLECULESS[15] = [True, 315338397] #actenoneext
#MOLECULESS[17] = [True, 3153380] #lactateext
#MOLECULESS[19] = [True, 5486881] #EtOHext
#MOLECULESS[21] = [True, 315338000] #NO3ext
#MOLECULESS[13] = [True, 261356] #NO2ext
MOLECULESS[77] = [True, 20000] #H202

REACTIONS = [([13], [14], 5e-3 ), ([15], [16], 1e-2), ([17], [18], 5e-3), ([19], [20], 1e-2), ([21], [22], 5e-3), ([23], [24], 5e-3), ([5, 77], [26], 1.15157e-005), ([26], [5, 77], 24), ([26, 4], [28], 7.77313e-006), ([28], [26, 4], 1), ([28], [5, 29], 240), ([7, 2], [32], 2.91492e-009 ), ([32], [7, 2], 0.0009), ([32], [7, 33], 0.009), ([6, 35], [36], 5.82985e-008), ([36], [6, 35], 15), ([36], [37, 77], 150), ([37], [6, 7], 10000), ([1, 2], [39], 9.50049e-012), ([39], [1, 2], 0.033), ([39, 16], [40], 1.08824e-006), ([40], [39, 16], 0.07), ([40], [41, 3], 0.7), ([41], [1, 35], 0.33), ([10, 3], [44], 1.8077e-006), ([44], [10, 3], 40), ([44, 14], [45], 9.71641e-009), ([45], [44, 14], 20), ([45], [46, 2], 200), ([46], [10, 47], 400), ([48], [49], 0.0019), ([11, 48], [51], 7.77313e-006), ([51], [11, 48], 0.2), ([51], [11, 52], 2), ([48, 52], [53], 0.00777313), ([53], [48, 52], 200), ([53], [54], 2000), ([12, 54], [56], 0.00777313), ([56], [12, 54], 200), ([56], [57], 2000), ([9, 24], [59], 1.05042e-006), ([59], [9, 24], 0.2), ([9, 2], [60], 1.94328e-006), ([60], [9, 2], 0.2), ([59, 2], [61], 1.94328e-006), ([60, 24], [62],1.05042e-006), ([61], [59, 2], 0.2), ([62], [60, 24],0.2), ([61], [63, 48], 2), ([62], [63, 48], 2), ([63], [9, 3], 10000), ([9, 22], [64], 5.44119e-005), ([64], [9, 22], 21), ([64, 2], [65], 0.000204045), ([60, 22], [66], 5.44119e-005), ([65], [64, 2], 21), ([66], [60, 22], 21), ([65], [67, 24], 210), ([66], [67, 24], 210), ([67], [9, 3], 10000), ([10, 3], [68], 3.37961e-007), ([68], [10, 3], 40), ([68, 14], [69], 9.71641e-009), ([69], [68, 14], 20), ([69], [70, 2], 200), ([70], [10, 47], 400), ([2, 71], [72], 3.88656e-006), ([72], [2, 71], 10), ([72], [3, 0] ,100), ([74, 18], [75], 4.63083e-006), ([75], [74, 18], 23.83), ([75], [76, 77], 238.3), ([76], [74, 78], 10000), ([1, 20], [79], 4.8661e-007), ([79], [1, 20], 30.8), ([1, 3], [80], 0.000119706), ([80], [1, 3], 30.8), ([79, 3], [81], 0.000119706), ([80, 20], [82], 4.8661e-007), ([81], [79, 3], 30.8), ([82], [80, 20], 30.8), ([81], [83, 84], 308), ([82], [83, 84], 308), ([83], [1, 2], 10000), ([8, 0], [86], 0.000590758), ([86], [8, 0], 76), ([86, 77], [87], 1.64099e-005), ([87], [86, 77], 76), ([87], [8, 71], 760)]


def Num_Reaction_Actif(PRESENCE): #Fonction auxilliaire qui me donne le numéro des réactions à considérer vu les molécules présentes
    Inter_Mol = []
    for a in range(len(PRESENCE)):
        if PRESENCE[a][0]:
            Inter_Mol.append(a)
    
    Inter_Reac = []
    for i in range(len(REACTIONS)):
        n = len(REACTIONS[i][0])
        compteur = 0
        while compteur < n:
            a = REACTIONS[i][0][compteur]
            if PRESENCE[a][0] == False:
                break
            compteur = compteur + 1
            if compteur ==  n:
                Inter_Reac.append(i)
    return Inter_Reac

def Calc_Vitesse(PRESENCE, li):#Fonction auxilliaire qui me donne la vitesse des réactions à considérer.
    vitesse = []
    for k in li:
        v = REACTIONS[k][2]
        n = len(REACTIONS[k][0])
        compteur = 0
        while compteur < n:
            a = REACTIONS[k][0][compteur]
            v = v * PRESENCE[a][1]
            compteur = compteur + 1
        vitesse.append(v)
    return vitesse

def Cre_PresenceBis(PRESENCE,Avancement,Inter):#Fonction auxiliaire qui me crée la liste actualiser des REACTIFS seulements après avoir effectuer la réaction
    reponse = PRESENCE
    for k in range(len(Inter)):
        i = Inter[k]
        n = len(REACTIONS[i][0])
        compteur = 0
        while compteur < n:
            a = REACTIONS[i][0][compteur]
            reponse[a][1] = reponse [a][1] - Avancement[k]
            compteur = compteur + 1
    for i in range(len(reponse)):
        if reponse[i][1] < 0:
            reponse[i][0] = False
            reponse[i][1] = 0
    return reponse

def Cre_Produit(INTER, avancement):
    reponse = []
    for k in range(len(INTER)):
        i = INTER[k]
        n = len(REACTIONS[i][1])
        compteur = 0
        while compteur < n : 
            a = REACTIONS[i][1][compteur]
            reponse.append([a, avancement[k]])
            compteur = compteur + 1
    return reponse
    

def simulation(nombre_etape, precision):
    compteur = 0 #Compteur pour le nombre de pas
    minimum = min(MOLECULESS[k][1] for k in range(len(MOLECULESS)) if MOLECULESS[k][0])
    pas = minimum/precision #la quantité de matière que je veux créer a chaque étape dans le produit, (peut faire disparaitre +, problème de création de matière)
    PRESENCE = copy.deepcopy(MOLECULESS)
    while compteur < nombre_etape:
        Inter = Num_Reaction_Actif(PRESENCE)#Liste des réactions à considérer
        Vitesse_Inter = Calc_Vitesse(PRESENCE, Inter)#calcul des différentes vitesse à considérer
        if sum(Vitesse_Inter) < 0.0000001: #Je prépare ma règle de 3, si le systeme est à l'équilibre, alors arrêt de la fonction.
            print('nombre etape :' , compteur)
            return PRESENCE
        RegleDeTrois = pas/sum(Vitesse_Inter)
        Avancement_reel = [Vitesse_Inter[i] * RegleDeTrois for i in range(len(Vitesse_Inter))] #Avancement pour chaque réaction
        PRESENCE = Cre_PresenceBis(PRESENCE, Avancement_reel, Inter) #Je soustrais aux réactifs ce qui a été consommé
        Produit = Cre_Produit(Inter, Avancement_reel)#Je crée les produits
        for k in range(len(Produit)): 
            i = Produit[k][0]
            PRESENCE[i][0] = True
            PRESENCE[i][1] = PRESENCE[i][1] + Produit[k][1] #j'ajoute les produits à la liste PRESENCE
        compteur = compteur + 1
    return (PRESENCE[29])#Retourne la composition du milieu Ici le 29 c'est le resorufin, cohérent, on voit que si on ajoute ou non NAD on obtien soit 1 soit un chiffre enorme. 
    

simulation(500099, 1000)
