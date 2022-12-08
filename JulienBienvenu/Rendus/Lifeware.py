from abc import ABC, abstractmethod
import colorama
from colorama import Fore
import numpy as np
import os
import itertools
import pandas as pd
from pprint import pprint

adress = #yourworkenv
mail = #yourmail
pw = #yourpassword
os.chdir(adress)

### SOAP API
from zeep import Client
import hashlib

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256(pw.encode("utf-8")).hexdigest()
client = Client(wsdl)


def hashed_metabolite(name):
    parameters = (mail,password,name)
    resultString = client.service.getLigandStructureIdByCompoundName(*parameters)
    return resultString

### Operations on sets/listes


def findsets(s, n='full'):
    if n=='full':
        n=len(s)
    return list(itertools.permutations(s, n))

def obj_is_in_list(Obj,List):
    Listofcompare = [(Obj==Listj) for Listj in List]
    Bool = False
    for b in Listofcompare:
        Bool = Bool or b
    return Bool


def dictionnary_lengths(dic):
    Dic = dict()
    for key in dic:
        Dic[key] = len(dic[key])
    return Dic

def disjoint(a,b):
    return (len(set(a).intersection(b))==0)

def union(a,b):
    return list(set(a).union(set(b)))

def included(a,b):
    return set(a).issubset(set(b))

def diff(a,b):
    return list(set(a).difference(set(b)))


### BRENDA Reading tools

def comptability_lecture(Lecture):
    Metabolites = set()
    Enzymes = set()
    ECs = set()
    i = 0
    for l in Lecture :
        i+=1
        Metabolites = Metabolites.union(set(l[3]+l[4]))
        ECs = ECs.union({l[0]})
        Enzymes = Enzymes.union({l[1]})
        if i%50==0:
            s = "\r "+str(i)+" ECs : "+str(len(ECs))+", Enzymes : "+str(len(Enzymes))+", Metabolites : "+str(len(Metabolites))
            print(s, end="")
    s = "\r "+str(i)+" ECs : "+str(len(ECs))+", Enzymes : "+str(len(Enzymes))+", Metabolites : "+str(len(Metabolites))
    print(s, end="")

    crypte=open("Metabolites.tsv",'w')
    j=0
    for m in Metabolites :
        crypte.write(m+'\n')
        j+=1
        s = "\r Metabolites written : "+str(j)
        print(s, end="")
    crypte.close()

def redundancy_lecture(Lecture):
    EC_Enzymes = set()
    i = 0
    for l in Lecture :
        i+=1
        EC_Enzymes = EC_Enzymes.union({(l[0],l[1])})
        if i%50==0:
            s = "\r "+str(i)+" EC-Enzymes couples : "+str(len(EC_Enzymes))
            print(s, end="")
    s = "\r "+str(i)+" EC-Enzymes : "+str(len(EC_Enzymes))
    print(s, end="")
    L=[]
    for e in EC_Enzymes:
        L.append(e[1])
    values,counts = np.unique(L,return_counts=True)
    print('\n')
    i=0
    for c in counts :
        if c>1:
            v = values[i]
            for e in EC_Enzymes :
                if e[1]==v:
                    print(e)
        i+=1

def clean_lecture(Lecture,badfilename,aff=True):
    L = []
    crypte=open(badfilename,'r')
    Texte=crypte.readlines()
    BadLists = [int(i) for i in Texte]
    crypte.close()
    i,j = 0,0
    for l in Lecture:
        if not(i in BadLists):
            L.append(l)
        else :
            j+=1
        if aff :
            s = "\r Entries read "+str(i+1)+" Entries rejected : "+str(j)
            print(s, end="")
        i+=1
    return L

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def with_KM(Lecture):
    i,j=0,0
    for l in Lecture:
        i+=1
        if is_number(l[6]):
            j+=1
        s = "\r Entries read "+str(i)+" Entries with Km : "+str(j)
        print(s, end="")


def from_brenda_to_clean_lecture(source,badentries):
    crypte=open(source,'r')
    Texte=crypte.readlines()
    Lecture = [i.split("\t") for i in Texte]

    for i in range(len(Lecture)) :
        l = Lecture[i]
        l[3],l[4] = l[3].split("-+-"),l[4].split("-+-")
        Lecture[i] = l[:12]
    crypte.close()
    Lecture = clean_lecture(Lecture,badentries)
    return Lecture

## Brenda Cleaning and Lexicon building
def lexicon_metabolites(Lecture):
    Metabolites = set()
    Dictionnary = dict()
    Bad_Entries = []
    i,j,k = 0,0,0
    for l in Lecture :
        i+=1
        L = [m for m in l[3]+l[4] if not(m in Metabolites)]
        for m in L :
            m_id = hashed_metabolite(m)
            if not(m_id==''):
                j+=1
                if m_id in Dictionnary :
                    Dictionnary[m_id].append(m)
                else :
                    Dictionnary[m_id] = [m]
            else :
                k+=1
                Bad_Entries.append(m)
            Metabolites.add(m)
        s = "\r Entries read "+str(i)+", Metabolites not recognized : "+str(k)+", Recognized : "+str(j)
        print(s, end="")
    return Dictionnary, Bad_Entries

def Write_Metabolites(Lecture):
    Metabolites = set()
    crypte=open("Metabolites.tsv",'w')
    i,j=0,0
    LL = Lecture.copy()
    LL.pop(0)
    for l in LL :
        i+=1
        L = [m for m in l[3]+l[4] if not(m in Metabolites)]
        for m in L :
            j+=1
            crypte.write(m+'\n')
            Metabolites.add(m)
        s = "\r Entries read "+str(i)+", Metabolites written : "+str(j)
        print(s, end="")
    crypte.close()

def Write_Metabolites_Id(filename):
    crypte=open("Metabolites.tsv",'r')
    Texte=crypte.read().split('\n')
    crypte.close()
    i=11396
    for m in Texte[i:] :
        i+=1
        crypte=open(filename,'a')
        if m.isprintable():
            m_id = hashed_metabolite(m)
        else :
            m_id ='000 not unicode'
        if m_id=='':
            m_id='000'
        if m_id==None:
            m_id='000'
        crypte.write(m+'\t'+m_id+'\n')
        crypte.close()
        s = "\r Entries read and identified "+str(i)
        print(s, end="")

def Write_Metabolites_Id_with_prefilling(filename,prefilename):
    crypte=open(prefilename,'r')
    Texte=crypte.read().split('\n')
    Lecture = [l.split('\t') for l in Texte]
    crypte.close()
    i,j=4379,0
    for u in Lecture[i:] :
        i+=1
        m,m_i = u[0],u[1]
        crypte=open(filename,'a')
        if m_i=='?':
            j+=1
            if m.isprintable():
                m_id = hashed_metabolite(m)
            else :
                m_id ='000 not unicode'
            if m_id=='':
                m_id='000'
            if m_id==None:
                m_id='000'

        else :
            m_id = m_i
        crypte.write(m+'\t'+m_id+'\n')
        crypte.close()
        s = "\r Entries read and identified "+str(i)+" previously unindentified : "+str(j)+" with id : "+m_i+'|'+m_id
        print(s, end="")





### Lexique exploitation

def Lexiques_building(filename):
    crypte=open(filename,'r')
    Texte=crypte.read().split('\n')
    crypte.close()
    meta_to_key = dict()
    key_to_metas = dict()
    Table = [l.split('\t') for l in Texte]
    Table.pop(-1)
    i=0
    for u in Table:
        i+=1
        if len(u)<2:
            print(i,u)
        meta,key = u[0],u[1]
        meta_to_key[meta]=key
        if key in key_to_metas :
            key_to_metas[key].append(meta)
        else :
            key_to_metas[key] = [meta]
    return meta_to_key,key_to_metas


def Write_Metabolites_Id_with_archive(filename,archivename):
    crypte=open("Metabolites.tsv",'r')
    Texte=crypte.read().split('\n')
    crypte.close()
    Lex1,Lex2 = Lexiques_building(archivename)
    i,j = 0,0
    for m in Texte :
        i+=1
        crypte=open(filename,'a')
        if m in Lex1:
            m_id = Lex1[m]
        else :
            m_id = '?'
            j+=1
        crypte.write(m+'\t'+m_id+'\n')
        crypte.close()
        s = "\r Entries read : "+str(i)+", unidentified :"+str(j)
        print(s, end="")

#Write_Metabolites_Id_with_archive("Lexique_test.tsv","LexiqueArchive.csv")

#
# for meta in Unknowns :
#     print(meta,hashed_metabolite(meta))


def Numbering(Reading,Lexique1):
    NumberedLecture = []
    i=0
    for l in Reading :
        n = l.copy()
        n[3]=[Lexique1[s] for s in l[3]]
        n[4]=[Lexique1[s] for s in l[4]]
        NumberedLecture.append(n)
    return NumberedLecture


### Class building


class Metabolite:
    def __init__(self,name):
        self.name = name
        self.synonyms = [name]

    def add_synonyms(self,liste):
        self.synonyms = self.synonyms + list(set(liste) - set(self.synonyms))

    def sameas(self,metabolite):
        return obj_is_in_list(self.name,metabolite.synonyms) or obj_is_in_list(metabolite.name,self.synonyms)

# met1 = Metabolite("Glucose")
# met1.add_synonyms(["Sugar"])
#
# met2 = Metabolite("Sugar")
# met1.sameas(met2)

# class Reaction:
#     def

def in_and_out_formating(in_and_outs):
    values,counts = np.unique(in_and_outs,return_counts=True)
    dic = dict(zip(values,counts))
    for label in ['in','out','inh'] :
        if label not in dic :
            dic[label] = 0
    #dic.pop('-')
    return dic

class Abstract_Reaction:
    def __init__(self,in_and_outs,env,kin):
        self.config = in_and_out_formating(in_and_outs)
        self.environment = env
        self.kinetics = kin


class Concrete_Reaction:
    def __init__(self,l):
        self.liste = l
        self.ec = self.liste[0]
        self.enzyme_name = self.liste[1]
        self.organism = self.liste[2]
        self.reactants = self.liste[3]
        self.products = self.liste[4]
        Lr,Lp = self.liste[3],self.liste[4]
        S = 'm'+Lr[0]
        for s in range(len(Lr)-1):
            S+=' + m'+Lr[s+1]
        S += ' => m'+Lp[0]
        for s in range(len(Lp)-1):
            S+=' + m'+Lp[s+1]
        S+='. %'+self.ec+', '+self.enzyme_name+', in '+self.organism
        self.trace = S


    def soft_match(self,abs_react):
        nb_in, nb_out = len(self.reactants), len(self.products)
        Nb_in, Nb_out = abs_react.config['in'], abs_react.config['out']
        return ((nb_in == Nb_in) and (nb_out == Nb_out))

    def hard_match(self,abs_react):
        nb_in, nb_out = len(self.reactants), len(self.products)
        Nb_in, Nb_out = abs_react.config['in'], abs_react.config['out']
        return ((nb_in >= Nb_in) and (nb_out >= Nb_out))

class Concrete_Reaction_Library:
    def __init__ (self,big_liste):
        self.library = []
        for l in big_liste :
            self.library.append(Concrete_Reaction(l))

    def select_candidates(self,CRN_abstracts, soft_matching=True):
        Book = CRN_abstracts.reaction_book()
        Precandidates = dict()
        n = len(self.library)
        for key in Book :
            l = []
            for i in range(n) :
                react = self.library[i]
                if soft_matching and (react.soft_match(Book[key])):
                    l.append(i)
                if (not soft_matching) and (react.hard_match(Book[key])):
                    l.append(i)
            Precandidates[key] = l
        return Precandidates


def dictionnarize(liste):
    dic = dict()
    for i in range(len(liste)):
        dic[liste[i]] = i
    return dic

def coloring(i,CRN_blue):
    for j in range(len(CRN_blue[i])):
        color = CRN_blue[i,j]
        if (color =='vert' or color=='rouge') :
            CRN_blue[i,j] = 'vert'

    return CRN_blue

def score(j,CRN_blue):
    values,counts = np.unique(CRN_blue[:,j],return_counts=True)
    dic = dict(zip(values,counts))
    if 'vert' not in dic :
        dic['vert'] = 0
    if 'rouge' not in dic :
        dic['rouge'] = 0
    return dic

def nb_red_out(j,CRNe,CRN_blue):
    N, nb_m = 0,len(CRNe[:,0])
    for i in range(nb_m):
        color,type = CRN_blue[i,j],CRNe[i,j]
        if color == 'rouge' and type =='out' :
            N+=1
    return N


class CRNa :

    def __init__(self,liste):
        self.abstract = np.array(liste)
        self.extract = self.abstract[1::,1::]
        self.enzymes = self.abstract[0,1:]
        self.metabolites = self.abstract[1:,0]

    def reaction_book(self):
        EC_dic = dictionnarize(self.enzymes)
        CRNe = self.extract
        Book = dict()
        for key in dictionnarize(self.enzymes):
            Book[key] = Abstract_Reaction(CRNe[:,EC_dic[key]],"any","any")
        return Book

    def blueprint(self):
        return np.array([[(j=="-")*'gris' + (j=='out' or j=='in' or j=='inh')*'rouge' for j in listj] for listj in list(np.copy(self.extract))])

    def Delta(self,A,E):
        CRNe = self.extract
        EC_dic, Meta_dic = dictionnarize(self.enzymes), dictionnarize(self.metabolites)
        return CRNe[Meta_dic[A],EC_dic[E]]


    def parcours(self,init):
        CRNe = self.extract
        CRN_EC = self.enzymes
        CRN_Meta = self.metabolites
        CRN_blue = self.blueprint()
        nb_e, nb_m = len(CRN_EC),len(CRN_Meta)
        EC_dic, Meta_dic = dictionnarize(self.enzymes), dictionnarize(self.metabolites)

        for Meta_i in init:
            CRN_blue = coloring(Meta_dic[Meta_i],CRN_blue)

        Steps = [['Empty',init,'Empty']]
        To_do_list = [i for i in range(nb_e)]
        while len(To_do_list)>0 :
            J, G, R = -1,-1,-1
            for j in range(len(To_do_list)):
                dic_j = score(To_do_list[j],CRN_blue)
                g,r = dic_j['vert'],dic_j['rouge']
                if g>0 :
                    if (g==G) and (r==R):
                        n,N = nb_red_out(To_do_list[j],CRNe,CRN_blue),nb_red_out(To_do_list[J],CRNe,CRN_blue)
                        if n>N :
                            J,G,R = j,g,r
                    if (g==G) and (r<R) :
                        J,G,R = j,g,r
                    if g > G :
                        J,G,R = j,g,r
            J = To_do_list.pop(J)
            Fillable, Filled = [],[]
            for i in range(nb_m):
                if CRN_blue[i,J] == 'vert' :
                    Filled.append(CRN_Meta[i])
                if CRN_blue[i,J] == 'rouge' :
                    Fillable.append(CRN_Meta[i])
                    CRN_blue = coloring(i,CRN_blue)
            Steps.append([CRN_EC[J],Fillable,Filled])
        return Steps


def insertion(CRNtab,E,Motif):
    CRN = np.array(CRNtab.abstract)
    Cj = np.where(CRN[0,:]==E)[0][0]
    C = CRN[:,Cj]
    pCRN = np.delete(CRN, Cj, axis=1)
    n,m = pCRN.shape
    u,v = Motif.shape
    M = np.zeros((n+u-1,m+v-1),dtype="<U100")
    for i in range(n+u-1):
        for j in range(m+v-1):
            M[i,j]='-'
    for i in range(n):
        for j in range(m):
            M[i,j] = pCRN[i,j]
    for i in range(u-1):
        for j in range(v-1):
            M[n+i,m+j]=Motif[i+1,j+1]
    for k in range(v-1):
        M[0,m+k] = E+Motif[0,k+1]
    for k in range(u-1):
        M[n+k,0] = E+Motif[k+1,0]
    for i_in in np.where(C=='in'):
        M[i_in,m]='in'
    for i_out in np.where(C=='out'):
        M[i_out,-1]='out'
    return CRNa(M)

def simple_insert(CRNtab,E):
    Simple_Mo = np.array([
    ['-','e','f'],
    ['a','out','in']])
    return insertion(CRNtab,E,Simple_Mo)

def fork_insert(CRNtab,E):
    CRN = np.array(CRNtab.abstract)
    Cj = np.where(CRN[0,:]==E)[0][0]
    C = CRN[:,Cj]
    ins = np.where(C=='in')[0]
    pCRN = np.copy(CRN)
    n,m = pCRN.shape
    u,v = len(ins),len(ins)
    M = np.zeros((n+u,m+v),dtype="<U100")
    for i in range(n+u):
        for j in range(m+v):
            M[i,j]='-'
    for i in range(n):
        for j in range(m):
            M[i,j] = pCRN[i,j]
    for i in range(len(ins)):
        M[0,m+i] = E+pCRN[ins[i],0]
        M[n+i,0] = pCRN[ins[i],0]+'2'+E
        M[ins[i],Cj] = '-'
        M[ins[i],m+i] = 'in'
        M[n+i,Cj] = 'in'
        M[n+i,m+i] = 'out'
    return CRNa(M)


## SISO* searching (le plat de résistance)

def concatenate(dic):
    s = ''
    for key in dic.keys():
        s+=key+': '+dic[key]+'; '
    return s

def f(CRNabs,Library,Starter,Grey=False,One=False):
    Parcours = CRNabs.parcours(list(Starter.keys()))
    Precandidats = Library.select_candidates(CRNabs,not(Grey))
    Blacklist_EC = []
    D_EC = dict()
    Blacklist_M = list(Starter.values())
    D_M = Starter.copy()
    Greylist_M = []
    global S, S_E
    S = []
    S_E = []
    def aux(Steps,B_EC,B_M):
        if len(Steps)==0 :
            S.append([D_EC,D_M])
        else :
            E,A,X = Steps[0]

            A_in = [key for key in A if (CRNabs.Delta(key,E)=='in')]
            A_out = [key for key in A if (CRNabs.Delta(key,E)=='out')]
            x_in = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='in')]
            x_out = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='out')]
            PCandidats = [Library.library[i] for i in Precandidats[E]]
            Candidats = []
            for c in PCandidats :
                c_in,c_out = c.reactants,c.products
                if not(c.ec in B_EC):
                    if included(x_in,c_in) and included(x_out,c_out):
                        a_in, a_out = diff(c_in,x_in),diff(c_out,x_out)
                        if disjoint(B_M,a_in+a_out):
                            Candidats.append(c)
            if len(Candidats)!=0 :
                for c in Candidats :
                    e,ec = c.enzyme_name,c.ec
                    D_EC[E] = e
                    c_in,c_out = c.reactants,c.products
                    e_in, e_out = diff(c_in,x_in),diff(c_out,x_out)
                    permutations_in, permutations_out = findsets(e_in),findsets(e_out)
                    for pi in permutations_in:
                        for po in permutations_out:
                            for i in range(len(A_in)):
                                D_M[A_in[i]] = pi[i]
                            for j in range(len(A_out)):
                                D_M[A_out[j]] = po[j]
                            aux(Steps[1:],B_EC+[ec],B_M+e_in+e_out)


    def aux_grey(Steps,B_EC,B_M,G_M):
        if len(Steps)==0 :
            s_EC = concatenate(D_EC)
            if not(s_EC in S_E):
                S_E.append(s_EC)
                print(D_EC)
            S.append([D_EC,D_M,G_M])
        else :
            E,A,X = Steps[0]

            A_in = [key for key in A if (CRNabs.Delta(key,E)=='in')]
            n_in = len(A_in)
            A_out = [key for key in A if (CRNabs.Delta(key,E)=='out')]
            n_out = len(A_out)
            x_in = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='in')]
            x_out = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='out')]
            PCandidats = [Library.library[i] for i in Precandidats[E]]
            Candidats = []
            for c in PCandidats :
                c_in,c_out = c.reactants,c.products
                if not(c.ec in B_EC):
                    if included(x_in,c_in) and included(x_out,c_out):
                        a_in, a_out = diff(c_in,x_in),diff(c_out,x_out)
                        if disjoint(B_M,a_in+a_out):
                            Candidats.append(c)
            if len(Candidats)!=0 :
                for c in Candidats :
                    e,ec = c.enzyme_name,c.ec
                    D_EC[E] = e
                    c_in,c_out = c.reactants,c.products
                    e_in, e_out = diff(diff(c_in,x_in),G_M),diff(diff(c_out,x_out),G_M)
                    if (len(e_in)>= n_in)and(len(e_out)>= n_out):
                        permutations_in, permutations_out = findsets(e_in,n_in),findsets(e_out,n_out)
                        for pi in permutations_in:
                            for po in permutations_out:
                                for i in range(len(A_in)):
                                    D_M[A_in[i]] = pi[i]
                                for j in range(len(A_out)):
                                    D_M[A_out[j]] = po[j]
                                aux_grey(Steps[1:],B_EC+[ec],B_M+e_in+e_out,G_M+diff(e_in,pi)+diff(e_out,po))
    def aux_grey_uno(Steps,B_EC,B_M,G_M,Trace):
        if len(Steps)==0 :
            return [D_EC,D_M,G_M,Trace]

        else :
            E,A,X = Steps[0]
            A_in = [key for key in A if (CRNabs.Delta(key,E)=='in')]
            n_in = len(A_in)
            A_out = [key for key in A if (CRNabs.Delta(key,E)=='out')]
            n_out = len(A_out)
            x_in = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='in')]
            x_out = [D_M[key] for key in X if (CRNabs.Delta(key,E)=='out')]
            PCandidats = [Library.library[i] for i in Precandidats[E]]
            Candidats = []
            for c in PCandidats :
                c_in,c_out = c.reactants,c.products
                if not(c.ec in B_EC):
                    if included(x_in,c_in) and included(x_out,c_out):
                        a_in, a_out = diff(c_in,x_in),diff(c_out,x_out)
                        if disjoint(B_M,a_in+a_out):
                            Candidats.append(c)
            if len(Candidats)!=0 :
                for c in Candidats :
                    e,ec = c.enzyme_name,c.ec
                    D_EC[E] = e
                    c_in,c_out = c.reactants,c.products
                    e_in, e_out = diff(diff(c_in,x_in),G_M),diff(diff(c_out,x_out),G_M)
                    if (len(e_in)>= n_in)and(len(e_out)>= n_out):
                        permutations_in, permutations_out = findsets(e_in,n_in),findsets(e_out,n_out)
                        for pi in permutations_in:
                            for po in permutations_out:
                                for i in range(len(A_in)):
                                    D_M[A_in[i]] = pi[i]
                                for j in range(len(A_out)):
                                    D_M[A_out[j]] = po[j]
                                S = aux_grey_uno(Steps[1:],B_EC+[ec],B_M+e_in+e_out,G_M+diff(e_in,pi)+diff(e_out,po),Trace+c.trace+'\n')
                                if S != None:
                                    return S
    if One:
        return aux_grey_uno(Parcours[1:],Blacklist_EC,Blacklist_M,Greylist_M,'')
    else :
        if Grey :
            aux_grey(Parcours[1:],Blacklist_EC,Blacklist_M,Greylist_M)
        else :
            aux(Parcours[1:],Blacklist_EC,Blacklist_M)
    return S


def f_with_inserts(CRNabs,Library,Starter,Depth,log=''):
    if Depth == 0 :
        S = f(CRNabs,Library,Starter,Grey = True, One=True)
        return S
    else :
        for e in CRNabs.abstract[0,1:] :
            print(log+'\t'+e+'\t'+'i')
            iCRNabs = simple_insert(CRNabs,e)
            fCRNabs = fork_insert(CRNabs,e)
            Si = f(iCRNabs,Library,Starter,Grey = True, One=True)
            if Si != None :
                return Si
            else :
                Sii = f_with_inserts(iCRNabs,Library,Starter,Depth-1,log+'i')
                if Sii != None:
                    return Sii
            print(log+'\t'+e+'\t'+'f')
            Sf = f(fCRNabs,Library,Starter,Grey = True, One=True)
            if Sf != None :
                return Sf
            else :
                Sff = f_with_inserts(fCRNabs,Library,Starter,Depth-1,log+'f')
                if Sff != None:
                    return Sff



## Loading everything neatly from BRENDA

Lecture = from_brenda_to_clean_lecture("brendaReaction.tsv","brendaReactionBadEntries.tsv")
Lex1,Lex2 = Lexiques_building("Lexique2.tsv")
Unknowns = Lex2['000']
NumberedLecture = Numbering(Lecture,Lex1)
Bibliotheque = Concrete_Reaction_Library(NumberedLecture)
# Should say "Entries read 304848 Entries rejected : 129482"
#D = dictionnary_lengths(Bibliotheque.select_candidates(CRN_abs))

### Test space

liste_abs = [['-','E','F','G','H'],['A','in','-','-','-'],['B','out','in','in','-'],['C','out','in','-','-'],['D','out','out','-','-'],['I','inh','-','-','-'],['J','out','-','-','inh'],['K','-','out','in','-'],['L','-','out','-','out'],['M','-','-','out','in'],['N','-','-','out','in'],['O','-','-','-','out']]

CRN_abs2 = np.array([['-','E','F','G','H'],['A','in','-','-','-'],['B','in','out','-','-'],['C','out','-','in','-'],['D','out','-','in','-'],['I','inh','out','-','-'],['J','-','-','-','in'],['K','-','-','out','-'],['L','in','-','-','out'],['M','-','in','-','out']])


CRN_REFERENCE = CRNa([
['-','W','X','Y','Z'],
['A','-','in','out','-'],
['B','in','-','-','-'],
['C','in','out','-','-'],
['D','out','in','-','-'],
['E','out','-','-','-'],
['G','-','out','in','-'],
['H','-','-','out','in'],
['I','-','-','-','in'],
['J','-','-','-','out'],
])

CRN_REF_LONG = CRNa([
['-','W','X','Y','Z','ZZ'],
['A','-','in','out','-','-'],
['B','in','-','-','-','-'],
['C','in','out','-','-','-'],
['D','out','in','-','-','-'],
['E','out','-','-','-','-'],
['G','-','out','in','-','-'],
['H','-','-','out','in','-'],
['I','-','-','-','in','-'],
['J','-','-','-','-','out'],
['K','-','-','-','out','in']
])


CRN_abs = CRNa(liste_abs)

CRN_test = CRNa([['-','E'],['A','in'],['B','out']])

CRN_REF_COURT =CRNa([
['-','W','X','Y'],
['A','-','in','out'],
['B','in','-','-'],
['C','in','out','-'],
['D','out','in','-'],
['E','out','-','-'],
['G','-','out','in'],
['H','-','-','out'],
])

CRN_TESTFORK =CRNa([
['-','E'],
['A','in'],
['B','in'],
['C','in'],
['D','out'],
['E','out']
])
CRN_very_short=CRNa([
['-','Z','ZZ'],
['I','in','-'],
['J','-','out'],
['H','in','-'],
['K','out','in']
])

CRN_goal=CRNa([
['-','E'],
['A','in'],
['B','in'],
['J','in']])

a = Lex1['Acetone']
b = Lex1['D-glucose']
h = Lex1['H2O2']
j = Lex1['resorufin']
i = Lex1['resazurin']

Starter = {'A': a, 'B': b, 'J': j}


#f_with_inserts(CRN_goal,Bibliotheque,Starter,4)
#f(CRN_REF_LONG,Bibliotheque,Starter,Grey=True,One=True)


### Reading solutionshort to find a good subset of reactions to include for sepi search
import ast
crypte = open("SolutionsShort.tsv",'r')
Texte = crypte.readlines()
crypte.close()
J = 5
CleanedT,CleanedM = set(),set()
for t in Texte[:J]:
    l,m = '',''
    b,c,d,e = False,False,False,False
    for s in t:
        if s=='{':
            b=True
        if b and not(c):
            l+=s
        if c and not(d):
            if s=='{':
                e=True
            if e and not(d):
                m+=s
            if s=='}':
                d=True
        if s=='}':
            c=True

    d1,d2=ast.literal_eval(l),ast.literal_eval(m)
    CleanedT = CleanedT.union(set(d1.values()))
    CleanedM = CleanedM.union(set(d2.values()))

CleanedM = CleanedM.difference({'1','2'})
SetR = set()
i=0

def does_intersect(set1,list1,list2):
    set2 = set(list1).union(set(list2))
    set3 = set1.intersection(set2)
    return (len(set3) > 0)

for n in NumberedLecture :
    if n[1] in CleanedT and does_intersect(CleanedM,n[3],n[4]) :
        s = 'm'+n[3][0]
        for r in n[3][1:]:
            s+=' + m'+r
        s+=' => m'+n[4][0]
        for p in n[4][1:]:
            s+=' + m'+p
        s+='.'
        if not(s in SetR):
            i+=1
            print(i,n[1])

            crypte = open('BrendaBiochamReactionsSmart'+str(J)+'SubsetM.bc','a')
            crypte.write(s+'\n')
            crypte.close
            SetR = SetR.union({s})

## Answer to Mathieu's question on BRENDA's types of reactions
Glossaire = []
j=0
for l in Lecture :
    g = '('+str(len(l[3]))+','+str(len(l[4]))+')'
    if g == '(2,2)':
        if (l[3][0] in l[4]) or (l[3][1] in l[4]):
            j+=1
    Glossaire.append(g)

values,counts = np.unique(Glossaire,return_counts=True)
print('(nbre entrées,nbre sorties) -> nbre de réaction de ce type dans Brenda nettoyé')
for i in range(len(values)):
    v,c = values[i],counts[i]
    print(v+' -> '+str(c)+' '+str(100*c/len(Glossaire))+'%')


### Alternate cleaning method for BRENDA
i,j=0,0
for l in Lecture[1:] :
    i+=1
    metas = l[3]+l[4]
    if not(disjoint(metas,Unknowns)):
        j+=1
        crypte=open('brendaReactionBadEntries.tsv','a')
        crypte.write(str(i)+'\n')
        crypte.close()
    s = "\r Entries read "+str(i)+" Entries rejected : "+str(j)
    print(s, end="")


### Translating BRENDA in Biocham

SetR = set()
for n in NumberedLecture :
    s = 'm'+n[3][0]
    for r in n[3][1:]:
        s+=' + m'+r
    s+=' => m'+n[4][0]
    for p in n[4][1:]:
        s+=' + m'+p
    s+='.'
    if not(s in SetR):
        #print(s)
        crypte = open("BrendaBiochamReactionsM.bc",'a')
        crypte.write(s+'\n')
        crypte.close
        SetR = SetR.union({s})


### LOGICS

import itertools
import pandas as pd
from pprint import pprint
import numpy as np

def AND(list):
    if False in list :
        return False
    if "epsilon" in list :
        return "epsilon"
    return True

def OR(list):
    if True in list :
        return True
    if "epsilon" in list :
        return "epsilon"
    return False

def from_biocham_to_CRN(S):
    L = S.split('\n')
    L = [l.strip() for l in L]
    L2 = []
    LM = []
    for s in L :
        s2 = s.strip('.')
        l = s2.split('=>')
        l[0],l[1] = l[0].split('+'),l[1].split('+')
        l[0],l[1] = [i.strip() for i in l[0]],[i.strip() for i in l[1]]
        Ei = dict()
        for j in l[0]:
            Ei[j]='in'
        for j in l[1]:
            Ei[j]='out'
        for j in list(set(l[0]).intersection(set(l[1]))) :
            Ei[j]='mid'
        Keys = list(Ei.keys())
        LM = LM + Keys
        L2.append(Ei)
    LM = list(set(LM))
    LM.sort()
    DM = dict()
    for i in range(len(LM)):
        DM[LM[i]] = i
    M = np.zeros((len(LM)+1,len(L)+1),dtype="<U100")
    M[1:,0] = LM
    M[0,0] = '-'
    for i,j in itertools.product(range(len(LM)),range(len(L))):
        M[i+1,j+1]='-'
    for j in range(len(L)):
        M[0,j+1]='E'+str(j+1)
        dm = L2[j]
        for k in dm.keys():
            M[DM[k]+1,j+1] = dm[k]
    return M

def deltaE(j,CRN,meta_prod):
    R = list(set(np.where(CRN[1:,j+1]=='in')[0]).union(set(np.where(CRN[1:,j+1]=='mid')[0])))
    return AND([meta_prod[i] for i in R])

def next_state(CRN,meta_prod,meta_conso):
    DE = []
    for e in range(len(CRN[0])-1):
        DE.append(deltaE(e,CRN,meta_prod))

    meta_prod2,meta_conso2 = meta_prod.copy(),meta_conso.copy()
    for e in range(len(CRN[0])-1):
        Emid = set(np.where(CRN[1:,e+1]=='mid')[0])
        Ein = set(np.where(CRN[1:,e+1]=='in')[0])
        Eout = set(np.where(CRN[1:,e+1]=='out')[0])
        Eall = (Emid.union(Ein)).union(Eout)
        de = (DE[e]==True)
        dmid = (len(Emid)>0)
        Aeps = [meta_prod[j] for j in Eall]
        Lepsilon = []
        for j in range(len(Aeps)) :
            if Aeps[j]=='epsilon':
                Lepsilon.append(j)
        depsilon = (len(list(Lepsilon))>0)
        if de and dmid and depsilon :
            for b in Ein :
                meta_prod2[b] = False
                meta_conso2[b] = True
            for a in Emid :
                meta_prod2[a] = True
                meta_conso2[a]= False
            for c in Eout:
                meta_prod2[c] = True
                meta_conso2[c] = False
    for a in range(len(meta_prod)):
        Rin = list(set(np.where(CRN[a+1,1:]=='in')[0]))
        Rout = list(np.where(CRN[a+1,1:]=='out')[0])
        if meta_prod[a] =='epsilon' and len(Rout)>0 :
            meta_prod2[a] = OR([DE[i] for i in Rout])
        conso_a =  OR([DE[i] for i in Rin])
        if conso_a != 'epsilon' :
            meta_conso2[a] = conso_a
    return meta_prod2, meta_conso2

def boolrange(Aprod,Aconso):
    if Aprod=='epsilon' or Aconso=='epsilon':
        return 'epsilon'
    if not(Aprod):
        return False
    if not(Aconso):
        return True
    return 'eta'

def boole(s):
    if s=='True':
        return True
    if s=='False':
        return False

def clean_list_input(l,boolean=False):
    L = l.strip('[')
    L = L.strip(']')
    L = L.split(',')
    if boolean :
        L = [boole(i.strip()) for i in L]
    else :
        L = [i.strip() for i in L]
    return L

def boolt(A):
    if A==True :
        return 1
    if A=='epsilon':
        return '\u03B5'
    if A=='eta':
        return '\u03B7'
    return 0

def boolisation(CRN,inputs,constraint='in',outputs='all',details = True):
    Conso_init, Prod_init = [], []
    DM = dict()
    DS = dict()
    for a in inputs :
        DS[a]=[]
    for i in range(len(CRN)-1):
        DM[CRN[i+1,0]] = i
        DS[CRN[i+1,0]+"-prod"] = []
        DS[CRN[i+1,0]+"-conso"] = []
        DS[CRN[i+1,0]] = []
    if constraint == 'in':
        for a in range(len(CRN)-1):
            Rout = list(np.where(CRN[a+1,1:]=='out')[0])
            if len(Rout)==0:
                Prod_init.append(True)
            else :
                Prod_init.append('epsilon')
            Conso_init.append(False)
    else :
        Prod_init, Conso_init = ['epsilon' for _ in CRN[1:,0]], [False for _ in CRN[1:,0]]
        const_meta, const_const = constraint
        for i in range(len(const_meta)) :
            Prod_init[DM[const_meta[i]]] = const_const[i]

    for input_const in itertools.product([False,True], repeat=len(inputs)):
        conso,prod = Conso_init.copy(),Prod_init.copy()
        for i in range(len(input_const)):
            prod[DM[inputs[i]]] = input_const[i]
        p,c=[],[]
        while p+c != prod+conso :
            p2,c2 = next_state(CRN,prod,conso)
            p,c,prod,conso = prod,conso,p2,c2
        bol = [boolrange(p[i],c[i]) for i in range(len(p))]
        for k in DM.keys():
            DS[k].append(boolt(bol[DM[k]]))
            DS[k+'-prod'].append(boolt(p[DM[k]]))
            DS[k+'-conso'].append(boolt(c[DM[k]]))

    DC = dict()
    for i in inputs:
        DC[i]=DS[i+'-prod']

    if outputs == 'all' :
        out = set(CRN[1:,0]).difference(set(inputs))
    else :
        out = set(outputs)
    if details :
        for o in out:
            DC[o+'-prod'] = DS[o+'-prod']
            DC[o+'-conso'] = DS[o+'-conso']
    for o in out:
        DC[o] = DS[o]
    return pd.DataFrame(DC)


t,f,e = True, False,'epsilon'

sss = "a+b=> c+d.\n e+ f  =>j+a."
ss = 'a+b=>c.'
s = 'a=>d. \n d+b =>c.'
S = 'b+c => d+e. \n a+d => c+g. \n g => a+h. \n h+i => j.'

aCRN = from_biocham_to_CRN(sss)
aCRNeasy = from_biocham_to_CRN(ss)
aCRNArtifice = from_biocham_to_CRN(s)
aCRNShort = from_biocham_to_CRN(S)
aCRNneg = from_biocham_to_CRN('A+B => A+C.')


print('Pour CRNShort :\n')
print(boolisation(aCRNShort,['a','b'], constraint=(['c','i'],[t,t]), outputs=['d','j'], details=True))

print('\n Pour CRNArtifice :\n')
print(boolisation(aCRNArtifice,['a','b']))

print('\n Pour CRNneg :\n')
print(boolisation(aCRNneg,['A'],outputs='all',details=True))