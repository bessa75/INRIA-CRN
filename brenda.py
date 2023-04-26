#@title get data
import pandas as pd
import json
import time
from functools import reduce
from tqdm import tqdm

link_to_data = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQW5Udu9IvmmRpJdl4GfCGhy0ZEq-kNhKIuo1bGpQUpYchPNmDdYjm846DmKRB6UVWjkIgCXTO_ChiV/pub?output=csv'

class DisJointSets():
   def __init__(self,N):
      # Initially, all elements are single element subsets
      self._parents = [node for node in range(N)]
      self._ranks = [1 for _ in range(N)]
      self._reachable = [set([node]) for node in range(N)]

   def find(self, u):
      while u != self._parents[u]:
         # path compression technique
         self._parents[u] = self._parents[self._parents[u]]
         u = self._parents[u]
      return u

   def connected(self, u, v):
      return self.find(u) == self.find(v)

   def _union(self, u, v):
      # Union by rank optimization
      root_u, root_v = self.find(u), self.find(v)
      if root_u == root_v:
         return True
      if self._ranks[root_u] > self._ranks[root_v]:
         self._parents[root_v] = root_u
      elif self._ranks[root_v] > self._ranks[root_u]:
         self._parents[root_u] = root_v
      else:
         self._parents[root_u] = root_v
         self._ranks[root_v] += 1
      return False
  
   def union(self, *args):
      for arg in args[1:]:
         self._union(args[0], arg)

   def nbr_of_sets(self):
      return len(set(self._parents))
class Reachable:
   def __init__(self, size):
      self._parents = list(range(size))
      self._reachable = [set() for _ in range(size)]
      self._rank = [0 for _ in range(size)]
   # Ici on remonte jusqu'au parent en ajoutant à chaque fois la liste des produits dans produits atteignables
   def find(self, u, produits):
      while u != self._parents[u]:
         # path compression technique
         self._parents[u] = self._parents[self._parents[u]]
         u = self._parents[u]
      return u

   def union(self, produits, *args):
      parent = self.find(args[0], produits)
      for arg in args[1:]:
         self._parents[self.find(arg, produits)] = parent


def write_data(fichier):  
  print("Reading data from distant database...")
  tic = time.perf_counter()
  with tqdm(total=100, desc="Downloading dataset") as pbar:
    dataframe = pd.read_csv(link_to_data)
    pbar.update(100)
  tac = time.perf_counter()
  print(f"\nDone in {tac - tic} seconds.\nBuilding dictionnaries...") #12seconds
  tic = time.perf_counter()
  data = dataframe[['Substrats', 'Produits', 'KM', 'KCAT', 'PH']].drop_duplicates(subset=['Substrats', 'Produits'])

  # strip for removing spaces before spliting
  data['reactifs']  = data['Substrats'].map(lambda v: v.strip(' ').split('-+-'))
  data['produits']  = data['Produits' ].map(lambda v: v.strip(' ').split('-+-'))
  data['kms'      ] = data['KM'       ].map(lambda v: v.strip(' ').split(';'))
  data['kcats'    ] = data['KCAT'     ].map(lambda v: v.strip(' ').split(';'))
  data['phs'      ] = data['PH'       ].map(lambda v: v.strip(' ').split('-'))

  # on retient ce qui nous interesse
  useful = data[['reactifs', 'produits', 'kms', 'kcats', 'phs']]

  #identification des différentes especes
  molecules = set()
  for r in useful['reactifs'] + useful['produits']:
    molecules.update([e for e in r if e != ''])

  # équivalence nom - index  (dictionnaire/lexique)
  names_to_idx = {mol : idx for idx, mol in enumerate(molecules)}
  # reciproque
  idx_to_name  = {idx : mol for idx, mol in enumerate(molecules)}

  print (f'\nDone in {time.perf_counter() - tic} seconds.\nData processing {fichier}...')
  tic = time.perf_counter()

  # UnionFind pour identifier les molécules qui sont liées par des réactions
  clusters = DisJointSets(len(molecules))
  reachable = Reachable(len(molecules))
  # encodage
  """
  reactions = [
    {
      [liste des reactifs],
      [liste des produits],
      [kms],
      [kcats].
      [phs]  
    } autant de fois qu'il y a de réactions
  ]
  """
  reactions = []
  """
  reactions_per_mol = [
    {
      [liste des index des réactions dans lesquelles la molecule réagit]
    } autant de fois qu'il y a de molécules
  ]
  """
  reactions_per_mol = [ set() for _ in range(len(molecules))]
  reaction_per_mol_p = [ set() for _ in range(len(molecules))]
  add = reactions.append
  counter = 0
  inconsistant = 0
  with tqdm(total = len(useful.index)) as pbar:
    for _, d in useful.iterrows():
      
      kms      = [km   if km   != '?' else -1 for km   in d.kms   ] 
      kcats    = [kcat if kcat != '?' else -1 for kcat in d.kcats ]
      
      if d.phs == ['?'] :   phs = [-1, -1]
      elif len(d.phs) == 1: phs = [d.phs[0], d.phs[0]]
      else :                phs = [d.phs[0], d.phs[1]]

      reactifs = [names_to_idx[reac] for reac in d.reactifs if reac !=  '']
      produits = [names_to_idx[prod] for prod in d.produits if prod != '']
      
      # update progress bar
      pbar.update(1)

      # Normalement pas besoin de faire ça, mais on sait jamais
      # S'il n'y a pas de réactifs ou de produits ça veut tout simplement dire que la réaction n'est pas conistante.
      if len(reactifs) == 0 or len(produits) == 0: 
        inconsistant += 1
        continue
      # Sinon on peut ajouter la réaction

      add([reactifs, produits, kms, kcats, phs])
      for r in reactifs:
          reactions_per_mol[r].add(counter)
      for p in produits:
          reaction_per_mol_p[p].add(counter)

      
      # Unions des molécules qui sont liées par cette réaction
      clusters.union(*reactifs)
      reachable.union(produits, *reactifs)
      
      counter += 1

      

  print(f'\nDone in {time.perf_counter() - tic} seconds.\n Saving data...')
  tic = time.perf_counter()
  # Dictionnaire rassemblant les informations
  dictionnaire =  {
    'dic_name_to_idx': names_to_idx, 
    'dic_idx_to_name': idx_to_name, 
    'molecules_list': list(molecules), 
    'reactions': reactions, 
    'reactions_per_molecule': [list(rpm) for rpm in reactions_per_mol],
    'reactions_per_molecule_p': [list(rpm) for rpm in reaction_per_mol_p],
    'cluster': clusters._parents,
    'reachable': [list(reach) for reach in reachable._reachable]
  }

  with open(fichier, 'w') as f: json.dump(dictionnaire, f)
  
  print(f'\nDone in {time.perf_counter() - tic} seconds.\n')
  return dictionnaire

def stats(fichier):
  data = read_data(fichier)
  # Nombre maximum de reactivité
  molecularities_ = [ (i, len(l)) for i, l in enumerate(data['reactions_per_molecule']) ]

  molecularities = sorted(molecularities_, key = lambda x: x[1], reverse=True)
  molecules_correspondantes = [data['dic_idx_to_name'][str(i)] for i, _ in molecularities]
  reactifs = sorted([(i, len(r[0])) for i, r in enumerate(data['reactions'])], key = lambda x : x[1], reverse=True)

  return molecularities, molecules_correspondantes, reactifs


def read_data(fichier):
  with open(fichier, 'r') as f:
    dictionnaire = json.load(f)
  return dictionnaire

if __name__ == '__main__':
  write_data('data.json')
