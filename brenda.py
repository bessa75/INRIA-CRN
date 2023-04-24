#@title get data
import pandas as pd
import json
import time

def write_data(fichier):
    
  print("Reading data from database...")
  tic = time.perf_counter()
  dataframe = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vQW5Udu9IvmmRpJdl4GfCGhy0ZEq-kNhKIuo1bGpQUpYchPNmDdYjm846DmKRB6UVWjkIgCXTO_ChiV/pub?output=csv')
  tac = time.perf_counter()
  print(f"Done in {tac - tic} seconds.\nComputing useful data...") #12seconds
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
  idx_to_name = {v : k for k, v in names_to_idx.items()}

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
  add = reactions.append
  counter = 0
  for _, d in useful.iterrows():
    # print(d)
    # reactifs = [r  for r  in d.reactifs if r  != '']

    # produits = [p for p    in d.produits if r    != '']
    kms      = [km   if km   != '?' else -1 for km   in d.kms   ] 
    kcats    = [kcat if kcat != '?' else -1 for kcat in d.kcats ]
    
    if d.phs == ['?'] :   phs = [-1, -1]
    elif len(d.phs) == 1: phs = [d.phs[0], d.phs[0]]
    else :                phs = [d.phs[0], d.phs[1]]
    reactifs = [names_to_idx[reac] for reac in d.reactifs if reac !=  '']
    produits = [names_to_idx[prod] for prod in d.produits if prod != '']
    # Normalement pas besoin de faire ça, mais on sait jamais
    # S'il n'y a pas de réactifs ou de produits ça veut tout simplement dire que la réaction n'est pas conistante.
    if len(reactifs) == 0 or len(produits) == 0: continue
    # Sinon on peut ajouter la réaction
    add([reactifs, produits, kms, kcats, phs])
    for r in d.reactifs:
        reactions_per_mol[names_to_idx[r]].add(counter)
    counter += 1
  # Dictionnaire rassemblant les informations
  dictionnaire =  {
    'dic_name_to_idx': names_to_idx, 
    'dic_idx_to_name': idx_to_name, 
    'molecules_list': list(molecules), 
    'reactions': reactions, 
    'reactions_per_molecule': [list(rpm) for rpm in reactions_per_mol],
    'nombre_de_reactions': counter
  }
  with open(fichier, 'w') as f: json.dump(dictionnaire, f)
  
  return dictionnaire


def read_data(fichier):
  with open(fichier, 'r') as f:
    dictionnaire = json.load(f)
  return dictionnaire

if __name__ == '__main__':
  write_data('data.json')
