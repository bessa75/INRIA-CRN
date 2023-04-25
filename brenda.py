#@title get data
import pandas as pd
def get_data(df):
    data = df[['Substrats', 'Produits', 'KM', 'KCAT', 'PH']].drop_duplicates(subset=['Substrats', 'Produits'])

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
    for i, d in useful.iterrows():
      # print(d)
      # reactifs = [r  for r  in d.reactifs if r  != '']

      for r in d.reactifs:
          reactions_per_mol[names_to_idx[r]].add(i)

      # produits = [p for p    in d.produits if r    != '']
      kms      = [km   if km   != '?' else -1 for km   in d.kms   ] 
      kcats    = [kcat if kcat != '?' else -1 for kcat in d.kcats ]
      
      if d.phs == ['?'] :   phs = [-1, -1]
      elif len(d.phs) == 1: phs = [d.phs[0], d.phs[0]]
      else :                phs = [d.phs[0], d.phs[1]]

      add([d.reactifs, d.produits, kms, kcats, phs])

    return {
        'dic_name_to_idx': names_to_idx, 
        'dic_idx_to_name': idx_to_name, 
        'molecules_list': list(molecules), 
        'reactions': reactions, 
        'reactions_per_molecule': list(reactions_per_mol)
      }