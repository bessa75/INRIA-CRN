import pandas as pd

reactions = pd.read_csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vQW5Udu9IvmmRpJdl4GfCGhy0ZEq-kNhKIuo1bGpQUpYchPNmDdYjm846DmKRB6UVWjkIgCXTO_ChiV/pub?output=csv')

data = reactions[['Substrats', 'Produits']].drop_duplicates()

data['reactifs'] = data['Substrats'].map(lambda v: v.split('-+-'))
data['produits'] = data['Produits'].map(lambda v: v.split('-+-'))

simple_reactions = data[['reactifs', 'produits']]

molecules = set()
for r in simple_reactions['reactifs']:
  molecules.update(r)
for p in simple_reactions['produits']:
  molecules.update(p)


# équivalence nom - index
names_to_idx = {mol : idx for idx, mol in enumerate(molecules)}

molecules_encodage = simple_reactions.apply(
    lambda reaction: pd.Series([[names_to_idx[reactif] for reactif in reaction.reactifs], [names_to_idx[produit] for produit in reaction.produits]], index = ['reactifs', 'produits']), 
    axis=1
)

mol_enc_to_list1 = molecules_encodage.reactifs.tolist()
mol_enc_to_list2 = molecules_encodage.produits.tolist()

# [[reactifs, produits]]
mol_enc_to_list = [[u, v] for u, v in zip(mol_enc_to_list1, mol_enc_to_list2)]


