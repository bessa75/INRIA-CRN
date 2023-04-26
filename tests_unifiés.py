import recherche_chemin
import brenda
import pandas as pd
import time

def test_recherche_chemin(_params):
   
   ENZ, REAC, prod, n, imprime, liste_reaction_texte, MOL, REACTIONS, REACPARMOL, bool = _params
   # print(REACTIONS[:20])
   # print(REACTIONS[:20])
   # recherche_chemin.pluscourtchemin(ENZ,REAC,prod,SETM,SETR,n,imprime,liste_reaction_texte,MOL,REACTIONS,REACPARMOL,bool=False,bool2=False)
   return recherche_chemin.pluscourtchemin(ENZ, REAC, prod, n, imprime, liste_reaction_texte, MOL, REACTIONS, REACPARMOL, bool)

def brenda_test():
   data = brenda.read_data('data.json')
   # nombre_de_reactions = data['nombre_de_reactions']	
   dic_name_to_idx = data['dic_name_to_idx']
   dic_idx_to_name = data['dic_idx_to_name']
   reactions = data['reactions']
   reactions_per_molecule = data['reactions_per_molecule']
   clusters = data['cluster']

   a = data['dic_name_to_idx']['pyruvate'] #data['dic_name_to_idx']['beta-alanine']

   # print(f'Il y a { data["nombre_de_reactions"] } reactions akip')
   # print(f'Il y a {nombre_de_reactions} réactions dans la base de données.')
   # print(f'Il y a {len(dic_name_to_idx)} molécules dans la base de données.')
   # print(f'On a aussi {len(reactions)} Dans la base de données.')
   # print(f'le plus gros indice dans la liste des reactions est {max(reactions_per_molecule)}')

   count = 0
   friends = []
   add = friends.extend
   h20 = dic_name_to_idx['H2O']
   idx_l_ascorbic_tatata = data['dic_name_to_idx']['L-ascorbic acid alpha-D-glucoside']
   idx_H2O = data['dic_name_to_idx']['H2O']
   idx_maltose = data['dic_name_to_idx']['maltose']
   l = 0
   # print(data['dic_idx_to_name'][str(idx_l_ascorbic_tatata)])
   for r in reactions:
      if a in r[1]: 
         # print(r)
         
         # count = max(count, len(r[0]))
         # if len(r[0]) >= 4:
         #    l += 1
         print(([[data['dic_idx_to_name'][str(i)] for i in r[0]],  [data['dic_idx_to_name'][str(i)] for i in r[1]]]))
         print('-------------------')
      #   if idx_l_ascorbic_tatata in r[0]:
      #     l.append([[data['dic_idx_to_name'][str(i)] for i in r[0]],  [data['dic_idx_to_name'][str(i)] for i in r[1]]])
   
   return count, l


def main():
    print(f'Starting reading data...')
    tic = time.perf_counter()
    data = brenda.read_data('data.json')
    tac = time.perf_counter()
    print(f'Done in {tac - tic} seconds.')
    # Reaction L-ascorbic acid alpha-D-glucoside + H2O	-> maltose + L-ascorbicacid
    idx_l_ascorbic_tatata = data['dic_name_to_idx']['L-ascorbic acid alpha-D-glucoside']
    idx_H2O = data['dic_name_to_idx']['H2O']
    idx_maltose = data['dic_name_to_idx']['maltose']
    #  test = data['dic_name_to_idx'][''] # Pour asavoir s'il n'y a pas de chaines vides dans le dictionnaire.
    idx_l_ascorbic_acid = data['dic_name_to_idx']['L-ascorbicacid']

    print('Finding solutions...')
    tic = time.perf_counter()

    # Test de cohérence des données
    # print(max(data['dic_idx_to_name']))
    print(len(data['dic_name_to_idx']))
    print(len(data['molecules_list']))
    print(f'Done in {tac - tic} seconds.\n')
   
    mecanisme = test_recherche_chemin((
        [], 
        [idx_l_ascorbic_tatata, idx_H2O], 
        idx_maltose, 
        10,
        False,
        [],
        data['dic_idx_to_name'], 
        data['reactions'], 
        data['reactions_per_molecule'],
        False
    ))
    tac = time.perf_counter()
    
    print(mecanisme)
   # brenda_test()


if __name__ == '__main__':
   # main()
   r = brenda_test()
   print(r)
   # print('c =  ', c)
   # print('len = ', len(f))
   # print('cluster ', n)

   # print(len(set(f)))

   
   # print(f[:10])
   # molecularities, molecules = brenda.stats('data.json')
   # print(molecularities[:5])
   # print(molecules[:5])


# phosphoenolpyruvate + HPr -> pyruvate + phospho-HPr
# L-asparate -> beta-alanine + CO2
# beta-alanine + pyruvate -> oxaloacetate + alamine
#