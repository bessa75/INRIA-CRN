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
   data = brenda.write_data('data.json')
   nombre_de_reactions = data['nombre_de_reactions']	
   dic_name_to_idx = data['dic_name_to_idx']
   dic_idx_to_name = data['dic_idx_to_name']
   reactions = data['reactions']
   reactions_per_molecule = data['reactions_per_molecule']

   print(f'Il y a { data["nombre_de_reactions"] } reactions akip')
   print(f'Il y a {nombre_de_reactions} réactions dans la base de données.')
   print(f'Il y a {len(dic_name_to_idx)} molécules dans la base de données.')
   print(f'On a aussi {len(reactions)} Dans la base de données.')
   print(f'le plus gros indice dans la liste des reactions est {max(reactions_per_molecule)}')


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
   main()
   # molecularities, molecules = brenda.stats('data.json')
   # print(molecularities[:5])
   # print(molecules[:5])
