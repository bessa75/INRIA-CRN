# Projet Scientifique Collectif : détermination de mécanismes réactionnels correspondant à des équations logiques chimiques

HAZERA Gabriel, LI Romain, TCHOMBA-NGUEKO Ivan, DAGONNEAU Thomas, BESSA Swann, MONTOYA Yann


## Fichiers importants
Rapport.pdf : rapport final du projet

main_code :

* algo_resolution_vesicules_simples.py : application de l'algorithme sur 3 exemples de CRNs

* pathFinder.py : résolution de type Breadth First Search à partir d'un produit et d'une liste fixée de réactifs, pour une équation logique du type R1 +R2 +... + Rn = C

* negation.py : résolution de type Breadth First Search d'équations logiques du type A + nonB = C

* algoresolution_système.py : résolution de systèmes d'équations logique par combinaison des ensembles d'enzymes des équations individuelles

* tests_brenda.py : extraction et analyse de la base de données Brenda

```
cd existing_repo
git remote add origin https://gitlab.binets.fr/psc-compiling-math-functions-in-biochemical-reactions/psc.git
git branch -M main
git push -uf origin main
```

***
