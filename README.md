# PSC

Git privé pour échanger les fichiers.

Branches "main" et "edit-in-progress"


## TODO

Clarifier affichage mécanismes trouvés
Vérifier fonctionnement du NON
Implémenter le OU

## Ajouter des  files

- [ ] [Créer](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) ou [Téléverser](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) fichiers
- [ ] [Ajouter des fichiers en ligne de commande](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) ou push un dossier Git existant avec la commande:

```
cd existing_repo
git remote add origin https://gitlab.binets.fr/psc-compiling-math-functions-in-biochemical-reactions/psc.git
git branch -M main
git push -uf origin main
```

***
