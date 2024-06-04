# Detection automatique d'interactions sucres-protéiones dans les séquences de protéomes viraux et bactériens

Le but est d'appliquer des profils HMM sur des protéomes viraux et bactériens afin de tenter de caractériser des protéines impliquées dans les interactions sucres-protéines. L'application manuelle de profils est lente et fastidieuse, de même que la transformation des résultats en graphiques. C'est pourquoi nous automatisons ces processus grâce à un programme python. 


## Fonctionnement

-   Main.py : Coordonne l'appel de tous les modules

-	Profiles.py : Ce module gère l'exécution des recherches HMMER, utilisant les profils préétablis pour identifier les séquences correspondantes dans les protéomes.

-	Plotting.py : Après l'analyse, ce module génère une HeatMap par modèle qui visualise la distribution des hits ainsi que leur e-value respective. Il produit également des diagrammes représentant la localisation génomique de chaque protéine identifiée comme un hit pour un modèle donné.

-	Matrix.py : Il produit des matrices d'occurrences qui résument les résultats de l'analyse, offrant une vue d'ensemble de la présence relative des familles de protéines dans les échantillons étudiés.

-	Hits_list.py : Répertorie dans un fichier texte tous les profils ayant enregistré au moins un hit. Dans un autre fichier, il liste tous ceux qui n'ont présenté aucune occurrence.  

- TaxonTable.py : Crée un table d'occurrence pour qui traduit les numéros des protéomes en nom d'espèce


## Utilisation


<details>

<summary>Plusieurs librairies doivent être installées avant l'utilisation </summary>

<br>
-	Os<br>
-	Sys<br>
-	Shutil<br>
-	Glob<br>
-	Argparse<br>
-	Time<br>
-	Pandas<br>
-	Openpyxl.utils<br>
-	Bio<br>
-	Collections<br>
-	Matplotlib<br>
-	Seaborn<br>
-	Numpy<br>
-	Subprocess<br>
-	Csv<br>
-	Multiprocessing<br>
</details>
<br>

Ce programme est prévu pour être utiliser sur les protéomes téléchargés depuis la base de donnée du NCBI, car il dépend d'une architecture précise de dossiers qui doit être : 
```bash
data/<Banque_données>_<ID_Protéome>/
```
Chaque dossier <Banque_données>_<ID_Protéome> doit **comprendre les fichiers genomic.gbff et protein.faa** propre à la souche.<br>

Pour executer :
  
```bash
python3 Main.py Proteome_path Model_path Bank_name --evalue 
```

- Proteome_path : Chemin d'accès au répertoire du protéome qui doit pointer sur data/ (str).

- Model_path : Chemin d'accès au répertoire contenant les modèles (str).

- Bank : Banque à utiliser pour l'analyse (str).

- Evalue : Seuil de l'evalue minimum pour les correspondances (float), avec une valeur par défaut de 0,012.


Si on choisit de télécharger les données des banques **RefSeq et GenBank du NCBI en même temps**, on trouve dans le dossier Data/ les versions d’un protéome à la fois en GCA et en GCF.<br>

**Par exemple**, pour le protéome 000009685.1 de Clostridium perfringens :<br>
●	Version **GCA** : <br>
Clostridium_Perfringens_Complete_Genome/ncbi_dataset/data/GCA_000009685.1<br>
●	Version **GCF** :<br>
Clostridium_Perfringens_Complete_Genome/ncbi_dataset/data/GCF_000009685.1<br>

L'argument **Bank** est utilisé pour spécifier quelle banque de données analyser, afin d’éviter la redondance des résultats. Si on souhaite utiliser RefSeq ou GenBank nous avons simplement à écrire GCF ou GCA. Il est possible d'utiliser une autre banque tant que tous les dossiers commencent par le même prefix exemple : data_<ID_protéome>/, donc Bank_name = "data".

## Résultats 
Les résultats s'organisent selon cette architecture : 
```
.
└── Results_<Organism_name>*_<Bank>_*<Evalue>/
    ├── <Model_name_1>/
    │   ├── All_Profiles/
    │   │   └── <Model_Name>_<Assembly>_Profile.hmmer
    │   ├── Hits_Profiles/
    │   │   └── <Model_Name>_<Assembly>_Profile.hmmer
    │   ├── Hits_Profiles_Filtered/
    │   │   └── <Model_Name>_<Assembly>_Profile.hmmer
    │   ├── <Model_Name>_<Bank>_cluster
    │   ├── <Model_Name>_<Bank>_cluster.clstr
    │   ├── <Model_Name>_<Bank>_Fasta_file.fasta
    │   ├── <Model_Name>_<Bank>_Genomic_Location_Plot.png
    │   ├── <Model_Name>_<Bank>_HeatMap.png
    │   ├── <Model_Name>_<Bank>_Mean_Evalue_Matrix.csv
    │   └── <Model_Name>_Filtered_Sequence_<Bank>.csv
    ├── <Model_name_2>
    ├── <Taxon_Name>_<Bank>_Hits_matrix.xlsx
    ├── <Taxon_Name>_With_Hits.txt
    └── <Taxon_Name>_Without_Hits.txt
```
