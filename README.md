BOUAOUD Malik


#### Contenue Du dépôt :
```
├── Diapo
│   └── diapo.pdf
├── README.md
├── data
│   ├── MA0056.1.jaspar
│   ├── MA0057.1.jaspar
│   ├── MA0083.1.jaspar
│   ├── MA0114.1.jaspar
│   ├── MA0114.3.jaspar
│   ├── jaspar.txt
│   ├── sequence.fasta
│   └── sequence.gb
├── src
│   ├── putative_TFBS.py
│   ├── pwm.py
│   ├── scan_pwm.py
│   ├── tp2_1.py
│   ├── tp2_2.py
│   └── utils.py
└── testDir
    └── test.py
```
## Lancer le programme :

##### REMARQUE IMPORTANTE 
Les commandes suivantes se lance dans le  repertoire **src**, si vous voullez lancer les test
vous devez suivre les étapes décrites ci-dessous.

### Ensuite plusieurs options s'offre à vous:

#### Lancer putative_TFBS avec FULL parametrage:
>python3 putative_TFBS.py  <List id GenBank> -m <matrice jaspar> -t <sueil score > -l <taille du Promoteur> -w <taille fenetre> -s  <seuil fentre> -p <pseudo pods>

##### Exemple :
>  python3 putative_TFBS.py "NM_007389" "NM_079420" "NM_001267550" "NM_002470" "NM_003279" "NM_005159" -m MA0056.1.jaspar -t -2 -l 1000 -w 100 -s  -2 -p 0.1

#### Lancer putative_TFBS Not Full parametrage :

#### Exemple

> python3 putative_TFBS.py "NM_007389" "NM_079420" "NM_001267550" "NM_002470" "NM_003279" "NM_005159" -m MA0056.1.jaspar -t -2 -l 1000  -s  -2 

#### Lancer scan_pwm
> python3 scan_pwm.py MA0056.1.jaspar NM_007389 1000 -2
#### Lancer les Tests unitaires :
Se placer dans la racine 

##### Lancer  
>export PYTHONPATH=.
#### Puis
>python3 testDir/test.py
