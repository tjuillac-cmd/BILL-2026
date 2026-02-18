# VCF Treatment Tool (v1)

Outil de traitement, d'annotation et de comparaison de fichiers VCF (SNPs et SV) optimisé pour l'étude des populations du virus KHV.

---

## Installation et Dépendances

Assure-toi d'être dans ton environnement Conda `viral_vcf_env`.

```bash
# Installation des bibliothèques nécessaires
pip install vcfpy numpy matplotlib matplotlib-venn upsetplot venny4py
```

---

## Utilisation (Lignes de commande)

### 1. Analyse et Annotation (Un seul échantillon)

Filtre les variants, calcule les statistiques de population (VAF, Entropie) et génère le rapport d'impact SnpEff.

```bash
python3 vcfTreatmentV1.py --input P25-4.vcf --out P25_annotated.vcf --annotate --statistics
```

### 2. Comparaison et Venn (2 à 3 échantillons)

Fusionne les fichiers et crée un diagramme de Venn montrant les variants communs et spécifiques avec transparence.

```bash
python3 vcfTreatmentV1.py --input P15-4.vcf P25-4.vcf --out comparatif.vcf -merge --venn --statistics
```

### 3. Analyse Complexe (4 échantillons et plus)

Fusionne les fichiers et bascule automatiquement sur un UpSet Plot pour une meilleure lisibilité des intersections complexes.

```bash
python3 vcfTreatmentV1.py --input P15.vcf P25.vcf P50.vcf P90.vcf --out global.vcf -merge --venn
```

### 4. Comparaison entre deux origines d'un fichier mergé

Génère un barplot comparatif entre deux échantillons avec, optionnellement, la VAF moyenne en axe secondaire.

```bash
python3 vcfTreatmentV1.py --input merged.vcf --compare P25-4 P25-8 --add-vaf
```

---

## Options Disponibles

### Entrées / Sorties

| Option | Description |
|---|---|
| `-i` / `--input` | Liste des fichiers VCF d'entrée (espace entre chaque fichier). **Obligatoire.** |
| `-o` / `--out` | Nom du fichier VCF de sortie (utilisé aussi comme préfixe pour tous les graphiques générés) |

### Filtrage

| Option | Type | Description |
|---|---|---|
| `--filter-snp` | `float` | Seuil de qualité GQ minimum pour conserver un SNP |
| `--filter-sv` | `int` | Seuil de profondeur DV minimum pour conserver un SV |

### Annotation SnpEff

| Option | Description |
|---|---|
| `--annotate` | Lance l'annotation SnpEff sur le fichier de sortie. Construit automatiquement la base de données si absente. Nécessite `snpEff` installé via Conda. |
| `--genome` | Génome de référence à utiliser pour SnpEff. **Défaut : `KHV_U`** |

### Visualisation

| Option | Description |
|---|---|
| `--distrib-qual` | Histogramme de distribution de la qualité GQ *(fichiers SNP uniquement)* |
| `--distrib-dv` | Histogramme de distribution de la profondeur DV *(fichiers SV uniquement)* |
| `--distrib-vaf` | Histogramme de distribution de la VAF *(fichiers SV uniquement)* |
| `--venn` | Diagramme de Venn pour 2 ou 3 échantillons, basé sur le tag `ORIGIN` |

### Comparaison entre origines

| Option | Description |
|---|---|
| `--compare ORIGIN1 ORIGIN2` | Compare les variants de deux origines dans un fichier mergé : génère un barplot avec variants uniques à chaque origine et variants communs. En cas d'origines invalides, la liste des origines disponibles est affichée. |
| `--add-vaf` | Ajoute la VAF moyenne ± écart-type en axe secondaire sur le graphique `--compare` *(SV uniquement)* |

### Statistiques et Comptage

| Option | Description |
|---|---|
| `--statistics` | Calcule et affiche la VAF moyenne, l'écart-type et l'Entropie de Shannon sur l'ensemble des variants |
| `-count` | Affiche le nombre total de variants traités dans le terminal |

### Fusion

| Option | Description |
|---|---|
| `-merge` | Fusionne plusieurs fichiers VCF en un seul en ajoutant un tag `ORIGIN` dans le champ INFO (valeur = nom du fichier source sans extension) |

---

## Fichiers en sortie

Le programme utilise le nom du fichier `--out` comme préfixe pour générer :

| Fichier | Condition | Description |
|---|---|---|
| `prefix.vcf` | Toujours (si `--out`) | Fichier de données filtré / fusionné |
| `prefix_summary.txt` | `--annotate` | Rapport textuel complet : distribution des impacts, top gènes, détails HIGH/MODERATE |
| `prefix_impacts_pie.png` | `--annotate` | Camembert de répartition des impacts fonctionnels SnpEff |
| `prefix_venn.png` | `--venn` | Diagramme de Venn des intersections entre 2 ou 3 échantillons |
| `prefix_qual.png` | `--distrib-qual` | Histogramme de distribution GQ |
| `prefix_dv.png` | `--distrib-dv` | Histogramme de distribution DV |
| `prefix_vaf.png` | `--distrib-vaf` | Histogramme de distribution VAF |
| `compare_ORIGIN1_ORIGIN2.png` | `--compare` | Barplot comparatif entre deux origines |

---

## Notes Techniques (Reprise du travail)

- **Erreur NoneType** : Le script gère désormais les valeurs `GQ` ou `DV` absentes (remplacées par `0.0`) pour éviter les crashs sur les fichiers mergés.
- **Champs INFO** : Pour éviter les erreurs vcfpy `Expected list value`, le champ `ORIGIN` est toujours écrit sous forme de liste : `record.INFO['ORIGIN'] = [nom]`.
- **Calcul VAF** : Si `DV` est absent, le script tente de le récupérer via le champ `AD` ou met `0` par défaut.
- **Détection SNP / SV** : La distinction se fait automatiquement à la lecture du header (présence de champs `VAF`, `DV` ou de lignes `##ALT=<ID=INS/DEL/...>`), et non plus uniquement sur le nom de fichier.
- **Renommage d'échantillon** : Lors de l'écriture, le nom de l'échantillon est normalisé en `SAMPLE` dans le header et dans chaque record pour garantir la compatibilité du writer vcfpy.
- **VS Code** : Si les imports sont soulignés en rouge, vérifie que l'interpréteur Python (en bas à droite) est bien sur `viral_vcf_env`.