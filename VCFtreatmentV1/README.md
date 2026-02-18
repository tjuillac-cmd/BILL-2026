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
python3 vcfTreatmentV1.py --input P15-4.vcf P25-4.vcf --out comparatif.vcf --merge --venn --statistics
```

### 3. Analyse Complexe (4 échantillons et plus)

Fusionne les fichiers et bascule automatiquement sur un UpSet Plot pour une meilleure lisibilité des intersections complexes.

```bash
python3 vcfTreatmentV1.py --input P15.vcf P25.vcf P50.vcf P90.vcf --out global.vcf --merge --venn
```

---

## Options Disponibles

| Option | Description |
|---|---|
| `--input` | Liste des fichiers VCF d'entrée (espace entre chaque fichier) |
| `--out` | Nom du fichier VCF de sortie (préfixe pour les graphiques) |
| `--merge` | Fusionne les entrées en ajoutant un tag `ORIGIN` dans le champ INFO |
| `--annotate` | Parse le champ `ANN` (SnpEff) pour les graphiques et le résumé |
| `--statistics` | Calcule la VAF moyenne et l'Entropie de Shannon |
| `--venn` | Génère un diagramme de Venn (2-3 éch.) ou UpSet Plot (4+ éch.) |
| `--count` | Affiche simplement le nombre de variants traités dans le terminal |

---

## Fichiers en sortie

Le programme utilise le nom du fichier `--out` comme préfixe pour générer :

1. `prefix.vcf` : Le fichier de données filtré/fusionné.
2. `prefix_summary.txt` : Rapport textuel complet (Stats, top gènes, détails HIGH/MODERATE).
3. `prefix_impacts_pie.png` : Répartition visuelle des impacts SnpEff.
4. `prefix_venn.png` ou `prefix_upset.png` : Comparaison des échantillons.

---

## Notes Techniques (Reprise du travail)

- **Erreur NoneType** : Le script gère désormais les valeurs `GQ` ou `DV` absentes (remplacées par `0.0`) pour éviter les crashs sur les fichiers mergés.
- **Champs INFO** : Pour éviter les erreurs vcfpy `Expected list value`, le champ `ORIGIN` est toujours écrit sous forme de liste : `record.INFO['ORIGIN'] = [nom]`.
- **Calcul VAF** : Si `DV` est absent, le script tente de le récupérer via le champ `AD` ou met `0` par défaut.
- **VS Code** : Si les imports sont soulignés en rouge, vérifie que l'interpréteur Python (en bas à droite) est bien sur `viral_vcf_env`.