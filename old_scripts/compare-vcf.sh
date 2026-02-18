#!/bin/bash                                                                     

# Vérification des arguments                                                    
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 fichier1.vcf fichier2.vcf"
    exit 1
fi

VCF1=$1
VCF2=$2

# Nom de sortie basé sur les noms des fichiers d'entrée                         
OUT="diff_$(basename $VCF1 .vcf)_vs_$(basename $VCF2 .vcf)"

vcftools --vcf "$VCF1" --diff "$VCF2" --diff-site --out "$OUT"

echo "Fichier de sortie : ${OUT}.diff.sites_in_files"
