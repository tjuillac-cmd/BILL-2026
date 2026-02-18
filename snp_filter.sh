#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fichier.vcf> <seuil_GQ> <sortie.vcf>"
    echo "Exemple: $0 P25-4.vcf 30 P25-4_filtered.vcf"
    exit 1
fi

VCF=$1
SEUIL=$2
OUTPUT=$3

(grep '^#' "$VCF"; grep -v '^#' "$VCF" | awk -F'\t' -v seuil="$SEUIL" '{
    n = split($NF, a, ":")
    gq = a[n]
    if (gq+0 >= seuil) print
}') > "$OUTPUT"

echo "Filtrage terminé : $(grep -v '^#' "$OUTPUT" | wc -l) variants conservés (GQ >= $SEUIL)"