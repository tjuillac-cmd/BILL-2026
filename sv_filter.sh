#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fichier.vcf> <seuil_DV> <sortie.vcf>"
    echo "Exemple: $0 P25-4.vcf 100 P25-4_filtered.vcf"
    exit 1
fi

VCF=$1
SEUIL=$2
OUTPUT=$3

(grep '^#' "$VCF"; grep -v '^#' "$VCF" | awk -F'\t' -v seuil="$SEUIL" '{
    # Récupère la position de DV dans le FORMAT
    n = split($9, fmt, ":")
    dv_idx = -1
    for (i = 1; i <= n; i++) {
        if (fmt[i] == "DV") dv_idx = i
    }
    if (dv_idx == -1) next  # pas de champ DV, on saute

    split($10, sample, ":")
    dv = sample[dv_idx]
    if (dv+0 >= seuil) print
}') > "$OUTPUT"

echo "Filtrage terminé : $(grep -v '^#' "$OUTPUT" | wc -l) variants conservés (DV >= $SEUIL)"