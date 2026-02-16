#!/bin/bash
set -e

VCF_FILE="$1"
OUTPUT_DIR="${2:-./results}"

if [ ! -f "$VCF_FILE" ]; then
    echo "Erreur: Fichier VCF '$VCF_FILE' introuvable"
    exit 1
fi

if [ ! -f "$HOME/snpEff/snpEff.jar" ]; then
    echo "Erreur: snpEff.jar introuvable dans $HOME/snpEff/"
    exit 1
fi

if [ ! -f "summary_snpeff.py" ]; then
    echo "Erreur: summary_snpeff.py introuvable"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

BASENAME=$(basename "$VCF_FILE" .vcf)
ANNOTATED_VCF="$OUTPUT_DIR/${BASENAME}_annotated.vcf"
STATS_HTML="$OUTPUT_DIR/${BASENAME}_stats.html"

echo "Annotation SnpEff..."
java -jar ~/snpEff/snpEff.jar ann -stats "$STATS_HTML" KHV_U "$VCF_FILE" > "$ANNOTATED_VCF"

echo "Génération summary..."
python3 summary_snpeff.py "$ANNOTATED_VCF" "$OUTPUT_DIR"

echo "Terminé. Fichiers dans $OUTPUT_DIR"