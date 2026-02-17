#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

vcf_file = sys.argv[1]

dico = {}
with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        qual = parts[9].split(":")[1]  # colonne QUAL
        if qual not in dico:
            dico[qual] = 0
        dico[qual] += 1

# Tri par valeur de QUAL croissante
quals  = sorted(dico.keys(), key=float)
counts = [dico[q] for q in quals]

plt.bar(range(len(quals)), counts)
plt.xticks(range(len(quals)), [round(float(q), 1) for q in quals], rotation=90)
plt.xlabel("QUAL")
plt.ylabel("Nombre de variants")
plt.title("Distribution des scores QUAL")
plt.tight_layout()
plt.savefig("qual_distribution.png", dpi=150)
plt.show()

