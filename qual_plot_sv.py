#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt

vcf_file = sys.argv[1]

# Définition des intervalles
vaf_bins      = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
ratio_bins    = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

dico_vaf   = {}
dico_ratio = {}

# Initialisation des tranches
for i in range(len(vaf_bins) - 1):
    key = f"{vaf_bins[i]}-{vaf_bins[i+1]}"
    dico_vaf[key]   = 0
    dico_ratio[key] = 0

def get_bin(value, bins):
    for i in range(len(bins) - 1):
        if bins[i] <= value < bins[i+1]:
            return f"{bins[i]}-{bins[i+1]}"
    return f"{bins[-2]}-{bins[-1]}"  # dernière tranche (valeur == 1.0)

with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        infos = parts[7].split(";")  # tous les champs INFO

        support  = None
        coverage = None
        vaf      = None

        for info in infos:
            if info.startswith("SUPPORT="):
                support = int(info.split("=")[1])
            elif info.startswith("COVERAGE="):
                coverage = int(info.split("=")[1].split(",")[0])
            elif info.startswith("VAF="):
                vaf = float(info.split("=")[1])

        if vaf is not None:
            dico_vaf[get_bin(vaf, vaf_bins)] += 1

        if support is not None and coverage is not None and coverage > 0:
            ratio = support / coverage
            dico_ratio[get_bin(ratio, ratio_bins)] += 1

# ── Plot ──────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

labels = list(dico_vaf.keys())

ax1.bar(range(len(labels)), dico_vaf.values(), color="#4C72B0", edgecolor="white")
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, rotation=45, ha="right")
ax1.set_xlabel("VAF")
ax1.set_ylabel("Nombre de variants")
ax1.set_title("Distribution des VAF")

ax2.bar(range(len(labels)), dico_ratio.values(), color="#55A868", edgecolor="white")
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels, rotation=45, ha="right")
ax2.set_xlabel("SUPPORT / COVERAGE")
ax2.set_ylabel("Nombre de variants")
ax2.set_title("Distribution du ratio Support/Coverage")

plt.suptitle(vcf_file, fontsize=10, color="grey")
plt.tight_layout()
plt.savefig("vaf_ratio_distribution.png", dpi=150)
plt.show()