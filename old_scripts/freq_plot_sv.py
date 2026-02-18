#!/usr/bin/env python3
import sys
import math
import matplotlib.pyplot as plt

vcf_file = sys.argv[1]

vaf_values = []

with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 10:
            continue

        fields = parts[7].split(";")
        for field in fields:
            if field.startswith("VAF="):
                vaf = field.split("=")[1]
                vaf_values.append(float(vaf))
        
if not vaf_values:
    print("Aucune valeur VAF trouvée dans le fichier VCF.")
    sys.exit(1)

min_vaf  = min(vaf_values)
max_vaf  = max(vaf_values)
mean_vaf = sum(vaf_values) / len(vaf_values)

print(f"{len(vaf_values)} variants trouvés")
print(f"Profondeur min : {min_vaf}, max : {max_vaf}, moyenne : {mean_vaf:.1f}")

# ── Fenêtres ──────────────────────────────────────────────────
bins = [i * 0.05 for i in range(21)]  # 0.0, 0.05, 0.10, ..., 1.0
# ── Plot ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))

ax.hist(vaf_values, bins=bins, color="#4C72B0", edgecolor="white")
ax.axvline(mean_vaf, color="red", linestyle="--", linewidth=1.2, label=f"Moyenne : {mean_vaf:.1f}")

ax.set_xlabel("Fréquence allélique (vaf)")
ax.set_ylabel("Nombre de variants")
ax.set_title("Distribution de la fréquence allélique (vaf)")
ax.legend()

fig.suptitle(vcf_file, fontsize=9, color="grey")
fig.tight_layout()
fig.savefig(sys.argv[1] + "_vaf_distribution.png", dpi=150)
plt.show()