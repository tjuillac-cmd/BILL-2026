#!/usr/bin/env python3
import sys
import math
import matplotlib.pyplot as plt

vcf_file = sys.argv[1]

depth_values = []

with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 10:
            continue

        depth = parts[9].split(":")[3]
        depth_values.append(int(depth))
        
if not depth_values:
    print("Aucune valeur DP trouvée dans le fichier VCF.")
    sys.exit(1)

min_dp  = min(depth_values)
max_dp  = max(depth_values)
mean_dp = sum(depth_values) / len(depth_values)

print(f"{len(depth_values)} variants trouvés")
print(f"Profondeur min : {min_dp}, max : {max_dp}, moyenne : {mean_dp:.1f}")

# ── Fenêtres ──────────────────────────────────────────────────
n_bins = 10
step   = math.ceil((max_dp - min_dp) / n_bins)
bins   = list(range(min_dp, max_dp + step + 1, step))

# ── Plot ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))

ax.hist(depth_values, bins=bins, color="#4C72B0", edgecolor="white")
ax.axvline(mean_dp, color="red", linestyle="--", linewidth=1.2, label=f"Moyenne : {mean_dp:.1f}")

ax.set_xlabel("Profondeur (DP)")
ax.set_ylabel("Nombre de variants")
ax.set_title("Distribution de la profondeur (DP)")
ax.legend()

fig.suptitle(vcf_file, fontsize=9, color="grey")
fig.tight_layout()
fig.savefig("depth_distribution.png", dpi=150)
plt.show()