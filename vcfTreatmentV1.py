import vcfpy
import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

def plot_data(values, title, xlabel, output_path, is_vaf=False):
    """Génère un graphique basé sur tes styles (hist ou bar)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    
    if is_vaf:
        # Style de ton script freq_plot_sv.py
        bins = [i * 0.05 for i in range(21)]
        ax.hist(values, bins=bins, color="#4C72B0", edgecolor="white")
    else:
        # Style de ton script qual_plot_sv.py (auto-binning)
        ax.hist(values, bins=15, color="#4C72B0", edgecolor="white")

    mean_val = sum(values) / len(values) if values else 0
    ax.axvline(mean_val, color="red", linestyle="--", label=f"Moyenne : {mean_val:.2f}")
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Nombre de variants")
    ax.legend()
    fig.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Graphique généré : {output_path}")
    plt.close()

def main():

    parser = argparse.ArgumentParser(description="Analyseur VCF évolutif pour populations virales")
    
    # Entrées/Sorties
    parser.add_argument("-i", "--input", required=True, nargs='+', help="Fichier(s) VCF d'entrée (un ou plusieurs)")    
    parser.add_argument("-o", "--out", help="Fichier de sortie filtré")
    
    # Options de filtrage
    parser.add_argument("--filter-snp", type=float, help="Seuil de qualité GQ pour SNP")
    parser.add_argument("--filter-sv", type=int, help="Seuil de profondeur DV pour SV")
    
    # Options de plotting
    parser.add_argument("--distrib-qual", action="store_true", help="Distribution GQ (SNP uniquement)")
    parser.add_argument("--distrib-dv", action="store_true", help="Distribution DV (SV uniquement)")
    parser.add_argument("--distrib-vaf", action="store_true", help="Distribution VAF (SV uniquement)")

    # Option de counting
    parser.add_argument("-count", action="store_true", help="Affiche le nombre de variants")
    
    # Option de merging
    parser.add_argument("-merge", action="store_true", help="Fusionne les fichiers en un seul avec trace d'origine")

    # Option de comparaison entre fichiers
    parser.add_argument("--compare", nargs=2, metavar=('ORIGIN1', 'ORIGIN2'), help="Compare deux origines (ex: P25-4 P25-8)")
    parser.add_argument("--add-vaf", action="store_true", help="Affiche la VAF moyenne et l'écart-type sur le graphique de comparaison")

    args = parser.parse_args()

    if len(args.input) > 1 and args.out and not args.merge:
        print("Attention : Plusieurs fichiers fournis sans l'option -merge.")
        print("Veuillez activer -merge pour fusionner dans l'output, ou traiter les fichiers un par un.")
        sys.exit(1)
    
    # Ouverture d'un reader temporaire pour inspecter le contenu
    with vcfpy.Reader.from_path(args.input[0]) as inspect_reader:
        # On vérifie si des lignes ##ALT=<ID=INS...> existent dans le header
        alt_lines = inspect_reader.header.get_lines('ALT')
        has_sv_alt = any(line.id in ['INS', 'DEL', 'INV', 'DUP', 'BND'] for line in alt_lines)
        
        # On vérifie si le champ VAF ou DV existe dans les définitions INFO/FORMAT
        has_sv_fields = 'VAF' in inspect_reader.header.info_ids() or \
                        'DV' in inspect_reader.header.format_ids()
        
        is_sv = has_sv_alt or has_sv_fields
        is_snp = not is_sv

    # Sécurités demandées
    if args.distrib_qual and is_sv:
        print("Erreur : --distrib-qual nécessite un fichier SNP.")
        sys.exit(1)
    if (args.distrib_dv or args.distrib_vaf) and not is_sv:
        print("Erreur : Les options DV/VAF nécessitent un fichier SV.")
        sys.exit(1)

    # Init des compteurs et de la liste pour le plotting
    total_in, total_out = 0, 0
    data_to_plot = []
    all_records = []  # Pour stocker tous les records si comparaison demandée
    writer = None

    if args.out:
        try:
            with vcfpy.Reader.from_path(args.input[0]) as header_reader:
                header = header_reader.header
                # On force le nom de l'échantillon à "SAMPLE" pour le fichier de sortie
                header.samples.names = ['SAMPLE']
                if args.merge:
                    header.add_info_line({'ID': 'ORIGIN', 'Number': '1', 'Type': 'String', 'Description': 'Fichier source du variant'})
            writer = vcfpy.Writer.from_path(args.out, header)
        except Exception as e:
            print(f"Erreur lors de l'initialisation du fichier de sortie : {e}")
            sys.exit(1)

    for vcf_file in args.input:
        with vcfpy.Reader.from_path(vcf_file) as reader:
            origin_name = os.path.basename(vcf_file).split(".")[0]
            for record in reader:
                total_in += 1

                # Récupération des valeurs pour le plotting et le filtrage
                sample = record.calls[0].data
                gq = float(sample.get('GQ', 0))
                dv = int(sample.get('DV', 0)) if is_sv else 0
                dp = int(sample.get('DP', 1)) if is_sv else 1

                # --- Logique de filtrage ---
                keep = True
                if is_snp and args.filter_snp and gq < args.filter_snp: keep = False
                if is_sv and args.filter_sv and dv < args.filter_sv: keep = False
            
                # Récupération sécurisée de la VAF
                vaf_info = record.INFO.get('VAF', 0.0)
                if isinstance(vaf_info, list):
                    vaf = float(vaf_info[0])
                else:
                    vaf = float(vaf_info)

                record.temp_vaf = vaf # On stocke la valeur temporairement dans l'objet
                
                # Si VAF est à 0 dans INFO, on tente le calcul via DV/DP (pour Sniffles)
                if vaf == 0:
                    vaf = dv / dp if dp > 0 else 0.0

                # Remplissage des données selon l'option choisie
                if args.distrib_qual: data_to_plot.append(gq)
                if args.distrib_dv: data_to_plot.append(dv)
                if args.distrib_vaf: data_to_plot.append(vaf)

                if keep:
                    total_out += 1
                    all_records.append(record)
                    # On renomme l'échantillon dans le record pour correspondre au writer
                    old_sample_name = record.calls[0].sample
                    if old_sample_name != 'SAMPLE':
                        record.call_for_sample['SAMPLE'] = record.call_for_sample.pop(old_sample_name)
                        record.calls[0].sample = 'SAMPLE'

                    if writer:
                        if args.merge:
                            record.INFO['ORIGIN'] = origin_name
                        writer.write_record(record)

    if args.compare:
        o1, o2 = args.compare
        available_origins = sorted(list(set(r.INFO.get('ORIGIN') for r in all_records if 'ORIGIN' in r.INFO)))

        if o1 not in available_origins or o2 not in available_origins:
            print(f"Erreur : Origines invalides. Disponibles : {', '.join(available_origins)}")
            sys.exit(1)

        # Calcul des ensembles de signatures
        vafs_o1 = {f"{r.CHROM}_{r.POS}_{r.REF}_{r.ALT}": getattr(r, 'temp_vaf', 0) for r in all_records if r.INFO.get('ORIGIN') == o1}
        vafs_o2 = {f"{r.CHROM}_{r.POS}_{r.REF}_{r.ALT}": getattr(r, 'temp_vaf', 0) for r in all_records if r.INFO.get('ORIGIN') == o2}

        sigs1, sigs2 = set(vafs_o1.keys()), set(vafs_o2.keys())

        data_groups = {
            f"Unique {o1}": [vafs_o1[s] for s in (sigs1 - sigs2)],
            "Communs": [(vafs_o1[s] + vafs_o2[s])/2 for s in (sigs1 & sigs2)],
            f"Unique {o2}": [vafs_o2[s] for s in (sigs2 - sigs1)]
        }

        labels = list(data_groups.keys())
        counts = [len(v) for v in data_groups.values()]

        # Initialisation du plot
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # 1. Barres pour le nombre de variants (toujours présentes)
        color_bars = ['#ff9999','#66b3ff','#99ff99']
        ax1.bar(labels, counts, color=color_bars, alpha=0.6, edgecolor='black')
        ax1.set_ylabel("Nombre de variants")
        
        # Ajout des étiquettes de texte sur les barres
        for i, v in enumerate(counts):
            ax1.text(i, v + (max(counts)*0.01), str(v), ha='center', fontweight='bold')

        # 2. Ajout conditionnel de la VAF (Axe secondaire)
        if args.add_vaf:
            if not is_sv:
                print("Note : L'option --add-vaf est ignorée car ce ne sont pas des fichiers SV.")
            else:
                means = [np.mean(v) if v else 0 for v in data_groups.values()]
                stds = [np.std(v) if v else 0 for v in data_groups.values()]
                
                ax2 = ax1.twinx()
                ax2.errorbar(labels, means, yerr=stds, fmt='o', color='red', capsize=8, elinewidth=2, label="VAF Moyenne ± SD")
                ax2.set_ylabel("VAF (Fréquence Allélique)", color='red')
                ax2.set_ylim(0, 1.1)
                ax2.legend(loc='upper right')

        plt.title(f"Comparaison de populations : {o1} vs {o2}")
        output_name = f"compare_{o1}_{o2}.png"
        plt.savefig(output_name, dpi=150)
        print(f"Graphique généré : {output_name}")

    if data_to_plot:
        # Si un output est défini, on l'utilise comme préfixe, sinon on prend le premier input
        prefix = args.out.replace(".vcf", "") if args.out else args.input[0].replace(".vcf", "")
        
        if args.distrib_qual: 
            plot_data(data_to_plot, "Qualité (GQ)", "GQ", f"{prefix}_qual.png")
        if args.distrib_dv: 
            plot_data(data_to_plot, "Profondeur (DV)", "DV", f"{prefix}_dv.png")
        if args.distrib_vaf: 
            plot_data(data_to_plot, "VAF", "VAF", f"{prefix}_vaf.png", is_vaf=True)

    # Finalisation
    if writer:
        writer.close()
        print(f"Filtrage : {total_out}/{total_in} variants conservés dans {args.out}")

    if args.count:
        target = args.out if args.out else "fichiers d'entrée"
        count_val = total_out if args.out else total_in
        print(f"{count_val} variants identifiés dans {target}")


if __name__ == "__main__":
    main()