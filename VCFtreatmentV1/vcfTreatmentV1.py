import vcfpy
import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from collections import Counter

def run_snpeff(input_vcf, genome="KHV_U"):
    config_file = "local_snpeff.config"
    gff_file = "reference/genes.gff"
    fasta_file = "reference/sequences.fa"
    
    # 1. Vérifier si la base de données SnpEff existe déjà dans ./data/KHV_U
    db_path = f"./data/{genome}/snpEffectPredictor.bin"
    
    if not os.path.exists(db_path):
        print(f"--- Base de données {genome} introuvable. Construction en cours... ---")
        # Création de l'arborescence requise par SnpEff
        os.makedirs(f"data/{genome}", exist_ok=True)
        # SnpEff attend des noms précis : genes.gff et sequences.fa
        import shutil
        shutil.copy(gff_file, f"data/{genome}/genes.gff")
        shutil.copy(fasta_file, f"data/{genome}/sequences.fa")
        
        # Commande pour construire la base
        build_cmd = ["snpEff", "build", "-gff3", "-v", "-noCheckCds", "-noCheckProtein", "-c", config_file, genome]
        subprocess.run(build_cmd, check=True)

    # 2. Lancer l'annotation
    output_ann = input_vcf.replace(".vcf", ".ann.vcf")
    print(f"--- Annotation de {input_vcf} ---")
    ann_cmd = ["snpEff", "ann", "-noStats", "-c", config_file, genome, input_vcf]
    
    with open(output_ann, "w") as out_f:
        subprocess.run(ann_cmd, stdout=out_f, check=True)
    
    return output_ann

def get_detailed_annotations(all_records):
    """Extrait les détails des annotations pour le résumé textuel."""
    detailed_ann = []
    for r in all_records:
        if 'ANN' in r.INFO:
            for ann in r.INFO['ANN']:
                parts = ann.split('|')
                if len(parts) >= 11:
                    detailed_ann.append({
                        'pos': r.POS,
                        'ref': r.REF,
                        'alt': r.ALT[0].value if hasattr(r.ALT[0], 'value') else r.ALT[0],
                        'effect': parts[1],
                        'impact': parts[2],
                        'gene': parts[3],
                        'aa_change': parts[10]
                    })
    return detailed_ann

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
    plt.close()

def main():

    parser = argparse.ArgumentParser(description="Analyseur VCF évolutif pour populations virales")
    
    # Entrées/Sorties
    parser.add_argument("-i", "--input", required=True, nargs='+', help="Fichier(s) VCF d'entrée (un ou plusieurs)")    
    parser.add_argument("-o", "--out", help="Fichier de sortie filtré")
    
    # Options SnpEff
    parser.add_argument("--annotate", action="store_true", help="Lance l'annotation SnpEff (nécessite snpEff installé via conda)")
    parser.add_argument("--genome", default="KHV_U", help="Génome de référence pour snpEff (default: KHV_U)")

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

    # --- 1. DÉFINITION DU PRÉFIXE (Pour éviter UnboundLocalError) ---
    if args.out:
        prefix = args.out.replace(".vcf", "")
    else:
        prefix = os.path.basename(args.input[0]).split(".")[0]

    generated_files = [] # Liste pour stocker les noms des fichiers créés

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
        # Annotation SnpEff si demandé (on annotera avant de filtrer pour garder les impacts même sur les variants filtrés)
        if args.annotate:
            vcf_file = run_snpeff(vcf_file, genome=args.genome)

        with vcfpy.Reader.from_path(vcf_file) as reader:
            origin_name = os.path.basename(vcf_file).split(".")[0]

            # Ajout manuel de la définition ANN si manquante dans le header pour vcfpy
            if args.annotate and 'ANN' not in reader.header.info_ids():
                reader.header.add_info_line({'ID': 'ANN', 'Number': '.', 'Type': 'String', 'Description': 'Functional annotations'})

            for record in reader:
                total_in += 1

                # --- Logique d'Impact SnpEff ---
                record.impact = "UNKNOWN"
                if 'ANN' in record.INFO:
                    first_ann = record.INFO['ANN'][0].split('|')
                    if len(first_ann) > 2:
                        record.impact = first_ann[2] # HIGH, MODERATE, etc.

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
                
                # PARSING DES ANNOTATIONS SnpEff (Champ ANN)
                record.impact = "UNKNOWN"
                if 'ANN' in record.INFO:
                    # On prend le premier impact de la liste d'annotations
                    first_ann = record.INFO['ANN'][0].split('|')
                    if len(first_ann) > 2:
                        record.impact = first_ann[2] # HIGH, MODERATE, LOW, MODIFIER

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

    # Finalisation des opérations
    if writer:
        writer.close()
        generated_files.append(args.out)
        if args.filter_snp or args.filter_sv:
            print(f"Filtrage : {total_out}/{total_in} variants conservés dans {args.out}")

    # Graphiques et summary

    if args.annotate and all_records:
        # Utilisation de Counter (nécessite 'from collections import Counter' en haut du script)
        impact_counts = Counter([r.impact for r in all_records])
        if impact_counts:
            labels = impact_counts.keys()
            sizes = impact_counts.values()
            colors_map = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c', 'MODIFIER': '#1f77b4', 'UNKNOWN': '#7f7f7f'}
            
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, 
                   colors=[colors_map.get(l, '#7f7f7f') for l in labels])
            ax.set_title(f"Répartition des impacts fonctionnels\n({len(all_records)} variants)")
            
            pie_output = f"{prefix}_impacts_pie.png"
            plt.savefig(pie_output, dpi=150)
            generated_files.append(pie_output)

    # GÉNÉRATION DU RÉSUMÉ TEXTUEL
    if args.annotate and all_records:
        details = get_detailed_annotations(all_records)
        summary_file = f"{prefix}_summary.txt"
        
        # Calcul des compteurs
        eff_counts = Counter([d['effect'] for d in details])
        imp_counts = Counter([d['impact'] for d in details])
        gene_counts = Counter([d['gene'] for d in details if d['gene']])

        try:
            with open(summary_file, 'w') as out:
                out.write(f"ANALYSE DÉTAILLÉE DES VARIANTS : {prefix}\n\n")
                out.write(f"Total variants : {len(all_records)}\n")
                out.write(f"Total annotations : {len(details)}\n\n")

                out.write("DISTRIBUTION DES IMPACTS :\n")
                for imp, count in imp_counts.most_common():
                    out.write(f"{imp:15} : {count} ({count/len(details)*100:.1f}%)\n")

                out.write("\nGÈNES LES PLUS AFFECTÉS :\n")
                for gene, count in gene_counts.most_common(15):
                    out.write(f"{gene:25} : {count} variants\n")

                out.write("\nDÉTAILS DES VARIANTS (Impacts HIGH et MODERATE) :\n")
                out.write("Pos\tRef>Alt\tEffet\tImpact\tGène\tAA_Change\n")
                # On filtre pour ne pas surcharger le fichier avec les MODIFIER/LOW
                for d in details:
                    if d['impact'] in ['HIGH', 'MODERATE']:
                        out.write(f"{d['pos']}\t{d['ref']}>{d['alt']}\t{d['effect']}\t{d['impact']}\t{d['gene']}\t{d['aa_change']}\n")

            generated_files.append(summary_file)
        except Exception as e:
            print(f"Erreur lors de l'écriture du résumé : {e}")

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
        generated_files.append(output_name)

    if data_to_plot:
        # Si un output est défini, on l'utilise comme préfixe, sinon on prend le premier input
        prefix = args.out.replace(".vcf", "") if args.out else args.input[0].replace(".vcf", "")
        
        if args.distrib_qual: 
            plot_data(data_to_plot, "Qualité (GQ)", "GQ", f"{prefix}_qual.png")
        if args.distrib_dv: 
            plot_data(data_to_plot, "Profondeur (DV)", "DV", f"{prefix}_dv.png")
        if args.distrib_vaf: 
            plot_data(data_to_plot, "VAF", "VAF", f"{prefix}_vaf.png", is_vaf=True)

    if args.count:
        target = args.out if args.out else "fichiers d'entrée"
        count_val = total_out if args.out else total_in
        print(f"{count_val} variants identifiés dans {target}")

    # --- AFFICHAGE FINAL ---
    if generated_files:
        print("\nTraitement terminé avec succès")
        print("Fichiers générés :")
        for f in generated_files:
            print(f" - {f}")

if __name__ == "__main__":
    main()