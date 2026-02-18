import vcfpy
import argparse
import sys
import os
import matplotlib.pyplot as plt

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

    args = parser.parse_args()

    if len(args.input) > 1 and args.out and not args.merge:
        print("Attention : Plusieurs fichiers fournis sans l'option -merge.")
        print("Veuillez activer -merge pour fusionner dans l'output, ou traiter les fichiers un par un.")
        sys.exit(1)
    
    # Identification automatique basée sur ton format de fichier
    filename = os.path.basename(args.input[0]).lower()
    is_sv = "sv" in filename
    is_snp = "snp" in filename or not is_sv

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
            origin_name = os.path.basename(vcf_file).split(".vcf")[0]
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

                # Si VAF est à 0 dans INFO, on tente le calcul via DV/DP (pour Sniffles)
                if vaf == 0:
                    vaf = dv / dp if dp > 0 else 0.0

                # Remplissage des données selon l'option choisie
                if args.distrib_qual: data_to_plot.append(gq)
                if args.distrib_dv: data_to_plot.append(dv)
                if args.distrib_vaf: data_to_plot.append(vaf)

                if keep:
                    total_out += 1
                    # On renomme l'échantillon dans le record pour correspondre au writer
                    old_sample_name = record.calls[0].sample
                    if old_sample_name != 'SAMPLE':
                        record.call_for_sample['SAMPLE'] = record.call_for_sample.pop(old_sample_name)
                        record.calls[0].sample = 'SAMPLE'

                    if writer:
                        if args.merge:
                            record.INFO['ORIGIN'] = origin_name
                        writer.write_record(record)

    # Après la boucle "for vcf_file in args.input:"
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