#!/usr/bin/env python3
import sys
import os
from collections import Counter

def parse_snpeff_vcf(vcf_file):
    snps = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): 
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, ref, alt, info = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
            
            # Extraire les annotations SnpEff (ANN=...)
            annotations = []
            for part in info.split(';'):
                if part.startswith('ANN='):
                    ann_data = part[4:]  # Enlever 'ANN='
                    for ann in ann_data.split(','):
                        ann_parts = ann.split('|')
                        if len(ann_parts) >= 4:
                            annotations.append({
                                'effect': ann_parts[1],
                                'impact': ann_parts[2],
                                'gene': ann_parts[3],
                                'transcript': ann_parts[6] if len(ann_parts) > 6 else '',
                                'aa_change': ann_parts[10] if len(ann_parts) > 10 else ''
                            })
            
            snps.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'annotations': annotations
            })
    
    return snps

def analyze_effects(snps):
    effects = Counter()
    impacts = Counter()
    genes = Counter()
    
    for snp in snps:
        for ann in snp['annotations']:
            effects[ann['effect']] += 1
            impacts[ann['impact']] += 1
            if ann['gene']:
                genes[ann['gene']] += 1
    
    return effects, impacts, genes

def main():
    # Vérification des arguments
    if len(sys.argv) < 2:
        print("Usage: python3 snp_analyzer.py <fichier.vcf> [dossier_output]")
        print("Exemples:")
        print("  python3 snp_analyzer.py snp_25-4_annotated.vcf")
        print("  python3 snp_analyzer.py snp_25-4_annotated.vcf /path/to/results/")
        print("  python3 snp_analyzer.py data/sample.vcf ./output/")
        sys.exit(1)
    
    # Arguments d'entrée
    vcf_file = sys.argv[1]
    
    # Dossier de sortie
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        output_dir = "."  # Dossier courant par défaut
    
    # Créer le nom du fichier de sortie
    base_name = os.path.splitext(os.path.basename(vcf_file))[0]
    output_filename = f"summary_{base_name}.txt"
    output_file = os.path.join(output_dir, output_filename)
    
    # Vérifier que le fichier VCF existe
    if not os.path.exists(vcf_file):
        print(f"Erreur: Le fichier '{vcf_file}' n'existe pas.")
        sys.exit(1)
    
    # Créer le dossier de sortie s'il n'existe pas
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Dossier créé: {output_dir}")
        except Exception as e:
            print(f"Erreur lors de la création du dossier '{output_dir}': {e}")
            sys.exit(1)
    
    print(f"Analyse du fichier: {vcf_file}")
    print(f"Sortie: {output_file}")
    print()
    
    # Analyse
    try:
        snps = parse_snpeff_vcf(vcf_file)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier VCF: {e}")
        sys.exit(1)
    
    print("=" * 60)
    print("ANALYSE DES ANNOTATIONS SnpEff")
    print("=" * 60)
    print(f"Total SNPs: {len(snps)}")
    print()
    
    effects, impacts, genes = analyze_effects(snps)
    
    print("EFFETS DES VARIANTS:")
    print("-" * 30)
    for effect, count in effects.most_common(10):
        print(f"{effect:30} {count:6} ({count/len(snps)*100:.1f}%)")
    
    print(f"\nIMPACT DES VARIANTS:")
    print("-" * 30)
    for impact, count in impacts.most_common():
        print(f"{impact:15} {count:6} ({count/len(snps)*100:.1f}%)")
    
    print(f"\nGÈNES LES PLUS AFFECTÉS:")
    print("-" * 30)
    for gene, count in genes.most_common(10):
        if gene:  # Ignorer les gènes vides
            print(f"{gene:25} {count:6} variants")
    
    # Sauvegarder le résumé détaillé
    try:
        with open(output_file, 'w') as out:
            out.write("ANALYSE DÉTAILLÉE DES SNPs KHV-U\n")
            out.write("=" * 50 + "\n\n")
            
            out.write("RÉSUMÉ:\n")
            out.write(f"Total SNPs analysés: {len(snps)}\n")
            out.write(f"Total annotations: {sum(effects.values())}\n\n")
            
            out.write("DISTRIBUTION DES EFFETS:\n")
            for effect, count in effects.most_common():
                out.write(f"{effect}: {count} ({count/len(snps)*100:.1f}%)\n")
            
            out.write(f"\nDISTRIBUTION DES IMPACTS:\n")
            for impact, count in impacts.most_common():
                out.write(f"{impact}: {count} ({count/len(snps)*100:.1f}%)\n")
            
            out.write(f"\nGÈNES LES PLUS AFFECTÉS:\n")
            for gene, count in genes.most_common(15):
                if gene:
                    out.write(f"{gene}: {count} variants\n")
            
            out.write(f"\nDÉTAIL PAR SNP (premiers 50):\n")
            out.write("Pos\tRef>Alt\tEffet\tImpact\tGène\tChangement_AA\n")
            for snp in snps[:50]:
                for ann in snp['annotations']:
                    out.write(f"{snp['pos']}\t{snp['ref']}>{snp['alt']}\t{ann['effect']}\t{ann['impact']}\t{ann['gene']}\t{ann['aa_change']}\n")
        
        print(f"\nRésultats détaillés sauvegardés dans: {output_file}")
        
    except Exception as e:
        print(f"Erreur lors de l'écriture du fichier de sortie: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()