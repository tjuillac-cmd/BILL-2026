#!/usr/bin/env python3
"""
Variant Annotator - Final Enhanced Version
Auto-detects files from folder, includes quality info, and filters by quality
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import subprocess
import re


class FileDetector:
    """Detect VCF, GFF, and BAM files in a directory"""
    
    def __init__(self, input_dir: Path):
        self.input_dir = input_dir
        self.vcf_files = []
        self.bam_files = []
        self.gff_files = []
        
    def scan_directory(self):
        """Scan directory and find relevant files"""
        if not self.input_dir.exists():
            print(f"Error: Directory not found: {self.input_dir}", file=sys.stderr)
            sys.exit(1)
        
        for file_path in self.input_dir.iterdir():
            if not file_path.is_file():
                continue
            
            name_lower = file_path.name.lower()
            
            if name_lower.endswith('.vcf') or name_lower.endswith('.vcf.gz'):
                self.vcf_files.append(file_path)
            elif name_lower.endswith('.bam'):
                self.bam_files.append(file_path)
            elif name_lower.endswith(('.gff', '.gff3', '.gtf')):
                self.gff_files.append(file_path)
        
        # Sort for consistency
        self.vcf_files.sort()
        self.bam_files.sort()
        self.gff_files.sort()
    
    def print_summary(self):
        """Print detected files summary"""
        print("\n" + "=" * 60)
        print("FILE DETECTION")
        print("=" * 60)
        print(f"Directory: {self.input_dir.absolute()}\n")
        
        print(f"VCF files: {len(self.vcf_files)}")
        for f in self.vcf_files:
            print(f"  • {f.name}")
        
        print(f"\nGFF files: {len(self.gff_files)}")
        for f in self.gff_files:
            print(f"  • {f.name}")
        
        print(f"\nBAM files: {len(self.bam_files)}")
        for f in self.bam_files:
            print(f"  • {f.name}")
        
        print("=" * 60 + "\n")


class VariantAnnotator:
    """Annotate variants with gene location, quality, and coverage"""
    
    def __init__(self, vcf_file: Path, gff_file: Optional[Path] = None, 
                 bam_file: Optional[Path] = None, min_quality: float = 0):
        self.vcf_file = vcf_file
        self.gff_file = gff_file
        self.bam_file = bam_file
        self.min_quality = min_quality
        
        self.variants = []
        self.genes = []
        self.filtered_variants = []
        self.samtools_available = self._check_samtools()
        
        self.vcf_quality_column = None  # Will be set when parsing VCF
    
    def _check_samtools(self) -> bool:
        """Check if samtools is installed"""
        try:
            result = subprocess.run(['samtools', '--version'], 
                                   capture_output=True, timeout=2)
            return result.returncode == 0
        except:
            return False
    
    def parse_vcf_quality(self, fields: List[str]) -> Dict:
        """Parse quality information from VCF record"""
        quality_info = {
            'qual_score': None,
            'filter': None,
            'depth': None,
            'allele_freq': None,
            'quality_flag': 'PASS'
        }
        
        # QUAL column (index 5)
        if len(fields) > 5 and fields[5] != '.':
            try:
                quality_info['qual_score'] = float(fields[5])
            except:
                quality_info['qual_score'] = None
        
        # FILTER column (index 6)
        if len(fields) > 6:
            quality_info['filter'] = fields[6]
            if fields[6] not in ['.', 'PASS', 'pass']:
                quality_info['quality_flag'] = 'FAIL'
        
        # INFO column (index 7) - extract DP and AF if available
        if len(fields) > 7:
            info = fields[7]
            
            # Extract depth (DP)
            dp_match = re.search(r'DP=(\d+)', info)
            if dp_match:
                quality_info['depth'] = int(dp_match.group(1))
            
            # Extract allele frequency (AF)
            af_match = re.search(r'AF=([\d.]+)', info)
            if af_match:
                quality_info['allele_freq'] = float(af_match.group(1))
        
        return quality_info
    
    def load_variants(self):
        """Load variants from VCF file with quality information"""
        print(f"→ Loading variants from {self.vcf_file.name}...")
        
        total_count = 0
        passed_count = 0
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                total_count += 1
                
                # Parse quality information
                quality_info = self.parse_vcf_quality(fields)
                
                # Check if passes quality filter
                passes_filter = True
                if self.min_quality > 0 and quality_info['qual_score'] is not None:
                    if quality_info['qual_score'] < self.min_quality:
                        passes_filter = False
                
                variant = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'id': fields[2] if fields[2] != '.' else None,
                    'ref': fields[3],
                    'alt': fields[4],
                    'qual_score': quality_info['qual_score'],
                    'filter': quality_info['filter'],
                    'depth': quality_info['depth'],
                    'allele_freq': quality_info['allele_freq'],
                    'quality_flag': quality_info['quality_flag'],
                    'passes_quality_filter': passes_filter,
                    'genes': [],
                    'coverage': None
                }
                
                self.variants.append(variant)
                
                if passes_filter:
                    passed_count += 1
        
        print(f"  Loaded {total_count} variants")
        if self.min_quality > 0:
            print(f"  {passed_count} variants pass quality filter (QUAL >= {self.min_quality})")
            print(f"  {total_count - passed_count} variants filtered out")
    
    def load_genes(self):
        """Load gene annotations from GFF file"""
        if not self.gff_file or not self.gff_file.exists():
            print("  No GFF file - skipping gene annotation")
            return
        
        print(f"→ Loading gene annotations from {self.gff_file.name}...")
        
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Parse attributes
                attributes = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key.strip()] = value.strip()
                
                # Extract gene name
                gene_name = (attributes.get('Name') or 
                           attributes.get('gene_name') or 
                           attributes.get('gene') or
                           attributes.get('ID') or
                           attributes.get('locus_tag') or
                           'Unknown')
                
                gene = {
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'type': fields[2],
                    'strand': fields[6],
                    'name': gene_name,
                    'attributes': attributes
                }
                
                self.genes.append(gene)
        
        print(f"  Loaded {len(self.genes)} gene features")
    
    def get_coverage(self, chrom: str, pos: int) -> Optional[Dict]:
        """Get coverage at position using samtools"""
        if not self.bam_file or not self.bam_file.exists():
            return None
        
        if not self.samtools_available:
            return None
        
        try:
            cmd = ['samtools', 'depth', '-r', f'{chrom}:{pos}-{pos}', 
                   str(self.bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if result.returncode == 0 and result.stdout.strip():
                fields = result.stdout.strip().split('\t')
                if len(fields) >= 3:
                    return {
                        'depth': int(fields[2]),
                        'available': True
                    }
        except:
            pass
        
        return {'depth': 0, 'available': False}
    
    def annotate_variants(self):
        """Annotate each variant with gene location"""
        print(f"→ Annotating variants with gene locations...")
        
        annotated_count = 0
        intergenic_count = 0
        
        for variant in self.variants:
            # Find overlapping genes
            overlapping = []
            for gene in self.genes:
                if (gene['chrom'] == variant['chrom'] and
                    gene['start'] <= variant['pos'] <= gene['end']):
                    overlapping.append(gene)
            
            if overlapping:
                variant['genes'] = overlapping
                annotated_count += 1
            else:
                variant['genes'] = []
                intergenic_count += 1
        
        print(f"  {annotated_count} variants in genes")
        print(f"  {intergenic_count} intergenic variants")
    
    def add_coverage_info(self):
        """Add BAM coverage information to variants"""
        if not self.bam_file:
            print("  No BAM file - skipping coverage")
            return
        
        if not self.samtools_available:
            print("  Samtools not available - skipping coverage")
            print("  Install samtools: http://www.htslib.org/")
            return
        
        print(f"→ Extracting coverage from {self.bam_file.name}...")
        
        for i, variant in enumerate(self.variants):
            coverage = self.get_coverage(variant['chrom'], variant['pos'])
            variant['coverage'] = coverage
            
            if (i + 1) % 100 == 0:
                print(f"  Processed {i + 1}/{len(self.variants)} variants...")
        
        print(f"  Coverage extracted for {len(self.variants)} variants")
    
    def get_location_type(self, variant: Dict) -> str:
        """Determine location type from genes"""
        if not variant['genes']:
            return 'intergenic'
        
        types = [g['type'] for g in variant['genes']]
        if 'CDS' in types:
            return 'coding'
        elif 'exon' in types:
            return 'exon'
        elif 'gene' in types:
            return 'gene'
        else:
            return 'genic'
    
    def write_annotation_table(self, output_file: Path, include_filtered: bool = False):
        """Write detailed annotation table"""
        print(f"→ Writing annotation table to {output_file.name}...")
        
        # Determine which variants to include
        if include_filtered:
            variants_to_write = self.variants
        else:
            variants_to_write = [v for v in self.variants if v['passes_quality_filter']]
        
        with open(output_file, 'w') as f:
            # Header
            f.write("# Variant Annotation Report\n")
            f.write(f"# VCF: {self.vcf_file.name}\n")
            if self.gff_file:
                f.write(f"# GFF: {self.gff_file.name}\n")
            if self.bam_file:
                f.write(f"# BAM: {self.bam_file.name}\n")
            if self.min_quality > 0:
                f.write(f"# Quality filter: QUAL >= {self.min_quality}\n")
            f.write("#\n")
            
            # Column headers
            headers = [
                "Chromosome",
                "Position",
                "Ref",
                "Alt",
                "QUAL_Score",
                "Filter",
                "VCF_Depth",
                "Allele_Freq",
                "BAM_Coverage",
                "Gene_Name",
                "Gene_Type",
                "Strand",
                "Location_Type",
                "Quality_Flag",
                "Variant_ID"
            ]
            f.write('\t'.join(headers) + '\n')
            
            # Data rows
            for variant in variants_to_write:
                # Format gene information
                if variant['genes']:
                    gene_names = [g['name'] for g in variant['genes']]
                    gene_types = [g['type'] for g in variant['genes']]
                    strands = [g['strand'] for g in variant['genes']]
                    
                    gene_name = '; '.join(gene_names)
                    gene_type = '; '.join(gene_types)
                    strand = strands[0]
                else:
                    gene_name = "Intergenic"
                    gene_type = "N/A"
                    strand = "."
                
                location = self.get_location_type(variant)
                
                # Quality info
                qual_score = f"{variant['qual_score']:.1f}" if variant['qual_score'] is not None else "N/A"
                filter_val = variant['filter'] if variant['filter'] else "."
                vcf_depth = str(variant['depth']) if variant['depth'] is not None else "N/A"
                allele_freq = f"{variant['allele_freq']:.3f}" if variant['allele_freq'] is not None else "N/A"
                
                # BAM coverage
                if variant['coverage'] and variant['coverage'].get('available'):
                    bam_coverage = str(variant['coverage']['depth'])
                else:
                    bam_coverage = "N/A"
                
                # Quality flag
                quality_flag = variant['quality_flag']
                
                # Variant ID
                variant_id = f"{variant['chrom']}:{variant['pos']}_{variant['ref']}>{variant['alt']}"
                
                # Write row
                row = [
                    variant['chrom'],
                    str(variant['pos']),
                    variant['ref'],
                    variant['alt'],
                    qual_score,
                    filter_val,
                    vcf_depth,
                    allele_freq,
                    bam_coverage,
                    gene_name,
                    gene_type,
                    strand,
                    location,
                    quality_flag,
                    variant_id
                ]
                f.write('\t'.join(row) + '\n')
        
        print(f"✓ Created {output_file}")
        print(f"  Wrote {len(variants_to_write)} variants")
    
    def write_annotated_bed(self, output_file: Path, include_filtered: bool = False):
        """Write color-coded BED file"""
        print(f"→ Writing color-coded BED file to {output_file.name}...")
        
        # Filter variants if needed
        if include_filtered:
            variants_to_write = self.variants
        else:
            variants_to_write = [v for v in self.variants if v['passes_quality_filter']]
        
        with open(output_file, 'w') as f:
            # BED track header
            f.write('track name="Annotated Variants" ')
            f.write('description="Variants color-coded by gene location" ')
            f.write('itemRgb="On"\n')
            
            for variant in variants_to_write:
                chrom = variant['chrom']
                start = variant['pos'] - 1  # BED is 0-based
                end = variant['pos']
                
                # Create name with gene and quality info
                if variant['genes']:
                    gene_names = [g['name'] for g in variant['genes']]
                    name = f"{variant['ref']}>{variant['alt']}|{','.join(gene_names)}"
                else:
                    name = f"{variant['ref']}>{variant['alt']}|intergenic"
                
                # Add quality to name
                if variant['qual_score'] is not None:
                    name += f"|Q{variant['qual_score']:.0f}"
                
                # Determine location for coloring
                location = self.get_location_type(variant)
                
                # Score based on location and quality
                location_scores = {
                    'coding': 1000,
                    'exon': 800,
                    'gene': 600,
                    'genic': 400,
                    'intergenic': 200
                }
                base_score = location_scores.get(location, 500)
                
                # Adjust score by quality if available
                if variant['qual_score'] is not None and variant['qual_score'] > 0:
                    # Scale quality (assume QUAL 0-100+ range)
                    quality_factor = min(variant['qual_score'] / 100, 1.0)
                    score = int(base_score * (0.5 + 0.5 * quality_factor))
                else:
                    score = base_score // 2  # Lower score if no quality
                
                # Color coding
                color_map = {
                    'coding': '255,0,0',        # Red
                    'exon': '255,128,0',        # Orange
                    'gene': '0,128,255',        # Blue
                    'genic': '0,255,128',       # Green
                    'intergenic': '128,128,128' # Gray
                }
                rgb = color_map.get(location, '0,0,0')
                
                # Dim color if low quality
                if variant['quality_flag'] == 'FAIL':
                    # Make colors more gray for failed quality
                    rgb = '128,128,128'
                
                strand = variant['genes'][0]['strand'] if variant['genes'] else '.'
                
                # BED12 format
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t")
                f.write(f"{start}\t{end}\t{rgb}\n")
        
        print(f"✓ Created {output_file}")
        print(f"  Wrote {len(variants_to_write)} variants")
    
    def write_summary(self, output_file: Path):
        """Write summary statistics"""
        print(f"→ Writing summary to {output_file.name}...")
        
        # Calculate statistics for all variants
        total_variants = len(self.variants)
        passed_variants = len([v for v in self.variants if v['passes_quality_filter']])
        
        location_counts = {'coding': 0, 'exon': 0, 'gene': 0, 'genic': 0, 'intergenic': 0}
        location_counts_passed = {'coding': 0, 'exon': 0, 'gene': 0, 'genic': 0, 'intergenic': 0}
        
        quality_stats = {
            'with_qual': 0,
            'without_qual': 0,
            'passed_filter': 0,
            'failed_filter': 0,
            'avg_qual': 0,
            'min_qual': float('inf'),
            'max_qual': 0
        }
        
        qual_sum = 0
        qual_count = 0
        
        genes_with_variants = set()
        
        for variant in self.variants:
            location = self.get_location_type(variant)
            location_counts[location] += 1
            
            if variant['passes_quality_filter']:
                location_counts_passed[location] += 1
            
            if variant['qual_score'] is not None:
                quality_stats['with_qual'] += 1
                qual_sum += variant['qual_score']
                qual_count += 1
                quality_stats['min_qual'] = min(quality_stats['min_qual'], variant['qual_score'])
                quality_stats['max_qual'] = max(quality_stats['max_qual'], variant['qual_score'])
            else:
                quality_stats['without_qual'] += 1
            
            if variant['quality_flag'] == 'PASS':
                quality_stats['passed_filter'] += 1
            else:
                quality_stats['failed_filter'] += 1
            
            for gene in variant['genes']:
                genes_with_variants.add(gene['name'])
        
        if qual_count > 0:
            quality_stats['avg_qual'] = qual_sum / qual_count
        
        with open(output_file, 'w') as f:
            f.write("# VARIANT ANNOTATION SUMMARY\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"VCF File: {self.vcf_file.name}\n")
            if self.gff_file:
                f.write(f"GFF File: {self.gff_file.name}\n")
            if self.bam_file:
                f.write(f"BAM File: {self.bam_file.name}\n")
            f.write("\n")
            
            f.write("VARIANT COUNTS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Variants: {total_variants}\n")
            if self.min_quality > 0:
                f.write(f"Quality Filter: QUAL >= {self.min_quality}\n")
                f.write(f"Passed Filter: {passed_variants} ({passed_variants/total_variants*100:.1f}%)\n")
                f.write(f"Filtered Out: {total_variants - passed_variants} ({(total_variants-passed_variants)/total_variants*100:.1f}%)\n")
            f.write(f"Total Gene Features: {len(self.genes)}\n")
            f.write(f"Genes with Variants: {len(genes_with_variants)}\n")
            f.write("\n")
            
            f.write("QUALITY STATISTICS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Variants with QUAL score: {quality_stats['with_qual']}\n")
            f.write(f"Variants without QUAL: {quality_stats['without_qual']}\n")
            if qual_count > 0:
                f.write(f"Average QUAL: {quality_stats['avg_qual']:.2f}\n")
                f.write(f"Min QUAL: {quality_stats['min_qual']:.2f}\n")
                f.write(f"Max QUAL: {quality_stats['max_qual']:.2f}\n")
            f.write(f"FILTER=PASS: {quality_stats['passed_filter']}\n")
            f.write(f"FILTER=FAIL: {quality_stats['failed_filter']}\n")
            f.write("\n")
            
            f.write("VARIANT LOCATIONS (All Variants):\n")
            f.write("-" * 40 + "\n")
            for location, count in sorted(location_counts.items(), key=lambda x: x[1], reverse=True):
                percentage = (count / total_variants * 100) if total_variants > 0 else 0
                f.write(f"{location:15s}: {count:5d} ({percentage:5.1f}%)\n")
            f.write("\n")
            
            if self.min_quality > 0 and passed_variants > 0:
                f.write(f"VARIANT LOCATIONS (Quality Filtered, QUAL >= {self.min_quality}):\n")
                f.write("-" * 40 + "\n")
                for location, count in sorted(location_counts_passed.items(), key=lambda x: x[1], reverse=True):
                    percentage = (count / passed_variants * 100) if passed_variants > 0 else 0
                    f.write(f"{location:15s}: {count:5d} ({percentage:5.1f}%)\n")
                f.write("\n")
            
            f.write("COLOR CODING (in BED file):\n")
            f.write("-" * 40 + "\n")
            f.write("🔴 Red    = Coding sequence (CDS)\n")
            f.write("🟠 Orange = Exon\n")
            f.write("🔵 Blue   = Gene\n")
            f.write("🟢 Green  = Other genic region\n")
            f.write("⚫ Gray   = Intergenic or low quality\n")
            f.write("\n")
            
            if genes_with_variants:
                f.write(f"TOP GENES BY VARIANT COUNT:\n")
                f.write("-" * 40 + "\n")
                
                # Count variants per gene
                gene_var_counts = {}
                for variant in self.variants:
                    if variant['passes_quality_filter']:  # Only count quality-filtered
                        for gene in variant['genes']:
                            name = gene['name']
                            gene_var_counts[name] = gene_var_counts.get(name, 0) + 1
                
                # Sort and show top 20
                for gene, count in sorted(gene_var_counts.items(), 
                                        key=lambda x: x[1], reverse=True)[:20]:
                    f.write(f"{gene:30s}: {count:3d} variants\n")
        
        print(f"✓ Created {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Variant Annotator - Final Version with Folder Input",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
FOLDER-BASED WORKFLOW:

Simply put all files in one folder and run:

  python variant_annotator_final.py -i input_folder -o output_folder

The tool will automatically detect:
  • VCF files (.vcf, .vcf.gz)
  • GFF files (.gff, .gff3, .gtf)
  • BAM files (.bam)

QUALITY FILTERING:

Filter variants by QUAL score:
  -q 20    Keep only variants with QUAL >= 20
  -q 30    More stringent (QUAL >= 30)

Output files:
  • annotated_variants.tsv       - Full table with quality info
  • annotated_variants.bed        - Color-coded for IGV
  • annotated_variants_all.tsv    - Includes filtered variants
  • annotation_summary.txt        - Statistics

Example folder structure:
  input_data/
  ├── sample.vcf       ← Your variants
  ├── genes.gff        ← Gene annotations
  └── reads.bam        ← Read alignments

Run:
  python variant_annotator_final.py -i input_data -o results -q 20

View in IGV:
  1. Load reads.bam
  2. Load sample.vcf
  3. Load results/annotated_variants.bed (color-coded)
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True, type=Path,
                       help='Input directory with VCF, GFF, BAM files')
    parser.add_argument('-o', '--output-dir', type=Path, default=Path('annotation_output'),
                       help='Output directory (default: annotation_output)')
    parser.add_argument('-q', '--min-quality', type=float, default=0,
                       help='Minimum QUAL score to include variant (default: 0 = no filter)')
    parser.add_argument('--include-filtered', action='store_true',
                       help='Also save filtered variants to separate file')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("VARIANT ANNOTATOR - FINAL VERSION")
    print("=" * 60)
    
    # Detect files
    detector = FileDetector(args.input_dir)
    detector.scan_directory()
    detector.print_summary()
    
    # Check requirements
    if not detector.vcf_files:
        print("Error: No VCF files found in input directory!")
        sys.exit(1)
    
    # Use first detected file of each type
    vcf_file = detector.vcf_files[0]
    gff_file = detector.gff_files[0] if detector.gff_files else None
    bam_file = detector.bam_files[0] if detector.bam_files else None
    
    if not gff_file:
        print("Warning: No GFF file found - variants will not be annotated with genes")
    
    # Initialize annotator
    print("=" * 60)
    print("PROCESSING")
    print("=" * 60)
    print()
    
    annotator = VariantAnnotator(
        vcf_file=vcf_file,
        gff_file=gff_file,
        bam_file=bam_file,
        min_quality=args.min_quality
    )
    
    # Load data
    annotator.load_variants()
    annotator.load_genes()
    
    # Annotate
    annotator.annotate_variants()
    
    # Add coverage if BAM provided
    if bam_file:
        annotator.add_coverage_info()
    
    # Write outputs
    print()
    print("=" * 60)
    print("GENERATING OUTPUT FILES")
    print("=" * 60)
    print()
    
    # Main output (quality-filtered)
    annotator.write_annotation_table(
        args.output_dir / 'annotated_variants.tsv',
        include_filtered=False
    )
    
    # Optional: all variants including filtered
    if args.include_filtered:
        annotator.write_annotation_table(
            args.output_dir / 'annotated_variants_all.tsv',
            include_filtered=True
        )
    
    annotator.write_annotated_bed(
        args.output_dir / 'annotated_variants.bed',
        include_filtered=False
    )
    
    annotator.write_summary(args.output_dir / 'annotation_summary.txt')
    
    # Final summary
    print()
    print("=" * 60)
    print("✓ ANNOTATION COMPLETE!")
    print("=" * 60)
    print(f"\nOutput directory: {args.output_dir.absolute()}")
    print("\nGenerated files:")
    print("  annotated_variants.tsv  - Variants passing quality filter")
    if args.include_filtered:
        print("  annotated_variants_all.tsv - All variants (including filtered)")
    print("  annotated_variants.bed  - Color-coded for IGV")
    print("  annotation_summary.txt  - Statistics and summary")
    
    if args.min_quality > 0:
        total = len(annotator.variants)
        passed = len([v for v in annotator.variants if v['passes_quality_filter']])
        print(f"\nQuality filter: QUAL >= {args.min_quality}")
        print(f"  Passed: {passed}/{total} variants ({passed/total*100:.1f}%)")
    
    print("\nTo view in IGV:")
    print(f"  1. Load {vcf_file.name}")
    if bam_file:
        print(f"  2. Load {bam_file.name}")
    print(f"  3. Load annotated_variants.bed")
    print("\nVariants color-coded by location:")
    print("  🔴 Red = Coding (check these first!)")
    print("  🟠 Orange = Exon")
    print("  🔵 Blue = Gene")
    print("  ⚫ Gray = Intergenic or low quality")
    print("=" * 60)


if __name__ == '__main__':
    main()