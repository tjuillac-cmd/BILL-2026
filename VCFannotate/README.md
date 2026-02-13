# Variant Annotator - Final Version
## Complete Guide with Folder Input & Quality Filtering

## Overview

The **final enhanced version** combines everything into a simple folder-based workflow with quality filtering.

### What's New

✅ **Folder Input** - Auto-detects VCF, GFF, BAM files  
✅ **Quality Information** - Extracts QUAL, Filter, Depth, Allele Frequency from VCF  
✅ **Quality Filtering** - Filter variants by QUAL score  
✅ **BAM Coverage** - Correlates VCF quality with actual coverage  
✅ **Gene Annotation** - Shows which gene each mutation is in  
✅ **Color Coding** - Visual prioritization in IGV  

## Quick Start

```bash
# 1. Put all files in one folder
mkdir my_analysis
cp variants.vcf my_analysis/
cp genes.gff my_analysis/
cp reads.bam my_analysis/

# 2. Run with quality filter
python variant_annotator_final.py -i my_analysis -o results -q 20

# 3. Open in IGV
Load: results/annotated_variants.bed
```

## Complete Workflow

```
┌──────────────────────────────────────────────────────┐
│ STEP 1: ORGANIZE FILES IN FOLDER                    │
└──────────────────────────────────────────────────────┘

input_data/
├── sample.vcf      ← Variants with QUAL scores
├── genes.gff       ← Gene annotations
└── reads.bam       ← Read alignments (optional)

┌──────────────────────────────────────────────────────┐
│ STEP 2: RUN ANNOTATION WITH QUALITY FILTER          │
└──────────────────────────────────────────────────────┘

$ python variant_annotator_final.py \
    -i input_data \
    -o results \
    -q 20 \              ← Keep only QUAL >= 20
    --include-filtered   ← Also save filtered variants

┌──────────────────────────────────────────────────────┐
│ STEP 3: CHECK OUTPUT                                │
└──────────────────────────────────────────────────────┘

results/
├── annotated_variants.tsv       ← QUAL >= 20 only
├── annotated_variants_all.tsv   ← All variants
├── annotated_variants.bed       ← Color-coded (QUAL >= 20)
└── annotation_summary.txt       ← Statistics

┌──────────────────────────────────────────────────────┐
│ STEP 4: VIEW IN IGV                                 │
└──────────────────────────────────────────────────────┘

Load annotated_variants.bed:
  🔴 Red = High-quality coding variants (priority!)
  🟠 Orange = High-quality exon variants
  🔵 Blue = High-quality gene variants
  ⚫ Gray = Low quality or intergenic
```

## Output Files Explained

### 1. annotated_variants.tsv (Quality-Filtered)

Main output with all quality information:

```
Chromosome  Position  Ref  Alt  QUAL_Score  Filter  VCF_Depth  Allele_Freq  BAM_Coverage  Gene_Name  Location_Type
chr1        100000    A    G    50.0        PASS    100        0.500        150          GENE1      gene
chr17       7577000   G    A    70.0        PASS    150        0.950        180          TP53       coding
```

**Quality Columns:**
- **QUAL_Score** - Phred-scaled quality from VCF (higher = better)
- **Filter** - PASS or filter flags from VCF
- **VCF_Depth** - Read depth from VCF INFO field (DP)
- **Allele_Freq** - Allele frequency from VCF INFO field (AF)
- **BAM_Coverage** - Actual coverage from BAM file (if available)
- **Quality_Flag** - PASS or FAIL based on VCF FILTER column

**Gene Annotation Columns:**
- **Gene_Name** - Gene(s) overlapping variant
- **Gene_Type** - Feature type (CDS, exon, gene)
- **Strand** - Gene strand (+/-)
- **Location_Type** - coding, exon, gene, or intergenic

### 2. annotated_variants_all.tsv (All Variants)

Only created if `--include-filtered` option used.

Contains ALL variants including those below quality threshold.  
Useful for checking what was filtered out.

### 3. annotated_variants.bed (Color-Coded for IGV)

BED file with RGB colors based on:
1. **Location type** (coding > exon > gene > intergenic)
2. **Quality score** (higher QUAL = brighter color)

```
chr1   99999  100000  A>G|GENE1|Q50     600  +  99999  100000  0,128,255   (Blue)
chr17  7576999 7577000 G>A|TP53|Q70    1000  -  7576999 7577000 255,0,0     (Red)
```

Names include quality: `A>G|GENE1|Q50` means QUAL=50

### 4. annotation_summary.txt (Statistics)

```
VARIANT COUNTS:
Total Variants: 247
Quality Filter: QUAL >= 20
Passed Filter: 189 (76.5%)
Filtered Out: 58 (23.5%)

QUALITY STATISTICS:
Variants with QUAL score: 247
Average QUAL: 45.3
Min QUAL: 10.0
Max QUAL: 99.0
FILTER=PASS: 230
FILTER=FAIL: 17

VARIANT LOCATIONS (Quality Filtered):
coding        :   95 (50.3%)  ← High priority!
exon          :   57 (30.2%)
gene          :   25 (13.2%)
intergenic    :   12 ( 6.3%)

TOP GENES BY VARIANT COUNT:
TP53          :  15 variants
BRCA1         :  12 variants
```

## Command-Line Options

```bash
python variant_annotator_final.py [OPTIONS]

Required:
  -i, --input-dir DIR       Input directory with VCF, GFF, BAM files

Optional:
  -o, --output-dir DIR      Output directory (default: annotation_output)
  -q, --min-quality FLOAT   Minimum QUAL score (default: 0 = no filter)
  --include-filtered        Also save filtered variants to _all.tsv file
```

## Quality Filtering Guide

### Choosing Quality Threshold

| Quality | Recommended Use |
|---------|----------------|
| `-q 0` | No filtering (default) - see all variants |
| `-q 10` | Very permissive - include most variants |
| `-q 20` | Standard - balanced sensitivity/specificity |
| `-q 30` | Stringent - high confidence only |
| `-q 50` | Very stringent - very high confidence |


## Understanding Quality Scores

### VCF Quality Fields

**QUAL (column 6):**
- Phred-scaled quality score
- QUAL=20 → 99% confidence
- QUAL=30 → 99.9% confidence
- QUAL=50 → 99.999% confidence

**FILTER (column 7):**
- PASS = passed all filters
- Other values indicate specific failures
- Tool marks as Quality_Flag=FAIL if not PASS

**INFO DP (Depth):**
- Total read depth at position
- Extracted from INFO field
- Higher = more evidence

**INFO AF (Allele Frequency):**
- Frequency of alternate allele
- AF=0.5 → heterozygous
- AF=1.0 → homozygous alternate
- AF=0.1-0.3 → possible somatic/subclonal

### BAM vs VCF Coverage

**VCF Depth (from INFO):**
- Reported by variant caller
- May include filtered reads

**BAM Coverage (from samtools):**
- Actual read count at position
- Independent verification

**Comparison helps identify:**
- Discrepancies (caller issues)
- Read quality problems
- Coverage artifacts

## Integration with IGV

### Loading Files

```
1. Open IGV
2. Genomes → Select hg38 (or your genome)
3. File → Load from File:
   - Load reads.bam (alignments track)
   - Load variants.vcf (variant track)
   - Load annotated_variants.bed (color-coded track)
```

### Visual Interpretation

```
IGV Display:
┌────────────────────────────────────────────┐
│ Coverage Track (BAM)                       │
│ ▓▓▓▓▓▓▓▓▓▓▓▓▓▓ (150x)                     │
├────────────────────────────────────────────┤
│ Variants (VCF)                             │
│      ▼ G>A (QUAL=70)                       │
├────────────────────────────────────────────┤
│ Annotated (BED)                            │
│      🔴 G>A|TP53|Q70  ← Red = coding!      │
├────────────────────────────────────────────┤
│ Genes (GFF)                                │
│ ━━━━[████]━━━━ TP53                        │
│      exon/CDS                              │
└────────────────────────────────────────────┘

Interpretation:
✓ High coverage (150x)
✓ High quality (QUAL=70)
✓ In coding sequence (red)
✓ In TP53 (important tumor suppressor)
→ HIGH PRIORITY for validation
```

### Color Coding Priority

```
Priority Order in IGV:

1. 🔴 Red (Coding, High QUAL)
   └─ Most important - affects protein

2. 🟠 Orange (Exon, High QUAL)
   └─ Important - may affect splicing

3. 🔵 Blue (Gene, High QUAL)
   └─ Moderate - intronic regions

4. ⚫ Gray (Low QUAL or Intergenic)
   └─ Lower priority
```

## Quality Control Workflow

### Step 1: Run Without Filter

```bash
# See all variants first
python variant_annotator_final.py -i data -o qc_all -q 0
```

Check `annotation_summary.txt`:
- What's the QUAL range?
- How many FILTER=FAIL?
- Distribution by location?

### Step 2: Choose Filter

Based on summary, choose appropriate `-q` value:
- Low quality variants (avg QUAL < 30) → use `-q 10`
- Good quality (avg QUAL 30-50) → use `-q 20`
- High quality (avg QUAL > 50) → use `-q 30`

### Step 3: Run With Filter

```bash
python variant_annotator_final.py -i data -o filtered -q 20 --include-filtered
```

### Step 4: Compare

```bash
# Count variants by quality
echo "All variants:"
wc -l qc_all/annotated_variants.tsv

echo "Quality filtered:"
wc -l filtered/annotated_variants.tsv

# Check what was filtered
grep -v "^#" filtered/annotated_variants_all.tsv | \
  awk '$5 < 20' > low_quality.txt
```

### Step 5: Validate

In IGV, spot-check:
- High QUAL variants (should look good)
- Low QUAL variants (check if justified)
- Borderline variants (around threshold)

## Advanced Usage

### Extract Specific Subsets

```bash
# Only coding variants with QUAL >= 30
grep "coding" annotated_variants.tsv | awk '$5 >= 30' > coding_high_qual.txt

# Only TP53 variants
grep "TP53" annotated_variants.tsv > TP53_variants.txt

# High coverage variants (VCF depth > 100)
awk '$7 > 100' annotated_variants.tsv > high_coverage.txt
```

### Batch Processing

```bash
#!/bin/bash
# Process multiple samples

for sample_dir in sample_*/; do
    sample=$(basename $sample_dir)
    
    python variant_annotator_final.py \
        -i "$sample_dir" \
        -o "results/${sample}" \
        -q 20
    
    # Extract coding variants
    grep "coding" "results/${sample}/annotated_variants.tsv" \
        > "results/${sample}_coding.txt"
done

# Summarize all samples
echo "Sample,Total,Coding,TP53" > summary.csv
for sample_dir in results/*/; do
    sample=$(basename $sample_dir)
    total=$(grep -v "^#" "${sample_dir}/annotated_variants.tsv" | wc -l)
    coding=$(grep "coding" "${sample_dir}/annotated_variants.tsv" | wc -l)
    tp53=$(grep "TP53" "${sample_dir}/annotated_variants.tsv" | wc -l)
    echo "${sample},${total},${coding},${tp53}" >> summary.csv
done
```

### Statistical Analysis

```python
import pandas as pd

# Load annotations
df = pd.read_csv('annotated_variants.tsv', sep='\t', comment='#')

# Quality distribution
print(df['QUAL_Score'].describe())

# Variants by location and quality
print(df.groupby(['Location_Type']).agg({
    'QUAL_Score': ['count', 'mean', 'std'],
    'BAM_Coverage': 'mean'
}))

# Filter high-confidence coding variants
high_conf_coding = df[
    (df['Location_Type'] == 'coding') & 
    (df['QUAL_Score'] >= 30) & 
    (df['BAM_Coverage'] >= 50)
]

print(f"High-confidence coding variants: {len(high_conf_coding)}")
```

## Troubleshooting

### No Quality Scores

**Problem:** QUAL_Score column shows N/A

**Cause:** VCF file doesn't have QUAL column filled

**Solution:** Check variant caller settings to output quality scores

### All Variants Filtered

**Problem:** No variants pass quality filter

**Cause:** Quality threshold too high for your data

**Solution:**
```bash
# Check quality distribution
python variant_annotator_final.py -i data -o check -q 0
grep "Average QUAL" check/annotation_summary.txt

# Use lower threshold
python variant_annotator_final.py -i data -o results -q 10
```

### Coverage Mismatch

**Problem:** VCF_Depth and BAM_Coverage very different

**Possible causes:**
- Variant caller filtered some reads
- BAM file doesn't match VCF
- Different counting methods

**Investigation:**
```bash
# Check specific position in IGV
# Compare VCF INFO vs actual BAM coverage
# Look for read quality issues
```

## Summary

### Input Requirements
✓ One folder with VCF, GFF (BAM optional)

### Key Features
✓ Auto-detection of files  
✓ Quality filtering by QUAL score  
✓ Gene annotation from GFF  
✓ Coverage from BAM  
✓ Color-coded IGV visualization  

### Output Files
✓ Annotated table with quality info  
✓ Color-coded BED for IGV  
✓ Summary statistics  
✓ Optional: filtered variants file  

### Typical Workflow
1. Put files in folder
2. Run with `-q 20` (or appropriate threshold)
3. Check summary statistics
4. Load BED in IGV
5. Focus on red (coding) variants

The tool combines quality information from VCF, gene locations from GFF, and coverage from BAM to give you a complete picture of each mutation!
